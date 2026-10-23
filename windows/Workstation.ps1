[CmdletBinding()]
param(
    [ValidateSet('Install','Start','Check','Update')][string]$Action='Start',
    [string]$InstallRoot='',
    [string]$Distro='',
    [ValidateRange(1024,65535)][int]$Port=8767,
    [switch]$FunctionsOnly
)

function ConvertTo-NativeArgument {
    param([AllowEmptyString()][string]$Value)
    # WSL parses management options from the raw command line and retains quotes.
    # Keep flags/simple values bare; use CRT quoting only where required.
    if ($Value.Length -gt 0 -and $Value -notmatch '[\s"]') { return $Value }
    $escaped=$Value -replace '(\\*)"','$1$1\"'
    $escaped=$escaped -replace '(\\+)$','$1$1'
    return '"'+$escaped+'"'
}

function Get-WslExecutable {
    # Get-Command -CommandType Application can return multiple PATH matches.
    # Passing their .Source array to FileName joins them into an invalid filename.
    if ($env:OS -eq 'Windows_NT' -and $env:SystemRoot) {
        $folders=@('System32')
        if ([Environment]::Is64BitOperatingSystem -and ![Environment]::Is64BitProcess) { $folders=@('Sysnative','System32') }
        foreach ($folder in $folders) {
            $candidate=Join-Path $env:SystemRoot ($folder+'\wsl.exe')
            if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
        }
    }
    foreach ($command in @(Get-Command wsl.exe -CommandType Application -ErrorAction SilentlyContinue)) {
        if ($command.Path -and (Test-Path -LiteralPath $command.Path -PathType Leaf)) { return [string]$command.Path }
    }
    throw '未找到可执行的 wsl.exe。请检查 Windows 的 WSL 安装。'
}

function ConvertFrom-WslBytes {
    param([byte[]]$Bytes)
    # wsl --list emits UTF-16LE on some Windows versions; Linux tools use UTF-8.
    if ($Bytes.Length -ge 2 -and (($Bytes[0] -eq 255 -and $Bytes[1] -eq 254) -or $Bytes[1] -eq 0)) {
        return [Text.Encoding]::Unicode.GetString($Bytes).TrimStart([char]0xfeff)
    }
    return [Text.Encoding]::UTF8.GetString($Bytes).TrimStart([char]0xfeff)
}

function Invoke-Wsl {
    param([string[]]$WslArguments,[switch]$Capture,[switch]$AllowFailure,[string]$Step='WSL 调用')
    # WSL flags remain bare; arguments after --exec follow Windows CRT quoting.
    # Never merge stderr into paths/JSON: WSL startup warnings also use stderr.
    $info=New-Object Diagnostics.ProcessStartInfo
    $info.FileName=Get-WslExecutable
    $info.Arguments=($WslArguments | ForEach-Object { ConvertTo-NativeArgument $_ }) -join ' '
    $info.UseShellExecute=$false
    # Native WSL warnings otherwise use UTF-16 while Linux stdout/stderr use UTF-8.
    $info.EnvironmentVariables['WSL_UTF8']='1'
    $info.RedirectStandardOutput=[bool]$Capture
    $info.RedirectStandardError=[bool]$Capture
    $process=New-Object Diagnostics.Process
    $process.StartInfo=$info
    $stdout=$null; $stderr=$null; $errorText=''; $outputText=''
    try {
        try {
            if (!$process.Start()) { throw '无法启动 wsl.exe。' }
        } catch {
            $detail=$_.Exception.GetBaseException().Message
            if ($script:DiagnosticFile) {
                [pscustomobject]@{time=[DateTime]::Now.ToString('o');step=$Step;executable=$info.FileName;arguments=$WslArguments;start_error=$detail} | ConvertTo-Json -Compress | Add-Content -LiteralPath $script:DiagnosticFile -Encoding UTF8
            }
            throw "$Step 无法启动。可执行文件：$($info.FileName)。原因：$detail。诊断日志：$script:DiagnosticFile"
        }
        if ($Capture) {
            $stdout=New-Object IO.MemoryStream
            $stderr=New-Object IO.MemoryStream
            # Drain both pipes concurrently, including output larger than a pipe buffer.
            $outTask=$process.StandardOutput.BaseStream.CopyToAsync($stdout)
            $errTask=$process.StandardError.BaseStream.CopyToAsync($stderr)
        }
        $process.WaitForExit()
        $code=$process.ExitCode
        if ($Capture) {
            $null=$outTask.GetAwaiter().GetResult(); $null=$errTask.GetAwaiter().GetResult()
            $outputText=ConvertFrom-WslBytes $stdout.ToArray()
            $errorText=ConvertFrom-WslBytes $stderr.ToArray()
        }
    } finally {
        if ($stdout) { $stdout.Dispose() }; if ($stderr) { $stderr.Dispose() }; $process.Dispose()
    }
    if ($script:DiagnosticFile) {
        $record=[pscustomobject]@{time=[DateTime]::Now.ToString('o');step=$Step;executable=$info.FileName;command_line=$info.Arguments;arguments=$WslArguments;exit_code=$code;stderr=$errorText}
        $record | ConvertTo-Json -Compress | Add-Content -LiteralPath $script:DiagnosticFile -Encoding UTF8
    }
    if ($Capture -and $code -eq 0 -and $errorText.Trim() -and !$AllowFailure) {
        if ($errorText -match '(?i)localhost.*(?:代理|proxy)') {
            if (!$script:ProxyNoticeShown) {
                Write-Host 'WSL 提示 Windows localhost 代理未转发；提示已与路径/配置输出分离。如后续下载失败，请查看日志中的网络报错。' -ForegroundColor Yellow
                $script:ProxyNoticeShown=$true
            }
        } else { Write-Host $errorText.Trim() -ForegroundColor Yellow }
    }
    if (!$AllowFailure -and $code -ne 0) {
        if ($Capture -and $errorText.Trim()) { Write-Host $errorText.Trim() -ForegroundColor Red }
        throw "$Step 失败（退出码 $code）。WSL 命令：$($info.Arguments)。诊断日志：$script:DiagnosticFile"
    }
    if ($Capture) { return [pscustomobject]@{ExitCode=$code;Text=$outputText;Stderr=$errorText} }
}

function Select-Ubuntu {
    param([string]$Requested)
    $result=Invoke-Wsl -WslArguments @('--list','--quiet') -Capture
    $names=@(($result.Text -replace "`0",'') -split '\r?\n' | ForEach-Object { $_.Trim() } | Where-Object { $_ })
    if ($Requested) {
        if ($names -notcontains $Requested) { throw "找不到 Ubuntu 发行版 $Requested。已安装：$($names -join ', ')。" }
        return $Requested
    }
    foreach ($preferred in @('Ubuntu-24.04','Ubuntu')) { if ($names -contains $preferred) { return $preferred } }
    $ubuntu=@($names | Where-Object { $_ -match '^Ubuntu(?:-|$)' })
    if ($ubuntu.Count -gt 0) { return $ubuntu[0] }
    throw '没有找到 Ubuntu。请先在 Windows 安装 WSL2 Ubuntu 24.04，然后再次运行部署。'
}

function Get-WslPath {
    param([string]$WindowsPath,[string]$SelectedDistro)
    if ($WindowsPath -match '^\\\\wsl(?:\.localhost|\$)\\([^\\]+)\\(.*)$') {
        if ($Matches[1] -ne $SelectedDistro) { throw '源码位于另一 Ubuntu 发行版，请先将源码解压到普通 Windows 文件夹。' }
        return '/'+($Matches[2] -replace '\\','/')
    }
    $result=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','wslpath','-u',$WindowsPath.TrimEnd('\')) -Capture -Step '转换 Windows 路径'
    $path=$result.Text.Trim()
    if (!$path.StartsWith('/') -or $path -match '[\r\n\x00]') { throw "WSL 返回了无效路径，尚未调用分析脚本。Windows 路径：$WindowsPath。诊断日志：$script:DiagnosticFile" }
    return $path
}

function Ensure-Docker {
    param([string]$SelectedDistro)
    $probe=@('-d',$SelectedDistro,'--exec','docker','info','--format','{{.ServerVersion}}')
    $result=Invoke-Wsl -WslArguments $probe -Capture -AllowFailure
    if ($result.ExitCode -eq 0) { Write-Host "Docker 已就绪：$($result.Text.Trim())"; return }
    $desktop=$null
    foreach ($candidate in @((Join-Path $env:ProgramFiles 'Docker\Docker\Docker Desktop.exe'),
                              (Join-Path $env:LOCALAPPDATA 'Programs\DockerDesktop\Docker Desktop.exe'))) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) { $desktop=$candidate; break }
    }
    if ($desktop) {
        Write-Host '正在打开 Docker Desktop，请等待引擎启动…'
        Start-Process -FilePath $desktop
        for ($i=0;$i -lt 60;$i++) {
            $result=Invoke-Wsl -WslArguments $probe -Capture -AllowFailure
            if ($result.ExitCode -eq 0) { Write-Host 'Docker 已就绪。'; return }
            if ($i % 10 -eq 0) { Write-Host "等待 Docker…如未开启 WSL Integration，请在设置中启用 $SelectedDistro。" }
            Start-Sleep -Seconds 2
        }
        throw "Ubuntu 无法连接 Docker。请在 Docker Desktop → Settings → Resources → WSL Integration 中启用 $SelectedDistro，再双击启动。"
    }
    $engine=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','sh','-c','command -v dockerd >/dev/null 2>&1') -Capture -AllowFailure
    if ($engine.ExitCode -eq 0) {
        Write-Host '检测到 Ubuntu Docker Engine。启动服务可能需要你本机 Ubuntu 的 sudo 密码。'
        Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','sudo','service','docker','start')
        $result=Invoke-Wsl -WslArguments $probe -Capture -AllowFailure
        if ($result.ExitCode -eq 0) { Write-Host 'Docker Engine 已就绪。'; return }
        throw 'Docker Engine 已启动，但当前 Ubuntu 用户无法访问它。请检查 Docker 用户组权限，或改用 Docker Desktop 的 WSL Integration。'
    }
    throw '没有可用的 Docker。请安装并打开 Docker Desktop：https://www.docker.com/products/docker-desktop/，开启 Ubuntu 的 WSL Integration 后再次部署。'
}

function Invoke-Helper {
    param($Config,[string]$Operation,[switch]$Capture)
    $native=@('-d',$Config.distro,'--exec','/bin/bash','--',($Config.program_wsl+'/windows/wsl_entry.sh'),
              $Operation,$Config.root_wsl,$Config.program_wsl,$Config.venv_wsl,[string]$Config.port,$Config.root_windows)
    return Invoke-Wsl -WslArguments $native -Capture:$Capture -Step ("工作站操作："+$Operation)
}

function Get-ServiceState {
    param($Config,[string]$Operation='probe')
    $result=Invoke-Wsl -WslArguments @('-d',$Config.distro,'--exec','python3',
        ($Config.program_wsl+'/windows/station_runtime.py'),$Operation,$Config.root_wsl,[string]$Config.port) -Capture -Step ("检查服务："+$Operation)
    return $result.Text | ConvertFrom-Json
}

function Assert-WindowsService {
    param($Config)
    for ($i=0;$i -lt 5;$i++) {
        try {
            $health=Invoke-RestMethod -Uri "http://127.0.0.1:$($Config.port)/api/health" -TimeoutSec 4
            if ($health.service_id -ne 'mouse-bulk-rnaseq' -or $health.run_dir -ne ($Config.root_wsl+'/分析结果')) {
                throw 'Windows 本地端口返回的服务不属于此工作站。'
            }
            if (!$health.docker_image) { throw '工作站的分析镜像缺失，请重新部署。' }
            return
        } catch { if ($i -eq 4) { throw ('Ubuntu 服务已启动，但 Windows 检查未通过：'+$_.Exception.Message+' 请检查 WSL 的 localhost 转发和端口设置。') } }
        Start-Sleep -Seconds 1
    }
}

function Start-Station {
    param($Config)
    Ensure-Docker $Config.distro
    $state=Get-ServiceState $Config
    if ($state.state -eq 'occupied') { throw "端口 $($Config.port) 已被其他程序占用；没有停止其他程序。重新部署时可用 -Port 指定新端口。" }
    if ($state.state -eq 'running') {
        if ($state.health.program_dir -ne $Config.program_wsl) { throw '旧版本工作站仍在运行。请先关闭旧终端，或运行更新工作站.bat。' }
        Write-Host '工作站已经运行，打开现有页面。'
    } else {
        Invoke-Helper $Config 'check'
        $native=@('-d',$Config.distro,'--exec','/bin/bash','--',($Config.program_wsl+'/windows/wsl_entry.sh'),
                  'start',$Config.root_wsl,$Config.program_wsl,$Config.venv_wsl,[string]$Config.port,$Config.root_windows)
        $arguments=($native | ForEach-Object { ConvertTo-NativeArgument $_ }) -join ' '
        Start-Process -FilePath (Get-WslExecutable) -ArgumentList $arguments
        $null=Get-ServiceState $Config 'wait'
    }
    Assert-WindowsService $Config
    $url="http://localhost:$($Config.port)"
    Write-Host "工作站位置：$($Config.root_windows)"
    Write-Host "分析结果：$($Config.runs_windows)"
    Write-Host "网页：$url。请保持新打开的 Ubuntu 工作站终端运行。"
    Start-Process $url
}

function New-StationShortcut {
    param([string]$Root)
    try {
        $desktop=[Environment]::GetFolderPath('Desktop')
        $path=Join-Path $desktop '转录组分析工作站.lnk'
        $shell=New-Object -ComObject WScript.Shell
        $shortcut=$shell.CreateShortcut($path)
        $arguments='/c ""'+(Join-Path $Root '启动工作站.bat')+'""'
        if (Test-Path -LiteralPath $path) {
            if ($shortcut.Arguments -ne $arguments -or $shortcut.TargetPath -ne $env:ComSpec) {
                Write-Warning '桌面已有同名的其他快捷方式，已保留它。请从 D 盘双击启动工作站.bat。'
            } else { Write-Host '桌面已有本工作站快捷方式，已保留。' }
            return
        }
        $shortcut.TargetPath=$env:ComSpec
        $shortcut.Arguments=$arguments
        $shortcut.WorkingDirectory=$Root
        $shortcut.Description='本地小鼠 bulk RNA-seq 转录组分析工作站'
        $shortcut.Save()
        Write-Host "桌面快捷方式已创建：$path"
    } catch { Write-Warning '无法创建桌面快捷方式；仍可从 D 盘双击启动工作站.bat。' }
}

function Get-UpdateSource {
    param([string]$Root)
    $folder=Join-Path $Root ('日志\updates\'+[Guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $folder -Force | Out-Null
    $archive=Join-Path $folder 'source.zip'
    Write-Host '下载 GitHub 最新源码；保留当前版本和全部分析结果…'
    [Net.ServicePointManager]::SecurityProtocol=[Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
    Invoke-WebRequest -Uri 'https://github.com/kangweihou522-cmyk/H-/archive/refs/heads/main.zip' -OutFile $archive -UseBasicParsing
    (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash | Set-Content -LiteralPath (Join-Path $folder 'source.sha256') -Encoding ASCII
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip=[IO.Compression.ZipFile]::OpenRead($archive)
    try {
        foreach ($entry in $zip.Entries) {
            $name=$entry.FullName -replace '\\','/'
            if ($name -match '^/|^[A-Za-z]:|(^|/)\.\.(/|$)') { throw '源码 ZIP 含越界路径，未解压。' }
        }
    } finally { $zip.Dispose() }
    $expanded=Join-Path $folder 'source'
    Expand-Archive -LiteralPath $archive -DestinationPath $expanded
    $sources=@(Get-ChildItem -LiteralPath $expanded -Directory | Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'windows\Workstation.ps1') })
    if ($sources.Count -ne 1) { throw '下载的源码不完整，当前版本和数据保持原样。' }
    return $sources[0].FullName
}

function Install-Station {
    param([string]$Source,[string]$Root,[string]$SelectedDistro,[int]$SelectedPort,$Existing)
    Ensure-Docker $SelectedDistro
    $sourceWsl=Get-WslPath $Source $SelectedDistro
    $rootWsl=Get-WslPath $Root $SelectedDistro
    Write-Host "Ubuntu：$SelectedDistro；源码：$sourceWsl；工作站：$rootWsl"
    Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','test','-f',($sourceWsl+'/windows/wsl_entry.sh')) -Step '检查源码脚本路径'
    Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','/bin/bash','--',($sourceWsl+'/windows/wsl_entry.sh'),'prerequisites',$rootWsl,$sourceWsl) -Step '准备 Ubuntu 依赖'
    if ($Existing) { $null=Get-ServiceState $Existing 'stop' }
    $prepared=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','python3',($sourceWsl+'/windows/deploy.py'),
        'prepare',$sourceWsl,$rootWsl,$Root,$SelectedDistro,[string]$SelectedPort) -Capture -Step '准备工作站版本'
    $proposal=$prepared.Text | ConvertFrom-Json
    Invoke-Helper $proposal.config 'install'
    $activated=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','python3',
        ($proposal.config.program_wsl+'/windows/deploy.py'),'activate',$rootWsl,$proposal.plan_file) -Capture -Step '激活工作站版本'
    $config=$activated.Text | ConvertFrom-Json
    # Keep an Ubuntu command entry as well as the Windows desktop entry.
    $null=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'--exec','python3',
        ($config.program_wsl+'/windows/ubuntu_shortcut.py'),$config.root_wsl) -Capture
    New-StationShortcut $Root
    Write-Host '部署完成。程序、原始数据目录、参考文件目录和分析结果均在 D 盘工作站文件夹。'
    Write-Host '真实小鼠 GO 可导入离线注释；联网注释安装是后续可选步骤，不影响打开工作站。'
    Start-Station $config
}

function Invoke-Workstation {
    $ErrorActionPreference='Stop'
    [Console]::OutputEncoding=New-Object Text.UTF8Encoding($false)
    $OutputEncoding=[Console]::OutputEncoding
    Write-Host ('WSL：'+(Get-WslExecutable))
    $root=$InstallRoot
    if (!$root) { $root='D:\转录组分析工作站' }
    $root=[IO.Path]::GetFullPath($root)
    if (!(Test-Path -LiteralPath ([IO.Path]::GetPathRoot($root)))) { throw '没有找到目标磁盘，请检查 D 盘是否存在。' }
    $logDir=Join-Path $root '日志'
    New-Item -ItemType Directory -Path $logDir -Force | Out-Null
    $script:DiagnosticFile=Join-Path $logDir 'wsl-launcher.log'
    $file=Join-Path $root '配置\workstation.json'
    $existing=$null
    if (Test-Path -LiteralPath $file) {
        $existing=Get-Content -LiteralPath $file -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($existing.service_id -ne 'mouse-bulk-rnaseq' -or $existing.schema -ne 1) { throw '目标目录包含其他配置，未覆盖，请选择独立安装目录。' }
    }
    $requested=$Distro
    if (!$requested -and $existing) { $requested=$existing.distro }
    if (!$requested -and $env:WSL_DISTRO) { $requested=$env:WSL_DISTRO }
    $selected=Select-Ubuntu $requested
    $selectedPort=$Port
    if (!$script:ExplicitPort -and $existing) { $selectedPort=$existing.port }
    $hasher=[Security.Cryptography.SHA256]::Create()
    try { $key=([BitConverter]::ToString($hasher.ComputeHash([Text.Encoding]::UTF8.GetBytes($root.ToLowerInvariant())))).Replace('-','').Substring(0,16) } finally { $hasher.Dispose() }
    $mutex=New-Object Threading.Mutex($false,('Local\MouseRnaseqStation-'+$key))
    if (!$mutex.WaitOne(0)) { $mutex.Dispose(); throw '工作站正在部署或启动，请等待当前操作完成后重试。' }
    try {
    if ($Action -in @('Install','Update')) {
        if ($Action -eq 'Update' -and !$existing) { throw '尚未部署，请先运行部署到D盘.bat。' }
        New-Item -ItemType Directory -Path $root -Force | Out-Null
        $source=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
        if ($Action -eq 'Update') { $source=Get-UpdateSource $root }
        Install-Station $source $root $selected $selectedPort $existing
    } else {
        if (!$existing) { throw '尚未部署到 D 盘。请解压新版源码并双击部署到D盘.bat。' }
        $existing.distro=$selected
        if ($Action -eq 'Check') { Ensure-Docker $selected; Invoke-Helper $existing 'check'; Write-Host '环境检查通过。' }
        else { Start-Station $existing }
    }
    } finally { $mutex.ReleaseMutex(); $mutex.Dispose() }
}

$script:ExplicitPort=$PSBoundParameters.ContainsKey('Port')
if (!$FunctionsOnly) {
    try { Invoke-Workstation; exit 0 }
    catch { Write-Host $_.Exception.Message -ForegroundColor Red; exit 1 }
}

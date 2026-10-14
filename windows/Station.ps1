param([ValidateSet('Start','Check','Update')][string]$Action='Start')
$ErrorActionPreference='Stop'
try {
    $root=$PSScriptRoot
    $config=Get-Content -LiteralPath (Join-Path $root '配置\workstation.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($config.service_id -ne 'mouse-bulk-rnaseq' -or $config.schema -ne 1) { throw '此目录不是本工作站的安装目录。' }
    $driver=Join-Path $config.program_windows 'windows\Workstation.ps1'
    if (!(Test-Path -LiteralPath $driver -PathType Leaf)) { throw '程序文件缺失，请重新下载新版源码并运行部署到D盘.bat。' }
    & $driver -Action $Action -InstallRoot $root -Distro $config.distro -Port $config.port
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} catch {
    Write-Host $_.Exception.Message -ForegroundColor Red
    exit 1
}

param([string]$Driver=(Join-Path $PSScriptRoot '../windows/Workstation.ps1'))
# Linux-only executable shim. Exercises actual Process pipes/arguments, not a function mock.
# It does not constitute a Windows/WSL VM test.
$ErrorActionPreference='Stop'
. $Driver -FunctionsOnly
$script:count=0
function Assert-True { param($Value,[string]$Message); if (!$Value) { throw $Message }; $script:count++ }
$shim=(Get-Command wsl.exe -CommandType Application).Source
Assert-True ($shim -match 'rnaseq-wsl-native') 'Test must use an isolated fake wsl executable'
$result=Invoke-Wsl -WslArguments @('--list','--quiet') -Capture
Assert-True (@($result).Count -eq 1) 'No async task objects leak into function output'
Assert-True ($result.Text -ceq "Ubuntu-24.04`nUbuntu`n") 'UTF16LE list decoded'
Assert-True ((Select-Ubuntu '') -eq 'Ubuntu-24.04') 'Native distro selection despite stderr warning'
$path="D:\转录组 数据\O'Brien"
Assert-True ((Get-WslPath $path 'Ubuntu-24.04') -ceq '/mnt/d/转录组 数据') 'Proxy warning cannot enter path'
$wanted=@('', '转录组 工作站', 'a"b', 'C:\RNA seq\', "O'Brien", '$(touch WRONG)', '--')
$result=Invoke-Wsl -WslArguments (@('echo-argv')+$wanted) -Capture
$actual=@($result.Text | ConvertFrom-Json)
Assert-True ($actual.Count -eq $wanted.Count) 'Argument count retained'
for ($i=0;$i -lt $wanted.Count;$i++) { Assert-True ($actual[$i] -ceq $wanted[$i]) "Argument $i retained verbatim" }
$fixture=Join-Path (Split-Path $shim) '下载 源码/转录组 工作站/windows'
$null=New-Item -ItemType Directory -Path $fixture -Force
$entry=Join-Path $fixture 'wsl_entry.sh'
$body=@'
#!/usr/bin/env bash
python3 -c 'import json,sys;print(json.dumps(sys.argv[1:],ensure_ascii=False))' "$@"
'@
[IO.File]::WriteAllText($entry,$body,(New-Object Text.UTF8Encoding($false)))
$config=[pscustomobject]@{distro='Ubuntu-24.04';program_wsl=(Split-Path $fixture);root_wsl='/mnt/d/转录组 工作站';venv_wsl='/home/user/Python env';port=8767;root_windows="D:\转录组 工作站 O'Brien"}
$result=Invoke-Helper $config 'check' -Capture
$actual=@($result.Text | ConvertFrom-Json)
$expected=@('check',$config.root_wsl,$config.program_wsl,$config.venv_wsl,'8767',$config.root_windows)
Assert-True ($actual.Count -eq $expected.Count) 'Real bash script receives all helper arguments'
for ($i=0;$i -lt $expected.Count;$i++) { Assert-True ($actual[$i] -ceq $expected[$i]) "Bash argument $i retained verbatim" }
$script:DiagnosticFile=Join-Path ([IO.Path]::GetTempPath()) ('rnaseq-wsl-native-'+[Guid]::NewGuid().ToString('N')+'.log')
try {
    $result=Invoke-Wsl -WslArguments @('large-pipes') -Capture -AllowFailure
    Assert-True ($result.Text.Length -eq 200000 -and $result.Stderr.Length -eq 200000) 'Both pipes drained without deadlock'
    $result=Invoke-Wsl -WslArguments @('fail') -Capture -AllowFailure
    Assert-True ($result.ExitCode -eq 126 -and $result.Stderr -match 'cannot execute binary') 'Exit126 and stderr retained'
    $rejected=$false
    try { Invoke-Wsl -WslArguments @('fail') -Capture -Step '部署原生调用测试' } catch {
        $rejected=$_.Exception.Message -match '部署原生调用测试.*126.*fail.*rnaseq-wsl-native'
    }
    Assert-True $rejected 'Failure names phase, arguments and diagnostic file'
    $records=@(Get-Content -LiteralPath $script:DiagnosticFile | ForEach-Object { $_ | ConvertFrom-Json })
    Assert-True ($records.Count -eq 3 -and $records[-1].exit_code -eq 126) 'Structured diagnostic records'
    # Noncaptured invocations must inherit terminal handles for sudo/password prompts.
    $null=Invoke-Wsl -WslArguments @('interactive') -Step '终端继承测试'
} finally { Remove-Item -LiteralPath $script:DiagnosticFile -ErrorAction SilentlyContinue }
"PASS: $script:count real executable/pipe checks on Linux (no Windows or WSL VM)"

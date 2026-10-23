param([string]$Driver=(Join-Path $PSScriptRoot '..\windows\Workstation.ps1'))
$ErrorActionPreference='Stop'
. $Driver -FunctionsOnly
$script:count=0
function Assert-True { param($Value,[string]$Message); if (!$Value) { throw $Message }; $script:count++ }
Assert-True ((ConvertTo-NativeArgument '--list') -ceq '--list') 'WSL management flag must be bare'
Assert-True ((ConvertTo-NativeArgument 'Ubuntu-24.04') -ceq 'Ubuntu-24.04') 'WSL distribution must be bare'
Assert-True ((ConvertTo-NativeArgument '') -ceq '""') 'Empty argument'
Assert-True ((ConvertTo-NativeArgument 'a"b') -ceq '"a\"b"') 'Embedded quote'
Assert-True ((ConvertTo-NativeArgument 'C:\RNA seq\') -ceq '"C:\RNA seq\\"') 'Trailing backslash'
Assert-True ((ConvertTo-NativeArgument '转录组 工作站') -ceq '"转录组 工作站"') 'Unicode whitespace'

function Invoke-Wsl {
    param($WslArguments,[switch]$Capture,[switch]$AllowFailure,[string]$Step)
    $script:lastWslArgs=@($WslArguments)
    $text=''
    if ($WslArguments -contains '--list') { $text="Ubuntu-24.04`0`nUbuntu`0`ndocker-desktop`0" }
    elseif ($WslArguments -contains 'wslpath') { $text='/mnt/d/转录组 工作站' }
    if ($script:invalidPath) { $text="wsl: warning`n/mnt/d/root" }
    return [pscustomobject]@{Text=$text;Stderr='wsl: proxy warning';ExitCode=0}
}
Assert-True ((Select-Ubuntu '') -eq 'Ubuntu-24.04') 'Ubuntu selection'
Assert-True ((Select-Ubuntu 'Ubuntu') -eq 'Ubuntu') 'Explicit Ubuntu selection'
$rejected=$false;try { Select-Ubuntu 'Missing' } catch {$rejected=$true}
Assert-True $rejected 'Missing Ubuntu rejected'
Assert-True ((Get-WslPath '\\wsl.localhost\Ubuntu-24.04\home\user' 'Ubuntu-24.04') -eq '/home/user') 'WSL UNC conversion'
$path="D:\RNAseq data\O'Brien"
$null=Get-WslPath $path 'Ubuntu-24.04'
Assert-True ($script:lastWslArgs[-1] -ceq $path) 'Paths passed as a single native argument'
$script:invalidPath=$true;$rejected=$false
try { Get-WslPath $path 'Ubuntu-24.04' } catch {$rejected=$true}
Assert-True $rejected 'Multiline/contaminated path rejected before bash'
$script:invalidPath=$false

function Get-WslExecutable { return 'wsl.exe' }
function Ensure-Docker { param($SelectedDistro) }
function Invoke-Helper { param($Config,$Operation,[switch]$Capture); $script:helper=$Operation }
function Start-Process { param($FilePath,$ArgumentList); $script:starts+=@([pscustomobject]@{File=$FilePath;Arguments=$ArgumentList}) }
function Get-ServiceState { param($Config,$Operation='probe'); if($Operation -eq 'wait'){return @{state='running'}}; return $script:state }
function Invoke-RestMethod { param($Uri,$TimeoutSec); return @{service_id='mouse-bulk-rnaseq';run_dir='/mnt/d/转录组 工作站/分析结果';docker_image=$true} }
$config=[pscustomobject]@{distro='Ubuntu-24.04';program_wsl="/mnt/d/转录组 工作站/程序/version";root_wsl='/mnt/d/转录组 工作站';venv_wsl='/home/user/runtime/venv';port=8767;root_windows='D:\转录组 工作站';runs_windows='D:\转录组 工作站\分析结果'}
$script:state=@{state='occupied'};$script:starts=@();$rejected=$false
try { Start-Station $config } catch {$rejected=$true}
Assert-True ($rejected -and $script:starts.Count -eq 0) 'Other service never opened or stopped'
$script:state=@{state='running';health=@{program_dir=$config.program_wsl}};$script:starts=@()
Start-Station $config
Assert-True ($script:starts.Count -eq 1 -and $script:starts[0].File -eq 'http://localhost:8767') 'Existing service opens browser only'
$script:state=@{state='offline'};$script:starts=@()
Start-Station $config
Assert-True ($script:starts.Count -eq 2 -and $script:starts[0].File -eq 'wsl.exe') 'New service plus browser'
Assert-True ($script:starts[0].Arguments.Contains('"/mnt/d/转录组 工作站/程序/version/windows/wsl_entry.sh"')) 'Quoted Chinese/space script path'
Assert-True ($script:starts[0].Arguments.Contains('-d Ubuntu-24.04 --exec /bin/bash -- ')) 'Absolute bash with option terminator'
"PASS: $script:count PowerShell launcher checks (WSL/native process calls mocked; no Windows desktop VM)"

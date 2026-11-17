$ErrorActionPreference = 'Stop'
$helper = Join-Path $PSScriptRoot 'transfer_results.ps1'
$tokens=$null; $parseErrors=$null
$ast = [System.Management.Automation.Language.Parser]::ParseFile($helper,[ref]$tokens,[ref]$parseErrors)
if ($parseErrors.Count) { throw 'Receiver syntax errors.' }
# Load only declarations. Never execute the receiver loop or a network call.
$addType = $ast.FindAll({param($n) $n -is [System.Management.Automation.Language.CommandAst] -and $n.GetCommandName() -eq 'Add-Type'},$false)
Invoke-Expression $addType[0].Extent.Text
foreach ($fn in $ast.FindAll({param($n) $n -is [System.Management.Automation.Language.FunctionDefinitionAst]},$false)) {
    Invoke-Expression $fn.Extent.Text
}
$Utf8 = New-Object System.Text.UTF8Encoding($false)
$TaskId = 'sigma_parameters_structured_v1'
$testRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot ('transfer-test-' + [Guid]::NewGuid().ToString('N'))))
$workspace = [IO.Path]::GetFullPath($PSScriptRoot).TrimEnd('\') + '\'
if (-not $testRoot.StartsWith($workspace,[StringComparison]::OrdinalIgnoreCase)) { throw 'Test directory outside intended workspace.' }
$null = Assert-SafePath $testRoot
$null = New-Item -ItemType Directory -Path $testRoot
$LocalRoot = $testRoot
$EventLog = Join-Path $testRoot 'events.jsonl'
$DestinationRoot = Join-Path $testRoot 'destination'
$null = New-Item -ItemType Directory -Path $DestinationRoot
$model = 'Fornax999999'
$modelRoot = Join-Path $DestinationRoot $model
$null = New-Item -ItemType Directory -Path $modelRoot
$MockBytes = $null
function Invoke-TestScp {
    # Only this in-memory payload is written; no SSH/SCP executable is run.
    [IO.File]::WriteAllBytes([string]$args[-1],$script:MockBytes)
    $global:LASTEXITCODE = 0
}
$Scp = 'Invoke-TestScp'
$Server = 'unused-offline-test'
try {
    $state = Join-Path $testRoot 'state.json'
    Write-AtomicJson $state ([ordered]@{generation=1})
    Write-AtomicJson $state ([ordered]@{generation=2})
    if (([IO.File]::ReadAllText($state,$Utf8) | ConvertFrom-Json).generation -ne 2) { throw 'Two atomic JSON writes failed.' }

    $unknownName = $model + '_unknown.csv'
    $unknownPath = Join-Path $modelRoot $unknownName
    [IO.File]::WriteAllText($unknownPath,'preserve existing unrelated content',$Utf8)
    $before = (Get-FileHash -LiteralPath $unknownPath -Algorithm SHA256).Hash
    $unknownDelivery = [PSCustomObject]@{run_id=$model;model_directory=$model;files=@([PSCustomObject]@{path=('/travail/xhuang/test/'+$unknownName);name=$unknownName;bytes=3;sha256=('0'*64)})}
    $rejected=$false
    try { $null=Receive-Delivery $unknownDelivery } catch { $rejected=$_.Exception.Message -like '*without matching task ownership*' }
    if (-not $rejected -or (Get-FileHash -LiteralPath $unknownPath -Algorithm SHA256).Hash -ne $before) { throw 'Unowned existing-file protection failed.' }

    $ownedName = $model + '_owned.csv'
    foreach ($content in @('generation-one','generation-two')) {
        $script:MockBytes = $Utf8.GetBytes($content)
        $hash = [Security.Cryptography.SHA256]::Create()
        try { $digest=([BitConverter]::ToString($hash.ComputeHash($script:MockBytes))).Replace('-','').ToLowerInvariant() } finally { $hash.Dispose() }
        $delivery = [PSCustomObject]@{run_id=$model;model_directory=$model;files=@([PSCustomObject]@{path=('/travail/xhuang/test/'+$ownedName);name=$ownedName;bytes=$script:MockBytes.Length;sha256=$digest})}
        if ((Receive-Delivery $delivery) -ne 1) { throw 'Owned-file delivery did not report one file.' }
        if ([IO.File]::ReadAllText((Join-Path $modelRoot $ownedName),$Utf8) -ne $content) { throw 'Owned-file content differs after atomic publication.' }
    }
    if (@(Get-ChildItem -LiteralPath $testRoot -Recurse -File | Where-Object {$_.Name -like '*.partial-*' -or $_.Name -like '*.transfer-*'}).Count) { throw 'Successful publication left a temporary file.' }
    Write-Output 'PASS: two atomic JSON writes, two owned-file publications, unowned-file preservation; no network or real-model data accessed.'
} finally {
    # Remove only this unique test fixture, after absolute-path, ancestor, volume
    # and FileId checks. No wildcard or shell-built deletion is used.
    $resolved = Assert-SafePath $testRoot
    if (-not $resolved.StartsWith($workspace,[StringComparison]::OrdinalIgnoreCase)) { throw 'Refusing test cleanup outside workspace.' }
    foreach ($item in Get-ChildItem -LiteralPath $resolved -Recurse -Force) {
        $null = Assert-SafePath $item.FullName
        $null = [SigmaTransferIdentity]::Read($item.FullName)
    }
    $null = [SigmaTransferIdentity]::Read($resolved)
    Remove-Item -LiteralPath $resolved -Recurse -Force
}

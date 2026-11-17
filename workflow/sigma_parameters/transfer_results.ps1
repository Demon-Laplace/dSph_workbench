[CmdletBinding()]
param(
    [string]$Server = 'mesopsl-codex',
    [string]$ServerRoot = '/travail/xhuang/workspace/jobs/sigma-parameters',
    [Parameter(Mandatory=$true)]
    [ValidateNotNullOrEmpty()]
    [string]$DestinationRoot,
    [int]$PollSeconds = 300,
    [switch]$Once
)

# This helper watches only published data-delivery manifests. It never queries
# Slurm, simulation processes, or job logs, and never deletes server data.
$ErrorActionPreference = 'Stop'
$TaskId = 'sigma_parameters_structured_v1'
$LocalRoot = $PSScriptRoot
$StateFile = Join-Path $LocalRoot 'transfer-state.json'
$Utf8 = New-Object System.Text.UTF8Encoding($false)
[Console]::OutputEncoding = $Utf8
$OutputEncoding = $Utf8
$Ssh = (Get-Command ssh.exe -ErrorAction Stop).Source
$Scp = (Get-Command scp.exe -ErrorAction Stop).Source
if ($Server -notmatch '^[A-Za-z0-9_.-]+$' -or
    $ServerRoot -notmatch '^/travail/xhuang/[A-Za-z0-9_./-]+$' -or
    $ServerRoot.Contains('..')) { throw 'Invalid server or task path.' }
if ($PollSeconds -lt 15) { throw 'PollSeconds must be at least 15.' }

Add-Type -TypeDefinition @'
using System;
public static class SigmaTransferAtomic {
    // Windows PowerShell 5.1 coerces a directly supplied $null string argument
    // to String.Empty. Keep the actual null backup name inside typed C#.
    public static void Replace(string source, string destination) {
        System.IO.File.Replace(source, destination, null);
    }
}
'@

function Assert-SafePath([string]$Path) {
    $absolute = [IO.Path]::GetFullPath($Path)
    $cursor = $absolute
    while ($cursor) {
        if (Test-Path -LiteralPath $cursor) {
            $item = Get-Item -LiteralPath $cursor -Force
            if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
                throw "A path component is a reparse point: $cursor"
            }
        }
        $parent = [IO.Directory]::GetParent($cursor)
        if ($null -eq $parent) { break }
        $cursor = $parent.FullName
    }
    return $absolute
}

function Write-AtomicJson([string]$Path, $Value) {
    $pathAbs = Assert-SafePath $Path
    $temporary = "$pathAbs.partial-$([Guid]::NewGuid().ToString('N'))"
    [IO.File]::WriteAllText($temporary, ($Value | ConvertTo-Json -Depth 30), $Utf8)
    $null = Assert-SafePath $temporary
    if (Test-Path -LiteralPath $pathAbs) {
        [SigmaTransferAtomic]::Replace($temporary, $pathAbs)
    } else { [IO.File]::Move($temporary, $pathAbs) }
}

function Read-RemoteJson([string]$Name) {
    if ($Name -notmatch '^[a-z-]+\.json$') { throw 'Invalid manifest name.' }
    $remotePath = "$ServerRoot/$Name"
    $errors = Join-Path $LocalRoot ('ssh-errors-' + [Guid]::NewGuid().ToString('N') + '.txt')
    try {
        $output = @(& $Ssh -o BatchMode=yes -o UpdateHostKeys=no -o ConnectTimeout=20 $Server "if test -f '$remotePath'; then cat '$remotePath'; else exit 44; fi" 2> $errors)
        $code = $LASTEXITCODE
        if ($code -eq 44) { return $null }
        if ($code -ne 0) { throw "Cannot read $Name (SSH exit $code)." }
        $text = $output -join "`n"
        if (-not $text.Trim()) { throw "Empty delivery manifest: $Name" }
        return ($text | ConvertFrom-Json)
    } finally {
        if (Test-Path -LiteralPath $errors) {
            $null = Assert-SafePath $errors
            Remove-Item -LiteralPath $errors -Force
        }
    }
}

function Test-Terminal($Final) {
    if ($null -eq $Final) { return $false }
    $terminals = @('complete','completed','complete_with_failures','failed','cancelled','canceled','terminal','done','success','partial')
    $expected = 2
    if ($Final.PSObject.Properties.Name -contains 'expected_batches') { $expected = [int]$Final.expected_batches }
    if ($expected -lt 1 -or $expected -gt 2) { return $false }
    if ($Final.PSObject.Properties.Name -contains 'completed_batches') {
        if ([int]$Final.completed_batches -ge $expected -and ($Final.terminal -eq $true -or $Final.final -eq $true)) { return $true }
    }
    if ($Final.PSObject.Properties.Name -contains 'batches') {
        $batches = @($Final.batches)
        if ($batches.Count -lt $expected) { return $false }
        foreach ($batch in $batches) {
            $status = [string]$batch.status
            if (-not $status) { $status = [string]$batch.state }
            if ($terminals -notcontains $status.ToLowerInvariant()) { return $false }
        }
        return $true
    }
    return $false
}

function Read-Receipt([string]$RunId) {
    $receiptPath = Join-Path $LocalRoot ($RunId + '_transfer_receipt.json')
    if (-not (Test-Path -LiteralPath $receiptPath)) {
        return [PSCustomObject]@{ task_id = $TaskId; run_id = $RunId; files = @() }
    }
    $receipt = [IO.File]::ReadAllText($receiptPath, $Utf8) | ConvertFrom-Json
    if ($receipt.task_id -ne $TaskId -or $receipt.run_id -ne $RunId) { throw "Unexpected receipt owner for $RunId" }
    return $receipt
}

function Receive-Delivery($Delivery) {
    $runId = [string]$Delivery.run_id
    $model = [string]$Delivery.model_directory
    if ($runId -notmatch '^Fornax[0-9]+[A-Za-z0-9_.-]*$' -or $runId.Contains('..') -or
        $model -notmatch '^Fornax[0-9]+[A-Za-z0-9_.-]*$' -or $model.Contains('..')) { throw 'Invalid run or model identity.' }
    $rootAbs = Assert-SafePath $DestinationRoot
    if (-not (Test-Path -LiteralPath $rootAbs -PathType Container)) { throw "Destination root unavailable: $rootAbs" }
    $modelRoot = Assert-SafePath (Join-Path $rootAbs $model)
    if ([IO.Directory]::GetParent($modelRoot).FullName -ne $rootAbs.TrimEnd('\')) { throw 'Model directory escaped destination root.' }
    if (-not (Test-Path -LiteralPath $modelRoot)) { $null = New-Item -ItemType Directory -Path $modelRoot }
    $modelPath = Assert-SafePath (Join-Path $modelRoot '_support')
    if (-not (Test-Path -LiteralPath $modelPath)) { $null = New-Item -ItemType Directory -Path $modelPath }
    $receipt = Read-Receipt $runId
    $records = @($receipt.files | ForEach-Object { [PSCustomObject]@{name=$_.name;path=$_.path;source=$_.source;disposition='downloaded'} })
    $files = @($Delivery.files)
    if ($files.Count -eq 0) { throw "No declared result files for $runId" }
    foreach ($file in $files) {
        $name = [string]$file.name
        $remote = [string]$file.path
        if (-not $name.StartsWith($runId, [StringComparison]::Ordinal) -or
            $name -notmatch '^[A-Za-z0-9_.-]+\.(csv|json|jsonl|npz|npy|parquet|txt)$' -or $name.Contains('..') -or
            $remote -notmatch '^/travail/xhuang/[A-Za-z0-9_./-]+$' -or $remote.Contains('..')) { throw "Invalid file declaration for $runId" }
        $target = Assert-SafePath (Join-Path $modelPath $name)
        $previous = @($records | Where-Object { $_.name -eq $name }) | Select-Object -Last 1
        if (Test-Path -LiteralPath $target) {
            if ($null -ne $previous -and $previous.path -eq $target -and $previous.source -eq $remote) {
                continue
            }
            throw "Existing target conflict; preserved file without matching task ownership: $target"
        }
        $temporary = Assert-SafePath (Join-Path $modelPath ('.' + $name + '.transfer-' + [Guid]::NewGuid().ToString('N')))
        try {
            & $Scp -q -o BatchMode=yes -o UpdateHostKeys=no -o ConnectTimeout=20 ($Server + ':' + $remote) $temporary
            if ($LASTEXITCODE -ne 0) { throw "SCP failed for $name (exit $LASTEXITCODE)" }
            $null = Assert-SafePath $temporary
            $null = Assert-SafePath $target
            if (Test-Path -LiteralPath $target) {
                throw "Existing target conflict during transfer; preserved file: $target"
            }
            [IO.File]::Move($temporary, $target)
            $records = @($records | Where-Object { $_.name -ne $name }) + [PSCustomObject]@{ name=$name; path=$target; source=$remote; disposition='downloaded' }
        } finally {
            if (Test-Path -LiteralPath $temporary) {
                $null = Assert-SafePath $temporary
                Remove-Item -LiteralPath $temporary -Force
            }
        }
        # Persist ownership immediately so a later interrupted transfer can resume.
        $receipt = [PSCustomObject]@{task_id=$TaskId;run_id=$runId;model_directory=$model;files=$records;updated_utc=[DateTime]::UtcNow.ToString('o')}
        Write-AtomicJson (Join-Path $LocalRoot ($runId + '_transfer_receipt.json')) $receipt
    }
    $receipt = [PSCustomObject]@{task_id=$TaskId;run_id=$runId;model_directory=$model;files=$records;updated_utc=[DateTime]::UtcNow.ToString('o')}
    Write-AtomicJson (Join-Path $LocalRoot ($runId + '_transfer_receipt.json')) $receipt
    return $files.Count
}

$lockPath = Join-Path $LocalRoot 'transfer.lock'
$null = Assert-SafePath $lockPath
$lock = [IO.File]::Open($lockPath,[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try {
    while ($true) {
        $errors = @(); $delivered = 0; $fileCount = 0; $final = $null
        try {
            $manifest = Read-RemoteJson 'ready-deliveries.json'
            $deliveries = @()
            if ($null -ne $manifest) {
                if ($manifest.PSObject.Properties.Name -contains 'deliveries') { $deliveries = @($manifest.deliveries) }
                else { $deliveries = @($manifest) }
            }
            foreach ($delivery in $deliveries) {
                try { $fileCount += Receive-Delivery $delivery; $delivered++ }
                catch { $errors += [PSCustomObject]@{run_id=$delivery.run_id;error=$_.Exception.Message} }
            }
            $final = Read-RemoteJson 'transfer-final.json'
            $terminal = Test-Terminal $final
            $status = 'waiting_for_data'
            if ($deliveries.Count -gt 0) { $status = 'receiving' }
            if ($errors.Count -gt 0) { $status = 'retry_pending' }
            if ($terminal -and $errors.Count -eq 0 -and $null -ne $manifest) { $status = 'delivered_all_available' }
            Write-AtomicJson $StateFile ([ordered]@{task_id=$TaskId;status=$status;updated_utc=[DateTime]::UtcNow.ToString('o');destination_root=$DestinationRoot;declared_runs=$deliveries.Count;delivered_runs=$delivered;delivered_files=$fileCount;errors=$errors;server_final=$final})
            if ($status -eq 'delivered_all_available') { break }
        } catch {
            Write-AtomicJson $StateFile ([ordered]@{task_id=$TaskId;status='retry_pending';updated_utc=[DateTime]::UtcNow.ToString('o');destination_root=$DestinationRoot;error=$_.Exception.Message})
        }
        if ($Once) { break }
        Start-Sleep -Seconds $PollSeconds
    }
} finally { $lock.Dispose() }

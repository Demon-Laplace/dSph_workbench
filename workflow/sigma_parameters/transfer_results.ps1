[CmdletBinding()]
param(
    [string]$Server = 'mesopsl-codex',
    [string]$ServerRoot = '/travail/xhuang/workspace/jobs/sigma-parameters',
    [string]$DestinationRoot = ('F:\Fornax_data\Fornax' + [char]0x5904 + [char]0x7406),
    [int]$PollSeconds = 300,
    [switch]$Once
)

# This helper watches only published data-delivery manifests. It never queries
# Slurm, simulation processes, or job logs, and never deletes server data.
$ErrorActionPreference = 'Stop'
$TaskId = 'sigma_parameters_structured_v1'
$LocalRoot = $PSScriptRoot
$StateFile = Join-Path $LocalRoot 'transfer-state.json'
$EventLog = Join-Path $LocalRoot 'transfer-events.jsonl'
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
using System.ComponentModel;
using System.Runtime.InteropServices;
using Microsoft.Win32.SafeHandles;
public static class SigmaTransferIdentity {
    [StructLayout(LayoutKind.Sequential)] private struct Info {
        public uint attr, ca, cb, aa, ab, wa, wb, volume, high, low, links, indexHigh, indexLow;
    }
    [DllImport("kernel32.dll", CharSet=CharSet.Unicode, SetLastError=true)]
    private static extern SafeFileHandle CreateFile(string name, uint access, uint share, IntPtr sec, uint mode, uint flags, IntPtr template);
    [DllImport("kernel32.dll", SetLastError=true)]
    private static extern bool GetFileInformationByHandle(SafeFileHandle file, out Info info);
    public static string Read(string path) {
        using (SafeFileHandle file=CreateFile(path,0,7,IntPtr.Zero,3,0x02000000,IntPtr.Zero)) {
            if(file.IsInvalid) throw new Win32Exception(Marshal.GetLastWin32Error());
            Info i; if(!GetFileInformationByHandle(file,out i)) throw new Win32Exception(Marshal.GetLastWin32Error());
            return i.volume.ToString("X8")+":"+i.indexHigh.ToString("X8")+i.indexLow.ToString("X8");
        }
    }
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
    $tempId = [SigmaTransferIdentity]::Read($temporary)
    if (Test-Path -LiteralPath $pathAbs) {
        $targetId = [SigmaTransferIdentity]::Read($pathAbs)
        if ($targetId -eq $tempId) { throw 'Temporary file aliases the JSON destination.' }
        [SigmaTransferIdentity]::Replace($temporary, $pathAbs)
    } else { [IO.File]::Move($temporary, $pathAbs) }
}

function Write-Event([string]$Kind, [string]$Message, $Details = $null) {
    $record = [ordered]@{ time_utc = [DateTime]::UtcNow.ToString('o'); kind = $Kind; message = $Message; details = $Details }
    [IO.File]::AppendAllText($EventLog, (($record | ConvertTo-Json -Compress -Depth 10) + [Environment]::NewLine), $Utf8)
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
            $null = [SigmaTransferIdentity]::Read($errors)
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
    $modelPath = Assert-SafePath (Join-Path $rootAbs $model)
    if ([IO.Directory]::GetParent($modelPath).FullName -ne $rootAbs.TrimEnd('\')) { throw 'Model directory escaped destination root.' }
    if (-not (Test-Path -LiteralPath $modelPath)) { $null = New-Item -ItemType Directory -Path $modelPath }
    $modelVolume = ([SigmaTransferIdentity]::Read($modelPath)).Split(':')[0]
    $receipt = Read-Receipt $runId
    $records = @($receipt.files)
    $files = @($Delivery.files)
    if ($files.Count -eq 0) { throw "No declared result files for $runId" }
    foreach ($file in $files) {
        $name = [string]$file.name
        $remote = [string]$file.path
        $expectedSha = ([string]$file.sha256).ToLowerInvariant()
        $bytes = [long]$file.bytes
        if (-not $name.StartsWith($runId, [StringComparison]::Ordinal) -or
            $name -notmatch '^[A-Za-z0-9_.-]+\.(csv|json|jsonl|npz|npy|parquet|txt)$' -or $name.Contains('..') -or
            $remote -notmatch '^/travail/xhuang/[A-Za-z0-9_./-]+$' -or $remote.Contains('..') -or
            $expectedSha -notmatch '^[0-9a-f]{64}$' -or $bytes -lt 0 -or $bytes -gt 2147483648) { throw "Invalid file declaration for $runId" }
        $target = Assert-SafePath (Join-Path $modelPath $name)
        $previous = @($records | Where-Object { $_.name -eq $name }) | Select-Object -Last 1
        $currentSha = $null
        if (Test-Path -LiteralPath $target) {
            $targetInfo = Get-Item -LiteralPath $target
            $targetId = [SigmaTransferIdentity]::Read($target)
            $currentSha = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant()
            if ($currentSha -eq $expectedSha -and $targetInfo.Length -eq $bytes) {
                $records = @($records | Where-Object { $_.name -ne $name }) + [PSCustomObject]@{ name=$name; path=$target; sha256=$expectedSha; bytes=$bytes; source=$remote; time_utc=[DateTime]::UtcNow.ToString('o'); disposition='verified_existing' }
                continue
            }
            if ($null -eq $previous -or $previous.sha256 -ne $currentSha) { throw "Preserved an existing file without matching task ownership: $target" }
        }
        $temporary = Assert-SafePath (Join-Path $modelPath ('.' + $name + '.transfer-' + [Guid]::NewGuid().ToString('N')))
        try {
            & $Scp -q -o BatchMode=yes -o UpdateHostKeys=no -o ConnectTimeout=20 ($Server + ':' + $remote) $temporary
            if ($LASTEXITCODE -ne 0) { throw "SCP failed for $name (exit $LASTEXITCODE)" }
            $null = Assert-SafePath $temporary
            $tempInfo = Get-Item -LiteralPath $temporary
            $tempId = [SigmaTransferIdentity]::Read($temporary)
            if ($tempId.Split(':')[0] -ne $modelVolume) { throw 'Temporary and target directory volumes differ.' }
            $downloadSha = (Get-FileHash -LiteralPath $temporary -Algorithm SHA256).Hash.ToLowerInvariant()
            if ($tempInfo.Length -ne $bytes -or $downloadSha -ne $expectedSha) { throw "Size/SHA verification failed for $name" }
            $null = Assert-SafePath $target
            if (Test-Path -LiteralPath $target) {
                if ([SigmaTransferIdentity]::Read($target) -ne $targetId -or
                    (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $currentSha) { throw 'Destination changed during transfer.' }
                if ($targetId -eq $tempId) { throw 'Temporary file aliases the destination.' }
                [SigmaTransferIdentity]::Replace($temporary, $target)
            } else { [IO.File]::Move($temporary, $target) }
            if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expectedSha) { throw 'Final file SHA verification failed.' }
            $records = @($records | Where-Object { $_.name -ne $name }) + [PSCustomObject]@{ name=$name; path=$target; sha256=$expectedSha; bytes=$bytes; source=$remote; time_utc=[DateTime]::UtcNow.ToString('o'); disposition='downloaded_verified' }
            Write-Event 'file_delivered' $target @{run_id=$runId;sha256=$expectedSha;bytes=$bytes}
        } finally {
            if (Test-Path -LiteralPath $temporary) {
                $null = Assert-SafePath $temporary
                $null = [SigmaTransferIdentity]::Read($temporary)
                Remove-Item -LiteralPath $temporary -Force
            }
        }
        # Persist ownership immediately so a later interrupted transfer can resume.
        $receipt = [PSCustomObject]@{task_id=$TaskId;run_id=$runId;model_directory=$model;files=$records;updated_utc=[DateTime]::UtcNow.ToString('o')}
        Write-AtomicJson (Join-Path $LocalRoot ($runId + '_transfer_receipt.json')) $receipt
    }
    $receipt = [PSCustomObject]@{task_id=$TaskId;run_id=$runId;model_directory=$model;files=$records;updated_utc=[DateTime]::UtcNow.ToString('o')}
    Write-AtomicJson (Join-Path $LocalRoot ($runId + '_transfer_receipt.json')) $receipt
    $modelReceipt = Join-Path $modelPath ($runId + '_sigma_parameters_delivery.json')
    if (Test-Path -LiteralPath $modelReceipt) {
        $oldReceipt = [IO.File]::ReadAllText((Assert-SafePath $modelReceipt),$Utf8) | ConvertFrom-Json
        if ($oldReceipt.task_id -ne $TaskId -or $oldReceipt.run_id -ne $runId) { throw 'Preserved an unrelated model receipt.' }
    }
    Write-AtomicJson $modelReceipt $receipt
    return $files.Count
}

$lockPath = Join-Path $LocalRoot 'transfer.lock'
$null = Assert-SafePath $lockPath
$lock = [IO.File]::Open($lockPath,[IO.FileMode]::OpenOrCreate,[IO.FileAccess]::ReadWrite,[IO.FileShare]::None)
try {
    Write-Event 'started' 'Structured-data receiver started; no scheduler queries.' @{destination=$DestinationRoot;server=$Server}
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
                catch { $errors += [PSCustomObject]@{run_id=$delivery.run_id;error=$_.Exception.Message}; Write-Event 'retry_pending' $_.Exception.Message @{run_id=$delivery.run_id} }
            }
            $final = Read-RemoteJson 'transfer-final.json'
            $terminal = Test-Terminal $final
            $status = 'waiting_for_data'
            if ($deliveries.Count -gt 0) { $status = 'receiving' }
            if ($errors.Count -gt 0) { $status = 'retry_pending' }
            if ($terminal -and $errors.Count -eq 0 -and $null -ne $manifest) { $status = 'delivered_all_available' }
            Write-AtomicJson $StateFile ([ordered]@{task_id=$TaskId;status=$status;updated_utc=[DateTime]::UtcNow.ToString('o');destination_root=$DestinationRoot;declared_runs=$deliveries.Count;delivered_runs=$delivered;verified_files=$fileCount;errors=$errors;server_final=$final})
            if ($status -eq 'delivered_all_available') { Write-Event 'finished' 'All declared files were verified locally; server computation is terminal.'; break }
        } catch {
            Write-Event 'retry_pending' $_.Exception.Message
            Write-AtomicJson $StateFile ([ordered]@{task_id=$TaskId;status='retry_pending';updated_utc=[DateTime]::UtcNow.ToString('o');destination_root=$DestinationRoot;error=$_.Exception.Message})
        }
        if ($Once) { break }
        Start-Sleep -Seconds $PollSeconds
    }
} finally { $lock.Dispose() }

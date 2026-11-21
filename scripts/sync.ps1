#Requires -Version 5.1
<#
.SYNOPSIS
  Deterministic git pull / lab-only push / branch push. Do not wrap in extra git commands.
.EXAMPLE
  .\scripts\sync.ps1 pull
  .\scripts/sync.ps1 push-lab
  .\scripts/sync.ps1 push-branch
#>
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [ValidateSet("pull", "push-lab", "push-branch")]
    [string]$Command
)

$ErrorActionPreference = "Stop"
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $Root

function Assert-NoAviStaged {
    $staged = git diff --cached --name-only
    foreach ($p in $staged) {
        if ($p -match '\.avi$') {
            Write-Error "Refusing to continue: staged AVI path '$p'. Unstage it. Raw recordings must not be committed."
            exit 1
        }
    }
}

if ($Command -eq "pull") {
    git fetch
    git pull --ff-only
    if ($LASTEXITCODE -ne 0) {
        Write-Error "git pull --ff-only failed. Stop. Do not rebase or force. Resolve with Sai."
        exit 1
    }
    Write-Host "pull ok (fast-forward only)"
    exit 0
}

if ($Command -eq "push-branch") {
    $branch = "tracking-fixes"
    $remoteRef = git ls-remote --heads origin "refs/heads/$branch"
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Could not check whether origin/$branch exists. Stop."
        exit 1
    }
    if ($remoteRef) {
        Write-Error "Refusing to continue: origin/$branch already exists."
        exit 1
    }

    git push origin "HEAD:refs/heads/$branch"
    if ($LASTEXITCODE -ne 0) {
        Write-Error "git push failed. Stop. Do not force."
        exit 1
    }
    Write-Host "push-branch ok: origin/$branch"
    exit 0
}

# push-lab: notes and evidence only — not lab/tools Python
$paths = @()
Get-ChildItem -Path (Join-Path $Root "lab") -Filter "*.md" -File -ErrorAction SilentlyContinue | ForEach-Object { $paths += $_.FullName }
$evidence = Join-Path $Root "lab\evidence"
if (Test-Path $evidence) { $paths += $evidence }
$writeup = Join-Path $Root "writeup"
if (Test-Path $writeup) { $paths += $writeup }

if ($paths.Count -eq 0) {
    Write-Host "nothing to stage for push-lab"
    exit 0
}

git add -- $paths
Assert-NoAviStaged

$status = git status --porcelain
if (-not $status) {
    Write-Host "no lab-note changes to push"
    exit 0
}

# If the only staged files are still empty, skip
$cached = git diff --cached --name-only
if (-not $cached) {
    Write-Host "no lab-note changes staged"
    exit 0
}

git commit -m "lab: sync notebook"
if ($LASTEXITCODE -ne 0) {
    Write-Error "commit failed"
    exit 1
}
git push
if ($LASTEXITCODE -ne 0) {
    Write-Error "git push failed. Stop. Do not force."
    exit 1
}
Write-Host "push-lab ok"

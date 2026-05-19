param(
    [string]$Root = ".",
    [string]$OutFile = "project_dump.txt"
)

$ExcludeDirs = @(
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "build",
    "dist",
    "results",
    "result",
    "output",
    "data",
    ".idea",
    ".vscode"
)

$ExcludeFiles = @(
    "*.pyc",
    "*.pyo",
    "*.pyd",
    "*.log",
    "*.zip",
    "*.png",
    "*.jpg",
    "*.jpeg",
    "*.gif",
    "*.npy",
    "*.npz",
    "*.csv",
    "*.tsv",
    "*.pdb",
    "*.mol2",
    "*.sdf",
    "*.json"
)

$AllowedTextExt = @(
    ".py", ".txt", ".md", ".toml", ".cfg", ".ini", ".yaml", ".yml",
    ".json", ".js", ".ts", ".sh", ".ps1", ".bat", ".cmd", ".rst"
)

$Root = (Resolve-Path $Root).Path
$OutFile = Join-Path $Root $OutFile

function Should-ExcludeDir($path) {
    $name = Split-Path $path -Leaf
    return $ExcludeDirs -contains $name
}

function Should-ExcludeFile($path) {
    $name = Split-Path $path -Leaf
    foreach ($pat in $ExcludeFiles) {
        if ($name -like $pat) { return $true }
    }
    return $false
}

function Get-TreeLines($dir, $prefix = "") {
    $items = Get-ChildItem -LiteralPath $dir | Where-Object {
        if ($_.PSIsContainer) { -not (Should-ExcludeDir $_.FullName) } else { -not (Should-ExcludeFile $_.FullName) }
    } | Sort-Object @{Expression={$_.PSIsContainer}; Descending=$true}, Name

    $count = $items.Count
    for ($i = 0; $i -lt $count; $i++) {
        $item = $items[$i]
        $isLast = ($i -eq $count - 1)
        $branch = if ($isLast) { "└── " } else { "├── " }
        $nextPrefix = $prefix + (if ($isLast) { "    " } else { "│   " })

        if ($item.PSIsContainer) {
            "$prefix$branch$($item.Name)"
            Get-TreeLines $item.FullName $nextPrefix
        } else {
            "$prefix$branch$($item.Name)"
        }
    }
}

$sb = New-Object System.Text.StringBuilder
[void]$sb.AppendLine("PROJECT TREE")
[void]$sb.AppendLine("ROOT: $Root")
[void]$sb.AppendLine("")

$rootName = Split-Path $Root -Leaf
[void]$sb.AppendLine($rootName)
Get-TreeLines $Root "" | ForEach-Object { [void]$sb.AppendLine($_) }

[void]$sb.AppendLine("")
[void]$sb.AppendLine("=" * 80)
[void]$sb.AppendLine("FILE CONTENTS")
[void]$sb.AppendLine("=" * 80)
[void]$sb.AppendLine("")

$files = Get-ChildItem -Path $Root -Recurse -File | Where-Object {
    -not (Should-ExcludeFile $_.FullName) -and
    -not ($ExcludeDirs | ForEach-Object { $_ }) -contains $_.Directory.Name -and
    ($AllowedTextExt -contains $_.Extension.ToLower())
}

foreach ($file in $files) {
    [void]$sb.AppendLine("### FILE: $($file.FullName)")
    try {
        $content = Get-Content -LiteralPath $file.FullName -Raw -ErrorAction Stop
        [void]$sb.AppendLine($content)
    } catch {
        [void]$sb.AppendLine("[ERROR READING FILE: $($_.Exception.Message)]")
    }
    [void]$sb.AppendLine("")
    [void]$sb.AppendLine("---")
    [void]$sb.AppendLine("")
}

Set-Content -LiteralPath $OutFile -Value $sb.ToString() -Encoding UTF8
Write-Host "Written: $OutFile"
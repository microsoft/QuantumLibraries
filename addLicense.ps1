param($target = "C:\Users\martinro\Source\Repos\OSS\Quantum")

$header = "// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
"

function Write-Header ($file)
{
    $content = Get-Content $file
    $filename = Split-Path -Leaf $file
    $fileheader = $header
    Set-Content $file $fileheader
    Add-Content $file $content
}

Get-ChildItem $target -Recurse | ? { $_.Extension -like ".qb" } | % `
{
    Write-Header $_.PSPath.Split(":", 3)[2]
}


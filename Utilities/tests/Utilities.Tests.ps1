# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

# Make sure we're running on Pester 4 or later.
if ((Get-Module Pester | Select-Object -ExpandProperty Version | Select-Object -ExpandProperty Major) -lt 4) {
    Write-Error "Pester version 4 or later required, but $(Get-Module Pester | Select-Object -ExpandProperty Version) found.";
    return;
}

Get-Module Microsoft.Quantum.Utilities | Remove-Module -Force
Import-Module -Force (
    Join-Path `
        $PSScriptRoot `
        "../src/PowerShell/Microsoft.Quantum.Utilities/Microsoft.Quantum.Utilities.psd1"
)

Describe 'GetQdkVersion' {
    It 'Does not throw' {
        { Get-QdkVersion } | Should -Not -Throw;
    }

    It 'Should have all expected versions' {
        $version = Get-QdkVersion;
        $version.Key.Count | Should -Be 6;
    }
}

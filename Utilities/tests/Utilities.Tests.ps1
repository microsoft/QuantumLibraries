# Copyright (c) Microsoft Corporation.
# Licensed under the MIT License.

Get-Module Microsoft.Quantum.Utilities | Remove-Module -Force
Import-Module -Force (
    Join-Path `
        $PSScriptRoot `
        "../src/PowerShell/Microsoft.Quantum.Utilities/Microsoft.Quantum.Utilities.psd1"
)

Describe 'GetQdkVersion' {
    It 'Does not throw' {
        { Get-QdkVersion } | Should Not Throw;
    }

    It 'Should have all expected versions' {
        $version = Get-QdkVersion;
        $version.Key.Count | Should Be 6;
    }
}

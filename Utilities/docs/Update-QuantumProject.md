---
external help file: Microsoft.Quantum.Utilities-help.xml
Module Name: Microsoft.Quantum.Utilities
online version:
schema: 2.0.0
---

# Update-QuantumProject

## SYNOPSIS
This cmdlet updates the Quantum SDK version as well as the references to Microsoft Quantum packages
to a specified version of the QDK.

## SYNTAX

```
Update-QuantumProject [[-Version] <String>] [-Path] <FileInfo> [-Revert] [-NoRestore] [-WhatIf] [-Confirm]
 [<CommonParameters>]
```

## DESCRIPTION
This cmdlet updates the Quantum SDK version as well as the references to Microsoft Quantum packages
to a specified version of the QDK.
If no version is specified, then the QDK version of the PowerShell
Microsoft.Quantum.Utilities module is used.
The dependencies of the project are restored, unless 
parameter `-NoRestore` is specified.

## EXAMPLES

### EXAMPLE 1
```
Update-QuantumProject -Path QuantumFourierTransform.csproj
```
Updates a single project to the version of the Quantum Development Kit matching this module.

### EXAMPLE 2
```
Update-QuantumProject -Path QuantumFourierTransform.csproj -Version 0.11.2004.2825
```
Updates a single project to version 0.11.2004.2825 of the Quantum Development Kit.

### EXAMPLE 3
```
Update-QuantumProject -Path QuantumFourierTransform.csproj -Version 0.11.2004.2825 -NoRestore
```
Updates a single project to version 0.11.2004.2825 of the Quantum Development Kit, but does not restore NuGet packages after updating.

### EXAMPLE 4
```
Update-QuantumProject -Path QuantumFourierTransform.csproj -Revert
```
Reverts changes made by a previous call to `Update-QuantumProject`.

## PARAMETERS

### -Version
(Optional) Version of the QDK to update the project to.
If not specified, the version of this PowerShell
module will be used.
To display the version of this modue, use cmdlet 'Get-QdkVersion'.

```yaml
Type: String
Parameter Sets: (All)
Aliases:

Required: False
Position: 1
Default value: $MyInvocation.MyCommand.Module.Version
Accept pipeline input: False
Accept wildcard characters: False
```

### -Path
Path to the Q# project file to update.

```yaml
Type: FileInfo
Parameter Sets: (All)
Aliases:

Required: True
Position: 2
Default value: None
Accept pipeline input: True (ByValue)
Accept wildcard characters: False
```

### -Revert
Restores a the target project to a previous version if available.
Dependency restoration is not performed during a revert.

```yaml
Type: SwitchParameter
Parameter Sets: (All)
Aliases:

Required: False
Position: Named
Default value: False
Accept pipeline input: False
Accept wildcard characters: False
```

### -NoRestore
(Optional) Skips performing the .NET restore step on the project after updating the references.

```yaml
Type: SwitchParameter
Parameter Sets: (All)
Aliases:

Required: False
Position: Named
Default value: False
Accept pipeline input: False
Accept wildcard characters: False
```

### -WhatIf
Shows what would happen if the cmdlet runs.
The cmdlet is not run.

```yaml
Type: SwitchParameter
Parameter Sets: (All)
Aliases: wi

Required: False
Position: Named
Default value: None
Accept pipeline input: False
Accept wildcard characters: False
```

### -Confirm
Prompts you for confirmation before running the cmdlet.

```yaml
Type: SwitchParameter
Parameter Sets: (All)
Aliases: cf

Required: False
Position: Named
Default value: None
Accept pipeline input: False
Accept wildcard characters: False
```

### CommonParameters
This cmdlet supports the common parameters: -Debug, -ErrorAction, -ErrorVariable, -InformationAction, -InformationVariable, -OutVariable, -OutBuffer, -PipelineVariable, -Verbose, -WarningAction, and -WarningVariable. For more information, see [about_CommonParameters](http://go.microsoft.com/fwlink/?LinkID=113216).

## INPUTS

## OUTPUTS

## NOTES

## RELATED LINKS

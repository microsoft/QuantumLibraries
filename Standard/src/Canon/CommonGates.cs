// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Canon
{
    public partial class CX
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args) =>
            args.Value switch
            {
                (Qubit ctrl, Qubit target) => new RuntimeMetadata()
                {
                    Label = "X",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                },
                _ => throw new Exception($"Failed to retrieve runtime metadata for {this.ToString()}."),
            };
    }

    public partial class CY
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args) =>
            args.Value switch
            {
                (Qubit ctrl, Qubit target) => new RuntimeMetadata()
                {
                    Label = "Y",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                },
                _ => throw new Exception($"Failed to retrieve runtime metadata for {this.ToString()}."),
            };
    }

    public partial class CZ
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args) =>
            args.Value switch
            {
                (Qubit ctrl, Qubit target) => new RuntimeMetadata()
                {
                    Label = "Z",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                },
                _ => throw new Exception($"Failed to retrieve runtime metadata for {this.ToString()}."),
            };
    }
}

// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Canon
{
    public partial class CX
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args)
        {
            Debug.Assert(args.Value is ValueTuple<Qubit, Qubit>, $"Failed to retrieve runtime metadata for {this.ToString()}.");

            if (args.Value is ValueTuple<Qubit, Qubit> cxArgs)
            {
                var (ctrl, target) = cxArgs;
                return new RuntimeMetadata()
                {
                    Label = "X",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                };
            }

            return null;
        }
    }

    public partial class CY
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args)
        {
            Debug.Assert(args.Value is ValueTuple<Qubit, Qubit>, $"Failed to retrieve runtime metadata for {this.ToString()}.");

            if (args.Value is ValueTuple<Qubit, Qubit> cyArgs)
            {
                var (ctrl, target) = cyArgs;
                return new RuntimeMetadata()
                {
                    Label = "Y",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                };
            }

            return null;
        }
    }

    public partial class CZ
    {
        public override RuntimeMetadata? GetRuntimeMetadata(IApplyData args)
        {
            Debug.Assert(args.Value is ValueTuple<Qubit, Qubit>, $"Failed to retrieve runtime metadata for {this.ToString()}.");

            if (args.Value is ValueTuple<Qubit, Qubit> czArgs)
            {
                var (ctrl, target) = czArgs;
                return new RuntimeMetadata()
                {
                    Label = "Z",
                    IsControlled = true,
                    Controls = new List<Qubit>() { ctrl },
                    Targets = new List<Qubit>() { target },
                };
            }

            return null;
        }
    }
}

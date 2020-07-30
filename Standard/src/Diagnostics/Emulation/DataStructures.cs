// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Collections.Generic;
using System.Collections.Immutable;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Diagnostics.Emulation
{    
    public class DisplayableUnitaryOperator
    {
        public IList<Qubit>? Qubits { get; set; }
        public NumSharp.NDArray? Data { get; set; }
    }

    public struct CallSites
    {
        public string Subject { get; set; }
        public ImmutableList<ImmutableStack<string>> Sites { get; set; }
    }
}

// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Collections.Generic;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    public class FailureRecord<T>
    {
        [JsonProperty("actual")]
        public T Actual { get; set; }

        [JsonProperty("expected")]
        public T Expected { get; set; }

        [JsonProperty("message")]
        public string Message { get; set; } = "";
    }
    
    public class DisplayableUnitaryOperator
    {
        public IList<Qubit>? Qubits { get; set; }
        public NumSharp.NDArray? Data { get; set; }
    }
}

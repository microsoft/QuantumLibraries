// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Arrays
{
    public partial class ForEach<__T__, __U__>
    {
        public override RuntimeMetadata GetRuntimeMetadata(IApplyData args) =>
            new RuntimeMetadata() { IsComposite = true };
    }
}

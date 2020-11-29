// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;
using System;

namespace Microsoft.Quantum.Arrays
{
    public partial class ForEach<__T__, __U__>
    {
        public override RuntimeMetadata GetRuntimeMetadata(IApplyData args)
        {
            var metadata = base.GetRuntimeMetadata(args);
            if (metadata == null) throw new NullReferenceException($"Null RuntimeMetadata found for {this.ToString()}.");
            metadata.IsComposite = true;
            return metadata;
        }
    }
}

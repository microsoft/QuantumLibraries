// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;
using System;

namespace Microsoft.Quantum.Canon
{
    public partial class ApplyToEach<__T__>
    {
        public override RuntimeMetadata GetRuntimeMetadata(IApplyData args)
        {
            var metadata = base.GetRuntimeMetadata(args);
            if (metadata == null) throw new NullReferenceException($"Null RuntimeMetadata found for {this.ToString()}.");
            metadata.IsComposite = true;
            return metadata;
        }
    }

    public partial class ApplyToEachC<__T__>
    {
        public override RuntimeMetadata GetRuntimeMetadata(IApplyData args)
        {
            var metadata = base.GetRuntimeMetadata(args);
            if (metadata == null) throw new NullReferenceException($"Null RuntimeMetadata found for {this.ToString()}.");
            metadata.IsComposite = true;
            return metadata;
        }
    }

    public partial class ApplyToEachA<__T__>
    {
        public override RuntimeMetadata GetRuntimeMetadata(IApplyData args)
        {
            var metadata = base.GetRuntimeMetadata(args);
            if (metadata == null) throw new NullReferenceException($"Null RuntimeMetadata found for {this.ToString()}.");
            metadata.IsComposite = true;
            return metadata;
        }
    }

    public partial class ApplyToEachCA<__T__>
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

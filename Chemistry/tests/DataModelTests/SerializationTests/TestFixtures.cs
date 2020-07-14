// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Microsoft.Quantum.Chemistry.Broombridge;

using Newtonsoft.Json;

using Xunit;

using static Microsoft.Quantum.Chemistry.Tests.Extensions;

namespace Microsoft.Quantum.Chemistry.Tests
{
    public class FciDumpDataFixture
    {
        public ElectronicStructureProblem problem;

        public FciDumpDataFixture()
        {
            using var reader = File.OpenText(
                Path.Join(
                    "FciDump", "h2_anorccmb.fcidump"
                )
            );
            problem = FciDumpSerializer.Deserialize(reader).Single();
        }
    }

}

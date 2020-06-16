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

    public class IntegrationTests : IClassFixture<FciDumpDataFixture>
    {
        FciDumpDataFixture fixture;
        public IntegrationTests(FciDumpDataFixture fixture)
        {
            this.fixture = fixture;
        }

        [Fact]
        public void RoundtripFcidumpIsCorrect()
        {
            // Start by converting the Fcidump data deserialized
            // in the test fixture out to Broombridge.
            using var writer = new StringWriter();
            BroombridgeSerializer.Serialize(writer, new [] { fixture.problem });

            // Now deserialize back from Broombridge.
            using var reader = new StringReader(writer.ToString());
            var roundtripData = BroombridgeSerializer.Deserialize(reader).Single();

            // Next, check that the resulting problem is the same.
            AssertThat(fixture.problem.IsEqualTo(roundtripData));
        }
    }

}

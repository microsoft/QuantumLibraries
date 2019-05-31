// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

using Microsoft.Quantum.Chemistry.OrbitalIntegrals;

namespace Microsoft.Quantum.Chemistry.Tests
{

    using OrbitalIntegral = OrbitalIntegral;
 
    /// <summary>
    /// These are examples that we may use or are currently using as sample
    /// code online in the chemistry library prose documentation, or in research
    /// manuscripts. If these break, those examples are no longer valid and need updating.
    /// </summary>
    public class DocumentationExamplesTests
    {
        [Fact]
        public void HydrogenOrbitalIntegrals()
        {
            // These orbital integrals represent terms in molecular Hydrogen
            var orbitalIntegrals = new OrbitalIntegral[]
            {
                new OrbitalIntegral(new[] { 0,0 }, -1.252477495),
                new OrbitalIntegral(new[] { 1,1 }, -0.475934275),
                new OrbitalIntegral(new[] { 0,0,0,0 }, 0.674493166),
                new OrbitalIntegral(new[] { 0,1,0,1 }, 0.181287518),
                new OrbitalIntegral(new[] { 0,1,1,0 }, 0.663472101),
                new OrbitalIntegral(new[] { 1,1,1,1 }, 0.697398010)
            };
        }
        
    }
}
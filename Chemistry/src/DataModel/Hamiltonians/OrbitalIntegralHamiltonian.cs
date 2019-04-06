// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;

using System.Runtime.Serialization.Formatters.Binary;
using System.IO;
using System.IO.Compression;
using YamlDotNet.Serialization;
using Microsoft.Extensions.Logging;

namespace Microsoft.Quantum.Chemistry
{

    public class OrbitalIntegralHamiltonian : Hamiltonian<TermType.OrbitalIntegral, OrbitalIntegral>
    {
        public OrbitalIntegralHamiltonian() : base() { }
    }
}
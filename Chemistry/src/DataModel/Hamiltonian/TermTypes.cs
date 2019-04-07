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
    /// <summary>
    /// All Hamiltonian terms must implement this interface.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    public interface HamiltonianTerm <TermClassification>
    {
        TermClassification GetTermType();
    }

    /// <summary>
    /// Class containing a indices to a variety of term categories. 
    /// </summary>
    // Might want to distribute these to the separate term classes.
    public static class TermType
    {
        public enum OrbitalIntegral
        {
            Identity, OneBody, TwoBody
        }

        public enum Fermion
        {
            Identity = 0, PP = 1, PQ = 2, PQQP = 3, PQQR = 4, PQRS = 5
        }

        public enum PauliTerm
        {
            Identity = 0,   // This has contribution from PP and PQQP terms.
            Z = 1,          // This has contribution from PP and PQQP terms.
            ZZ = 2,         // This has contribution from PQQP terms.
            PQ = 3,         // This has contributions from PQ and PQQR terms.
            PQQR = 4,       // This has contributions from PQQR terms.
            v01234 = 5      // This has contributions from PQRS terms.
        }
    }
}
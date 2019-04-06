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
    public static class TermType
    {
        public enum OrbitalIntegral
        {
            Identity, OneBody, TwoBody
        }

        public enum Fermion
        {
            Identity, PP, PQ, PQQP, PQQR, PQRS
        }

        public enum JordanWigner
        {
            Identity, Z, ZZ, PQ, PQQR, PQRS
        }
    }
}
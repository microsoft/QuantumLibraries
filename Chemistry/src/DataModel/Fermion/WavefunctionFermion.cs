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

namespace Microsoft.Quantum.Chemistry.Fermion
{
    /// <summary>
    /// Enum over valid input state types.
    /// </summary>
    public enum StateType
    {
        Default = 0, SingleConfigurational = 1, SparseMultiConfigurational = 2, UnitaryCoupledCluster = 3
    }

    /// <summary>
    /// Representation of a general Fermion Hamiltonian. This Hamiltonian is
    /// assumed to be a linear combination of sequences of creation
    /// and annihilation operators. 
    /// </summary>
    public partial class Wavefunction
    {
        /// <summary>
        /// List of Hamiltonian input states.
        /// </summary>
        public Dictionary<string, InputState> InputStates = new Dictionary<string, InputState>();

        
        /// <summary>
        /// Data structure representing an input state.
        /// </summary>
        public struct InputState
        {
            public StateType type;
            public string Label;
            public Double Energy;
            public ((Double, Double) complexCoeff, FermionTerm term)[] Superposition;
        }
        
    }

}
// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

namespace Microsoft.Quantum.Chemistry
{
    /// <summary>
    /// Generic Hamiltonian class. This is the base class for any Hamiltonians,
    /// which are collections of categorized terms.
    /// </summary>
    /// <typeparam name="TermClassification">Index to categories of terms.</typeparam>
    /// <typeparam name="TermIndexing">Index to individual terms.</typeparam>
    public class GenericQuantumState<TermClassification, TermIndexing> 
        //where TermClassification: IEquatable<TermClassification>
        where TermIndexing: HamiltonianTerm<TermClassification>//, IEquatable<TermIndexing>
    {
        /// <summary>
        /// Container for all terms in a Hamiltonian.
        /// </summary>
        public Dictionary<TermClassification, Dictionary<TermIndexing, double>> terms;

        /// <summary>
        /// Indices to systems (e.g. fermions, qubits, or orbitals) the Hamiltonian acts on.
        /// </summary>
        public HashSet<int> systemIndices = new HashSet<int>();

        /// <summary>
        /// Constructor for empty Hamiltonian.
        /// </summary>
        public GenericQuantumState()
        {
            terms = new Dictionary<TermClassification, Dictionary<TermIndexing, double>>();
        }

        /// <summary>
        /// Constructor for copying a Hamiltonian.
        /// </summary>
        public GenericQuantumState(GenericQuantumState<TermClassification, TermIndexing> hamiltonian)
        {
            terms = hamiltonian.terms;
        }
    }
    
}
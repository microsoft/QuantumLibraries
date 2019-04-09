// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.Linq;
using Microsoft.Quantum.Chemistry.Generic;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{
    public class FermionHamiltonian : Hamiltonian<TermType.Fermion, HermitianFermionTerm, DoubleCoeff>
    {
        public FermionHamiltonian() : base() { }

        /// <summary>
        /// Convert a fermion Hamiltonian object into a serialization-friendly object.
        /// </summary>
        /// <returns>Representation of Hamiltonian in terms of primitive types.</returns>
        // system indices
        public (int[],Dictionary<TermType.Fermion,IEnumerable<(int[],double)>>) SerializationFormat()
        {
            var systemIndices = SystemIndices;
            var terms = Terms
                .ToDictionary(termType => termType.Key, termType => termType.Value
                    .Select(termIndex => (termIndex.Key.Indices().ToArray(), ((double)termIndex.Key.Coefficient) * termIndex.Value))
                    .Select(term => (term.Item1, term.Item2)));

            return (systemIndices.ToArray(), terms);
        }

        /// <summary>
        /// Create a fermion Hamiltonian object from a its serialization.
        /// </summary>
        /// <returns>Deserialized fermion Hamiltonian.</returns>
        // system indices
        public FermionHamiltonian((int[], Dictionary<TermType.Fermion, (int[], double)[]>) serialization)
        {
            SystemIndices = new HashSet<int>(serialization.Item1);
            Terms = serialization.Item2.ToDictionary(termType => termType.Key, termType => termType.Value
                .ToDictionary(kv => new HermitianFermionTerm(kv.Item1), kv => kv.Item2.ToDoubleCoeff()));
        }

        /// <summary>
        /// Method for collecting all distinct system (orbital) indices.
        /// </summary>
        /// <param name="orbitalIntegral">Collate orbital indices from this orbital integral.</param>
        public override void AddToSystemIndices(HermitianFermionTerm index)
        {
            foreach (var idx in index.Sequence)
            {
                SystemIndices.Add(idx.Index);
            }
        }
    }
    
}
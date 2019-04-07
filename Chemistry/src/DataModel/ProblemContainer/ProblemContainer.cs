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
    public class Convenience
    {
        /// <summary>
        /// Sample implementation of end-to-end electronic structure problem simulation. 
        /// </summary>
        /// <param name="filename"></param>
        public void SampleWorkflow(string filename)
        {
            // Deserialize Broombridge from file.
            Broombridge.Current.Data broombridge = Broombridge.DeserializeBroombridge(filename);

            // A single file can contain multiple problem descriptions. Let us pick the first one.
            Broombridge.Current.ProblemDescription problemData = broombridge.ProblemDescriptions.First();

            #region Create electronic structure Hamiltonian
            // Electronic structure Hamiltonians are usually represented compactly by orbital integrals. Let us construct
            // such a Hamiltonian from broombridge.
            OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = problemData.CreateOrbitalIntegralHamiltonian();

            // We can obtain the full fermion Hamiltonian from the more compact orbital integral representation.
            // This transformation requires us to pick a convention for converting a spin-orbital index to a single integer.
            // Let us pick one according to the formula `integer = 2 * orbitalIndex + spinIndex`.
            FermionHamiltonian fermionHamiltonian = orbitalIntegralHamiltonian.ToFermionHamiltonian(SpinOrbital.IndexConvention.UpDown);

            // We target a qubit quantum computer, which requires a Pauli representation of the fermion Hamiltonian.
            // A number of mappings from fermions to qubits are possible. Let us choose the Jordan--Wigner encoding.
            PauliHamiltonian pauliHamiltonian = fermionHamiltonian.ToPauliHamiltonian(PauliTerm.Encoding.JordanWigner);
            #endregion

            #region Create wavefunction Ansatz

            #endregion

            #region Pipe to QSharp and simulate

            #endregion
        }
    }

    public class ProblemContainer
    {
        public class Config { }
        //public readonly Config.IndexConvention.Type IndexConvention;

        // For now, this only has a few options.
        public OrbitalIntegralHamiltonian orbitalIntegralHamiltonian = new OrbitalIntegralHamiltonian();
        public FermionHamiltonian fermionHamiltonian = new FermionHamiltonian();
        //public PauliHamiltonian =

        // For now, do not process input states much.
        public Dictionary<string, InputState> InputStates = new Dictionary<string, InputState>();

        // Additional data
        public Int64 NOrbitals = 0;
        public Int64 NElectrons = 0;
        public string MiscellaneousInformation;

        //public Hamiltonian hamiltonian;

        /// <summary>
        ///     A label for this particular Hamiltonian.
        ///     Can be used to identify the Hamiltonian out of set
        ///     loaded from the same file.
        /// </summary>
        public string Name { get; set; } = "<unknown>";

        /*
                 /// <summary>
        /// Constructor for <see cref="FermionHamiltonian"/>.
        /// </summary>
        /// <param name="fermionTerms">Dictionary of all fermion terms.</param>
        /// <param name="nOrbitals">Total number of distinct orbital indices.</param>
        public FermionHamiltonian(Dictionary<FermionTermType, List<FermionTerm>> fermionTerms, Int64 nOrbitals, Int64 nElectrons = 0, Double energyOffset = 0.0)
        {
            FermionTerms = fermionTerms;
            NOrbitals = nOrbitals;
            NElectrons = nElectrons;
            EnergyOffset = energyOffset;
            SortAndAccumulate();
        }
        */


        /// <summary>
        /// Checks that all terms in a <see cref="FermionHamiltonian"/>
        /// are in canonical order. 
        /// </summary>
        /// <returns>
        /// Returns <c>true</c> if all terms are in canonical order, and <c>false</c><
        /// otherwise.
        /// </returns>
        //public bool VerifyFermionTerms()
        // VerifyHamiltonian();
    }


    /// <summary>
    /// Representation of a general Fermion Hamiltonian. This Hamiltonian is
    /// assumed to be a linear combination of sequences of creation
    /// and annihilation operators. 
    /// </summary>
   // public partial class FermionHamiltonian
   // {
        ///public Config Configuration = Config.Default();

        // This will sort terms and accumulate terms in a canonical format.
        // This format is sorted by:
        // 1) List {1,1,...,0,...} of ones for each `a^\dag` followed by zeros for each `a` 
        // 2) Number of unique spin-orbital indices
        // 3) spin-orbital index for a^\dag terms in ascending order, then spin-orbital index for a term in descending order.

        // A Fermionic term has orbital and spin indices.
        // A Fermionic term e.g. a^{\dag}_{4,0} a^{\dag}_{7,1} a_{2,0} a^{\dag}_{9,1} is
        // represented by {{1,1,0,1},{{4,0},{7,1},{2,0},{9,1}}}
        // This is indexes by QArray<(QArray<Int64>)>

        // Every Fermionic term is classified into a term type based on:
        // 1) The list of ones and zeroes representing creation and annihilation operators.
        // 2) The number of distinct spin-orbital indices.
      /*  
        /// <summary>
        ///  Converts a <see cref="FermionHamiltonian"/> to canonical order. This generates
        ///  new terms and modifies the coefficient as needed.
        /// </summary>
        public void ToCanonicalOrder()
        {
            // For now, this assumes that the creation and annihilation operators are in canonical order.
            // Thus anticommutation here only introduces a minus sign.
            // This will loop over FermionTerm.ToCanonicalOrder();
            // This will then call SortAndAccumulate();
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Assuming that all <see cref="FermionTerm"/>s of a <see cref="FermionHamiltonian"/> are
        /// in canonical order, this sorts all terms according to <see cref="Comparers.FermionTermIComparer"/>
        /// and combines duplicates term types.
        /// </summary>
        public void SortAndAccumulate()
        {
            foreach (var termType in FermionTerms)
            {
                // This assumes that Fermion terms are in canonical order.
                var NOnes = termType.Key.type.Item2.Where(o => o == 1).ToArray().Length;
                termType.Value.Sort(new Comparers.FermionTermIComparer(nOnesIn: NOnes));
            }
            foreach (var termType in new HashSet<FermionTermType>(FermionTerms.Keys))
            {
                FermionTerms[termType] = FermionTerms[termType].AccumulateFermionTerm();
            }
        }

        /// <summary>
        /// String representation of Fermion Hamiltonian.
        /// </summary>
        /// <returns>String representation of Fermion Hamiltonian.</returns>
        public override string ToString()
        {
            var output = "";
            foreach (var termType in FermionTerms)
            {
                output += $"(FermionTermType, Number of entries): ({termType.Key}, {termType.Value.Count()}).\n";
                foreach(var term in termType.Value)
                {
                    output += $"{term}\n";
                }
            }
            return output;
        }

    */
  //  }

}
 
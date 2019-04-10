// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using System.Linq;
using System.Text.RegularExpressions;
using YamlDotNet.Core;
using YamlDotNet.Serialization;

namespace Microsoft.Quantum.Chemistry.Broombridge
{
    /// <summary>
    /// Functionality to validate and parse Broombridge
    /// </summary>
    public static partial class DataStructures
    {

        /// <summary>
        /// Broombridge v0.2 format.
        /// 
        /// Changes from v0.1:
        /// - `format -> version` replaced with `version`.
        /// - `integral_sets` replaced with `problem_description`.
        /// - `metadata` is now optional.
        /// - Initial state suggestion `state` key removed. All values in this key are moved up one level.
        /// </summary>
        #region Broombridge v0.2 format
        public class V0_2
        {
            public static class Strings
            {
                public static string SparseMultiConfigurational = "sparse_multi_configurational";
                public static string UnitaryCoupledCluster = "unitary_coupled_cluster";
                public static string VersionNumber = "0.2";
            }

            // Root of Broombridge data structure
            public struct Data
            {

                [YamlMember(Alias = "$schema", ApplyNamingConventions = false)]
                public string Schema { get; set; }

                [YamlMember(Alias = "format", ApplyNamingConventions = false)]
                public DataStructure.Format Format { get; set; }

                [YamlMember(Alias = "generator", ApplyNamingConventions = false)]
                public DataStructure.Generator Generator { get; set; }

                [YamlMember(Alias = "bibliography", ApplyNamingConventions = false)]
                public List<DataStructure.BibliographyItem> Bibliography { get; set; }

                [YamlMember(Alias = "problem_description", ApplyNamingConventions = false)]
                public List<ProblemDescription> ProblemDescriptions { get; set; }

            }

            public struct ProblemDescription
            {
                [YamlMember(Alias = "metadata", ApplyNamingConventions = false)]
                public Dictionary<string, object> Metadata { get; set; }

                [YamlMember(Alias = "basis_set", ApplyNamingConventions = false)]
                public DataStructure.BasisSet BasisSet { get; set; }

                [YamlMember(Alias = "geometry", ApplyNamingConventions = false)]
                public DataStructure.Geometry Geometry { get; set; }

                [YamlMember(Alias = "coulomb_repulsion", ApplyNamingConventions = false)]
                public DataStructure.SimpleQuantity CoulombRepulsion { get; set; }

                [YamlMember(Alias = "scf_energy", ApplyNamingConventions = false)]
                public DataStructure.SimpleQuantity ScfEnergy { get; set; }

                [YamlMember(Alias = "scf_energy_offset", ApplyNamingConventions = false)]
                public DataStructure.SimpleQuantity ScfEnergyOffset { get; set; }

                [YamlMember(Alias = "fci_energy", ApplyNamingConventions = false)]
                public DataStructure.BoundedQuantity FciEnergy { get; set; }

                [YamlMember(Alias = "n_orbitals", ApplyNamingConventions = false)]
                public int NOrbitals { get; set; }

                [YamlMember(Alias = "n_electrons", ApplyNamingConventions = false)]
                public int NElectrons { get; set; }

                [YamlMember(Alias = "energy_offset", ApplyNamingConventions = false)]
                public DataStructure.SimpleQuantity EnergyOffset { get; set; }

                [YamlMember(Alias = "hamiltonian", ApplyNamingConventions = false)]
                public DataStructure.HamiltonianData Hamiltonian { get; set; }

                // FIXME: actually specify what initial_state_suggestions looks like.
                //[YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
                //public List<Dictionary<string, object>> InitialStateSuggestions { get; set; }

                [YamlMember(Alias = "initial_state_suggestions", ApplyNamingConventions = false)]
                public List<State> InitialStates { get; set; }
            }



            public struct State
            {
                [YamlMember(Alias = "label", ApplyNamingConventions = false)]
                public string Label { get; set; }

                [YamlMember(Alias = "method", ApplyNamingConventions = false)]
                public string Method { get; set; }

                [YamlMember(Alias = "energy", ApplyNamingConventions = false)]
                public DataStructure.SimpleQuantity Energy { get; set; }

                /// <summary>
                /// Sparse multi-configurational data
                /// </summary>
                [YamlMember(Alias = "superposition", ApplyNamingConventions = false)]
                public List<List<object>> Superposition { get; set; }
                
                /// <summary>
                /// Coupled-cluster operator
                /// </summary>
                [YamlMember(Alias = "cluster_operator", ApplyNamingConventions = false)]
                public DataStructure.ClusterOperator ClusterOperator { get; set; }
            }


        }
        #endregion

    }
    
}
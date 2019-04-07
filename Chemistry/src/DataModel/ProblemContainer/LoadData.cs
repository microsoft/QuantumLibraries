// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;
using YamlDotNet;
using static System.Math;
using static System.Convert;

namespace Microsoft.Quantum.Chemistry
{
    
    /// <summary>
    ///     Represents the possible formats that can be used to represent integral
    ///     data sets.
    /// </summary>
    public enum IntegralDataFormat
    {
        Liquid, Broombridge
    }

    /// <summary>
    /// Methods for loading Hamiltonian data from standard formats
    /// into a <see cref="FermionHamiltonian"/>.
    /// </summary>
    public partial class LoadData
    {
        /// <summary>
        /// Only spin 1/2 electrons currently supported.
        /// </summary>
        private const Int64 NSpins = 2;

    }


    public partial class FermionHamiltonian
    {
        /// <summary>
        ///      Loads a Hamiltonian from integral data represented
        ///      in Broombridge format.
        ///      Please see the <a href="https://docs.microsoft.com/quantum/libraries/chemistry/schema/spec">
        ///      for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      An instance of <see cref="FermionHamiltonian"/> representing the
        ///      data contained in <paramref name="filename"/>.
        /// </returns>
        public static IEnumerable<FermionHamiltonian> LoadFromBroombridge(string filename, Chemistry.Config configuration)
        {
            var broombridgeData = Broombridge.Deserialize.Source(filename);
            IEnumerable<BroombridgeTyped> broombridgeDataTyped = broombridgeData.ProblemDescription.Select(o => new BroombridgeTyped(o, configuration.indexConvention));
            return broombridgeDataTyped.Select(o => LoadData.LoadIntegralData(o, configuration));
        }

        public static IEnumerable<FermionHamiltonian> LoadFromBroombridge(string filename)
        {
            var configuration = Microsoft.Quantum.Chemistry.Config.Default();
            return LoadFromBroombridge(filename, configuration);
        }

    }


    public partial class LoadData
    {
        internal static FermionHamiltonian LoadIntegralData(BroombridgeTyped schemaInstance)
        {
            return LoadIntegralData(schemaInstance, Config.Default());
        }

        internal static FermionHamiltonian LoadIntegralData(BroombridgeTyped schemaInstance, Config configuration)
        {
            var hamiltonian = new FermionHamiltonian(configuration)
            {

                NOrbitals = schemaInstance.NOrbitals,
                NElectrons = schemaInstance.NElectrons,
                EnergyOffset = schemaInstance.IdentityTerm,
                // TODO: optionally look at metadata to decide a better label.
                Name = $"hamiltonian"
            };

            foreach (var orbitalIntegral in schemaInstance.OneBodyTerms)
            {
                if (Abs(orbitalIntegral.Coefficient) >= configuration.truncationThreshold)
                {
                    hamiltonian.AddFermionTerm(orbitalIntegral);
                }
            }

            foreach (var orbitalIntegral in schemaInstance.TwoBodyTerms)
            {
                if (Abs(orbitalIntegral.Coefficient) >= configuration.truncationThreshold)
                {
                    hamiltonian.AddFermionTerm(orbitalIntegral);
                }
            }

            hamiltonian.SortAndAccumulate();

            hamiltonian.InputStates = schemaInstance.InitialStates;

            return hamiltonian;
        }



        
    }
}
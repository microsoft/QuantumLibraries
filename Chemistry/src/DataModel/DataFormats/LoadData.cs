// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;
using YamlDotNet;

namespace Microsoft.Quantum.Chemistry
{
    
    /// <summary>
    ///     Represents the possible formats that can be used to represent integral
    ///     data sets.
    /// </summary>
    public enum IntegralDataFormat
    {
        Liquid, YAML
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
        ///      Please see the <a href="https://raw.githubusercontent.com/Microsoft/Quantum/master/Chemistry/Schema/broombridge-0.1.schema.json">
        ///      for further details about the
        ///      format parsed by this method.
        /// </summary>
        /// <param name="filename">The name of the file to be loaded.</param>
        /// <returns>
        ///      An instance of <see cref="FermionHamiltonian"/> representing the
        ///      data contained in <paramref name="filename"/>.
        /// </returns>
        public static IEnumerable<FermionHamiltonian> LoadFromBroombridge(string filename)
        {
            var broombridgeData = Broombridge.Deserialize.Source(filename);

            return LoadData.LoadIntegralData(broombridgeData);
        }

    }


    public partial class LoadData
    {
        internal static IEnumerable<FermionHamiltonian> LoadIntegralData(Broombridge.V0_2.Data schemaInstance, Double threshold = 1e-8)
        {
            return schemaInstance.ProblemDescription.Select(
                (hamiltonianData, index) =>
                {
                    var hamiltonian = new FermionHamiltonian()
                    {

                        NOrbitals = hamiltonianData.NOrbitals,
                        EnergyOffset = 0 * hamiltonianData.EnergyOffset.Value + hamiltonianData.CoulombRepulsion.Value,
                        // TODO: optionally look at metadata to decide a better label.
                        Name = $"hamiltonian_{index}"
                    };

                    var fileHIJTerms = new Dictionary<Int64[], Double>(new Extensions.IntArrayIEqualityComparer());
                    var fileHIJKLTerms = new Dictionary<Int64[], Double>(new Extensions.IntArrayIEqualityComparer());

                    fileHIJTerms = hamiltonianData.Hamiltonian.OneElectronIntegrals.Values.ToDictionary(entry => entry.Item1, entry => entry.Item2);
                    fileHIJKLTerms = hamiltonianData.Hamiltonian.TwoElectronIntegrals.Values.ToDictionary(entry => entry.Item1, entry => entry.Item2);
                    hamiltonian.NOrbitals = hamiltonianData.NOrbitals;
                    hamiltonian.NElectrons = hamiltonianData.NElectrons;
                    HashSet<Int64[]> recordedIJTerms = new HashSet<Int64[]>(new Extensions.IntArrayIEqualityComparer());
                    foreach (var ijTerm in fileHIJTerms)
                    {
                        // If equivalent overlaps are specified, only record the first instance.
                        var orderedTerm = ijTerm.Key[0] > ijTerm.Key[1] ? ijTerm.Key.Reverse().ToArray() : ijTerm.Key;
                        if (recordedIJTerms.Contains(orderedTerm) == false)
                        {
                            recordedIJTerms.Add(orderedTerm);
                            if (Math.Abs(ijTerm.Value) >= threshold)
                            {
                                hamiltonian.AddFermionTerm(new OrbitalIntegral(ijTerm.Key.Select(o => o - 1), ijTerm.Value, OrbitalIntegral.Convention.Mulliken));
                            }
                        }
                    }
                    foreach (var ijklTerm in fileHIJKLTerms)
                    {
                        var tmp = new OrbitalIntegral(ijklTerm.Key.Select(o => o - 1), ijklTerm.Value, OrbitalIntegral.Convention.Mulliken);
                        if (Math.Abs(ijklTerm.Value) >= threshold)
                        {
                            hamiltonian.AddFermionTerm(new OrbitalIntegral(ijklTerm.Key.Select(o => o - 1), ijklTerm.Value, OrbitalIntegral.Convention.Mulliken));
                        }
                    }
                    hamiltonian.SortAndAccumulate();

                    if (!(hamiltonianData.InitialStates == null))
                    {

                        foreach (var state in hamiltonianData.InitialStates)
                        {
                            var label = state.Label;
                            var energy = state.Energy?.Value ?? 0.0;
                            var superpositionRaw = state.Superposition;
                            var stringData = superpositionRaw.Select(o => o.Select(k => k.ToString()).ToList());
                            var superposition = stringData.Select(o => BroombridgeTyped.ParseInputState(o)).ToArray();
                            //Only have terms with non-zero amplitudes.
                            superposition = superposition.Where(o => Math.Abs(o.Item1.Item1) >= threshold).ToArray();
                            if (superposition.Count() > 0)
                            {
                                hamiltonian.InputStates.Add(label,
                                    new FermionHamiltonian.InputState {  Energy = energy, Superposition = superposition }
                                    );
                            }
                        }
                    }
                    return hamiltonian;
                }
            );
        }



        
    }
}
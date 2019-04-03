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
            var broombridgeData = Broombridge.Deserialize.v0_1(filename);

            return LoadData.LoadIntegralData(broombridgeData);
        }

    }


    public partial class LoadData
    {
        internal static IEnumerable<FermionHamiltonian> LoadIntegralData(Broombridge.V0_1.Data schemaInstance, Double threshold = 1e-8)
        {
            return schemaInstance.IntegralSets.Select(
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

                    fileHIJTerms = hamiltonianData.Hamiltonian.OneElectronIntegrals.Values;
                    fileHIJKLTerms = hamiltonianData.Hamiltonian.TwoElectronIntegrals.Values;
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

                    if (!(hamiltonianData.SuggestedState == null))
                    {

                        foreach (var state in hamiltonianData.SuggestedState)
                        {
                            var label = state.SuggestedStateData.Label;
                            var energy = state.SuggestedStateData.Energy?.Value ?? 0.0;
                            var superpositionRaw = state.SuggestedStateData.Superposition;
                            var stringData = superpositionRaw.Select(o => o.Select(k => k.ToString()).ToList());
                            var superposition = stringData.Select(o => ParseInputState(o)).ToArray();
                            //Only have terms with non-zero amplitudes.
                            superposition = superposition.Where(o => Math.Abs(o.Item1.Item1) >= threshold).ToArray();
                            if (superposition.Count() > 0)
                            {
                                hamiltonian.InputStates.Add(
                                    new FermionHamiltonian.InputState { Label = label, Energy = energy, Superposition = superposition }
                                    );
                            }
                        }
                    }
                    return hamiltonian;
                }
            );
        }



        public static ((Double, Double), FermionTerm) ParseInputState(List<string> superpositionElement)
        {
            var amplitude = Double.Parse(superpositionElement.First(), System.Globalization.CultureInfo.InvariantCulture);
            var initialState = superpositionElement.Last();
            var ca = new List<Int64>();
            var so = new List<SpinOrbital>();

            for (int i = 1; i < superpositionElement.Count() - 1; i++)
            {
                FermionTerm singleTerm = ParsePolishNotation(superpositionElement[i]);
                ca.Add(singleTerm.CreationAnnihilationIndices.First());
                so.Add(singleTerm.SpinOrbitalIndices.First());
            }
            FermionTerm term = new FermionTerm(ca.ToArray(), so.ToArray(), amplitude);

            var canonicalOrder = term.ToCanonicalOrder();
            var noAnnihilationTerm = canonicalOrder.Where(o => !(o.CreationAnnihilationIndices.Contains(0)));
            var finalAmplitude = 0.0;
            // check if there are any valid terms.
            if (noAnnihilationTerm.Count() == 1)
            {
                term = noAnnihilationTerm.Single();
                finalAmplitude = term.coeff;
            }
            term.coeff = 1.0;
            return ((finalAmplitude, 0.0), term);
        }

        public static FermionTerm ParsePolishNotation(string input)
        {
            // Regex match examples: (1a)+ (2a)+ (3a)+ (4a)+ (5a)+ (6a)+ (1b)+ (2b)- (3b)+
            Regex regex = new Regex(@"(\((?<orbital>\d+)(?<spin>[ab])\)(?<operator>\+*))");
            Match match = regex.Match(input);
            if (match.Success)
            {
                var orbital = Int64.Parse(match.Groups["orbital"].ToString()) - 1;
                var spin = match.Groups["spin"].ToString() == "a" ? Spin.u : Spin.d;
                var conjugate = match.Groups["operator"].ToString() == "+" ? 1 : 0;
                return new FermionTerm(new long[] { conjugate }, new SpinOrbital[] { new SpinOrbital(orbital, spin) }, 1.0);
            }
            else
            {
                throw new System.ArgumentException($"{input} is not valid Polish notation");
            }
        }
    }
}
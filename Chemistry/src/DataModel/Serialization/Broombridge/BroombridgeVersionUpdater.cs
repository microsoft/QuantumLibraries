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

    internal static partial class DataStructures
    {
        public static V0_3.Data Update(V0_2.Data input) =>
            new V0_3.Data
            {
                Bibliography = input.Bibliography,
                Format = input.Format,
                Generator = input.Generator,
                Schema = V0_3.SchemaUrl,
                ProblemDescriptions = input
                    .ProblemDescriptions
                    .Select(problem => new V0_3.ProblemDescription
                    {
                        BasisSet = problem.BasisSet,
                        CoulombRepulsion = problem.CoulombRepulsion,
                        EnergyOffset = problem.EnergyOffset,
                        FciEnergy = problem.FciEnergy,
                        Geometry = problem.Geometry,
                        InitialStates = problem.InitialStates,
                        Metadata = problem.Metadata,
                        NElectrons = problem.NElectrons,
                        NOrbitals = problem.NOrbitals,
                        ScfEnergy = problem.ScfEnergy,
                        ScfEnergyOffset = problem.ScfEnergyOffset,
                        Hamiltonian = problem.Hamiltonian.ToBroombridgeV0_3()
                    })
                    .ToList()
            };

        /// <summary>
        /// Converts v0.1 Broombridge to v0.2.
        /// </summary>
        /// <param name="input">Source Broombridge in v0.1 format.</param>
        /// <returns>Converted Broombridge in v0.2 format.</returns>
        public static V0_3.Data Update(V0_1.Data input)
        {
            var output = new V0_3.Data()
            {
                Schema = input.Schema,
                Format = input.Format,
                Generator = input.Generator,
                Bibliography = input.Bibliography,
                ProblemDescriptions = new List<V0_3.ProblemDescription>()
            };

            foreach (var integralSet in input.IntegralSets)
            {
                var problemDescription = new V0_3.ProblemDescription()
                {
                    Metadata = integralSet.Metadata,
                    BasisSet = integralSet.BasisSet,
                    Geometry = integralSet.Geometry,
                    CoulombRepulsion = integralSet.CoulombRepulsion,
                    ScfEnergy = integralSet.ScfEnergy,
                    ScfEnergyOffset = integralSet.ScfEnergyOffset,
                    FciEnergy = integralSet.FciEnergy,
                    NOrbitals = integralSet.NOrbitals,
                    NElectrons = integralSet.NElectrons,
                    EnergyOffset = integralSet.EnergyOffset,
                    Hamiltonian = integralSet.Hamiltonian.ToBroombridgeV0_3(),
                    InitialStates = new List<V0_2.State>()
                };
            

                
                if (integralSet.SuggestedState != null)
                {
                    foreach (var sourceInitialState in integralSet.SuggestedState)
                    {
                        var initialState = new V0_2.State()
                        {
                            Label = sourceInitialState.SuggestedStateData.Label,
                            Energy = sourceInitialState.SuggestedStateData.Energy,
                            Method = V0_2.UpdaterStrings.SparseMultiConfigurational,
                            Superposition = sourceInitialState.SuggestedStateData.Superposition
                        };

                        problemDescription.InitialStates.Add(initialState);
                    }
                }

                output.ProblemDescriptions.Add(problemDescription);
            }

            return output;
        }
    }

    
}
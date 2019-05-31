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

    public static partial class DataStructures
    {
        
        /// <summary>
        /// Converts v0.1 Broombridge to v0.2.
        /// </summary>
        /// <param name="input">Source Broombridge in v0.1 format.</param>
        /// <returns>Converted Broombridge in v0.2 format.</returns>
        public static V0_2.Data Update(V0_1.Data input)
        {
            var output = new V0_2.Data()
            {
                Schema = input.Schema,
                Format = input.Format,
                Generator = input.Generator,
                Bibliography = input.Bibliography,
                ProblemDescriptions = new List<V0_2.ProblemDescription>()
            };

            foreach (var integralSet in input.IntegralSets)
            {
                var problemDescription = new V0_2.ProblemDescription()
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
                    Hamiltonian = integralSet.Hamiltonian,
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
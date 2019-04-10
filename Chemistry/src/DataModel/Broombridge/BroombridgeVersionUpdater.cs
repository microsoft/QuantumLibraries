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
            var output = new V0_2.Data();
            output.Schema = input.Schema;
            output.Format = input.Format;
            output.Generator = input.Generator;
            output.Bibliography = input.Bibliography;
            output.ProblemDescriptions = new List<V0_2.ProblemDescription>();

            foreach (var integralSet in input.IntegralSets)
            {
                var problemDescription = new V0_2.ProblemDescription();
                problemDescription.Metadata = integralSet.Metadata;
                problemDescription.BasisSet = integralSet.BasisSet;
                problemDescription.Geometry = integralSet.Geometry;
                problemDescription.CoulombRepulsion = integralSet.CoulombRepulsion;
                problemDescription.ScfEnergy = integralSet.ScfEnergy;
                problemDescription.ScfEnergyOffset = integralSet.ScfEnergyOffset;
                problemDescription.FciEnergy = integralSet.FciEnergy;
                problemDescription.NOrbitals = integralSet.NOrbitals;
                problemDescription.NElectrons = integralSet.NElectrons;
                problemDescription.EnergyOffset = integralSet.EnergyOffset;
                problemDescription.Hamiltonian = integralSet.Hamiltonian;

                problemDescription.InitialStates = new List<V0_2.State>();
                foreach (var sourceInitialState in integralSet.SuggestedState)
                {
                    var initialState = new V0_2.State();
                    initialState.Label = sourceInitialState.SuggestedStateData.Label;
                    initialState.Energy = sourceInitialState.SuggestedStateData.Energy;
                    initialState.Method = V0_2.UpdaterStrings.SparseMultiConfigurational;
                    initialState.Superposition = sourceInitialState.SuggestedStateData.Superposition;

                    problemDescription.InitialStates.Add(initialState);
                }

                output.ProblemDescriptions.Add(problemDescription);
            }

            return output;
        }
    }

    
}
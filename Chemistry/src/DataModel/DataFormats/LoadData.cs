// Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
// Microsoft Software License Terms for Microsoft Quantum Simulation Library (Preview).
// See LICENSE.md in the project root for license information.


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
    
    
}
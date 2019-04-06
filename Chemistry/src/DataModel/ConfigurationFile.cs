// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{

    
    /// <summary>
    /// Configuration settings for modifying chemistry library behavior.
    /// </summary>
    public class Config
    {

        /// <summary>
        /// Construct default configuration.
        /// </summary>
        /// <returns>Default configuration class.</returns>
        public static Config Default()
        {
            return new Config();
        }

        /// <summary>
        /// Default configuration constructor;
        /// </summary>
        public Config()
        {
            indexConvention = IndexConvention.Default;
            truncationThreshold = TruncationThreshold.Default;
        }

        /// <summary>
        /// Choose indexing convention from spin-orbital index to an integer.
        /// </summary>
        public IndexConvention.Type IndexConvention;

        /// <summary>
        /// Available indexing convention from spin-orbital index to an integer.
        /// </summary>
        public class IndexConvention : SpinOrbital.Config.IndexConvention { }

        /// <summary>
        /// Threshold below which to truncate Hamiltonian coefficients.
        /// </summary>
        public double truncationThreshold;

        public static class TruncationThreshold
        {
            public const double Default = 1e-8;
        }
        
    }

  
}


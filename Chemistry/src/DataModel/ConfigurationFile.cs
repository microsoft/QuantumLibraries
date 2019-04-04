// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{
    public class Configuration
    {
        /// <summary>
        /// Default configuration constructor;
        /// </summary>
        public Configuration()
        {
            IndexConvention = IndexConventionClass.Default();
            TruncationThreshold = 1e-8;
        }

        /// <summary>
        /// Choose indexing convention from spin-orbital index to an integer.
        /// </summary>
        public IndexConventionClass.Type IndexConvention = IndexConventionClass.Default();

        /// <summary>
        /// Available indexing convention from spin-orbital index to an integer.
        /// </summary>
        public static class IndexConventionClass {
            public enum Type
            {
                UpDown, HalfUp
            }
            public static Type Default()
            {
                return Type.HalfUp;
            }                
        }

        /// <summary>
        /// Threshold hold below which to truncate Hamiltonian coefficients.
        /// </summary>
        public double TruncationThreshold;
        
    }
}


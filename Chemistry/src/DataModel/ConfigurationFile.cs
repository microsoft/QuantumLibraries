// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.



using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;

namespace Microsoft.Quantum.Chemistry
{
    
    public class Config
    {
        /// <summary>
        /// Default configuration constructor;
        /// </summary>
        public Config()
        {
            indexConvention = IndexConvention.Default();
            truncationThreshold = TruncationThreshold.Default;
        }

        /// <summary>
        /// Choose indexing convention from spin-orbital index to an integer.
        /// </summary>
        public IndexConvention.Type indexConvention;

        /// <summary>
        /// Available indexing convention from spin-orbital index to an integer.
        /// </summary>
        public static class IndexConvention {
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
        public double truncationThreshold;

        public static class TruncationThreshold
        {
            public static double Default = 1e-8;
        }
        
    }

  
}


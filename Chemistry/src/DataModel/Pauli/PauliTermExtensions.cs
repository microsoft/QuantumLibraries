// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;


namespace Microsoft.Quantum.Chemistry.Pauli
{

    /// <summary>
    /// Extensions For Pauli objects.
    /// </summary>
    public static class Extensions
    {
        /// <summary>
        /// Construct PauliTermValue that implements the ITermValue interface. 
        /// </summary>
        /// <param name="x">Input double.</param>
        /// <returns>PauliTermValue representing the input double.</returns>
        public static PauliTermValue ToPauliTermValue(this double x)
        {
            return new PauliTermValue(x);
        }

        /// <summary>
        /// Construct PauliTermValue that implements the ITermValue interface. 
        /// </summary>
        /// <param name="x">Input double  sequence.</param>
        /// <returns>PauliTermValue representing the input double sequence.</returns>
        public static PauliTermValue ToPauliTermValue(this IEnumerable<double> x)
        {
            return new PauliTermValue(x);
        }
    }
}


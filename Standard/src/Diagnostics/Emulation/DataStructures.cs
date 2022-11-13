// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Collections.Generic;
using System.Collections.Immutable;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    /// <summary>
    ///     Represents a complex-vaued square matrix using a one-dimensional list of `double` values.
    /// </summary>
    public class ComplexSquareMatrix
    {
        public int Dimension { get; set; } = 0;
        private double[] entries;

        public ComplexSquareMatrix(int dim)
        {
            Dimension = dim;
            entries = new double[dim * dim * 2];
        }

        /// <summary>
        /// Accesses an entry by providing row and column and part (0 = real, 1 = imag) of complex value.
        /// </summary>
        public double this[int row, int col, int part]
        {
            get
            {
                var index = row * (Dimension * 2) + col * 2 + part;
                return entries[index];
            }
            set
            {
                var index = row * (Dimension * 2) + col * 2 + part;
                entries[index] = value;
            }
        }

        public override string ToString()
        {
            var realRows = new string[Dimension];
            var imagRows = new string[Dimension];
            var realEntries = new double[Dimension];
            var imagEntries = new double[Dimension];

            var index = 0;
            for (var row = 0; row < Dimension; ++row)
            {
                for (var col = 0; col < Dimension; ++col)
                {
                    realEntries[col] = entries[index++];
                    imagEntries[col] = entries[index++];
                }

                realRows[row] = $"[{string.Join(", ", realEntries)}]";
                imagRows[row] = $"[{string.Join(", ", imagEntries)}]";
            }

            return $"Real:\n[{string.Join(", \n", realRows)}]\nImag:\n[{string.Join(", \n", imagRows)}]";
        }
    }

    /// <summary>
    ///     Represents a unitary operator intended for use as a diagnostic
    ///     display.
    /// </summary>
    public class DisplayableUnitaryOperator
    {
        /// <summary>
        ///     The qubits on which the represented operator acts, or
        ///     <c>null</c> if there is no specific register associated with
        ///     this operator.
        /// </summary>
        public IList<Qubit>? Qubits { get; set; }
        
        /// <summary>
        ///     An array of matrix elements for the given unitary operator.
        /// </summary>
        /// <remarks>
        ///     For ease of use, this matrix has stores <c>double</c> values in a
        ///     shape <c>(dim, dim, 2)</c>, where <c>dim</c> is the dimension
        ///     of the represented unitary operator (i.e.: 2^nQubits), and
        ///     where the last axis represents the real (index 0) and
        ///     imaginary parts, respectively.
        /// </remarks>
        public ComplexSquareMatrix? Data { get; set; }

        public override string ToString() => Data == null ? "" : Data!.ToString();
    }

    /// <summary>
    ///     A diagnostic record of sites where a given operation or function
    ///     was called.
    /// </summary>
    public struct CallSites
    {
        /// <summary>
        ///     The name of the operation or function whose calls are
        ///     represented by this record.
        /// </summary>
        public string Subject { get; set; }
        
        /// <summary>
        ///     A collection of calls to the given operation or function, each
        ///     represented as a call stack at the point where the subject was
        ///     called.
        /// </summary>
        public ImmutableList<ImmutableStack<string>> Sites { get; set; }
    }
}

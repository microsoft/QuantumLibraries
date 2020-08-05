// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Collections.Generic;
using System.Collections.Immutable;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Newtonsoft.Json;

namespace Microsoft.Quantum.Diagnostics.Emulation
{    
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
        ///     For ease of use, this array has a dtype of <c>double</c>, and
        ///     shape <c>(dim, dim, 2)</c>, where <c>dim</c> is the dimension
        ///     of the represented unitary operator (i.e.: 2^nQubits), and
        ///     where the last axis represents the real (index 0) and
        ///     imaginary parts, respectively.
        /// </remarks>
        public NumSharp.NDArray? Data { get; set; }
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

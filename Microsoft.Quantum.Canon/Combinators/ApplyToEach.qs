// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// <seealso cref="ApplyToEach" />
    operation ApplyToEachAC(singleQubitOperation : (Qubit => ():Adjoint,Controlled), register : Qubit[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// <seealso cref="ApplyToEach" />
    operation ApplyToEachA(singleQubitOperation : (Qubit => ():Adjoint), register : Qubit[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        adjoint auto
    }

    /// <seealso cref="ApplyToEach" />
    operation ApplyToEachC(singleQubitOperation : (Qubit => ():Controlled), register : Qubit[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }

        controlled auto
    }

    /// <summary>
    ///     Applies a single-qubit operation to each qubit in a register.
    /// </summary>
    /// <param name="singleQubitOperation">Operation to apply to each qubit.</param>
    /// <param name="register">Array of qubits on which to apply the given operation.</param>
    /// <example>
    ///     Prepare a three-qubit |+〉 state:
    ///     <c>
    ///         using (register = Qubit[3]) {
    ///             ApplyToEach(H, register)
    ///         }
    ///     </c>
    /// </example>
    /// <seealso cref="ApplyToEachA" />
    /// <seealso cref="ApplyToEachC" />
    /// <seealso cref="ApplyToEachAC" />
    operation ApplyToEach(singleQubitOperation : (Qubit => ()), register : Qubit[])  : ()
    {
        body {
            for (idxQubit in 0..(Length(register) - 1)) {
                singleQubitOperation(register[idxQubit]);
            }
        }
    }

}

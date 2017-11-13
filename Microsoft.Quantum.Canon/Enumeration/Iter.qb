// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    // Needed for range accessors.
    open Microsoft.Quantum.Primitive

    /// <seealso cref="Iter" />
    //operation IterCA<'T>(op : ('T => () : Controlled, Adjoint), array : 'T[])  : ()
    {
        Body {
            Iter(op, array)
        }

        Adjoint {
            (Adjoint IterA)(op, array)
        }
        Controlled (controls) {
            (Controlled IterC)(controls, (op, array))
        }
        Controlled Adjoint (controls) {
            (Controlled IterC)(controls, (Adjoint op, array))
        }
    }

    /// <seealso cref="Iter" />
    operation IterA<'T>(op : ('T => () : Adjoint), array : 'T[])  : ()
    {
        Body {
            Iter(op, array)
        }

        Adjoint {
            Iter(Adjoint op, array)
        }
    }

    /// <seealso cref="Iter" />
    operation IterC<'T>(op : ('T => () : Controlled), array : 'T[])  : ()
    {
        Body {
            Iter(op, array)
        }

        Controlled (controls) {
            Iter((Controlled op)(controls, _), array)
        }
    }

    /// <summary>
    ///     Applies an operation to each element of an array.
    /// </summary>
    /// <param name="singleQubitOperation">Operation to apply to each qubit.</param>
    /// <param name="register">Array of qubits on which to apply the given operation.</param>
    /// <example>
    ///     Prepare a three-qubit |+〉 state:
    ///     <c>
    ///         using (register = Qubit[3]) {
    ///             Iter(H, register)
    ///         }
    ///     </c>
    /// </example>
    /// <seealso cref="IterA" />
    /// <seealso cref="IterC" />
    /// <seealso cref="IterCA" />
    operation Iter<'T>(op : ('T => ()), array : 'T[])  : ()
    {
        Body {
            let nElements = Length(array)
            for (idxElement in 0..(nElements - 1)) {
                op(array[idxElement])
            }
        }
    }

    operation IterIndex<'T>(op : ((Int, 'T) => ()), array : 'T[]) : () {
        Body {
            let nElements = Length(array)
            for (idxElement in 0..(nElements - 1)) {
                op(idxElement, array[idxElement])
            }
        }
    } 

    // TODO: A, C, CA variants.

    operation IterRange(op : (Int => ()), range : Range) : () {
        Body {
            for (idx in range) {
                op(idx)
            }
        }
    }

    operation IterRangeCA(op : (Int => () : Adjoint, Controlled), range : Range) : () {
        Body {
            IterRange(op, range)
        }
        Adjoint {
            // FIXME: add a ReverseRange fn.
            //IterRange(Adjoint op, Stop(range)..(-Step(range))..Start(range))
            // For now using just range, which is wrong, but allows it to parse
            IterRange(Adjoint op, range)
        }
        Controlled (controlRegister) {
            IterRange((Controlled op)(controlRegister, _), range)
        }
        Controlled Adjoint (controlRegister) {
        // For now using just range, which is wrong, but allows it to parse
           // IterRange((Controlled Adjoint op)(controlRegister, _), Stop(range)..(-Step(range))..Start(range))
           IterRange((Controlled Adjoint op)(controlRegister, _), range)
        }
    }


    // TODO: A, C variants

}

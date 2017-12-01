// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    operation ApplyToSubregister(op : (Qubit[] => ()), idxs : Int[], target : Qubit[]) : () {
        body {
            let subregister = Subarray(idxs, target);
            op(subregister);
        }
    }

    operation ApplyToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
    }

    operation ApplyToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }

    operation ApplyToSubregisterCA(op : (Qubit[] => () : Controlled, Adjoint), idxs : Int[], target : Qubit[]) : () {
        body {
            ApplyToSubregister(op, idxs, target);
        }
        adjoint {
            ApplyToSubregister(Adjoint op, idxs, target);
        }
        controlled (controls) {
            let cop = (Controlled op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
        controlled adjoint (controls) {
            let cop = (Controlled Adjoint op);
            ApplyToSubregister(cop(controls, _), idxs, target);
        }
    }

    function RestrictToSubregister(op : (Qubit[] => ()), idxs : Int[]) : (Qubit[] => ()) {
        return ApplyToSubregister(op, idxs, _);
    }
    function RestrictToSubregisterA(op : (Qubit[] => () : Adjoint), idxs : Int[]) : (Qubit[] => () : Adjoint) {
        return ApplyToSubregisterA(op, idxs, _);
    }
    function RestrictToSubregisterC(op : (Qubit[] => () : Controlled), idxs : Int[]) : (Qubit[] => () : Controlled) {
        return ApplyToSubregisterC(op, idxs, _);
    }
    function RestrictToSubregisterCA(op : (Qubit[] => () : Adjoint, Controlled), idxs : Int[]) : (Qubit[] => () : Adjoint, Controlled) {
        return ApplyToSubregisterCA(op, idxs, _);
    }

}

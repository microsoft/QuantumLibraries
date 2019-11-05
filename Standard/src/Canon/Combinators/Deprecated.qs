// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWith".
    @Deprecated("Microsoft.Quantum.Canon.ApplyWith")
    operation With<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit), target : 'T) : Unit {
        ApplyWith(outerOperation, innerOperation, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithA".
    @Deprecated("Microsoft.Quantum.Canon.ApplyWithA")
    operation WithA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj), target : 'T) : Unit {
        body (...) {
            ApplyWithA(outerOperation, innerOperation, target);
        }
        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithC".
    @Deprecated("Microsoft.Quantum.Canon.ApplyWithC")
    operation WithC<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Ctl), target : 'T) : Unit {
        body (...) {
            ApplyWithC(outerOperation, innerOperation, target);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithCA".
    @Deprecated("Microsoft.Quantum.Canon.ApplyWithCA")
    operation WithCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl), target : 'T) : Unit {
        body (...) {
            ApplyWithCA(outerOperation, innerOperation, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregister".
    @Deprecated("Microsoft.Quantum.Canon.RestrictedToSubregister")
    function RestrictToSubregister(op : (Qubit[] => Unit), idxs : Int[]) : (Qubit[] => Unit) {
        return RestrictedToSubregister(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregistera".
    @Deprecated("Microsoft.Quantum.Canon.RestrictedToSubregisterA")
    function RestrictToSubregisterA(op : (Qubit[] => Unit is Adj), idxs : Int[]) : (Qubit[] => Unit) {
        return RestrictedToSubregisterA(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterc".
    @Deprecated("Microsoft.Quantum.Canon.RestrictedToSubregisterC")
    function RestrictToSubregisterC(op : (Qubit[] => Unit is Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        return RestrictedToSubregisterC(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterca".
    @Deprecated("Microsoft.Quantum.Canon.RestrictedToSubregisterCA")
    function RestrictToSubregisterCA(op : (Qubit[] => Unit is Adj + Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        return RestrictedToSubregisterCA(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.curriedop".
    @Deprecated("Microsoft.Quantum.Canon.CurriedOp")
    function CurryOp<'T, 'U> (op : (('T, 'U) => Unit)) : ('T -> ('U => Unit)) {
        return CurriedOp(op);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedop".
    @Deprecated("Microsoft.Quantum.Canon.UncurriedOp")
    function UncurryOp<'T, 'U> (curriedOp : ('T -> ('U => Unit))) : (('T, 'U) => Unit) {
        return UncurriedOp(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopa".
    @Deprecated("Microsoft.Quantum.Canon.UncurriedOpA")
    function UncurryOpA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj))) : (('T, 'U) => Unit is Adj) {
        return UncurriedOpA(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopc".
    @Deprecated("Microsoft.Quantum.Canon.UncurriedOpC")
    function UncurryOpC<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl))) : (('T, 'U) => Unit is Ctl) {
        return UncurriedOpC(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopca".
    @Deprecated("Microsoft.Quantum.Canon.UncurriedOpCA")
    function UncurryOpCA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj + Ctl))) : (('T, 'U) => Unit is Adj + Ctl) {
        return UncurriedOpCA(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.bound".
    @Deprecated("Microsoft.Quantum.Canon.Bound")
    function Bind<'T> (operations : ('T => Unit)[]) : ('T => Unit) {
        return Bound(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.bounda".
    @Deprecated("Microsoft.Quantum.Canon.BoundA")
    function BindA<'T> (operations : ('T => Unit is Adj)[]) : ('T => Unit is Adj) {
        return BoundA(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.boundc".
    @Deprecated("Microsoft.Quantum.Canon.BoundC")
    function BindC<'T> (operations : ('T => Unit is Ctl)[]) : ('T => Unit is Ctl) {
        return BoundC(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.boundca".
    @Deprecated("Microsoft.Quantum.Canon.BoundCA")
    function BindCA<'T> (operations : ('T => Unit is Adj + Ctl)[]) : ('T => Unit is Adj + Ctl) {
        return BoundCA(operations);
    }

}

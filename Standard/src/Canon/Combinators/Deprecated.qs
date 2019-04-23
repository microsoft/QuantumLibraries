// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Warnings;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWith".
    operation With<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit), target : 'T) : Unit {
        _Renamed("Microsoft.Quantum.Canon.With", "Microsoft.Quantum.Canon.ApplyWith");
        ApplyWith(outerOperation, innerOperation, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithA".
    operation WithA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj), target : 'T) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.WithA", "Microsoft.Quantum.Canon.ApplyWithA");
            ApplyWithA(outerOperation, innerOperation, target);
        }
        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithC".
    operation WithC<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Ctl), target : 'T) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.WithC", "Microsoft.Quantum.Canon.ApplyWithC");
            ApplyWithC(outerOperation, innerOperation, target);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithCA".
    operation WithCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl), target : 'T) : Unit {
        body (...) {
            _Renamed("Microsoft.Quantum.Canon.WithCA", "Microsoft.Quantum.Canon.ApplyWithCA");
            ApplyWithCA(outerOperation, innerOperation, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregister".
    function RestrictToSubregister(op : (Qubit[] => Unit), idxs : Int[]) : (Qubit[] => Unit) {
        _Renamed("Microsoft.Quantum.Canon.RestrictToSubregister", "Microsoft.Quantum.Canon.RestrictedToSubregister");
        return RestrictedToSubregister(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregistera".
    function RestrictToSubregisterA(op : (Qubit[] => Unit is Adj), idxs : Int[]) : (Qubit[] => Unit) {
        _Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterA", "Microsoft.Quantum.Canon.RestrictedToSubregisterA");
        return RestrictedToSubregisterA(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterc".
    function RestrictToSubregisterC(op : (Qubit[] => Unit is Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        _Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterC", "Microsoft.Quantum.Canon.RestrictedToSubregisterC");
        return RestrictedToSubregisterC(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterca".
    function RestrictToSubregisterCA(op : (Qubit[] => Unit is Adj + Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        _Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterCA", "Microsoft.Quantum.Canon.RestrictedToSubregisterCA");
        return RestrictedToSubregisterCA(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.curriedop".
    function CurryOp<'T, 'U> (op : (('T, 'U) => Unit)) : ('T -> ('U => Unit)) {
        _Renamed("Microsoft.Quantum.Canon.CurryOp", "Microsoft.Quantum.Canon.CurriedOp");
        return CurriedOp(op);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedop".
    function UncurryOp<'T, 'U> (curriedOp : ('T -> ('U => Unit))) : (('T, 'U) => Unit) {
        _Renamed("Microsoft.Quantum.Canon.UncurryOp", "Microsoft.Quantum.Canon.UncurriedOp");
        return UncurriedOp(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopa".
    function UncurryOpA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj))) : (('T, 'U) => Unit is Adj) {
        _Renamed("Microsoft.Quantum.Canon.UncurryOpA", "Microsoft.Quantum.Canon.UncurriedOpA");
        return UncurriedOpA(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopc".
    function UncurryOpC<'T, 'U> (curriedOp : ('T -> ('U => Unit is Ctl))) : (('T, 'U) => Unit is Ctl) {
        _Renamed("Microsoft.Quantum.Canon.UncurryOpC", "Microsoft.Quantum.Canon.UncurriedOpC");
        return UncurriedOpC(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.uncurriedopca".
    function UncurryOpCA<'T, 'U> (curriedOp : ('T -> ('U => Unit is Adj + Ctl))) : (('T, 'U) => Unit is Adj + Ctl) {
        _Renamed("Microsoft.Quantum.Canon.UncurryOpCA", "Microsoft.Quantum.Canon.UncurriedOpCA");
        return UncurriedOpCA(curriedOp);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.bound".
    function Bind<'T> (operations : ('T => Unit)[]) : ('T => Unit) {
        _Renamed("Microsoft.Quantum.Canon.Bind", "Microsoft.Quantum.Canon.Bound");
        return Bound(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.bounda".
    function BindA<'T> (operations : ('T => Unit is Adj)[]) : ('T => Unit is Adj) {
        _Renamed("Microsoft.Quantum.Canon.BindA", "Microsoft.Quantum.Canon.BoundA");
        return BoundA(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.boundc".
    function BindC<'T> (operations : ('T => Unit is Ctl)[]) : ('T => Unit is Ctl) {
        _Renamed("Microsoft.Quantum.Canon.BindC", "Microsoft.Quantum.Canon.BoundC");
        return BoundC(operations);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.boundca".
    function BindCA<'T> (operations : ('T => Unit is Adj + Ctl)[]) : ('T => Unit is Adj + Ctl) {
        _Renamed("Microsoft.Quantum.Canon.BindCA", "Microsoft.Quantum.Canon.BoundCA");
        return BoundCA(operations);
    }

}

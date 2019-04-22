// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWith".
    operation With<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit), target : 'T) : Unit {
        Renamed("Microsoft.Quantum.Canon.With", "Microsoft.Quantum.Canon.ApplyWith");
        ApplyWith(outerOperation, innerOperation, target);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithA".
    operation WithA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj), target : 'T) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Canon.WithA", "Microsoft.Quantum.Canon.ApplyWithA");
            ApplyWithA(outerOperation, innerOperation, target);
        }
        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithC".
    operation WithC<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Ctl), target : 'T) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Canon.WithC", "Microsoft.Quantum.Canon.ApplyWithC");
            ApplyWithC(outerOperation, innerOperation, target);
        }
        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Canon.ApplyWithCA".
    operation WithCA<'T>(outerOperation : ('T => Unit is Adj), innerOperation : ('T => Unit is Adj + Ctl), target : 'T) : Unit {
        body (...) {
            Renamed("Microsoft.Quantum.Canon.WithCA", "Microsoft.Quantum.Canon.ApplyWithCA");
            ApplyWithCA(outerOperation, innerOperation, target);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregister".
    function RestrictToSubregister(op : (Qubit[] => Unit), idxs : Int[]) : (Qubit[] => Unit) {
        Renamed("Microsoft.Quantum.Canon.RestrictToSubregister", "Microsoft.Quantum.Canon.RestrictedToSubregister");
        return RestrictedToSubregister(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregistera".
    function RestrictToSubregisterA(op : (Qubit[] => Unit is Adj), idxs : Int[]) : (Qubit[] => Unit) {
        Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterA", "Microsoft.Quantum.Canon.RestrictedToSubregisterA");
        return RestrictedToSubregisterA(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterc".
    function RestrictToSubregisterC(op : (Qubit[] => Unit is Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterC", "Microsoft.Quantum.Canon.RestrictedToSubregisterC");
        return RestrictedToSubregisterC(op, idxs);
    }

    /// # Deprecated
    /// Please use @"microsoft.quantum.canon.restrictedtosubregisterca".
    function RestrictToSubregisterCA(op : (Qubit[] => Unit is Adj + Ctl), idxs : Int[]) : (Qubit[] => Unit) {
        Renamed("Microsoft.Quantum.Canon.RestrictToSubregisterCA", "Microsoft.Quantum.Canon.RestrictedToSubregisterCA");
        return RestrictedToSubregisterCA(op, idxs);
    }

}

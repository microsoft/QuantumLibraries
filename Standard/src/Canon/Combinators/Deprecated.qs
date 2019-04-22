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

}

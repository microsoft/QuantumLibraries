// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;

    @Test("ToffoliSimulator")
    operation TestSinglyControlled() : Unit {
        use ctl = Qubit();
        use target = Qubit[3];
        let targetLE = LittleEndian(target);
        let value = 5;

        for singlyControlled in [SinglyControlled, SinglyControlledA] {
            for enableControl in [false, true] {
                within {
                    ApplyIfA(enableControl, X, ctl);
                } apply {
                    singlyControlled(ApplyXorInPlace(value, _))(ctl, targetLE);
                    EqualityFactI(MeasureInteger(targetLE), enableControl ? value | 0, "Unexpected measurement result for SinglyControlled");
                }
            }
        }
    }
}

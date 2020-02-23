// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Diagnostics;

    operation ManuallyReflectAboutFive(register : Qubit[]) : Unit is Adj + Ctl {
        within {
            X(register[1]);
        } apply {
            Controlled Z(register[0..1], register[2]);
        }
    }

    operation ReflectAboutFiveUsingLibrary(register : Qubit[]) : Unit is Adj + Ctl {
        let littleEndian = LittleEndian(register);
        ReflectAboutInteger(5, littleEndian);
    }

    operation ReflectAboutIntegerTest() : Unit {
        AssertOperationsEqualReferenced(3,
            ReflectAboutFiveUsingLibrary,
            ManuallyReflectAboutFive
        );
    }

}

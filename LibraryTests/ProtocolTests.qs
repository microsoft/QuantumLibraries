// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    // given that the Teleportation circuit is correct this operation must be an identity
    operation TeleportationIdentityTestHelper (arg : Qubit[]) : () {
        body {
            AssertIntEqual(Length(arg), 1, "Helper is defined only of single qubit input");
            using (anc = Qubit[1])
            {
                Teleportation(arg[0], anc[0]);
                SWAP(arg[0], anc[0]);
            }
        }
    }

    operation TeleportationTest () : () {
        body {
            // given that there is randomness involved in the Teleportation, repeat the tests several times.
            for(idxIteration in 1 .. 8)
            {
                AssertOperationsEqualInPlace(TeleportationIdentityTestHelper, IdentityTestHelper, 1);
                AssertOperationsEqualReferenced(TeleportationIdentityTestHelper, IdentityTestHelper, 1);
            }
        }
    }

}

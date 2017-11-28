// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    // NB: the following are stubs to allow this file to at least parse, even without
    //     defining the correct action of the operations.

    operation LEStub( target : LittleEndian)  : ()
    {
        body { }
    }

    operation BEStub( target : BigEndian)  : ()
    {
        body { }
    }

    // Design notes:
    //
    //     We want a single operation that either performs a cyclic shift
    //     or drops a qubit, depending on the value of an argument.
    //     The most natural means of doing so might be to have a Boolean
    //     argument, but that's not currently in the language (see Solid #406).
    //
    //     As a next best thing, we could use global constants to kluge together
    //     an enum-like concept, but that's also not in the language at the moment.
    //
    //     In lieu, we take an Int argument instead, and document the values we
    //     accept.

    /// # Summary
    /// Given an operation acting on a little-endian register
    ///     |q0 q1 . q?>,
    /// returns a new operation acting on the big end-shifted register
    ///     |q? q0 q1 . q??1>,
    /// with qubit q? being dropped if the cyclic argument is 0.
    ///
    /// # Input
    /// ## cyclic
    /// Specifies if the shift is cyclic (1) or drops
    /// shifted qubits (0).
    /// ## op
    /// Operation whose action is to be shifted.
    operation BigShiftOpLE(cyclic : Int, op : (LittleEndian => ()))  : (LittleEndian => ())
    {
        body {
            return LEStub;
        }
    }

    /// # Summary
    /// Given an operation acting on a little-endian register
    ///     |q0 q1 . q?>,
    /// returns a new operation acting on the little end-shifted register
    ///     |q1 . q?>,
    /// with qubit q0 being dropped being dropped if the cyclic argument is 0.
    operation LittleShiftOpLE(cyclic : Int, op : (LittleEndian => ()))  : (LittleEndian => ())
    {
        body {
            return LEStub;
        }
    }

    /// # Summary
    /// Given an operation acting on a little-endian register
    ///     |q? q??1 . q0>,
    /// returns a new operation acting on the big end-shifted register
    ///     |q??1 . q0>,
    /// with qubit q? being dropped being dropped if the cyclic argument is 0.
    operation BigShiftOpBE(cyclic : Int, op : (BigEndian => ()))  : (BigEndian => ())
    {
        body {
            return BEStub;
        }
    }

    /// # Summary
    /// Given an operation acting on a big-endian register
    ///     |q? q??1 . q0>,
    /// returns a new operation acting on the little end-shifted register
    ///     |q? q??1 . q1>,
    /// with qubit q0 being dropped being dropped if the cyclic argument is 0.
    operation LittleShiftOpBE(cyclic : Int, op : (BigEndian => ()))  : (BigEndian => ())
    {
        body {
            return BEStub;
        }
    }

    // Design notes:
    //     We want to expose a signature of the form
    //         ((LittleEndian) => ()) => ((LittleEndian) => ())
    //     for each operation. With partial application, this is easiest to
    //     do by having an "impl" operation of the form
    //         (((LittleEndian) => ()), LittleEndian)
    //     that we partially apply on the second argument.
    //
    //     There are only two shifts (left and right),
    //     but they need to be exposed using both BE and LE types,
    //     each of which has different naming conventions.
    //     Depending on the variance rules in Solid #511, we may need to split
    //     the two "impl" operations into four to remove UDT labels.
    //     For example:

    operation RightShiftOpImpl(cyclic : Int, op : (LittleEndian => ()),  target : LittleEndian)  : ()
    {
        body {
            // TODO
        }
    }

}

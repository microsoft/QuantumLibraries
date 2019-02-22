// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Applies an operation that takes little-endian input to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a little-endian register.
    /// ## register
    /// A big-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpLEA
    /// - ApplyReversedOpLEC
    /// - ApplyReversedOpLECA
    operation ApplyReversedOpLE(op : (LittleEndian => Unit), register : BigEndian) : Unit {
        ApplyWithInputTransformation(BigEndianAsLittleEndian, op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpLE".
    operation ApplyReversedOpLittleEndian(op : (LittleEndian => Unit), register : BigEndian) : Unit {
        Renamed("ApplyReversedOpLittleEndian", "ApplyReversedOpLE");
        ApplyReversedOpLE(op, register);
    }

    /// # Summary
    /// Given an operation that takes a little-endian input, returns a new
    /// operation that takes a big-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a big-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpLE
    /// - ReversedOpLEA
    /// - ReversedOpLEC
    /// - ReversedOpLECA
    function ReversedOpLE(op : (LittleEndian => Unit)) : (BigEndian => Unit) {
        return ApplyReversedOpLE(op, _);
    }

    /// # Summary
    /// Applies an operation that takes little-endian input to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a little-endian register.
    /// ## register
    /// A big-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpLE
    /// - ApplyReversedOpLEC
    /// - ApplyReversedOpLECA
    operation ApplyReversedOpLEA (op : (LittleEndian => Unit : Adjoint), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationA(BigEndianAsLittleEndian, op, register);
        }

        adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpLEA".
    operation ApplyReversedOpLittleEndianA(op : (LittleEndian => Unit : Adjoint), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianA", "ApplyReversedOpLEA");
            ApplyReversedOpLEA(op, register);
        }
        adjoint auto;
    }

    /// # Summary
    /// Given an operation that takes a little-endian input, returns a new
    /// operation that takes a big-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a big-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpLEA
    /// - ReversedOpLE
    /// - ReversedOpLEC
    /// - ReversedOpLECA
    function ReversedOpLEA(op : (LittleEndian => Unit : Adjoint)) : (BigEndian => Unit : Adjoint) {
        return ApplyReversedOpLEA(op, _);
    }

    /// # Summary
    /// Applies an operation that takes little-endian input to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a little-endian register.
    /// ## register
    /// A big-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpLE
    /// - ApplyReversedOpLEA
    /// - ApplyReversedOpLECA
    operation ApplyReversedOpLEC(op : (LittleEndian => Unit : Controlled), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationC(BigEndianAsLittleEndian, op, register);
        }

        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianC(op : (LittleEndian => Unit : Controlled), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianC", "ApplyReversedOpLEC");
            ApplyReversedOpLEC(op, register);
        }
        controlled auto;
    }

    /// # Summary
    /// Given an operation that takes a little-endian input, returns a new
    /// operation that takes a big-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a big-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpLEC
    /// - ReversedOpLE
    /// - ReversedOpLEA
    /// - ReversedOpLECA
    function ReversedOpLEC(op : (LittleEndian => Unit : Controlled)) : (BigEndian => Unit : Controlled) {
        return ApplyReversedOpLEC(op, _);
    }

    /// # Summary
    /// Applies an operation that takes little-endian input to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a little-endian register.
    /// ## register
    /// A big-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpLE
    /// - ApplyReversedOpLEA
    /// - ApplyReversedOpLEC
    operation ApplyReversedOpLECA(op : (LittleEndian => Unit : Controlled, Adjoint), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationCA(BigEndianAsLittleEndian, op, register);
        }

        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpLEC".
    operation ApplyReversedOpLittleEndianCA(op : (LittleEndian => Unit : Adjoint, Controlled), register : BigEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpLittleEndianCA", "ApplyReversedOpLECA");
            ApplyReversedOpLECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Summary
    /// Given an operation that takes a little-endian input, returns a new
    /// operation that takes a big-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a big-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpLECA
    /// - ReversedOpLE
    /// - ReversedOpLEA
    /// - ReversedOpLEC
    function ReversedOpLECA(op : (LittleEndian => Unit : Adjoint, Controlled)) : (BigEndian => Unit : Adjoint, Controlled) {
        return ApplyReversedOpLECA(op, _);
    }

    /// # Summary
    /// Applies an operation that takes big-endian input to a register encoding
    /// an unsigned integer using little-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a big-endian register.
    /// ## register
    /// A little-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpBEA
    /// - ApplyReversedOpBEC
    /// - ApplyReversedOpBECA
    operation ApplyReversedOpBE(op : (BigEndian => Unit), register : LittleEndian) : Unit {
        ApplyWithInputTransformation(LittleEndianAsBigEndian, op, register);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpBE".
    operation ApplyReversedOpBigEndian(op : (BigEndian => Unit), register : LittleEndian) : Unit {
        Renamed("ApplyReversedOpBigEndian", "ApplyReversedOpBE");
        ApplyReversedOpBE(op, register);
    }

    /// # Summary
    /// Given an operation that takes a big-endian input, returns a new
    /// operation that takes a little-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a little-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpBE
    /// - ReversedOpBEA
    /// - ReversedOpBEC
    /// - ReversedOpBECA
    function ReversedOpBE(op : (BigEndian => Unit)) : (LittleEndian => Unit) {
        return ApplyReversedOpBE(op, _);
    }

    /// # Summary
    /// Applies an operation that takes big-endian input to a register encoding
    /// an unsigned integer using little-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a big-endian register.
    /// ## register
    /// A little-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpBE
    /// - ApplyReversedOpBEC
    /// - ApplyReversedOpBECA
    operation ApplyReversedOpBEA (op : (BigEndian => Unit : Adjoint), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationA(LittleEndianAsBigEndian, op, register);
        }

        adjoint invert;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpBEA".
    operation ApplyReversedOpBigEndianA(op : (BigEndian => Unit : Adjoint), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianA", "ApplyReversedOpBEA");
            ApplyReversedOpBEA(op, register);
        }
        adjoint auto;
    }

    /// # Summary
    /// Given an operation that takes a big-endian input, returns a new
    /// operation that takes a little-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a little-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpBEA
    /// - ReversedOpBE
    /// - ReversedOpBEC
    /// - ReversedOpBECA
    function ReversedOpBEA(op : (BigEndian => Unit : Adjoint)) : (LittleEndian => Unit : Adjoint) {
        return ApplyReversedOpBEA(op, _);
    }

    /// # Summary
    /// Applies an operation that takes big-endian input to a register encoding
    /// an unsigned integer using little-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a big-endian register.
    /// ## register
    /// A little-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpBE
    /// - ApplyReversedOpBEA
    /// - ApplyReversedOpBECA
    operation ApplyReversedOpBEC(op : (BigEndian => Unit : Controlled), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationC(LittleEndianAsBigEndian, op, register);
        }

        controlled auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianC(op : (BigEndian => Unit : Controlled), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianC", "ApplyReversedOpBEC");
            ApplyReversedOpBEC(op, register);
        }
        controlled auto;
    }

    /// # Summary
    /// Given an operation that takes a big-endian input, returns a new
    /// operation that takes a little-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a little-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpBEC
    /// - ReversedOpBE
    /// - ReversedOpBEA
    /// - ReversedOpBECA
    function ReversedOpBEC(op : (BigEndian => Unit : Controlled)) : (LittleEndian => Unit : Controlled) {
        return ApplyReversedOpBEC(op, _);
    }

    /// # Summary
    /// Applies an operation that takes big-endian input to a register encoding
    /// an unsigned integer using little-endian format.
    ///
    /// # Input
    /// ## op
    /// Operation that acts on a big-endian register.
    /// ## register
    /// A little-endian register to be transformed.
    ///
    /// # See Also
    /// - ApplyReversedOpBE
    /// - ApplyReversedOpBEA
    /// - ApplyReversedOpBEC
    operation ApplyReversedOpBECA(op : (BigEndian => Unit : Controlled, Adjoint), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationCA(LittleEndianAsBigEndian, op, register);
        }

        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.ApplyReversedOpBEC".
    operation ApplyReversedOpBigEndianCA(op : (BigEndian => Unit : Adjoint, Controlled), register : LittleEndian) : Unit {
        body (...) {
            Renamed("ApplyReversedOpBigEndianCA", "ApplyReversedOpBECA");
            ApplyReversedOpBECA(op, register);
        }
        adjoint auto;
        controlled auto;
        controlled adjoint auto;
    }

    /// # Summary
    /// Given an operation that takes a big-endian input, returns a new
    /// operation that takes a little-endian input.
    ///
    /// # Input
    /// ## op
    /// The operation whose input is to be reversed.
    ///
    /// # Output
    /// A new operation that accepts its input as a little-endian register.
    ///
    /// # See Also
    /// - ApplyReversedOpBECA
    /// - ReversedOpBE
    /// - ReversedOpBEA
    /// - ReversedOpBEC
    function ReversedOpBECA(op : (BigEndian => Unit : Adjoint, Controlled)) : (LittleEndian => Unit : Adjoint, Controlled) {
        return ApplyReversedOpBECA(op, _);
    }

}

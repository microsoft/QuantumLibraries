// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLECA
    operation ApplyReversedOpLE(op : (LittleEndian => Unit), register : BigEndian) : Unit {
        ApplyWithInputTransformation(BigEndianAsLittleEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLECA
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLECA
    operation ApplyReversedOpLEA (op : (LittleEndian => Unit is Adj), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationA(BigEndianAsLittleEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLECA
    function ReversedOpLEA(op : (LittleEndian => Unit is Adj)) : (BigEndian => Unit is Adj) {
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLECA
    operation ApplyReversedOpLEC(op : (LittleEndian => Unit is Ctl), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationC(BigEndianAsLittleEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLECA
    function ReversedOpLEC(op : (LittleEndian => Unit is Ctl)) : (BigEndian => Unit is Ctl) {
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLEC
    operation ApplyReversedOpLECA(op : (LittleEndian => Unit is Ctl + Adj), register : BigEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationCA(BigEndianAsLittleEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpLECA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpLEC
    function ReversedOpLECA(op : (LittleEndian => Unit is Adj + Ctl)) : (BigEndian => Unit is Adj + Ctl) {
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBECA
    operation ApplyReversedOpBE(op : (BigEndian => Unit), register : LittleEndian) : Unit {
        ApplyWithInputTransformation(LittleEndianAsBigEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBECA
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBECA
    operation ApplyReversedOpBEA (op : (BigEndian => Unit is Adj), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationA(LittleEndianAsBigEndian, op, register);
        }

        adjoint invert;
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBECA
    function ReversedOpBEA(op : (BigEndian => Unit is Adj)) : (LittleEndian => Unit is Adj) {
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBECA
    operation ApplyReversedOpBEC(op : (BigEndian => Unit is Ctl), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationC(LittleEndianAsBigEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBECA
    function ReversedOpBEC(op : (BigEndian => Unit is Ctl)) : (LittleEndian => Unit is Ctl) {
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBEC
    operation ApplyReversedOpBECA(op : (BigEndian => Unit is Ctl + Adj), register : LittleEndian) : Unit {
        body (...) {
            ApplyWithInputTransformationCA(LittleEndianAsBigEndian, op, register);
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
    /// - Microsoft.Quantum.Arithmetic.ApplyReversedOpBECA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBE
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEA
    /// - Microsoft.Quantum.Arithmetic.ReversedOpBEC
    function ReversedOpBECA(op : (BigEndian => Unit is Adj + Ctl)) : (LittleEndian => Unit is Adj + Ctl) {
        return ApplyReversedOpBECA(op, _);
    }

}

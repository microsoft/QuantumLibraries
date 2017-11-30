// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Reorders the qubits in a register to obtain a new
    /// qubit register.
    /// 
    /// # Input 
    /// ## register
    /// The register from which a new register with the reversed order of qubits being built
    function ReverseRegister(register : Qubit[])  : Qubit[]
    {
        mutable reversed = new Qubit[Length(register)];

        for (idxQubit in 0..Length(register) - 1) {
            set reversed[idxQubit] = register[Length(register) - 1 - idxQubit ];
        }

        return reversed;
    }

    /// # Summary
    /// Uses SWAP gates to reverse the order of the qubits in
    /// a register.
    /// 
    /// # Input 
    /// ## register
    /// The qubits order of which should be reversed using SWAP gates
    operation SwapReverseRegister(register : Qubit[])  : ()
    {
        body {
            let totalQubits = Length(register);
            let halfTotal = totalQubits / 2;
            for( i in 0 .. halfTotal - 1 ) {
                SWAP(register[i],register[ totalQubits - i - 1 ]);
            }
        }
        adjoint self
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Register that encodes unsigned integer in little-endian order
    /// qubit with index 0 encodes the lowest bit of an unsigned integer
    newtype LittleEndian = Qubit[];

    /// # Summary
    /// Register that encodes unsigned integer in big-endian order
    /// qubit with index 0 encodes the highest bit of an unsigned integer
    newtype BigEndian = Qubit[];

    /// # Summary 
    /// Applies an operation that takes little-endian input to a register encoding 
    /// an unsigned integer using big-endian format.
    /// 
    /// # Input 
    /// ## op 
    /// Operation that acts on little-endian register
    /// ## register
    /// big-endian register to be transformed
    /// 
    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedoplittleendiana"
    /// - @"microsoft.quantum.canon.applyreversedoplittleendianc"
    /// - @"microsoft.quantum.canon.applyreversedoplittleendianca"
    operation ApplyReversedOpLittleEndian( 
              op : (LittleEndian => ()),
              register : BigEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianA( 
              op : (LittleEndian => () : Adjoint),  
              register : BigEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        adjoint auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianC(
              op : (LittleEndian => () : Controlled),
              register : BigEndian )  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        controlled auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianCA(
              op : (LittleEndian => () : Controlled, Adjoint),  
              register : BigEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary 
    /// Applies an operation that takes big-endian input to a register encoding 
    /// an unsigned integer using little-endian format.
    ///
    /// # Input 
    /// ## op 
    /// Operation that acts on big-endian register
    /// ## register
    /// little-endian register to be transformed
    /// 
    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedopbigendiana"
    /// - @"microsoft.quantum.canon.applyreversedopbigendianc"
    /// - @"microsoft.quantum.canon.applyreversedopbigendianca"
    operation ApplyReversedOpBigEndian(
              op : (BigEndian => ()),
              register : LittleEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianA(
              op : (BigEndian => () : Adjoint),
              register : LittleEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        adjoint auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianC(
              op : (BigEndian => () : Controlled),
              register : LittleEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        controlled auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianCA(
              op : (BigEndian => () : Controlled, Adjoint),
              register : LittleEndian)  : ()
    {
        body {
            let bareReversed = ReverseRegister(register);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }
}

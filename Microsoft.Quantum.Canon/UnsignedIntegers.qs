// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Canon {

    open Microsoft.Quantum.Primitive;

    /// # Summary
    /// Register that encodes unsigned integer in little-endian order
    /// qubit with index 0 encodes the lowest bit of an unsigned integer
    newtype LittleEndian = Qubit[];

    /// # Summary
    /// Register that encodes unsigned integer in big-endian order
    /// qubit with index 0 encodes the highest bit of an unsigned integer
    newtype BigEndian = Qubit[];

    /// # Summary 
    /// Little-endian unsigned integers in QFT basis. 
    /// For example, if |x⟩ is little-endian encoding of integer x in computational basis, 
    /// then QFTLE|x⟩ is encoding of x in QFT basis. 
    newtype PhaseLittleEndian = (Qubit[]);

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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
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
            let bareReversed = Reverse(AsQubitArray(register));
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
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
    /// Applies an operation that takes PhaseLittleEndian input on LittleEndian target
    ///
    /// # Input
    /// ## op
    /// Operation to be applied
    /// ## target
    /// The register to which the operation is applied
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applyphaseleoperationonlea"
    /// - @"microsoft.quantum.canon.applyphaseleoperationonlec"
    /// - @"microsoft.quantum.canon.applyphaseleoperationonleca"
    operation ApplyPhaseLEOperationOnLE( op : (PhaseLittleEndian => ()), target : LittleEndian ) : () {
        body {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            op(phaseLE);
            (Adjoint QFTLE)(target);
        }
    }
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLEA( op : (PhaseLittleEndian => () : Adjoint), target : LittleEndian ) : () {
        body {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            op(phaseLE);
            (Adjoint QFTLE)(target);
        }
        adjoint auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLEC( op : (PhaseLittleEndian => () : Controlled), target : LittleEndian ) : () {
        body {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            op(phaseLE);
            (Adjoint QFTLE)(target);
        }
        controlled( controls ) {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            (Controlled op)(controls, phaseLE);
            (Adjoint QFTLE)(target);
        }
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLECA( op : (PhaseLittleEndian => () : Controlled, Adjoint), target : LittleEndian ) : () {
        body {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            op(phaseLE);
            (Adjoint QFTLE)(target);
        }
        adjoint auto
        controlled( controls ) {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target);
            (Controlled op)(controls, phaseLE);
            (Adjoint QFTLE)(target);
        }
        controlled adjoint auto
    }

    /// # Summary
    /// Applies an operation that takes PhaseLittleEndian input on LittleEndian target
    ///
    /// # Input
    /// ## op
    /// Operation to be applied
    /// ## target
    /// The register to which the operation is applied
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applyleoperationonphaselea"
    /// - @"microsoft.quantum.canon.applyleoperationonphaselec"
    /// - @"microsoft.quantum.canon.applyleoperationonphaseleca"
    operation ApplyLEOperationOnPhaseLE( op : (LittleEndian => ()), target : PhaseLittleEndian ) : () {
        body {
            let targetLE = LittleEndian(target);
            With(Adjoint(QFTLE),op,targetLE);
        }
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyleoperationonphasele"
    operation ApplyLEOperationOnPhaseLEA( op : (LittleEndian => () : Adjoint), target : PhaseLittleEndian ) : () {
        body {
            let targetLE = LittleEndian(target);
            WithA(Adjoint(QFTLE),op,targetLE);
        }
        adjoint auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyleoperationonphasele"
    operation ApplyLEOperationOnPhaseLEC( op : (LittleEndian => () : Controlled), target : PhaseLittleEndian ) : () {
        body {
            let targetLE = LittleEndian(target);
            WithC(Adjoint(QFTLE),op,targetLE);
        }
        controlled auto
    }

    /// # See Also 
    /// - @"microsoft.quantum.canon.applyleoperationonphasele"
    operation ApplyLEOperationOnPhaseLECA( op : (LittleEndian => () : Controlled, Adjoint ), target : PhaseLittleEndian ) : () {
        body {
            let targetLE = LittleEndian(target);
            WithCA(Adjoint(QFTLE),op,targetLE);
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }
}

// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Canon
{
    
    open Microsoft.Quantum.Primitive;
    
    
    /// # Summary
    /// Register that encodes an unsigned integer in little-endian order. The
    /// qubit with index `0` encodes the lowest bit of an unsigned integer
    ///
    /// # Remarks
    /// We abbreviate `LittleEndian` as `LE` in the documentation.
    newtype LittleEndian = Qubit[];
    
    /// # Summary
    /// Register that encodes an unsigned integer in big-endian order. The
    /// qubit with index `0` encodes the highest bit of an unsigned integer
    ///
    /// # Remarks
    /// We abbreviate `BigEndian` as `BE` in the documentation.
    newtype BigEndian = Qubit[];
    
    /// # Summary
    /// Little-endian unsigned integers in QFT basis.
	///
    /// For example, if `|x⟩` is the little-endian encoding of the integer `x` in the computational basis,
    /// then `QFTLE|x⟩` is the encoding of `x` in the QFT basis.
    ///
    /// # Remarks
    /// We abbreviate `LittleEndian` as `LE` in the documentation.
    ///
    /// # See Also
    /// - microsoft.quantum.canon.qft
    /// - microsoft.quantum.canon.qftle
    newtype PhaseLittleEndian = Qubit[];
    
    
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
    /// - Microsoft.Quantum.Canon.ApplyReversedOpLittlEndianA
    /// - Microsoft.Quantum.Canon.ApplyReversedOpLittlEndianC
    /// - Microsoft.Quantum.Canon.ApplyReversedOpLittlEndianCA
    operation ApplyReversedOpLittleEndian (op : (LittleEndian => Unit), register : BigEndian) : Unit
    {
        let bareReversed = Reverse(register!);
        let reversed = LittleEndian(bareReversed);
        op(reversed);
    }
    
    
    /// # Summary
    /// Applies an operation that takes little-endian input and that supports
    /// the adjoint functor to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianA (op : (LittleEndian => Unit : Adjoint), register : BigEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        
        adjoint invert;
    }
    
    
    /// # Summary
    /// Applies an operation that takes little-endian input and that supports
    /// the controlled functor to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianC (op : (LittleEndian => Unit : Controlled), register : BigEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        
        controlled distribute;
    }
    
    
    /// # Summary
    /// Applies an operation that takes little-endian input and that supports
    /// the controlled and adjoint functors to a register encoding
    /// an unsigned integer using big-endian format.
    ///
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedoplittleendian"
    operation ApplyReversedOpLittleEndianCA (op : (LittleEndian => Unit : Controlled, Adjoint), register : BigEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = LittleEndian(bareReversed);
            op(reversed);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
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
    operation ApplyReversedOpBigEndian (op : (BigEndian => Unit), register : LittleEndian) : Unit
    {
        let bareReversed = Reverse(register!);
        let reversed = BigEndian(bareReversed);
        op(reversed);
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianA (op : (BigEndian => Unit : Adjoint), register : LittleEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        
        adjoint invert;
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianC (op : (BigEndian => Unit : Controlled), register : LittleEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        
        controlled distribute;
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyreversedopbigendian"
    operation ApplyReversedOpBigEndianCA (op : (BigEndian => Unit : Controlled, Adjoint), register : LittleEndian) : Unit
    {
        body (...)
        {
            let bareReversed = Reverse(register!);
            let reversed = BigEndian(bareReversed);
            op(reversed);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    /// # Summary
    /// Uses SWAP gates to reverse the order of the qubits in
    /// a register.
    ///
    /// # Input
    /// ## register
    /// The qubits order of which should be reversed using SWAP gates
    operation SwapReverseRegister (register : Qubit[]) : Unit
    {
        body (...)
        {
            let totalQubits = Length(register);
            let halfTotal = totalQubits / 2;
            
            for (i in 0 .. halfTotal - 1)
            {
                SWAP(register[i], register[(totalQubits - i) - 1]);
            }
        }
        
        adjoint self;
        controlled distribute;
        controlled adjoint self;
    }
    
    
    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.canon.littleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.canon.phaselittleendian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `PhaseLittleEndian` by the use of
    /// <xref:microsoft.quantum.canon.qftle> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationonLEA
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationonLEA
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationonLECA
    operation ApplyPhaseLEOperationOnLE (op : (PhaseLittleEndian => Unit), target : LittleEndian) : Unit
    {
        QFTLE(target);
        let phaseLE = PhaseLittleEndian(target!);
        op(phaseLE);
        Adjoint QFTLE(target);
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLEA (op : (PhaseLittleEndian => Unit : Adjoint), target : LittleEndian) : Unit
    {
        body (...)
        {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
            Adjoint QFTLE(target);
        }
        
        adjoint invert;
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLEC (op : (PhaseLittleEndian => Unit : Controlled), target : LittleEndian) : Unit
    {
        body (...)
        {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
            Adjoint QFTLE(target);
        }
        
        controlled (controls, ...)
        {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target!);
            Controlled op(controls, phaseLE);
            Adjoint QFTLE(target);
        }
    }
    
    
    /// # See Also
    /// - @"microsoft.quantum.canon.applyphaseleoperationonle"
    operation ApplyPhaseLEOperationOnLECA (op : (PhaseLittleEndian => Unit : Controlled, Adjoint), target : LittleEndian) : Unit
    {
        body (...)
        {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
            Adjoint QFTLE(target);
        }
        
        adjoint invert;
        
        controlled (controls, ...)
        {
            QFTLE(target);
            let phaseLE = PhaseLittleEndian(target!);
            Controlled op(controls, phaseLE);
            Adjoint QFTLE(target);
        }
        
        controlled adjoint invert;
    }
    
    
    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.canon.phaselittleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.canon.littleendian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `LittleEndian` by the use of
    /// <xref:microsoft.quantum.canon.qftle> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEA
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEA
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLECA
    operation ApplyLEOperationOnPhaseLE (op : (LittleEndian => Unit), target : PhaseLittleEndian) : Unit
    {
        let targetLE = LittleEndian(target!);
        With(Adjoint QFTLE, op, targetLE);
    }
    
    
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    operation ApplyLEOperationOnPhaseLEA (op : (LittleEndian => Unit : Adjoint), target : PhaseLittleEndian) : Unit
    {
        body (...)
        {
            let targetLE = LittleEndian(target!);
            WithA(Adjoint QFTLE, op, targetLE);
        }
        
        adjoint invert;
    }
    
    
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    operation ApplyLEOperationOnPhaseLEC (op : (LittleEndian => Unit : Controlled), target : PhaseLittleEndian) : Unit
    {
        body (...)
        {
            let targetLE = LittleEndian(target!);
            WithC(Adjoint QFTLE, op, targetLE);
        }
        
        controlled distribute;
    }
    
    
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    operation ApplyLEOperationOnPhaseLECA (op : (LittleEndian => Unit : Controlled, Adjoint), target : PhaseLittleEndian) : Unit
    {
        body (...)
        {
            let targetLE = LittleEndian(target!);
            WithCA(Adjoint QFTLE, op, targetLE);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
}



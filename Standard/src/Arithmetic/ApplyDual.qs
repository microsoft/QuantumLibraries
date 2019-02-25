// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;

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

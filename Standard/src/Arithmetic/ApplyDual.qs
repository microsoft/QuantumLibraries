// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

///////////////////////////////////////////////////////////////////////////////////////////
// Types and supporting functions for representing unsigned integers in arrays of qubits //
///////////////////////////////////////////////////////////////////////////////////////////

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;

    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.arithmetic.littleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.arithmetic.phaselittleendian>.
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
    operation ApplyPhaseLEOperationOnLE (op : (PhaseLittleEndian => Unit), target : LittleEndian) : Unit {
        QFTLE(target);
        let phaseLE = PhaseLittleEndian(target!);
        op(phaseLE);
        Adjoint QFTLE(target);
    }

    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationOnLE
    operation ApplyPhaseLEOperationOnLEA (op : (PhaseLittleEndian => Unit is Adj), target : LittleEndian) : Unit
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
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationOnLE
    operation ApplyPhaseLEOperationOnLEC (op : (PhaseLittleEndian => Unit is Ctl), target : LittleEndian) : Unit
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
    /// - Microsoft.Quantum.Canon.ApplyPhaseLEOperationOnLE
    operation ApplyPhaseLEOperationOnLECA (op : (PhaseLittleEndian => Unit is Adj + Ctl), target : LittleEndian) : Unit
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
    /// <xref:microsoft.quantum.arithmetic.phaselittleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.arithmetic.littleendian>.
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
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEC
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLECA
    operation ApplyLEOperationOnPhaseLE (op : (LittleEndian => Unit), target : PhaseLittleEndian) : Unit
    {
        let targetLE = LittleEndian(target!);
        ApplyWith(Adjoint QFTLE, op, targetLE);
    }



    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.arithmetic.phaselittleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.arithmetic.littleendian>.
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
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEC
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLECA
    operation ApplyLEOperationOnPhaseLEA (op : (LittleEndian => Unit is Adj), target : PhaseLittleEndian)
    : Unit is Adj {
        let targetLE = LittleEndian(target!);
        ApplyWithA(Adjoint QFTLE, op, targetLE);
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.arithmetic.phaselittleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.arithmetic.littleendian>.
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
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEA
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLECA
    operation ApplyLEOperationOnPhaseLEC (op : (LittleEndian => Unit is Ctl), target : PhaseLittleEndian)
    : Unit is Ctl {
        let targetLE = LittleEndian(target!);
        ApplyWithC(Adjoint QFTLE, op, targetLE);
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:microsoft.quantum.arithmetic.phaselittleendian> register as input
    /// on a target register of type <xref:microsoft.quantum.arithmetic.littleendian>.
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
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLE
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEA
    /// - Microsoft.Quantum.Canon.ApplyLEOperationonPhaseLEC
    operation ApplyLEOperationOnPhaseLECA(op : (LittleEndian => Unit is Adj + Ctl), target : PhaseLittleEndian)
    : Unit is Adj + Ctl {
        let targetLE = LittleEndian(target!);
        ApplyWithCA(Adjoint QFTLE, op, targetLE);
    }

}

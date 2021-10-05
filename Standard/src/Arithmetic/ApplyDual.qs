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
    /// <xref:Microsoft.Quantum.Arithmetic.LittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `PhaseLittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEC
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLECA
    operation ApplyPhaseLEOperationOnLE (op : (PhaseLittleEndian => Unit), target : LittleEndian) : Unit {
        within {
            QFTLE(target);
        } apply {
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
        }
    }

    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.LittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `PhaseLittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLE
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEC
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLECA
    operation ApplyPhaseLEOperationOnLEA (op : (PhaseLittleEndian => Unit is Adj), target : LittleEndian)
    : Unit is Adj {
        within {
            QFTLE(target);
        } apply {
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
        }
    }


      /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.LittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `PhaseLittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLE
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLECA
    operation ApplyPhaseLEOperationOnLEC (op : (PhaseLittleEndian => Unit is Ctl), target : LittleEndian)
    : Unit is Ctl {
        within {
            QFTLE(target);
        } apply {
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
        }
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.LittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `PhaseLittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLE
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyPhaseLEOperationOnLEC
    operation ApplyPhaseLEOperationOnLECA (op : (PhaseLittleEndian => Unit is Adj + Ctl), target : LittleEndian)
    : Unit is Adj + Ctl {
        within {
            QFTLE(target);
        } apply {
            let phaseLE = PhaseLittleEndian(target!);
            op(phaseLE);
        }
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.LittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `LittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyLEOperationOnPhaseLEA
    /// - Microsoft.Quantum.Canon.ApplyLEOperationOnPhaseLEC
    /// - Microsoft.Quantum.Canon.ApplyLEOperationOnPhaseLECA
    operation ApplyLEOperationOnPhaseLE (op : (LittleEndian => Unit), target : PhaseLittleEndian)
    : Unit {
        let targetLE = LittleEndian(target!);
        ApplyWith(Adjoint QFTLE, op, targetLE);
    }



    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.LittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `LittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLE
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLEC
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLECA
    operation ApplyLEOperationOnPhaseLEA (op : (LittleEndian => Unit is Adj), target : PhaseLittleEndian)
    : Unit is Adj {
        let targetLE = LittleEndian(target!);
        ApplyWithA(Adjoint QFTLE, op, targetLE);
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.LittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `LittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLE
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLECA
    operation ApplyLEOperationOnPhaseLEC (op : (LittleEndian => Unit is Ctl), target : PhaseLittleEndian)
    : Unit is Ctl {
        let targetLE = LittleEndian(target!);
        ApplyWithC(Adjoint QFTLE, op, targetLE);
    }


    /// # Summary
    /// Applies an operation that takes a
    /// <xref:Microsoft.Quantum.Arithmetic.PhaseLittleEndian> register as input
    /// on a target register of type <xref:Microsoft.Quantum.Arithmetic.LittleEndian>.
    ///
    /// # Input
    /// ## op
    /// The operation to be applied.
    /// ## target
    /// The register to which the operation is applied.
    ///
    /// # Remarks
    /// The register is transformed to `LittleEndian` by the use of
    /// <xref:Microsoft.Quantum.Canon.QFTLE> and is then returned to
    /// its original representation after application of `op`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLE
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLEA
    /// - Microsoft.Quantum.Arithmetic.ApplyLEOperationOnPhaseLEC
    operation ApplyLEOperationOnPhaseLECA(op : (LittleEndian => Unit is Adj + Ctl), target : PhaseLittleEndian)
    : Unit is Adj + Ctl {
        let targetLE = LittleEndian(target!);
        ApplyWithCA(Adjoint QFTLE, op, targetLE);
    }

}

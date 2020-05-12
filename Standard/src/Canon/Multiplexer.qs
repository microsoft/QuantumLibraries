// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Applies a multiply-controlled unitary operation $U$ that applies a
    /// unitary $V_j$ when controlled by n-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{N-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaryGenerator
    /// A tuple where the first element `Int` is the number of unitaries $N$,
    /// and the second element `(Int -> ('T => () is Adj + Ctl))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$.
    ///
    /// ## index
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Generic qubit register that $V_j$ acts on.
    ///
    /// # Remarks
    /// `coefficients` will be padded with identity elements if
    /// fewer than $2^n$ are specified. This implementation uses
    /// $n-1$ auxiliary qubits.
    ///
    /// # References
    /// - [ *Andrew M. Childs, Dmitri Maslov, Yunseong Nam, Neil J. Ross, Yuan Su*,
    ///      arXiv:1711.10980](https://arxiv.org/abs/1711.10980)
    operation MultiplexOperationsFromGenerator<'T>(unitaryGenerator : (Int, (Int -> ('T => Unit is Adj + Ctl))), index: LittleEndian, target: 'T) : Unit is Ctl + Adj {
        let (nUnitaries, unitaryFunction) = unitaryGenerator;
        let unitaryGeneratorWithOffset = (nUnitaries, 0, unitaryFunction);
        if (Length(index!) == 0) {
            fail "MultiplexOperations failed. Number of index qubits must be greater than 0.";
        }
        if (nUnitaries > 0) {
            let auxiliary = new Qubit[0];
            Adjoint _MultiplexOperationsFromGenerator(unitaryGeneratorWithOffset, auxiliary, index, target);
        }
    }

    /// # Summary
    /// Implementation step of `MultiplexOperationsFromGenerator`.
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexOperationsFromGenerator
    operation _MultiplexOperationsFromGenerator<'T>(unitaryGenerator : (Int, Int, (Int -> ('T => Unit is Adj + Ctl))), auxiliary: Qubit[], index: LittleEndian, target: 'T)
    : Unit {
        body (...) {
            let nIndex = Length(index!);
            let nStates = 2^nIndex;

            let (nUnitaries, unitaryOffset, unitaryFunction) = unitaryGenerator;

            let nUnitariesLeft = MinI(nUnitaries, nStates/2);
            let nUnitariesRight = MinI(nUnitaries, nStates);

            let leftUnitaries = (nUnitariesLeft, unitaryOffset, unitaryFunction);
            let rightUnitaries = (nUnitariesRight-nUnitariesLeft, unitaryOffset + nUnitariesLeft, unitaryFunction);

            let newControls = LittleEndian(index![0..nIndex - 2]);

            if (nUnitaries > 0) {
                if (Length(auxiliary) == 1 and nIndex == 0) {
                    // Termination case

                    (Controlled Adjoint (unitaryFunction(unitaryOffset)))(auxiliary, target);
                } elif (Length(auxiliary) == 0 and nIndex >= 1) {
                    // Start case
                    let newauxiliary = [index![Length(index!) - 1]];
                    if(nUnitariesRight > 0){
                        _MultiplexOperationsFromGenerator(rightUnitaries, newauxiliary, newControls, target);
                    }
                    within {
                        X(newauxiliary[0]);
                    } apply {
                        _MultiplexOperationsFromGenerator(leftUnitaries, newauxiliary, newControls, target);
                    }
                } else {
                    // Recursion that reduces nIndex by 1 & sets Length(auxiliary) to 1.
                    using (newauxiliary = Qubit[1]) {
                        let op = LogicalANDMeasAndFix(_, _);
                        // Naive measurement-free approach uses 4x more T gates with
                        // let op = (Controlled X);
                        op([index![Length(index!) - 1]] + auxiliary, newauxiliary[0]);
                        if (nUnitariesRight > 0) {
                            _MultiplexOperationsFromGenerator(rightUnitaries, newauxiliary, newControls, target);
                        }
                        within {
                            (Controlled X)(auxiliary, newauxiliary[0]);
                        } apply {
                            _MultiplexOperationsFromGenerator(leftUnitaries, newauxiliary, newControls, target);
                        }
                        (Adjoint op)([index![Length(index!) - 1]] + auxiliary, newauxiliary[0]);
                    }
                }
            }
        }
        adjoint auto;
        controlled (controlRegister, (...)) {
            _MultiplexOperationsFromGenerator(unitaryGenerator, auxiliary + controlRegister, index, target);
        }
        adjoint controlled auto;
    }

    /// # Summary
    /// Applies multiply-controlled unitary operation $U$ that applies a
    /// unitary $V_j$ when controlled by n-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{N-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaryGenerator
    /// A tuple where the first element `Int` is the number of unitaries $N$,
    /// and the second element `(Int -> ('T => () is Adj + Ctl))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$.
    ///
    /// ## index
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// little-endian format.
    ///
    /// ## target
    /// Generic qubit register that $V_j$ acts on.
    ///
    /// # Remarks
    /// `coefficients` will be padded with identity elements if
    /// fewer than $2^n$ are specified. This version is implemented
    /// directly by looping through n-controlled unitary operators.
    operation MultiplexOperationsBruteForceFromGenerator<'T>(unitaryGenerator : (Int, (Int -> ('T => Unit is Adj + Ctl))), index: LittleEndian, target: 'T)
    : Unit is Adj + Ctl {
        let nIndex = Length(index!);
        let nStates = 2^nIndex;
        let (nUnitaries, unitaryFunction) = unitaryGenerator;
        for (idxOp in 0..MinI(nStates,nUnitaries) - 1){
            (ControlledOnInt(idxOp, unitaryFunction(idxOp)))(index!, target);
        }
    }

    /// # Summary
    /// Returns a multiply-controlled unitary operation $U$ that applies a
    /// unitary $V_j$ when controlled by n-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{2^n-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaryGenerator
    /// A tuple where the first element `Int` is the number of unitaries $N$,
    /// and the second element `(Int -> ('T => () is Adj + Ctl))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$.
    ///
    /// # Output
    /// A multiply-controlled unitary operation $U$ that applies unitaries
    /// described by `unitaryGenerator`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexOperationsFromGenerator
    function MultiplexerFromGenerator(unitaryGenerator : (Int, (Int -> (Qubit[] => Unit is Adj + Ctl)))) : ((LittleEndian, Qubit[]) => Unit is Adj + Ctl) {
        return MultiplexOperationsFromGenerator(unitaryGenerator, _, _);
    }

    /// # Summary
    /// Returns a multiply-controlled unitary operation $U$ that applies a
    /// unitary $V_j$ when controlled by n-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{2^n-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaryGenerator
    /// A tuple where the first element `Int` is the number of unitaries $N$,
    /// and the second element `(Int -> ('T => () is Adj + Ctl))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$.
    ///
    /// # Output
    /// A multiply-controlled unitary operation $U$ that applies unitaries
    /// described by `unitaryGenerator`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexerBruteForceFromGenerator
    function MultiplexerBruteForceFromGenerator(unitaryGenerator : (Int, (Int -> (Qubit[] => Unit is Adj + Ctl)))) : ((LittleEndian, Qubit[]) => Unit is Adj + Ctl) {
        return MultiplexOperationsBruteForceFromGenerator(unitaryGenerator, _, _);
    }

    /// # Summary
    /// Computes the logical AND of multiple qubits.
    /// # Input
    /// ## ctrlRegister
    /// Qubits input to the multiple-input AND gate.
    /// ## target
    /// Qubit on which output of AND is written to. This is
    /// assumed to always start in the $\ket{0}$ state.
    /// # Remarks
    /// When `ctrlRegister` has exactly $2$ qubits,
    /// this is equivalent to the `CCNOT` operation but phase with a phase
    /// $e^{i\Pi/2}$ on $\ket{001}$ and $-e^{i\Pi/2}$ on $\ket{101}$ and $\ket{011}$.
    /// The Adjoint is also measure-and-fixup approach that is the inverse
    /// of the original operation only in special cases (see references),
    /// but uses fewer T-gates.
    ///
    /// # References
    /// - [ *Craig Gidney*, 1709.06648](https://arxiv.org/abs/1709.06648)
    internal operation LogicalANDMeasAndFix(ctrlRegister : Qubit[], target : Qubit)
    : Unit {
        body (...) {
            if(Length(ctrlRegister) == 2){
                let c1 = ctrlRegister[0];
                let c2 = ctrlRegister[1];
                H(target);
                T(target);
                CNOT(c1,target);
                CNOT(c2,target);
                CNOT(target,c1);
                CNOT(target,c2);
                (Adjoint T)(c1);
                (Adjoint T)(c2);
                T(target);
                CNOT(target,c2);
                CNOT(target,c1);
                H(target);
                S(target);
            } else {
                (Controlled X)(ctrlRegister, target);
            }
        }
        adjoint (...)  {
            if(Length(ctrlRegister) == 2){
                let c1 = ctrlRegister[0];
                let c2 = ctrlRegister[1];
                H(target);
                let Meas = M(target);
                if (Meas == One) {
                    within {
                        H(c2);
                    } apply {
                        CNOT(c1,c2);
                    }
                    X(target);
                }
            } else {
                 (Controlled X)(ctrlRegister, target);
            }
        }
    }
}

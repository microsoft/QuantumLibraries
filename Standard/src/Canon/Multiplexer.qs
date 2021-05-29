// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;
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
        if Length(index!) == 0 {
            fail "MultiplexOperations failed. Number of index qubits must be greater than 0.";
        }
        if nUnitaries > 0 {
            let auxiliary = [];
            Adjoint MultiplexOperationsFromGeneratorImpl(unitaryGeneratorWithOffset, auxiliary, index, target);
        }
    }

    /// # Summary
    /// Implementation step of `MultiplexOperationsFromGenerator`.
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexOperationsFromGenerator
    internal operation MultiplexOperationsFromGeneratorImpl<'T>(unitaryGenerator : (Int, Int, (Int -> ('T => Unit is Adj + Ctl))), auxiliary: Qubit[], index: LittleEndian, target: 'T)
    : Unit {
        body (...) {
            let nIndex = Length(index!);
            let nStates = 2^nIndex;

            let (nUnitaries, unitaryOffset, unitaryFunction) = unitaryGenerator;

            let nUnitariesLeft = MinI(nUnitaries, nStates / 2);
            let nUnitariesRight = MinI(nUnitaries, nStates);

            let leftUnitaries = (nUnitariesLeft, unitaryOffset, unitaryFunction);
            let rightUnitaries = (nUnitariesRight - nUnitariesLeft, unitaryOffset + nUnitariesLeft, unitaryFunction);

            let newControls = LittleEndian(Most(index!));

            if nUnitaries > 0 {
                if Length(auxiliary) == 1 and nIndex == 0 {
                    // Termination case

                    (Controlled Adjoint (unitaryFunction(unitaryOffset)))(auxiliary, target);
                } elif Length(auxiliary) == 0 and nIndex >= 1 {
                    // Start case
                    let newauxiliary = Tail(index!);
                    if nUnitariesRight > 0 {
                        MultiplexOperationsFromGeneratorImpl(rightUnitaries, [newauxiliary], newControls, target);
                    }
                    within {
                        X(newauxiliary);
                    } apply {
                        MultiplexOperationsFromGeneratorImpl(leftUnitaries, [newauxiliary], newControls, target);
                    }
                } else {
                    // Recursion that reduces nIndex by 1 and sets Length(auxiliary) to 1.
                    let controls = [Tail(index!)] + auxiliary;
                    use newauxiliary = Qubit();
                    use andauxiliary = Qubit[MaxI(0, Length(controls) - 2)];
                    within {
                        ApplyAndChain(andauxiliary, controls, newauxiliary);
                    } apply {
                        if nUnitariesRight > 0 {
                            MultiplexOperationsFromGeneratorImpl(rightUnitaries, [newauxiliary], newControls, target);
                        }
                        within {
                            (Controlled X)(auxiliary, newauxiliary);
                        } apply {
                            MultiplexOperationsFromGeneratorImpl(leftUnitaries, [newauxiliary], newControls, target);
                        }
                    }
                }
            }
        }
        adjoint auto;
        controlled (controlRegister, (...)) {
            MultiplexOperationsFromGeneratorImpl(unitaryGenerator, auxiliary + controlRegister, index, target);
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
        for idxOp in 0..MinI(nStates,nUnitaries) - 1 {
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
    /// Computes a chain of AND gates
    ///
    /// # Description
    /// The auxiliary qubits to compute temporary results must be specified explicitly.
    /// The length of that register is `Length(ctrlRegister) - 2`, if there are at least
    /// two controls, otherwise the length is 0.
    internal operation ApplyAndChain(auxRegister : Qubit[], ctrlRegister : Qubit[], target : Qubit)
    : Unit is Adj {
        if Length(ctrlRegister) == 0 {
            X(target);
        } elif Length(ctrlRegister) == 1 {
            CNOT(Head(ctrlRegister), target);
        } else {
            EqualityFactI(Length(auxRegister), Length(ctrlRegister) - 2, "Unexpected number of auxiliary qubits");
            let controls1 = ctrlRegister[0..0] + auxRegister;
            let controls2 = Rest(ctrlRegister);
            let targets = auxRegister + [target];
            ApplyToEachA(ApplyAnd, Zipped3(controls1, controls2, targets));
        }
    }
}

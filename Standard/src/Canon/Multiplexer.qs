// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Applies a multiply-controlled unitary operation $U$ that applies a 
    /// unitary $V_j$ when controlled by n-qubit number state $\ket{j}$.
    ///
    /// $U = \sum^{N-1}_{j=0}\ket{j}\bra{j}\otimes V_j$.
    ///
    /// # Input
    /// ## unitaryGenerator
    /// A tuple where the first element `Int` is the number of unitaries $N$,
    /// and the second element `(Int -> ('T => () : Adjoint, Controlled))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$. 
    ///
    /// ## index
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// big-endian format.
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
    operation MultiplexOperationsFromGenerator<'T>(unitaryGenerator : (Int, (Int -> ('T => Unit : Adjoint, Controlled))), index: BigEndian, target: 'T) : Unit {
        body (...) {
            let (nUnitaries, unitaryFunction) = unitaryGenerator;
            let unitaryGeneratorWithOffset = (nUnitaries, 0, unitaryFunction);
            if (Length(index!) == 0) {
                fail "MultiplexOperations failed. Number of index qubits must be greater than 0.";
            }
            if (nUnitaries > 0) {
                let auxiliary = new Qubit[0];
                Adjoint MultiplexOperationsFromGenerator_(unitaryGeneratorWithOffset, auxiliary, index, target);
            }
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
    }

    /// # Summary
    /// Implementation step of `MultiplexOperationsFromGenerator`.
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexOperationsFromGenerator
    operation MultiplexOperationsFromGenerator_<'T>(unitaryGenerator : (Int, Int, (Int -> ('T => Unit : Adjoint, Controlled))), auxiliary: Qubit[], index: BigEndian, target: 'T) : Unit {
        body (...) {
            let nIndex = Length(index!);
            let nStates = 2^nIndex;
            
            let (nUnitaries, unitaryOffset, unitaryFunction) = unitaryGenerator;

            let nUnitariesLeft = Microsoft.Quantum.Extensions.Math.MinI(nUnitaries, nStates/2);
            let nUnitariesRight = Microsoft.Quantum.Extensions.Math.MinI(nUnitaries, nStates);

            let leftUnitaries = (nUnitariesLeft, unitaryOffset, unitaryFunction);
            let rightUnitaries = (nUnitariesRight-nUnitariesLeft, unitaryOffset + nUnitariesLeft, unitaryFunction);
            
            let newControls = BigEndian(index![1..nIndex-1]);

            if(nUnitaries > 0){
                if(Length(auxiliary) == 1 and nIndex==0){
                    // Termination case
                    
                    (Controlled Adjoint (unitaryFunction(unitaryOffset)))(auxiliary, target);
                }
                elif(Length(auxiliary) == 0 and nIndex>=1){
                    // Start case
                    let newauxiliary = [index![0]];
                    if(nUnitariesRight > 0){
                        MultiplexOperationsFromGenerator_(rightUnitaries, newauxiliary, newControls, target);
                    }
                    X(newauxiliary[0]);
                    MultiplexOperationsFromGenerator_(leftUnitaries, newauxiliary, newControls, target);
                    X(newauxiliary[0]);
                }
                else{
                    // Recursion that reduces nIndex by 1 & sets Length(auxiliary) to 1.
                    using(newauxiliary = Qubit[1]){
                        let op = LogicalANDMeasAndFix_(_, _);
                        // Naive measurement-free approach uses 4x more T gates with 
                        // let op = (Controlled X);
                        op(auxiliary + [index![0]], newauxiliary[0]);
                        if(nUnitariesRight > 0){
                            MultiplexOperationsFromGenerator_(rightUnitaries, newauxiliary, newControls, target);
                        }
                        (Controlled X)(auxiliary, newauxiliary[0]);
                        MultiplexOperationsFromGenerator_(leftUnitaries, newauxiliary, newControls, target);
                        (Controlled X)(auxiliary, newauxiliary[0]);
                        (Adjoint op)(auxiliary + [index![0]], newauxiliary[0]);
                    }
                }
            }
        }
        adjoint auto;
        controlled (controlRegister, (...)) {
            MultiplexOperationsFromGenerator_(unitaryGenerator, auxiliary + controlRegister, index, target);
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
    /// and the second element `(Int -> ('T => () : Adjoint, Controlled))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$. 
    ///
    /// ## index
    /// $n$-qubit control register that encodes number states $\ket{j}$ in
    /// big-endian format.
    ///
    /// ## target
    /// Generic qubit register that $V_j$ acts on.
    ///
    /// # Remarks
    /// `coefficients` will be padded with identity elements if 
    /// fewer than $2^n$ are specified. This version is implemented 
    /// directly by looping through n-controlled unitary operators.
    operation MultiplexOperationsBruteForceFromGenerator<'T>(unitaryGenerator : (Int, (Int -> ('T => Unit : Adjoint, Controlled))), index: BigEndian, target: 'T) : Unit {
        body (...) {
            let nIndex = Length(index!);
            let nStates = 2^nIndex;
            let (nUnitaries, unitaryFunction) = unitaryGenerator;
            for(idxOp in 0..Microsoft.Quantum.Extensions.Math.MinI(nStates,nUnitaries)-1){
                (ControlledOnInt(idxOp, unitaryFunction(idxOp)))(Reversed(index!),target);
            }
        }
        adjoint auto;
        controlled auto;
        adjoint controlled auto;
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
    /// and the second element `(Int -> ('T => () : Adjoint, Controlled))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$. 
    ///
    /// # Output
    /// A multiply-controlled unitary operation $U$ that applies unitaries
    /// described by `unitaryGenerator`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexOperationsFromGenerator
    function MultiplexerFromGenerator(unitaryGenerator : (Int, (Int -> (Qubit[] => Unit : Adjoint, Controlled)))) : ((BigEndian, Qubit[]) => Unit : Adjoint, Controlled) {
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
    /// and the second element `(Int -> ('T => () : Adjoint, Controlled))`
    /// is a function that takes an integer $j$ in $[0,N-1]$ and outputs the unitary
    /// operation $V_j$. 
    ///
    /// # Output
    /// A multiply-controlled unitary operation $U$ that applies unitaries
    /// described by `unitaryGenerator`.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.MultiplexerBruteForceFromGenerator
    function MultiplexerBruteForceFromGenerator(unitaryGenerator : (Int, (Int -> (Qubit[] => Unit : Adjoint, Controlled)))) : ((BigEndian, Qubit[]) => Unit : Adjoint, Controlled) {
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
    operation LogicalANDMeasAndFix_ (ctrlRegister: Qubit[], target: Qubit) : Unit
    {
        body (...) 
        {
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
            }
            else{
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
                    H(c2);
                    CNOT(c1,c2);
                    H(c2);
                    X(target);
                }
            }
            else{
                 (Controlled X)(ctrlRegister, target);
            }
        }
    }
}

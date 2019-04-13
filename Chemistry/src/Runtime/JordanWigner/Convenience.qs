// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Chemistry.JordanWigner {
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Chemistry;
    open Microsoft.Quantum.Math;
    
    
    // Convenience functions for performing simulation.
    
    /// # Summary
    /// Returns Trotter step operation and the parameters necessary to run it.
    ///
    /// # Input
    /// ## qSharpData
    /// Hamiltonian described by `JordanWignerEncodingData` format.
    /// ## trotterStepSize
    /// Step size of Trotter integrator.
    /// ## trotterOrder
    /// Order of Trotter integrator.
    ///
    /// # Output
    /// A tuple where: `Int` is the number of qubits allocated,
    /// `Double` is `1.0/trotterStepSize`, and the operation
    /// is the Trotter step.
    function TrotterStepOracle (qSharpData : JordanWignerEncodingData, trotterStepSize : Double, trotterOrder : Int) : (Int, (Double, (Qubit[] => Unit : Adjoint, Controlled))) {
        
        let (nSpinOrbitals, data, statePrepData, energyShift) = qSharpData!;
        let generatorSystem = JordanWignerGeneratorSystem(data);
        let evolutionGenerator = EvolutionGenerator(JordanWignerFermionEvolutionSet(), generatorSystem);
        let simulationAlgorithm = (TrotterSimulationAlgorithm(trotterStepSize, trotterOrder))!;
        let oracle = simulationAlgorithm(trotterStepSize, evolutionGenerator, _);
        let nTargetRegisterQubits = nSpinOrbitals;
        let rescaleFactor = 1.0 / trotterStepSize;
        return (nTargetRegisterQubits, (rescaleFactor, oracle));
    }
    
    
    function _QubitizationOracleSeperatedRegisters (qSharpData : JordanWignerEncodingData) : ((Int, Int), (Double, ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled))) {
        
        let (nSpinOrbitals, data, statePrepData, energyShift) = qSharpData!;
        let generatorSystem = JordanWignerBlockEncodingGeneratorSystem(data);
        let (nTerms, genIdxFunction) = generatorSystem!;
        let (oneNorm, blockEncodingReflection) = PauliBlockEncoding(generatorSystem);
        let nTargetRegisterQubits = nSpinOrbitals;
        let nCtrlRegisterQubits = Microsoft.Quantum.Extensions.Math.Ceiling(Lg(ToDouble(nTerms)));
        return ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, QuantumWalkByQubitization(blockEncodingReflection)));
    }
    
    
    /// # Summary
    /// Returns Qubitization operation and the parameters necessary to run it.
    ///
    /// # Input
    /// ## qSharpData
    /// Hamiltonian described by `JordanWignerEncodingData` format.
    ///
    /// # Output
    /// A tuple where: `Int` is the number of qubits allocated,
    /// `Double` is the one-norm of Hamiltonian coefficients, and the operation
    /// is the Quantum walk created by Qubitization.
    function QubitizationOracle (qSharpData : JordanWignerEncodingData) : (Int, (Double, (Qubit[] => Unit : Adjoint, Controlled))) {

        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, oracle)) = (_QubitizationOracleSeperatedRegisters(qSharpData));
        let nQubits = nCtrlRegisterQubits + nTargetRegisterQubits;
        return (nQubits, (oneNorm, _MergeTwoRegisters_(oracle, nTargetRegisterQubits, _)));
    }
    
    
    operation _MergeTwoRegisters_ (oracle : ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled), nSystemQubits : Int, allQubits : Qubit[]) : Unit {
        
        body (...) {
            oracle(allQubits[nSystemQubits .. Length(allQubits) - 1], allQubits[0 .. nSystemQubits - 1]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    function _OptimizedQubitizationOracleSeperatedRegisters_ (qSharpData : JordanWignerEncodingData, targetError : Double) : ((Int, Int), (Double, ((Qubit[], Qubit[]) => Unit : Adjoint, Controlled))) {
        
        let (nSpinOrbitals, data, statePrepData, energyShift) = qSharpData!;
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, blockEncodingReflection)) = (_JordanWignerOptimizedBlockEncoding_(targetError, data, nSpinOrbitals));
        return ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, QuantumWalkByQubitization(blockEncodingReflection)));
    }
    
    
    /// # Summary
    /// Returns T-count optimized Qubitization operation
    /// and the parameters necessary to run it.
    ///
    /// # Input
    /// ## qSharpData
    /// Hamiltonian described by `JordanWignerEncodingData` format.
    /// ## targetError
    /// Error of auxillary state preparation step.
    ///
    /// # Output
    /// A tuple where: `Int` is the number of qubits allocated,
    /// `Double` is the one-norm of Hamiltonian coefficients, and the operation
    /// is the Quantum walk created by Qubitization.
    function OptimizedQubitizationOracle (qSharpData : JordanWignerEncodingData, targetError : Double) : (Int, (Double, (Qubit[] => Unit : Adjoint, Controlled))) {
        
        let ((nCtrlRegisterQubits, nTargetRegisterQubits), (oneNorm, oracle)) = (_OptimizedQubitizationOracleSeperatedRegisters_(qSharpData, targetError));
        let nQubits = nCtrlRegisterQubits + nTargetRegisterQubits;
        return (nQubits, (oneNorm, _MergeTwoRegisters_(oracle, nTargetRegisterQubits, _)));
    }
    
}



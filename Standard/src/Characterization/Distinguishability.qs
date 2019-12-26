// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Arrays;


    operation _ApplyHadamardTest(
        phaseShift : Bool,
        commonPreparation : (Qubit[] => Unit is Adj),
        preparation1 : (Qubit[] => Unit is Adj + Ctl),
        preparation2 : (Qubit[] => Unit is Adj + Ctl),
        control : Qubit,
        target : Qubit[]
    )
    : Unit is Adj
    {
        within {
            H(control);
        } apply {
            commonPreparation(target);
            Controlled preparation1([control], target);
            within { X(control); }
            apply { Controlled preparation2([control], target); }
        }

        (phaseShift ? S | I)(control);
    }

    operation _ApplyHadamardTestOnSingleRegister(
        phaseShift : Bool,
        commonPreparation : (Qubit[] => Unit is Adj),
        preparation1 : (Qubit[] => Unit is Adj + Ctl),
        preparation2 : (Qubit[] => Unit is Adj + Ctl),
        register : Qubit[]
    )
    : Unit is Adj
    {
        let control = Head(register);
        let target = Rest(register);
        _ApplyHadamardTest(
            phaseShift,
            commonPreparation,
            preparation1, preparation2,
            control, target
        );
    }

    /// # Summary
    /// Given two operations which each prepare copies of a state, estimates
    /// the real part of the overlap between the states prepared by each
    /// operation.
    ///
    /// # Input
    /// ## commonPreparation
    /// An operation that prepares a fixed input state.
    /// ## preparation1
    /// The first of the two state preparation operations to be compared.
    /// ## preparation2
    /// The second of the two state preparation operations to be compared.
    /// ## nQubits
    /// The number of qubits on which `commonPreparation`, `preparation1`, and
    /// `preparation2` all act.
    /// ## nMeasurements
    /// The number of measurements to use in estimating the overlap.
    ///
    /// # Remarks
    /// This operation uses the Hadamard test to find the real part of
    /// $$
    /// \begin{align}
    ///     \braket{\psi | V^{\dagger} U | \psi}
    /// \end{align}
    /// $$
    /// where $\ket{\psi}$ is the state prepared by `commonPreparation`,
    /// $U$ is the unitary representation of the action of `preparation1`,
    /// and where $V$ corresponds to `preparation2`.
    ///
    /// # References
    /// TODO
    ///
    /// # See Also
    /// - Microsoft.Quantum.Characterization.EstimateImagOverlapBetweenStates
    /// - Microsoft.Quantum.Characterization.EstimateOverlapBetweenStates
    operation EstimateRealOverlapBetweenStates(
        commonPreparation : (Qubit[] => Unit is Adj),
        preparation1 : (Qubit[] => Unit is Adj + Ctl),
        preparation2 : (Qubit[] => Unit is Adj + Ctl),
        nQubits : Int, nMeasurements : Int
    )
    : Double {
        return 2.0 * EstimateFrequencyA(
            _ApplyHadamardTestOnSingleRegister(false, commonPreparation, preparation1, preparation2, _),
            _HeadMeasurement(nQubits + 1),
            nQubits + 1, nMeasurements
        ) - 1.0;
    }

    operation EstimateImagOverlapBetweenStates(
        commonPreparation : (Qubit[] => Unit is Adj),
        preparation1 : (Qubit[] => Unit is Adj + Ctl),
        preparation2 : (Qubit[] => Unit is Adj + Ctl),
        nQubits : Int, nMeasurements : Int
    )
    : Double {
        return 2.0 * EstimateFrequencyA(
            _ApplyHadamardTestOnSingleRegister(true, commonPreparation, preparation1, preparation2, _),
            _HeadMeasurement(nQubits + 1),
            nQubits + 1, nMeasurements
        ) - 1.0;
    }

    operation _ApplySwapTest(
        preparation1 : (Qubit[] => Unit is Adj),
        preparation2 : (Qubit[] => Unit is Adj),
        control : Qubit,
        target1 : Qubit[],
        target2 : Qubit[]
    )
    : Unit is Adj {
        within {
            H(control);
        } apply {
            preparation1(target1);
            preparation2(target2);
            ApplyToEachCA(Controlled SWAP([control], _), Zip(target1, target2));
        }
    }

    operation _ApplySwapTestOnSingleRegister(
        preparation1 : (Qubit[] => Unit is Adj),
        preparation2 : (Qubit[] => Unit is Adj),
        register : Qubit[]
    )
    : Unit is Adj {
        let control = Head(register);
        let targets = Rest(register);
        _ApplySwapTest(
            preparation1, preparation2,
            control,
            targets[...Length(targets) / 2 - 1],
            targets[Length(targets) / 2...]
        );
    }

    function _HeadMeasurement(nQubits : Int) : (Qubit[] => Result) {
        return Measure(
            ConstantArray(nQubits, PauliI) w/ 0 <- PauliZ,
            _
        );
    }

    operation EstimateOverlapBetweenStates(
        preparation1 : (Qubit[] => Unit is Adj),
        preparation2 : (Qubit[] => Unit is Adj),
        nQubits : Int, nMeasurements : Int
    )
    : Double {
        let nTotalQubits = 2 * nQubits + 1;
        return 2.0 * EstimateFrequencyA(
            _ApplySwapTestOnSingleRegister(preparation1, preparation2, _),
            _HeadMeasurement(nTotalQubits),
            nTotalQubits, nMeasurements
        ) - 1.0;
    }

}
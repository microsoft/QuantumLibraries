// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Characterization {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    

    /// # Summary
    /// Given a preparation and measurement, estimates the frequency
    /// with which that measurement succeeds (returns `Zero`) by
    /// performing a given number of trials.
    ///
    /// # Input
    /// ## preparation
    /// An operation $P$ that prepares a given state $\rho$ on
    /// its input register.
    /// ## measurement
    /// An operation $M$ representing the measurement of interest.
    /// ## nQubits
    /// The number of qubits on which the preparation and measurement
    /// each act.
    /// ## nMeasurements
    /// The number of times that the measurement should be performed
    /// in order to estimate the frequency of interest.
    ///
    /// # Output
    /// An estimate $\hat{p}$ of the frequency with which
    /// $M(P(\ket{00 \cdots 0}\bra{00 \cdots 0}))$ returns `Zero`,
    /// obtained using the unbiased binomial estimator $\hat{p} =
    /// n\_{\uparrow} / n\_{\text{measurements}}$, where $n\_{\uparrow}$ is
    /// the number of `Zero` results observed.
    ///
    /// This is particularly important on target machines which respect
    /// physical limitations, such that probabilities cannot be measured.
    operation EstimateFrequency (preparation : (Qubit[] => Unit), measurement : (Qubit[] => Result), nQubits : Int, nMeasurements : Int) : Double
    {
        mutable nUp = 0;

        for (idxMeasurement in 0 .. nMeasurements - 1) {
            using (register = Qubit[nQubits]) {
                preparation(register);
                let result = measurement(register);

                if (result == Zero) {
                    // NB!!!!! This reverses Zero and One to use conventions
                    //         common in the QCVV community. That is confusing
                    //         but is confusing with an actual purpose.
                    set nUp = nUp + 1;
                }

                // NB: We absolutely must reset here, since preparation()
                //     and measurement() can each use randomness internally.
                ApplyToEach(Reset, register);
            }
        }

        return IntAsDouble(nUp) / IntAsDouble(nMeasurements);
    }

    /// # Summary
    /// Given a preparation that is adjointable and measurement, estimates the frequency
    /// with which that measurement succeeds (returns `Zero`) by
    /// performing a given number of trials.
    ///
    /// # Input
    /// ## preparation
    /// An adjointable operation $P$ that prepares a given state $\rho$ on
    /// its input register.
    /// ## measurement
    /// An operation $M$ representing the measurement of interest.
    /// ## nQubits
    /// The number of qubits on which the preparation and measurement
    /// each act.
    /// ## nMeasurements
    /// The number of times that the measurement should be performed
    /// in order to estimate the frequency of interest.
    ///
    /// # Output
    /// An estimate $\hat{p}$ of the frequency with which
    /// $M(P(\ket{00 \cdots 0}\bra{00 \cdots 0}))$ returns `Zero`,
    /// obtained using the unbiased binomial estimator $\hat{p} =
    /// n\_{\uparrow} / n\_{\text{measurements}}$, where $n\_{\uparrow}$ is
    /// the number of `Zero` results observed.
    ///
    /// For adjointable operations, certain assumptions can be made such like
    /// calling the operation will prepare the qubits to exactly the same state,
    /// which allow target machines to make some performance optimizations.
    operation EstimateFrequencyA (preparation : (Qubit[] => Unit is Adj), measurement : (Qubit[] => Result), nQubits : Int, nMeasurements : Int) : Double
    {
        return EstimateFrequency(preparation, measurement, nQubits, nMeasurements);
    }
}

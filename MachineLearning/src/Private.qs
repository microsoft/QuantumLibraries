// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.MachineLearning {
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;

    internal function AllNearlyEqualD(v1 : Double[], v2 : Double[]) : Bool {
        return Length(v1) == Length(v2) and All(NearlyEqualD, Zipped(v1, v2));
    }

    internal function TailMeasurement(nQubits : Int) : (Qubit[] => Result) {
        let paulis = [PauliI, size = nQubits] w/ (nQubits - 1) <- PauliZ;
        return Measure(paulis, _);
    }

}

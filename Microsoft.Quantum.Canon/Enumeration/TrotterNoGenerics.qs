// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Implementation of first-order Trotter-Suzuki integrator.
    ///
    /// # Input
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated. See @"microsoft.quantum.canon.evolutiongenerator" for input format.
    /// ## stepSize
    /// Multiplier on size of each step of the simulation.
    operation Trotter1ImplCA(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), stepSize : Double, target : Qubit[]) : () {
        body {
            let (nSteps, op) = evolutionGenerator;
            for(idx in 0..nSteps-1){
                op(idx, stepSize, target);
            }
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Implementation of second-order Trotter-Suzuki integrator.
    ///
    /// # Input
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated. See @"microsoft.quantum.canon.evolutiongenerator" for input format.
    /// ## stepSize
    /// Multiplier on size of each step of the simulation.
    operation Trotter2ImplCA(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), stepSize : Double, target : Qubit[]) : () {
        body {
            let (nSteps, op) = evolutionGenerator;
            for(idx in 0..nSteps-1){
                op(idx, stepSize * 0.5, target);
            }
            for(idx in (nSteps-1)..(-1)..0){
                op(idx, stepSize * 0.5, target);
            }
        }
        adjoint auto
        controlled auto
        controlled adjoint auto
    }

    /// # Summary
    /// Function that returns unitary implementation of Trotter-Suzuki integrator.
    ///
    /// # Input
    /// ## evolutionGenerator
    /// A complete description of the system to be simulated. See @"microsoft.quantum.canon.evolutiongenerator" for input format.
    /// ## trotterOrder
    /// Selects order of Trotter-Suzuki integrator. Order 1 and 2 are currently supported.
    ///
    /// # Output
    /// Returns a unitary implementing the Trotter-Suzuki integrator, where the first parameter `Double` is the integration step size, and the second parameter `Qubit[]` is the register acted upon.
    function DecomposeIntoTimeStepsCA(evolutionGenerator : (Int, ((Int, Double, Qubit[]) => () : Adjoint, Controlled)), trotterOrder : Int) : ((Double, Qubit[]) => () : Adjoint, Controlled) {
        if (trotterOrder == 1) {
            return Trotter1ImplCA(evolutionGenerator, _, _);
        } elif (trotterOrder == 2) {
            return Trotter2ImplCA(evolutionGenerator, _, _);
        } else {
            fail "Order $order not yet supported.";
        }

        // Needed so we have a return value of the right type in all cases, but
        //        this line is unreachable.
        return Trotter1ImplCA(evolutionGenerator, _, _);
    }
}

// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Implementation of the first-order Trotter–Suzuki integrator.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type which each time step should act upon; typically, either
    /// `Qubit[]` or `Qubit`.
    ///
    /// # Input
    /// ### nSteps
    /// The number of operations to be decomposed into time steps.
    /// ### op
    /// An operation which accepts an index input (type `Int`) and a time
    /// input (type `Double`) and a quantum register (type `'T`) for decomposition.
    /// ## stepSize
    /// Multiplier on size of each step of the simulation.
    /// ## target
    /// A quantum register on which the operations act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// op(0, deltaT, target);
    /// op(1, deltaT, target);
    /// ```
	/// and
	/// ```qsharp
    /// _Trotter1ImplCA((2, op), deltaT, target);
    /// ```
    operation _Trotter1ImplCA<'T> ((nSteps : Int, op : ((Int, Double, 'T) => Unit is Adj + Ctl)), stepSize : Double, target : 'T) : Unit is Adj + Ctl {
        for (idx in 0 .. nSteps - 1)
        {
            op(idx, stepSize, target);
        }
    }


    /// # Summary
    /// Implementation of the second-order Trotter–Suzuki integrator.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type which each time step should act upon; typically, either
    /// `Qubit[]` or `Qubit`.
    ///
    /// # Input
    /// ### nSteps
    /// The number of operations to be decomposed into time steps.
    /// ### op
    /// An operation which accepts an index input (type `Int`) and a time
    /// input (type `Double`) and a quantum register (type `'T`) for decomposition.
    /// ## stepSize
    /// Multiplier on size of each step of the simulation.
    /// ## target
    /// A quantum register on which the operations act.
    ///
    /// # Remarks
    /// ## Example
    /// The following are equivalent:
    /// ```qsharp
    /// op(0, deltaT / 2.0, target);
    /// op(1, deltaT / 2.0, target);
    /// op(1, deltaT / 2.0, target);
    /// op(0, deltaT / 2.0, target);
    /// ```
	/// and
	/// ```qsharp
    /// _Trotter2ImplCA((2, op), deltaT, target);
    /// ```
    operation _Trotter2ImplCA<'T> ((nSteps : Int, op : ((Int, Double, 'T) => Unit is Adj + Ctl)), stepSize : Double, target : 'T) : Unit is Adj + Ctl
    {
        for (idx in 0 .. nSteps - 1)
        {
            op(idx, stepSize * 0.5, target);
        }
        
        for (idx in nSteps - 1 .. -1 .. 0)
        {
            op(idx, stepSize * 0.5, target);
        }
    }

    /// # Summary
    /// Recursive implementation of even-order Trotter–Suzuki integrator.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type which each time step should act upon; typically, either
    /// `Qubit[]` or `Qubit`.
    ///
    /// # Input
    /// ### Order
    /// Order of Trotter-Suzuki integrator.
    /// ### nSteps
    /// The number of operations to be decomposed into time steps.
    /// ### op
    /// An operation which accepts an index input (type `Int`) and a time
    /// input (type `Double`) and a quantum register (type `'T`) for decomposition.
    /// ## stepSize
    /// Multiplier on size of each step of the simulation.
    /// ## target
    /// A quantum register on which the operations act.
    operation _TrotterArbitraryImplCA<'T> (order:Int, (nSteps : Int, op : ((Int, Double, 'T) => Unit is Adj + Ctl)), stepSize : Double, target : 'T)
    : Unit is Adj + Ctl {
        if(order > 2){
            let stepSizeOuter = _TrotterStepSize(order);
            let stepSizeInner = 1.0 - 4.0 * stepSizeOuter;
            _TrotterArbitraryImplCA(order -2, (nSteps, op), stepSizeOuter * stepSize, target);
            _TrotterArbitraryImplCA(order -2, (nSteps, op), stepSizeOuter * stepSize, target);
            _TrotterArbitraryImplCA(order -2, (nSteps, op), stepSizeInner * stepSize, target);
            _TrotterArbitraryImplCA(order -2, (nSteps, op), stepSizeOuter * stepSize, target);
            _TrotterArbitraryImplCA(order -2, (nSteps, op), stepSizeOuter * stepSize, target);
        }
        elif(order == 2){
            _Trotter2ImplCA((nSteps, op), stepSize, target);
        }
        else{
            _Trotter1ImplCA((nSteps, op), stepSize, target);
        }
    }

    /// # Summary
    /// Computes Trotter step size in recursive implementation of
    /// Trotter simulation algorithm.
    function _TrotterStepSize (order: Int) : Double {
        return 1.0 / (4.0 - PowD(4.0, 1.0 / (2.0 * IntAsDouble(order) - 1.0)));
    }


    /// # Summary
    /// Returns an operation implementing the Trotter–Suzuki integrator for
    /// a given operation.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type which each time step should act upon; typically, either
    /// `Qubit[]` or `Qubit`.
    ///
    /// # Input
    /// ## nSteps
    /// The number of operations to be decomposed into time steps.
    /// ## op
    /// An operation which accepts an index input (type `Int`) and a time
    /// input (type `Double`) for decomposition.
    /// ## trotterOrder
    /// Selects the order of the Trotter–Suzuki integrator to be used.
    /// Order 1 and even orders 2,4,6,... are currently supported.
    ///
    /// # Output
    /// Returns a unitary implementing the Trotter–Suzuki integrator, where
    /// the first parameter `Double` is the integration step size, and the
    /// second parameter is the target acted upon.
    ///
    /// # References
    /// We use the implementation in
    /// - [ *D. W. Berry, G. Ahokas, R. Cleve, B. C. Sanders* ](https://arxiv.org/abs/quant-ph/0508139)
    function DecomposeIntoTimeStepsCA<'T> ((nSteps : Int, op : ((Int, Double, 'T) => Unit is Adj + Ctl)), trotterOrder : Int) : ((Double, 'T) => Unit is Adj + Ctl)
    {
        if (trotterOrder == 1)
        {
            return _Trotter1ImplCA((nSteps, op), _, _);
        }
        elif (trotterOrder == 2)
        {
            return _Trotter2ImplCA((nSteps, op), _, _);
        }
        elif(trotterOrder % 2 ==0){
            return _TrotterArbitraryImplCA(trotterOrder, (nSteps, op), _, _);
        }
        else
        {
            fail $"Odd order {trotterOrder} not yet supported.";
        }
    }

}



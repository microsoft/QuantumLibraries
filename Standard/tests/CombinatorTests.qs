// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    
    
    function ComposeTest () : Unit {
        let target = [3, 17, 2];
        EqualityFactI((Compose(Modulus(_, 14), Max))(target), 3, $"Compose(Modulus(_, 14), Max) did not return expected result.");
    }
    
    
    operation WithTest () : Unit {
        
        let actual = ApplyWith(H, X, _);
        let expected = Z;
        AssertOperationsEqualReferenced(4, ApplyToEach(actual, _), ApplyToEachA(expected, _));
    }
    
    
    // Make sure that if CurryTest fails, it's because of Curry and not
    // something else.
    operation CurryPreTest () : Unit {
        
        AssertOperationsEqualInPlace(Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _), 1);
        AssertOperationsEqualReferenced(1, Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _));
    }
    
    
    operation CurryTest () : Unit {
        
        let curried = CurryOp(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(curried(1.7), Exp([PauliZ], 1.7, _), 1);
        AssertOperationsEqualReferenced(1, curried(1.7), Exp([PauliZ], 1.7, _));
    }
    
    
    operation BindTest () : Unit {
        
        let bound = BindCA([H, X, H]);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(Z, _));
    }
    
    
    function StripControlled<'T> (op : ('T => Unit : Adjoint, Controlled)) : ('T => Unit : Adjoint) {
        
        return op;
    }
    
    
    operation BindATest () : Unit {
        
        let bound = BindA(Mapped(StripControlled<Qubit>, [T, T]));
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        AssertOperationsEqualReferenced(3, ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _));
    }
    
    
    operation BindCTestHelper0 (op : (Qubit => Unit : Controlled), qubits : Qubit[]) : Unit {
        
        Controlled op([qubits[0]], qubits[1]);
    }
    
    
    operation BindCTestHelper1 (op : (Qubit => Unit : Adjoint, Controlled), qubits : Qubit[]) : Unit {
        
        body (...) {
            Controlled op([qubits[0]], qubits[1]);
        }
        
        adjoint invert;
    }
    
    
    function StripAdjoint<'T> (op : ('T => Unit : Adjoint, Controlled)) : ('T => Unit : Controlled) {
        
        return op;
    }
    
    
    operation BindCTest () : Unit {
        
        let stripped = Mapped(StripAdjoint<Qubit>, [T, T]);
        let bound = BindC(stripped);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        let op = BindCTestHelper0(bound, _);
        let target = BindCTestHelper1(S, _);
        AssertOperationsEqualReferenced(6, op, target);
    }
    
    
    operation BindCATest () : Unit {
        let bound = BindCA([T, T]);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        AssertOperationsEqualReferenced(3, ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _));
        let op = BindCTestHelper0(Adjoint bound, _);
        let target = BindCTestHelper1(Adjoint S, _);
        AssertOperationsEqualReferenced(4, op, target);
    }
    
    
    operation OperationPowTest () : Unit {
        
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(H, 2), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(Z, 2), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(S, 4), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(T, 8), _), NoOp<Qubit[]>);
    }
    
    
    operation ApplyToSubregisterTest () : Unit {
        
        let bigOp = ApplyPauli([PauliI, PauliX, PauliY, PauliZ, PauliI], _);
        let smallOp = ApplyPauli([PauliX, PauliY, PauliZ], _);
        AssertOperationsEqualReferenced(5, ApplyToSubregister(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictToSubregister(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterC(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictToSubregisterC(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterA(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictToSubregisterA(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterCA(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictToSubregisterCA(smallOp, [1, 2, 3]), bigOp);
    }
    
    
    operation CControlledExpected (op : (Qubit => Unit : Adjoint, Controlled), target : Qubit[]) : Unit {
        
        body (...) {
            op(target[0]);
            op(target[2]);
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    operation CControlledActual (op : (Qubit => Unit), target : Qubit[]) : Unit {
        
        ApplyToEach(CControlled(op), Zip([true, false, true], target));
    }
    
    
    operation CControlledActualC (op : (Qubit => Unit : Controlled), target : Qubit[]) : Unit {
        
        body (...) {
            ApplyToEachC(CControlledC(op), Zip([true, false, true], target));
        }
        
        controlled distribute;
    }
    
    
    operation CControlledActualA (op : (Qubit => Unit : Adjoint), target : Qubit[]) : Unit {
        
        body (...) {
            ApplyToEachA(CControlledA(op), Zip([true, false, true], target));
        }
        
        adjoint invert;
    }
    
    
    operation CControlledActualCA (op : (Qubit => Unit : Adjoint, Controlled), target : Qubit[]) : Unit {
        
        body (...) {
            ApplyToEachCA(CControlledCA(op), Zip([true, false, true], target));
        }
        
        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }
    
    
    operation CControlledTest () : Unit {
        
        AssertOperationsEqualReferenced(3, CControlledActual(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActual(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActual(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActual(T, _), CControlledExpected(T, _));
    }
    
    
    operation CControlledTestC () : Unit {
        
        AssertOperationsEqualReferenced(3, CControlledActualC(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(T, _), CControlledExpected(T, _));
    }
    
    
    operation CControlledTestA () : Unit {
        
        AssertOperationsEqualReferenced(3, CControlledActualA(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(T, _), CControlledExpected(T, _));
    }
    
    
    operation CControlledTestCA () : Unit {
        
        AssertOperationsEqualReferenced(3, CControlledActualCA(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(T, _), CControlledExpected(T, _));
    }
    
}



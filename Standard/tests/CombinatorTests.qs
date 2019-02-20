// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Microsoft.Quantum.Tests {
    
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Testing;
    
    
    function ComposeTest () : Unit {
        
        let target = [3, 17, 2];
        AssertIntEqual((Compose(Modulus(_, 14), Max))(target), 3, $"Compose(Modulus(_, 14), Max) did not return expected result.");
    }
    
    
    operation WithTest () : Unit {
        
        let actual = With(H, X, _);
        let expected = Z;
        AssertOperationsEqualReferenced(ApplyToEach(actual, _), ApplyToEachA(expected, _), 4);
    }
    
    
    // Make sure that if CurryTest fails, it's because of Curry and not
    // something else.
    operation CurryPreTest () : Unit {
        
        AssertOperationsEqualInPlace(Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _), 1);
        AssertOperationsEqualReferenced(Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _), 1);
    }
    
    
    operation CurryTest () : Unit {
        
        let curried = CurryOp(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(curried(1.7), Exp([PauliZ], 1.7, _), 1);
        AssertOperationsEqualReferenced(curried(1.7), Exp([PauliZ], 1.7, _), 1);
    }
    
    
    operation BindTest () : Unit {
        
        let bound = BindCA([H, X, H]);
        AssertOperationsEqualReferenced(ApplyToEach(bound, _), ApplyToEachA(Z, _), 3);
    }
    
    
    function StripControlled<'T> (op : ('T => Unit : Adjoint, Controlled)) : ('T => Unit : Adjoint) {
        
        return op;
    }
    
    
    operation BindATest () : Unit {
        
        let bound = BindA(Map(StripControlled<Qubit>, [T, T]));
        AssertOperationsEqualReferenced(ApplyToEach(bound, _), ApplyToEachA(S, _), 3);
        AssertOperationsEqualReferenced(ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _), 3);
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
        
        let stripped = Map(StripAdjoint<Qubit>, [T, T]);
        let bound = BindC(stripped);
        AssertOperationsEqualReferenced(ApplyToEach(bound, _), ApplyToEachA(S, _), 3);
        let op = BindCTestHelper0(bound, _);
        let target = BindCTestHelper1(S, _);
        AssertOperationsEqualReferenced(op, target, 6);
    }
    
    
    operation BindCATest () : Unit {
        
        let bound = BindCA([T, T]);
        AssertOperationsEqualReferenced(ApplyToEach(bound, _), ApplyToEachA(S, _), 3);
        AssertOperationsEqualReferenced(ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _), 3);
        let op = BindCTestHelper0(Adjoint bound, _);
        let target = BindCTestHelper1(Adjoint S, _);
        AssertOperationsEqualReferenced(op, target, 4);
    }
    
    
    operation OperationPowTest () : Unit {
        
        AssertOperationsEqualReferenced(ApplyToEach(OperationPow(H, 2), _), NoOp<Qubit[]>, 3);
        AssertOperationsEqualReferenced(ApplyToEach(OperationPow(Z, 2), _), NoOp<Qubit[]>, 3);
        AssertOperationsEqualReferenced(ApplyToEach(OperationPow(S, 4), _), NoOp<Qubit[]>, 3);
        AssertOperationsEqualReferenced(ApplyToEach(OperationPow(T, 8), _), NoOp<Qubit[]>, 3);
    }
    
    
    operation ApplyToSubregisterTest () : Unit {
        
        let bigOp = ApplyPauli([PauliI, PauliX, PauliY, PauliZ, PauliI], _);
        let smallOp = ApplyPauli([PauliX, PauliY, PauliZ], _);
        AssertOperationsEqualReferenced(ApplyToSubregister(smallOp, [1, 2, 3], _), bigOp, 5);
        AssertOperationsEqualReferenced(RestrictToSubregister(smallOp, [1, 2, 3]), bigOp, 5);
        AssertOperationsEqualReferenced(ApplyToSubregisterC(smallOp, [1, 2, 3], _), bigOp, 5);
        AssertOperationsEqualReferenced(RestrictToSubregisterC(smallOp, [1, 2, 3]), bigOp, 5);
        AssertOperationsEqualReferenced(ApplyToSubregisterA(smallOp, [1, 2, 3], _), bigOp, 5);
        AssertOperationsEqualReferenced(RestrictToSubregisterA(smallOp, [1, 2, 3]), bigOp, 5);
        AssertOperationsEqualReferenced(ApplyToSubregisterCA(smallOp, [1, 2, 3], _), bigOp, 5);
        AssertOperationsEqualReferenced(RestrictToSubregisterCA(smallOp, [1, 2, 3]), bigOp, 5);
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
        
        AssertOperationsEqualReferenced(CControlledActual(H, _), CControlledExpected(H, _), 3);
        AssertOperationsEqualReferenced(CControlledActual(Z, _), CControlledExpected(Z, _), 3);
        AssertOperationsEqualReferenced(CControlledActual(S, _), CControlledExpected(S, _), 3);
        AssertOperationsEqualReferenced(CControlledActual(T, _), CControlledExpected(T, _), 3);
    }
    
    
    operation CControlledTestC () : Unit {
        
        AssertOperationsEqualReferenced(CControlledActualC(H, _), CControlledExpected(H, _), 3);
        AssertOperationsEqualReferenced(CControlledActualC(Z, _), CControlledExpected(Z, _), 3);
        AssertOperationsEqualReferenced(CControlledActualC(S, _), CControlledExpected(S, _), 3);
        AssertOperationsEqualReferenced(CControlledActualC(T, _), CControlledExpected(T, _), 3);
    }
    
    
    operation CControlledTestA () : Unit {
        
        AssertOperationsEqualReferenced(CControlledActualA(H, _), CControlledExpected(H, _), 3);
        AssertOperationsEqualReferenced(CControlledActualA(Z, _), CControlledExpected(Z, _), 3);
        AssertOperationsEqualReferenced(CControlledActualA(S, _), CControlledExpected(S, _), 3);
        AssertOperationsEqualReferenced(CControlledActualA(T, _), CControlledExpected(T, _), 3);
    }
    
    
    operation CControlledTestCA () : Unit {
        
        AssertOperationsEqualReferenced(CControlledActualCA(H, _), CControlledExpected(H, _), 3);
        AssertOperationsEqualReferenced(CControlledActualCA(Z, _), CControlledExpected(Z, _), 3);
        AssertOperationsEqualReferenced(CControlledActualCA(S, _), CControlledExpected(S, _), 3);
        AssertOperationsEqualReferenced(CControlledActualCA(T, _), CControlledExpected(T, _), 3);
    }
    
}



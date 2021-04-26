// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Tests {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    @Test("QuantumSimulator")
    function ComposeIsCorrect() : Unit {
        let target = [3, 17, 2];
        EqualityFactI((Compose(ModulusI(_, 14), Max))(target), 3, $"Compose(ModulusI(_, 14), Max) did not return expected result.");
    }

    @Test("QuantumSimulator")
    operation TestWith() : Unit {

        let actual = ApplyWith(H, X, _);
        let expected = Z;
        AssertOperationsEqualReferenced(4, ApplyToEach(actual, _), ApplyToEachA(expected, _));
    }


    // Make sure that if CurryTest fails, it's because of Curry and not
    // something else.
    @Test("QuantumSimulator")
    operation PreTestCurried() : Unit {

        AssertOperationsEqualInPlace(1, Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _));
        AssertOperationsEqualReferenced(1, Exp([PauliZ], 1.7, _), Exp([PauliZ], 1.7, _));
    }

    @Test("QuantumSimulator")
    operation TestCurried() : Unit {
        let curried = CurriedOp(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(1, curried(1.7), Exp([PauliZ], 1.7, _));
        AssertOperationsEqualReferenced(1, curried(1.7), Exp([PauliZ], 1.7, _));
    }

    @Test("QuantumSimulator")
    operation TestCurriedA() : Unit {
        let curried = CurriedOpA(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(1, curried(1.7), Exp([PauliZ], 1.7, _));
        AssertOperationsEqualReferenced(1, curried(1.7), Exp([PauliZ], 1.7, _));
    }

    @Test("QuantumSimulator")
    operation TestCurriedC() : Unit {
        let curried = CurriedOpC(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(1, curried(1.7), Exp([PauliZ], 1.7, _));
        AssertOperationsEqualReferenced(1, curried(1.7), Exp([PauliZ], 1.7, _));
    }

    @Test("QuantumSimulator")
    operation TestCurriedCA() : Unit {
        let curried = CurriedOpCA(Exp([PauliZ], _, _));
        AssertOperationsEqualInPlace(1, curried(1.7), Exp([PauliZ], 1.7, _));
        AssertOperationsEqualReferenced(1, curried(1.7), Exp([PauliZ], 1.7, _));
    }

    @Test("QuantumSimulator")
    operation TestBind() : Unit {

        let bound = BoundCA([H, X, H]);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(Z, _));
    }


    function StripControlled<'T> (op : ('T => Unit is Adj + Ctl)) : ('T => Unit is Adj) {

        return op;
    }

    @Test("QuantumSimulator")
    operation TestBindA() : Unit {

        let bound = BoundA(Mapped(StripControlled<Qubit>, [T, T]));
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        AssertOperationsEqualReferenced(3, ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _));
    }


    operation BindCTestHelper0 (op : (Qubit => Unit is Ctl), qubits : Qubit[]) : Unit {

        Controlled op([qubits[0]], qubits[1]);
    }


    operation BindCTestHelper1 (op : (Qubit => Unit is Adj + Ctl), qubits : Qubit[]) : Unit {
        body (...) {
            Controlled op([qubits[0]], qubits[1]);
        }

        adjoint invert;
    }


    function StripAdjoint<'T> (op : ('T => Unit is Adj + Ctl)) : ('T => Unit is Ctl) {
        return op;
    }

    @Test("QuantumSimulator")
    operation TestBindC() : Unit {

        let stripped = Mapped(StripAdjoint<Qubit>, [T, T]);
        let bound = BoundC(stripped);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        let op = BindCTestHelper0(bound, _);
        let target = BindCTestHelper1(S, _);
        AssertOperationsEqualReferenced(6, op, target);
    }

    @Test("QuantumSimulator")
    operation TestBindCA() : Unit {
        let bound = BoundCA([T, T]);
        AssertOperationsEqualReferenced(3, ApplyToEach(bound, _), ApplyToEachA(S, _));
        AssertOperationsEqualReferenced(3, ApplyToEach(Adjoint bound, _), ApplyToEachA(Adjoint S, _));
        let op = BindCTestHelper0(Adjoint bound, _);
        let target = BindCTestHelper1(Adjoint S, _);
        AssertOperationsEqualReferenced(4, op, target);
    }

    @Test("QuantumSimulator")
    operation TestOperationPow() : Unit {

        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(H, 2), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(Z, 2), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(S, 4), _), NoOp<Qubit[]>);
        AssertOperationsEqualReferenced(3, ApplyToEach(OperationPow(T, 8), _), NoOp<Qubit[]>);
    }

    @Test("QuantumSimulator")
    operation TestApplyToSubregister() : Unit {

        let bigOp = ApplyPauli([PauliI, PauliX, PauliY, PauliZ, PauliI], _);
        let smallOp = ApplyPauli([PauliX, PauliY, PauliZ], _);
        AssertOperationsEqualReferenced(5, ApplyToSubregister(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictedToSubregister(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterC(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictedToSubregisterC(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterA(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictedToSubregisterA(smallOp, [1, 2, 3]), bigOp);
        AssertOperationsEqualReferenced(5, ApplyToSubregisterCA(smallOp, [1, 2, 3], _), bigOp);
        AssertOperationsEqualReferenced(5, RestrictedToSubregisterCA(smallOp, [1, 2, 3]), bigOp);
    }


    operation CControlledExpected (op : (Qubit => Unit is Adj + Ctl), target : Qubit[]) : Unit is Adj + Ctl {
        op(target[0]);
        op(target[2]);
    }


    operation CControlledActual (op : (Qubit => Unit), target : Qubit[]) : Unit {

        ApplyToEach(CControlled(op), Zipped([true, false, true], target));
    }


    operation CControlledActualC (op : (Qubit => Unit is Ctl), target : Qubit[]) : Unit {

        body (...) {
            ApplyToEachC(CControlledC(op), Zipped([true, false, true], target));
        }

        controlled distribute;
    }


    operation CControlledActualA (op : (Qubit => Unit is Adj), target : Qubit[]) : Unit {

        body (...) {
            ApplyToEachA(CControlledA(op), Zipped([true, false, true], target));
        }

        adjoint invert;
    }


    operation CControlledActualCA (op : (Qubit => Unit is Adj + Ctl), target : Qubit[]) : Unit {

        body (...) {
            ApplyToEachCA(CControlledCA(op), Zipped([true, false, true], target));
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    @Test("QuantumSimulator")
    operation TestCControlled() : Unit {

        AssertOperationsEqualReferenced(3, CControlledActual(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActual(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActual(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActual(T, _), CControlledExpected(T, _));
    }

    @Test("QuantumSimulator")
    operation CControlledTestC () : Unit {

        AssertOperationsEqualReferenced(3, CControlledActualC(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualC(T, _), CControlledExpected(T, _));
    }

    @Test("QuantumSimulator")
    operation CControlledTestA () : Unit {

        AssertOperationsEqualReferenced(3, CControlledActualA(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualA(T, _), CControlledExpected(T, _));
    }

    @Test("QuantumSimulator")
    operation CControlledTestCA () : Unit {

        AssertOperationsEqualReferenced(3, CControlledActualCA(H, _), CControlledExpected(H, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(Z, _), CControlledExpected(Z, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(S, _), CControlledExpected(S, _));
        AssertOperationsEqualReferenced(3, CControlledActualCA(T, _), CControlledExpected(T, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfZero() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfZero(One, (ApplyToEach(H, _), _)), ApplyToEachA(I, _));
        AssertOperationsEqualReferenced(2, ApplyIfZero(Zero, (ApplyToEach(H, _), _)), ApplyToEachA(H, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfOne() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfOne(One, (ApplyToEach(H, _), _)), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfOne(Zero, (ApplyToEach(H, _), _)), ApplyToEachA(I, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfZeroC() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfZeroC(One, (ApplyToEachC(H, _), _)), ApplyToEachA(I, _));
        AssertOperationsEqualReferenced(2, ApplyIfZeroC(Zero, (ApplyToEachC(H, _), _)), ApplyToEachA(H, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfOneC() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfOneC(One, (ApplyToEachC(H, _), _)), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfOneC(Zero, (ApplyToEachC(H, _), _)), ApplyToEachA(I, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfZeroCA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfZeroCA(One, (ApplyToEachCA(H, _), _)), ApplyToEachA(I, _));
        AssertOperationsEqualReferenced(2, ApplyIfZeroCA(Zero, (ApplyToEachCA(H, _), _)), ApplyToEachA(H, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfOneCA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfOneCA(One, (ApplyToEachCA(H, _), _)), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfOneCA(Zero, (ApplyToEachCA(H, _), _)), ApplyToEachA(I, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfZeroA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfZeroA(One, (ApplyToEachA(H, _), _)), ApplyToEachA(I, _));
        AssertOperationsEqualReferenced(2, ApplyIfZeroA(Zero, (ApplyToEachA(H, _), _)), ApplyToEachA(H, _));
    }

    @Test("QuantumSimulator")
    operation TestApplyIfOneA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfOneA(One, (ApplyToEachA(H, _), _)), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfOneA(Zero, (ApplyToEachA(H, _), _)), ApplyToEachA(I, _));
    }

    operation ApplyIfElseRCase(result : Result, register : Qubit[]) : Unit {
        ApplyIfElseR(
            result,
            (ApplyToEach(H, _), register),
            (ApplyToEach(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseR() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseRCase(Zero, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseRCase(One, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseRACase(result : Result, register : Qubit[]) : Unit is Adj {
        ApplyIfElseRA(
            result,
            (ApplyToEachA(H, _), register),
            (ApplyToEachA(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseRA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseRACase(Zero, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseRACase(One, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseRCCase(result : Result, register : Qubit[]) : Unit is Ctl {
        ApplyIfElseRC(
            result,
            (ApplyToEachC(H, _), register),
            (ApplyToEachC(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseRC() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseRCCase(Zero, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseRCCase(One, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseRCACase(result : Result, register : Qubit[]) : Unit is Adj + Ctl {
        ApplyIfElseRCA(
            result,
            (ApplyToEachCA(H, _), register),
            (ApplyToEachCA(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseRCA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseRCACase(Zero, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseRCACase(One, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseBCase(bit : Bool, register : Qubit[]) : Unit {
        ApplyIfElseB(
            bit,
            (ApplyToEach(H, _), register),
            (ApplyToEach(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseB() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseBCase(true, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseBCase(false, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseBACase(bit : Bool, register : Qubit[]) : Unit is Adj {
        ApplyIfElseBA(
            bit,
            (ApplyToEachA(H, _), register),
            (ApplyToEachA(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseBA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseBACase(true, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseBACase(false, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseBCCase(bit : Bool, register : Qubit[]) : Unit is Ctl {
        ApplyIfElseBC(
            bit,
            (ApplyToEachC(H, _), register),
            (ApplyToEachC(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseBC() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseBCCase(true, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseBCCase(false, _), ApplyToEachA(X, _));
    }

    operation ApplyIfElseBCACase(bit : Bool, register : Qubit[]) : Unit is Adj + Ctl {
        ApplyIfElseBCA(
            bit,
            (ApplyToEachCA(H, _), register),
            (ApplyToEachCA(X, _), register)
        );
    }

    @Test("QuantumSimulator")
    operation TestApplyIfElseBCA() : Unit {
        AssertOperationsEqualReferenced(2, ApplyIfElseBCACase(true, _), ApplyToEachA(H, _));
        AssertOperationsEqualReferenced(2, ApplyIfElseBCACase(false, _), ApplyToEachA(X, _));
    }

    operation ApplyXToSecondQubit(qubits : Qubit[]) : Unit is Adj + Ctl {
        X(qubits[1]);
    }

    @Test("QuantumSimulator")
    operation TestApplyToElement() : Unit {
        AssertOperationsEqualReferenced(3,
            ApplyToElement(X, 1, _),
            ApplyXToSecondQubit
        );
        
        AssertOperationsEqualReferenced(3,
            ApplyToElementC(X, 1, _),
            ApplyXToSecondQubit
        );
        
        AssertOperationsEqualReferenced(3,
            ApplyToElementA(X, 1, _),
            ApplyXToSecondQubit
        );
        
        AssertOperationsEqualReferenced(3,
            ApplyToElementCA(X, 1, _),
            ApplyXToSecondQubit
        );
    }

}

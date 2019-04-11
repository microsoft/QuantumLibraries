// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arithmetic {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    /// # Summary
    /// Performs a modular increment of a qubit register by an integer constant.
    ///
    /// Let us denote `increment` by a, `modulus` by N and integer encoded in `target` by y
    /// Then the operation performs the following transformation:
    /// \begin{align}
    ///     \ket{y} \mapsto \ket{y + 1 \operatorname{mod} N}
    /// \end{align}
    /// Integers are encoded in little-endian format.
    ///
    /// # Input
    /// ## increment
    /// Integer increment a to be added to y.
    /// ## modulus
    /// Integer N that mods y + a.
    /// ## target
    /// Integer y in `LittleEndian` format that `increment` a is added to.
    ///
    /// # See Also
    /// - IncrementPhaseByModularInteger
    ///
    /// # Remarks
    /// Assumes that the value of target is less than N. Note that
    /// <xref:microsoft.quantum.arithmetic.incrementphasebymodularinteger> implements
    /// the same operation, but in the `PhaseLittleEndian` basis.
    operation IncrementByModularInteger(increment : Int, modulus : Int, target : LittleEndian) : Unit {
        body (...) {
            let inner = IncrementPhaseByModularInteger(increment, modulus, _);

            using (extraZeroBit = Qubit()) {
                ApplyPhaseLEOperationOnLECA(inner, LittleEndian(target! + [extraZeroBit]));
            }
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// Performs a modular increment of a qubit register by an integer constant.
    ///
    /// Let us denote `increment` by a, `modulus` by N and integer encoded in `target` by y
    /// Then the operation performs the following transformation:
    /// |y⟩ ↦ |y+a (mod N)⟩
    /// Integers are encoded in little-endian format in QFT basis
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ModularIncrementLE
    ///
    /// # Remarks
    /// Assumes that `target` has the highest bit set to 0.
    /// Also assumes that the value of target is less than N.
    ///
    /// For the circuit diagram and explanation see Figure 5 on [Page 5
    /// of arXiv:quant-ph/0205095v3](https://arxiv.org/pdf/quant-ph/0205095v3.pdf#page=5).
    operation IncrementPhaseByModularInteger(increment : Int, modulus : Int, target : PhaseLittleEndian) : Unit {
        body (...) {
            Controlled IncrementPhaseByModularInteger(new Qubit[0], (increment, modulus, target));
        }
        adjoint auto;

        controlled (controls, ...) {
            Fact(modulus <= 2 ^ (Length(target!) - 1), $"`multiplier` must be big enough to fit integers modulo `modulus`" + $"with highest bit set to 0");

            if (_EnableExtraAssertsForArithmetic()) {
                // assert that the highest bit is zero, by switching to computational basis
                ApplyLEOperationOnPhaseLEA(AssertHighestBit(Zero, _), target);

                // check that the input is less than modulus
                AssertPhaseLessThan(modulus, target);
            }

            // note that controlled version is correct only under the assumption
            // that the value of target is less than modulus
            using (lessThanModulusFlag = Qubit()) {
                let copyMostSignificantBitPhaseLE = ApplyLEOperationOnPhaseLEA(CopyMostSignificantBit(_, lessThanModulusFlag), _);

                // lets track the state of target register through the computation
                Controlled IncrementPhaseByInteger(controls, (increment, target));

                // the state is |x+a⟩ in QFT basis
                Adjoint IncrementPhaseByInteger(modulus, target);

                // the state is |x+a-N⟩ in QFT basis
                copyMostSignificantBitPhaseLE(target);

                // lessThanModulusFlag is set to 1 if x+a < N
                Controlled IncrementPhaseByInteger([lessThanModulusFlag], (modulus, target));

                // the state is |x+a (mod N)⟩ in QFT basis
                // Now let us restore the lessThanModulusFlag qubit back to zero
                Controlled (Adjoint IncrementPhaseByInteger)(controls, (increment, target));
                X(lessThanModulusFlag);
                copyMostSignificantBitPhaseLE(target);
                Controlled IncrementPhaseByInteger(controls, (increment, target));
            }
        }

        controlled adjoint auto;
    }

    /// # Summary
    /// Performs a modular multiply-and-add by integer constants on a qubit register.
    ///
    /// Implements the map
    /// $$
    /// \begin{align}
    ///     \ket{x} \ket{b} \mapsto \ket{x} \ket{b + a \cdot x \operatorname{mod} N}
    /// \end{align}
    /// $$
    /// for a given modulus $N$, constant multiplier $a$, and summand $y$.
    ///
    /// # Input
    /// ## constantMultiplier
    /// An integer $a$ to be added to each basis state label.
    /// ## modulus
    /// The modulus $N$ which addition and multiplication is taken with respect to.
    /// ## multiplier
    /// A quantum register representing an unsigned integer whose value is to
    /// be added to each basis state label of `summand`.
    /// ## summand
    /// A quantum register representing an unsigned integer to use as the target
    /// for this operation.
    ///
    /// # Remarks
    /// - For the circuit diagram and explanation see Figure 6 on [Page 7
    ///        of arXiv:quant-ph/0205095v3](https://arxiv.org/pdf/quant-ph/0205095v3.pdf#page=7)
    /// - This operation corresponds to CMULT(a)MOD(N) in
    ///   [arXiv:quant-ph/0205095v3](https://arxiv.org/pdf/quant-ph/0205095v3.pdf)
    operation MultiplyAndAddByModularInteger(constMultiplier : Int, modulus : Int, multiplier : LittleEndian, summand : LittleEndian) : Unit {
        body (...) {
            let inner = MultiplyAndAddPhaseByModularInteger(constMultiplier, modulus, multiplier, _);

            using (extraZeroBit = Qubit()) {
                ApplyPhaseLEOperationOnLECA(inner, LittleEndian(summand! + [extraZeroBit]));
            }
        }

        adjoint invert;
        controlled distribute;
        controlled adjoint distribute;
    }

    /// # Summary
    /// The same as ModularAddProductLE, but assumes that summand encodes
    /// integers in QFT basis
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ModularAddProductLE
    ///
    /// # Remarks
    /// Assumes that `phaseSummand` has the highest bit set to 0.
    /// Also assumes that the value of `phaseSummand` is less than N.
    operation MultiplyAndAddPhaseByModularInteger(constMultiplier : Int, modulus : Int, multiplier : LittleEndian, phaseSummand : PhaseLittleEndian) : Unit is Adj + Ctl {
        EqualityFactB(modulus <= 2 ^ (Length(phaseSummand!) - 1), true, $"`multiplier` must be big enough to fit integers modulo `modulus`" + $"with highest bit set to 0");
        EqualityFactB(constMultiplier >= 0 and constMultiplier < modulus, true, $"`constMultiplier` must be between 0 and `modulus`-1");

        if (_EnableExtraAssertsForArithmetic()) {
            // assert that the highest bit is zero, by switching to computational basis
            ApplyLEOperationOnPhaseLECA(AssertHighestBit(Zero, _), phaseSummand);

            // check that the input is less than modulus
            AssertPhaseLessThan(modulus, phaseSummand);
        }

        for (i in 0 .. Length(multiplier!) - 1) {
            let summand = (ExpMod(2, i, modulus) * constMultiplier) % modulus;
            Controlled IncrementPhaseByModularInteger([(multiplier!)[i]], (summand, modulus, phaseSummand));
        }
    }

    /// # Summary
    /// Performs modular multiplication by an integer constant on a qubit register.
    ///
    /// Let us denote modulus by N and constMultiplier by a
    /// then this operation implements a unitary defined by the following map on
    /// computational basis:
    /// |y⟩ ↦ |a⋅y (mod N) ⟩, for all y between 0 and N - 1
    ///
    /// # Input
    /// ## constMultiplier
    /// Constant by which multiplier is being multiplied. Must be co-prime to modulus.
    /// ## modulus
    /// The multiplication operation is performed modulo `modulus`
    /// ## multiplier
    /// The number being multiplied by a constant.
    /// This is an array of qubits representing integer in little-endian bit order.
    ///
    /// # Remarks
    /// - For the circuit diagram and explanation see Figure 7 on [Page 8
    ///        of arXiv:quant-ph/0205095v3](https://arxiv.org/pdf/quant-ph/0205095v3.pdf#page=8)
    /// - This operation corresponds to Uₐ in
    ///   [arXiv:quant-ph/0205095v3](https://arxiv.org/pdf/quant-ph/0205095v3.pdf)
    operation MultiplyByModularInteger(constMultiplier : Int, modulus : Int, multiplier : LittleEndian) : Unit is Adj + Ctl {
        // Check the preconditions using Microsoft.Quantum.Canon.EqualityFactB
        EqualityFactB(constMultiplier >= 0 and constMultiplier < modulus, true, $"`constMultiplier` must be between 0 and `modulus`");
        EqualityFactB(modulus <= 2 ^ Length(multiplier!), true, $"`multiplier` must be big enough to fit integers modulo `modulus`");
        EqualityFactB(IsCoprime(constMultiplier, modulus), true, $"`constMultiplier` and `modulus` must be co-prime");

        using (summand = Qubit[Length(multiplier!)]) {
            // recall that newly allocated qubits are all in 0 state
            // and therefore summandLE encodes 0.
            let summandLE = LittleEndian(summand);

            // Let us look at what is the result of operations below assuming
            // multiplier is in computational basis and encodes x
            // Currently the joint state of multiplier and summandLE is
            // |x⟩|0⟩
            MultiplyAndAddByModularInteger(constMultiplier, modulus, multiplier, summandLE);

            // now the joint state is |x⟩|x⋅a(mod N)⟩
            ApplyToEachCA(SWAP, Zip(summandLE!, multiplier!));

            // now the joint state is |x⋅a(mod N)⟩|x⟩
            let inverseMod = InverseMod(constMultiplier, modulus);

            // note that the operation below implements the following map:
            // |x⟩|y⟩ ↦ |x⟩|y - a⁻¹⋅x (mod N)⟩
            Adjoint MultiplyAndAddByModularInteger(inverseMod, modulus, multiplier, summandLE);
            // now the joint state is |x⋅a(mod N)⟩|x - a⁻¹⋅x⋅a (mod N)⟩ = |x⋅a(mod N)⟩|0⟩
        }
    }

}

// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Math;

    /// # Summary
    /// Syndrome measurement and the inverse of embedding.
    /// $X$- and $Z$-stabilizers are not treated equally,
    /// which is due to the particular choice of the encoding circuit.
    /// This asymmetry leads to a different syndrome extraction routine.
    /// One could measure the syndrome by measuring multi-qubit Pauli operator
    /// directly on the code state, but for the distillation purpose
    /// the logical qubit is returned into a single qubit,
    /// in course of which the syndrome measurements can be done without further ancillas.
    ///
    /// # Output
    /// The logical qubit and a pair of integers for $X$-syndrome and $Z$-syndrome.
    /// They represent the index of the code qubit on which a single $X$- or $Z$-error
    /// would have caused the measured syndrome.
    ///
    /// # Remarks
    /// > [!WARNING]
    /// > This routine is tailored 
    /// > to a particular encoding circuit for Steane's 7 qubit code;
    /// > if the encoding circuit is modified then the syndrome outcome
    /// > might have to be interpreted differently.
    operation _ExtractLogicalQubitFromSteaneCode(code: LogicalRegister) : (Qubit, Int, Int)
    {
        body {
            (Adjoint SteaneCodeEncoderImpl)(code[0..0], code[1..6]);

            let x0 = M( code[6] );
            let x1 = M( code[1] );
            let x2 = M( code[3] );

            mutable xsyn = 0;
            if( x0 == One ) { set xsyn = xsyn ^^^ 1; }
            if( x1 == One ) { set xsyn = xsyn ^^^ 2; }
            if( x2 == One ) { set xsyn = xsyn ^^^ 4; }
            set xsyn = xsyn - 1;
            // xsyn contains the qubit index (0..6) at which a single Z-error would
            // produce the given syndrome.

            let z0 = M(code[5]);
            let z1 = M(code[2]);
            let z2 = M(code[4]);

            mutable zsyn = 0;
            if( z0 == One ) { set zsyn = zsyn ^^^ 1; }
            if( z1 == One ) { set zsyn = zsyn ^^^ 2; }
            if( z2 == One ) { set zsyn = zsyn ^^^ 5; }
            set zsyn = zsyn - 1;
            // zsyn contains the qubit index (0..6) at which a single X-error would
            // produce the given syndrome.

            return (code[0], xsyn, zsyn);
        }
    }

    /// # Summary
    /// InjectPi4YRotation needs a magic state 
    /// $\cos\frac{\pi}{8} \ket{0} + \sin \frac{\pi}{8} \ket{1}$.
    /// The magic state required is the same for the Adjoint,
    /// that implements $-\pi/4$ $Y$-rotation.
    operation InjectPi4YRotation(data: Qubit, magic: Qubit) : ()
    {
        body {
            (Adjoint S)(data);
            CNOT(magic, data);
            S(data);
            let r = Measure([PauliY], [magic]);
            if ( r == One ) {
                // The following five gates is equal to	Ry( Pi()/2.0, data)
                // up to global phase.
                S(data);
                H(data);
                (Adjoint S)(data);
                H(data);
                (Adjoint S)(data);
            }
        }
        adjoint {
            (Adjoint S)(data);
            CNOT(magic, data);
            S(data);
            let r = Measure([PauliY], [magic]);
            if ( r == Zero ) {
                S(data);
                H(data);
                S(data);
                H(data);
                (Adjoint S)(data);
            }
        }
    }

    operation Pi4YInjectionTest():()
    {
        body {
            using (anc = Qubit[2]) {
                // magic state in anc[1]
                Ry( PI() / 4.0, anc[1]);
                InjectPi4YRotation( anc[0], anc[1] );

                // For test, we bring the data qubit anc[0] to the original state.
                Ry(-PI() / 4.0, anc[0]);

                // So, the data must be in Zero.
                AssertQubit(Zero, anc[0] );
                ResetAll(anc);
            }
            using (anc2 = Qubit[2]) {
                // magic state in anc[1]
                Ry( PI() / 4.0, anc2[1]);
                (Adjoint InjectPi4YRotation)( anc2[0], anc2[1] );

                // For test, we bring the data qubit anc[0] to the original state.
                Ry( PI() / 4.0, anc2[0]);

                // So, the data must be in Zero.
                AssertQubit(Zero, anc2[0] );
                ResetAll(anc2);
            }
        }
    }

    // TODO: define magic states here.
    // TODO: define what happens to the other fourteen qubits.
    /// # Summary
    /// Given 15 approximate copies of a magic state, yields
    /// one higher-quality copy.
    ///
    /// # Input
    /// ## roughMagic
    /// A register of fifteen qubits containing approximate copies
    /// of a magic state. Following application of this distillation
    /// procedure, ``roughMagic[0]` will contain one higher-quality
    /// copy.
    ///
    /// # Output
    /// If `true`, then the procedure succeeded and the higher-quality
    /// copy should be accepted. If `false`, the procedure failed, and
    /// the state of the register should be considered undefined.
    ///
    /// # Remarks
    /// We follow the algorithm of Knill.
    /// However, the present implementation is far from being optimal,
    /// as it uses too many qubits.
    /// The magic states are injected in this routine,
    /// in which case there are better protocols.
    ///
    /// # References
    /// - [Knill](https://arxiv.org/abs/quant-ph/0402171)
    operation KnillDistill( roughMagic : Qubit[] ) : ( Bool )
    {
        body {
            mutable accept = false;
            using (scratch = Qubit[8]) {
                let anc = scratch[7];
                let code = scratch[0..6];
                InjectPi4YRotation( code[0], roughMagic[14] );
                SteaneCodeEncoderImpl( code[0..0], code[1..6] );
                for ( idx in 0..6 ) {
                    (Adjoint InjectPi4YRotation)( code[idx], roughMagic[idx] );
                    CNOT(code[idx], anc);
                    InjectPi4YRotation(code[idx], roughMagic[idx + 7]);
                }
                let (logicalQubit, xsyn, zsyn) = 
                    _ExtractLogicalQubitFromSteaneCode(LogicalRegister(code));
                let m = M(anc);
                if( xsyn == -1 && zsyn == -1 && m == Zero ) {
                    SWAP(logicalQubit, roughMagic[0]);
                    set accept = true;
                }

                ResetAll(scratch);
            }
            return accept;
        }
    }

    
}

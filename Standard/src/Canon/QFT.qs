// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Apply the Approximate Quantum Fourier Transform (AQFT) to a quantum register.
    ///
    /// # Input
    /// ## a
    /// approximation parameter which determines at which level the controlled Z-rotations that
    /// occur in the QFT circuit are pruned.
    ///
    /// The approximation parameter a determines the pruning level of the Z-rotations, i.e.,
    /// a ∈ {0..n} and all Z-rotations 2π/2ᵏ where k>a are
    /// removed from the QFT circuit. It is known that for k >= log₂(n)+log₂(1/ε)+3
    /// one can bound ||QFT-AQFT||<ε.
    ///
    /// ## qs
    /// quantum register of n qubits to which the Approximate Quantum Fourier Transform is applied.
    ///
    /// # Remarks
    /// AQFT requires Z-rotation gates of the form 2π/2ᵏ and Hadamard gates.
    ///
    /// The input and output are assumed to be encoded in big endian encoding.
    ///
    ///
    /// # References
    /// - [ *M. Roetteler, Th. Beth*,
    ///      Appl. Algebra Eng. Commun. Comput.
    ///      19(3): 177-193 (2008) ](http://doi.org/10.1007/s00200-008-0072-2)
    /// - [ *D. Coppersmith* arXiv:quant-ph/0201067v1 ](https://arxiv.org/abs/quant-ph/0201067)
    operation ApproximateQFT (a : Int, qs : BigEndian) : Unit is Adj + Ctl {
        let nQubits = Length(qs!);
        Fact(nQubits > 0, $"`Length(qs)` must be least 1");
        Fact(a > 0 and a <= nQubits, $"`a` must be positive and less than `Length(qs)`");

        for (i in 0 .. nQubits - 1) {
            for (j in 0 .. i - 1) {
                if (i - j < a) {
                    Controlled R1Frac([(qs!)[i]], (1, i - j, (qs!)[j]));
                }
            }

            H((qs!)[i]);
        }

        // Apply the bit reversal permutation to the quantum register as
        // a side effect, such that we enforce the invariants specified
        // by the BigEndian UDT.
        SwapReverseRegister(qs!);
    }

    /// # Summary
    /// Performs the Quantum Fourier Transform on a quantum register containing an
    /// integer in the big-endian representation.
    ///
    /// # Input
    /// ## qs
    /// Quantum register to which the Quantum Fourier Transform is applied
    ///
    /// # Remarks
    /// The input and output are assumed to be in big endian encoding.
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApproximateQFT
    /// - Microsoft.Quantum.Canon.QFTLE
    internal operation ApplyQuantumFourierTransformBE(qs : BigEndian) : Unit is Adj + Ctl {
        ApproximateQFT(Length(qs!), qs);
    }

    /// # Summary
    /// Performs the Quantum Fourier Transform on a quantum register containing an
    /// integer in the little-endian representation.
    ///
    /// # Input
    /// ## qs
    /// Quantum register to which the Quantum Fourier Transform is applied
    ///
    /// # Remarks
    /// The input and output are assumed to be in little endian encoding.
    ///
    ///
    /// # See Also
    /// - Microsoft.Quantum.Canon.ApplyQuantumFourierTransformBE
    operation ApplyQuantumFourierTransform(qs : LittleEndian) : Unit is Adj + Ctl {
        ApplyReversedOpBECA(ApplyQuantumFourierTransformBE, qs);
    }

}

// Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the 
// Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries 
// and Samples. See LICENSE in the project root for license information.

namespace Microsoft.Quantum.Canon {
    // Including the namespace Primitive gives access to basic operations such as the 
	// Hadamard gates, CNOT gates, etc. that are useful for defining circuits. The 
	// implementation of these operations is dependent on the targeted machine. 
	open Microsoft.Quantum.Primitive;
	// The canon namespace contains many useful library functions for creating 
	// larger circuits, combinators, and utility functions. The implementation of 
	// the operations in the canon is machine independent as they are built on 
	// top of the primitive operations. 
	open Microsoft.Quantum.Canon;

	//////////////////////////////////////////////////////////////////////////
	// Introduction //////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	/// This sample contains serveral simple quantum algorithms coded in Q#. The 
	/// intent is to highlight the expressive capabilities of the language that 
	/// allow to express quantum algirthm that consist of a short quantum part and 
	/// classical post-processing that is simple, or in some cases, trivial.

	// First, note that every Q# function needs to have a namespace. We define 
	// a new one for this purpose. 

    //////////////////////////////////////////////////////////////////////////
    // Bernstein-Vazirani Fouier Sampling Quantum Algorithm //////////////////
    //////////////////////////////////////////////////////////////////////////

	/// # Summary 
	/// ParityViaFourierSampling implements the Bernstein-Vazirani quantum algorith. 
	/// This algorithm computes for a given Boolean function that is promised to be
	/// a parity $f(x_0, \ldots, x_{n-1}) = \sum_i r_i x_i$ a result in form of  
	/// a bit vector $(r_0, \ldots, r_{n-1})$ corresponding to the parity function. 
	/// Note that it is promised that the function is actually a parity function. 
	///
	/// # Input
	/// ## Uf
	/// A quantum circuit that implements $\ket{x}\ket{y}\mapsto\ket{x}\ket{y\oplus f(x)}$, 
	/// where $f$ is a Boolean function that implements a parity $\sum_i r_i x_i$. 
	/// ## n 
	/// The number of bits of the input register $\ket{x}$.
	///
	/// # Output
	/// An array of type `Bool[]` that contains the parity $(r_0, \ldots, r_{n-1})$. 
	///
	/// # See Also
    /// - For details see Section 1.4.3 of Nielsen & Chuang
    ///
    /// # References
    /// - [ *Ethan Bernstein and Umesh Vazirani*, 
	///     SIAM J. Comput., 26(5), 1411–1473, 1997 ](https://doi.org/10.1137/S0097539796300921)
	operation ParityViaFourierSampling(Uf : (Qubit[] => ()), n : Int) : Bool[] { 
		body {
			// we first create an array of size n which will hold the final result.
			mutable resultArray = new Result[n];
			// now, we allocate n+1 clean qubits. Note that the function Uf is defined
			// on inputs of the form (x, y), where x has n bits and y has 1 bit.
			using(qubits=Qubit[n+1]) {				
				// the last qubit needs to be flipped so that the function will 
				// actually be computed into the phase when Uf is applied. 
				X(qubits[n]);
				// now, a Hadamard transform is applied to each of the qubits.
				ApplyToEach(H, qubits);
				// we now apply Uf to the n+1 qubits, computing |x,y〉 -> |x, y+f(x)〉
				Uf(qubits);
				// as the last step before the measurement, a Hadamard transform is 
				// but the very last one. We could apply the Hadamard transform to 
				// the last qubit also, but this would not affect the final outcome. 
				ApplyToEach(H, qubits[0..(n-1)]); 
				// the following for-loop measures all qubits and resets them to 
				// zero so that they can be safely returned at the end of the 
				// using-block.
				for (idx in 0..(n-1)) {
					set resultArray[idx] = MResetZ(qubits[idx]);
				}
				// finally, the last qubit, which held the y-register, is reset. 
				Reset(qubits[n]);							
			}	
			// the result is already contained in resultArray and not further 
			// post-processing is necessary.
			Message($"measured: {resultArray}");
			return BoolArrFromResultArr(resultArray);
		 }
    }

	//////////////////////////////////////////////////////////////////////////
    // Deusch-Jozsa Quantum Algorithm ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

	/// # Summary 
	/// Deutsch-Jozsa is a quantum algorithm that decides whether a given Boolean function 
	/// $f$ that is promised to either be constant or to be balanced---i.e., taking the 
	/// values 0 and 1 the exact same number of times---is actually constant or balanced. 
	/// The operation `IsConstantBooleanFunction` answers this question by returning the 
	/// Boolean value `true` if the function is constant and `false` if it is not. Note 
	/// that the promise that the function is either constant or balanced is assumed.  
	///
	/// # Input
	/// ## Uf
	/// A quantum circuit that implements $\ket{x}\ket{y}\mapsto\ket{x}\ket{y\oplus f(x)}$, 
	/// where $f$ is a Boolean function, $x$ is an $n$ bit register and $y$ is a single qubit. 
	/// ## n 
	/// The number of bits of the input register $\ket{x}$.
	///
	/// # Output
	/// A boolean value `true` that indicates that the function is constant and `false` 
	/// that indicates that the function is balanced.
    ///
	/// # See Also
    /// - For details see Section 1.4.3 of Nielsen & Chuang
    ///
    /// # References
    /// - [ *Michael A. Nielsen , Isaac L. Chuang*,
    ///     Quantum Computation and Quantum Information ](http://doi.org/10.1017/CBO9780511976667)
	operation IsConstantBooleanFunction(Uf : (Qubit[] => ()), n : Int) : Bool { 
		body {
			// we first create an array of size n from which we compute the final result. 
			mutable resultArray = new Result[n];
			// now, we allocate n+1 clean qubits. Note that the function Uf is defined
			// on inputs of the form (x, y), where x has n bits and y has 1 bit.
			using(qubits=Qubit[n+1]) {				
				// the last qubit needs to be flipped so that the function will 
				// actually be computed into the phase when Uf is applied. 
				X(qubits[n]);
				// now, a Hadamard transform is applied to each of the qubits.
				ApplyToEach(H, qubits);
				// we now apply Uf to the n+1 qubits, computing |x,y〉 -> |x, y+f(x)〉
				Uf(qubits);
				// as the last step before the measurement, a Hadamard transform is 
				// but the very last one. We could apply the Hadamard transform to 
				// the last qubit also, but this would not affect the final outcome. 
				ApplyToEach(H, qubits[0..(n-1)]); 
				// the following for-loop measures all qubits and resets them to 
				// zero so that they can be safely returned at the end of the 
				// using-block.
				for (idx in 0..(n-1)) {
					set resultArray[idx] = MResetZ(qubits[idx]);
				}
				// finally, the last qubit, which held the y-register, is reset. 
				Reset(qubits[n]);							
			}	
			// we use the predicte `IsResultZero` from Microsoft.Quantum.Canon.Utils
			// (Predicates.qs) and compose it with the ForAll function from 
			// Microsoft.Quantum.Canon.Enumeration (ForAll.qs). This will return 
			// `true` if the all zero string has been measured, i.e., if the function 
			// was a constant function and `false` if not, which according to the 
			// promise on f means that it must have been balanced. 
			return ForAll(IsResultZero, resultArray);
		 }
    }
	
	///////////////////////////////////////////////////////////////////////////
    // Hidden Shift Quantum Algorithm /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

	/// # Summary 
	/// Correlation-based algorithm to solve the hidden shift problem for bent functions. 
	/// The problem is to identify an unknown shift $s$ of the arguments of two Boolean functions 
	/// $f$ and $g$ that are promised to satisfy the relation $g(x) = f(x \oplus s)$ for all $x$. 
	/// Note that the promise about the functions $f$ and $g$ to be bent functions is assumed, 
	/// i.e., they both have a flat Fourier (Walsh-Hadamard) spectra. Input to this algorithm 
	/// are implementations $U_g$ of the Boolean function $g$ and $U_f^*$, an implementation of 
	/// dual bent function of the function $f$. Both functions are given via phase encoding.
	/// 
	/// # Input
	/// ## Ufstar
	/// A quantum circuit that implements $U_f^*:\ket{x}\mapsto (-1)^{f^*(x)} \ket{x}$, 
	/// where $f^*$ is a Boolean function, $x$ is an $n$ bit register and $y$ is a single qubit. 
	/// ## Ug 
	/// A quantum circuit that implements $U_g:\ket{x}\mapsto (-1)^{g(x)} \ket{x}$, 
	/// where $g$ is a Boolean function that is shifted by unknown $s$ from $f$, and $x$ is 
	/// an $n$ bit register.
	/// ## n 
	/// The number of bits of the input register $\ket{x}$.
	///
	/// # Output
	/// An array of type `Bool[]` which encodes the bit representation of the hidden shift.
	/// 
	/// # References
	/// - [*Martin Roetteler*, 
	///    Proc. SODA 2010, ACM, pp. 448-457, 2010](https://doi.org/10.1137/1.9781611973075.37)
	operation HiddenShiftBentCorrelation (Ufstar : (Qubit[] => ()), Ug : (Qubit[] => ()), n : Int) : Bool[] { 
		body {
			// we first create an array of size n from which we compute the final result. 
			mutable resultArray = new Result[n];
			// now, we allocate n clean qubits. Note that the function Ufstar and Ug are 
			// unitary operations on n qubits defined via phase encoding.
			using(qubits=Qubit[n]) {				
				// first, a Hadamard transform is applied to each of the qubits.
				ApplyToEach(H, qubits);
				// we now apply the shifted function Ug to the n qubits, computing 
				// |x〉 -> (-1)^{g(x)} |x〉. 
				Ug(qubits);
				// now, a Hadamard transform is applied to each of the n qubits.
				ApplyToEach(H, qubits);
				// we now apply the dual function of the unshifted function, i.e., Ufstar, 
				// to the n qubits, computing |x〉 -> (-1)^{fstar(x)} |x〉.
				Ufstar(qubits);
				// now, a Hadamard transform is applied to each of the n qubits.
				ApplyToEach(H, qubits);
				// the following for-loop measures the n qubits and resets them to 
				// zero so that they can be safely returned at the end of the 
				// using-block.
				for (idx in 0..(n-1)) {
					set resultArray[idx] = MResetZ(qubits[idx]);
				}
			}	
			// the result is already contained in resultArray and not further 
			// post-processing is necessary except for a conversion from Result[] to 
			// Bool[] for which we use a canon function (from TypeConversion.qs).
			Message($"measured: {resultArray}");
			return BoolArrFromResultArr(resultArray);
		 }
    }

}


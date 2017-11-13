// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Samples.Ising {
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    //////////////////////////////////////////////////////////////////////////
    // Introduction //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // In this sample, we demonstrate using the generator representation
    // functionality offered by the Q# canon to represent Ising model
    // Hamiltonians.
    // In later samples, we will extend these techniques to represent the
    // 1D Heisenberg XXZ model.
    
    // We will begin by constructing a representation of the 1D transverse
    // Ising model Hamiltonian,
    //     H = - J0 Z0 Z1 + J1 Z1 Z2 + … - h? (X0 + X1 + …),
    // where {J?} are nearest-neighbor couplings, and where h? is a
    // transverse field.

    // Since this Hamiltonian is naturally expressed in the Pauli basis,
    // we will use the PauliEvolutionSet() function to obtain a simulatable
    // basis to use in representing H. Thus, we begin by defining our
    // indices with respect to the Pauli basis. In doing so, we will define
    // helper functions to return single-site and two-site generator indices.


    // FIXME: rename to make the role of this clear and to expand abbreviations.
    // FIXME: rename hC to be more clear.
    // TODO: consider generalizing (that is, remove idxPauli notation) and promoting
    //       back to Canon as a useful utility function.
    /// # Summary
    /// Returns a generator index that is supported on a single site.
    ///
    /// # Input
    /// ## idxPauli
    /// Index of the Pauli operator to be represented, where `1` denotes `PauliX`, `2` denotes
    /// `PauliY` and `3` denotes `PauliZ`.
    /// ## idxQubit
    /// Index of the qubit that represented term will act upon.
    /// ## hC
    /// Function returning coefficients for each site. E.g.: `hC(3)` should return the coefficient
    /// for the index at `idxQubit = 3`.
    ///
    /// # Output
    /// A `GeneratorIndex` representing the term $h(i) \sigma_{\mu}^{(i)}$, where $h(i)$ is the
    /// function `hC` evaluated at the site index `i`, and where
    /// $\sigma_{\mu}^{(i)}$ is the $\mu$th Pauli operator acting at the site
    /// index `i`.
    function OneSiteGenIdx(idxPauli: Int, idxQubit: Int, hC: (Int -> Double)) : GeneratorIndex {
        let coeff = - 1.0 * hC(idxQubit);
        let idxPauliString = [idxPauli];
        let idxQubits = [idxQubit];
        return GeneratorIndex((idxPauliString, [coeff]), idxQubits);
    }
    // TODO: document as the above.
    function OneSiteGenSys(idxPauli: Int, nQubits: Int, hC: (Int -> Double)) : GeneratorSystem {
        return GeneratorSystem(nQubits, OneSiteGenIdx(idxPauli, _, hC));
    }
   
    /// We define the coupling Hamiltonian - (J0 Z0 Z1 + J1 Z1 Z2 + ...)
    function TwoSiteGenIdx(idxPauli: Int, nSites: Int, idxQubit : Int , jC: (Int -> Double)) : GeneratorIndex
    {
        /// when idxQubit is in [0, nSites - 1], apply Ising couplings jC(idxQubit)
        let idx = idxQubit;
        let coeff = - 1.0 * jC(idx);
        let idxPauliString = [idxPauli; idxPauli];
        let idxQubits = [idx; (idx + 1) % nSites];
        let generatorIndex = GeneratorIndex((idxPauliString,[coeff]),idxQubits);
        if (idx >= nSites) {
            fail "Qubit index must be smaller than number of sites.";
        }
        return generatorIndex;
    }
    /// We define the coupling Hamiltonian - (J0 Z0 Z1 + J1 Z1 Z2 + ...)
    function TwoSiteGenSys(idxPauli: Int, nSites: Int, jC: (Int -> Double)) : GeneratorSystem
    {
        return GeneratorSystem(nSites,  TwoSiteGenIdx(idxPauli, nSites, _, jC));
    }

    /// We now add the transverse and coupling Hamiltonians
    function IsingGenSys(nSites: Int, hX: (Int -> Double), jC: (Int -> Double)) : GeneratorSystem {
        let XGenSys = OneSiteGenSys(1, nSites, hX);
        let ZZGenSys = TwoSiteGenSys(3, nSites, jC);
        return AddGeneratorSystems(XGenSys, ZZGenSys);
    }

    /// This provides a full description of the Hamiltonian -- a map from a 1D index to terms in the Hamiltonian, and then to unitaries
    function Ising1DEvolutionGenerator(nSites : Int, hX: (Int -> Double),  jC: (Int -> Double)) : EvolutionGenerator {
        let generatorSystem = IsingGenSys(nSites, hX, jC);
        let evolutionSet = PauliEvolutionSet();
        return EvolutionGenerator(evolutionSet, generatorSystem);
    }

    /// We now define functions for the coefficients
    /// Uniform transvers field couplings
    function GenerateUniformHCoupling( amplitude : Double, idxQubit : Int) : Double
    {
        return 1.0 * amplitude;
    }

    /// Uniform coupling with open boundary conditions
    /// Q - ... - Q - Q - ... - Q
    function GenerateUniform1DJCoupling(nSites: Int, amplitude: Double, idxQubit: Int): Double {
        mutable coeff = amplitude;
        if(idxQubit == nSites - 1){
            set coeff = ToDouble(0);
        }
        return coeff;
    }

    /// It is straightforward to generalize this to the Heisenberg XXZ model
    /// Hamiltonian = -( JXY0 (X0 X1  + Y0 Y1) JXY1 (X1 X2 + ...) + ... ) - hz (X0 + X1 + ...)
    function HeisenbergXXZGenSys(nSites: Int, hZ: (Int -> Double), jC: (Int -> Double)) : GeneratorSystem {
        let ZGenSys = OneSiteGenSys(3, nSites, hZ);
        let XXGenSys = TwoSiteGenSys(3, nSites, jC);
        let YYGenSys = TwoSiteGenSys(3, nSites, jC);
        /// This multiplies all coefficients in a generator system
        let jZZmultiplier = 0.5;
        let ZZGenSys = MultiplyGeneratorSystem( jZZmultiplier, TwoSiteGenSys(3, nSites, jC) );
        return AddGeneratorSystems(AddGeneratorSystems(ZGenSys, ZZGenSys),AddGeneratorSystems(YYGenSys, XXGenSys));
    }    /// We now add the transverse and coupling Hamiltonians

}

// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Xunit;
using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;

using System.Text.RegularExpressions;
using System;
using System.Linq;
using System.Collections.Generic;

namespace Microsoft.Quantum.Chemistry.Tests
{
    //using static FermionTermType.Common;
    using FermionTerm = FermionTerm;
    //using FermionTermType = FermionTermType;
    using SpinOrbital = SpinOrbital;
    using OrbitalIntegral = OrbitalIntegral;
    using Spin = Spin;
    /*
public class FermionHamiltonianTests
{

   public FermionHamiltonian GenerateTestHamiltonian()
    {
        Int64 nOrbitals = 6;
        Dictionary<FermionTermType, List<FermionTerm>> fermionTerms = new Dictionary<FermionTermType, List<FermionTerm>>
        {
            { PPTermType, new List<FermionTerm>() },
            { PQTermType, new List<FermionTerm>() },
            { PQQPTermType, new List<FermionTerm>() },
            { PQQRTermType, new List<FermionTerm>() },
            { PQRSTermType, new List<FermionTerm>() }
        };

        fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 0 }, 1.0));
        fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 1 }, 1.0));
        fermionTerms[PPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 2 }, 1.0));

        fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 0, 2 }, 1.0));
        fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 1, 3 }, 1.0));
        fermionTerms[PQTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 0 }, new Int64[] { 2, 6 }, 1.0));

        fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2,2 ,0 }, 1.0));
        fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3,1 }, 1.0));
        fermionTerms[PQQPTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 2 }, 1.0));

        fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 2, 1 }, 1.0));
        fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 3, 3, 2 }, 1.0));
        fermionTerms[PQQRTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 6, 6, 5 }, 1.0));

        fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 0, 2, 4, 3 }, 1.0));
        fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 1, 4, 3, 2 }, 1.0));
        fermionTerms[PQRSTermType].Add(new FermionTerm(nOrbitals, new Int64[] { 1, 1, 0, 0 }, new Int64[] { 2, 4,5, 3 }, 1.0));
        return new FermionHamiltonian(fermionTerms, nOrbitals);
    }

    [Fact]
    public void VerifyFermionTermsTest()
    {
        var ham = GenerateTestHamiltonian();
        Assert.True(ham.VerifyFermionTerms());
    }


}

public class FermionTermTypeTests
{

    [Theory]
    [InlineData(true, 1L, new Int64[] { 0, 0, 0, 0 })]
    [InlineData(true, 1L, new Int64[] { 1, 1, 0, 0 })]
    [InlineData(true, 2L, new Int64[] { 1, 0, })]
    [InlineData(true, 2L, new Int64[] { 1, 1, })]
    [InlineData(true, 0L, new Int64[] { })]
    [InlineData(false, 1L, new Int64[] { 0, 0, 1, 0 })]
    [InlineData(false, 1L, new Int64[] { 1, 0, 1, 0 })]
    [InlineData(false, 3L, new Int64[] { 1, 0, })]
    [InlineData(false, 1L, new Int64[] { })]
    public void IsInCanonicalOrderTest(bool pass, Int64 uniqueTerms, Int64[] type)
    {
        var tmp = new FermionTermType(uniqueTerms, type);
        if (pass)
        {
            Assert.True(tmp.IsInCanonicalOrder());
        }
        else
        {
            Assert.False(tmp.IsInCanonicalOrder());
        }
    }


    [Fact]
    public void IsInCanonicalOrderCommonTypesTest()
    {
        var tmp = new FermionTermType[] { IdentityTermType,
        PPTermType,
        PQTermType,
        PQQPTermType,
        PQQRTermType,
        PQRSTermType};
        foreach(var item in tmp){
            Assert.True(item.IsInCanonicalOrder());            
        }
    }
}
*/

}
// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.Quantum.Chemistry;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.Simulators.QCTraceSimulators;
using Xunit;

namespace SystemTestsOptimizedBlockEncoding
{

    /*
    using static FermionTermType.Common;

    public class TraceConfig
    {
        public static QCTraceSimulatorConfiguration config = new QCTraceSimulatorConfiguration()
        {
            usePrimitiveOperationsCounter = true,
            throwOnUnconstraintMeasurement = false
        };
    }

    public class SystemTestFromGeneralHamiltonian
    {
        public Double targetError = 0.5;

        [Fact]
        public void PPTermFromGeneralHamiltonianTest()
        {
            
            var generalHamiltonian = new FermionHamiltonian(nOrbitals: 1, nElectrons: 1);
            generalHamiltonian.AddFermionTerm(PPTermType, new Int64[] { 0, 0 }, 1.0);
            generalHamiltonian.AddFermionTerm(PPTermType, new Int64[] { 1, 1 }, -1.0);
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }



        [Fact]
        public void PQTermABCFromGeneralHamiltonianTest()
        {
            
            var generalHamiltonian = new FermionHamiltonian(nOrbitals: 3, nElectrons: 1);
            generalHamiltonian.AddFermionTerm(PQTermType, new Int64[] { 0, 2 }, 2.0);
            generalHamiltonian.AddFermionTerm(PQTermType, new Int64[] { 0, 1 }, -1.0);
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQQPTermFromGeneralHamiltonianTest()
        {
            
            var generalHamiltonian = new FermionHamiltonian(nOrbitals: 3, nElectrons: 1);
            generalHamiltonian.AddFermionTerm(PQQPTermType, new Int64[] { 0, 1, 1, 0 }, 1.0);
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQQRTermFromGeneralHamiltonianTest()
        {
            
            var generalHamiltonian = new FermionHamiltonian(nOrbitals: 2, nElectrons: 1);
            generalHamiltonian.AddFermionTerm(PQQRTermType, new Int64[] { 0, 1, 2, 0 }, 1.0);
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQRSTermFromGeneralHamiltonianTest()
        {
            
            var generalHamiltonian = new FermionHamiltonian(nOrbitals: 2, nElectrons: 1);
            generalHamiltonian.AddFermionTerm(PQRSTermType, new Int64[] { 0, 1, 3, 2 }, 2.0);
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }
    }

    public class SystemTestFromLiquidOrbital
    {
        public Double targetError = 0.5;

        [Fact]
        public void PPTermFromLiquidOrbitalTest()
        {
            
            string orbitals = "0,0=1.0";
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQTermABFromLiquidOrbitalTest()
        {
            
            string orbitals = "0,1=1.0";
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQTermACFromLiquidOrbitalTest()
        {
            
            string orbitals = "0,2=1.0";
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQQPTermFromLiquidOrbitalTest()
        {
            
            string orbitals = "0,1,1,0=1.0";
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var identityCoefficient = jwEvolutionSetData.energyOffset;
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQQRTermFromLiquidOrbitalTest()
        {
            
            string orbitals =  "0,0,1,0=1.0" ;
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

        [Fact]
        public void PQRSTermFromLiquidOrbitalTest()
        {
            
            string orbitals = "0,1,0,1=1.0" ;
            var generalHamiltonian = LoadData.LoadFromLiquid(orbitals);
            generalHamiltonian.NElectrons = 1L;
            var jwEvolutionSetData = JordanWignerEncoding.Create(generalHamiltonian);
            var termData = jwEvolutionSetData.Terms;
            var qsim = new QCTraceSimulator(TraceConfig.config);
            RunOptimizedBlockEncoding.Run(qsim, generalHamiltonian.NOrbitals * 2, termData, targetError).Wait();
        }

    
    }*/
}

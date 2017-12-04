// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Simulation.Core;

namespace IntegerFactorizationWithCanon
{
    class Program
    {
        static void Main(string[] args)
        {
            long numberToFactor = 15;
            long nTrials = 100;
            bool useRobustPhaseEstimation = true;

            // Parse the arguments provided in command line
            if( args.Length >= 1 )
            {
                // first argument is the number to factor
                Int64.TryParse(args[0], out numberToFactor);
            }

            if (args.Length >= 2 )
            {
                // the second is the number of trials 
                Int64.TryParse(args[1], out nTrials);
            }

            if (args.Length >= 3)
            {
                // the third on indicates if Robust or Quantum Phase Estimation 
                // should be used
                bool.TryParse(args[2], out useRobustPhaseEstimation);
            }

            for (int i = 0; i < nTrials; ++i)
            {
                try
                {
                    using (QuantumSimulator sim = new QuantumSimulator())
                    {
                        Console.WriteLine($"==========================================");
                        Console.WriteLine($"Factoring {numberToFactor}");
                        (long factor1, long factor2) = 
                            ShorWithCanon.Run(
                                sim, numberToFactor, useRobustPhaseEstimation).Result;

                        Console.WriteLine($"Factors are {factor1} and {factor2}");
                    }
                }
                catch(AggregateException e )
                {
                    Console.WriteLine($"This run of Shor's algorithm failed:");
                    foreach ( Exception eInner in e.InnerExceptions )
                    {
                        ExecutionFailException failException = eInner as ExecutionFailException;
                        Console.WriteLine($"   {failException.Message}");
                    }
                }
            }
        }
    }
}
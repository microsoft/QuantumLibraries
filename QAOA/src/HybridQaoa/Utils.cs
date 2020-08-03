using Microsoft.Quantum.Simulation.Core;
using System;
using System.Collections.Generic;
using System.Text;

namespace Quantum.QAOA
{
    public class Utils
    {

        public struct FreeParamsVector
        {
            public Double[] beta;
            public Double[] gamma;
        }

        /// # Summary
        /// Returns a vector of random doubles in a range from 0 to maximum.
        ///
        /// # Input
        /// ## length
        /// Length of a random vector.
        /// ## maximum
        /// Maximum value of a random double.
        /// 
        /// # Output
        /// A random vector of doubles.
        
        public static double[] GetRandomVector(int length, double maximumValue)
        {
            var rand = new Random();
            double[] randomVector = new double[length];
            for (int i = 0; i < length; i++)
            {
                randomVector[i] = maximumValue * rand.NextDouble();
            }

            return randomVector;
        }

        /// # Summary
        /// Return the most common boolean string from a list of boolean values.
        ///
        /// # Input
        /// ## list
        /// List of boolean values.
        /// 
        /// # Output
        /// The most common boolean string.

        public static String GetModeFromBoolList(List<bool[]> list)
        {
            Dictionary<string, int> counter = new Dictionary<string, int>();
            foreach (bool[] boolArray in list)
            {
                String boolString = GetBoolStringFromBoolArray(boolArray);
                if (counter.ContainsKey(boolString))
                {
                    counter[boolString] += 1;
                }
                else
                {
                    counter[boolString] = 1;
                }

            }
            int maximum = 0;
            String result = null;
            foreach (string key in counter.Keys)
            {
                if (counter[key] > maximum)
                {
                    maximum = counter[key];
                    result = key;
                }
            }

            return result;
        }

        /// # Summary
        /// Converts an array of bools to a boolean string.
        ///
        /// # Input
        /// ## boolArray
        /// An array of bools.
        /// 
        /// # Output
        /// A boolean string.

        public static string GetBoolStringFromBoolArray(bool[] boolArray)
        {
            System.Text.StringBuilder sb = new StringBuilder();
            foreach (bool b in boolArray)
            {
                sb.Append(b ? "1" : "0");
            }
            return sb.ToString();
        }

        /// # Summary
        /// Converts concatenated beta and gamma vectors into separate beta and gamma vector.
        ///
        /// # Input
        /// ## bigfreeParamsVector
        /// Concatenated beta and gamma vectors.
        ///
        /// # Output
        /// FreeParamsVector that contains beta and gamma vectors.
        ///
        /// # Remarks
        /// Useful for getting beta and gamma vectors from a concatenated vector inside the optimized function.

        public static FreeParamsVector ConvertVectorIntoHalves(double[] bigfreeParamsVector)
        {
            int size = bigfreeParamsVector.Length;
            int vectorTermsNumber = size / 2;
            FreeParamsVector freeParamsVector = new FreeParamsVector
            {

                beta = bigfreeParamsVector[0..vectorTermsNumber],
                gamma = bigfreeParamsVector[vectorTermsNumber..(2*vectorTermsNumber)],

            };

            return freeParamsVector;
        }

        /// # Summary
        /// Prints current values of beta and gamma vectors as optimization is being performed.
        ///
        /// # Input
        /// ## beta
        /// Beta vector of coefficients.
        /// ## gamma
        /// Gamma vector of coefficients.
        public static void PrintCurrentBetaGamma(QArray<Double> beta, QArray<Double> gamma)
        {
            Console.WriteLine("Current beta vector:");
            Console.WriteLine(beta);

            Console.WriteLine("Current gamma vector:");
            Console.WriteLine(gamma);

        }

        /// # Summary
        /// Prints current values of the best fidelity and the best solution vector as optimization is being performed.
        ///
        /// # Input
        /// ## bestHamiltonian
        /// Best value of a Hamiltonian so far.
        /// ## bestVector
        /// Best solution vector that generates the above value of a Hamiltonian  so far.
        public static void PrintCurrentBestSolution(double bestHamiltonian, String bestVector)
        {
            Console.WriteLine("Current best fidelity");
            Console.WriteLine(bestHamiltonian);
            Console.WriteLine("Current best string");
            Console.WriteLine(bestVector);
        }

        /// # Summary
        /// Prints whether an optimization finished successfully.
        ///
        /// # Input
        /// ## success
        /// A flag that indiciates whether an optimization finished successfully.
        public static void PrintSuccess(bool success)
        {
            Console.WriteLine("Was optimization successful?");
            Console.WriteLine(success);
            Console.WriteLine("##################################");
        }

    }
}
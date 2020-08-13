// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;

    public class Utils
    {

        public struct FreeParameters
        {
            public double[] Beta;
            public double[] Gamma;
        }

        /// <summary>
        /// Returns a vector of random doubles in a range from 0 to maximum.
        /// </summary>
        /// <param name="length">
        /// Length of a random vector.
        /// </param>
        /// <param name="maximumValue">
        /// Maximum value of a random double.
        /// </param>
        /// <returns>
        /// A random vector of doubles.
        /// </returns>
        public static double[] GetRandomVector(int length, double maximumValue)
        {
            var rand = new Random();
            double[] randomVector = new double[length];
            for (var i = 0; i < length; i++)
            {
                randomVector[i] = maximumValue * rand.NextDouble();
            }

            return randomVector;
        }

        /// <summary>
        /// Return the most common boolean string from a list of boolean values.
        /// </summary>
        /// <param name="list">
        /// List of boolean values.
        /// </param>
        /// <returns>
        /// The most common boolean string.
        /// </returns>
        public static bool[] GetModeFromBoolList(List<bool[]> list)
        {
            var counter = new Dictionary<string, int>();
            foreach (var boolString in list.Select(GetBoolStringFromBoolArray))
            {
                if (counter.ContainsKey(boolString))
                {
                    counter[boolString] += 1;
                }
                else
                {
                    counter[boolString] = 1;
                }
            }

            var maximum = 0;
            string result = null;
            foreach (var key in counter.Keys.Where(key => counter[key] > maximum))
            {
                maximum = counter[key];
                result = key;
            }

            return result.Select(chr => chr == '1').ToArray();
        }

        /// <summary>
        /// Converts an array of bools to a boolean string.
        /// </summary>
        /// <param name="boolArray">
        /// An array of bools.
        /// </param>
        /// <returns>
        /// A boolean string.
        /// </returns>
        public static string GetBoolStringFromBoolArray(bool[] boolArray)
        {
            var sb = new StringBuilder();
            foreach (var b in boolArray)
            {
                sb.Append(b ? "1" : "0");
            }

            return sb.ToString();
        }

        /// <summary>
        /// Converts concatenated beta and gamma vectors into separate beta and gamma vector.
        /// </summary>
        /// <param name="bigFreeParamsVector">
        /// Concatenated beta and gamma vectors.
        /// </param>
        /// <returns>
        /// FreeParameters that contains beta and gamma vectors.
        /// </returns>
        /// <remarks>
        /// Useful for getting beta and gamma vectors from a concatenated vector inside the optimized function.
        /// </remarks>
        public static FreeParameters ConvertVectorIntoHalves(double[] bigFreeParamsVector)
        {
            var size = bigFreeParamsVector.Length;
            var vectorTermsNumber = size / 2;
            var freeParamsVector = new FreeParameters
            {

                Beta = bigFreeParamsVector[0..vectorTermsNumber],
                Gamma = bigFreeParamsVector[vectorTermsNumber..(2 * vectorTermsNumber)],

            };

            return freeParamsVector;
        }
    }
}
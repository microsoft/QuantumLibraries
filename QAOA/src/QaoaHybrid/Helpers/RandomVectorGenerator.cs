// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.QAOA.QaoaHybrid.Helpers
{
    using System;

    public class RandomVectorGenerator
    {
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
        public static double[] GenerateRandomVector(int length, double maximumValue)
        {
            var rand = new Random();
            double[] randomVector = new double[length];
            for (var i = 0; i < length; i++)
            {
                randomVector[i] = maximumValue * rand.NextDouble();
            }

            return randomVector;
        }
    }
}

namespace QAOA.QaoaHybrid
{
    using System;
    using System.Collections.Generic;
    using System.Text;

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
    }
}
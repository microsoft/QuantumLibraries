// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.QaoaHybrid.Helpers
{
    using System.Text;

    public class ArrayToStringConverter
    {
        /// <summary>
        /// Converts an array of bools to a boolean string.
        /// </summary>
        /// <param name="boolArray">
        /// An array of bools.
        /// </param>
        /// <returns>
        /// A boolean string.
        /// </returns>
        public static string ConvertBoolArrayToString(bool[] boolArray)
        {
            var sb = new StringBuilder();
            foreach (var b in boolArray)
            {
                sb.Append(b ? "1" : "0");
            }

            return sb.ToString();
        }
    }
}

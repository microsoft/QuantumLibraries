// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Qaoa.QaoaHybrid
{
    using System.Linq;
    using Newtonsoft.Json;
    using Microsoft.Quantum.Qaoa.QaoaHybrid.Helpers;

    public class QaoaParameters
    {
        public double[] Betas { get; set; }
        public double[] Gammas { get; set; }

        public QaoaParameters(double[] concatenatedQaoaParameters)
        {
            var size = concatenatedQaoaParameters.Length;
            var vectorTermsNumber = size / 2;
            Betas = concatenatedQaoaParameters[0..vectorTermsNumber];
            Gammas = concatenatedQaoaParameters[vectorTermsNumber..(2 * vectorTermsNumber)];
        }

        [JsonConstructor]
        public QaoaParameters(double[] betas, double[] gammas)
        {
            Betas = betas;
            Gammas = gammas;
        }

        public QaoaParameters(int p)
        {
           Betas = RandomVectorGenerator.GenerateRandomVector(p, System.Math.PI);
           Gammas = RandomVectorGenerator.GenerateRandomVector(p, 2 * System.Math.PI);
        }

        /// <summary>
        /// Gets betas and gammas vectors as a concatenated vector.
        /// </summary>
        /// <returns>
        /// Array of concatenated betas and gammas arrays.
        /// </returns>
        public double[] ConcatenatedQaoaParameters =>
            Betas.Concat(Gammas).ToArray();
    }
}

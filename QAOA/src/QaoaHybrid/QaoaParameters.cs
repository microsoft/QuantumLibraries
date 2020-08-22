using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Microsoft.Quantum.Math;
using Newtonsoft.Json;

namespace Microsoft.Quantum.QAOA.QaoaHybrid
{
    public class QaoaParameters
    {
        public double[] Beta { get; }
        public double[] Gamma { get; }

        public QaoaParameters(double[] concatenatedQaoaParameters)
        {
            var size = concatenatedQaoaParameters.Length;
            var vectorTermsNumber = size / 2;
            Beta = concatenatedQaoaParameters[0..vectorTermsNumber];
            Gamma = concatenatedQaoaParameters[vectorTermsNumber..(2 * vectorTermsNumber)];
        }

        [JsonConstructor]
        public QaoaParameters(double[] beta, double[] gamma)
        {
            Beta = beta;
            Gamma = gamma;
        }

        public QaoaParameters(int p)
        {
           Beta = Utils.GetRandomVector(p, System.Math.PI);
           Gamma = Utils.GetRandomVector(p, 2 * System.Math.PI);
        }

        /// <summary>
        /// Converts beta and gamma vectors into a concatenated vector.
        /// </summary>
        /// <returns>
        /// Array of concatenated beta and gamma arrays.
        /// </returns>
        public double[] getConcatenatedQaoaParameters()
        {
            return Beta.Concat(Gamma).ToArray();
        }
    }
}

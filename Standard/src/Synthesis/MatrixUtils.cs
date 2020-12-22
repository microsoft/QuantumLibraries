using System.Diagnostics;
using System.Numerics;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Synthesis
{
    internal class MatrixUtils {
        // Checks whether given matrix is unitary.
        public static bool isMatrixUnitary(Complex[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n) return false; // Unitary matrix must be square.

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Complex dotProduct = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        dotProduct += matrix[i, k] * Complex.Conjugate(matrix[j, k]);
                    }
                    Complex expectedDotProduct = (i == j) ? 1.0 : 0.0;
                    if ((dotProduct - expectedDotProduct).Magnitude > 1e-9)
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        // Converts matrix from C# array to Q# array.
        public static QArray<QArray<Quantum.Math.Complex>> matrixToQs(Complex[,] b)
        {
            long n1 = b.GetLength(0);
            long n2 = b.GetLength(1);
            var a = new QArray<Quantum.Math.Complex>[n1];
            for (long i = 0; i < n1; i++)
            {
                var row = new Quantum.Math.Complex[n2];
                for (int j = 0; j < n2; j++)
                {
                    row[j] = new Quantum.Math.Complex((b[i, j].Real, b[i, j].Imaginary));
                }
                a[i] = new QArray<Quantum.Math.Complex>(row);
            }
            return new QArray<QArray<Quantum.Math.Complex>>(a);
        }

        // Converts square matrix from Q# array to C#.
        public static Complex[,] squareMatrixFromQs(IQArray<IQArray<Quantum.Math.Complex>> a)
        {
            long n = a.Length;
            var b = new Complex[n, n];
            for (long i = 0; i < n; i++)
            {
                Debug.Assert(a[i].Length == n, "Matrix is not square");
                for (long j = 0; j < n; j++)
                {
                    b[i, j] = new Complex(a[i][j].Real, a[i][j].Imag);
                }
            }
            return b;
        }
    }
}
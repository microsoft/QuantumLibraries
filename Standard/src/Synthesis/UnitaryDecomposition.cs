// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
//using static System.Math;
//using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
//using System.Diagnostics.CodeAnalysis;
//using Microsoft.Quantum.Diagnostics.Emulation;
//using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
//using Microsoft.Quantum.Simulation.Simulators;

namespace Microsoft.Quantum.Synthesis
{

    internal partial class TwoLevelDecomposition
    {
        public class Native : TwoLevelDecomposition
        {
            public Native(IOperationFactory m) : base(m) { }

            // Represents a square matrix which is an identity matrix with elements on positions
            // (i1, i1), (i1, i2), (i2, i1), (i2, i2) replaced with elements from mx. 
            internal class TwoLevelUnitary
            {
                private Complex[,] mx;   // 2x2 non-trivial principal submatrix.
                private int i1, i2;      // Indices of non-trivial submatrix.

                public TwoLevelUnitary(Complex[,] mx, int i1, int i2)
                {
                    this.mx = mx;
                    this.i1 = i1;
                    this.i2 = i2;
                }

                public void applyPermutation(int[] perm)
                {
                    i1 = perm[i1];
                    i2 = perm[i2];
                }

                // Ensures that index1 < index2.
                public void orderIndices()
                {
                    if (this.i1 > this.i2)
                    {
                        (i1, i2) = (i2, i1);
                        (mx[0, 0], mx[1, 1]) = (mx[1, 1], mx[0, 0]);
                        (mx[0, 1], mx[1, 0]) = (mx[1, 0], mx[0, 1]);
                    }
                }

                // Equivalent to inversion for unitary matrix.
                public void conjugateTranspose()
                {
                    mx[0, 0] = Complex.Conjugate(mx[0, 0]);
                    mx[1, 1] = Complex.Conjugate(mx[1, 1]);
                    (mx[0, 1], mx[1, 0]) = (Complex.Conjugate(mx[1, 0]), Complex.Conjugate(mx[0, 1]));
                }

                // Applies A = A * M, where M is this matrix.
                public void multiplyRight(Complex[,] A)
                {
                    int n = A.GetLength(0);
                    for (int i = 0; i < n; i++)
                    {
                        (A[i, i1], A[i, i2]) = (A[i, i1] * mx[0, 0] + A[i, i2] * mx[1, 0], A[i, i1] * mx[0, 1] + A[i, i2] * mx[1, 1]);
                    }
                }

                public bool isIdentity()
                {
                    // TODO: implement.
                    return false;
                }

                // Converts square matrix from C# to Q#.
                // TODO: can inline in toQsharp() and get rid of loops.
                private static QArray<QArray<Quantum.Math.Complex>> matrixToQs(Complex[,] b)
                {
                    long n = b.GetLength(0);
                    var a = new QArray<Quantum.Math.Complex>[n];
                    for (long i = 0; i < n; i++)
                    {
                        var row = new Quantum.Math.Complex[n];
                        for (int j = 0; j < n; j++)
                        {
                            row[j] = new Quantum.Math.Complex((b[i, j].Real, b[i, j].Imaginary));
                        }
                        a[i] = new QArray<Quantum.Math.Complex>(row);
                    }
                    return new QArray<QArray<Quantum.Math.Complex>>(a);
                }

                // Converts to tuple to be passed to Q#.
                public (IQArray<IQArray<Quantum.Math.Complex>>, long, long) toQsharp()
                {
                    return (matrixToQs(this.mx), this.i1, this.i2);
                }
            }

            // Returns unitary 2x2 matrix U, s.t. [a, b] U = [c, 0].
            private static Complex[,] makeEliminatingMatrix(Complex a, Complex b)
            {
                Debug.Assert(a.Magnitude > 1e-9 && b.Magnitude > 1e-9);
                double theta = System.Math.Atan((b / a).Magnitude);
                double lmbda = -a.Phase;
                double mu = System.Math.PI + b.Phase - a.Phase - lmbda;
                return new Complex[,] {{
                    System.Math.Cos(theta) * Complex.FromPolarCoordinates(1.0, lmbda),
                    System.Math.Sin(theta) * Complex.FromPolarCoordinates(1.0, mu)
                }, {
                    -System.Math.Sin(theta) * Complex.FromPolarCoordinates(1.0, -mu),
                    System.Math.Cos(theta) * Complex.FromPolarCoordinates(1.0, -lmbda)
                }};
            }

            // Returns list of two-level unitary matrices, which multiply to A.
            //
            // Matrices are listed in application order.
            // Every matrix has indices differring in 1.
            // A is modified as result of this function.
            private static List<TwoLevelUnitary> twoLevelDecompose(Complex[,] A)
            {
                int n = A.GetLength(0);
                var result = new List<TwoLevelUnitary>();

                for (int i = 0; i < n - 2; i++)
                {
                    for (int j = n - 1; j > i; j--)
                    {
                        if (A[i, j].Magnitude < 1e-9)
                        {
                            // Element is already zero, skipping.
                        }
                        else
                        {
                            Complex[,] mx;
                            if (A[i, j - 1].Magnitude < 1e-9)
                            {
                                // Just swap columns with Pauli X matrix.
                                mx = new Complex[,] { { 0.0, 1.0 }, { 1.0, 0.0 } };
                            }
                            else
                            {
                                mx = makeEliminatingMatrix(A[i, j - 1], A[i, j]);
                            }
                            var u2x2 = new TwoLevelUnitary(mx, j - 1, j);
                            u2x2.multiplyRight(A);
                            u2x2.conjugateTranspose();
                            result.Add(u2x2);
                        }
                    }
                }

                var lastMatrix = new TwoLevelUnitary(new Complex[,] { { A[n - 2, n - 2], A[n - 2, n - 1] }, { A[n - 1, n - 2], A[n - 1, n - 1] } }, n - 2, n - 1);
                if (!lastMatrix.isIdentity())
                {
                    result.Add(lastMatrix);
                }
                return result;
            }

            // Builds binary-reflected gray code.
            // Returns permutation of numbers from 0 to n-1 such as any two consequent number differ in exactly one bit.
            // n must be power of 2.
            private static int[] grayCode(int n)
            {
                var result = new int[n];
                for (int i = 0; i < n; i++)
                {
                    result[i] = i ^ (i / 2);
                }
                return result;
            }

            // Applies permutation to columns and rows of matrix.
            private static void permuteMatrix(Complex[,] matrix, int[] perm)
            {
                int n = perm.Length;
                Debug.Assert(matrix.GetLength(0) == n);
                Debug.Assert(matrix.GetLength(1) == n);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        matrix[i, j] = matrix[i, perm[j]];
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        matrix[j, i] = matrix[perm[j], i];
                    }
                }
            }

            private static bool isPowerOfTwo(int x)
            {
                // TODO: implement.
                return true;
            }

            private static bool isUnitary(Complex[,] matrix)
            {
                // TODO: implement.
                return true;
            }


            // Returns list of two-level unitary matrices, which multiply to A.
            //
            // Every matrix has indices differing in exactly 1 bit (this is achieved with Gray code).
            private static List<TwoLevelUnitary> twoLevelDecomposeGray(Complex[,] A)
            {
                int n = A.GetLength(0);
                Debug.Assert(isPowerOfTwo(n), "Matrix sizeis not power of two.");
                Debug.Assert(A.GetLength(1) == n, "Matrix is not square.");
                Debug.Assert(isUnitary(A), "Matrix is not unitary.");

                int[] perm = grayCode(n);
                permuteMatrix(A, perm);
                List<TwoLevelUnitary> result = twoLevelDecompose(A);
                foreach (TwoLevelUnitary matrix in result)
                {
                    matrix.applyPermutation(perm);
                }
                return result;
            }


            // Converts square matrix from Q# to C#.
            private static Complex[,] matrixFromQs(IQArray<IQArray<Quantum.Math.Complex>> a)
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

            private IQArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)> DoThis(IQArray<IQArray<Quantum.Math.Complex>> unitary)
            {
                Complex[,] a = matrixFromQs(unitary);
                List<TwoLevelUnitary> matrices = twoLevelDecomposeGray(a);
                var result = new List<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>();
                foreach (TwoLevelUnitary matrix in matrices)
                {
                    matrix.orderIndices();
                    result.Add(matrix.toQsharp());
                }
                return new QArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>(result);
            }

            // Override __Body__ property to use C# function
            public override Func<IQArray<IQArray<Quantum.Math.Complex>>, IQArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>> __Body__ => DoThis;
        }
    }
}

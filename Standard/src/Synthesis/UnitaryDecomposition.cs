// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Synthesis
{
    internal partial class TwoLevelDecomposition
    {
        public class Native : TwoLevelDecomposition
        {
            public Native(IOperationFactory m) : base(m) { }

            // Returns special unitary 2x2 matrix U, s.t. [a, b] U = [c, 0].
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

                var lastMatrix = new TwoLevelUnitary(new Complex[,] {
                    { A[n - 2, n - 2], A[n - 2, n - 1] },
                    { A[n - 1, n - 2], A[n - 1, n - 1] } }, n - 2, n - 1);
                if (!lastMatrix.isIdentity())
                {
                    result.Add(lastMatrix);
                }
                return result;
            }

            // Returns permutation of numbers from 0 to n-1 such as any two consequent number 
            // differ in exactly one bit.
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
            private static Complex[,] permuteMatrix(Complex[,] matrix, int[] perm)
            {
                int n = perm.Length;
                Debug.Assert(matrix.GetLength(0) == n);
                Debug.Assert(matrix.GetLength(1) == n);
                var result = new Complex[n, n];
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        result[i, j] = matrix[perm[i], perm[j]];
                    }
                }
                return result;
            }

            // Returns list of two-level unitary matrices, which multiply to A.
            // Every matrix has indices differing in exactly 1 bit.
            private static List<TwoLevelUnitary> twoLevelDecomposeGray(Complex[,] A)
            {
                int n = A.GetLength(0);
                Debug.Assert(A.GetLength(1) == n, "Matrix is not square.");
                Debug.Assert(MatrixUtils.isMatrixUnitary(A), "Matrix is not unitary.");

                int[] perm = grayCode(n);
                A = permuteMatrix(A, perm);
                List<TwoLevelUnitary> result = twoLevelDecompose(A);
                foreach (TwoLevelUnitary matrix in result)
                {
                    matrix.applyPermutation(perm);
                }
                return result;
            }

            // Decomposes unitary matrix into product of 2-level unitary matrices.
            private IQArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)> Decompose(
                IQArray<IQArray<Quantum.Math.Complex>> unitary)
            {
                Complex[,] a = MatrixUtils.squareMatrixFromQs(unitary);
                List<TwoLevelUnitary> matrices = twoLevelDecomposeGray(a);
                var result = new List<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>();
                foreach (TwoLevelUnitary matrix in matrices)
                {
                    matrix.orderIndices();
                    result.Add(matrix.toQsharp());
                }
                return new QArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>(result);
            }

            public override Func<IQArray<IQArray<Quantum.Math.Complex>>,
                IQArray<(IQArray<IQArray<Quantum.Math.Complex>>, long, long)>> __Body__ =>
                Decompose;
        }
    }
}

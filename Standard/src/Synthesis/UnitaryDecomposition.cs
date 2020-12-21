// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

using System;
//using static System.Math;
// TODO: trim not needed stuff.
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using Microsoft.Quantum.Diagnostics.Emulation;
using Microsoft.Quantum.Simulation.Common;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Math;

namespace Microsoft.Quantum.Synthesis {
    public partial class TwoLevelDecomposition {
        public class Native : TwoLevelDecomposition {
            public Native(IOperationFactory m) : base(m) {}



            private Complex[,] matrixFromQs(IQArray<IQArray<Complex>> a) {
                long n = a.Length;
                Complex[,] b = new Complex[n, n];
                for(long i=0;i<n;i++) {
                    Debug.Assert(a[i].Length == n, "Matrix is not square");
                    for (long j=0;j<n;j++) {
                        b[i,j] = a[i][j];
                    }
                }
                return b;
            }

            private QArray<QArray<Complex>> matrixToQs(Complex[,] b) {
                long n = b.GetLength(0);
                QArray<Complex>[] a = new QArray<Complex>[n];
                for(long i=0;i<n;i++) {
                    Complex[] row = new Complex[n];
                    for(int j = 0;j<n;j++) {
                        row[j] = b[i,j];
                    }
                    a[i] = new QArray<Complex>(row);
                }
                return new QArray<QArray<Complex>>(a);
            }

            private QArray<QArray<QArray<Complex>>> matricesToQs(List<Complex[,]> matrices) {
                int n = (int)matrices.Count;
                QArray<QArray<Complex>>[] result = new QArray<QArray<Complex>>[n];
                for(int i=0;i<n;i++) {
                    result[i] = matrixToQs(matrices[i]);
                }
                return new QArray<QArray<QArray<Complex>>>(result);
            }

            private (IQArray<IQArray<IQArray<Complex>>>, IQArray<IQArray<long>>) TwoLevelDecomposeGray(IQArray<IQArray<Complex>> unitary) {
                Complex[,] a = matrixFromQs(unitary);
                List<Complex[,]> r1 = new List<Complex[,]>();
                r1.Add(a);
                
                QArray<QArray<QArray<Complex>>> r2 = matricesToQs(r1);

                QArray<long> r3 = new QArray<long> ( new long[] {0, 1});
                QArray<QArray<long>> r4 = new QArray<QArray<long>>(new QArray<long>[] {r3});

                return (r2, r4);
            }

            // Override __Body__ property to use C# function
            public override Func<IQArray<IQArray<Complex>>, (IQArray<IQArray<IQArray<Complex>>>, IQArray<IQArray<long>>)> __Body__ => TwoLevelDecomposeGray;
        }
    }
}
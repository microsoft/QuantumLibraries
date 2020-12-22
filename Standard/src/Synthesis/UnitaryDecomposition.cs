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
    internal partial class TwoLevelDecomposition {
        public class Native : TwoLevelDecomposition {
            public Native(IOperationFactory m) : base(m) {}


            private Complex[,] matrixFromQs(IQArray<IQArray<Complex>> a) {
                long n = a.Length;
                var b = new Complex[n, n];
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
                var a = new QArray<Complex>[n];
                for(long i=0;i<n;i++) {
                    var row = new Complex[n];
                    for(int j = 0;j<n;j++) {
                        row[j] = b[i,j];
                    }
                    a[i] = new QArray<Complex>(row);
                }
                return new QArray<QArray<Complex>>(a);
            }

            private QArray<QArray<QArray<Complex>>> matricesToQs(List<Complex[,]> matrices) {
                int n = matrices.Count;
                var result = new QArray<QArray<Complex>>[n];
                for(int i=0;i<n;i++) {
                    result[i] = matrixToQs(matrices[i]);
                }
                return new QArray<QArray<QArray<Complex>>>(result);
            }

            private IQArray<(IQArray<IQArray<Complex>>, long, long)> TwoLevelDecomposeGray(IQArray<IQArray<Complex>> unitary) {
                Complex[,] a = matrixFromQs(unitary);
                var result = new List<(IQArray<IQArray<Complex>>, long, long)>(); 
                result.Add((matrixToQs(a), 0, 1));
                
                return new QArray<(IQArray<IQArray<Complex>>, long, long)>(result);
            }

            // Override __Body__ property to use C# function
            public override Func<IQArray<IQArray<Complex>>, IQArray<(IQArray<IQArray<Complex>>, long, long)>> __Body__ => TwoLevelDecomposeGray;
        }
    }
}

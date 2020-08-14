// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

// Portions of this file are derived from NumPy's binomial
// sampling implementation, under the following license:

/* Copyright 2005 Robert Kern (robert.kern@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

using System;
using static System.Math;

namespace Microsoft.Quantum.Standard.Emulation
{
    /// <summary>
    ///      Represents the binomial distribution $\operatorname{Bin}(N, p)$,
    ///      and provides samples according to this distribution.
    /// </summary>
    internal class BinomialDistribution
    {
        public System.Random RandomState { get; }
        public long NSamples { get; }
        public double SuccessProbability { get; }

        public readonly double FailureProbability;
        public readonly double Expectation;
        public readonly double Variance;
        public readonly double StandardDeviation;

        // Some useful numeric parameters to have precomputed.
        private readonly double qn;
        private readonly long bound;

        public BinomialDistribution ComplementaryDistribution =>
            new BinomialDistribution(NSamples, FailureProbability, RandomState);

        public BinomialDistribution(long nSamples, double successProbability, System.Random randomState)
        {
            NSamples = nSamples;
            SuccessProbability = successProbability;
            RandomState = randomState;

            FailureProbability = 1.0 - SuccessProbability;
            Expectation = NSamples * SuccessProbability;
            Variance = NSamples * SuccessProbability * FailureProbability;
            StandardDeviation = Sqrt(Variance);

            qn = Exp(NSamples * Log(FailureProbability));
            bound = (long)Min(NSamples, Expectation + 10.0 * Sqrt(Variance + 1));
        }

        public BinomialDistribution(long nSamples, double successProbability)
        : this(nSamples, successProbability, new Random())
        { }

        public long NextSample()
        {
            if (SuccessProbability > 0.5)
            {
                return NSamples - ComplementaryDistribution.NextSample();
            }
            else if (Expectation <= 30.0)
            {
                return SampleUsingInversion();
            }
            else
            {
                return SampleUsingBTPE();
            }
        }

        private long SampleUsingInversion()
        {
            var px = qn;
            var U = RandomState.NextDouble();
            var X = 0L;

            while (U > px)
            {
                X++;
                if (X > bound)
                {
                    X = 0;
                    px = qn;
                    U = RandomState.NextDouble();
                }
                else
                {
                    U -= px;
                    px  = ((NSamples - X + 1) * SuccessProbability * px) / (X * FailureProbability);
                }
            }
            return X;
        }

        private long SampleUsingBTPE()
        {
            double r,q,fm,p1,xm,xl,xr,c,laml,lamr,p2,p3,p4;
            double a,u,v,s,F,rho,t,A,nrq,x1,x2,f1,f2,z,z2,w,w2,x;
            long m,y,k,i;

            /* initialize */
            r = Min(SuccessProbability, FailureProbability);
            q = 1.0 - r;
            fm = NSamples * r + r;
            m = (long)Floor(fm);
            p1 = Floor(2.195 * Sqrt(NSamples * r * q) - 4.6 * q) + 0.5;
            xm = m + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5/(15.3 + m);
            a = (fm - xl)/(fm-xl*r);
            laml = a*(1.0 + a/2.0);
            a = (xr - fm)/(xr*q);
            lamr = a*(1.0 + a/2.0);
            p2 = p1*(1.0 + 2.0*c);
            p3 = p2 + c/laml;
            p4 = p3 + c/lamr;


        /* sigh ... */
        Step10:
            nrq = NSamples * r * q;
            u = RandomState.NextDouble() * p4;
            v = RandomState.NextDouble();
            if (u > p1) goto Step20;
            y = (long)Floor(xm - p1*v + u);
            goto Step60;

        Step20:
            if (u > p2) goto Step30;
            x = xl + (u - p1)/c;
            v = v*c + 1.0 - Abs(m - x + 0.5)/p1;
            if (v > 1.0) goto Step10;
            y = (long)Floor(x);
            goto Step50;

        Step30:
            if (u > p3) goto Step40;
            y = (long)Floor(xl + Log(v)/laml);
            /* Reject if v == 0.0 since cast of inf not well defined */
            if ((y < 0) || (v == 0.0)) goto Step10;
            v = v*(u-p2)*laml;
            goto Step50;

        Step40:
            y = (long)Floor(xr - Log(v)/lamr);
            /* Reject if v == 0.0 since cast of inf not well defined */
            if ((y > NSamples) || (v == 0.0)) goto Step10;
            v = v*(u-p3)*lamr;

        Step50:
            k = Abs(y - m);
            if ((k > 20) && (k < ((nrq)/2.0 - 1))) goto Step52;

            s = r/q;
            a = s*(NSamples+1);
            F = 1.0;
            if (m < y)
            {
                for (i=m+1; i<=y; i++)
                {
                    F *= (a/i - s);
                }
            }
            else if (m > y)
            {
                for (i=y+1; i<=m; i++)
                {
                    F /= (a/i - s);
                }
            }
            if (v > F) goto Step10;
            goto Step60;

            Step52:
            rho = (k/(nrq))*((k*(k/3.0 + 0.625) + 0.16666666666666666)/nrq + 0.5);
            t = -k*k/(2*nrq);
            A = Log(v);
            if (A < (t - rho)) goto Step60;
            if (A > (t + rho)) goto Step10;

            x1 = y+1;
            f1 = m+1;
            z = NSamples + 1-m;
            w = NSamples - y+1;
            x2 = x1*x1;
            f2 = f1*f1;
            z2 = z*z;
            w2 = w*w;
            if (A > (xm*Log(f1/x1)
                + (NSamples - m+0.5)*Log(z/w)
                + (y-m)*Log(w*r/(x1*q))
                + (13680.0-(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0
                + (13680.0-(462.0-(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0
                + (13680.0-(462.0-(132.0-(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0
                + (13680.0-(462.0-(132.0-(99.0-140.0/w2)/w2)/w2)/w2)/w/166320.0))
            {
                goto Step10;
            }

        Step60:
            if (SuccessProbability > 0.5)
            {
                y = NSamples - y;
            }

            return y;
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.Quantum.Simulation.Core;

namespace Microsoft.Quantum.Canon.Native
{
    public class Sin : Microsoft.Quantum.Canon.Sin
    {

        public Sin(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Sin(theta);
        
    }

    public class Cos : Microsoft.Quantum.Canon.Cos
    {

        public Cos(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Cos(theta);

    }

    public class Tan : Microsoft.Quantum.Canon.Tan
    {

        public Tan(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Tan(theta);

    }

    public class ArcSin : Microsoft.Quantum.Canon.ArcSin
    {
        
        public ArcSin(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Asin(theta);

    }

    public class ArcCos : Microsoft.Quantum.Canon.ArcCos
    {
        
        public ArcCos(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Acos(theta);

    }

    public class ArcTan : Microsoft.Quantum.Canon.ArcTan
    {

        public ArcTan(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Atan(theta);

    }

    public class ArcTan2 : Microsoft.Quantum.Canon.ArcTan2
    {

        public ArcTan2(IOperationFactory m) : base(m) { }

        public override Func<(double, double), double> Body =>
            (args) =>
                System.Math.Atan2(args.Item1, args.Item2);

    }

    public class Log : Microsoft.Quantum.Canon.Log
    {

        public Log(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Log(theta);

    }

    public class Ceiling : Microsoft.Quantum.Canon.Ceiling
    {

        public Ceiling(IOperationFactory m) : base(m) { }

        public override Func<double, long> Body =>
            (theta) =>
                // FIXME: add range checking!
                (long) System.Math.Ceiling(theta);

    }

    public class Sqrt : Microsoft.Quantum.Canon.Sqrt
    {

        public Sqrt(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Sqrt(theta);

    }

    public class Cosh : Microsoft.Quantum.Canon.Cosh
    {

        public Cosh(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Cosh(theta);

    }

    public class Sinh : Microsoft.Quantum.Canon.Sinh
    {

        public Sinh(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Sinh(theta);

    }

    public class Tanh : Microsoft.Quantum.Canon.Tanh
    {

        public Tanh(IOperationFactory m) : base(m) { }

        public override Func<double, double> Body =>
            (theta) =>
                System.Math.Tanh(theta);

    }

    public class Floor : Microsoft.Quantum.Canon.Floor
    {

        public Floor(IOperationFactory m) : base(m) { }

        public override Func<double, long> Body =>
            (theta) =>
                // FIXME: add range checking!
                (long) System.Math.Floor(theta);

    }

    public class ToDouble : Microsoft.Quantum.Canon.ToDouble
    {

        public ToDouble(IOperationFactory m) : base(m) { }

        public override Func<long, double> Body =>
            (arg) =>
                (double) arg;

    }
}

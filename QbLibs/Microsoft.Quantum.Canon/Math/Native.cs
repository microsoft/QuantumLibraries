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
        public Sin(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<double, double> Body
        {
            get => (theta) =>
                System.Math.Sin(theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, Double theta)
        {
            return m.Run<Sin, double, double>(theta);
        }
    }

    public class Cos : Microsoft.Quantum.Canon.Cos
    {
        public Cos(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<double, double> Body
        {
            get => (theta) =>
                System.Math.Cos(theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, Double theta)
        {
            return m.Run<Cos, double, double>(theta);
        }
    }

    public class ArcSin : Microsoft.Quantum.Canon.ArcSin
    {
        public ArcSin(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<double, double> Body
        {
            get => (theta) =>
                System.Math.Asin(theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, Double theta)
        {
            return m.Run<ArcSin, double, double>(theta);
        }
    }

    public class ArcCos : Microsoft.Quantum.Canon.ArcCos
    {
        public ArcCos(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<double, double> Body
        {
            get => (theta) =>
                System.Math.Acos(theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, Double theta)
        {
            return m.Run<ArcCos, double, double>(theta);
        }
    }

    public class Sqrt : Microsoft.Quantum.Canon.Sqrt
    {
        public Sqrt(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<double, double> Body
        {
            get => (theta) =>
                System.Math.Sqrt(theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, Double theta)
        {
            return m.Run<Sqrt, double, double>(theta);
        }
    }

    public class ToDouble : Microsoft.Quantum.Canon.ToDouble
    {
        public ToDouble(IOperationFactory m) : base(m)
        {
            this.Dependencies = new Type[] { };
        }

        public override Type[] Dependencies
        {
            get;
        }

        public override Func<long, double> Body
        {
            get => (theta) => ((double) theta);
        }

        public static System.Threading.Tasks.Task<double> Run(IOperationFactory m, long theta)
        {
            return m.Run<ToDouble, long, double>(theta);
        }
    }


}

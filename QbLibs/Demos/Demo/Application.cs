using System;
using System.Collections.Generic;
using System.Text;
using Microsoft.Quantum.Simulation.Core;
using Microsoft.Quantum.Simulation.Simulators;
using Microsoft.Quantum.Canon;

using System.Windows.Forms;

using OxyPlot;
using OxyPlot.Series;

namespace Microsoft.Quantum.Examples.H2ExampleApplication
{

    class MainViewModel
    {
        public MainViewModel()
        {
            this.Title = "Example 2";
            this.Points = new List<DataPoint>
                              {
                                  new DataPoint(0, 4),
                                  new DataPoint(10, 13),
                                  new DataPoint(20, 15),
                                  new DataPoint(30, 16),
                                  new DataPoint(40, 12),
                                  new DataPoint(50, 12)
                              };
        }

        public string Title { get; private set; }

        public IList<DataPoint> Points { get; private set; }

    }

    class Program
    {
        static void Main(string[] args)
        {
            var qsim = new QuantumSimulator();
            qsim.Register(typeof(Microsoft.Quantum.Canon.Ceiling), typeof(Microsoft.Quantum.Canon.Native.Ceiling));
            qsim.Register(typeof(Microsoft.Quantum.Canon.ArcTan2), typeof(Microsoft.Quantum.Canon.Native.ArcTan2));
            qsim.Register(typeof(Microsoft.Quantum.Canon.ToDouble), typeof(Microsoft.Quantum.Canon.Native.ToDouble));

            //H2EstimateEnergy(idxBondLength: Int, trotterStepSize: Double, bitsPrecision: Int)
            for (int idxBondLength = 0; idxBondLength < 54; idxBondLength++)
            {
                //var phaseSet = 0.34545;
                //Repeat 3 times, take lowest energy
                var phaseEst = (double)0.0;
                for (int rep = 0; rep < 3; rep++)
                {
                    phaseEst =  Math.Min(phaseEst, H2EstimateEnergy.Run(qsim, idxBondLength, 0.5, 8).Result);
                }
                Console.WriteLine("Estimated energy in Hartrees is : {0}", phaseEst);
            }
            

            //PlotModel.InvalidatePlot(true)
            Console.ReadLine();
        }
    }
}

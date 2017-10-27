using Microsoft.Quantum.Simulation.Simulators;
using System.Linq;

namespace Microsoft.Quantum.Examples.TeleportDemo
{
    class Program
    {
        static void Main(string[] args)
        {
            var sim = new QuantumSimulator();
            var rand = new System.Random();

            foreach (var idxRun in Enumerable.Range(0, 10)) {
                var message = rand.Next(2) == 0;
                var task = Teleport.Run(sim, message);
                System.Console.WriteLine($"{message} =====> {task.GetAwaiter().GetResult()}\n");
            }

            System.Console.ReadLine();
        }
    }
}

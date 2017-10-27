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

            foreach (var idxRun in Enumerable.Range(0, 8)) {
                var sent = rand.Next(2) == 0;
                var recieved = Teleport.Run(sim, sent).Result;
                System.Console.WriteLine($"Round {idxRun}:\tSent {sent},\tgot {recieved}.");
                System.Console.WriteLine(sent == recieved ? "Teleportation successful!!\n" : "\n");
            }

            System.Console.ReadLine();
        }
    }
}

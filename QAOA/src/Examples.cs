namespace Quantum.QAOA
{
    using global::QAOA.QaoaHybrid;
    using System;


    class Examples
    {
        static void Main(string[] args)
        {
            //PARAMETERS
            int numberOfIterations = 50;
            int p = 3;
            int numberOfRandomStartingPoints = 1;

            //EXAMPLES

            //Quantum Santa (http://quantumalgorithmzoo.org/traveling_santa/)
            double[] dtx = { 0.619193, 0.742566, 0.060035, -1.568955, 0.045490 };
            double[] dtz = { 3.182203, -1.139045, 0.221082, 0.537753, -0.417222 };
            double[] segmentCosts = { 4.70, 9.09, 9.03, 5.70, 8.02, 1.71 };
            double[] dh = { 4 * 20 - 0.5 * 4.7, 4 * 20 - 0.5 * 9.09, 4 * 20 - 0.5 * 9.03, 4 * 20 - 0.5 * 5.70, 4 * 20 - 0.5 * 8.02, 4 * 20 - 0.5 * 1.71 };
            double[] dJ = { 40.0,40.0,20.0,40.0,40.0,40.0,
                            40.0,40.0,40.0,20.0,40.0,40.0,
                            40.0,40.0,40.0,40.0,40.0,40.0,
                            40.0,40.0,40.0,40.0,40.0,40.0,
                            40.0,40.0,40.0,40.0,40.0,20.0,
                            40.0,40.0,40.0,40.0,40.0,40.0};
            ProblemInstance quantumSanta = new ProblemInstance(dh, dJ);


            //MaxCut (medium.com/mdr-inc/qaoa-maxcut-using-blueqat-aaf33038f46e)
            dh = new Double[] { 0,0,0,0,0};
            dJ = new Double[]{ 0,1,0,1,0,
                               0,0,1,0,0,
                               0,0,0,1,1,
                               0,0,0,0,1,
                               0,0,0,0,0};
            ProblemInstance maxCut1 = new ProblemInstance(dh, dJ);
            

            //Rigetti MaxCut unit tests
            dh = new Double[]{-0.5,0,-1,0.5};
            dJ = new Double[]{0,1,2,0,
                              0,0,0.5,0,
                              0,0,0,2.5,
                              0,0,0,0};
            ProblemInstance maxCut2 = new ProblemInstance(dh, dJ);

            
            dh = new Double[] { 0.8, -0.5 };
            dJ = new Double[]{ 0, -1,
                               0, 0};
            ProblemInstance maxCut3 = new ProblemInstance(dh, dJ);

            dh = new Double[] {0, 0 };
            dJ = new Double[]{ 0, 1,
                               0, 0};
            ProblemInstance maxCut4 = new ProblemInstance(dh, dJ);

            //END EXAMPLES

            HybridQaoa cop = new HybridQaoa(numberOfIterations, p, maxCut4, numberOfRandomStartingPoints, true);

            OptimalSolution res = cop.RunOptimization();
            Console.WriteLine(res.optimalVector);

            }
    }
}
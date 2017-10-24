namespace Microsoft.Quantum.Examples

open System
open System.Windows
open System.Windows.Controls
open System.Windows.Markup
open FsXaml
open OxyPlot
open System.ComponentModel

type MainViewModel() = class
    let mutable title = "H₂ Simulation"
    let mutable points = [
        DataPoint(0.0, 4.0)
        DataPoint(10.0, 13.0)
        DataPoint(20.0, 15.0)
        DataPoint(30.0, 16.0)
        DataPoint(40.0, 12.0)
        DataPoint(50.0, 12.0)
    ]
    let propertyChanged = new Event<_, _>()
    interface INotifyPropertyChanged with
        [<CLIEvent>]
        member this.PropertyChanged = propertyChanged.Publish

    member this.Title
        with get() =
            title
        and private set(value) =
            title <- value
            propertyChanged.Trigger(this, PropertyChangedEventArgs("Title"))

    member this.Points
        with get() =
            points
        and set(value) =
            points <- value
            propertyChanged.Trigger(this, PropertyChangedEventArgs("Points"))

        
end


module H2PlottingDemo =
    open Microsoft.Quantum.Simulation.Simulators
    open Microsoft.Quantum.Canon

    type MainWindow = XAML<"MainWindow.xaml">

    let estimateEnergies =
        let qsim = QuantumSimulator()
        qsim.Register(typeof<Microsoft.Quantum.Canon.Ceiling>, typeof<Microsoft.Quantum.Canon.Native.Ceiling>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ArcTan2>, typeof<Microsoft.Quantum.Canon.Native.ArcTan2>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ToDouble>, typeof<Microsoft.Quantum.Canon.Native.ToDouble>);

        let H2EstimateEnergy = qsim.Get(typeof<H2EstimateEnergy>) :?> H2EstimateEnergy

        let estAtBondLength idx =
            // TODO: repeat and take lowest.
            H2EstimateEnergy.Body.Invoke (struct (idx, 0.5, int64 8))

        [0..10] // fixme: change to 53
            |> Seq.map (int64 >> estAtBondLength)
            |> Seq.mapi (fun idx est -> DataPoint(float idx, est))
            |> List.ofSeq

            ////H2EstimateEnergy(idxBondLength: Int, trotterStepSize: Double, bitsPrecision: Int)
            //for (int idxBondLength = 0; idxBondLength < 54; idxBondLength++)
            //{
            //    //var phaseSet = 0.34545;
            //    //Repeat 3 times, take lowest energy
            //    var phaseEst = (double)0.0;
            //    for (int rep = 0; rep < 3; rep++)
            //    {
            //        phaseEst =  Math.Min(phaseEst, H2EstimateEnergy.Run(qsim, idxBondLength, 0.5, 8).Result);
            //    }
            //    Console.WriteLine("Estimated energy in Hartrees is : {0}", phaseEst);
            //}
            

            ////PlotModel.InvalidatePlot(true)
            //Console.ReadLine();

    [<STAThread>]
    [<EntryPoint>]
    let main argv = 
        let mainWindow = MainWindow()
 
        let ctx = mainWindow.DataContext
        //ctx.
        let application = new Application()
        application.Run(mainWindow)


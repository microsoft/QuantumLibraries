namespace Microsoft.Quantum.Examples

open System
open System.Windows.Forms
open System.Windows.Controls
open System.Windows.Markup

module H2PlottingDemo =
    open Microsoft.Quantum.Simulation.Simulators
    open Microsoft.Quantum.Canon
    open FSharp.Charting
    open FSharp.Control
    open FSharp.Charting
    open Microsoft.Quantum.Simulation.Core

    let estimateEnergies =
        let qsim = QuantumSimulator()
        qsim.Register(typeof<Microsoft.Quantum.Canon.Ceiling>, typeof<Microsoft.Quantum.Canon.Native.Ceiling>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ArcTan2>, typeof<Microsoft.Quantum.Canon.Native.ArcTan2>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ToDouble>, typeof<Microsoft.Quantum.Canon.Native.ToDouble>);

        let H2EstimateEnergy = qsim.Get(typeof<H2EstimateEnergy>) :?> H2EstimateEnergy
        let H2BondLengths = qsim.Get(typeof<H2BondLengths>) :?> H2BondLengths

        let estAtBondLength idx =
            // TODO: repeat and take lowest.
            H2EstimateEnergy.Body.Invoke (struct (idx, 0.5, int64 8))

        let bondLengths =
            H2BondLengths.Body.Invoke QVoid.Instance

        let data =
            [0..53] // fixme: change to 53
            |> Seq.map (int64 >> estAtBondLength)
            |> Seq.zip bondLengths

        data

    [<STAThread>]
    [<EntryPoint>]
    let main argv = 
        Application.EnableVisualStyles()
        Application.SetCompatibleTextRenderingDefault false
        
        let chart =
            estimateEnergies
            |> Chart.Line
        chart.ShowChart()
            |> Application.Run
        0

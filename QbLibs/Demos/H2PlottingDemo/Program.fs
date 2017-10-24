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

        let H2EstimateEnergyRPE = qsim.Get(typeof<H2EstimateEnergyRPE>) :?> H2EstimateEnergyRPE
        let H2BondLengths = qsim.Get(typeof<H2BondLengths>) :?> H2BondLengths

        let estAtBondLength idx =
            // TODO: repeat and take lowest.
            [0..3]
            |> Seq.map (fun idxRep ->
                    H2EstimateEnergyRPE.Body.Invoke (struct (idx, int64 6, float 1))
                )
            |> Seq.min

        let bondLengths =
            H2BondLengths.Body.Invoke QVoid.Instance

        let data =
            [0..53] // fixme: change to 53
            |> Seq.map (int64 >> estAtBondLength)
            |> Seq.zip bondLengths
            |> Seq.map(fun (bondLen, est) -> 
                    printfn "%A %A" bondLen est;
                    est
                )

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

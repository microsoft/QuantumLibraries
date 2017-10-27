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
    open FSharp.Control.Reactive
    open Microsoft.FSharp.Core
    open Microsoft.FSharp.Collections
    open System.Reactive
    open System.Windows
    open FSharp.Charting.ChartTypes

    let registerCallables (qsim : AbstractFactory<Operation>) =        
        qsim.Register(typeof<Microsoft.Quantum.Canon.Ceiling>, typeof<Microsoft.Quantum.Canon.Native.Ceiling>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ArcTan2>, typeof<Microsoft.Quantum.Canon.Native.ArcTan2>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ToDouble>, typeof<Microsoft.Quantum.Canon.Native.ToDouble>);

    

    let estimateEnergies =
        let qsim = QuantumSimulator()
        registerCallables qsim

        let H2EstimateEnergyRPE = qsim.Get(typeof<H2EstimateEnergyRPE>) :?> H2EstimateEnergyRPE
        let H2BondLengths = qsim.Get(typeof<H2BondLengths>) :?> H2BondLengths

        let estAtBondLength idx =
            [0..3]
            |> Seq.map (fun idxRep ->
                    H2EstimateEnergyRPE.Body.Invoke (struct (idx, int64 8, float 1.0))
                )
            |> Seq.min

        let bondLengths =
            H2BondLengths.Body.Invoke QVoid.Instance

        let data =
            asyncSeq {
                for idxBond in [0..53] do
                yield bondLengths.[idxBond], (int64 >> estAtBondLength) idxBond
            }
            
        data
        
    [<STAThread>]
    [<EntryPoint>]
    let main argv = 
        let window = Window()
        window.Loaded.Add <| fun eventArgs ->
            
            let chart =
                estimateEnergies
                |> AsyncSeq.toObservable
                |> LiveChart.LineIncremental
                |> fun chart ->
                    chart
                        .WithXAxis(Title = "Bond Length")
                        .WithYAxis(Title = "Est. GS Energy (Ha)")
                |> ChartControl
            
            let integrationHost = new Forms.Integration.WindowsFormsHost(Child = chart)
            
            window.Content <- integrationHost

        Application().Run(window)
            |> ignore
        0

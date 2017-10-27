namespace Microsoft.Quantum.Examples

open System
open System.Windows
open System.Windows.Controls

open Microsoft.Quantum.Simulation.Simulators
open Microsoft.Quantum.Canon
open Microsoft.Quantum.Simulation.Core

open Microsoft.FSharp.Core
open Microsoft.FSharp.Collections

open FSharp.Charting
open FSharp.Control
open FSharp.Charting.ChartTypes

module H2PlottingDemo =

    let registerCallables (qsim : AbstractFactory<Operation>) =        
        qsim.Register(typeof<Microsoft.Quantum.Canon.Ceiling>, typeof<Microsoft.Quantum.Canon.Native.Ceiling>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ArcTan2>, typeof<Microsoft.Quantum.Canon.Native.ArcTan2>)
        qsim.Register(typeof<Microsoft.Quantum.Canon.ToDouble>, typeof<Microsoft.Quantum.Canon.Native.ToDouble>);

    let qsim = QuantumSimulator()
    registerCallables qsim
    
    let H2EstimateEnergyRPE = qsim.Get(typeof<H2EstimateEnergyRPE>) :?> H2EstimateEnergyRPE
    let H2BondLengths = qsim.Get(typeof<H2BondLengths>) :?> H2BondLengths
    
    let bondLengths =
        H2BondLengths.Body.Invoke QVoid.Instance        
    
    let estAtBondLength idx =
        [0..3]
        |> Seq.map (fun idxRep ->
                H2EstimateEnergyRPE.Body.Invoke (struct (idx, int64 6, float 1))
            )
        |> Seq.min

    let estimateEnergies =
        asyncSeq {
            for idxBond in [0..53] do
            yield bondLengths.[idxBond], (int64 >> estAtBondLength) idxBond
        }
        
        
    [<STAThread>]
    [<EntryPoint>]
    let main argv = 
        let window = Window()
        window.Loaded.Add <| fun eventArgs ->
            
            let estEnergyChart =
                estimateEnergies
                |> AsyncSeq.toObservable
                |> fun data ->
                    LiveChart
                        .LineIncremental(data, Name="Estimated")
                        .WithMarkers(Size=12, Style=MarkerStyle.Circle)

            let theoryChart =
                [
                    0.14421
                    -0.323939
                    -0.612975
                    -0.80051
                    -0.92526
                    -1.00901
                    -1.06539
                    -1.10233
                    -1.12559
                    -1.13894
                    -1.14496
                    -1.1456
                    -1.14268
                    -1.13663
                    -1.12856
                    -1.1193
                    -1.10892
                    -1.09802
                    -1.08684
                    -1.07537
                    -1.06424
                    -1.05344
                    -1.043
                    -1.03293
                    -1.02358
                    -1.01482
                    -1.00665
                    -0.999025
                    -0.992226
                    -0.985805
                    -0.980147
                    -0.975156
                    -0.970807
                    -0.966831
                    -0.963298
                    -0.960356
                    -0.957615
                    -0.95529
                    -0.953451
                    -0.951604
                    -0.950183
                    -0.949016
                    -0.947872
                    -0.946982
                    -0.946219
                    -0.945464
                    -0.944887
                    -0.944566
                    -0.94415
                    -0.943861
                    -0.943664
                    -0.943238
                    -0.943172
                    -0.942973
                ]
                |> Seq.zip bondLengths
                |> fun data ->
                    Chart
                        .Line(data, Name="Theory")
                        

            let chart =
                Chart.Combine([estEnergyChart; theoryChart])
                |> fun chart ->
                    chart
                        .WithXAxis(Title = "BOND LENGTH", TitleFontName="Segoe UI Semibold", TitleFontSize = 24.0)
                        .WithYAxis(Title = "ENERGY", TitleFontName="Segoe UI Semibold", TitleFontSize = 24.0)
                        .WithLegend(FontSize = 24.0, FontName = "Segoe UI")
                        .WithTitle("H₂", FontName="Segoe UI", FontSize = 24.0)
                |> ChartControl
            
            let integrationHost = new Forms.Integration.WindowsFormsHost(Child = chart)
            
            window.Content <- integrationHost

        Application().Run(window)
            |> ignore
        0

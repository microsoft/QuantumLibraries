module Plots

open System
open System.IO
open System.Drawing
open System.Windows.Forms
open System.Windows.Forms.DataVisualization.Charting
open System.Windows.Forms.DataVisualization.Charting.Utilities

/////////////////////// Plotting utils from Dave Wecker ///////////////////////////////////

type pltT() =
    inherit Chart()
    
    static let axisFont = new Font("Calibri", 16.f, FontStyle.Bold)
    static let axisColor = Color.Blue
    static let titleFont = new Font("Calibri", 20.f, FontStyle.Bold)
    static let titleColor = Color.DarkBlue
    static let legendFont = new Font("Calibri", 12.f, FontStyle.Bold)
    
//    static let axisFont = new Font("Times New Roman", 16.f, FontStyle.Bold)
//    static let axisColor = Color.Blue
//    static let titleFont = new Font("Times New Roman", 20.f, FontStyle.Bold)
//    static let titleColor = Color.DarkBlue
//    static let legendFont = new Font("Times New Roman", 12.f, FontStyle.Bold)
//    
    static let copyChart (c1:pltT) =
        let copyTitle (t1:Title) =
            let t = new Title()
            t.Text <- t1.Text
            t.Font <- t1.Font
            t.ForeColor <- t1.ForeColor
            t.BorderColor <- t1.BorderColor
            t.BackColor <- t1.BackColor
            t.Alignment <- t1.Alignment
            t
        let c = new pltT()
        c1.Titles |> Seq.iter (fun t -> c.Titles.Add(copyTitle t))
        c
        
    static let copyArea (a1:ChartArea) =
        let a = new ChartArea()
        a.AxisX.MajorGrid.Enabled <- a1.AxisX.MajorGrid.Enabled
        a.AxisY.MajorGrid.Enabled <- a1.AxisY.MajorGrid.Enabled
        a.AxisX.Title <- a1.AxisX.Title
        a.AxisX.TitleFont <- a1.AxisX.TitleFont 
        a.AxisX.TitleForeColor <- a1.AxisX.TitleForeColor       
        a.AxisY.Title <- a1.AxisY.Title
        a.AxisY.TitleFont <- a1.AxisY.TitleFont
        a.AxisY.TitleForeColor <- a1.AxisY.TitleForeColor
        a        

    static let rec copySeries (s1:Series) (c:pltT) a: Series =
        let findLastNotEmpty (seriesCol:SeriesCollection) =
            Array.FindLast(seriesCol |> Seq.toArray, fun s -> not(s.Points.Count = 0) && not(s.ChartType = SeriesChartType.BoxPlot))    
        if s1.Points.Count = 0 && not(s1.ChartType = SeriesChartType.BoxPlot) then
            let s = copySeries (c.Series |> findLastNotEmpty) c a
            s.ChartType <- s1.ChartType
            s.ChartArea <- a
            s
        else
            let copyPoint (p1:DataPoint) =
                let p = new DataPoint(p1.XValue, p1.YValues |> Array.copy)
                p
            let s = new Series()
            s1.Points |> Seq.iter (fun p -> s.Points.Add(copyPoint p))
            s.ChartType <- s1.ChartType
            if s.ChartType = SeriesChartType.BoxPlot then
                let lastSeries = c.Series |> findLastNotEmpty
                s.["BoxPlotSeries"] <- (lastSeries).Name
                lastSeries.Enabled <- false
            s.IsValueShownAsLabel <- s1.IsValueShownAsLabel
            s.MarkerSize <- s1.MarkerSize
            s.MarkerStyle <- s1.MarkerStyle
            s.Color <- s1.Color
            s.Enabled <- s1.Enabled
            s.["DrawingStyle"] <- s1.["DrawingStyle"]
            s.ChartArea <- a
            s
       
    static let addAreaAndSeries (c:pltT) oldArea (seriesCol:SeriesCollection) =
        let newArea = copyArea oldArea
        c.ChartAreas.Add(newArea)
        let relatedSeries = seriesCol |> Seq.filter (fun s -> s.ChartArea = oldArea.Name)
        relatedSeries |> Seq.iter (fun s -> c.Series.Add(copySeries s c newArea.Name))         
        
    static let Create (chartType, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title) =
        let c = new pltT()
        let a = new ChartArea()
        let s = new Series()
        s.ChartType <- chartType
        c.ChartAreas.Add(a)
        c.Series.Add(s)
        match x, y with
            | Some(x), None     -> failwith "You cannot pass only x to a chart drawing function"
            | Some(x), Some(y)  -> s.Points.DataBindXY(x, [|y|])
            | None, Some(y)     -> s.Points.DataBindY([|y|])
            | None, None        -> ()
        s.IsValueShownAsLabel <- defaultArg isValueShownAsLabel s.IsValueShownAsLabel
        s.MarkerSize <- defaultArg markerSize s.MarkerSize
        s.MarkerStyle <- defaultArg markerStyle s.MarkerStyle
        s.Color <- defaultArg color s.Color
        
        a.AxisX.MajorGrid.Enabled <- false
        a.AxisY.MajorGrid.Enabled <- false
        
        match xname with
        | Some(xname) ->
            a.AxisX.Title <- xname
            a.AxisX.TitleFont <- axisFont
            a.AxisX.TitleForeColor <- axisColor
        | _ -> ()
        match yname with
        | Some(yname) ->
            a.AxisY.Title <- yname
            a.AxisY.TitleFont <- axisFont
            a.AxisY.TitleForeColor <- axisColor
        | _ -> ()
        match seriesName with
        | Some(seriesName) -> s.Name <- seriesName
        | _ -> ()
        match title with
        | Some(title) ->
            let t = c.Titles.Add(title: string)
            t.Font <- titleFont
            t.ForeColor <- titleColor
        | _ -> ()
        c
    
    member c.Title(title:string) =
        let t = c.Titles.Add(title: string)
        t.Font <- titleFont
        t.ForeColor <- titleColor

    static member (+) (c1:pltT, c2:pltT) =    
        let c = copyChart(c1)   
        c1.ChartAreas |> Seq.iter (fun a -> addAreaAndSeries c a c1.Series)
        let lastArea = c.ChartAreas |> Seq.nth ((c.ChartAreas |> Seq.length) - 1)
        c2.Series |> Seq.iter(fun s -> c.Series.Add(copySeries s c lastArea.Name))
        let l = c.Legends.Add("")
        l.Font <- legendFont
        c
       
    static member (++) (c1:pltT, c2:pltT) =
        let c = copyChart(c1)   
        c1.ChartAreas |> Seq.iter (fun a -> addAreaAndSeries c a c1.Series)
        let lastArea = c.ChartAreas |> Seq.nth ((c.ChartAreas |> Seq.length) - 1)
        addAreaAndSeries c c2.ChartAreas.[0] c2.Series
        let firstArea = c.ChartAreas |> Seq.nth ((c.ChartAreas |> Seq.length) - 1)
        c2.ChartAreas |> Seq.skip 1 |> Seq.iter (fun a -> addAreaAndSeries c a c2.Series)
        c    
        
    // Here are the interesting public static construction methods
    // There is one method for every type of chart that we support
    
    static member scatter (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title) = 
        Create (SeriesChartType.Point, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
    static member line (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title) = 
        Create (SeriesChartType.Line, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
    static member spline (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title) = 
        Create (SeriesChartType.Spline, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
    static member stepline (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title) = 
        Create (SeriesChartType.StepLine, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
    static member bar (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title, ?drawingStyle) = 
        let c = Create (SeriesChartType.Bar, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
        c.Series.[0].["DrawingStyle"] <- defaultArg drawingStyle (c.Series.[0].["DrawingStyle"])
        c
        
    // This is the F# wrapper for the .NET Bubble Chart.
        
    static member bubble (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title, ?drawingStyle) = 
        let c = Create (SeriesChartType.Bubble, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
        c.Series.[0].["DrawingStyle"] <- defaultArg drawingStyle (c.Series.[0].["DrawingStyle"])
        c
    static member column (?y,?x, ?isValueShownAsLabel, ?markerSize, ?markerStyle, ?color, ?xname, ?yname, ?seriesName, ?title, ?drawingStyle) = 
        let c = Create (SeriesChartType.Column, x, y, isValueShownAsLabel, markerSize, markerStyle, color, xname, yname, seriesName, title)
        c.Series.[0].["DrawingStyle"] <- defaultArg drawingStyle (c.Series.[0].["DrawingStyle"])
        c
    static member boxplot (?y, ?color, ?xname, ?yname, ?seriesName, ?title, ?whiskerPercentile, ?percentile, ?showAverage, ?showMedian, ?showUnusualValues) = 
        let c = Create (SeriesChartType.Spline, None, y, None, None, None, color, xname, yname, seriesName, title)
        c.Series.[0].Enabled <- false
        let s = new Series()
        c.Series.Add(s)
        s.ChartType <- SeriesChartType.BoxPlot
        s.["BoxPlotSeries"] <- c.Series.[0].Name
        whiskerPercentile |> Option.iter (fun x -> s.["BoxPlotWhiskerPercentile"] <- string x) 
        percentile |> Option.iter (fun x -> s.["BoxPlotPercentile"] <- string x) 
        showAverage |> Option.iter (fun x -> s.["BoxPlotShowAverage"] <- string x) 
        showMedian |> Option.iter (fun x -> s.["BoxPlotShowMedian"] <- string x)
        s.["BoxPlotShowUnusualValues"] <- "true" 
        showUnusualValues |> Option.iter (fun x -> s.["BoxPlotShowUnusualValues"] <- string x)
        c
    static member legend(?title) =
        let l = new Legend()
        title |> Option.iter (fun t -> l.Title <- string title)
        l 


    // Display the plot and optionally don't wait
    member c.display(?wait,?width,?height,?offset) =
        let wait    = defaultArg wait true
        let width   = defaultArg width 800
        let height  = defaultArg height 600
        let offset  = defaultArg offset 100
    
        let copy () =
            let stream = new MemoryStream()
            c.SaveImage(stream, Imaging.ImageFormat.Bmp)
            let bmp = new Bitmap(stream) 
            Clipboard.SetDataObject(bmp) 
       
        c.KeyDown.Add(fun e -> if e.KeyCode = Keys.A then copy () )
        // c.KeyDown.Add(fun e -> if e.Control = true && e.KeyCode = Keys.C then copy ())
        let pressToCopy = "(press CTRL+C to copy)"     
        let name = if c.Titles.Count = 0 then sprintf "%s %s " "pltT" pressToCopy else sprintf "%s %s " c.Titles.[0].Text  pressToCopy
        let f = new Form(Text = name, Size = new Size(width,height), 
                            Location = new Point(offset,offset), 
                            StartPosition=FormStartPosition.Manual,
                            TopMost=false)
        c.Dock <- DockStyle.Fill
        f.Controls.Add(c)
        if wait then f.ShowDialog() |> ignore 
        else 
            f.Show()
            Application.DoEvents()

    // ---------------------------------------------------------------------------
    // Sample tests for the API above.
    // ---------------------------------------------------------------------------

    static member debug() =
        let x = [1.;2.5;3.1;4.;4.8;6.0;7.5;8.;9.1;15.]
        let y = [1.6;2.1;1.4;4.;2.3;1.9;2.4;1.4;5.;2.9]

        pltT.scatter(x, y).display()
        pltT.scatter(x = x, y = y, markerSize = 10, markerStyle = MarkerStyle.Diamond,
            xname = "Players", yname = "Ratings", title = "Players' Ratings").display()
        pltT.line(y = y, markerSize = 10, markerStyle = MarkerStyle.Diamond, xname = "Players", yname = "Ratings", title = "Players' Ratings", isValueShownAsLabel = true,
            color = Color.Red).display()
        pltT.spline(x = x, y = y, markerSize = 10, markerStyle = MarkerStyle.Diamond, xname = "Players", yname = "Ratings",
            title = "Players' Ratings", isValueShownAsLabel = true, color = Color.Red).display()
        pltT.stepline(x = x, y = y, markerSize = 10, markerStyle = MarkerStyle.Diamond, xname = "Players", yname = "Ratings",
            title = "Players' Ratings", isValueShownAsLabel = true, color = Color.Red).display()
        pltT.bar(y = y, xname = "Players", yname = "Ratings", title = "Players' Ratings", isValueShownAsLabel = true,
            drawingStyle = "Emboss").display()
        pltT.column(y = y, xname = "Players", yname = "Ratings", title = "Players' Ratings",
            isValueShownAsLabel = true, drawingStyle = "Cylinder").display()
        pltT.boxplot(y = y, xname = "Players", yname = "Ratings", title = "Players' Ratings", color = Color.Blue, whiskerPercentile = 5, percentile = 30,
            showAverage = false, showMedian = false, showUnusualValues = true).display()
        (pltT.bubble(y = y, xname = "Players", yname = "Ratings", title = "Players' Ratings", color = Color.Blue)).display()
        (pltT.scatter(y) + pltT.line(y)).display()
        (pltT.scatter(y, markerSize = 10) + pltT.column() + pltT.boxplot() + pltT.line()  + pltT.column(x) + pltT.boxplot()).display()
        (pltT.scatter(y, markerSize = 10) + pltT.column() + (pltT.line(x)  + pltT.column()) + pltT.scatter(markerSize = 20)).display()

        // This sample will plot a combination of graphs
        // We plot a line and column graph justiposed on top of each other,
        // combined with a line graph at the bottom.
        // Then we pipe it to the display

        let h = [1.;2.5;3.1;4.;4.8;6.0;7.5;8.;9.1;15.]
        let w = h |> List.map (fun h -> h * 1.2)

        (pltT.line(h) + pltT.column() ++ pltT.line(w) ++ pltT.bubble(w)).display()



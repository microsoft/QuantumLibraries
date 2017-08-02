namespace Microsoft.Research.Liquid

open System
open System.Text

open Microsoft.Research.Liquid
open Util
open Operations
open Shor
open Tests

///////////////////////////// Allow installation as a service /////////////////////

/// <summary>
/// This allows us to be run as a service (install with InstallUtil.exe from .Net Framework)
/// Start with: sc start "LiquidService"
/// Stop  with: sc stop  "LiquidService"
/// </summary>
[<ComponentModel.RunInstaller(true)>]
type LiquidInstaller() =
    inherit Configuration.Install.Installer()
    do LiquidWindowsService.Install base.Installers

module Tests =
    [<LQD>]
    let __doDebug() =
        let runProc (cmd:string) (args:string) =
            show "=============== RUNNING: %s %s =========================" cmd args
            let pi  = Diagnostics.ProcessStartInfo(cmd,args)
            pi.RedirectStandardOutput  <- true
            pi.UseShellExecute         <- false
            let p                       = Diagnostics.Process.Start pi
            show "%s" (p.StandardOutput.ReadToEnd())
            show ""
            p.WaitForExit() 

#if FALSE
        let binPath = "\"" + Environment.CurrentDirectory + "\\Liquid.exe\""
        runProc "sc.exe" (@"\\localhost create LiquidService start= auto binPath= " + binPath)
        Threading.Thread.Sleep 2000
        runProc "sc.exe" @"\\localhost start LiquidService"
        Threading.Thread.Sleep 2000
        runProc "sc.exe" @"\\localhost query LiquidService"
#endif
        try
            let url = sprintf "tcp://localhost:%d/LiquidService" Parser.Port
            let svc = Activator.GetObject(typeof<ILiquidService>, url) :?> ILiquidService
            show "Sending command"
            svc.SendCommand "-e __doTeleport"
            show "Getting output"
            let rec loop() =
                Threading.Thread.Sleep 1000
                let running,str     = svc.GetOutput()
                show "%s" str
                if running then loop()
            loop()
        with e ->
            show "Error talking to service: %s" e.Message
            show "%s" e.StackTrace
#if FALSE
        runProc "sc.exe" @"\\localhost stop LiquidService"
        Threading.Thread.Sleep 2000
        runProc "sc.exe" @"\\localhost query LiquidService"
#endif

module App =

    /// <summary>
    /// Main entry point for application (both interactive and service)
    /// </summary>
    [<EntryPoint>]
    let Main _ =
        // TCP service port to use for all sends and receives
        Parser.Port        <- 28028

        // Skip the first argument (the program name)
        let args    = 
            Environment.GetCommandLineArgs()
            |> Seq.skip 1
            |> Seq.toList

        // Need to get rid of any shell escaping (can't use '^' in an argument)
        let args    = List.map (fun (arg:string) -> arg.Replace("^","")) args

        // Figures out if we need to fake being interactive for HPC
        let p       = Parser(args)
        let las     = p.CommandArgs()

        // Are we supposed to run as a service?
        if las.svcInst then
            LiquidWindowsService.Main Parser.Port
        else
            // Run the command line
            p.CommandRun las


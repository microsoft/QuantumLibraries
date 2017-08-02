module Microsoft.Research.Liquid.SignalProcessing

open System
open Util
open Operations
open Montgomery


///////////////////////////////////////////////////////////////////////////
//
// Quantum discrete cosine transforms
//
///////////////////////////////////////////////////////////////////////////

/// <summary>
/// Controlled rotation gate
/// </summary>
/// <param name="k">Rotation by 2^k</param>
/// <param name="qs">Qubits ([0]=Control [1]=rotated)</param>
let CR (k:int) (qs:Qubits) =
    let gate (k:int) (qs:Qubits)    =
        Gate.Build("CR_" + k.ToString() ,fun () ->
            new Gate(
                Qubits  = qs.Length,
                Name    = "CR",
                Help    = "Controlled R gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{R%d}" k,
                Op      = WrapOp (fun (qs:Qubits) -> Cgate (R k) qs)
        ))
    (gate k qs).Run qs

/// <summary>
/// Take QFT
/// </summary>
/// <param name="qs">Qubits to take QFT of</param>
let QFTdave (qs:Qubits) =
    let gate (qs:Qubits)    =
        let op (qs:Qubits) =
            let n   = qs.Length
            
            for i in 0..(n/2-1) do 
                CNOT [qs.[i]; qs.[n-i-1]]
                CNOT [qs.[n-i-1]; qs.[i]]
                CNOT [qs.[i]; qs.[n-i-1]]
            for aIdx in n-1..-1..0 do
                let a   = qs.[aIdx]
                H [a]
                for k in 2..aIdx+1 do
                    let c   = qs.[aIdx-(k-1)]
                    CR k [c;a]
            
        Gate.Build("QFT_" + qs.Length.ToString() ,fun () ->
            new Gate(
                Qubits  = qs.Length,
                Name    = "QFT",
                Help    = "QFT",
                Draw    = sprintf "\\multigate{#%d}{QFT}" (qs.Length-1),
                Op      = WrapOp op
        ))
    (gate qs).Run qs

/// <summary>
/// Take QFT
/// </summary>
/// <param name="qs">Qubits to take QFT of</param>
let QFT (qs:Qubits) =
    let gate (qs:Qubits)    =
        let op (qs:Qubits) =
            let n   = qs.Length
            for i in 0..(n/2-1) do 
                CNOT [qs.[i]; qs.[n-i-1]]
                CNOT [qs.[n-i-1]; qs.[i]]
                CNOT [qs.[i]; qs.[n-i-1]]
            for i in 0..(n-1) do
                for j in 0..(i-1) do
                    CR (i-j+1) [qs.[i];qs.[j]]
                H [qs.[i]]
                
        Gate.Build("QFT_" + qs.Length.ToString() ,fun () ->
            new Gate(
                Qubits  = qs.Length,
                Name    = "QFT",
                Help    = "QFT",
                Draw    = sprintf "\\multigate{#%d}{QFT}" (qs.Length-1),
                Op      = WrapOp op
        ))
    (gate qs).Run qs

/// <summary>
/// B gate
/// </summary>
let B (qs:Qubits) =
        let gate =
            let nam     = "B"
            new Gate(
                Name    = nam,
                Help    = sprintf "B gate",
                Mat     = (
                    let sq = (1.0 / (Math.Sqrt 2.0))
                    CSMat(2,[0,0,sq,0.;0,1,0.,sq;1,0,sq,0.;1,1,0.,-sq])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// Badj gate
/// </summary>
let Badj (qs:Qubits) =
        let gate =
            let nam     = "Badj"
            new Gate(
                Name    = nam,
                Help    = sprintf "Badj gate",
                Mat     = (
                    let sq = (1.0 / (Math.Sqrt 2.0))
                    CSMat(2,[0,0,sq,0.;0,1,sq,0.;1,0,0.,-sq;1,1,0.,sq])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// Bbar gate
/// </summary>
let Bbar (qs:Qubits) =
        let gate =
            let nam     = "Bbar"
            new Gate(
                Name    = nam,
                Help    = sprintf "Bbar gate",
                Mat     = (
                    let sq = (1.0 / (Math.Sqrt 2.0))
                    CSMat(2,[0,0,sq,0.;0,1,0.,-sq;1,0,sq,0.;1,1,0.,sq])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// controlled B gate
/// </summary>
let cB (qs:Qubits) =
        let gate =
            new Gate(
                Qubits = qs.Length,
                Name    = "cB",
                Help    = sprintf "controlled B gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{B}",
                Op      = WrapOp (fun (qs:Qubits) -> Cgate B qs)
            )
        gate.Run qs

/// <summary>
/// controlled Badj gate
/// </summary>
let cBadj (qs:Qubits) =
        let gate =
            new Gate(
                Qubits = qs.Length,
                Name    = "cBadj",
                Help    = sprintf "controlled Badj gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{Badj}",
                Op      = WrapOp (fun (qs:Qubits) -> Cgate Badj qs)
            )
        gate.Run qs

/// <summary>
/// J gate
/// </summary>
let J (qs:Qubits) =
        let gate =
            let nam     = "J"
            new Gate(
                Name    = nam,
                Help    = sprintf "J gate",
                Mat     = (
                    let sq = (1.0 / (Math.Sqrt 2.0))
                    CSMat(2,[0,0,sq,0.;0,1,0.,-sq;1,0,0.,-sq;1,1,sq,0.])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// controlled J gate
/// </summary>
let cJ (qs:Qubits) =
        let gate =
            new Gate(
                Qubits = qs.Length,
                Name    = "cJ",
                Help    = sprintf "controlled J gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{J}",
                Op      = WrapOp (fun (qs:Qubits) -> Cgate J qs)
            )
        gate.Run qs

/// <summary>
/// MM gate
/// </summary>
let MM (k:int) (qs:Qubits) =
        let gate =
            let nam     = "MM"
            new Gate(
                Name    = nam,
                Help    = sprintf "MM gate",
                Mat     = (
                    let phi     = (2.0 * Math.PI)/(8.0 * (pown 2.0 k))
                    let phiR    = Math.Cos phi
                    let phiI    = Math.Sin phi
                    CSMat(2,[(0,0,phiR,phiI);(1,1,phiR,phiI)])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// C gate
/// </summary>
let C (k:int) (qs:Qubits) =
        let gate =
            let nam     = "C"
            new Gate(
                Name    = nam,
                Help    = sprintf "C gate",
                Mat     = (
                    let phi     = -(2.0*Math.PI)/(4.0 * (pown 2.0 k)) 
                    let phiR    = Math.Cos phi
                    let phiI    = Math.Sin phi
                    CSMat(2,[(0,0,1.,0.);(1,1,phiR,phiI)])),
                Draw    = "\\gate{" + nam + "}"
                )
        gate.Run qs

/// <summary>
/// K gate
/// </summary>
/// <param name="k"> order of root of unity</param>
/// <param name="i"> index of K gate</param>
let K (k:int) (i:int) (qs:Qubits) =
    let gate (k:int) (i:int) =
        Gate.Build("K_" + k.ToString() + "_" + i.ToString(),fun () ->
            new Gate(
                Name    = sprintf "K_%d_%d" k i,
                Help    = sprintf "(2pi/2^k)^(2^i): %d %d" k i,
                Mat     = (
                    let phi     = -(2.0*Math.PI)*(pown 2.0 i)/(4.0 * (pown 2.0 k))
                    let phiR    = Math.Cos phi
                    let phiI    = Math.Sin phi
                    CSMat(2,[(0,0,1.,0.);(1,1,phiR,phiI)])),
                    //CSMat(2,[(0,0,phiR,phiI);(1,1,1.,0.)])),
                Draw    = "\\gate{" + "K_" + (k.ToString()) + "_" + (i.ToString()) + "}"
               ))
    (gate k i).Run qs

/// <summary>
/// controlled K gate
/// </summary>
let cK (k:int) (i:int) (qs:Qubits) =
        let gate =
            new Gate(
                Qubits = qs.Length,
                Name    = "cK",
                Help    = sprintf "controlled K gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{K_" + (k.ToString()) + "_" + (i.ToString()) + "}",
                Op      = WrapOp (fun (qs:Qubits) -> Cgate (K k i) qs)
            )
        gate.Run qs

/// <summary>
/// L gate
/// </summary>
/// <param name="k"> order of root of unity</param>
/// <param name="i"> index of K gate</param>
let L (k:int) (i:int) (qs:Qubits) =
    let gate (k:int) (i:int) =
        Gate.Build("L_" + k.ToString() + "_" + i.ToString(),fun () ->
            new Gate(
                Name    = sprintf "L_%d_%d" k i,
                Help    = sprintf "(2pi/2^k)^(2^i): %d %d" k i,
                Mat     = (
                    let phi     = (2.0*Math.PI)*(pown 2.0 i)/(4.0 * (pown 2.0 k))
                    let phiR    = Math.Cos phi
                    let phiI    = Math.Sin phi
                    CSMat(2,[(0,0,1.,0.);(1,1,phiR,phiI)])),
                Draw    = "\\gate{" + "L_" + (k.ToString()) + "_" + (i.ToString()) + "}"
               ))
    (gate k i).Run qs
        
/// <summary>
/// controlled L gate
/// </summary>
let cL (k:int) (i:int) (qs:Qubits) =
        let gate =
            new Gate(
                Qubits = qs.Length,
                Name    = "cL",
                Help    = sprintf "controlled L gate",
                Draw    = sprintf "\\ctrl{#1}\\go[#1]\\gate{L_" + (k.ToString()) + "_" + (i.ToString()) + "}",
                Op      = WrapOp (fun (qs:Qubits) -> Cgate (L k i) qs)
            )
        gate.Run qs
        
let DCT_I (qs:Qubits) = 
    let n = qs.Length-1 // last qubit is the ancilla, transform operates on qubits 0..(n-2)
    
    let CtrlPerm (qs:Qubits) = 
        for i in (n-2)..(-1)..0 do 
            BuildMultiplyControlledNOT((slice qs [0..(i-1)]) @ [qs.[n-1]]) [qs.[i]] [qs.[n]]

    let CtrlPermAdj (qs:Qubits) = 
        for i in 0..(n-2) do 
            BuildMultiplyControlledNOT((slice qs [0..(i-1)]) @ [qs.[n-1]]) [qs.[i]] [qs.[n]]
        
    let CtrlBgate (qs:Qubits) = 
        for i in 0..(n-2) do 
            X [qs.[i]]
        BuildMultiplyControlledNOT (slice qs [0..(n-2)]) [qs.[n]] [qs.[n-1]]
        cB [qs.[n]; qs.[n-1]]
        BuildMultiplyControlledNOT (slice qs [0..(n-2)]) [qs.[n]] [qs.[n-1]]
        for i in 0..(n-2) do 
            X [qs.[i]]   
    
    let CtrlBadjgate (qs:Qubits) = 
        for i in 0..(n-2) do 
            X [qs.[i]]
        BuildMultiplyControlledNOT (slice qs [0..(n-2)]) [qs.[n]] [qs.[n-1]]
        cBadj [qs.[n]; qs.[n-1]]
        BuildMultiplyControlledNOT (slice qs [0..(n-2)]) [qs.[n]] [qs.[n-1]]
        for i in 0..(n-2) do 
            X [qs.[i]]   

    // now apply the actual functions to compute the DCT_II
    B [qs.[n-1]]
    CtrlBadjgate qs
    for i in 0..(n-2) do 
        CNOT [qs.[n-1]; qs.[i]]
    CtrlPerm qs
    QFT (slice qs [0..(n-1)])
    CtrlPermAdj qs
    for i in 0..(n-2) do 
        CNOT [qs.[n-1]; qs.[i]]
    CtrlBgate qs
    Badj [qs.[n-1]]

let DCT_IV (qs:Qubits) = 
    let n = qs.Length 
    
    // now apply the actual functions to compute the DCT_II
    Bbar [qs.[n-1]]
    for i in 0..(n-2) do 
        (cK (n-1) i) [qs.[n-1]; qs.[i]]
    X [qs.[n-1]] 
    for i in 0..(n-2) do 
        (cL (n-1) i) [qs.[n-1]; qs.[i]]
    X [qs.[n-1]]
    for i in 0..(n-2) do 
        CNOT [qs.[n-1]; qs.[i]]
    (C (n-1)) [qs.[n-1]]
    QFT qs
    (C (n-1)) [qs.[n-1]]
    for i in 0..(n-2) do 
        CNOT [qs.[n-1]; qs.[i]]
    for i in 0..(n-2) do 
        (cK (n-1) i) [qs.[n-1]; qs.[i]]
    X [qs.[n-1]] 
    for i in 0..(n-2) do 
        (cL (n-1) i) [qs.[n-1]; qs.[i]]
    X [qs.[n-1]]
    Badj [qs.[n-1]]
    (MM (n-1)) [qs.[n-1]]

[<LQD>]
/// <summary> Compute discrete cosine transforms based on Klappenecker/Roetteler quant-ph/0111038 </summary>
/// <param name="name">DCT type (I, II, III, or IV) </param>
/// <param name="n">number of qubits </param>
let __DCT (name:string) (n:int) =
    // for DCT_I need 2 ancillas: 1 for embedding DCT_(2^n) + DST_(2^n) -> DFT_(2^(n+1)) and 1 extra ancilla
    // for DCT_IV need 1 ancillas: 1 for embedding DCT_(2^n) + DST_(2^n) -> DFT_(2^(n+1))
    let numQubits = function 
        | "I" -> n+2 
        | "IV" -> n+1
        | _ -> failwith "unsupported DCT type"
    let m = numQubits name
    let k = Ket(m)
    let v = k.Single()
    let qs = k.Qubits
    
    let dctfun = 
        match name with 
        | "I" -> DCT_I
        | "IV" -> DCT_IV
        | _ -> failwith "unsupported DCT type"
    
    let circ = Circuit.Compile dctfun qs
    //circ.Dump()

    let BoolInt n x = [ for i in [0..(n-1)] do yield (int ((x >>> i) % 2)) ] 
    
    // dump the entire unitary into the logfile
    for i in 0..((pown 2 m)-1) do 
        let L = BoolInt m i 
        qs |> List.iteri (fun j x -> M [x]; x.ReAnimate Zero; if (L.[j]=1) then X [x];) 
        circ.Run qs
        k.Dump(showLogInd,0,false,true)
    
    show "Done"

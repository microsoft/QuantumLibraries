module Microsoft.Research.Liquid.UserSample

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions
open System.Drawing
open System.Windows.Forms
open System.Windows.Forms.DataVisualization.Charting
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open Shor               // factoring using semi-classical qft
open Plots              // nifty .net based plot lib written by Dave


///////////////////////////////////////////////////////////////////////////
//
// Tools for slicing and to convert between bigints and list of bools
//
///////////////////////////////////////////////////////////////////////////

let BoolInt n x = 
    [ for i in [0..(n-1)] do yield (int ((x >>> i) % 2I)) ] 

let IntBool n (x: int array) = 
    List.fold (+) 0I <| List.map (fun i -> (1I <<< i) * (bigint x.[i])) [0..x.Length-1]

let PrepBool (qs:Qubits) offset step ls = 
    List.mapi (fun i x -> if x <> 0 then X !!(qs,(offset+i*step))) ls 

let slice (qs:Qubits) ls = List.map (fun n-> qs.[n]) ls

[<LQD>]
let __QuickTest () = 
    // add stuff here for quick test of F# functionality
    for i in 10..(-1)..1 do 
        printf "%A\n" i

///////////////////////////////////////////////////////////////////////////
//
// Tools for circuit simulation of Toffoli networks and pretty printing 
//
///////////////////////////////////////////////////////////////////////////

type TofGate = 
    | MyNOT of int 
    | MyCNOT of int * int 
    | MyTOFF of int * int * int 

let rec MyCircuitDump (file:string) (c:Circuit) = 
     // pretty prints circuit to a file. Code patterned after member function c.Dump() in Circuits.fs  
    let MyPrinter (file:string) (s:string) = 
        let a = [| for m in Regex.Matches(s, "[a-zA-Z0-9\- ]+") do yield (string m) |]
        File.AppendAllLines(file, a)

    // (incomplete) match over the discriminated union of Circuit type. Todo: implement missing cases later
    match c with
    | Empty -> MyPrinter file "EOF"
    | Seq cs       ->
        for c in cs do MyCircuitDump file c
    | Par cs       ->
        for c in cs do MyCircuitDump file c
    | Apply(g,ws)   ->
        match g.Name with 
        | "X" -> MyPrinter file (String.concat "" [| "tof"; " "; (ws.[0].ToString()) |])
        | "CNOT" -> MyPrinter file (String.concat "" 
                        [| "tof"; " "; (ws.[0].ToString()); " "; (ws.[1].ToString()); |])
        | "CCNOT" -> MyPrinter file (String.concat ""
                        [| "tof"; " "; (ws.[0].ToString()); " "; (ws.[1].ToString()); " "; (ws.[2].ToString()) |])
        | _ -> MyPrinter file "unrecognized gate" 
    | _ ->  MyPrinter file "unrecognized construct or gate"

let rec MyCircuitExport (c:Circuit) = 
    // converts circuit to list of NOT, CNOT, and Toffoli gates. 
    match c with
    | Empty -> []
    | Seq cs       -> List.concat (List.map MyCircuitExport cs)
    | Par cs       -> List.concat (List.map MyCircuitExport cs)
    | Apply(g,ws)   ->
        match g.Name with 
        | "X" -> [ MyNOT ws.[0] ]
        | "CNOT" | "CNOT'" -> [ MyCNOT (ws.[0], ws.[1]) ]
        | "CCNOT" | "CCNOT'" -> [ MyTOFF (ws.[0], ws.[1], ws.[2]) ]
        | _ -> failwithf "Expo: unrecognized gate %s" g.Name 
    | _ ->  
        failwith "Expo: unrecognized construct or gate"

let MyCircuitSimulate c (b:int []) = 
    // simulates a Toffoli network by applying it to an input bit pattern 
    let rec MyApply a (b:int []) = 
        match a with 
        | [] -> b
        | (MyNOT i)::gs -> 
            MyApply gs (Array.mapi (fun ind x -> if ind=i then (-x+1) else x) b)
        | (MyCNOT (i, j))::gs -> 
            MyApply gs (Array.mapi (fun ind x -> if ind=j then (x+b.[i]-2*x*b.[i]) else x) b)
        | (MyTOFF (i, j, k))::gs -> 
            MyApply gs (Array.mapi (fun ind x -> if ind=k then (x+(b.[i]*b.[j])-2*x*b.[i]*b.[j]) else x) b)         
    MyApply c b

///////////////////////////////////////////////////////////////////////////
//
// Decomposing multiply controlled NOT gates with dirty ancillas
//
///////////////////////////////////////////////////////////////////////////

let rec BuildMultiplyControlledNOT (cs:Qubits) (t:Qubits) (a:Qubits) = 
    let n = cs.Length
    let slice (qs:Qubits) ls = List.map (fun n-> qs.[n]) ls
    match n with 
    | 0 -> X [t.[0]]
    | 1 -> CNOT [cs.[0]; t.[0]] 
    | 2 -> CCNOT [cs.[0]; cs.[1]; t.[0]]
    | n -> BuildMultiplyControlledNOT (slice cs [0..n-2]) a [cs.[n-1]]        
           CCNOT [cs.[n-1]; a.[0]; t.[0]]
           BuildMultiplyControlledNOT (slice cs [0..n-2]) a [cs.[n-1]]        
           CCNOT [cs.[n-1]; a.[0]; t.[0]]
       
// MultiplyControlledNOT uses the method given in [Barenco et al, quant-ph/9503016] 
// to implement an n-fold controlled NOT using O(n) Toffoli gates, provided at least 
// one qubit is available as an ancilla. This qubit, which can be dirty is given 
// as value <a>, the control qubits are given as list <cs> and the target is <t>
let MultiplyControlledNOT (qs:Qubits) = 
    let n = qs.Length-2
    let cs = slice qs [0..(n-1)]
    let t = [ qs.[n] ]
    let a = [ qs.[n+1] ]
    BuildMultiplyControlledNOT cs t a


///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: integer addition
//
///////////////////////////////////////////////////////////////////////////

/// Cuccaro adder (inplace adder described in [Cuccaro et al, quant-ph/0410184])

let BuildCuccaroAdder (xs:Qubits) (ys:Qubits) (z:Qubits) (a:Qubits) = 
    if (xs.Length <> ys.Length) || (z.Length <> 1) then 
        failwith "Error: Cuccaro adder requires length(x) = length(y) and 1 qubit for carry"
    if xs.Length < 4 then 
        failwith "Error: Cuccaro adder requires n >= 4"
    let n = xs.Length

    for i in 1..(n-1) do 
        CNOT [xs.[i]; ys.[i]]
    CNOT [xs.[1]; a.[0]]
    CCNOT [xs.[0]; ys.[0]; a.[0]]
    CNOT [xs.[2]; xs.[1]]
    CCNOT [a.[0]; ys.[1]; xs.[1]]
    CNOT [xs.[3]; xs.[2]]
    for i in 2..(n-3) do 
        CCNOT [xs.[i-1]; ys.[i]; xs.[i]]
        CNOT [xs.[i+2]; xs.[i+1]]
    CCNOT [xs.[n-3]; ys.[n-2]; xs.[n-2]]
    CNOT [xs.[n-1]; z.[0]]
    CCNOT [xs.[n-2]; ys.[n-1]; z.[0]]
    for i in 1..(n-2) do
        X [ys.[i]]
    CNOT [a.[0]; ys.[1]] 
    for i in 2..(n-1) do 
        CNOT [xs.[i-1]; ys.[i]]
    CCNOT [xs.[n-3]; ys.[n-2]; xs.[n-2]]
    for i in (n-3)..(-1)..2 do 
        CCNOT [xs.[i-1]; ys.[i]; xs.[i]]
        CNOT [xs.[i+2]; xs.[i+1]]
        X [ys.[i+1]]
    CCNOT [a.[0]; ys.[1]; xs.[1]]
    CNOT [xs.[3]; xs.[2]]
    X [ys.[2]]
    CCNOT [xs.[0]; ys.[0]; a.[0]]
    CNOT [xs.[2]; xs.[1]]
    X [ys.[1]]
    CNOT [xs.[1]; a.[0]]
    for i in 0..(n-1) do
        CNOT [xs.[i]; ys.[i]]

let CuccaroAdder (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ys = slice qs [n..2*n-1]
    let z  = slice qs [2*n]
    let a  = slice qs [2*n+1]
    BuildCuccaroAdder xs ys z a

/// Simple ripple adder (adder described e.g. in [Draper, quant-ph/0008033])

let BuildRippleAdder (xs:Qubits) (ys:Qubits) (z:Qubits) (cs:Qubits)  = 
    if (xs.Length <> ys.Length) || (z.Length <> 1) then 
        failwith "Error: Ripple adder requires length(x) = length(y) and 1 qubit for carry"
    if xs.Length = 0 then 
        failwith "Error: Ripple adder requires n > 0"
    let n = xs.Length 
    
    let Carry (qs:Qubits) =
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        CNOT [qs.[1]; qs.[2]]
        CCNOT [qs.[0]; qs.[2]; qs.[3]]

    let CarryInv (qs:Qubits) =
        CCNOT [qs.[0]; qs.[2]; qs.[3]]
        CNOT [qs.[1]; qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
    
    let Sum (qs:Qubits) =
        CNOT [qs.[1]; qs.[2]]
        CNOT [qs.[0]; qs.[2]]
        
    for i in 0..(n-2) do 
        Carry [cs.[i]; xs.[i]; ys.[i]; cs.[i+1]]
    Carry [cs.[n-1]; xs.[n-1]; ys.[n-1]; z.[0]]
    CNOT [xs.[n-1]; ys.[n-1]]
    Sum [cs.[n-1]; xs.[n-1]; ys.[n-1]]
    for i in (n-2)..(-1)..0 do 
        CarryInv [cs.[i]; xs.[i]; ys.[i]; cs.[i+1]]
        Sum [cs.[i]; xs.[i]; ys.[i]]
    
let RippleAdder (qs:Qubits) = 
    let n = (qs.Length-1)/3
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let z  = slice qs [2*n]
    let cs = slice qs [2*n+1..3*n]
    BuildRippleAdder xs ys z cs

/// Takahashi adder (inplace adder described e.g. in [Takahasi et al, arxiv:0910.2530])

let BuildTakahashiAdder (xs:Qubits) (ys:Qubits) = 
    if (xs.Length+1 <> ys.Length) then 
        failwith "Error: Takahashi adder requires length(x)+1 = length(y)"
    if xs.Length < 3 then 
        failwith "Error: Takahashi adder requires n >= 3"
    let n = xs.Length 
    
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    for i in (n-1)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-1) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in (n-1)..(-1)..1 do 
        CNOT [ys.[i]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildTakahashiAdderInverse (xs:Qubits) (ys:Qubits) = 
    if (xs.Length+1 <> ys.Length) then 
        failwith "Error: Takahashi adder requires length(x)+1 = length(y)"
    if xs.Length < 3 then 
        failwith "Error: Takahashi adder requires n >= 3"
    let n = xs.Length 
    
    for i in 0..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        CNOT [ys.[i]; xs.[i]]
    for i in (n-1)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
   
let BuildTakahashiModAdder (xs:Qubits) (ys:Qubits) = 
    if (xs.Length <> ys.Length) then 
        failwith "Error: Takahashi mod 2^n adder requires length(x) = length(y)"
    if xs.Length < 3 then 
        failwith "Error: Takahashi adder requires n >= 3"
    let n = xs.Length 
    
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-2) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in (n-1)..(-1)..1 do 
        CNOT [ys.[i]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildCtrlTakahashiAdder (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    
    //following the same pattern as Takahashi, using partial reflection symmetry to save controls
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-2) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-1)..(-1)..1 do 
        CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildCtrlTakahashiModAdder (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    
    //following the same pattern as Takahashi, using partial reflection symmetry to save controls
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-2) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in (n-1)..(-1)..1 do 
        CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

/// The PlusOne adders specialize the Takahashi adder to add the constant 1, i.e. increment an integer by 1. 
/// They were obtained by setting ys.[0] = 1 and ys.[i] = 0 for i in 1..(n-1), but still need qubits for the
/// carry bits.

let BuildTakahashiPlusOneAdder (xs:Qubits) (cs:Qubits) = 
    let n = xs.Length 
    
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-1) do 
        CCNOT [cs.[i-1]; xs.[i]; cs.[i]]
    for i in (n-1)..(-1)..2 do 
        CNOT [cs.[i-1]; xs.[i]]
        CCNOT [cs.[i-2]; xs.[i-1]; cs.[i-1]]
    CNOT [cs.[0]; xs.[1]]
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CNOT [cs.[i-1]; cs.[i]]
    X [xs.[0]]
    for i in 1..(n-1) do 
        CNOT [cs.[i-1]; xs.[i]]
 
let BuildTakahashiPlusOneModAdder (xs:Qubits) (cs:Qubits) = 
    let n = xs.Length 
    
    CNOT [xs.[0]; cs.[0]];
    for i in (n-1)..(-1)..2 do 
        CNOT [cs.[i-1]; xs.[i]]
        CCNOT [cs.[i-2]; xs.[i-1]; cs.[i-1]]
    CNOT [cs.[0]; xs.[1]]
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CNOT [cs.[i-1]; cs.[i]]
    X [xs.[0]]
    for i in 1..(n-1) do 
        CNOT [cs.[i-1]; xs.[i]]   

let BuildCtrlTakahashiPlusOneAdder (xs:Qubits) (cs:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 

    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CCNOT [cs.[i-1]; xs.[i]; cs.[i]]
    MultiplyControlledNOT [cs.[n-2]; xs.[n-1]; ctrl.[0]; cs.[n-1]; cs.[n-3]] // needs one dirty ancilla taken to be cs.[n-3]
    for i in (n-1)..(-1)..2 do 
        CCNOT [cs.[i-1]; ctrl.[0]; xs.[i]]
        CCNOT [cs.[i-2]; xs.[i-1]; cs.[i-1]]
    CCNOT [cs.[0]; ctrl.[0]; xs.[1]]
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CNOT [cs.[i-1]; cs.[i]]
    CNOT [ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [cs.[i-1]; xs.[i]]


let BuildCtrlTakahashiPlusOneModAdder (xs:Qubits) (cs:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CCNOT [cs.[i-1]; xs.[i]; cs.[i]]
    for i in (n-1)..(-1)..2 do 
        CCNOT [cs.[i-1]; ctrl.[0]; xs.[i]]
        CCNOT [cs.[i-2]; xs.[i-1]; cs.[i-1]]
    CCNOT [cs.[0]; ctrl.[0]; xs.[1]]
    CNOT [xs.[0]; cs.[0]]
    for i in 1..(n-2) do 
        CNOT [cs.[i-1]; cs.[i]]
    CNOT [ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [cs.[i-1]; xs.[i]]


// unpackaging from a single qubit registers as required for the addition         
let TakahashiAdder (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs  ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiAdder xs ys 

let TakahashiAdderInverse (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs  ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiAdderInverse xs ys 
    
let TakahashiModAdder (qs:Qubits) = 
    let n = (qs.Length)/2
    let ys = slice qs  ([0..(n-1)]) //
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiModAdder xs ys 

let CtrlTakahashiAdder (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n+1)]
    BuildCtrlTakahashiAdder xs ys ctrl

let CtrlTakahashiModAdder (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n)]
    BuildCtrlTakahashiModAdder xs ys ctrl


let TakahashiPlusOneAdder (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs  ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiPlusOneAdder xs ys 

let TakahashiPlusOneModAdder (qs:Qubits) = 
    let n = (qs.Length)/2
    let ys = slice qs  ([0..(n-1)]) 
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiPlusOneModAdder xs ys 

let CtrlTakahashiPlusOneAdder (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n+1)]
    BuildCtrlTakahashiPlusOneAdder xs ys ctrl

let CtrlTakahashiPlusOneModAdder (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n)]
    BuildCtrlTakahashiPlusOneModAdder xs ys ctrl

/// <summary> Apply an integer adder to a register of qubits </summary>
/// <param name="name">Name of the adder. Currently implemented are: Cuccaro, simple carry ripple, Takahashi</param>
/// <param name="n">Input bit size</param> <param name="s1">n bit input (first addend)</param>
/// <param name="s1">n+1 bit input (second addend). Gets overwritten</param>
/// <param name="ctrl">a bit that (coherently) indicates whether to apply adder or not</param>

let RunAdder (name:string) (n:int) (s1:bigint) (s2:bigint) (ctrl:bigint) = 
    // Prepare the initial state. Dispatching different cases as most adders
    // in the literature have different default way to label inputs and outputs. 
    let arrangeAdderInputs =   
        match name with 
        | "Cuccaro" -> 
            let k = Ket(2*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "Ripple" -> 
            let k = Ket(3*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "Takahashi" -> 
            let k = Ket(2*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "Ctrl-Takahashi" -> 
            let k = Ket(2*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (2*n+1) 1 |> ignore
            qs
        | _ -> failwith "Unknown adder."
              
    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        match name with 
        | "Cuccaro" -> 
            let ADD = CuccaroAdder            // constructing the adder circuit
            let rslt = Array.create (2*n+2) 0 // used to store final result vector after measurement 
            let rsind = [n..2*n]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Ripple" -> 
            let ADD = RippleAdder             // constructing the adder circuit
            let rslt = Array.create (3*n+1) 0 // used to store final result vector after measurement 
            let rsind = [n..2*n]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Takahashi" -> 
            let ADD = TakahashiAdder          // constructing the adder circuit
            let rslt = Array.create (2*n+1) 0 // used to store final result vector after measurement 
            let rsind = [n..2*n]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Ctrl-Takahashi" -> 
            let ADD = CtrlTakahashiAdder      // constructing the adder circuit
            let rslt = Array.create (2*n+2) 0 // used to store final result vector after measurement 
            let rsind = [n..2*n]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | _ -> failwith "Unknown adder."

    let ADD, rslt, rsind = arrangeAdderCircuit
    ADD qs                                    // run the circuit (full quantum sim)
    
    // measure the register 
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1] 
    for i in 0..(qs.Length-1) do 
        rslt.[i] <- qs.[i].Bit.v
    
    let ps      = procStats(true)
    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt

    let sInt = rsind 
               |> List.map (fun i -> rslt.[i])
               |> List.toArray |> IntBool n // have to reverse order to interpret as int
    show "RES: %s adder n s1 s2 rslt %A %A %A %A" name n s1 s2 sInt



let RunPlusOneAdder (name:string) (n:int) (s1:bigint) (ctrl:bigint) = 
    // Prepare the initial state. Dispatching different cases as most adders
    // in the literature have different default way to label inputs and outputs. 
    let arrangeAdderInputs =   
        match name with 
        | "Takahashi-PlusOne" -> 
            let k = Ket(2*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "Takahashi-PlusOneMod" -> 
            let k = Ket(2*n)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "Ctrl-Takahashi-PlusOne" -> 
            let k = Ket(2*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs n 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (2*n+1) 1 |> ignore
            qs
        | "Ctrl-Takahashi-PlusOneMod" -> 
            let k = Ket(2*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs n 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (2*n) 1 |> ignore
            qs
        | _ -> failwith "Unknown adder."
              
    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        match name with 
        | "Takahashi-PlusOne" -> 
            let ADD = TakahashiPlusOneAdder // constructing the adder circuit
            let rslt = Array.create (2*n+1) 0 // used to store final result vector after measurement 
            let rsind = [n..(2*n-1)]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Takahashi-PlusOneMod" -> 
            let ADD = TakahashiPlusOneModAdder // constructing the adder circuit
            let rslt = Array.create (2*n) 0 // used to store final result vector after measurement 
            let rsind = [n..(2*n-1)]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Ctrl-Takahashi-PlusOne" -> 
            let ADD = CtrlTakahashiPlusOneAdder // constructing the adder circuit
            let rslt = Array.create (2*n+2) 0 // used to store final result vector after measurement 
            let rsind = [n..2*n]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | "Ctrl-Takahashi-PlusOneMod" -> 
            let ADD = CtrlTakahashiPlusOneModAdder // constructing the adder circuit
            let rslt = Array.create (2*n+1) 0 // used to store final result vector after measurement 
            let rsind = [n..(2*n-1)]              // picks up the labels for the n+1 bit result
            ADD, rslt, rsind
        | _ -> failwith "Unknown adder."

    let ADD, rslt, rsind = arrangeAdderCircuit
    ADD qs                                    // run the circuit (full quantum sim)
    
    // measure the register 
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1] 
    for i in 0..(qs.Length-1) do 
        rslt.[i] <- qs.[i].Bit.v
    
    let ps      = procStats(true)
    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt

    let sInt = rsind 
               |> List.map (fun i -> rslt.[i])
               |> List.toArray |> IntBool n // have to reverse order to interpret as int
    show "RES: %s adder n s1 rslt %A %A %A" name n s1 sInt

///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular addition
//
///////////////////////////////////////////////////////////////////////////

let BuildModularAdder (xs:Qubits) (ys:Qubits) (ms:Qubits) (tmp:Qubits) =
    let n = xs.Length

    BuildTakahashiAdder xs ys
    BuildTakahashiAdderInverse xs (ms @ [ys.[n]]) 
    X [ys.[n]]
    CNOT [ys.[n]; tmp.[0]]
    X [ys.[n]]
    X [tmp.[0]]
    BuildCtrlTakahashiAdder xs (ms @ [ys.[n]]) [tmp.[0]]
    X [tmp.[0]]
    BuildTakahashiAdderInverse xs ys
    CNOT [ys.[n]; tmp.[0]]
    BuildTakahashiAdder xs ys

let ModularAdder (qs:Qubits) =
    let n = (qs.Length-2)/3
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ms = slice qs [(2*n+1)..(3*n)]
    let tmp = slice qs [(3*n+1)]
    BuildModularAdder xs ys ms tmp

let ModADD (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModularAdder qs 
    if verbose then  
        [a] |> List.map (fun x -> x.GateCount(false,(fun x -> x.Name="CCNOT")))
                        |> List.fold (+) 0 
                        |> show "Number of Toffoli gates = %A" 
    [a] |> List.map (fun x -> MyCircuitExport x) |> List.concat


let RunModularAdder (n:int) (s1:bigint) (s2:bigint) (m:bigint) (verbosity) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    // Prepare the initial state. Dispatching different cases as most adders
    // in the literature have different default way to label inputs and outputs. 
    let arrangeAdderInputs =   
        let k = Ket(3*n+2)
        let qs = k.Qubits
        s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
        s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
        m |> BoolInt n |> PrepBool qs (2*n+1) 1 |> ignore
        qs

    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        let ADD = ModADD verbose1
        let rslt = Array.create (3*n+2) 0   // used to store final result vector after measurement 
        ADD, rslt
  
    let ADD, rslt = arrangeAdderCircuit
               // run the circuit (full quantum sim)
    let ModularAdderCircuit = ADD qs

    let initialState = Array.zeroCreate (3*n+2)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulate ModularAdderCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[n..(2*n-1)] 
    // printfn "Final result = %A" res
    let resInt = (IntBool res.Length res)
    // printfn "As a number mod p this is = %A" resInt
    resInt

///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular integer multiplication
//
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
// Montgomery multiplication
///////////////////////////////////////////////////////////////////////////

let BuildMontgomeryMultiplierForward (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (mg:Qubits) = 
    let n = xs.Length
    
    let AdderRound (a:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (b:Qubits) = 
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ys @ [acc.[n+1]]) a // invariant: MSB of acc is equal to 0
        CNOT [acc.[0]; b.[0]] 
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ms @ [acc.[n+1]]) b // invariant: LSB of acc is equal to 0

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+1]) @ (slice acc [0..(i-1)]) // reindex the accumulator so we never have to rewire
        AdderRound [xs.[i]] ys acci ms [mg.[i]] 
   
let MontgomeryMultiplierForward (qs:Qubits) = 
    let n = (qs.Length-4)/6
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n)]
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let ms = slice qs [(3*n+3)..(4*n+3)]
    let mg = slice qs [(4*n+4)..(5*n+3)]
    BuildMontgomeryMultiplierForward xs ys acc ms mg

let MontgomeryConditionalSubtraction (qs:Qubits) = 
    let n = (qs.Length-4)/6
    let ms  = slice qs [(3*n+3)..(4*n+3)]
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    // check if final result in the accumulator overflows and if so, then subtract the modulus once
    BuildTakahashiAdderInverse (slice acci [0..n]) (ms @ [acci.[n+1]]) 
    BuildCtrlTakahashiModAdder (slice acci [0..n]) ms [acci.[n+1]]

let CopyCircuit (qs:Qubits) = 
    let n = (qs.Length-4)/6
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let res = slice qs [(5*n+4)..(6*n+3)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acci res
 
let MontgomeryMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MontgomeryMultiplierForward qs
    let b = Circuit.Compile MontgomeryConditionalSubtraction qs
    let c = Circuit.Compile CopyCircuit qs
    let d = b.Reverse() // Run the final modular reduction in reverse
    let e = a.Reverse() // Run the sequence of Montgomery steps in reverse
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        [a; b; c; d; e] |> List.map (fun x -> x.GateCount(false,(fun x -> x.Name="CCNOT")))
                        |> List.fold (+) 0 
                        |> show "Number of Toffoli gates = %A" 
    [a; b; c; d; e] |> List.map (fun x -> MyCircuitExport x) |> List.concat


///////////////////////////////////////////////////////////////////////////
// Multiplication modulo a Mersenne number 2^n-1
///////////////////////////////////////////////////////////////////////////

let BuildMersenneMultiplierForward (xs:Qubits) (ys:Qubits) (acc:Qubits) (cs:Qubits) = 
    let n = xs.Length    

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+i])
        BuildCtrlTakahashiAdder acci (ys @ (slice acc [2*n+2])) [xs.[i]] 

    let acc0 = (slice acc [0..(n-1)]) @ (slice acc [2*n+2])
    let acc1 = slice acc [n..(2*n+1)]
    BuildTakahashiAdder acc0 acc1 
    BuildCtrlTakahashiPlusOneModAdder acc0 cs (slice acc [2*n+2])

let MersenneMultiplierForward (qs:Qubits) = 
    let n = (qs.Length-4)/6
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n)]
    let acc = slice qs [(2*n+1)..(4*n+3)]
    let cs = slice qs [(4*n+4)..(5*n+3)]
    BuildMersenneMultiplierForward xs ys acc cs

let MersenneCopyCircuit (qs:Qubits) = 
    let n = (qs.Length-4)/6
    let acc = slice qs [(2*n+1)..(4*n+3)]
    let res = slice qs [(5*n+4)..(6*n+3)]
    let acc0 = slice acc [0..(n-1)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acc0 res
 
let MersenneMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MersenneMultiplierForward qs
    let b = Circuit.Compile MersenneCopyCircuit qs
    let c = a.Reverse() // Run the sequence of multiplication steps in reverse
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        [a; b; c] |> List.map (fun x -> x.GateCount(false,(fun x -> x.Name="CCNOT")))
                        |> List.fold (+) 0 
                        |> show "Number of Toffoli gates = %A" 
    [a; b; c] |> List.map (fun x -> MyCircuitExport x) |> List.concat


///////////////////////////////////////////////////////////////////////////
// Multiplication using modular adder for each bit as in Proos/Zalka
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// Function to run multipliers
///////////////////////////////////////////////////////////////////////////

let RunMultiplier (name:string) (n:int) (s1:bigint) (s2:bigint) (m:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    // Prepare the initial state
    let arrangeMultiplierInputs =   
        match name with 
        | "Montgomery" -> 
            let k = Ket(6*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
            m  |> BoolInt (n+1) |> PrepBool qs (3*n+3) 1 |> ignore
            qs
        | "Mersenne" -> 
            let k = Ket(6*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
            qs
        | _ -> failwith "Unknown multiplier."
                      
    let qs = arrangeMultiplierInputs          // arranges inputs in pattern needed for the multiplier
    
    let arrangeMultiplierCircuit = 
        match name with 
        | "Montgomery" -> 
            let MUL = MontgomeryMultiplier verbose1 // constructing the multiplier circuit
            let rslt = Array.create (6*n+4) 0 // will be used to store final result 
            MUL, rslt
        | "Mersenne" -> 
            let MUL = MersenneMultiplier verbose1 // constructing the multiplier circuit
            let rslt = Array.create (6*n+4) 0 // will be used to store final result 
            MUL, rslt
        | _ -> failwith "This multiplier is not implemented yet."

    let MUL, rslt = arrangeMultiplierCircuit
    let MultiplierCircuit = MUL qs 
    if verbose1 then 
        show "Number of qubits = %A" (6*n+4)

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
//    let qs2 = arrangeMultiplierInputs 
//    let initialState = Array.zeroCreate (6*n+4)
//    List.iter (fun i -> M !!(qs2,i)) [0..qs2.Length-1]
//    for i in 0..(qs2.Length-1) do 
//        initialState.[i] <- qs2.[i].Bit.v

    let initialState = Array.zeroCreate (6*n+4)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulate MultiplierCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[(5*n+4)..(6*n+3)] 
    // printfn "Final result = %A" res
    let resInt = (IntBool res.Length res)
    // printfn "As a number mod p this is = %A" resInt
    resInt




///////////////////////////////////////////////////////////////////////////
//
// Some testing of the adder modules and the Montgomery multiplier module
//
///////////////////////////////////////////////////////////////////////////

[<LQD>]
let __RunSmallAdderTests () = 
    RunAdder "Cuccaro" 4 3I 5I 1I
    RunAdder "Ripple" 4 3I 5I 1I
    RunAdder "Takahashi" 4 3I 5I 1I
    RunAdder "Ctrl-Takahashi" 4 3I 5I 0I
    RunAdder "Ctrl-Takahashi" 4 3I 5I 1I

[<LQD>]
let __RunSmallPlusOneAdderTests () = 
    RunPlusOneAdder "Takahashi-PlusOne" 4 7I 0I
    RunPlusOneAdder "Takahashi-PlusOneMod" 4 5I 0I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOne" 4 7I 0I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOne" 4 7I 1I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOne" 7 112I 0I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOne" 7 112I 1I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOneMod" 3 3I 0I
    RunPlusOneAdder "Ctrl-Takahashi-PlusOneMod" 3 3I 1I
    //RunPlusOneAdder "Ctrl-Takahashi-PlusOneMod" 7 112I 1I

let RunModularAdderTest n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose")
    if verbose then 
        show "Running modular adder circuit test for %A + %A mod %A" s1 s2 p
        show "====================================================="
    let result =  RunModularAdder n s1 s2 p verbosity
    if (s1 + s2) % p = result then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result 
        0 // return integer code zero in case of failure


[<LQD>]
let __RunModularAdderTests () = 
    RunModularAdderTest 4 3I 5I 11I "all" |> ignore
    RunModularAdderTest 4 3I 8I 11I "all" |> ignore
    RunModularAdderTest 4 5I 7I 11I "all" |> ignore
    RunModularAdderTest 4 3I 5I 11I "verbose" |> ignore
    RunModularAdderTest 4 3I 8I 11I "verbose" |> ignore
    RunModularAdderTest 4 5I 7I 11I "verbose" |> ignore


let RunMultiplierTest n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose")
    if verbose then 
        show "Running Montgomery circuit test for %A * %A mod %A" s1 s2 p
        show "====================================================="
    let result =  RunMultiplier "Montgomery" n s1 s2 p verbosity
    let decres =  RunMultiplier "Montgomery" n result 1I p "default"
    let decode1 = RunMultiplier "Montgomery" n s1 1I p "default"
    let decode2 = RunMultiplier "Montgomery" n s2 1I p "default"
    if (decode1 * decode2) % p = decres then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Results: product decprod decs1 decs2 %A %A %A %A\n" n s1 s2 p result decres decode1 decode2
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A\n" n s1 s2 p 
        0 // return integer code zero in case of failure
   
let RunMersenneMultiplierTest n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose")
    if verbose then 
        show "Running Mersenne prime multiplier circuit test for %A * %A mod 2^%A-1" s1 s2 n
        show "====================================================="
    let result =  RunMultiplier "Mersenne" n s1 s2 p verbosity
    if (s1 * s2) % p = result then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Results: product %A\n" n s1 s2 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A\n" n s1 s2 p 
        0 // return integer code zero in case of failure

[<LQD>]
let __RunCheckAllInputsModularArithmeticTest (name:string) (n:int) (p:int) = 
    let P = bigint p
    let SelectedFunction =
        match name with 
        | "Montgomery" -> RunMultiplierTest 
        | "Mersenne" -> RunMersenneMultiplierTest
        | "ModADD" -> RunModularAdderTest
        | _ -> failwith "This function is not implemented yet."
    if name = "Mersenne" && not (P = 2I**n - 1I) then failwith "Modulus is not the right power of 2."
    
    let L = [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (SelectedFunction n i j P "default") ] ]
            |> List.concat
    List.fold (*) 1 L 
    |> function 
        | 1 -> show "Test passed for all inputs: n p = %A %A" n p 
        | 0 -> let ind = List.findIndex (fun i -> (i=0)) L
               show "Test failed for some inputs: n p i j = %A %A %A %A" n p (ind/p) (ind%p) 
        | _ -> show "this should never happen"


        
[<LQD>]
let __RunSmallMultiplierTests () = 
//    RunMultiplierTest 4 12I 11I 13I "verbose" |> ignore
//    RunMultiplierTest 4 11I 3I 13I "verbose"  |> ignore
//    RunMultiplierTest 4 7I 5I 13I "verbose"   |> ignore
//    RunMultiplierTest 5 13I 7I 29I "verbose"  |> ignore
//    RunMultiplierTest 6 7I 9I 63I "verbose"   |> ignore
//    show "\n\n\n\n"
//    RunMultiplierTest 5 5I 6I 31I "verbose"   |> ignore
//    RunMersenneMultiplierTest 5 5I 6I 31I "verbose"   |> ignore
//    RunMultiplierTest 5 23I 12I 31I "verbose"   |> ignore
//    RunMersenneMultiplierTest 5 23I 12I 31I "verbose"   |> ignore
//    RunMultiplierTest 6 7I 43I 63I "verbose"   |> ignore
//    RunMersenneMultiplierTest 6 7I 43I 63I "verbose"   |> ignore
//    RunMultiplierTest 6 7I 11I 63I "all"   |> ignore
//    RunMersenneMultiplierTest 6 7I 11I 63I "all"   |> ignore
    RunMultiplierTest 6 7I 9I 63I "all"   |> ignore
    RunMersenneMultiplierTest 6 7I 9I 63I "all"   |> ignore
    RunMultiplierTest 3 5I 2I 7I "all"   |> ignore
    RunMersenneMultiplierTest 3 5I 2I 7I "all"   |> ignore
    RunMultiplierTest 3 1I 3I 7I "all"   |> ignore
    RunMersenneMultiplierTest 3 1I 3I 7I "all"   |> ignore
    

[<LQD>]
let __RunMediumMultiplierTests () = 
    RunMultiplierTest 10 542I 7I 1021I "verbose"       |> ignore
    RunMultiplierTest 20 500I 7I 1000003I "verbose"    |> ignore
    RunMultiplierTest 32 500I 7I 2000000011I "verbose" |> ignore
    show "\n\n\n\n"
    RunMultiplierTest 10 788I 534I 1023I "verbose"   |> ignore
    RunMersenneMultiplierTest 10 788I 534I 1023I "verbose"   |> ignore

[<LQD>]
let __RunLargeMultiplierTests () = // following are all the primes from the NIST ECC standard (as of April 2010)
    let p192 = 6277101735386680763835789423207666416083908700390324961279I
    let p224 = 26959946667150639794667015087019630673557916260026308143510066298881I
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I
    let p384 = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319I
    let p521 = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151I 
    RunMultiplierTest 192 500I 7I p192 "verbose" |> ignore
    RunMultiplierTest 224 500I 7I p224 "verbose" |> ignore
    RunMultiplierTest 256 500I 7I p256 "verbose" |> ignore
    RunMultiplierTest 384 500I 7I p384 "verbose" |> ignore
    RunMultiplierTest 521 500I 7I p521 "verbose" |> ignore
           
[<LQD>]
let __RunSweepMultiplierTests (startSize:int) (sweepLen:int) = 
    for i in startSize..(startSize+sweepLen) do 
        let p = (pown 2I i)-1I // note: in general this is not prime
        RunMultiplierTest i 500I 7I p "verbose" |> ignore


///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular integer inverse
//
///////////////////////////////////////////////////////////////////////////

let StrictlyLargerThanComparator (qs:Qubits) = 
    let a = Circuit.Compile TakahashiAdder qs
    a.Reverse() |> MyCircuitExport

let MultiplyControlledNOTgate (qs:Qubits) = 
    let a = Circuit.Compile MultiplyControlledNOT qs 
    MyCircuitDump "multCNOT.qc" a
    MyCircuitExport a

[<LQD>]
let __TestMultiplyNOT (n:int) (s:int) =
    // Tools to convert between bigints and list of bools
    let BoolInt n x = 
        [ for i in [0..(n-1)] do yield (int ((x >>> i) % 2)) ] 
    let IntBool n (x: int array) = 
        List.fold (+) 0 <| List.map (fun i -> (1 <<< i) * (x.[i])) [0..x.Length-1]
    let PrepBool (qs:Qubits) offset step ls = 
        List.mapi (fun i x -> if x <> 0 then X !!(qs,(offset+i*step))) ls 

    // Prepare the initial state
    let k = Ket(n+2)
    let qs = k.Qubits
    s |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0 |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    0 |> BoolInt 1 |> PrepBool qs (n+1) 1 |> ignore
    
    let MultiplyControlledNOTCircuit = MultiplyControlledNOTgate qs

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+2)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulate MultiplyControlledNOTCircuit initialState
    printf "The initial state is   %A\n" initialState
    printf "And the final state is %A\n" finalState

    let res = finalState.[0..n+1] 
    printfn "Final result = %A" res
    show "done"

[<LQD>]
let __TestComparator (n:int) (s1:int) (s2:int) =
    // Tools to convert between bigints and list of bools
    let BoolInt n x = 
        [ for i in [0..(n-1)] do yield (int ((x >>> i) % 2)) ] 
    let IntBool n (x: int array) = 
        List.fold (+) 0 <| List.map (fun i -> (1 <<< i) * (x.[i])) [0..x.Length-1]
    let PrepBool (qs:Qubits) offset step ls = 
        List.mapi (fun i x -> if x <> 0 then X !!(qs,(offset+i*step))) ls 

    // Prepare the initial state
    let k = Ket(2*n+1)
    let qs = k.Qubits
    s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
    s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
    0 |> BoolInt 1 |> PrepBool qs (2*n) 1 |> ignore
    
    let ComparatorCircuit = StrictlyLargerThanComparator qs 
        
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (2*n+1)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulate ComparatorCircuit initialState
    printf "The initial state is   %A\n" initialState
    printf "And the final state is %A\n" finalState

    let res = finalState.[0..(2*n)] 
    printfn "Final result = %A" res
//    let resInt = (IntBool res.Length res)
//    printfn "As a number mod p this is = %A" resInt
//    resInt
    show "done"



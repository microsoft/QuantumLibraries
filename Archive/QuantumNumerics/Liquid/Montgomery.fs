// Quantum circuits for Montgomery arithmetic and elliptic curve point additions
// =============================================================================
// M. Roetteler, M. Naehrig, MSR, May 2016

// Todos: 
// - implement in place modular adder: DONE
// - implement out of place integer adder: DONE
// - implement out of place modular adder: DONE
// - introduce type Adder and type Multiplier 
// - introduce type Bitvector: DONE
// - add static methods for adder, using optional arguments: Cuccaro, Ripple, Takahasi
// - add further optional arguments: "mod 2^n", "controlled (by a list of qubits)", "inverse", "inplace"

// - refactor inverse and multiplier to make it input, output, ancilla format
// - rewrite simulator to adapt it to binary ops: DONE
// - write function for circuit "fold" to compute the depth
// - implement mapper for general multipliers: DONE
// - implement entire affine point addition: DONE
// - fix problem with gate counts for Mersenne: DONE

module Microsoft.Research.Liquid.Montgomery

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

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

let BigintHex (s:string) = 
    let dispatch = function
        | '0' -> 0I | '1' -> 1I | '2' -> 2I  | '3' -> 3I | '4' -> 4I | '5' -> 5I | '6' -> 6I | '7' -> 7I 
        | '8' -> 8I | '9' -> 9I | 'A' | 'a' -> 10I | 'B' | 'b' -> 11I | 'C' | 'c' -> 12I | 'D' | 'd' -> 13I  
        | 'E' | 'e' -> 14I | 'F' | 'f' -> 15I | _ -> -1I

    s.ToCharArray()
    |> Array.map dispatch
    |> Array.filter (fun x -> x <> -1I) 
    |> Array.rev
    |> Array.mapi (fun i x -> (pown 16I i) * x)
    |> Array.fold (+) 0I

let HexBigint x = 
    let digitsInt (n:BigInteger) = 
        let rec digitExt x acc = 
            match x with 
            | a when (a=0I) -> acc
            | _ -> digitExt (x/16I) (List.append acc [int (x%16I)])
        digitExt n [] |> List.rev
     
    let dispatch = function
        | 0 -> '0' | 1 -> '1' | 2 -> '2' | 3 -> '3' | 4 -> '4' | 5 -> '5' | 6 -> '6' | 7 -> '7' 
        | 8 -> '8' | 9 -> '9' | 10 -> 'a' | 11 -> 'b' | 12 -> 'c' | 13 -> 'd' | 14 -> 'e' | 15 -> 'f' 
        | _ -> failwith "this should never happen"
        
    digitsInt x 
    |> List.map dispatch 
    |> List.map string 
    |> List.fold (+) "" 

let slice (qs:Qubits) ls = List.map (fun n-> qs.[n]) ls

///////////////////////////////////////////////////////////////////////////
//
// Tools for circuit simulation of Toffoli networks and pretty printing 
//
///////////////////////////////////////////////////////////////////////////

type TofGate = 
    | MyNOT of int 
    | MyCNOT of int * int 
    | MyTOFF of int * int * int 

type eccpoint = 
    | AffinePoint of bigint * bigint
    | ProjectivePoint of bigint * bigint * bigint

 
// bit vector datastructure for fast simulation of Toffoli networks

type BitVector public(n:int) = 
    let groups = (n >>> 5) + 1 // encoding 32 bits into each integer
    let mask   = 0x1F // bit mask 111110...0 to strip off the 5 least significant bits
    let state  = Array.create groups 0u
    let pow2s  = Array.init 32 (fun i -> 1u <<< i) 
    
    let pick(group,element) = (state.[group] &&& element) <> 0u

    member this.State = state

    member this.Init(inputstate:int []) = 
        if inputstate.Length <> n then 
            failwith "Bitvector length and length of initializing array must have same length"
        let iterbound = function 
            | a when (a = (n >>> 5)) -> (min 31 ((n-1) % 32)) 
            | _ -> 31
        for i in 0..((n-1) >>> 5) do 
            state.[i] <- List.fold (+) 0u 
                      <| List.map (fun j -> pow2s.[j] * (uint32 inputstate.[32*i+j])) [0..(iterbound i)]
                
    member this.PrettyPrint() = 
        for i in 0..(n-1) do 
            let group = (i >>> 5)            // i div 32
            let element = pow2s.[i &&& mask] // i mod 32
            match pick(group,element) with 
            | true -> printf "1"
            | false -> printf "0"
        printf "\n"
        () // return a unit

    member this.Run(toffCirc:TofGate list) =
        let rec apply toffCirc = 
            match toffCirc with 
            | [] -> () 
            | (MyNOT i)::gs -> 
                let igroup = (i >>> 5)
                let ibit = pow2s.[i &&& mask]
                state.[igroup] <- state.[igroup] ^^^ ibit                                
                apply gs                
            | (MyCNOT (i, j))::gs -> 
                let igroup = (i >>> 5)
                let jgroup = (j >>> 5)
                let ibit = pow2s.[i &&& mask]
                let jbit = pow2s.[j &&& mask]
                let ival = pick(igroup,ibit) 
                if ival then 
                    state.[jgroup] <- state.[jgroup] ^^^ jbit
                apply gs                
            | (MyTOFF (i, j, k))::gs -> 
                let igroup = (i >>> 5)
                let jgroup = (j >>> 5)
                let kgroup = (k >>> 5)
                let ibit = pow2s.[i &&& mask]
                let jbit = pow2s.[j &&& mask]
                let kbit = pow2s.[k &&& mask]
                let ival = pick(igroup,ibit)
                let jval = pick(jgroup,jbit)
                if (ival && jval) then 
                    state.[kgroup] <- state.[kgroup] ^^^ kbit                
                apply gs
        apply toffCirc

let rec MyCircuitDump (file:string) (c:Circuit) = 
     // pretty prints circuit to a file
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
    | _ ->  failwith "Expo: unrecognized construct or gate"

let MyCircuitSimulateSlow c (b:int []) = 
    // simulates a Toffoli network by applying it to an input bit pattern encodes into ints. 
    // This is very slow as it relies on integer arithemtics.
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

let MyCircuitSimulateFast c (b:int []) = 
    // simulates a Toffoli network by applying it to an input bit pattern encodes into bools. 
    // This is much faster than the integer arithmetic simulator. 
    let n = b.Length 
    let vec = BitVector(n)
    vec.Init(b)
    vec.Run(c) 
    let res = vec.State
    let x = List.concat [ for i in res do yield [ for j in [0..31] do yield (int ((i >>> j) % 2u)) ] ]
    let y = List.toArray x
    y.[0..(n-1)]

///////////////////////////////////////////////////////////////////////////
//
// Decomposing multiply controlled NOT gates with dirty ancillas
//
///////////////////////////////////////////////////////////////////////////

// BuildMultiplyControlledNOT dispatches on small cases first and for n >= 6 it uses the paper 
// [Barenco et al, quant-ph/9503016] for a O(n) implementation of n-fold controlled NOT gates. 
// The constrution uses n control qubits, 1 target qubit and 1 additional ancilla which can be dirty. 
// The dirty is given as list <a>, the control qubits are given as list <cs> and the target is <t>
// Overall cost: for n = 1..5 we obtain 0, 1, 4, 10, 16 and for n >= 6 we obtain 8n-24 Toffoli gates. 

let rec BuildMultiplyControlledNOT (cs:Qubits) (t:Qubits) (a:Qubits) = 
    let n = cs.Length

    let BuildCascade (xs:Qubits) (ys:Qubits) =
        let len = ys.Length-1
        for i in 1..len do 
            CCNOT [ys.[i]; xs.[i-1]; ys.[i-1]]
        CCNOT [xs.[len]; xs.[len+1]; ys.[len]]
        //CCNOT [xs.[len]; xs.[len+1]; xs.[len-1]] 
        for i in len..(-1)..1 do 
            CCNOT [ys.[i]; xs.[i-1]; ys.[i-1]]
        for i in 2..len do 
            CCNOT [ys.[i]; xs.[i-1]; ys.[i-1]]
        CCNOT [xs.[len]; xs.[len+1]; ys.[len]]
        //CCNOT [xs.[len]; xs.[len+1]; xs.[len-1]]
        for i in len..(-1)..2 do 
            CCNOT [ys.[i]; xs.[i-1]; ys.[i-1]]

    match n with 
    | 0 -> X [t.[0]]
    | 1 -> CNOT [cs.[0]; t.[0]] 
    | 2 -> CCNOT [cs.[0]; cs.[1]; t.[0]]
    | 3 | 4 -> BuildMultiplyControlledNOT (slice cs [0..(n-2)]) a [cs.[n-1]]        
               CCNOT [cs.[n-1]; a.[0]; t.[0]]
               BuildMultiplyControlledNOT (slice cs [0..(n-2)]) a [cs.[n-1]]        
               CCNOT [cs.[n-1]; a.[0]; t.[0]]
    | 5 -> BuildMultiplyControlledNOT (slice cs [0..(n-3)]) a [cs.[n-1]]
           BuildMultiplyControlledNOT ((slice cs [(n-2)..(n-1)]) @ a) t [cs.[0]]
           BuildMultiplyControlledNOT (slice cs [0..(n-3)]) a [cs.[n-1]]
           BuildMultiplyControlledNOT ((slice cs [(n-2)..(n-1)]) @ a) t [cs.[0]]           
    | n -> // apply Barenco et al's technique for linear synthesis if n >= 6 
        let n1 = double n |> fun x -> x/2.0 |> ceil |> int |> fun x -> x+1
        let n2 = double n |> fun x -> x/2.0 |> floor |> int 
        BuildCascade (slice cs [0..(n1-1)]) (a @ t @ (slice cs [n1..(2*n1-4)])) 
        match n2 with 
        | n2 when (n2 < 4) -> BuildMultiplyControlledNOT ((slice cs [n1..(n1+n2-2)]) @ a)  t [cs.[0]] 
        | _ -> BuildCascade ((slice cs [n1..(n1+n2-2)]) @ a) (t @ slice cs [0..(n2-3)])
        BuildCascade (slice cs [0..(n1-1)]) (a @ t @ (slice cs [n1..(2*n1-4)])) 
        match n2 with 
        | n2 when (n2 < 4) -> BuildMultiplyControlledNOT ((slice cs [n1..(n1+n2-2)]) @ a)  t [cs.[0]] 
        | _ -> BuildCascade ((slice cs [n1..(n1+n2-2)]) @ a) (t @ slice cs [0..(n2-3)])
               
// MultiplyControlledNOT is a wrapper function for the n-fold controlled NOT using O(n) Toffoli gates. 
let MultiplyControlledNOT (qs:Qubits) = 
    let n = qs.Length-2
    let cs = slice qs [0..(n-1)]
    let t = [ qs.[n] ]
    let a = [ qs.[n+1] ]
    BuildMultiplyControlledNOT cs t a

//// TODO: implement the following tool: take arbitrary circuit as input, then add 1 control
//let rec SingleControlledToffoliCircuit (c:Circuit) = 
//    // converts circuit to list of NOT, CNOT, and Toffoli gates. 
//    match c with
//    | Empty -> []
//    | Seq cs       -> List.concat (List.map SingleControlledToffoliCircuit cs)
//        //[ for c in cs do yield (MyCircuitExport c) ] |> List.concat 
//    | Par cs       -> List.concat (List.map SingleControlledToffoliCircuit cs)
//        //[ for c in cs do yield (MyCircuitExport c) ] |> List.concat 
//    | Apply(g,ws)   ->
//        match g.Name with 
//        | "X" -> [ MyNOT ws.[0] ]
//        | "CNOT" | "CNOT'" -> [ MyCNOT (ws.[0], ws.[1]) ]
//        | "CCNOT" | "CCNOT'" -> [ MyTOFF (ws.[0], ws.[1], ws.[2]) ]
//        | _ -> failwithf "Expo: unrecognized gate %s" g.Name 
//    | _ ->  
//        failwith "Expo: unrecognized construct or gate"


///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: integer addition
//
///////////////////////////////////////////////////////////////////////////

// TODO: more refactoring to be done... (use adder/multiplier/inverter type w/option args)
//       change naming convention from "BuildXYZ" + Wrapper to "XYZ" and move wrapper to tests. 
// TODO: combine all adders into one adder type "QuantumAdder" (ongoing work)

#if FALSE 
type ToffoliNetwork =
    // this type is not used yet; it will be used for 'stitching' of circuits
    { C : TofGate list; M : Map<int,int> }
    
type QuantumAdder (name:string,ctrl:bool,cs:Qubits,mod2n:bool,inv:bool,ancs:Qubits) =  
    // this type is not used yet; it will be used to simplify the various adder types
    // properties
    member this.name  = name 
    member this.ctrl  = ctrl
    member this.cs    = cs
    member this.mod2n = mod2n
    member this.inv   = inv
    member this.ancs  = ancs
    // methods
    member this.finalCarry = 
        match mod2n with 
        | true -> CNOT
        | _    -> I
    // constructors
    new(name,?ctrl,?cs,?mod2n,?inv,?ancs) = 
        let ctrl    = defaultArg ctrl false 
        let cs      = defaultArg cs []
        let mod2n   = defaultArg mod2n false
        let inv     = defaultArg inv false
        let ancs    = defaultArg ancs []
        QuantumAdder(name,ctrl,cs,mod2n,inv,ancs)
            
//let BuildAdder (t:QuantumAdder) (xs:Qubits) (ys:Qubits) = 
//    match t.adderName with 
//    | "Cuccaro" -> 
//        match adderInv with 
//        | false -> 
                
#endif

// case-by-case definitions of adders 

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

// CuccaroAdder is a wrapper function for the Cuccaro in-place addition circuit
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
    
// RippleAdder is a wrapper function for the carry ripple out-of-place adder
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
    if xs.Length < 2 then 
        CCNOT [ys.[0];xs.[0];ys.[1]]
        CNOT [ys.[0];xs.[0]]
    else
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
    if xs.Length < 2 then
        CNOT [ys.[0];xs.[0]]
        CCNOT [ys.[0];xs.[0];ys.[1]]    
    else
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
       
let BuildCtrlTakahashiAdderInverse (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
    MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildMultiCtrlTakahashiAdderInverse (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    MultiplyControlledNOT ([ys.[0]] @ ctrl @ [xs.[0]; ys.[1]]) // ys.[1] taken as dirty ancilla.
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        MultiplyControlledNOT ([ys.[i]] @ ctrl @ [xs.[i]; ys.[0]]) // ys.[0] taken as dirty ancilla.
    MultiplyControlledNOT ([ys.[n-1]; xs.[n-1]] @ ctrl @ [ys.[n]; ys.[n-2]]) // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    MultiplyControlledNOT ([ys.[n-1]] @ ctrl @ [ys.[n]; ys.[0]]) // ys.[0] taken as dirty ancilla. 
    for i in 1..(n-1) do 
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
    MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2] //diff
    for i in (n-1)..(-1)..1 do 
        CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

// Works for small sizes as well but then makes use of an ancilla qubit a
// to construct the multiply-controlled NOT
let BuildCtrlTakahashiAdderInverseSmall (xs:Qubits) (ys:Qubits) (ctrl:Qubits) (a:Qubits) = 
    let n = xs.Length 
    if n >= 1 then
        for i in 1..(n-1) do 
            CNOT [ys.[i]; xs.[i]]
        CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
        for i in (n-2)..(-1)..1 do 
            CNOT [ys.[i]; ys.[i+1]]
        for i in 1..(n-1) do 
            CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
            CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
        if n >= 2 then
            MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2]
        else
            MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; a.[0]]
        for i in (n-2)..(-1)..0 do 
            CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
        for i in 1..(n-2) do 
            CNOT [ys.[i]; ys.[i+1]]
        CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]]
        for i in 1..(n-1) do 
            CNOT [ys.[i]; xs.[i]]

// Works for small sizes as well but then makes use of an ancilla qubit a to 
// construct the multiply-controlled NOT
let BuildCtrlTakahashiAdderSmall (xs:Qubits) (ys:Qubits) (ctrl:Qubits) (a:Qubits) = 
    let n = xs.Length 
    //following the same pattern as Takahashi, using partial reflection symmetry to save controls
    if n >= 1 then
        for i in 1..(n-1) do 
            CNOT [ys.[i]; xs.[i]]
        CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]] 
        for i in (n-2)..(-1)..1 do 
            CNOT [ys.[i]; ys.[i+1]]
        for i in 0..(n-2) do 
            CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
        if n >= 2 then
            MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2] //diff
        else
            MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; a.[0]]
        for i in (n-1)..(-1)..1 do 
            CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
            CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        for i in 1..(n-2) do 
            CNOT [ys.[i]; ys.[i+1]]
        CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
        for i in 1..(n-1) do 
            CNOT [ys.[i]; xs.[i]]

let BuildMultiCtrlTakahashiAdder (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    //following the same pattern as Takahashi, using partial reflection symmetry to save controls
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    MultiplyControlledNOT ([ys.[n-1]] @ ctrl @ [ys.[n]; ys.[0]]) // ys.[0] taken as dirty ancilla. 
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-2) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    MultiplyControlledNOT ([ys.[n-1]; xs.[n-1]] @ ctrl @ [ys.[n]; ys.[n-2]]) // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-1)..(-1)..1 do 
        MultiplyControlledNOT ([ys.[i]] @ ctrl @ [xs.[i]; ys.[0]]) // ys.[0] taken as dirty ancilla.
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    MultiplyControlledNOT ([ys.[0]] @ ctrl @ [xs.[0]; ys.[1]]) // ys.[1] taken as dirty ancilla.
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
   
let BuildCtrlTakahashiModAdder (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
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

let BuildTakahashiModAdderInverse (xs:Qubits) (ys:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    CNOT [ys.[0]; xs.[0]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        CNOT [ys.[i]; xs.[i]]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildCtrlTakahashiModAdderInverse (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildMultiCtrlTakahashiModAdder (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 0..(n-2) do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in (n-1)..(-1)..1 do 
        BuildMultiplyControlledNOT ([ys.[i]] @ ctrl) [xs.[i]] [ys.[0]] // dirty ancilla = ys.[0]
        //CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    BuildMultiplyControlledNOT ([ys.[0]] @ ctrl) [xs.[0]] [ys.[1]] // dirty ancilla = ys.[1]
    //CCNOT [ys.[0]; ctrl.[0]; xs.[0]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]

let BuildMultiCtrlTakahashiModAdderInverse (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    BuildMultiplyControlledNOT ([ys.[0]] @ ctrl) [xs.[0]] [ys.[1]] // dirty ancilla = ys.[1]
    //CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]]
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        BuildMultiplyControlledNOT ([ys.[i]] @ ctrl) [xs.[i]] [ys.[0]] // dirty ancilla = ys.[0]
        //CCNOT [ys.[i]; ctrl.[0]; xs.[i]]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
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

let CtrlTakahashiAdderInverse (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n+1)]
    BuildCtrlTakahashiAdderInverse xs ys ctrl

let CtrlTakahashiModAdder (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n)]
    BuildCtrlTakahashiModAdder xs ys ctrl

let TakahashiModAdderInverse (qs:Qubits) = 
    let n = (qs.Length)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    BuildTakahashiModAdderInverse xs ys

let CtrlTakahashiModAdderInverse (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n)]
    BuildCtrlTakahashiModAdderInverse xs ys ctrl

let MultiCtrlTakahashiModAdderInverse (numCtrl:int) (qs:Qubits) = 
    let n = (qs.Length-numCtrl)/2
    let ys = slice qs ([0..(n-1)])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n)..qs.Length-1]
    BuildMultiCtrlTakahashiModAdderInverse xs ys ctrl
   
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

let MultiplyControlledNOTgate (qs:Qubits) = 
    let a = Circuit.Compile MultiplyControlledNOT qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
//    if verbose then 
//        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
//    MyCircuitDump "multCNOT.qc" a
    gates                         
    
let CtrlStrictlyLargerThanComparator (qs:Qubits) = 
    CtrlTakahashiAdderInverse qs 
    
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
// Quantum counters
//
///////////////////////////////////////////////////////////////////////////

// table of irreducible sparse polynomials over GF(2) of degrees 2..100
let IRRED_TABLE_GF2 = 
    [ [ 1 ]; [ 1 ]; [ 1 ]; [ 2 ]; [ 1 ]; [ 1 ]; [ 1; 3; 4 ]; [ 1 ]; [ 3 ]; [ 2 ]; [ 3 ]; [ 1; 3; 4 ]; [ 5 ]; [ 1 ]; [ 1; 3; 5 ]; [ 3 ]; [ 3 ]; [ 1; 2; 5 ]; [ 3 ]; [ 2 ]; 
      [ 1 ]; [ 5 ]; [ 1; 3; 4 ]; [ 3 ]; [ 1; 3; 4 ]; [ 1; 2; 5 ]; [ 1 ]; [ 2 ]; [ 1 ]; [ 3 ]; [ 2; 3; 7 ]; [ 10 ]; [ 7 ]; [ 2 ]; [ 9 ]; [ 1; 4; 6 ]; [ 1; 5; 6 ]; [ 4 ]; 
      [ 3; 4; 5 ]; [ 3 ]; [ 7 ]; [ 3; 4; 6 ]; [ 5 ]; [ 1; 3; 4 ]; [ 1 ]; [ 5 ]; [ 2; 3; 5 ]; [ 9 ]; [ 2; 3; 4 ]; [ 1; 3; 6 ]; [ 3 ]; [ 1; 2; 6 ]; [ 9 ]; [ 7 ]; [ 2; 4; 7 ]; 
      [ 4 ]; [ 19 ]; [ 2; 4; 7 ]; [ 1 ]; [ 1; 2; 5 ]; [ 29 ]; [ 1 ]; [ 1; 3; 4 ]; [ 18 ]; [ 3 ]; [ 1; 2; 5 ]; [ 9 ]; [ 2; 5; 6 ]; [ 1; 3; 5 ]; [ 6 ]; [ 3; 9; 10 ]; [ 25 ]; 
      [ 35 ]; [ 1; 3; 6 ]; [ 21 ]; [ 2; 5; 6 ]; [ 3; 5; 6 ]; [ 9 ]; [ 2; 4; 9 ]; [ 4 ]; [ 1; 3; 8 ]; [ 2; 4; 7 ]; [ 5 ]; [ 1; 2; 8 ]; [ 21 ]; [ 13 ]; [ 2; 6; 7 ]; [ 38 ]; 
      [ 27 ]; [ 1; 5; 8 ]; [ 21 ]; [ 2 ]; [ 21 ]; [ 11 ]; [ 6; 9; 10 ]; [ 6 ]; [ 11 ]; [ 1; 3; 6 ]; [ 15 ] ]

let PRIMITIVE_TABLE_GF2 = 
    [ [ 1 ]; [ 1 ]; [ 1 ]; [ 2 ]; [ 1 ]; [ 1 ]; [ 2; 3; 4 ]; [ 4 ]; [ 3 ]; [ 2 ]; [ 1; 3; 5; 6; 7 ]; [ 1; 3; 4 ]; [ 3; 5; 7 ]; [ 1 ]; [ 2; 3; 5 ]; [ 3 ]; [ 1; 10; 12 ]; 
      [ 1; 2; 5 ]; [ 3 ]; [ 2 ]; [ 1 ]; [ 5 ]; [ 1; 3; 4 ]; [ 3 ]; [ 1; 4; 6; 7; 8; 10; 14 ]; [ 1; 2; 5 ]; [ 2; 5; 6; 7; 13 ]; [ 2 ]; [ 1; 2; 3; 5; 7; 11; 13; 16; 17 ]; 
      [ 3 ]; [ 3; 4; 7; 9; 15 ]; [ 3; 6; 8; 10; 11; 12; 13 ]; [ 1; 2; 4; 5; 6; 7; 8; 11; 12; 15; 16 ]; [ 2 ]; [ 1; 5; 6; 8; 13; 14; 17; 19; 20; 22; 23 ]; [ 1; 4; 6 ]; 
      [ 1; 5; 6 ]; [ 4 ]; [ 3; 4; 5 ]; [ 3 ]; [ 1; 2; 5; 6; 9; 11; 12; 18; 20; 24; 25; 26; 30 ]; [ 3; 4; 6 ]; [ 1; 3; 4; 16; 17; 19; 24 ]; [ 1; 3; 4 ]; [ 14; 17; 20; 21; 23 ]; 
      [ 5 ]; [ 3; 7; 8; 10; 11; 12; 17; 23; 25 ]; [ 9 ]; [ 2; 3; 4 ]; [ 1; 3; 6 ]; [ 3 ]; [ 1; 2; 6 ]; [ 1; 2; 4; 7; 13; 15; 16; 17; 18; 21; 25; 27; 29; 30; 31; 32; 34 ]; 
      [ 4; 7; 9; 10; 11 ]; [ 2; 4; 7 ]; [ 1; 2; 3; 4; 5; 6; 8; 10; 11; 13; 16; 19; 21 ]; [ 19 ]; [ 2; 4; 7 ]; [ 1 ]; [ 1; 2; 5 ]; 
      [ 1; 6; 12; 13; 14; 16; 17; 18; 19; 20; 21; 24; 25; 26; 27; 28; 29; 30; 32  ]; [ 1 ]; [ 1; 3; 4 ]; [ 18 ]; [ 2; 4; 5; 6; 7; 11; 13; 16; 19; 21; 22; 23; 24; 27; 29; 30; 31; 32; 33; 34; 38; 40; 42; 45; 46 ]; 
      [ 1; 2; 5 ]; [ 9 ]; [ 2; 5; 6 ]; [ 1; 3; 5 ]; [ 6 ]; [ 3; 9; 10 ]; [ 25 ]; [ 3; 8; 11; 12; 13; 16; 17; 21; 24; 26; 27; 28; 29; 30; 31; 32; 33; 34; 35;  36; 37 ]; [ 1; 3; 6 ]; 
      [ 1; 2; 5; 14; 15; 19; 20; 23; 24; 25; 27; 29; 31; 33; 34; 35; 36; 37; 38 ]; [ 2; 5; 6 ]; [ 1; 4; 6; 7; 8; 9; 11; 13; 14; 15; 17; 18; 19; 21; 22; 24; 25; 27; 28; 32;  33; 34; 35; 37; 39; 42; 43; 44; 46 ]; 
      [ 9 ]; [ 2; 4; 9 ]; [ 4 ]; [ 1; 2; 3; 4; 6; 7; 9; 11; 12; 13; 15; 16; 17; 21; 22; 24; 28; 32; 33; 34;  35; 36; 37 ]; [ 2; 4; 7 ]; [ 2; 4; 5; 7; 8; 9; 11; 12; 13; 16; 17; 20; 21; 24; 25; 26; 28; 31; 32; 37;  39; 41; 46; 47; 48; 52; 57 ];
      [ 1; 2; 8 ]; [ 1; 2; 5; 6; 7; 8; 9; 10; 12; 14; 16; 21; 22; 25; 27; 31; 38; 40; 43 ]; [ 13 ]; [ 1; 3; 4; 5; 7; 8; 9; 11; 12; 13; 14; 15; 18; 19; 21; 24; 27; 28; 30; 31;  33; 36; 38; 41; 44; 45; 47 ]; [ 38 ]; [ 2; 4; 7; 10; 11; 12; 14; 16; 18; 24; 25; 26; 28; 30; 32; 34; 36; 37; 38;  39; 40; 42; 44; 45; 46; 47; 48; 50; 51; 53; 55; 56; 57; 58; 60; 61; 62; 63;  64 ]; [ 1; 5; 8 ]; [ 2; 4; 5; 9; 12; 14; 15; 16; 17; 18; 21; 22; 23; 24; 30; 32; 33; 37; 38;  40; 41; 42; 43; 44; 48 ]; [ 2 ]; [ 21 ]; [ 11 ]; [ 6; 9; 10 ]; [ 6 ]; [ 11 ]; [ 1; 2; 4; 6; 32; 33; 34; 35; 64; 65; 66; 67; 96; 97; 98 ]; 
      [ 3; 5; 6; 8; 9; 11; 15; 16; 19; 20; 22; 24; 25; 27; 30; 31; 34; 35; 36; 37; 41; 43; 44; 45; 46; 47; 48; 52; 55; 56; 57 ]]

let ShiftRegisterCounter (qs:Qubits) = 
    let n = qs.Length
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 0..(L.Length-1) do
        CNOT [qs.[n-1]; qs.[L.[i]-1]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]

// ControlShiftRegisterCounter is a counter that cycles through the states of a linear shift register 
// given by a primitive polynomial p in Galois form. The initial state can be an arbitrary non-zero state.
// The Toffoli count is given by the sparsity of p plus n, where n is the length of the LFSR. 
// Assumed format for qs: n bits for LFSR, 1 bit for control
let ControlShiftRegisterCounter (qs:Qubits) =
    let n = qs.Length-1
    let cs = [qs.[n]]
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 0..(L.Length-1) do
        CCNOT [cs.[0]; qs.[n-1]; qs.[L.[i]-1]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CCNOT [cs.[0]; qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CCNOT [cs.[0]; qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]

// ControlShiftRegisterCounterSpecial is a counter that cycles through the states of a linear shift register 
// given by a primitive polynomial p in Galois form. The initial state must be the all ones vector. 
// The Toffoli count is given by the sparsity of p. 
// Assumed format for qs: n bits for LFSR, 1 bit for control
let BuildControlShiftRegisterCounterSpecial (qs:Qubits) (cs:Qubits) =
    let n = qs.Length
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 0..(L.Length-1) do
        CCNOT [cs.[0]; qs.[n-1]; qs.[L.[i]-1]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]

let BuildShiftRegisterCounterSpecial (qs:Qubits) =
    let n = qs.Length
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 0..(L.Length-1) do
        CNOT [qs.[n-1]; qs.[L.[i]-1]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]

let BuildShiftRegisterCounterSpecialInverse (qs:Qubits) =
    let n = qs.Length
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 0..(L.Length-1) do
        CNOT [qs.[n-1]; qs.[L.[i]-1]]

let BuildControlShiftRegisterCounterSpecialInverse (qs:Qubits) (cs:Qubits) =
    let n = qs.Length
    if n<2 then failwith "n must be at least 2"
    //let L = IRRED_TABLE_GF2.[n-2] // table starts at degree n=2
    let L = PRIMITIVE_TABLE_GF2.[n-2] // table starts at degree n=2
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 0..(L.Length-1) do
        CCNOT [cs.[0]; qs.[n-1]; qs.[L.[i]-1]]
    
let ControlShiftRegisterCounterSpecial (qs:Qubits) =
    let n = qs.Length
    let xs = slice qs [0..(n-2)]
    let cs = slice qs [(n-1)]
    BuildControlShiftRegisterCounterSpecial xs cs 

let BuildBinaryDoubling (qs:Qubits) = 
    let n = qs.Length
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]]
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]] 
        CNOT [qs.[i]; qs.[n-i]]

let BuildBinaryHalfing (qs:Qubits) = 
    let n = qs.Length
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        CNOT [qs.[n-i]; qs.[i]] 
        CNOT [qs.[i]; qs.[n-i]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        CNOT [qs.[n-i-1]; qs.[i]] 
        CNOT [qs.[i]; qs.[n-i-1]]

let BuildCtrlBinaryDoubling (qs:Qubits) (cs:Qubits) (a:Qubits) = 
    let n = qs.Length
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        BuildMultiplyControlledNOT (cs @ [qs.[n-i-1]]) [qs.[i]] a
        CNOT [qs.[i]; qs.[n-i-1]]
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        BuildMultiplyControlledNOT (cs @ [qs.[n-i]]) [qs.[i]] a
        CNOT [qs.[i]; qs.[n-i]]

let BuildCtrlBinaryHalfing (qs:Qubits) (cs:Qubits) (a:Qubits) = 
    let n = qs.Length
    for i in 1..(n-1)/2 do 
        CNOT [qs.[i]; qs.[n-i]]
        BuildMultiplyControlledNOT (cs @ [qs.[n-i]]) [qs.[i]] a
        CNOT [qs.[i]; qs.[n-i]]
    for i in 0..(n/2-1) do 
        CNOT [qs.[i]; qs.[n-i-1]]
        BuildMultiplyControlledNOT (cs @ [qs.[n-i-1]]) [qs.[i]] a
        CNOT [qs.[i]; qs.[n-i-1]]

let CtrlBinaryDoubling (n:int) (qs:Qubits) = 
    let xs = slice qs [0..(n-1)]
    let cs = slice qs [n..qs.Length-2]
    BuildCtrlBinaryDoubling xs cs [qs.[qs.Length-1]]

let CtrlBinaryHalfing (n:int) (qs:Qubits) = 
    let xs = slice qs [0..(n-1)]
    let cs = slice qs [n..qs.Length-2]
    BuildCtrlBinaryHalfing xs cs [qs.[qs.Length-1]]
    
[<LQD>]
let __RunShiftRegisterCounter (n:int) (s:int) (rounds:int) = 
    let k = Ket(n)
    let qs = k.Qubits
    s |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
     
    let ShiftRegisterCounterCircuit = 
        Circuit.Compile ShiftRegisterCounter qs 
        |> MyCircuitExport

    let currentState = Array.zeroCreate n
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        currentState.[i] <- qs.[i].Bit.v
    
    let nextState = Array.zeroCreate n
    for i in 0..(rounds-1) do 
        let nextState = MyCircuitSimulateFast ShiftRegisterCounterCircuit currentState
        show "The current state is   %A" currentState
        for i in 0..(qs.Length-1) do currentState.[i] <- nextState.[i]

    show "done"

[<LQD>]
let __FindPeriodShiftRegister (n:int) (s:int) (c:int) = 
    let b  = match c with | 0 -> 1 | _ -> double (c+1) |> log |> fun x -> x/(log 2.0) |> ceil |> int 
    if b > 1 then failwith "ShiftRegisterCounterCircuit currently only implemented for 1 control"
    let k = Ket(n+b)
    let qs = k.Qubits
    s |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    c |> fun x -> (bigint x) |> BoolInt b |> PrepBool qs n 1 |> ignore

//    let ShiftRegisterCounterCircuit = 
//        Circuit.Compile ControlShiftRegisterCounter qs
//        |> MyCircuitExport

    let ShiftRegisterCounterCircuitSpecial = 
        Circuit.Compile ControlShiftRegisterCounterSpecial qs 
        |> MyCircuitExport

    let initialState = (BoolInt n (bigint s) @ BoolInt b (bigint c)) |> List.toArray 

    let currentState = Array.zeroCreate (n+b)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        currentState.[i] <- qs.[i].Bit.v
        
    let currentState = MyCircuitSimulateFast ShiftRegisterCounterCircuitSpecial currentState
    let mutable per = 1 
    while currentState <> initialState do 
        per <- per + 1
        let nextState = MyCircuitSimulateFast ShiftRegisterCounterCircuitSpecial currentState
        show "The current state is   %A" currentState
        for i in 0..(qs.Length-1) do currentState.[i] <- nextState.[i]

    show "The period of %A is %A\n" s per
    show "done"

let CoupledLFSRCircuit (qs:Qubits) =
    let n = (qs.Length-1)/2
    let xs = slice qs [0..(n-1)]
    let fs = slice qs [n]
    let ys = slice qs [(n+1)..(2*n)]
    BuildMultiplyControlledNOT xs fs [ys.[0]] // using ys.[0] as dirty ancilla
    X fs
    BuildMultiplyControlledNOT ys fs [xs.[0]] // using xs.[0] as dirty ancilla
    BuildControlShiftRegisterCounterSpecial xs fs 
    X fs 
    BuildControlShiftRegisterCounterSpecial ys fs 
    X fs 
    
[<LQD>]
let __RunTwoCoupledLFSRs (n:int) (s:int) = 
    let initCounter1 = (bigint s)
    let initCounter2 = (pown 2I n) - 1I
    let k = Ket(2*n+1)
    let qs = k.Qubits
    initCounter1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
    1I           |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    initCounter2 |> BoolInt n |> PrepBool qs (n+1) 1 |> ignore

    let LFSRCircuit = 
        Circuit.Compile CoupledLFSRCircuit qs
        |> MyCircuitExport

    let currentState = Array.zeroCreate (2*n+1)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        currentState.[i] <- qs.[i].Bit.v
        
    for i in 0..((pown 2 (n+1))-1) do 
        show "The current state is   %A" currentState
        let nextState = MyCircuitSimulateFast LFSRCircuit currentState
        for i in 0..(qs.Length-1) do currentState.[i] <- nextState.[i]

    show "done"

  
///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular addition
//
///////////////////////////////////////////////////////////////////////////

let BuildModularAdder (xs:Qubits) (ys:Qubits) (ms:Qubits) (tmp:Qubits) =
    // Computes the sum x + y mod m. The circuit is almost identical to 
    // the one in Fig. 4 of http://arxiv.org/abs/quant-ph/9511018v1.
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

let BuildCtrlModularAdder (xs:Qubits) (ys:Qubits) (ms:Qubits) (tmp:Qubits) (c:Qubits) =
    let n = xs.Length

    BuildTakahashiAdder xs ys  // (1) If c=0, cancels with (1) below
    BuildCtrlTakahashiAdderInverse xs (ms @ [ys.[n]]) c // Controlled inverse
    X [ys.[n]] // (3) If c=0, cancels with (3) below
    CCNOT [c.[0]; ys.[n]; tmp.[0]] // Controlled
    X [ys.[n]] // (3)
    X [tmp.[0]] // (6) If c=0, cancels with (6) below 
    CCNOT [c.[0]; tmp.[0]; tmp.[1]] // (7) Multiply controlling the next adder
    BuildCtrlTakahashiAdder xs (ms @ [ys.[n]]) [tmp.[1]] // Controlled by c and tmp.[0]
    CCNOT [c.[0]; tmp.[0]; tmp.[1]] // Undoing (7)
    X [tmp.[0]] // (6)
    BuildTakahashiAdderInverse xs ys // (1)
    CCNOT [c.[0]; ys.[n]; tmp.[0]] // Controlled
    BuildCtrlTakahashiAdder xs ys c // Controlled

let BuildCtrlNegator (xs:Qubits) (ms:Qubits) (c:Qubits) =
    let n = xs.Length
    
    BuildCtrlTakahashiModAdderInverse ms xs c
    BuildCtrlTakahashiModAdder xs ms c
    for i in 0..(n-1) do 
        CNOT  [xs.[i]; ms.[i]]
        CCNOT [c.[0]; ms.[i]; xs.[i]]
        CNOT  [xs.[i]; ms.[i]]
     
let ModularAdder (qs:Qubits) =
    let n = (qs.Length-1)/3
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ms = slice qs [(2*n+1)..(3*n)]
    let tmp = slice qs [(3*n+1)]
    BuildModularAdder xs ys ms tmp

let CtrlModularAdder (qs:Qubits) =
    let n = (qs.Length-4)/3
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ms = slice qs [(2*n+1)..(3*n)]
    let tmp = slice qs [(3*n+1)..(3*n+2)]
    let c = slice qs [(3*n+3)]
    BuildCtrlModularAdder xs ys ms tmp c

let CtrlNegator (qs:Qubits) =
    let n = (qs.Length-1)/2
    let xs = slice qs [0..(n-1)]
    let ms = slice qs [n..(2*n-1)]
    let c = slice qs [(2*n)]
    BuildCtrlNegator xs ms c

let ModADD (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModularAdder qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let CtrlModADD (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile CtrlModularAdder qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let CtrlNeg (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile CtrlNegator qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let RunModularAdder (name:string) (n:int) (s1:bigint) (s2:bigint) (m:bigint) (ctrl:bigint) (verbosity) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    // Prepare the initial state. Dispatching different cases as most adders
    // in the literature have different default way to label inputs and outputs. 
    let arrangeAdderInputs =   
         match name with 
         | "ModADD" -> 
            // name: || xs | ys | tmp | modulus | tmp ||
            // bits: || n  | n  | 1   | n       | 1   ||
            // init: || x  | y  | 0   | p       | 0   ||
            let k = Ket(3*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            m |> BoolInt n |> PrepBool qs (2*n+1) 1 |> ignore
            qs
         | "CtrlModADD" -> 
            // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | tmp | modulus | tmp | ctrl ||
            // bits: || n  | n  | 2   | n       | 1   | 1    ||
            // init: || x  | y  | 00  | p       | 0   | c    ||
            let k = Ket(3*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            m |> BoolInt n |> PrepBool qs (2*n+1) 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (3*n+3) 1 |> ignore
            qs
         | "CtrlNeg" -> 
            let k = Ket(2*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            m |> BoolInt n |> PrepBool qs n 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (2*n) 1 |> ignore
            qs
         | _ -> failwith "Unknown modular adder."
    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        match name with 
        | "ModADD" -> 
            let ADD = ModADD verbose1
            let rslt = Array.create (3*n+2) 0   // used to store final result vector after measurement 
            if verbose1 then 
                show "Number of qubits = %A" (3*n+2)
            ADD, rslt
        | "CtrlModADD" ->
            let ADD = CtrlModADD verbose1
            let rslt = Array.create (3*n+4) 0   // used to store final result vector after measurement 
            if verbose1 then 
                show "Number of qubits = %A" (3*n+4)
            ADD, rslt
        | "CtrlNeg" ->
            let ADD = CtrlNeg verbose1
            let rslt = Array.create (2*n+1) 0   // used to store final result vector after measurement 
            if verbose1 then 
                show "Number of qubits = %A" (2*n+1)
            ADD, rslt
        | _ -> failwith "Unknown modular adder."
  
    let ADD, rslt = arrangeAdderCircuit
               // run the circuit (full quantum sim)
    let ModularAdderCircuit = ADD qs

    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast ModularAdderCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = 
        match name with 
        | "ModADD" | "CtrlModADD" -> finalState.[n..(2*n-1)] 
        | "CtrlNeg" -> finalState.[0..(n-1)]
        | _ -> failwith "Unknown modular adder type"
    // printfn "Final result = %A" res
    let resInt = (IntBool res.Length res)
    // printfn "As a number mod p this is = %A" resInt
    resInt

let CompileModularAdder (name:string) (n:int) (s1:bigint) (s2:bigint) (m:bigint) (ctrl:bigint) (verbosity) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    // Prepare the initial state. Dispatching different cases as most adders
    // in the literature have different default way to label inputs and outputs. 
    let arrangeAdderInputs =   
         match name with 
         | "ModADD" -> 
            // Prepare the initial state for modular adder. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | tmp | modulus | tmp ||
            // bits: || n  | n  | 1   | n       | 1   ||
            // init: || x  | y  | 0   | p       | 0   ||
            let k = Ket(3*n+2)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (3*n+2)
            qs
         | "CtrlModADD" -> 
            // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | tmp | modulus | tmp | ctrl ||
            // bits: || n  | n  | 2   | n       | 1   | 1    ||
            // init: || x  | y  | 00  | p       | 0   | c    ||
            let k = Ket(3*n+4)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (3*n+4)
            qs
         | "CtrlNeg" -> 
            // Prepare the initial state for controlled modular negator. The partitioning of the quantum register is as follows: 
            // name: || xs | modulus | ctrl ||
            // bits: || n  | n       | 1    ||
            // init: || x  | p       | c    ||
            let k = Ket(2*n+1)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (2*n+1)
            qs
         | _ -> failwith "Unknown modular adder."
    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        match name with 
        | "ModADD" -> 
            let ADD = ModADD verbose1
            let rslt = Array.create (3*n+2) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlModADD" ->
            let ADD = CtrlModADD verbose1
            let rslt = Array.create (3*n+4) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlNeg" ->
            let ADD = CtrlNeg verbose1
            let rslt = Array.create (2*n+1) 0   // used to store final result vector after measurement 
            ADD, rslt
        | _ -> failwith "Unknown modular adder."
  
    let ADD, rslt = arrangeAdderCircuit           
    let ModularAdderCircuit = ADD qs
    ModularAdderCircuit

///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular doubling 
//
///////////////////////////////////////////////////////////////////////////

// BuildModularDoubler computes x -> 2*x mod m for odd m and replaces xs with the result. The implementation 
// resembles the prose description in Section 4.3.2 of Proos/Zalka http://arxiv.org/abs/quant-ph/0301141 and 
// only cleans up the ancillas correctly if m is odd. 
let BuildModularDoubler (xs:Qubits) (ms:Qubits) (a:Qubits) =
    let n = xs.Length
    BuildBinaryDoubling (xs @ [a.[0]]) // Double input by bit shift.
    // Subtract m to see whether modular reduction is necessary.
    BuildTakahashiAdderInverse xs (ms @ [a.[0]])
    CNOT [a.[0]; a.[1]] // Copy out control bit.
    // Conditional addition of m, if result of subtraction was negative.
    BuildCtrlTakahashiAdder xs (ms @ [a.[0]]) [a.[1]]
    // BuildCtrlTakahashiAdder (slice xtmp [0..(n-1)]) (ms @ [xtmp.[n]]) [a.[1]]
    CNOT [xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    X [a.[1]]

let BuildModularHalver (xs:Qubits) (ms:Qubits) (a:Qubits) =
    let n = xs.Length
    X [a.[1]]
    CNOT [xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    BuildCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) [a.[1]]
    CNOT [a.[0]; a.[1]] // Copy out control bit.
    BuildTakahashiAdder xs (ms @ [a.[0]])
    BuildBinaryHalfing (xs @ [a.[0]]) // Halving of input by bit shift.

let BuildCtrlModularDoubler (xs:Qubits) (ms:Qubits) (cs:Qubits) (a:Qubits) =
    let n = xs.Length
    BuildCtrlBinaryDoubling (xs @ [a.[0]]) cs [a.[1]] // Halving of input by bit shift. Using a.[1] as dirty ancilla.
    BuildCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) cs
    CCNOT [cs.[0]; a.[0]; a.[1]] // Copy out control bit.
    BuildMultiCtrlTakahashiAdder xs (ms @ [a.[0]]) [cs.[0]; a.[1]]
    CCNOT [cs.[0]; xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    CNOT [cs.[0]; a.[1]]

let BuildCtrlModularHalver (xs:Qubits) (ms:Qubits) (cs:Qubits) (a:Qubits) =
    let n = xs.Length
    CNOT [cs.[0]; a.[1]]
    CCNOT [cs.[0]; xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    BuildMultiCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) [cs.[0]; a.[1]]
    CCNOT [cs.[0]; a.[0]; a.[1]] // Copy out control bit.
    BuildCtrlTakahashiAdder xs (ms @ [a.[0]]) cs
    BuildCtrlBinaryHalfing (xs @ [a.[0]]) cs [a.[1]] // Halving of input by bit shift. Using a.[1] as dirty ancilla.

let ModularDoubler (qs:Qubits) =
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ms = slice qs [n..(2*n-1)]
    let a = slice qs [2*n..(2*n+1)]
    BuildModularDoubler xs ms a

let ModularHalver (qs:Qubits) =
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ms = slice qs [n..(2*n-1)]
    let a = slice qs [2*n..(2*n+1)]
    BuildModularHalver xs ms a

let ModDBL (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModularDoubler qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let RunModularDoubler (name:string) (n:int) (s:bigint) (m:bigint) (verbosity) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    let arrangeDoublerInputs =   
         match name with 
         | "ModDBL" -> 
            let k = Ket(2*n+2)
            let qs = k.Qubits
            s |> BoolInt n |> PrepBool qs 0 1 |> ignore
            m |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
         | _ -> failwith "Unknown modular adder."
    let qs = arrangeDoublerInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeDoublerCircuit = 
        match name with 
        | "ModDBL" -> 
            let DBL = ModDBL verbose1
            let rslt = Array.create (2*n+2) 0   // used to store final result vector after measurement 
            DBL, rslt
        | _ -> failwith "Unknown modular adder."
  
    let DBL, rslt = arrangeDoublerCircuit
               // run the circuit (full quantum sim)
    
    let ModularDoublerCircuit = DBL qs 
    // note: modular halfing can be implemented by running this in reverse: 
    // let ModularDoublerCircuit = DBL qs |> List.rev 

    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast ModularDoublerCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[0..(n-1)] 
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
    //BuildCtrlTakahashiAdder (slice acci [0..n]) ms [acci.[n+1]]

let CopyCircuit (qs:Qubits) = 
    let n, flag = 
        match ((qs.Length-5) % 6) with 
        | 0 -> (qs.Length-5)/6, true // this is a control copy circuit conditioned on the last qubit
        | _ -> (qs.Length-4)/6, false// this is a regular copy circuit
    let cs = slice qs [(qs.Length-1)..(qs.Length-1)]
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let res = slice qs [(5*n+4)..(6*n+3)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) (cs:Qubits) = 
        for i in 0..(n-1) do 
           match flag with 
           | true -> CCNOT [cs.[0]; acc.[i]; res.[i]]        
           | _    ->  CNOT [acc.[i]; res.[i]]        
    CopyRegs acci res cs
 
let MontgomeryMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MontgomeryMultiplierForward qs
    let b = Circuit.Compile MontgomeryConditionalSubtraction qs
    let c = Circuit.Compile CopyCircuit qs
    let d = b.Reverse() // Run the final modular reduction in reverse
    let e = a.Reverse() // Run the sequence of Montgomery steps in reverse
    //let f = Circuit.Seq [a; b; c; d; e]
    //f.Render("ModMont.htm", "svg", 1, 50.0)
    let gates = 
        [a; b; c; d; e] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let CtrlMontgomeryMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MontgomeryMultiplierForward qs
    let b = Circuit.Compile MontgomeryConditionalSubtraction qs
    let c = Circuit.Compile CopyCircuit qs
    let d = b.Reverse() // Run the final modular reduction in reverse
    let e = a.Reverse() // Run the sequence of Montgomery steps in reverse
    //let f = Circuit.Seq [a; b; c; d; e]
    //f.Render("ModMont.htm", "svg", 1, 50.0)
    let gates = 
        [a; b; c; d; e] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         


let MontgomeryMultiplierFast (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MontgomeryMultiplierForward qs |> MyCircuitExport
    let b = Circuit.Compile MontgomeryConditionalSubtraction qs |> MyCircuitExport
    let c = Circuit.Compile CopyCircuit qs |> MyCircuitExport
    let d = List.rev b // Run the final modular reduction in reverse
    let e = List.rev a // Run the sequence of Montgomery steps in reverse
    let gates = [a; b; c; d; e] |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
        
///////////////////////////////////////////////////////////////////////////
// Montgomery squarer
///////////////////////////////////////////////////////////////////////////

let BuildMontgomerySquarerForward (xs:Qubits) (acc:Qubits) (ms:Qubits) (mg:Qubits) = 
    let n = xs.Length-2
    
    let AdderRound (a:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (b:Qubits) = 
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ys @ [acc.[n+1]]) a // invariant: MSB of acc is equal to 0
        CNOT [acc.[0]; b.[0]] 
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ms @ [acc.[n+1]]) b // invariant: LSB of acc is equal to 0

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+1]) @ (slice acc [0..(i-1)]) // reindex the accumulator so we never have to rewire
        CNOT [xs.[i]; xs.[n+1]]
        AdderRound [xs.[n+1]] (slice xs [0..n]) acci ms [mg.[i]] 
        CNOT [xs.[i]; xs.[n+1]]

let MontgomerySquarerForward (qs:Qubits) = 
    let n = (qs.Length-5)/5
    let xs = slice qs [0..(n+1)] // n bits of input x and 2 ancilla bits
    let acc = slice qs [(n+2)..(2*n+3)]
    let ms = slice qs [(2*n+4)..(3*n+4)]
    let mg = slice qs [(3*n+5)..(4*n+4)]
    BuildMontgomerySquarerForward xs acc ms mg
           
let MontgomerySquarerConditionalSubtraction (qs:Qubits) = 
    let n = (qs.Length-5)/5
    let ms  = slice qs [(2*n+4)..(3*n+4)]
    let acc = slice qs [(n+2)..(2*n+3)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    // check if final result in the accumulator overflows and if so, then subtract the modulus once
    BuildTakahashiAdderInverse (slice acci [0..n]) (ms @ [acci.[n+1]]) 
    BuildCtrlTakahashiModAdder (slice acci [0..n]) ms [acci.[n+1]]
    //BuildCtrlTakahashiAdder (slice acci [0..n]) ms [acci.[n+1]]

let CopySquarerCircuit (qs:Qubits) = 
    let n = (qs.Length-5)/5
    let acc = slice qs [(n+2)..(2*n+3)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let res = slice qs [(4*n+5)..(5*n+4)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acci res
 
let MontgomerySquarer (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MontgomerySquarerForward qs
    let b = Circuit.Compile MontgomerySquarerConditionalSubtraction qs
    let c = Circuit.Compile CopySquarerCircuit qs
    let d = b.Reverse() // Run the final modular reduction in reverse
    let e = a.Reverse() // Run the sequence of Montgomery steps in reverse
    //let f = Circuit.Seq [a; b; c; d; e]
    //f.Render("ModMont.htm", "svg", 1, 50.0)
    // Compute the number of Toffoli gates in the circuit and write to console
    let gates = 
        [a; b; c; d; e] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
///////////////////////////////////////////////////////////////////////////
// Multiplication modulo a Mersenne number 2^n-1
///////////////////////////////////////////////////////////////////////////

let BuildMersenneMultiplierForward (xs:Qubits) (ys:Qubits) (acc:Qubits) (cs:Qubits) = 
    let n = xs.Length    

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+i])
        BuildCtrlTakahashiAdder (slice acci [0..(n-1)]) (ys @ [acci.[n]]) [xs.[i]] 

    let acc0 = slice acc [0..(n-1)]
    let acc1 = slice acc [n..(2*n)]
    BuildTakahashiAdder acc0 acc1 
    BuildCtrlTakahashiPlusOneModAdder acc0 cs [acc.[2*n]]

let MersenneMultiplierForward (qs:Qubits) = 
    let n = (qs.Length-1)/6
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(4*n)]
    let cs = slice qs [(4*n+1)..(5*n)]
    BuildMersenneMultiplierForward xs ys acc cs

let MersenneCopyCircuit (qs:Qubits) = 
    let n = (qs.Length-1)/6
    let acc = slice qs [(2*n)..(4*n)]
    let res = slice qs [(5*n+1)..(6*n)]
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
    let gates = 
        [a; b; c] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates        

///////////////////////////////////////////////////////////////////////////
// Squaring modulo a Mersenne number 2^n-1
///////////////////////////////////////////////////////////////////////////

let BuildMersenneSquarerForward (xs:Qubits) (acc:Qubits) (cs:Qubits) (tmp:Qubits) = 
    let n = xs.Length    

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+i])
        CNOT [xs.[i]; tmp.[0]]
        BuildCtrlTakahashiAdder (slice acci [0..(n-1)]) (xs @ [acci.[n]]) [tmp.[0]] 
        CNOT [xs.[i]; tmp.[0]]

    let acc0 = slice acc [0..(n-1)]
    let acc1 = slice acc [n..(2*n)]
    BuildTakahashiAdder acc0 acc1 
    BuildCtrlTakahashiPlusOneModAdder acc0 cs [acc.[2*n]]

let MersenneSquarerForward (qs:Qubits) = 
    let n = (qs.Length-2)/5
    let xs = slice qs [0..(n-1)]
    let acc = slice qs [(n)..(3*n)]
    let cs = slice qs [(3*n+1)..(4*n)]
    let tmp = slice qs [(4*n+1)]
    BuildMersenneSquarerForward xs acc cs tmp

let MersenneSquarerCopyCircuit (qs:Qubits) = 
    let n = (qs.Length-2)/5
    let acc = slice qs [(n)..(3*n)]
    let res = slice qs [(4*n+2)..(5*n+1)]
    let acc0 = slice acc [0..(n-1)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acc0 res

let MersenneSquarer (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile MersenneSquarerForward qs
    let b = Circuit.Compile MersenneSquarerCopyCircuit qs
    let c = a.Reverse() // Run the sequence of multiplication steps in reverse
    let gates = 
        [a; b; c] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
///////////////////////////////////////////////////////////////////////////
// Multiplication modulo a Pseudo-Mersenne number 2^n-c, c small
///////////////////////////////////////////////////////////////////////////

let BuildPseudoMersenneMultiplierForward (xs:Qubits) (ys:Qubits) (acc:Qubits) (cs:Qubits) (C:int) = 
    let n = xs.Length
    let mutable c = C  
    let mutable cbits = Array.empty
    while c > 0 do
        cbits <- (Array.append cbits [| c % 2 |])
        c <- c/2
    //printf "\n%A\n\n" cbits

    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+i])
        BuildCtrlTakahashiAdder (slice acci [0..(n-1)]) (ys @ [acci.[n]]) [xs.[i]] 

    let acc0 = slice acc [0..(n-1)]
    let acc1 = slice acc [n..(2*n)]
    for i in [0..(cbits.Length-1)] do
        if cbits.[i] = 1 then
            BuildTakahashiAdder acc0 acc1 
            // Reduce acc0 depending on acc.[2*n]
        // Modular doubling of acc1

let SetPseudoMersenneMultiplierForward (C:int) =
    let PseudoMersenneMultiplierForward (qs:Qubits) = 
        let n = (qs.Length-1)/6
        let xs = slice qs [0..(n-1)]
        let ys = slice qs [n..(2*n-1)]
        let acc = slice qs [(2*n)..(4*n)]
        let cs = slice qs [(4*n+1)..(5*n)]
        BuildPseudoMersenneMultiplierForward xs ys acc cs C
    PseudoMersenneMultiplierForward

let PseudoMersenneCopyCircuit (qs:Qubits) = 
    let n = (qs.Length-1)/6
    let acc = slice qs [(2*n)..(4*n)]
    let res = slice qs [(5*n+1)..(6*n)]
    let acc0 = slice acc [0..(n-1)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acc0 res

let SetPseudoMersenneMultiplier (C:int) =
    let PseudoMersenneMultiplier (verbose:bool) (qs:Qubits) = 
        let a = Circuit.Compile (SetPseudoMersenneMultiplierForward C) qs
        let b = Circuit.Compile PseudoMersenneCopyCircuit qs
        let c = a.Reverse() // Run the sequence of multiplication steps in reverse
        // Compute the number of Toffoli gates in the circuit and write to console
        let gates = 
            [a; b; c] 
            |> List.map (fun x -> MyCircuitExport x) 
            |> List.concat
        // Compute the number of Toffoli gates in the circuit and write to console
        if verbose then 
            show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
        gates                         
    PseudoMersenneMultiplier

///////////////////////////////////////////////////////////////////////////
// Multiplication using modular DBL-ADD for each bit as in Proos/Zalka
///////////////////////////////////////////////////////////////////////////

let BuildModDBLADDMultiplier (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..1 do 
        BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ slice acc [n]) ms tmp [xs.[i]]
        BuildModularDoubler (slice acc [0..(n-1)]) ms ([acc.[n]] @ tmp)
    BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ [acc.[n]]) ms tmp [xs.[0]] 
    
let ModDblAddMultiplier (qs:Qubits) = 
    let n = (qs.Length-3)/4
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(3*n)]
    let ms = slice qs [(3*n+1)..(4*n)]
    let tmp = slice qs [(4*n+1)..(4*n+2)]
    BuildModDBLADDMultiplier xs ys acc ms tmp

let ModDBLADDMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModDblAddMultiplier qs
    //a.Render("ModDBLADD.htm", "svg", 1, 50.0)
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let BuildCtrlModDBLADDMultiplier (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) (ctrl:Qubits) = 
    let n = xs.Length
    // Making the ModDBLADD multiplier controlled. tmp needs one more qubit for multiple control
    for i in (n-1)..(-1)..1 do
        CCNOT [ctrl.[0]; xs.[i]; tmp.[2]] 
        BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ slice acc [n]) ms tmp [tmp.[2]]
        CCNOT [ctrl.[0]; xs.[i]; tmp.[2]]
        BuildCtrlModularDoubler (slice acc [0..(n-1)]) ms [ctrl.[0]] ([acc.[n]] @ tmp) 
    CCNOT [ctrl.[0]; xs.[0]; tmp.[2]]
    BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ [acc.[n]]) ms tmp [tmp.[2]]
    CCNOT [ctrl.[0]; xs.[0]; tmp.[2]]

//let BuildCtrlModularDoubler (xs:Qubits) (ms:Qubits) (cs:Qubits) (a:Qubits) =
//let BuildModularDoubler (xs:Qubits) (ms:Qubits) (a:Qubits) =
//        BuildModularDoubler (slice acc [0..(n-1)]) ms ([acc.[n]] @ tmp)

let CtrlModDblAddMultiplier (qs:Qubits) = 
    let n = (qs.Length-5)/4
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(3*n)]
    let ms = slice qs [(3*n+1)..(4*n)]
    let tmp = slice qs [(4*n+1)..(4*n+3)]
    let ctrl = slice qs [(4*n+4)]
    BuildCtrlModDBLADDMultiplier xs ys acc ms tmp ctrl

let CtrlModDBLADDMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile CtrlModDblAddMultiplier qs
    //a.Render("ModDBLADD.htm", "svg", 1, 50.0)
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates    


#if FALSE
let BuildModDBLADDMultiplier (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..1 do 
        BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ slice acc [n]) ms tmp [xs.[i]]
        BuildModularDoubler (slice acc [0..(n-1)]) ms ([acc.[n]] @ tmp)
    BuildCtrlModularAdder (slice acc [0..(n-1)]) (ys @ [acc.[n]]) ms tmp [xs.[0]] 
    
//let BuildModDBLADDMultiplier (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) = 
//    let n = xs.Length
//
//    for i in (n-1)..(-1)..1 do 
//        let acci = (slice acc [(i+2)..n]) @ (slice acc [0..(i+1)]) // reindex the accumulator 
//        BuildCtrlModularAdder (slice acci [0..(n-1)]) (ys @ slice acci [n]) ms tmp [xs.[i]]
//        BuildModularDoubler (slice acci [0..(n-1)]) ms ([acci.[n]] @ tmp)
//    let acc0 = (slice acc [2..n]) @ (slice acc [0..1])
//    BuildCtrlModularAdder (slice acc0 [0..(n-1)]) (ys @ slice acc0 [n]) ms tmp [xs.[0]]

let ModDblAddMultiplier (qs:Qubits) = 
    let n = (qs.Length-3)/4
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(3*n)]
    let ms = slice qs [(3*n+1)..(4*n)]
    let tmp = slice qs [(4*n+1)..(4*n+2)]
    BuildModDBLADDMultiplier xs ys acc ms tmp

let ModDBLADDMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModDblAddMultiplier qs
    //a.Render("ModDBLADD.htm", "svg", 1, 50.0)
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
#endif 

// Specializing to Mersenne numbers.
let BuildModDBLADDMersenneMultiplier (xs:Qubits) (ys:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..0 do 
        let acci = (slice acc [(i+1)..(n-1)]) @ (slice acc [0..i]) // reindex the accumulator 
        BuildCtrlModularAdder (slice acci [0..(n-1)]) (ys @ [tmp.[2]]) ms tmp [xs.[i]]
  
let ModDblAddMersenneMultiplier (qs:Qubits) = 
    let n = (qs.Length-3)/4
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(3*n-1)]
    let ms = slice qs [(3*n)..(4*n-1)]
    let tmp = slice qs [(4*n)..(4*n+2)]
    BuildModDBLADDMersenneMultiplier xs ys acc ms tmp

let ModDBLADDMersenneMultiplier (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModDblAddMersenneMultiplier qs
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

///////////////////////////////////////////////////////////////////////////
// Squaring specialized from modular DBL-ADD multiplication (Proos/Zalka)
///////////////////////////////////////////////////////////////////////////

let BuildModDBLADDSquarer (xs:Qubits) (acc:Qubits) (ms:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..1 do 
        let acci = (slice acc [(i+2)..n]) @ (slice acc [0..(i+1)]) // reindex the accumulator 
        CNOT [xs.[i]; tmp.[2]]
        BuildCtrlModularAdder (slice acci [0..(n-1)]) (xs @ slice acci [n]) ms tmp [tmp.[2]]
        CNOT [xs.[i]; tmp.[2]]
        BuildModularDoubler acci ms tmp
    let acc0 = (slice acc [2..n]) @ (slice acc [0..1])
    CNOT [xs.[0]; tmp.[2]]
    BuildCtrlModularAdder (slice acc0 [0..(n-1)]) (xs @ slice acc0 [n]) ms tmp [tmp.[2]]
    CNOT [xs.[0]; tmp.[2]]

let ModDblAddSquarer (qs:Qubits) = 
    let n = (qs.Length-4)/3
    let xs = slice qs [0..(n-1)]
    let acc = slice qs [(n)..(2*n)]
    let ms = slice qs [(2*n+1)..(3*n)]
    let tmp = slice qs [(3*n+1)..(3*n+3)]
    BuildModDBLADDSquarer xs acc ms tmp

let ModDBLADDSquarer (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModDblAddSquarer qs
    //a.Render("ModDBLADDSqu.htm", "svg", 1, 50.0)

    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

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
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | modulus | ancilla | mg-rounds | result || 
            // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      ||
            // init: || x  | y  | 0       | 0    | p       | 0       | 0..0      | 0      || 
            let k = Ket(6*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
            m  |> BoolInt (n+1) |> PrepBool qs (3*n+3) 1 |> ignore
            qs
        | "MontgomerySquarer" -> 
            // Prepare the initial state for Montgomery squarer. The partitioning of the quantum register is as follows: 
            // name: || xs | ancilla | acc  | modulus | ancilla | mg-rounds | result || 
            // bits: || n  | 2       | n+2  | n       | 1       | n         | n      ||
            // init: || x  | 0       | 0    | p       | 0       | 0..0      | 0      || 
            let k = Ket(5*n+5)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            m  |> BoolInt (n+1) |> PrepBool qs (2*n+4) 1 |> ignore
            qs
        | "Mersenne" -> 
            // Prepare the initial state for Mersenne. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(6*n+1)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            qs
        | "MersenneSqu" -> 
            // Prepare the initial state for Mersenne squaring. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(5*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            qs
//        | "PseudoMersenne" -> 
//            let k = Ket(6*n+1)
//            let qs = k.Qubits
//            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
//            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
//            qs
        | "DoubleAdd" -> 
            // Prepare the initial state for double&add. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(4*n+3)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            m  |> BoolInt n |> PrepBool qs (3*n+1) 1 |> ignore
            qs
        | "DoubleAddSqu" -> 
            // Prepare the initial state for double&add squaring. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(3*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            m  |> BoolInt n |> PrepBool qs (2*n+1) 1 |> ignore
            qs
        | "DoubleAddMersenne" -> 
            // Prepare the initial state for double&add Mersenne. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(4*n+3)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            m  |> BoolInt n |> PrepBool qs (3*n) 1 |> ignore
            qs
        | _ -> failwith "Unknown multiplier."
                      
    let qs = arrangeMultiplierInputs   // arranges inputs in pattern needed for the multiplier
    
    let arrangeMultiplierCircuit = // constructing the multiplier circuit
        match name with 
        | "Montgomery" -> 
            //let MUL = MontgomeryMultiplier verbose1 
            let MUL = MontgomeryMultiplierFast verbose1 
            MUL
        | "MontgomerySquarer" -> 
            let MUL = MontgomerySquarer verbose1 
            MUL
        | "Mersenne" -> 
            let MUL = MersenneMultiplier verbose1 
            MUL
        | "MersenneSqu" -> 
            let MUL = MersenneSquarer verbose1 
            MUL
//        | "PseudoMersenne" -> 
//            let MUL = (SetPseudoMersenneMultiplier 11) verbose1 
//            MUL
        | "DoubleAdd" -> 
            let MUL = ModDBLADDMultiplier verbose1 
            MUL
        | "DoubleAddSqu" -> 
            let MUL = ModDBLADDSquarer verbose1 
            MUL
        | "DoubleAddMersenne" -> 
            let MUL = ModDBLADDMersenneMultiplier verbose1 
            MUL
        | _ -> failwith "This multiplier is not implemented yet."

    let MUL = arrangeMultiplierCircuit
    let MultiplierCircuit = MUL qs 
    if verbose1 then 
        show "Number of qubits = %A" qs.Length

    let initialState = Array.zeroCreate qs.Length
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MultiplierCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = 
        match name with
        | "Montgomery" -> 
            finalState.[(5*n+4)..(6*n+3)]
        | "MontgomerySquarer" -> 
            finalState.[(4*n+5)..(5*n+4)]
        | "Mersenne" -> 
            finalState.[(5*n+1)..(6*n)]
        | "MersenneSqu" -> 
            finalState.[(4*n+2)..(5*n+1)]
//        | "PseudoMersenne" -> 
//            finalState.[(5*n+1)..(6*n)]
        | "DoubleAdd" -> 
            finalState.[(2*n)..(3*n-1)] 
            //Array.append finalState.[(2*n+2)..(3*n)] [|finalState.[2*n]|] // Result is shifted by 2 bits.
        | "DoubleAddSqu" -> 
            Array.append finalState.[(n+2)..(2*n)] [|finalState.[n]|] // Result is shifted by 2 bits.
        | "DoubleAddMersenne" -> 
            Array.append finalState.[(2*n+1)..(3*n-1)] [|finalState.[2*n]|] // Result is shifted by 1 bits.
        | _ -> failwith "Wrong multiplier name in reading result from final state."
    let resInt = (IntBool res.Length res)
    resInt


let CompileModularMultiplier (name:string) (n:int) (s1:bigint) (s2:bigint) (m:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    
    // Prepare the initial state
    let arrangeMultiplierInputs =   
        match name with 
        | "Montgomery" -> 
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | modulus | ancilla | mg-rounds | result || 
            // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      ||
            // init: || x  | y  | 0       | 0    | p       | 0       | 0..0      | 0      || 
            let k = Ket(6*n+4)
            let qs = k.Qubits
            qs
        | "CtrlMontgomery" -> 
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | modulus | ancilla | mg-rounds | result | ctrl || 
            // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      | 1    ||
            // init: || x  | y  | 0       | 0    | p       | 0       | 0..0      | 0      | ctrl || 
            let k = Ket(6*n+5)
            let qs = k.Qubits
            qs
        | "MontgomerySquarer" -> 
            // Prepare the initial state for Montgomery squarer. The partitioning of the quantum register is as follows: 
            // name: || xs | ancilla | acc  | modulus | ancilla | mg-rounds | result || 
            // bits: || n  | 2       | n+2  | n       | 1       | n         | n      ||
            // init: || x  | 0       | 0    | p       | 0       | 0..0      | 0      || 
            let k = Ket(5*n+5)
            let qs = k.Qubits
            qs
        | "Mersenne" -> 
            // Prepare the initial state for Mersenne. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(6*n+1)
            let qs = k.Qubits
            qs
        | "MersenneSqu" -> 
            // Prepare the initial state for Mersenne squaring. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(5*n+2)
            let qs = k.Qubits
            qs
//        | "PseudoMersenne" -> 
//            let k = Ket(6*n+1)
//            let qs = k.Qubits
//            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
//            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
//            qs
        | "DoubleAdd" -> 
            // Prepare the initial state for double&add. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | anc || 
            // bits: || n  | n  | n+1  | n       | 1   ||
            // init: || x  | y  | 0    | p       | 0   || 
            // note that the final result is stored in the lower n bits of acc. 
            let k = Ket(4*n+3)
            let qs = k.Qubits
            qs
        | "CtrlDoubleAdd" -> 
            // Prepare the initial state for double&add. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | anc | ctrl || 
            // bits: || n  | n  | n+1  | n       | 3   | 1    ||
            // init: || x  | y  | 0    | p       | 000 | ctrl || 
            // note that the final result is stored in the lower n bits of acc. 
            let k = Ket(4*n+5)
            let qs = k.Qubits
            qs
        | "DoubleAddSqu" -> 
            // Prepare the initial state for double&add squaring. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(3*n+4)
            let qs = k.Qubits
            qs
        | "DoubleAddMersenne" -> 
            // Prepare the initial state for double&add Mersenne. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | mg-rounds | result || 
            // bits: || n  | n  | n+2  | n       | n         | n      ||
            // init: || x  | y  | 0    | p       | 0..0      | 0      || 
            let k = Ket(4*n+3)
            let qs = k.Qubits
            qs
        | _ -> failwith "Unknown multiplier."
                      
    let qs = arrangeMultiplierInputs   // arranges inputs in pattern needed for the multiplier
    
    let arrangeMultiplierCircuit = // constructing the multiplier circuit
        match name with 
        | "Montgomery" -> 
            let MUL = MontgomeryMultiplier verbose1 
            MUL
        | "CtrlMontgomery" -> 
            let MUL = CtrlMontgomeryMultiplier verbose1 
            MUL
        | "MontgomerySquarer" -> 
            let MUL = MontgomerySquarer verbose1 
            MUL
        | "Mersenne" -> 
            let MUL = MersenneMultiplier verbose1 
            MUL
        | "MersenneSqu" -> 
            let MUL = MersenneSquarer verbose1 
            MUL
//        | "PseudoMersenne" -> 
//            let MUL = (SetPseudoMersenneMultiplier 11) verbose1 
//            MUL
        | "CtrlDoubleAdd" -> 
            let MUL = CtrlModDBLADDMultiplier verbose1 
            MUL
        | "DoubleAdd" -> 
            let MUL = ModDBLADDMultiplier verbose1 
            MUL
        | "DoubleAddSqu" -> 
            let MUL = ModDBLADDSquarer verbose1 
            MUL
        | "DoubleAddMersenne" -> 
            let MUL = ModDBLADDMersenneMultiplier verbose1 
            MUL
        | _ -> failwith "This multiplier is not implemented yet."

    let MUL = arrangeMultiplierCircuit
    let MultiplierCircuit = MUL qs 
    if verbose2 then 
        show "Number of qubits = %A" qs.Length
    MultiplierCircuit  

///////////////////////////////////////////////////////////////////////////
//
// Quantum arithmetic: modular integer inversion
//
///////////////////////////////////////////////////////////////////////////

// BuildMontgomeryInverseForward implements the inverse mod p of numbers that are given in Montgomery foom x * 2^n mod p. The output is 
// again in Montogomery form. The implementation follows closely the Kalinski paper which in turn is a binary extended Euclidean algorithm. 
// The implementation uses a total of 14+6+2*ceil(log(n)) qubits. The registers are: 
// us = intially equal to modulus p, vs = intitially equal to input x, ss = intially 1, rs = intiially 0, ms = modulus p, mg = 8n qubits 
// holding the results of the dispatch predicates for each round, fs = flag that determines whether we are in compute or in counter mode, 
// ks = a register of qubits to implement a counter from 0 up to worst case iteration bound of 2n, and a = 1 ancilla. 

let BuildMontgomeryInverseForward (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (mg:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits) = 
    let EuclidRound (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (round:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits)  = 
        // step 1: check termination condition: is vs the all zero string
        X >< vs
        BuildMultiplyControlledNOT vs [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X >< vs
        // step 2: check counter state and flip flag
        BuildMultiplyControlledNOT ks [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X fs
        // step 3: compute match statement to detect next step in binary Euclid
        // step 3a: check whether u is even 
        X us
        CCNOT [fs.[0]; us.[0]; round.[0]] 
        X us 
        // step 3b: or else, check whether v is even
        X vs 
        X [round.[0]]
        BuildMultiplyControlledNOT [fs.[0]; round.[0]; vs.[0]] [round.[1]] [ss.[0]]  
        X [round.[0]]
        X vs
        // step 3c: or else, check whether u > v
        BuildCtrlTakahashiAdderInverse vs (us @ a) [fs.[0]] 
        X >< [round.[0]; round.[1]]
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; round.[0]; round.[1]] [round.[2]] [ss.[0]]
        // step 3d: else, we have that u <= v
        X a
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; round.[0]; round.[1]] [round.[3]] [ss.[0]]
        X a
        X >< [round.[0]; round.[1]]
        BuildCtrlTakahashiAdder vs (us @ a) [fs.[0]] // clean up the ancilla a
        // step 4: dispatch 4 mutually exclusive cases
        // case 4a: u even 
        BuildCtrlBinaryDoubling ss [round.[0]; fs.[0]] [a.[0]] // using a as dirty ancilla
        BuildCtrlBinaryHalfing  us [round.[0]; fs.[0]] [a.[0]]
        // case 4b: v even 
        BuildCtrlBinaryDoubling rs [round.[1]; fs.[0]] [a.[0]]
        BuildCtrlBinaryHalfing  vs [round.[1]; fs.[0]] [a.[0]]
        // case 4c: u odd, v odd, u > v 
        BuildMultiCtrlTakahashiModAdderInverse us vs [round.[2]; fs.[0]] // replace (u,v) -> (u-v,v)
        BuildCtrlBinaryHalfing  us [round.[2]; fs.[0]] [a.[0]] // replace (u-v) -> (u-v)/2
        BuildMultiCtrlTakahashiModAdder rs ss [round.[2]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling ss [round.[2]; fs.[0]] [a.[0]] // replace s -> 2s        
        // case 4d: u odd, v odd, u <= v 
        BuildMultiCtrlTakahashiModAdderInverse vs us [round.[3]; fs.[0]] // replace (u,v) -> (u,v-u)
        BuildCtrlBinaryHalfing  vs [round.[3]; fs.[0]] [a.[0]] // replace (v-u) -> (v-u)/2
        BuildMultiCtrlTakahashiModAdder ss rs [round.[3]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling rs [round.[3]; fs.[0]] [a.[0]] // replace r -> 2r      
        // step 5: increase counter if required
        X fs
        BuildControlShiftRegisterCounterSpecial ks fs 
        X fs
    
    let n = us.Length
    for i in 0..(2*n-1) do 
        let round = (slice mg [4*i..4*i+3]) // get 4 bits per round for match statement on (u, v)
        EuclidRound us vs ss rs round fs a ks 
   
let InverterModularReduction (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/14
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let rs = slice qs [(3*n+1)..(4*n+1)] // n+1 bits for the non-reduced register holding r
    // check if final result in r overflows and if so, then subtract the modulus once
    BuildTakahashiAdderInverse (slice rs [0..(n-1)]) (ms @ [rs.[n]]) 
    BuildCtrlTakahashiModAdder (slice rs [0..(n-1)]) ms [rs.[n]]

let InverterFixAlmostInverse (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/14
    let ks = slice qs [14*n+4..14*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [14*n+6+counterSize..14*n+5+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [13*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [14*n+4+counterSize..14*n+5+counterSize] // 2 clean ancillas 
    // replace counter register content k -> k -n 
    for i in 0..(n-1) do 
        BuildModularHalver rs ms a 
        //BuildShiftRegisterCounterSpecialInverse ks        
    // count down from k-n and in each step take half of the r register, mapping -x^-1 2^k mod p to -x^-1 2^k 2^(k-n) = -x^-1 2^n mod p 
    for i in 0..(n-1) do 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        BuildCtrlModularDoubler rs ms fs a 
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1 2^n)=x^-1 2^n mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    X >< ms // one's complement
    BuildTakahashiModAdder rs ms
    X >< ms // return modulus to original state
    X >< rs // one's complement of output   
         
let InverterCopyCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/14
    let rs = slice qs [(3*n+1)..(4*n)] // n bits for the reduced register holding r
    let res = slice qs [13*n+4..14*n+3] // n bits for the final result
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (rs:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [rs.[i]; res.[i]]        
    CopyRegs rs res
 
let MontgomeryInverseForward (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/14
    let us = slice qs [0..(n-1)] // u (n bits)
    let vs = slice qs [n..(2*n-1)] // v (n bits)
    let ss = slice qs [(2*n)..(3*n)] // s (n+1 bits)
    let rs = slice qs [(3*n+1)..(4*n+1)] // r (n+1 bits); holds the almost MG inverse
    let mg = slice qs [(5*n+2)..(13*n+1)] // round status indicators
    let fs = slice qs [13*n+2] // flag needed for 1st while loop
    let a  = slice qs [14*n+4+counterSize] // ancilla
    let ks = slice qs [14*n+4..14*n+3+counterSize] // shift register counter 
    BuildMontgomeryInverseForward us vs ss rs mg fs a ks
 
let MontgomeryInverse (counterSize:int) (verbose:bool) (qs:Qubits) =     
    let a = Circuit.Compile (MontgomeryInverseForward counterSize) qs
    let b = Circuit.Compile (InverterModularReduction counterSize) qs
    let c = Circuit.Compile (InverterFixAlmostInverse counterSize) qs
    let d = Circuit.Compile (InverterCopyCircuit counterSize) qs
    let e = c.Reverse() // Run the correction to fix almost inverse in reverse
    let f = b.Reverse() // Run the final modular reduction in reverse
    let g = a.Reverse() // Run the sequence of Montgomery steps in reverse
    // Compute the number of Toffoli gates in the circuit and write to console
    let gates = 
        [a; b; c; d; e; f; g] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let RunInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    let initCounter =  (pown 2I counterSize) - 1I

    // Prepare the initial state. The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | modulus | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | n       | 8n        | 2     | n      | log(n)   | 2        | log(n)   ||
    // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 00       | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "Montgomery" -> 
            let k = Ket(14*n+6+2*counterSize)
            let qs = k.Qubits
            m  |> BoolInt n |> PrepBool qs 0 1 |> ignore            // u register: n bits
            x  |> BoolInt n |> PrepBool qs n 1 |> ignore            // v register: n bits
            1I |> BoolInt (n+1) |> PrepBool qs (2*n) 1 |> ignore    // s register: n+1 bits
            0I |> BoolInt (n+1) |> PrepBool qs (3*n+1) 1 |> ignore  // r register: n+1 bits
            m  |> BoolInt n |> PrepBool qs (4*n+2) 1 |> ignore      // m register: n bits (the modulus)
            0I |> BoolInt (8*n) |> PrepBool qs (5*n+2) 1 |> ignore  // information to unravel the rounds: 8*n bits
            3I |> BoolInt 2 |> PrepBool qs (13*n+2) 1 |> ignore     // flags for counter 1 and counter 2: 2 bit
            0I |> BoolInt n |> PrepBool qs (13*n+4) 1 |> ignore  // final result: n bits
            initCounter |> BoolInt counterSize |> PrepBool qs (14*n+4) 1 |> ignore // counter: ceil(log(n))+1 bits
            0I |> BoolInt 2 |> PrepBool qs (14*n+4+counterSize) 1 |> ignore  // 2 ancillas for comparator, inplace doubler, and MultiplyCNOTs
            initCounter |> BoolInt counterSize |> PrepBool qs (14*n+6+counterSize) 1 |> ignore // counter: ceil(log(n))+1 bits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier
    
    let arrangeInverterCircuit = 
        match name with 
        | "Montgomery" -> 
            let INV = MontgomeryInverse counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (14*n+6+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    watch.Reset()
    watch.Start()   
    let INV, rslt = arrangeInverterCircuit
    let MontgomeryInverterCircuit = INV qs 
    watch.Stop()
    show "Time to compile %f" watch.Elapsed.TotalSeconds
    
    if verbose1 then 
        show "Number of qubits = %A" (14*n+6+2*counterSize)

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let qs2 = arrangeInverterInputs 
    let initialState = Array.zeroCreate (14*n+6+2*counterSize)
    List.iter (fun i -> M !!(qs2,i)) [0..qs2.Length-1]
    for i in 0..(qs2.Length-1) do 
        initialState.[i] <- qs2.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MontgomeryInverterCircuit initialState
    if verbose2 then 
        show "The initial state is   %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[(13*n+4)..(14*n+3)] 
    show "Final result = %A" res
    let resInt = (IntBool res.Length res)
    show "As a number mod p this is = %A" resInt
    resInt
    

///////////////////////////////////////////////////////////////////////////
//
// Space efficient quantum arithmetic: modular integer inversion with 10n qubits
//
///////////////////////////////////////////////////////////////////////////

// BuildMontgomeryInverseForward implements the inverse mod p of numbers that are given in Montgomery foom x * 2^n mod p. The output is 
// again in Montogomery form. The implementation follows closely the Kalinski paper which in turn is a binary extended Euclidean algorithm. 
// The implementation uses a total of 14+6+2*ceil(log(n)) qubits. The registers are: 
// us = intially equal to modulus p, vs = intitially equal to input x, ss = intially 1, rs = intiially 0, ms = modulus p, mg = 8n qubits 
// holding the results of the dispatch predicates for each round, fs = flag that determines whether we are in compute or in counter mode, 
// ks = a register of qubits to implement a counter from 0 up to worst case iteration bound of 2n, and a = 1 ancilla. 

let BuildEfficientMontgomeryInverseForward (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (mg:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits) = 
    let EuclidRound (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (round:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits)  = 
        // step 1: check termination condition: is vs the all zero string
        X >< vs
        BuildMultiplyControlledNOT vs [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X >< vs
        // step 2: check counter state and flip flag
        BuildMultiplyControlledNOT ks [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X fs
        // step 3: compute match statement to detect next step in binary Euclid
        // step 3a: check whether u is even 
        X us
        // now, use the more efficient encoding: ue -> 01, ve -> 10, ug -> 11, ul -> 00 which can 
        // be computed with 4 Toffoli gates and 1 ancilla which is returned clean after the case statement.
        CCNOT [fs.[0]; us.[0]; round.[0]] // if ue=1 then round.[0] round.[1] = 10
        X us 
        // step 3b: or else, check whether v is even
        X vs 
        X [round.[0]]
        BuildMultiplyControlledNOT [fs.[0]; round.[0]; vs.[0]] [round.[1]] [ss.[0]] // if ve=1 and ue = 0 then round.[0] round.[1] = 01  
        X [round.[0]]
        X vs
        // step 3c: or else, check whether u > v
        BuildCtrlTakahashiAdderInverse vs (us @ a) [fs.[0]] 
        CNOT [round.[0]; a.[1]]
        CNOT [round.[1]; a.[1]] // compute partity of round.[0] round.[1] into ancilla a.[1]
        X [a.[1]] // apply following gate only if ancilla is 0
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [round.[0]] [ss.[0]]
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [round.[1]] [ss.[0]] // if ue=ve=0 and ug=1 then round.[0] round.[1]= 11
        X [a.[1]] 
        CNOT [round.[0]; a.[1]]
        CNOT [round.[1]; a.[1]] // clean up ancilla a.[1] after each round. Note that parity is preseved in all cases        
        // step 3d: else, we have that u <= v
        // note: no computation needed in this case: ue=ve=ug=0 then round.[0] round.[1] = 00
        BuildCtrlTakahashiAdder vs (us @ a) [fs.[0]] // clean up the ancilla a.[0] after each round
        // step 4: dispatch 4 mutually exclusive cases
        // case 4a: u even: conditioning on pattern 10 
        X [round.[1]]
        BuildCtrlBinaryDoubling ss [round.[0]; round.[1]; fs.[0]] [a.[0]] // using a as dirty ancilla
        BuildCtrlBinaryHalfing  us [round.[0]; round.[1]; fs.[0]] [a.[0]]
        X [round.[1]]
        // case 4b: v even: conditioning on pattern 01 
        X [round.[0]]
        BuildCtrlBinaryDoubling rs [round.[0]; round.[1]; fs.[0]] [a.[0]]  
        BuildCtrlBinaryHalfing  vs [round.[0]; round.[1]; fs.[0]] [a.[0]]
        X [round.[0]]
        // case 4c: u odd, v odd, u > v: conditioning on pattern 11
        BuildMultiCtrlTakahashiModAdderInverse us vs [round.[0]; round.[1]; fs.[0]] // replace (u,v) -> (u-v,v)
        BuildCtrlBinaryHalfing  us [round.[0]; round.[1]; fs.[0]] [a.[0]] // replace (u-v) -> (u-v)/2
        BuildMultiCtrlTakahashiModAdder rs ss [round.[0]; round.[1]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling ss [round.[0]; round.[1]; fs.[0]] [a.[0]] // replace s -> 2s        
        // case 4d: u odd, v odd, u <= v: conditioning on pattern 00
        X [round.[0]]
        X [round.[1]]
        BuildMultiCtrlTakahashiModAdderInverse vs us [round.[0]; round.[1]; fs.[0]] // replace (u,v) -> (u,v-u)
        BuildCtrlBinaryHalfing  vs [round.[0]; round.[1]; fs.[0]] [a.[0]] // replace (v-u) -> (v-u)/2
        BuildMultiCtrlTakahashiModAdder ss rs [round.[0]; round.[1]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling rs [round.[0]; round.[1]; fs.[0]] [a.[0]] // replace r -> 2r      
        X [round.[0]]
        X [round.[1]]
        // step 5: increase counter if required
        X fs
        BuildControlShiftRegisterCounterSpecial ks fs 
        X fs
    
    let n = us.Length
    for i in 0..(2*n-1) do 
        let round = (slice mg [2*i..2*i+1]) // get 4 bits per round for match statement on (u, v)
        EuclidRound us vs ss rs round fs a ks 
   
let InverterEfficientModularReduction (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let rs = slice qs [(3*n+1)..(4*n+1)] // n+1 bits for the non-reduced register holding r
    // check if final result in r overflows and if so, then subtract the modulus once
    BuildTakahashiAdderInverse (slice rs [0..(n-1)]) (ms @ [rs.[n]]) 
    BuildCtrlTakahashiModAdder (slice rs [0..(n-1)]) ms [rs.[n]]

let InverterEfficientFixAlmostInverseStandardToMonty (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let ks = slice qs [10*n+4..10*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [10*n+6+counterSize..10*n+5+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [9*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [10*n+4+counterSize..10*n+5+counterSize] // 2 clean ancillas 
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times
    for i in 0..(n-1) do 
         // mapping -x^-1 2^k mod p to -x^-1 2^k 2^-n mod p 
         BuildModularHalver rs ms a         
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -x^-1 2^k 2^-n mod p to -x^-1 2^k 2^-n 2^(2n-k) = -x^-1 2^n mod p 
    for i in 0..(n-1) do 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        BuildCtrlModularDoubler rs ms fs a 
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1 2^n)=x^-1 2^n mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    // BUGBUG: check why this dooesn't work: BuildTakahashiModAdderInverse ms rs 
    X >< ms // one's complement
    BuildTakahashiModAdder rs ms
    X >< ms // return modulus to original state
    X >< rs // one's complement of output   

let InverterEfficientFixAlmostInverseMontyToMonty (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let ks = slice qs [10*n+4..10*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [10*n+6+counterSize..10*n+5+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [9*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [10*n+4+counterSize..10*n+5+counterSize] // 2 clean ancillas 
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times.
    for i in 0..(n-1) do 
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -(x^-1 2^-n) 2^k mod p to -x^-1 2^k 2^-n 2^(2n-k) = -x^-1 2^n mod p 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        BuildCtrlModularDoubler rs ms fs a 
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1 2^n)=x^-1 2^n mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    // BUGBUG: check why this dooesn't work: BuildTakahashiModAdderInverse ms rs 
    X >< ms // one's complement
    BuildTakahashiModAdder rs ms
    X >< ms // return modulus to original state
    X >< rs // one's complement of output   

let InverterEfficientFixAlmostInverseStandardToStandard (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let ks = slice qs [10*n+4..10*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [10*n+6+counterSize..10*n+5+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [9*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [10*n+4+counterSize..10*n+5+counterSize] // 2 clean ancillas 
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times
    for i in 0..(2*n-1) do 
         // mapping -x^-1 2^k mod p to -x^-1 2^k 2^-2n mod p 
         BuildModularHalver rs ms a         
    for i in 0..(n-1) do 
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -x^-1 2^k 2^-2n mod p to -x^-1 2^k 2^-2n 2^(2n-k) = -x^-1 mod p 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        BuildCtrlModularDoubler rs ms fs a 
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1)=x^-1 mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    // BUGBUG: check why this dooesn't work: BuildTakahashiModAdderInverse ms rs 
    X >< ms // one's complement
    BuildTakahashiModAdder rs ms
    X >< ms // return modulus to original state
    X >< rs // one's complement of output   
         
let InverterEfficientPrepCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    let initCounter =  (pown 2I counterSize) - 1I
    
    let ss = slice qs [(2*n)..(3*n)] // prepare initial state on s register: n+1 bits
    let fs = slice qs [(9*n+2)..(9*n+3)] // prepare flags for counter 1 and counter 2: 2 bit
    let cs1 = slice qs [(10*n+4)..(10*n+counterSize+3)]
    let cs2 = slice qs [(10*n+counterSize+6)..(10*n+2*counterSize+5)]
    let us = slice qs [0..(n-1)] // n bits for the register holding the us
    let ms = slice qs [4*n+2..5*n+1] // n bits for the register holding the modulus
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (ms:Qubits) (us:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [ms.[i]; us.[i]]        
    CopyRegs ms us 
    X ss; X >< fs; X >< cs1; X >< cs2 // set first bit in ss and all remaining qubits in fs, cs1, and cs2 to 1
    
let InverterEfficientCopyCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let rs = slice qs [(3*n+1)..(4*n)] // n bits for the reduced register holding r
    let res = slice qs [9*n+4..10*n+3] // n bits for the final result
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (rs:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [rs.[i]; res.[i]]        
    CopyRegs rs res
 
let MontgomeryEfficientInverseForward (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-6)/10
    let us = slice qs [0..(n-1)] // u (n bits)
    let vs = slice qs [n..(2*n-1)] // v (n bits)
    let ss = slice qs [(2*n)..(3*n)] // s (n+1 bits)
    let rs = slice qs [(3*n+1)..(4*n+1)] // r (n+1 bits); holds the almost MG inverse
    let mg = slice qs [(5*n+2)..(9*n+1)] // round status indicators
    let fs = slice qs [9*n+2] // flag needed for 1st while loop
    let a  = slice qs [10*n+4+counterSize..10*n+5+counterSize] // need 2 ancillas
    let ks = slice qs [10*n+4..10*n+3+counterSize] // shift register counter 
    BuildEfficientMontgomeryInverseForward us vs ss rs mg fs a ks
 
let MontgomeryEfficientInverse (name:string) (counterSize:int) (verbose:bool) (qs:Qubits) =     
    let watch = Diagnostics.Stopwatch() 
    
    watch.Reset()
    watch.Start()   
    let a = Circuit.Compile (InverterEfficientPrepCircuit counterSize) qs
    watch.Stop()
    show "INVERTER: Time for step a: %f" watch.Elapsed.TotalSeconds

    watch.Reset()
    watch.Start()   
    let b = Circuit.Compile (MontgomeryEfficientInverseForward counterSize) qs
    watch.Stop()
    show "INVERTER: Time for step b: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let c = Circuit.Compile (InverterEfficientModularReduction counterSize) qs
    watch.Stop()
    show "INVERTER: Time for step c: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let d = match name with 
            | "MontgomerySM" -> Circuit.Compile (InverterEfficientFixAlmostInverseStandardToMonty counterSize) qs
            | "MontgomeryMM" -> Circuit.Compile (InverterEfficientFixAlmostInverseMontyToMonty counterSize) qs
            | "MontgomerySS" -> Circuit.Compile (InverterEfficientFixAlmostInverseStandardToStandard counterSize) qs
            | _ -> failwith "MG inverse: this shoud never happen"             
    watch.Stop()
    show "INVERTER: Time for step d: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let e = Circuit.Compile (InverterEfficientCopyCircuit counterSize) qs
    watch.Stop()
    show "INVERTER: Time for step e: %f" watch.Elapsed.TotalSeconds
    
    let f = d.Reverse() // Run the correction to fix almost inverse in reverse
    watch.Stop()
    show "INVERTER: Time for step f: %f" watch.Elapsed.TotalSeconds
    
    let g = c.Reverse() // Run the final modular reduction in reverse
    watch.Stop()
    show "INVERTER: Time for step g: %f" watch.Elapsed.TotalSeconds
    
    let h = b.Reverse() // Run the sequence of Montgomery steps in reverse
    watch.Stop()
    show "INVERTER: Time for step h: %f" watch.Elapsed.TotalSeconds
    
    let i = a.Reverse() // Run the sequence of state preparation of the us in reverse
    watch.Stop()
    show "INVERTER: Time for step i: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let gates = 
        [a; b; c; d; e; f; g; h; i] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    watch.Stop()
    show "INVERTER: Time to contatenate: %f" watch.Elapsed.TotalSeconds
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let MontgomeryEfficientInverseManualReverse (name:string) (counterSize:int) (verbose:bool) (qs:Qubits) =     
    let watch = Diagnostics.Stopwatch() 
    
    watch.Reset()
    watch.Start()   
    let a = Circuit.Compile (InverterEfficientPrepCircuit counterSize) qs |> MyCircuitExport
    watch.Stop()
    show "INVERTER: Time for step a: %f" watch.Elapsed.TotalSeconds

    watch.Reset()
    watch.Start()   
    let b = Circuit.Compile (MontgomeryEfficientInverseForward counterSize) qs |> MyCircuitExport
    watch.Stop()
    show "INVERTER: Time for step b: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let c = Circuit.Compile (InverterEfficientModularReduction counterSize) qs |> MyCircuitExport
    watch.Stop()
    show "INVERTER: Time for step c: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let d = match name with 
            | "MontgomerySM" -> Circuit.Compile (InverterEfficientFixAlmostInverseStandardToMonty counterSize) qs 
            | "MontgomeryMM" -> Circuit.Compile (InverterEfficientFixAlmostInverseMontyToMonty counterSize) qs
            | "MontgomerySS" -> Circuit.Compile (InverterEfficientFixAlmostInverseStandardToStandard counterSize) qs
            | _ -> failwith "MG inverse: this shoud never happen"          
            |> MyCircuitExport   
    watch.Stop()
    show "INVERTER: Time for step d: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let e = Circuit.Compile (InverterEfficientCopyCircuit counterSize) qs |> MyCircuitExport
    watch.Stop()
    show "INVERTER: Time for step e: %f" watch.Elapsed.TotalSeconds
    
    let f = List.rev d // Run the correction to fix almost inverse in reverse
    watch.Stop()
    show "INVERTER: Time for step f: %f" watch.Elapsed.TotalSeconds
    
    let g = List.rev c // Run the final modular reduction in reverse
    watch.Stop()
    show "INVERTER: Time for step g: %f" watch.Elapsed.TotalSeconds
    
    let h = List.rev b // Run the sequence of Montgomery steps in reverse
    watch.Stop()
    show "INVERTER: Time for step h: %f" watch.Elapsed.TotalSeconds
    
    let i = List.rev a // Run the sequence of state preparation of the us in reverse
    watch.Stop()
    show "INVERTER: Time for step i: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let gates = 
        [a; b; c; d; e; f; g; h; i] 
        |> List.concat
    watch.Stop()
    show "INVERTER: Time to contatenate: %f" watch.Elapsed.TotalSeconds
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let RunEfficientInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    
    // Prepare the initial state. The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | modulus | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | n       | 4n        | 2     | n      | log(n)   | 2        | log(n)   ||
    // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 00       | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let k = Ket(10*n+6+2*counterSize)
            let qs = k.Qubits
            0I  |> BoolInt n |> PrepBool qs 0 1 |> ignore            // u register: n bits
            x  |> BoolInt n |> PrepBool qs n 1 |> ignore            // v register: n bits
            0I |> BoolInt (n+1) |> PrepBool qs (2*n) 1 |> ignore    // s register: n+1 bits
            0I |> BoolInt (n+1) |> PrepBool qs (3*n+1) 1 |> ignore  // r register: n+1 bits
            m  |> BoolInt n |> PrepBool qs (4*n+2) 1 |> ignore      // m register: n bits (the modulus)
            0I |> BoolInt (4*n) |> PrepBool qs (5*n+2) 1 |> ignore  // information to unravel the compressed rounds: 4*n bits
            0I |> BoolInt 2 |> PrepBool qs (9*n+2) 1 |> ignore     // flags for counter 1 and counter 2: 2 bit
            0I |> BoolInt n |> PrepBool qs (9*n+4) 1 |> ignore  // final result: n bits
            0I |> BoolInt counterSize |> PrepBool qs (10*n+4) 1 |> ignore // counter: ceil(log(n))+1 bits
            0I |> BoolInt 2 |> PrepBool qs (10*n+4+counterSize) 1 |> ignore  // 2 ancillas for comparator, inplace doubler, and MultiplyCNOTs
            0I |> BoolInt counterSize |> PrepBool qs (10*n+6+counterSize) 1 |> ignore // counter: ceil(log(n))+1 bits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier
    
    let arrangeInverterCircuit = 
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let INV = MontgomeryEfficientInverse name counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (10*n+6+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    watch.Reset()
    watch.Start()   
    let INV, rslt = arrangeInverterCircuit
    let MontgomeryEfficientInverterCircuit = INV qs 
    watch.Stop()
    show "Time to compile %f" watch.Elapsed.TotalSeconds
    
    if verbose1 then 
        show "Number of qubits = %A" (10*n+6+2*counterSize)

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let qs2 = arrangeInverterInputs 
    let initialState = Array.zeroCreate (10*n+6+2*counterSize)
    List.iter (fun i -> M !!(qs2,i)) [0..qs2.Length-1]
    for i in 0..(qs2.Length-1) do 
        initialState.[i] <- qs2.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MontgomeryEfficientInverterCircuit initialState
    if verbose2 then 
        show "The initial state is %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[(9*n+4)..(10*n+3)] 
    show "Final result = %A" res
    let resInt = (IntBool res.Length res)
    show "As a number mod p this is = %A" resInt
    resInt

let CompileEfficientModularInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    let initCounter =  (pown 2I counterSize) - 1I

    // The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | modulus | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | n       | 4n        | 2     | n      | log(n)   | 2        | log(n)   ||
    // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 00       | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "MontgomerySM" | "MontgomeryMM" | "MontgomerySS" -> 
            let k = Ket(10*n+6+2*counterSize)
            let qs = k.Qubits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier

    if verbose2 then 
        show "Number of qubits = %A" (10*n+6+2*counterSize)

    let arrangeInverterCircuit = 
        match name with 
        | "MontgomerySM" | "MontgomeryMM" | "MontgomerySS" -> 
            let INV = MontgomeryEfficientInverseManualReverse name counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (10*n+6+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    let INV, rslt = arrangeInverterCircuit
    let MontgomeryEfficientInverterCircuit = INV qs 
    
    MontgomeryEfficientInverterCircuit


///////////////////////////////////////////////////////////////////////////
//
// Tools for stitching of Toffoli gate lists into larger circuits 
//
///////////////////////////////////////////////////////////////////////////

let MapToffoli (T : TofGate)  (embed:Map<int,int>) = 
    match T with 
    | MyNOT(a) ->  MyNOT(embed.[a]) 
    | MyCNOT(a,b) -> MyCNOT(embed.[a], embed.[b]) 
    | MyTOFF(a,b,c) -> MyTOFF(embed.[a], embed.[b], embed.[c]) 
  
let RemapToffoli (T : TofGate)  (embed:Dictionary<int,int>) = 
    match T with 
    | MyNOT(a) ->  MyNOT(embed.[a]) 
    | MyCNOT(a,b) -> MyCNOT(embed.[a], embed.[b]) 
    | MyTOFF(a,b,c) -> MyTOFF(embed.[a], embed.[b], embed.[c]) 
  
let MapToffoliCircuit (T: TofGate list) (embed:Map<int,int>) = 
    T |> List.map (fun x -> MapToffoli x embed)

let RemapToffoliCircuit (T: TofGate list) (embed:Dictionary<int,int>) = 
    T |> List.map (fun x -> RemapToffoli x embed)

let difference (L1:list<'a>) (L2:list<'a>) =
        let cache = HashSet<'a>(L2, HashIdentity.Structural)
        L1 |> List.filter (fun x -> not (cache.Contains x))

let rec concat (L:list<list<'a>>) =
    match L with 
    | [] -> []
    | x::[] -> x
    | x::xs -> List.append x (concat xs)


///////////////////////////////////////////////////////////////////////////
//
// Application of modular arithmetic: Shor for factoring
//
///////////////////////////////////////////////////////////////////////////

let factoring (name:string) (n:int) (N:bigint) verbosity = 
    let verbose1 = (verbosity = "verbose" || verbosity = "all")
    let verbose2 = (verbosity = "all")

    if verbose2 then printf "Creating circuit: "
    let watch = Diagnostics.Stopwatch()        
    watch.Start()
    let MUL = 
        match name with 
        | "Montgomery" -> 
            let MULmon = CompileModularMultiplier "CtrlMontgomery" n 1I 1I N verbosity
            if verbose2 then printf "MUL MM done. "
            watch.Stop()
            show "CSV: Time to construct MUL MM %f" watch.Elapsed.TotalSeconds            
            MULmon 
        | "DoubleAdd"  ->
            let MULdadd = CompileModularMultiplier "CtrlDoubleAdd" n 1I 1I N verbosity
            if verbose2 then printf "MUL DA done. "
            watch.Stop()
            show "CSV: Time to construct MUL DA %f" watch.Elapsed.TotalSeconds               
            MULdadd
        | _ -> failwith "unsupported multiplier type for Shor's factoring algorithm"
    watch.Reset()
    MUL 

let RunShorFactoring (name:string) (n:int) (N:bigint) (x:bigint) (m:bigint) (ctrl:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states and other info
    
    let Ns = BoolInt n N;  
    let cs = BoolInt n ctrl;
    let xs = BoolInt n x; 
    let ms = BoolInt n m; 
            
    // Prepare the initial state
    let mutable factoring_initialstate =   
        match name with 
        | "Montgomery" -> 
            let init = Array.zeroCreate (6*n+5) // total #qubits = 6n+5
            // name: || xs | ms | ancilla | acc  | modulus | ancilla | mg-rounds | result | ctrl || 
            // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      | 1    ||
            // init: || x  | m  | 0       | 0    | N       | 0       | 0..0      | 0      | ctrl || 
            for i in [0..(n-1)] do 
                init.[i]       <- xs.[i]
                init.[n+i]     <- ms.[i]
                init.[3*n+3+i] <- Ns.[i]            
            init.[(6*n+4)]     <- cs.[0]
            init
        | "DoubleAdd" -> 
            let init = Array.zeroCreate (4*n+5) // total #qubits = 4n+5
            // Prepare the initial state for double&add. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | acc  | modulus | anc | ctrl || 
            // bits: || n  | n  | n+1  | n       | 3   | 1    ||
            // init: || x  | m  | 0    | N       | 000 | ctrl || 
            // note that the final result is stored in the lower n bits of acc. 
            for i in [0..(n-1)] do 
                init.[i]       <- xs.[i]
                init.[n+i]     <- ms.[i]
                init.[3*n+1+i] <- Ns.[i]            
            init.[(4*n+4)]     <- cs.[0]
            init
    
    let factoring_circuit = // constructing the factoring circuit
        factoring name n N verbosity
        
    if verbose2 then 
        show "Number of qubits = %A" factoring_initialstate.Length
        factoring_circuit  |> List.filter (fun x -> match x with | MyTOFF(a,b,c) -> true | _ -> false) 
                           |> List.length
                           |> show "Number of Toffoli gates = %A" 
        
    let watch = Diagnostics.Stopwatch()        
    watch.Start()
    let mutable factoring_finalstate = 
        match name with 
        | "Montgomery" -> Array.zeroCreate (6*n+5) 
        | "DoubleAdd"  -> Array.zeroCreate (4*n+5) 
    
    factoring_finalstate <- MyCircuitSimulateFast factoring_circuit factoring_initialstate
    watch.Stop()
    show "CSV: Time to simulate with fast simulator %f" watch.Elapsed.TotalSeconds            
    
    if verbose2 then 
        show "The initial state is   %A" factoring_initialstate
        show "And the final state is %A" factoring_finalstate
                
    let xdump = factoring_finalstate.[0..(n-1)] |> IntBool n
    let mdump = factoring_finalstate.[n..(2*n-1)] |> IntBool n 
    let Ndump =
        match name with 
        | "Montgomery" -> factoring_finalstate.[(3*n+3)..(4*n+2)] |> IntBool n
        | "DoubleAdd"  -> factoring_finalstate.[(3*n+1)..(4*n)] |> IntBool n
    let resdump = 
        match name with 
        | "Montgomery" -> factoring_finalstate.[(5*n+4)..(6*n+3)] |> IntBool n
        | "DoubleAdd"  -> factoring_finalstate.[(2*n)..(3*n-1)] |> IntBool n
    if verbose2 then 
        show "x=%A m=%A N=%A res=%A" xdump mdump Ndump resdump 
            
    resdump

///////////////////////////////////////////////////////////////////////////
//
// Application of modular arithmetic: Shor for dlog in ECC curves
//
///////////////////////////////////////////////////////////////////////////

let ec_add_affine (n:int) (p:bigint) verbosity = 
    let verbose1 = (verbosity = "verbose" || verbosity = "all")
    let verbose2 = (verbosity = "all")
    let m = double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 // counter size
    
    // Construction of the point addition formula for ECC in affine Weierstrass form, 
    // first point P=(X1,Y1), second point Q=(x2,y2), where Q is classically known 
    // and P is overwritten by P+Q. The entire point addition is controlled by a bit <ctrl>. 
    let qubits = [|0..(14*n+2*m+6)|] // total #qubits = 14n+2m+7
        // name: || ctrl | p  | X1 | Y1 | x2 | y2 | t0 | lam | ancillas || 
        // bits: || 1    | n  | n  | n  | n  | n  | n  | n   | 7n+2m+6 ||
        // maximum ancillas required: n + (n+1) + (n+1) + 4n + 2 + log(n) + 2 + log(n) = 7n + 2m + 6
    let ctrl = qubits.[0..0]
    let ps   = qubits.[1..n]
    let X1   = qubits.[(n+1)..(2*n)]
    let Y1   = qubits.[(2*n+1)..(3*n)]
    let x2   = qubits.[(3*n+1)..(4*n)]
    let y2   = qubits.[(4*n+1)..(5*n)]
    let t0   = qubits.[(5*n+1)..(6*n)]
    let lam  = qubits.[(6*n+1)..(7*n)]
    let anc  = qubits.[(7*n+1)..(14*n+2*m+6)] // BUGBUG : check this!!!

    // Contruct all the arithmetic functions needed to implement the ECC point additions reversibly. 
    // First the circuits are constructed as Toffoli networks operating on qubits in a range [0..(R-1)], 
    // where R is the largest qubit index used by the gates. Then these circuits are embedded using a map. 
    
    let watch = Diagnostics.Stopwatch()        
    
    if verbose2 then printf "Creating circuits: "
    
    watch.Start()
    let NEG = CompileModularAdder "CtrlNeg" n 1I 1I p 1I verbosity
    if verbose2 then printf "cNEG done. "
    watch.Stop()
    show "CSV: Time to construct cNEG %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let ADD = CompileModularAdder "ModADD" n 1I 1I p 1I verbosity
    if verbose2 then printf "ADD done. "
    watch.Stop()
    show "CSV: Time to construct ADD %f" watch.Elapsed.TotalSeconds            
    watch.Reset()
       
    watch.Start()
    let cADD = CompileModularAdder "CtrlModADD" n 1I 1I p 1I verbosity
    if verbose2 then printf "cADD done. "
    watch.Stop()
    show "CSV: Time to construct cADD %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let SUB  = List.rev ADD
    if verbose2 then printf "SUB done. "
    watch.Stop()
    show "CSV: Time to construct SUB %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let cSUB = List.rev cADD
    if verbose2 then printf "cSUB done. "
    watch.Stop()
    show "CSV: Time to construct sSUB %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let MUL = CompileModularMultiplier "Montgomery" n 1I 1I p verbosity
    if verbose2 then printf "MUL MM done. "
    watch.Stop()
    show "CSV: Time to construct MUL MM %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let MUL2 = CompileModularMultiplier "Mersenne" n 1I 1I p verbosity
    if verbose2 then printf "MUL Mers done. "
    watch.Stop()
    show "CSV: Time to construct MUL Mers %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let MUL3 = CompileModularMultiplier "DoubleAdd" n 1I 1I p verbosity
    if verbose2 then printf "MUL DA done. "
    watch.Stop()
    show "CSV: Time to construct MUL DA %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let MUL4 = CompileModularMultiplier "DoubleAddMersenne" n 1I 1I p verbosity
    if verbose2 then printf "MUL DAM done. "
    watch.Stop()
    show "CSV: Time to construct MUL DAM %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let SQU = CompileModularMultiplier "MontgomerySquarer" n 1I 1I p verbosity
    if verbose2 then printf "SQU done. "
    watch.Stop()
    show "CSV: Time to construct SQU: %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let INV = CompileEfficientModularInverter "MontgomeryMM" n p 1I verbosity
    if verbose2 then printf "INV done.\n"
    watch.Stop()
    show "CSV: Time to construct INV: %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    // Mapping step: in the following register allotments, the following conventions are used: 
    // "-" stands for an input that is returned unchanged, "*" stands for an input that is overwritten and 
    // carries the output all other bits and registers are considered ancillas.

    let ctrl_neg_modp (ps:int []) (xs:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 2n+1 qubits
        // name: ||*xs*|-modulus-|-ctrl-||
        // bits: || n  | n       | 1    ||
        // init: || x  | p       | c    ||   
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, xs.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, ps.[i]) ] @
//                          [ for i in [0..0] do yield (2*n+i, cs.[i]) ] )
//        MapToffoliCircuit NEG embed
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]
            embed.[n+i] <- ps.[i]
        embed.[2*n] <- cs.[0]
        RemapToffoliCircuit NEG embed

    let add_modp (ps:int []) (xs:int []) (ys:int []) = // result is stored in xs
        // circuit operates nontrivially on 3n+2 qubits
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp ||
        // bits: || n  | n  | 1   | n       | 1   ||
        // init: || x  | y  | 0   | p       | 0   ||
        let ds = concat [ [2*n]; [3*n+1] ] // pick up all ancillas       
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, ys.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, xs.[i]) ] @ // note: addressing chosen so that result stored in xs
//                          [ for i in [0..(n-1)] do yield (2*n+1+i, ps.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )                      
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]
            embed.[2*n+1+i] <- ps.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit ADD embed

    let ctrl_add_modp (ps:int []) (xs:int []) (ys:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 3n+4 qubits
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        let ds = concat [ [2*n]; [(3*n+1)..(3*n+2)] ] // pick up all ancillas
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, ys.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, xs.[i]) ] @ // note: addressing chosen so that result stored in xs
//                          [ for i in [0..(n-1)] do yield (2*n+1+i, ps.[i]) ] @
//                          [ for i in [0..0] do yield (3*n+3, cs.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )                      
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]
            embed.[2*n+1+i] <- ps.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[3*n+3] <- cs.[0]
        RemapToffoliCircuit cADD embed

    let sub_modp (ps:int []) (xs:int []) (ys:int []) = // result is stored in xs
        // circuit operates nontrivially on 3n+2 qubits
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp ||
        // bits: || n  | n  | 1   | n       | 1   ||
        // init: || x  | y  | 0   | p       | 0   ||
        let ds = concat [ [2*n]; [3*n+1] ] // pick up all ancillas
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, ys.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, xs.[i]) ] @ // note: addressing chosen so that result stored in xs
//                          [ for i in [0..(n-1)] do yield (2*n+1+i, ps.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )                      
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]
            embed.[2*n+1+i] <- ps.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit SUB embed

    let ctrl_sub_modp (ps:int []) (xs:int []) (ys:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 3n+4 qubits
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        let ds = concat [ [2*n]; [(3*n+1)..(3*n+2)] ] // pick up all ancillas
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, ys.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, xs.[i]) ] @ // note: addressing chosen so that result stored in xs
//                          [ for i in [0..(n-1)] do yield (2*n+1+i, ps.[i]) ] @
//                          [ for i in [0..0] do yield (3*n+3, cs.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )    
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]
            embed.[2*n+1+i] <- ps.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[3*n+3] <- cs.[0]
        RemapToffoliCircuit cSUB embed

    let mul_modp (ps:int []) (xs:int []) (ys:int []) (res:int []) = // result is stored in res
        // circuit operates nontrivially on 6n+4 qubits
        // name: ||-xs-|-ys-| ancilla | acc  |-modulus-| ancilla | mg-rounds |*result*|| 
        // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      ||
        // init: || x  | y  | 0       | 0    | p       | 0       | 0..0      | 0      || 
        let ds = concat [ [2*n]; [2*n..(3*n+2)]; [4*n+3]; [4*n+4..5*n+3] ] // pick up all ancillas
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, xs.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (n+i, ys.[i]) ] @      
//                          [ for i in [0..(n-1)] do yield (5*n+4+i, res.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (3*n+3+i, ps.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )                      
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]
            embed.[n+i] <- ys.[i]
            embed.[5*n+4+i] <- res.[i]
            embed.[3*n+3+i] <- ps.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit MUL embed

    let squ_modp (ps:int []) (xs:int []) (res:int []) = // result is stored in res
        // circuit operates nontrivially on 5n+5 qubits
        // name: ||-xs-| ancilla | acc  |-modulus-| ancilla | mg-rounds |*result*|| 
        // bits: || n  | 2       | n+2  | n       | 1       | n         | n      ||
        // init: || x  | 0       | 0    | p       | 0       | 0..0      | 0      ||         
        let ds = concat [ [n..(n+1)]; [(n+2)..(2*n+3)]; [3*n+4]; [(3*n+5)..(4*n+4)] ] // pick up all ancillas
//        let embed = Map ( [ for i in [0..(n-1)] do yield (i, xs.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (4*n+5+i, res.[i]) ] @
//                          [ for i in [0..(n-1)] do yield (2*n+4+i, ps.[i]) ] @
//                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )                      
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i] 
            embed.[4*n+5+i] <- res.[i] 
            embed.[2*n+4+i] <- ps.[i] 
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]

        RemapToffoliCircuit SQU embed
        
    let inv_modp (ps:int []) (xs:int []) (res:int []) =  // result is stored in res
        // circuit operates nontrivially on 10n+2m+6 qubits
        // name: || us |-vs-| ss  | rs  |-modulus-| mg-rounds | flags |*result*| counter1 | ancillas | counter2 || 
        // bits: || n  | n  | n+1 | n+1 | n       | 4n        | 2     | n      | log(n)   | 2        | log(n)   ||
        // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 00       | 1..1     || 
        let ds = concat [ [0..(n-1)]; [(2*n)..(3*n)]; [(3*n+1)..(4*n+1)]; [(5*n+2)..(9*n+1)]; [(9*n+2)..(9*n+3)]; 
                          [(10*n+4)..(10*n+m+3)]; [(10*n+m+4)..(10*n+m+5)]; [(10*n+m+6)..(10*n+2*m+5)] ] // pick up all ancillas
        let embedOld = Map ( [ for i in [0..(n-1)] do yield (n+i, xs.[i]) ] @
                          [ for i in [0..(n-1)] do yield (4*n+2+i, ps.[i]) ] @
                          [ for i in [0..(n-1)] do yield (9*n+4+i, res.[i]) ] @
                          [ for i in [0..(ds.Length-1)] do yield (ds.[i], anc.[i]) ] )          
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[n+i] <- xs.[i]
            embed.[4*n+2+i] <- ps.[i]
            embed.[9*n+4+i] <- res.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]

//        show "I'm here again %A" xs 
//        //List.iter (fun x -> printf "%A" (x, embed.[x])) [0..(10*n+2*m+5)]; printf "\n"
//        watch.Start()
//        MapToffoliCircuit INV embedOld |> ignore
//        watch.Stop()
//        show "CSV: Time to remap construct INV (old): %f" watch.Elapsed.TotalSeconds            
//        watch.Reset()
//        
        watch.Start()
        let a = RemapToffoliCircuit INV embed    
        watch.Stop()
        show "CSV: Time to remap construct INV: %f" watch.Elapsed.TotalSeconds            
        watch.Reset()
        a
    
    let opList = [
        sub_modp ps X1 x2;              // x1 <- x1 - x2
        ctrl_sub_modp ps Y1 y2 ctrl;    // y1 <- y1 - y2    // CTRL 
        inv_modp ps X1 t0;              // t0 <- 1/(x1 - x2) 
        mul_modp ps Y1 t0 lam;          // lam <- (y1 - y2)/(x1 - x2)
        mul_modp ps lam X1 Y1;          // y1 <- 0
        inv_modp ps X1 t0;              // t0 <- 0 
        squ_modp ps lam t0;             // t0 <- ((y1 - y2)/(x1 - x2))^2 = lam^2
        ctrl_sub_modp ps X1 t0 ctrl;    // x1 <- x1 - x2 - ((y1 - y2)/(x1 - x2))^2 // CTRL
        ctrl_add_modp ps X1 x2 ctrl;    // CTRL
        ctrl_add_modp ps X1 x2 ctrl;    // CTRL
        ctrl_add_modp ps X1 x2 ctrl;    // x1 <- x1 - x2 - ((y1 - y2)/(x1 - x2))^2 + 3*x2 = -(x3 - x2) // CTRL
        squ_modp ps lam t0;             // t0 <- 0
        mul_modp ps lam X1 Y1;          // y1 <- y3 + y2
        inv_modp ps X1 t0;              // t0 <- -1/(x3 - x2)
        mul_modp ps t0 Y1 lam;          // lam <- 0
        inv_modp ps X1 t0;              // t0 <- 0
        ctrl_neg_modp ps X1 ctrl;       // x1 <- x3 - x2 // CTRL
        ctrl_sub_modp ps Y1 y2 ctrl;    // y3 done // CTRL
        add_modp ps X1 x2;              // x3 done
    ]
    // output basic profiling of number of Toffoli gates per each instruction
    if verbose2 then 
        printf "Toffoli gates: "
        opList |> List.iter (fun x -> printf "%A " (x |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length) )
        printf "\n"

//    watch.Start()
//    let a = opList |> List.concat
//    watch.Stop()
//    show "CSV: Time to put the concat together: %f" watch.Elapsed.TotalSeconds            
//    watch.Reset()
            
    //opList |> List.concat
    opList 

let RunShorEllipticCurvePointAddition (name:string) (n:int) (p:bigint) (P:eccpoint) (Q:eccpoint) (ctrl:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states and other info
    
    let m = double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 // counter size
   
    // Prepare the initial state
    let mutable ecc_initialstate =   
        match name with 
        | "AffineWeierstrass" -> 
            let init = Array.zeroCreate (14*n+2*m+7) // total #qubits = 14n+2m+7
            // name: || ctrl | p  | X1 | Y1 | x2 | y2 | t0 | lam | ancillas || 
            // bits: || 1    | n  | n  | n  | n  | n  | n  | n   | 7n+2m+6  ||
            // maximum ancillas required: n + (n+1) + (n+1) + 4n + 2 + log(n) + 2 + log(n) = 7n + 2m + 6
            match P, Q with 
            | AffinePoint(x1,y1), AffinePoint(x2,y2) ->
                let ps = BoolInt n p;  
                let cs = BoolInt n ctrl;
                let X1 = BoolInt n x1; 
                let Y1 = BoolInt n y1;
                let X2 = BoolInt n x2; 
                let Y2 = BoolInt n y2;
                init.[0] <- cs.[0]
                for i in [0..(n-1)] do 
                    init.[i+1]     <- ps.[i]
                    init.[n+i+1]   <- X1.[i]
                    init.[2*n+i+1] <- Y1.[i]
                    init.[3*n+i+1] <- X2.[i]
                    init.[4*n+i+1] <- Y2.[i]  
            | _ -> failwith "wrong ECC point format; need P, Q both to be affine." 
            init
        | _ -> failwith "unknown ECC type."
   
    let ecc_circuit = // constructing the multiplier circuit
        match name with 
        | "AffineWeierstrass" -> ec_add_affine n p verbosity
        | _ -> failwith "unknown ECC type"
        
//    if verbose2 then 
//        show "Number of qubits = %A" ecc_initialstate.Length
//        ecc_circuit  |> List.filter (fun x -> match x with | MyTOFF(a,b,c) -> true | _ -> false) 
//                     |> List.length
//                     |> show "Number of Toffoli gates = %A" 
// Rewrite this to make it compatible with piecewise networks
    
    // Measure timings of slow vs fast simulator
//    let watch = Diagnostics.Stopwatch()        
//    if verbose2 then 
//        watch.Start()
//        let ecc_finalstateOld = MyCircuitSimulateFast ecc_circuit ecc_initialstate
//        watch.Stop()
//        show "CSV: Time to simulate with slow simulator %f" watch.Elapsed.TotalSeconds            
//        watch.Reset()
//    watch.Start()
//    let ecc_finalstate = MyCircuitSimulateFast ecc_circuit ecc_initialstate
//    watch.Stop()
//    show "CSV: Time to simulate with fast simulator %f" watch.Elapsed.TotalSeconds            
    
    let watch = Diagnostics.Stopwatch()        
    watch.Start()
    let mutable ecc_finalstate = Array.zeroCreate (14*n+2*m+7) // total #qubits = 14n+2m+7
    for op in ecc_circuit do 
        ecc_finalstate <-  MyCircuitSimulateFast op ecc_initialstate
        ecc_initialstate <- ecc_finalstate
    watch.Stop()
    show "CSV: Time to simulate with fast simulator %f" watch.Elapsed.TotalSeconds            
    
    if verbose2 then 
        show "The initial state is   %A" ecc_initialstate
        show "And the final state is %A" ecc_finalstate        
        let pdump = ecc_finalstate.[1..n] |> IntBool n
        let x1dump = ecc_finalstate.[(n+1)..(2*n)] |> IntBool n
        let y1dump = ecc_finalstate.[(2*n+1)..(3*n)] |> IntBool n
        let x2dump = ecc_finalstate.[(3*n+1)..(4*n)] |> IntBool n
        let y2dump = ecc_finalstate.[(4*n+1)..(5*n)] |> IntBool n
        let t0dump = ecc_finalstate.[(5*n+1)..(6*n)] |> IntBool n
        let lamdump = ecc_finalstate.[(6*n+1)..(7*n)] |> IntBool n
        let ancdump = ecc_finalstate.[(7*n+1)..(14*n+2*m+6)] |> IntBool n
        show "p=%A x1=%A y1=%A x2=%A y2=%A t0=%A lam=%A anc=%A" pdump x1dump y1dump x2dump y2dump t0dump lamdump ancdump

    let res = 
        match name with
        | "AffineWeierstrass" -> 
            AffinePoint( (IntBool n ecc_finalstate.[(n+1)..(2*n)]), (IntBool n ecc_finalstate.[(2*n+1)..(3*n)]))
        | _ -> failwith "unknown ECC type."
    res

///////////////////////////////////////////////////////////////////////////
//
// Some testing of adder and the Montgomery multiplier/inversion modules
//
///////////////////////////////////////////////////////////////////////////

[<LQD>]
let __TestHexConverter () = 
    let b521hex = "5ac635d8 aa3a93e7 b3ebbd55 769886bc 651d06b0 cc53b0f6 3bce3c3e 27d2604b"
    let b521dec = BigintHex b521hex 
    printf "Hex = %A\n" b521hex 
    printf "Dec = %A\n" b521dec
    
    let b521hex2 = HexBigint b521dec    
    printf "Now the other way around\n"
    printf "Dec = %A\n" b521dec    
    printf "Hex = %A\n" b521hex2
    show "done."

[<LQD>]
let __TestFastCircuitSimulate() = 
    let gates = [ MyCNOT(0,1); MyNOT(2); MyTOFF(0,2,1); MyTOFF(10,9,2); MyCNOT(8,4); MyCNOT(2,7); MyTOFF(7,1,2); MyNOT(5); MyNOT(4)] 
    let input = [| 1; 1; 1; 0; 1; 0; 0; 1; 1; 1; 1 |]
    let b = BitVector(11)
    b.Init(input)
    b.PrettyPrint()
    b.Run(gates)
    b.PrettyPrint()

[<LQD>]
let __TestSlowCircuitSimulate() =     
    let gates = [ MyCNOT(0,1); MyNOT(2); MyTOFF(0,2,1); MyTOFF(10,9,2); MyCNOT(8,4); MyCNOT(2,7); MyTOFF(7,1,2); MyNOT(5); MyNOT(4)] 
    let input = [| 1; 1; 1; 0; 1; 0; 0; 1; 1; 1; 1 |]
    let b = MyCircuitSimulateSlow gates input
    for i in 0..(b.Length-1) do printf "%A" b.[i]
    printf "\n"

[<LQD>]
let __RunSmallAdderTests () = 
    RunAdder "Cuccaro" 4 3I 5I 1I
    RunAdder "Ripple" 4 3I 5I 1I
    RunAdder "Takahashi" 4 3I 5I 1I
    RunAdder "Ctrl-Takahashi" 4 3I 5I 0I
    RunAdder "Ctrl-Takahashi" 4 3I 5I 1I

[<LQD>]
let __TestMultiplyNOT (n:int) (s:int) =
    // Prepare the initial state
    let k = Ket(n+2)
    let qs = k.Qubits
    s |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0I |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    1I |> BoolInt 1 |> PrepBool qs (n+1) 1 |> ignore
    
    let MultiplyControlledNOTCircuit = MultiplyControlledNOTgate qs

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+2)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MultiplyControlledNOTCircuit initialState
    printf "The initial state is   %A\n" initialState
    printf "And the final state is %A\n" finalState
    show "done"

let TestComparator (n:int) (s1:int) (s2:int) (c:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let b  = double (c+1) |> log |> fun x -> x/(log 2.0) |> ceil |> int 
    let k = Ket(2*n+1+b)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    s2 |> fun x -> (bigint x) |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
    c  |> fun x -> (bigint x) |> BoolInt b |> PrepBool qs (2*n+1) 1 |> ignore
    
    let ComparatorCircuit = Circuit.Compile CtrlStrictlyLargerThanComparator qs
                            |> MyCircuitExport
                       
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (2*n+1+b)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast ComparatorCircuit initialState
    (finalState.[2*n] = 1)

[<LQD>]
let __TestComparatorAllInputs (n:int) (s:int) (c:int) = 
    for i in 0..(pown 2 n)-1 do 
        let res = TestComparator n s i c
        show "Is %i > %i? : %b" s i res 

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
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running modular adder circuit test for %A + %A mod %A" s1 s2 p
        show "====================================================="
    let result =  RunModularAdder "ModADD" n s1 s2 p 1I verbosity
    if (s1 + s2) % p = result then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result 
        0 // return integer code zero in case of failure

let RunModularNegatorTest n s p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running controlled modular negator circuit test for - %A mod %A" s p
        show "===================================================================="
    let result =  RunModularAdder "CtrlNeg" n s s p 1I verbosity
    if (s + result) % p = 0I then 
        if verbose then 
            show "Test passed: n s p = %A %A %A. Result: inverse sum %A %A\n" n s p result ((s + result) % p)
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s p = %A %A %A. Result: inverse sum %A %A\n" n s p result ((s + result) % p)
        0 // return integer code zero in case of failure

let RunCtrlModularAdderTest n s1 s2 p ctrl verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running controlled modular adder circuit test for %A + %A mod %A" s1 s2 p
        show "Control bit = %A" ctrl
        show "====================================================="
    let result =  RunModularAdder "CtrlModADD" n s1 s2 p ctrl verbosity
    if (ctrl*s1 + s2) % p = result then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result 
        0 // return integer code zero in case of failure

[<LQD>]
let __RunModularNegatorTests () = 
    RunModularNegatorTest 4 3I 11I "all" |> ignore
    RunModularNegatorTest 4 8I 11I "all" |> ignore
    RunModularNegatorTest 5 8I 29I "all" |> ignore
    RunModularNegatorTest 6 8I 41I "all" |> ignore
    
[<LQD>]
let __RunModularAdderTests () = 
    RunModularAdderTest 4 3I 5I 11I "all" |> ignore
    RunModularAdderTest 4 3I 8I 11I "all" |> ignore
    RunModularAdderTest 4 5I 7I 11I "all" |> ignore
    RunModularAdderTest 4 3I 5I 11I "verbose" |> ignore
    RunModularAdderTest 4 3I 8I 11I "verbose" |> ignore
    RunModularAdderTest 4 5I 7I 11I "verbose" |> ignore
    
    RunCtrlModularAdderTest 4 3I 5I 11I 0I "all" |> ignore
    RunCtrlModularAdderTest 4 3I 5I 11I 1I "all" |> ignore
    RunCtrlModularAdderTest 4 3I 8I 11I 0I "verbose" |> ignore
    RunCtrlModularAdderTest 4 3I 8I 11I 1I "verbose" |> ignore
    RunCtrlModularAdderTest 4 5I 7I 11I 0I "verbose" |> ignore
    RunCtrlModularAdderTest 4 5I 7I 11I 1I "verbose" |> ignore

[<LQD>]
let __RunBinaryDoubler (n:int) (s:int) (c:int) = 
    let k = Ket(n+2)
    let qs = k.Qubits
    s |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    c |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs n 1 |> ignore

    let BinaryDoublerCircuit = 
        Circuit.Compile (CtrlBinaryDoubling n) qs 
        //Circuit.Compile (CtrlBinaryHalfing n) qs 
        |> MyCircuitExport

    let currentState = Array.zeroCreate (n+2)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        currentState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast BinaryDoublerCircuit currentState
    show "The final state is   %A" finalState        
    show "done"
    
[<LQD>]
let __RunModularDoublerTest () = 
    for i in 1I..12I do 
        RunModularDoubler "ModDBL" 4 i 13I "all" |> ignore
        
let RunModularDoublerTest n s1 p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running modular doubler circuit test for 2*%A mod %A" s1 p
        show "====================================================="
    let result =  RunModularDoubler "ModDBL" n s1 p verbosity
    if (2I*s1) % p = result then 
        if verbose then 
            show "Test passed: n s1 p = %A %A %A. Result: double %A\n" n s1 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 p = %A %A %A. Result: double %A\n" n s1 p result 
        0 // return integer code zero in case of failure

[<LQD>]
let __RunModularDoublerTests () = 
    RunModularDoublerTest 4 3I 11I "all" |> ignore
    RunModularDoublerTest 4 8I 11I "all" |> ignore
    RunModularDoublerTest 3 2I 7I "all" |> ignore
    RunModularDoublerTest 3 6I 7I "all" |> ignore


let RunMultiplierTest (name:string) n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running %A circuit test for %A * %A mod %A" name s1 s2 p
        show "====================================================="
    let result =  RunMultiplier name n s1 s2 p verbosity
    if name = "Montgomery" then 
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
    else
        if (s1 * s2) % p = result then 
            if verbose then 
                show "Test passed: n s1 s2 p = %A %A %A %A. Results: product %A\n" n s1 s2 p result
            1 // return integer code one in case of success
        else 
            if verbose then 
                show "Test failed: n s1 s2 p = %A %A %A %A\n" n s1 s2 p 
            0 // return integer code zero in case of failure

//let RunMultiplierTest n s1 s2 p verbosity =  
//    let verbose = (verbosity = "verbose")
//    if verbose then 
//        show "Running Montgomery circuit test for %A * %A mod %A" s1 s2 p
//        show "====================================================="
//    let result =  RunMultiplier "Montgomery" n s1 s2 p verbosity
//    let decres =  RunMultiplier "Montgomery" n result 1I p "default"
//    let decode1 = RunMultiplier "Montgomery" n s1 1I p "default"
//    let decode2 = RunMultiplier "Montgomery" n s2 1I p "default"
//    if (decode1 * decode2) % p = decres then 
//        if verbose then 
//            show "Test passed: n s1 s2 p = %A %A %A %A. Results: product decprod decs1 decs2 %A %A %A %A\n" n s1 s2 p result decres decode1 decode2
//        1 // return integer code one in case of success
//    else 
//        if verbose then 
//            show "Test failed: n s1 s2 p = %A %A %A %A\n" n s1 s2 p 
//        0 // return integer code zero in case of failure


let RunSquarerTest (name:string) n s1 p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running %A circuit test for %A^2 mod %A" name s1 p
        show "====================================================="
    let result =  RunMultiplier name n s1 0I p verbosity
    if (s1 * s1) % p = result then 
        if verbose then 
            show "Test passed: n s1 p = %A %A %A. Results: product %A\n" n s1 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 p = %A %A %A\n" n s1 p 
        0 // return integer code zero in case of failure
        
//let RunSquarerTest n s p verbosity =  
//    let verbose = (verbosity = "verbose")
//    if verbose then 
//        show "Running Montgomery squarer circuit test for %A ^2 mod %A" s p
//        show "====================================================="
//    let result =  RunMultiplier "MontgomerySquarer" n s s p verbosity
//    let decres =  RunMultiplier "Montgomery" n result 1I p "default"
//    let decode = RunMultiplier "Montgomery" n s 1I p "default"
//    if (decode * decode) % p = decres then 
//        if verbose then 
//            show "Test passed: n s p = %A %A %A. Results: product decprod decs %A %A %A \n" n s p result decres decode
//        1 // return integer code one in case of success
//    else 
//        if verbose then 
//            show "Test failed: n s p = %A %A %A\n" n s p 
//        0 // return integer code zero in case of failure

[<LQD>]
let __RunSeqMultiplierTests () = 
    //RunInverter "Montgomery" 4 13I 8I "all" |> ignore
    let watch = Diagnostics.Stopwatch()
    for i in 3..120 do
        watch.Reset()
        watch.Start()
        let p = (pown 2I i)+1I
        RunMultiplierTest "Montgomery" i 2I 3I p "verbose" |> ignore
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
  
[<LQD>]
let __RunCheckAllInputsMultiplierTest (n:int) (p:int) = 
    let P = bigint p
    let L = [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunMultiplierTest "Montgomery" n i j P "default") ] ]
            // |> fun x -> printf "%A" x; x
            |> List.concat
    List.fold (*) 1 L 
    |> function 
        | 1 -> show "Test passed for all inputs: n p = %A %A" n p 
        | 0 -> let ind = List.findIndex (fun i -> (i=0)) L
               show "Test failed for some inputs: n p i j = %A %A %A %A" n p (ind/p) (ind%p) 
        | _ -> show "this should never happen"
        
[<LQD>]
let __RunSmallSquarerTests () = 
    RunSquarerTest "Montgomery" 4 12I 13I "verbose" |> ignore
    RunSquarerTest "Montgomery" 4 11I 13I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 4 7I 13I "verbose"   |> ignore
    RunSquarerTest "Montgomery" 5 13I 29I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 5 10I 17I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 6 7I 63I "verbose"   |> ignore

[<LQD>]
let __RunSmallMultiplierTests () = 
//    show "\n"
    RunMultiplierTest "Montgomery" 3 5I 2I 7I "all"   |> ignore
    RunMultiplierTest "Mersenne" 3 5I 2I 7I "all"   |> ignore
    RunMultiplierTest "DoubleAdd" 3 5I 2I 7I "all"   |> ignore
//    RunMultiplierTest "DoubleAddMersenne" 3 5I 2I 7I "all"   |> ignore
//    show "\n\n"
//    RunMultiplierTest "Montgomery" 4 5I 4I 15I "all"   |> ignore
//    RunMultiplierTest "Mersenne" 4 5I 4I 15I "all"   |> ignore
//    RunMultiplierTest "DoubleAdd" 4 5I 4I 15I "all"   |> ignore
//    RunMultiplierTest "DoubleAddMersenne" 4 5I 4I 15I "all"   |> ignore
//    show "\n\n"
//    RunMultiplierTest "Montgomery" 5 25I 19I 31I "all"   |> ignore
//    RunMultiplierTest "Mersenne" 5 25I 19I 31I "all"   |> ignore
//    RunMultiplierTest "DoubleAdd" 5 25I 19I 31I "all"   |> ignore
//    RunMultiplierTest "DoubleAddMersenne" 5 25I 19I 31I "all"   |> ignore
//    show "\n\n"
//    RunMultiplierTest "Montgomery" 6 53I 44I 63I "all"   |> ignore
//    RunMultiplierTest "Mersenne" 6 53I 44I 63I "all"   |> ignore
//    RunMultiplierTest "DoubleAdd" 6 53I 44I 63I "all"   |> ignore
//    RunMultiplierTest "DoubleAddMersenne" 6 53I 44I 63I "all"   |> ignore
//    show "\n\n"
//    RunMultiplierTest "Montgomery" 6 53I 44I 101I "all"   |> ignore
//    RunMultiplierTest "Mersenne" 6 53I 44I 101I "all"   |> ignore
//    RunMultiplierTest "DoubleAdd" 6 53I 44I 101I "all"   |> ignore
//    RunMultiplierTest "DoubleAddMersenne" 6 53I 44I 101I "all"   |> ignore
//    RunSquarerTest "DoubleAddSqu" 3 5I 7I "all"   |> ignore
//    RunSquarerTest "DoubleAddSqu" 3 2I 7I "all"   |> ignore
//    RunSquarerTest "DoubleAddSqu" 3 6I 7I "all"   |> ignore
//    RunSquarerTest "DoubleAddSqu" 3 1I 7I "all"   |> ignore
//    RunMultiplierTest "DoubleAdd" 7 88I 88I 127I "all"   |> ignore
//    RunSquarerTest "DoubleAddSqu" 7 88I 127I "all"   |> ignore
//    RunMultiplierTest "Mersenne" 7 88I 88I 127I "all"   |> ignore
//    RunSquarerTest "MersenneSqu" 7 88I 127I "all"   |> ignore

[<LQD>]
let __RunMediumMultiplierTests () = 
    RunMultiplierTest "Montgomery" 10 542I 7I 1021I "verbose"       |> ignore
    RunMultiplierTest "Montgomery" 20 500I 7I 1000003I "verbose"    |> ignore
    RunMultiplierTest "Montgomery" 32 500I 7I 2000000011I "verbose" |> ignore
    RunMultiplierTest "DoubleAdd" 32 500I 7I 2000000011I "verbose" |> ignore
    show "\n\n\n\n"
    RunMultiplierTest "Montgomery" 10 788I 534I 1023I "verbose"   |> ignore
    RunMultiplierTest "Mersenne" 10 788I 534I 1023I "verbose"   |> ignore

[<LQD>]
let __RunLargeMultiplierTests () = // following are all the primes from the NIST ECC standard (as of April 2010)
    let p192 = 6277101735386680763835789423207666416083908700390324961279I
    let p224 = 26959946667150639794667015087019630673557916260026308143510066298881I
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I
    let p384 = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319I
    let p521 = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151I 
    RunMultiplierTest "Montgomery" 192 500I 7I p192 "verbose" |> ignore
    RunMultiplierTest "Montgomery" 224 500I 7I p224 "verbose" |> ignore
    RunMultiplierTest "Montgomery" 256 500I 7I p256 "verbose" |> ignore
    RunMultiplierTest "Montgomery" 384 500I 7I p384 "verbose" |> ignore
    RunMultiplierTest "Montgomery" 521 500I 7I p521 "verbose" |> ignore
           
[<LQD>]
let __RunSweepMultiplierTests (startSize:int) (sweepLen:int) = 
    for i in startSize..(startSize+sweepLen) do 
        let p = (pown 2I i)-1I // note: in general this is not prime
        RunMultiplierTest "Montgomery" i 500I 7I p "verbose" |> ignore

[<LQD>]
let __RunCheckAllInputsModularArithmeticTest (name:string) (n:int) (p:int) = 
    let P = bigint p
    if (name = "Mersenne" || name = "MersenneSqu") && not (P = 2I**n - 1I) then failwith "Modulus is not the right power of 2."
    
    let L =
        match name with
        | "ModADD" ->     
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunModularAdderTest n i j P "default") ] ]
                |> List.concat
        | "CtrlModADD" ->
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunCtrlModularAdderTest n i j P 0I "default") ] ]
            @ [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunCtrlModularAdderTest n i j P 1I "default") ] ]
            |> List.concat
        | "ModDBL" ->     
            [ for i in 0I..(P-1I) do yield (RunModularDoublerTest n i P "default") ]
            //    |> List.concat
        | "Montgomery" -> 
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunMultiplierTest name n i j P "default") ] ]
                |> List.concat
        | "Mersenne" -> 
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunMultiplierTest name n i j P "default") ] ]
                |> List.concat
        | "MersenneSqu" ->     
            [ for i in 0I..(P-1I) do yield (RunSquarerTest name n i P "default") ]
        | "DoubleAdd" -> 
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunMultiplierTest name n i j P "default") ] ]
                |> List.concat
        | "DoubleAddSqu" ->     
            [ for i in 0I..(P-1I) do yield (RunSquarerTest name n i P "default") ]
        | "DoubleAddMersenne" -> 
            [ for i in 0I..(P-1I) do yield [ for j in 0I..(P-1I) do yield (RunMultiplierTest name n i j P "default") ] ]
                |> List.concat
        | _ -> failwith "This function is not implemented yet."

    List.fold (*) 1 L 
    |> function 
        | 1 -> show "Test passed for all inputs: n p = %A %A" n p 
        | 0 -> let ind = List.findIndex (fun i -> (i=0)) L
               if name = "ModDBL" then
                   show "Test failed for some inputs: n p i = %A %A %A" n p (ind%p) 
               else    
                   show "Test failed for some inputs: n p i j = %A %A %A %A" n p (ind/p) (ind%p) 
        | _ -> show "this should never happen"

let RunInverterTest n p s verbosity =  
    let verbose = (verbosity = "verbose")
    if verbose then 
        show "Running Montgomery inverter circuit test for %A mod %A" s p
        show "======================================================"
    let result =  RunInverter "Montgomery" n p s verbosity
    let decres =  RunMultiplier "Montgomery" n result 1I p "default"
    let decode = RunMultiplier "Montgomery" n s 1I p "default"
    if (decres * s) % p = 1I then 
        if verbose then 
            show "Test passed: n s p = %A %A %A. Results: inverse decinv decinput product %A %A %A %A\n" n s p result decres decode ((decres * decode) % p)
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s p = %A %A %A. Results: inverse decinv decinput product %A %A %A %A\n" n s p result decres decode ((decres * decode) % p)
        0 // return integer code zero in case of failure

let RunEfficientInverterTest n p s verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running efficient Montgomery inverter circuit test for %A mod %A" s p
        show "==============================================================="
    let resultSM =  RunEfficientInverter "MontgomerySM" n p s verbosity
    let resultMM =  RunEfficientInverter "MontgomeryMM" n p s verbosity
    let resultSS =  RunEfficientInverter "MontgomerySS" n p s verbosity   
    let decresSM =  RunMultiplier "Montgomery" n resultSM 1I p "default"
    let decresMM =  RunMultiplier "Montgomery" n resultMM 1I p "default"
    let decode = RunMultiplier "Montgomery" n s 1I p "default"
    if verbose then 
        match (decresSM * s) % p = 1I with 
        | true -> show "Test Standard -> Montgomery passed." 
        | _ -> show "Test Standard -> Montgomery failed." 
        show "n s p = %A %A %A. Results: inverse decinv input product %A %A %A %A\n" n s p resultSM decresSM s ((decresSM * s) % p)
        match (decresMM * decode) % p = 1I with 
        | true -> show "Test Montgomery -> Montgomery passed." 
        | _ -> show "Test Montgomery -> Montgomery failed." 
        show "n s p = %A %A %A. Results: inverse decinv decinput product %A %A %A %A\n" n s p resultMM decresMM decode ((decresMM * decode) % p)     
        match (resultSS * s) % p = 1I with 
        | true -> show "Test Standard -> Standard passed." 
        | _ -> show "Test Standard -> Standard failed." 
        show "n s p = %A %A %A. Results: inverse input product %A %A %A\n" n s p resultSS s ((resultSS * s) % p)     
    0 // return an integer code

[<LQD>]
let __RunSmallInverterTests () = 
    RunInverterTest 4 11I 10I "all" |> ignore
//    RunInverterTest 4 13I 8I "verbose" |> ignore
//    RunInverterTest 5 17I 10I "verbose" |> ignore
//    RunInverterTest 20 1000003I 500I "verbose" |> ignore    
//    RunInverter "Montgomery" 4 13I 8I "all" |> ignore
//    RunInverter "Montgomery" 5 17I 10I "all" |> ignore
//    RunInverter "Montgomery" 20 1000003I 500I "all"    |> ignore
////    for i in 1I..16I do 
//        RunInverter "Montgomery" 6 61I i "verbose"       |> ignore
//    let watch = Diagnostics.Stopwatch()
//    for i in 3..120 do
//        watch.Reset()
//        watch.Start()
//        let p = (pown 2I i)+1I
//        RunInverter "Montgomery" i p 2I "verbose" |> ignore
//        watch.Stop()
//        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
    
[<LQD>]
let __RunSmallEfficientInverterTests () = 
    RunEfficientInverterTest 4 11I 6I "all" |> ignore
//    RunEfficientInverterTest 4 13I 8I "verbose" |> ignore
//    RunEfficientInverterTest 5 17I 10I "verbose" |> ignore
//    RunEfficientInverterTest 20 1000003I 500I "verbose" |> ignore    

[<LQD>]
let __RunSmallEfficientInverterSweepTests () = 
//    RunInverter "Montgomery" 4 13I 8I "all" |> ignore
//    RunEfficientInverter "Montgomery" 4 13I 8I "all" |> ignore
//    RunInverter "Montgomery" 5 17I 10I "all" |> ignore
//    RunEfficientInverter "Montgomery" 5 17I 10I "all" |> ignore
//    RunInverter "Montgomery" 20 1000003I 500I "all"    |> ignore
//    RunEfficientInverter "Montgomery" 20 1000003I 500I "all"    |> ignore
//    for i in 1I..16I do 
//        RunEfficientInverter "Montgomery" 6 61I i "verbose"       |> ignore
    let watch = Diagnostics.Stopwatch()
    for i in 3..120 do
        watch.Reset()
        watch.Start()
        let p = (pown 2I i)-1I // note: in general not a prime; included for performance testing only
        RunEfficientInverterTest i p 3I "verbose" |> ignore
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
    
[<LQD>]
let __RunLargeEfficientInverterTests () = 
    let p192 = 6277101735386680763835789423207666416083908700390324961279I
    let p224 = 26959946667150639794667015087019630673557916260026308143510066298881I
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I
    let p384 = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319I
    let p521 = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151I 
    
    let watch = Diagnostics.Stopwatch()
    watch.Reset()
    watch.Start()
    RunEfficientInverter "Montgomery" 192 p192 2I "verbose" |> ignore
    watch.Stop()
    show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
    
//    RunInverter "Montgomery" 224 p224 2I "verbose" |> ignore
//    RunInverter "Montgomery" 256 p256 2I "verbose" |> ignore
//    RunInverter "Montgomery" 384 p384 2I "verbose" |> ignore
//    RunInverter "Montgomery" 521 p521 2I "verbose" |> ignore
//    

let RunECCPointAdditionTest name (n:int) P Q R (p:bigint) (b:bigint) verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running ECC point addition for affine Weierstrass curve y^2 = x^3-3x+%A over GF(%A)" b p 
        show "====================================================================================="
        show "Adding the points %A + %A = %A." P Q R
    // let result0 =  RunShorEllipticCurvePointAddition name n p P Q 0I verbosity
    let result1 =  RunShorEllipticCurvePointAddition name n p P Q 1I verbosity
    if verbose then 
        //match (result0 = P) && (result1 = R) with 
        match (result1 = R) with
        | true -> show "Test passed." 
        | _ -> show "Test failed." 
        show "n=%A p=%A P=%A Q=%A." n p P Q
        //show "result0=%A result1=%A R=%A.\n" result0 result1 R        
        show "result1=%A R=%A.\n" result1 R        
    0 // return an integer code
    
[<LQD>]
let __RunSmallECCPointAdditionTests () = 
    let p = 11I 
    let b = 5I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=5 over GF(11)
    let P = AffinePoint(0I,2I) // MG representation corresponding to (0 : 7 : 1) 
    let Q = AffinePoint(5I,8I) // MG representation corresponding to (1 : 6 : 1) 
    let R = AffinePoint(0I,9I) // MG representation corresponding to (0 : 4 : 1) 
    RunECCPointAdditionTest "AffineWeierstrass" 4 P Q R p b "verbose" |> ignore
    
    let p = 1000003I
    let b = 5I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=5 over GF(1000003)
    let P = AffinePoint(34236I,177368I) // MG representation corresponding to (563093 : 779719 : 1)
    let Q = AffinePoint(481427I,185548I) // MG representation corresponding to (159317 : 146176 : 1)
    let R = AffinePoint(39120I,600367I) // MG representation corresponding to (187760 : 407071 : 1)
    RunECCPointAdditionTest "AffineWeierstrass" 20 P Q R p b "verbose" |> ignore

[<LQD>]
let __RunMediumECCPointAdditionTests () = 
    let p = 1125899906842679I // 51 bit prime number
    let b = 111I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=111 over GF(112589990684267)
    let P = AffinePoint(451254016350526I,119715075055576I) // MG representation corresponding to (36839505554729 : 449271643872930 : 1)
    let Q = AffinePoint(102913185642772I,40875630210012I) // MG representation corresponding to (1002138888044907 : 490930180893078 : 1)
    let R = AffinePoint(939594529930895I,740730201703437I) // MG representation corresponding to (247344574010328 : 331036061128227 : 1)
    RunECCPointAdditionTest "AffineWeierstrass" 51 P Q R p b "verbose" |> ignore

    let p = 1267650600228229401496703205653I // 101 bit prime number
    let b = 500I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=500 over GF(1267650600228229401496703205653)
    let P = AffinePoint(1050973377032356303253151680466I,901012499748817281889417696786I) // MG representation corresponding to (332176921040233098041995516868 : 1211107952384860587565637691876 : 1)
    let Q = AffinePoint(1193382548062512834590307231000I,282766732732562641760008066570I) // MG representation corresponding to (784979104567596392202663529808 : 777470103510587425752922530425 : 1)
    let R = AffinePoint(52487571062535416645933108150I,549778904875717537058134673149I) // MG representation corresponding to (1025008991572534722840933398961 : 1218606481979730240903762949350 : 1)
    RunECCPointAdditionTest "AffineWeierstrass" 101 P Q R p b "verbose" |> ignore
    
[<LQD>]
let __RunLargeECCPointAdditionTests () = 
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I // 256 bit prime number
    let b = 41058363725152142129326129780047268409114441015993725554835256314039467401291I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b as above over GF(p256)
    // MG representation corresponding to (66790887510078570501074792547226042371050454650016343615879109835932523443065 : 106377730366318381478632396320403129066886107237681426316144016021673156344273 : 1)
    let P = AffinePoint(71528713545620684423536968445225241465408169212990955778170977604457805197600I,103275792837057906611974512458293333482266004537501318958510683525869453045362I) 
    // MG representation corresponding to (91776155225397256558779214231279668060857253502872526815652263470579408313795 : 103169371820847966113293214366645798727644724219171959581593614101799275769904 : 1)
    let Q = AffinePoint(55060438472976370102375399605978382023362429340763975384115074388233082886021I,31367148643791062976004566319732149799763757143358924179768669057101515390805I) 
    // MG representation corresponding to (455438157303235004258818524917930849403739804230436973770526759698198641623 : 107400709691385187225095324036499254562775224598883627815051692242731812890812 : 1)    
    let R = AffinePoint(76441926148553224235917927366185924393481358162275007948038451322418981861496I,83037667562094125842330603399066966639949874188752994804195456661992373332472I)
    RunECCPointAdditionTest "AffineWeierstrass" 256 P Q R p256 b "verbose" |> ignore

    let p521 = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151I // 521 bit prime number = (2^521-1)
    //let bhex = "051 953eb961 8e1c9a1f 929a21a0 b68540ee a2da725b 99b315f3 b8b48991 8ef109e1 56193951 ec7e937b 1652c0bd 3bb1bf07 3573df88 3d2c34f1 ef451fd4 6b503f00"
    let b = 1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b as above over GF(p521)
    let P = AffinePoint(2019402704654183728462351680707615222529965139422490780562105682799695337147007446363284789008754467299020228258354671619002960455124613248567210880508767296I,1458022328452401795217394740146225190614131805280093502448120620964052932115555877014656712347738560998519864959329334000186165661257935266396355500650605593I)
    // MG representation corresponding to 
    // AffinePoint(2019402704654183728462351680707615222529965139422490780562105682799695337147007446363284789008754467299020228258354671619002960455124613248567210880508767296 : 1458022328452401795217394740146225190614131805280093502448120620964052932115555877014656712347738560998519864959329334000186165661257935266396355500650605593 : 1) 
    let Q = AffinePoint(3461116442272996866124164996493831743691438070006772927809532245952862225431756771739852593264067723426836921091230666794367271350434171851828402338313755720I,1188953449582389783122107743290325019324372418987639141663599524115622781004635527902389228732004684685851076733855776288999119147282493981522107172719269490I)
    // MG representation corresponding to 
    // AffinePoint(3461116442272996866124164996493831743691438070006772927809532245952862225431756771739852593264067723426836921091230666794367271350434171851828402338313755720 : 1188953449582389783122107743290325019324372418987639141663599524115622781004635527902389228732004684685851076733855776288999119147282493981522107172719269490 : 1) 
    let R = AffinePoint(1520763616772063944347337203057795069608949511323431077489007758463057259475338376804177701424277256377356420277707785637897056096371105632036219405839402879I,5734643871458554093752783078290558101442887861644131534745261210924045539944388171451689122970667215696907182918233222773762151844628132090314629632726195757I)
    // MG representation corresponding to 
    // AffinePoint(1520763616772063944347337203057795069608949511323431077489007758463057259475338376804177701424277256377356420277707785637897056096371105632036219405839402879 : 5734643871458554093752783078290558101442887861644131534745261210924045539944388171451689122970667215696907182918233222773762151844628132090314629632726195757 : 1)
    RunECCPointAdditionTest "AffineWeierstrass" 521 P Q R p521 b "verbose" |> ignore

let RunFactoringTest name (n:int) (x:bigint) (m:bigint) (N:bigint) verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running factoring test %A * %A over Z/%A Z" x m N
        show "================================================="        
    let result1 =  RunShorFactoring name n N x m 1I verbosity

    if verbose then 
        //match (result0 = P) && (result1 = R) with 
        match name with 
        | "DoubleAdd" -> 
            match (result1 = ((x * m) % N)) with
            | true -> show "Test passed." 
            | _ -> show "Test failed." 
            show "n=%A x=%A m=%A N=%A result=%A." n x m N result1
        | "Montgomery" -> 
            match (((result1*(pown 2I n)) % N) = ((x * m) % N)) with 
            | true -> show "Test passed." 
            | _ -> show "Test failed." 
            show "n=%A x=%A m=%A N=%A result=%A (in MG encoding)." n x m N result1
    0 // return an integer code
    
[<LQD>]
let __RunSmallFactoringTests () = 
    let N = 3I * 5I 
    let x = 5I
    let m = 2I
    RunFactoringTest "Montgomery" 4 x m N "verbose" |> ignore
    
    let N = 17I * 19I
    let x = 29I
    let m = 200I
    RunFactoringTest "Montgomery" 9 x m N "verbose" |> ignore

    let N = 17I * 19I
    let x = 29I
    let m = 200I
    RunFactoringTest "DoubleAdd" 9 x m N "verbose" |> ignore

[<LQD>]
let __RunMediumFactoringTests () = 
    let N = 1048583I * 1048613I 
    let x = 5500005500I
    let m = 212121I
    RunFactoringTest "DoubleAdd" 41 x m N "verbose" |> ignore
    RunFactoringTest "Montgomery" 41 x m N "verbose" |> ignore

[<LQD>]
let __RunLargeFactoringTests () = 
    let n576  = 188198812920607963838697239461650439807163563379417382700763356422988859715234665485319060606504743045317388011303396716199692321205734031879550656996221305168759307650257059I
    let n1024 = 135066410865995223349603216278805969938881475605667027524485143851526510604859533833940287150571909441798207282164471551373680419703964191743046496589274256239341020864383202110372958725762358509643110564073501508187510676594629205563685529475213500852879416377328533906109750544334999811150056977236890927563I
    let n1536 = 1847699703211741474306835620200164403018549338663410171471785774910651696711161249859337684305435744585616061544571794052229717732524660960646946071249623720442022269756756687378427562389508764678440933285157496578843415088475528298186726451339863364931908084671990431874381283363502795470282653297802934916155811881049844908319545009848393775227257052578591944993870073695755688436933812779613089230392569695253261620823676490316036551371447913932347169566988069I
    let n2048 = 25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357I
    // primes of bit-sizes 576, 1024, 1536, and 2048 were taken from the RSA factoring challenge. 
  
    let x = n576 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    RunFactoringTest "DoubleAdd" 576 x m n576 "verbose" |> ignore
    RunFactoringTest "Montgomery" 576 x m n576 "verbose" |> ignore

    let x = n1024 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    RunFactoringTest "DoubleAdd" 1024 x m n1024 "verbose" |> ignore
    RunFactoringTest "Montgomery" 1024 x m n1024 "verbose" |> ignore

    let x = n1536 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    RunFactoringTest "DoubleAdd" 1536 x m n1536 "verbose" |> ignore
    RunFactoringTest "Montgomery" 1536 x m n1536 "verbose" |> ignore

    let x = n2048 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    RunFactoringTest "DoubleAdd" 2048 x m n2048 "verbose" |> ignore
    RunFactoringTest "Montgomery" 2048 x m n2048 "verbose" |> ignore
              
[<LQD>]
let __RunVeryLargeFactoringTests () = 
    let n4096 = 1044388881413152506691752710716624382579964249047383780384233483283953907971557456848826811934997558340890106714439262837987573438185793607263236087851365277945956976543709998340361590134383718314428070011855946226376318839397712745672334684344586617496807908705803704071284048740118609114467977783598029006686938976881787785946905630190260940599579453432823469303026696443059025015972399867714215541693835559885291486318237914434496734087811872639496475100189041349008417061675093668333850551032972088269550769983616369411933015213796825837188091833656751221318492846368125550225998300412344784862595674492194616443716246933213029777899798818461970932131191348804759825573468593610446932516843547031871082533988701912296616359089685540536703762087577894961935770084452075260024141226911707933733334829279200075988660067382055419272677487938772582789085572793356060472402630868147530739655123992685056875065481628365716808346058511893983996604674482947443330804864016159772117403184411098879310084911985820576707955476414019490642324082779864275147623753677862525374215759379112292614276897779676360046934875731128797682090481735253923263326681166800529934324565723067016719014885465039728219454090779142376015916232389820650894246737I
    // n4096 is the product p*q where p, q are 2048-bit primes given by
    // p = 32317006071311007300714876688669951960444102669715484032130345427524655138867890893197201411522913463688717960921898019494119559150490921095088152386448283120630877367300996091750197750389652106796057638384067568276792218642619756161838094338476170470581645852036305042887575891541065808607552399123930385521914333389668342420684974786564569494856176035326322058077805659331026192708460314150258592864177116725943603718461857357598351152301645904403697613233287231227125684710820209725157101726931323469678542580656697935045997268352998638215525166389437335543602135433229604645318478604952148193555853611059596222149I
    // q = 32317006071311007300714876688669951960444102669715484032130345427524655138867890893197201411522913463688717960921898019494119559150490921095088152386448283120630877367300996091750197750389652106796057638384067568276792218642619756161838094338476170470581645852036305042887575891541065808607552399123930385521914333389668342420684974786564569494856176035326322058077805659331026192708460314150258592864177116725943603718461857357598351152301645904403697613233287231227125684710820209725157101726931323469678542580656697935045997268352998638215525166389437335543602135433229604645318478604952148193555853611059596221213I

    let n8192 = 1090748135619415929462984244733782862448264161996232692431832786189721331849119295216264234525201987223957291796157025273109870820177184063610979765077554799078906298842192989538609825228048205159696851613591638196771886542609324560121290553901886301017900252535799917200010079600026535836800905297805880952350501630195475653911005312364560014847426035293551245843928918752768696279344088055617515694349945406677825140814900616105920256438504578013326493565836047242407382442812245131517757519164899226365743722432277368075027627883045206501792761700945699168497257879683851737049996900961120515655050115561271491492515342105748966629547032786321505730828430221664970324396138635251626409516168005427623435996308921691446181187406395310665404885739434832877428167407495370993511868756359970390117021823616749458620969857006263612082706715408157066575137281027022310927564910276759160520878304632411049364568754920967322982459184763427383790272448438018526977764941072715611580434690827459339991961414242741410599117426060556483763756314527611362658628383368621157993638020878537675545336789915694234433955666315070087213535470255670312004130725495834508357439653828936077080978550578912967907352780054935621561090795845172954115955492451891216570613296084365855602561987825554164828983776100070639535058975031421325993557963092582246894596377877165664134176626579977819203711571822341535408979557604056707123186626577676416188241336594881408068968961998717359787691335649357729429863437279636352583621433720199268128212734394658767583155547333863583775463140939229166786852053411359347506292273580759659221329615535915374296921312179266278677936047679491864303619458615535439684608832386143593871684717434990189386528537845712625304664823679920329431996145195217962866020168617725228731562929836932081530970550958194455848937712504788617035029734883083588942942311417818515459029999559814744757117521882790786148595159691476349244130338981981174851473059507740147241521672189074910330952319180638318199504090691810250168776309204618406888213015285193048095145396995645340002985115979602262321310662883088594437748318499511605309954705034226005722298162852349163242408574091397080892968895451060657757201199454907101458809408532429510189097387093136564036047797365062857395373088494023754833694236076111826663219393138920039495403085305056656519919814219844108778484235361412730008772810730108656701948930183929294294950489931276241738768497131921385470975219731917721I
    // n8192 is the product p*q where p, q are 4096-bit primes given by
    // p = 1044388881413152506691752710716624382579964249047383780384233483283953907971557456848826811934997558340890106714439262837987573438185793607263236087851365277945956976543709998340361590134383718314428070011855946226376318839397712745672334684344586617496807908705803704071284048740118609114467977783598029006686938976881787785946905630190260940599579453432823469303026696443059025015972399867714215541693835559885291486318237914434496734087811872639496475100189041349008417061675093668333850551032972088269550769983616369411933015213796825837188091833656751221318492846368125550225998300412344784862595674492194617023806505913245610825731835380087608622102834270197698202313169017678006675195485079921636419370285375124784014907159135459982790513399611551794271106831134090584272884279791554849782954323534517065223269061394905987693002122963395687782878948440616007412945674919823050571642377154816321380631045902916136926708342856440730447899971901781465763473223850267253059899795996090799469201774624817718449867455659250178329070473119433165550807568221846571746373296884912819520317457002440926616910874148385078411929804522981857338977648103126085903001302413467189726673216491511131602920781738033436090243804708340403154182269I
    // q = 1044388881413152506691752710716624382579964249047383780384233483283953907971557456848826811934997558340890106714439262837987573438185793607263236087851365277945956976543709998340361590134383718314428070011855946226376318839397712745672334684344586617496807908705803704071284048740118609114467977783598029006686938976881787785946905630190260940599579453432823469303026696443059025015972399867714215541693835559885291486318237914434496734087811872639496475100189041349008417061675093668333850551032972088269550769983616369411933015213796825837188091833656751221318492846368125550225998300412344784862595674492194617023806505913245610825731835380087608622102834270197698202313169017678006675195485079921636419370285375124784014907159135459982790513399611551794271106831134090584272884279791554849782954323534517065223269061394905987693002122963395687782878948440616007412945674919823050571642377154816321380631045902916136926708342856440730447899971901781465763473223850267253059899795996090799469201774624817718449867455659250178329070473119433165550807568221846571746373296884912819520317457002440926616910874148385078411929804522981857338977648103126085903001302413467189726673216491511131602920781738033436090243804708340403154181709I

    let n16384 = 1189731495357231765085759326628007130763444687096510237472674821233261358180483686904488595472612039915115437484839309258897667381308687426274524698341565006080871634366004897522143251619531446845952345709482135847036647464830984784714280967845614138476044338404886122905286855313236158695999885790106357018120815363320780964323712757164290613406875202417365323950267880089067517372270610835647545755780793431622213451903817859630690311343850657539360649645193283178291767658965405285113556134369793281725888015908414675289832538063419234888599898980623114025121674472051872439321323198402942705341366951274739014593816898288994445173400364617928377138074411345791848573595077170437644191743889644885377684738322240608239079061399475675334739784016491742621485229014847672335977897158397334226349734811441653077758250988926030894789604676153104257260141806823027588003441951455327701598071281589597169413965608439504983171255062282026626200048042149808200002060993433681237623857880627479727072877482838438705048034164633337013385405998040701908662387301605018188262573723766279240798931717708807901740265407930976419648877869604017517691938687988088008944251258826969688364194133945780157844364946052713655454906327187428531895100278695119323496808703630436193927592692344820812834297364478686862064169042458555136532055050508189891866846863799917647547291371573500701015197559097453040033031520683518216494195636696077748110598284901343611469214274121810495077979275556645164983850062051066517084647369464036640569339464837172183352956873912042640003611618789278195710052094562761306703551840330110645101995435167626688669627763820604342480357906415354212732946756073006907088870496125050068156659252761297664065498347492661798824062312210409274584565587264846417650160123175874034726261957289081466197651553830744424709698634753627770356227126145052549125229448040149114795681359875968512808575244271871455454084894986155020794806980939215658055319165641681105966454159951476908583129721503298816585142073061480888021769818338417129396878371459575846052583142928447249703698548125295775920936450022651427249949580708203966082847550921891152133321048011973883636577825533325988852156325439335021315312134081390451021255363707903495916963125924201167877190108935255914539488216897117943269373608639074472792751116715127106396425081353553137213552890539802602978645319795100976432939091924660228878912900654210118287298298707382159717184569540515403029173307292433820279730892035211089568317921472966412729818117888454622835234533481549572743196734571634129403963684506083367368909845783940774676061371524391659800174583099528253482739757522802678131797463892371605031905393218989649729585777477658297715183959931786937290802047985884915659668363557899564591609695996145630027327988546172359576120674838045740053640396189717133917288042702062394781105995775027475794051243596930032470884082214287907917283972886566217666099928929736659305895401188482643507350834137594417765254338122252480578715215887377442751706844188139483627322937448198918074963942373224651327246325315257389993715375779389735367127951944954922262846304494110166728898867798152486928838574921534377364723691418182036282103573740854602704670142191797776310269453779790001271690513762843505679799117285781202701919548686629332950416374405421584385147988271411879919984987621119528596275512991867785554210230568216704790576492500881265464921474895559459290185727309851736982251106571217022721811129816003904051940819258133219624777452655172398944653308802365249212729736452165680043648457346007468020497859675817059619455164828406637438597547809499476599841667405756636926594683316183551527318347414506716327282224104675912411296799047208452645924316981738234697893596112182781775518514959248984587435401557122123845271399302125402393337922749394455393324480804314744560744187982171228782269878349258607387176290025943618757718706454898247726013046950133059199119527299386145131246957734198972549900572364145426997855548761897904248475449664993993862428294060559878742860518378058716440614430083917881823389319994271944195984564613910857996369511365798097068122437154330466861277932794002897625255843879189558691957373003353663905696047447765667201581975169298500940792736093264986662256899062338573645314621899957140466596072241410164464110856818238166950911474536753400753329335262575148713324594079410166286488150567351475933291175562137155949959416056516881232905680654150013798790046699053524795019980299093858029924914578907051418730788232093062110578135615360678690695102511579902265905593484061796350040671437474498455768307744249390537909776724456222613234554195683081115022999893977896277499918028773373149365321042211436343920206922932005094285083739085355321384854749817569151104942327971701097722576153091798048279364363438424522407742110033167468367004558861384900599732150535646653656240721316539419706881676263475309286637885522950176418670198965931I
    // n16384 is the product p*q where p, q are 8192-bit primes given by
    // p = 1090748135619415929462984244733782862448264161996232692431832786189721331849119295216264234525201987223957291796157025273109870820177184063610979765077554799078906298842192989538609825228048205159696851613591638196771886542609324560121290553901886301017900252535799917200010079600026535836800905297805880952350501630195475653911005312364560014847426035293551245843928918752768696279344088055617515694349945406677825140814900616105920256438504578013326493565836047242407382442812245131517757519164899226365743722432277368075027627883045206501792761700945699168497257879683851737049996900961120515655050115561271491492515342105748966629547032786321505730828430221664970324396138635251626409516168005427623435996308921691446181187406395310665404885739434832877428167407495370993511868756359970390117021823616749458620969857006263612082706715408157066575137281027022310927564910276759160520878304632411049364568754920967322982459184763427383790272448438018526977764941072715611580434690827459339991961414242741410599117426060556483763756314527611362658628383368621157993638020878537675545336789915694234433955666315070087213535470255670312004130725495834508357439653828936077080978550578912967907352780054935621561090795845172954115972927479877527738560008204118558930004777748727761853813510493840581861598652211605960308356405941821189714037868726219481498727603653616298856174822413033485438785324024751419417183012281078209729303537372804574372095228703622776363945290869806258422355148507571039619387449629866808188769662815778153079393179093143648340761738581819563002994422790754955061288818308430079648693232179158765918035565216157115402992120276155607873107937477466841528362987708699450152031231862594203085693838944657061346236704234026821102958954951197087076546186622796294536451620756509351018906023773821539532776208676978589731966330308893304665169436185078350641568336944530051437491311298834367265238595404904273455928723949525227184617404367854754610474377019768025576605881038077270707717942221977090385438585844095492116099852538903974655703943973086090930596963360767529964938414598185705963754561497355827813623833288906309004288017321424808663962671333528009232758350873059614118723781422101460198615747386855096896089189180441339558524822867541113212638793675567650340362970031930023397828465318547238244232028015189689660418822976000815437610652254270163595650875433851147123214227266605403581781469090806576468950587661997186505665475715783551I
    // q = 1090748135619415929462984244733782862448264161996232692431832786189721331849119295216264234525201987223957291796157025273109870820177184063610979765077554799078906298842192989538609825228048205159696851613591638196771886542609324560121290553901886301017900252535799917200010079600026535836800905297805880952350501630195475653911005312364560014847426035293551245843928918752768696279344088055617515694349945406677825140814900616105920256438504578013326493565836047242407382442812245131517757519164899226365743722432277368075027627883045206501792761700945699168497257879683851737049996900961120515655050115561271491492515342105748966629547032786321505730828430221664970324396138635251626409516168005427623435996308921691446181187406395310665404885739434832877428167407495370993511868756359970390117021823616749458620969857006263612082706715408157066575137281027022310927564910276759160520878304632411049364568754920967322982459184763427383790272448438018526977764941072715611580434690827459339991961414242741410599117426060556483763756314527611362658628383368621157993638020878537675545336789915694234433955666315070087213535470255670312004130725495834508357439653828936077080978550578912967907352780054935621561090795845172954115972927479877527738560008204118558930004777748727761853813510493840581861598652211605960308356405941821189714037868726219481498727603653616298856174822413033485438785324024751419417183012281078209729303537372804574372095228703622776363945290869806258422355148507571039619387449629866808188769662815778153079393179093143648340761738581819563002994422790754955061288818308430079648693232179158765918035565216157115402992120276155607873107937477466841528362987708699450152031231862594203085693838944657061346236704234026821102958954951197087076546186622796294536451620756509351018906023773821539532776208676978589731966330308893304665169436185078350641568336944530051437491311298834367265238595404904273455928723949525227184617404367854754610474377019768025576605881038077270707717942221977090385438585844095492116099852538903974655703943973086090930596963360767529964938414598185705963754561497355827813623833288906309004288017321424808663962671333528009232758350873059614118723781422101460198615747386855096896089189180441339558524822867541113212638793675567650340362970031930023397828465318547238244232028015189689660418822976000815437610652254270163595650875433851147123214227266605403581781469090806576468950587661997186505665475715783381I
    
    let x = n4096 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    //RunFactoringTest "DoubleAdd" 4096 x m n4096 "verbose" |> ignore
    RunFactoringTest "Montgomery" 4096 x m n4096 "verbose" |> ignore

    let x = n8192 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    //RunFactoringTest "DoubleAdd" 8192 x m n8192 "verbose" |> ignore
    RunFactoringTest "Montgomery" 8192 x m n8192 "verbose" |> ignore

    let x = n16384 - (pown 2I 100)
    let m = (pown 2I 200) + 59I
    //RunFactoringTest "DoubleAdd" 16384 x m n16384 "verbose" |> ignore
    RunFactoringTest "Montgomery" 16384 x m n16384 "verbose" |> ignore
              
///////////////////////////////////////////////////////////////////////////
//
// Read ECC data file with tests
//
///////////////////////////////////////////////////////////////////////////

let pathECC = "data/ecc_50.dat"

// eccinstance encapsulates an instance of a point addition on an elliptic curve in affine 
// Weierstrass form y^2 = x^3-ax+b over GF(p). Initialization is from a format produced by a 
// Magma script. Running yields timing info, #qubits, and #gates which is all logged to a file.
type eccinstance public () = 
    let mutable p = 0I
    let mutable n = 0I 
    let mutable a = 0I
    let mutable b = 0I
    let mutable ps = AffinePoint(0I,0I) // input point P in standard encoding
    let mutable qs = AffinePoint(0I,0I) // input point Q in standard encoding
    let mutable rs = AffinePoint(0I,0I) // result point P+Q in standard encoding
    let mutable pm = AffinePoint(0I,0I) // input point P in Montgomery encoding
    let mutable qm = AffinePoint(0I,0I) // input point Q in Montgomery encoding
    let mutable rm = AffinePoint(0I,0I) // result point P+Q in Montgomery encoding
    member this.Prime = p
    member this.Bitsize = n
    
    member this.Init(s:string) =         
        let L = s.Split(' ') |> Array.map ((fun y -> (BigInteger.TryParse y)) >> snd) 
        p <- L.[0]; n <- L.[1]; a <- L.[2]; b <- L.[3]
        pm <- AffinePoint(L.[10],L.[11]); qm <- AffinePoint(L.[12],L.[13]); rm <- AffinePoint(L.[14],L.[15])
        
    member this.Run() =
        let watch = Diagnostics.Stopwatch()
        watch.Start()
        show "Running ECC point addition for affine Weierstrass curve y^2 = x^3-%Ax+%A over GF(%A)" a b p 
        show "==========================================================================================="
        show "Adding the points %A + %A = %A." pm qm rm
        let result1 =  RunShorEllipticCurvePointAddition "AffineWeierstrass" (int n) p pm qm 1I "none"
        match (result1 = rm) with
        | true -> show "Test passed." 
        | _ -> show "Test failed: p=%A n=%A a=%A n=%A p=%A q=%A r=%A" p n a b pm qm rm   
        //show "result0=%A result1=%A R=%A.\n" result0 result1 R        
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
    
[<LQD>]
let __RunMeasureECCPointAdditionBatch (path:string) =
    let lines = File.ReadAllLines path
    for data in lines do 
        let T = new eccinstance()
        T.Init(data)
        T.Run()
    show "done"


///////////////////////////////////////////////////////////////////////////
//
// Read RSA data file with tests
//
///////////////////////////////////////////////////////////////////////////

let pathRSA = "data/rsa_50.dat"

// RSAinstance encapsulates an instance of a modular multiplication modula an RSA number 
// Running yields timing info, #qubits, and #gates which is all logged to a file.
type RSAinstance public () = 
    let mutable n = 0I
    let mutable p = 0I
    let mutable q = 0I
    let mutable N = 0I
    let mutable x = 0I 
    let mutable m = 0I
    let mutable r = 0I
    let mutable xm = 0I
    let mutable mm = 0I
    let mutable rm = 0I
    member this.Modulus = (p*q)
    member this.Bitsize = n
    
    member this.Init(s:string) =         
        let L = s.Split(' ') |> Array.map ((fun y -> (BigInteger.TryParse y)) >> snd) 
        p <- L.[0]; q <- L.[1]; x <- L.[2]; m <- L.[3]; r <- L.[4]
        xm <- L.[5]; mm <- L.[6]; rm <- L.[7]
        N <- (p*q)
        
    member this.Run() =
        let watch = Diagnostics.Stopwatch()
        watch.Start()
        
        show "Running factoring test %A * %A over Z/%A Z" x m N
        show "================================================="        
        let result1 =  RunShorFactoring "DoubleAdd" (int n) N x m 1I "none"

        match (result1 = ((x * m) % N)) with
        | true -> show "Test passed." 
        | _ -> show "Test failed." 
        show "n=%A x=%A m=%A N=%A result=%A." n x m N result1      
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds

[<LQD>]
let __RunMeasureShorFactoringBatch (path:string) =
    let lines = File.ReadAllLines path
    for data in lines do 
        let T = new RSAinstance()
        T.Init(data)
        T.Run()
    show "done"


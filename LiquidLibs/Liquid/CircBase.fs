module Microsoft.Research.Liquid.CircBase

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open HamiltonianGates


///////////////////////////////////////////////////////////////////////////
//
// Some useful tools 
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Conversions and slicing 
///////////////////////////////////////////////////////////////////////////

let BoolInt n x = 
    [ for i in [0..(n-1)] do yield (int ((x >>> i) % 2I)) ] 

let IntBool n (x: int array) = 
    List.fold (+) 0I <| List.map (fun i -> (1I <<< i) * (bigint x.[i])) [0..x.Length-1]

let PrepBool (qs:Qubits) offset step ls = 
    List.mapi (fun i x -> if x <> 0 then X !!(qs,(offset+i*step))) ls 

let slice (qs:Qubits) ls = List.map (fun n-> qs.[n]) ls

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

// (Classical algorithm) determines a^-1 mod N 
// This is used for the modular multiplication to uncompute x after having computed y=ax mod N
let ModularInverse (a:bigint) (n:bigint) =
    let mutable t = bigint 0
    let mutable r = n
    let mutable newt = bigint 1
    let mutable newr = a

    while (newr = bigint 0) = false do
        let q = r / newr
        let oldt = newt
        newt <- (t - q*newt)
        t <- oldt
        let oldr = newr
        newr <- (r - q*newr)
        r <- oldr
    let t = (t+n)%n
    t

///////////////////////////////////////////////////////////////////////////
// Toffoli circuit simulation and pretty printing 
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
// Tools for stitching of Toffoli gate lists into larger circuits 
///////////////////////////////////////////////////////////////////////////

let MapToffoli (T:TofGate)  (embed:Map<int,int>) = 
    match T with 
    | MyNOT(a) ->  MyNOT(embed.[a]) 
    | MyCNOT(a,b) -> MyCNOT(embed.[a], embed.[b]) 
    | MyTOFF(a,b,c) -> MyTOFF(embed.[a], embed.[b], embed.[c]) 
  
let MapToffoliToLiquid (T:TofGate) (qs:Qubits) = 
    match T with 
    | MyNOT(a) ->  X [qs.[a]]
    | MyCNOT(a,b) -> CNOT [qs.[a]; qs.[b]]
    | MyTOFF(a,b,c) -> CCNOT [qs.[a]; qs.[b]; qs.[c]] 
   
let RemapToffoli (T : TofGate)  (embed:Dictionary<int,int>) = 
    match T with 
    | MyNOT(a) ->  MyNOT(embed.[a]) 
    | MyCNOT(a,b) -> MyCNOT(embed.[a], embed.[b]) 
    | MyTOFF(a,b,c) -> MyTOFF(embed.[a], embed.[b], embed.[c]) 
  
let MapToffoliCircuit (T: TofGate list) (embed:Map<int,int>) = 
    T |> List.map (fun x -> MapToffoli x embed)

let RemapToffoliCircuitOld (T: TofGate list) (embed:Dictionary<int,int>) = 
    T |> List.map (fun x -> RemapToffoli x embed)

let RemapToffoliCircuit (T: TofGate list) (embed:Dictionary<int,int>) = 
    T |> List.map (function
                       | MyNOT(a) ->  MyNOT(embed.[a]) 
                       | MyCNOT(a,b) -> MyCNOT(embed.[a], embed.[b]) 
                       | MyTOFF(a,b,c) -> MyTOFF(embed.[a], embed.[b], embed.[c]) 
                   )

let difference (L1:list<'a>) (L2:list<'a>) =
        let cache = HashSet<'a>(L2, HashIdentity.Structural)
        L1 |> List.filter (fun x -> not (cache.Contains x))

let rec concat (L:list<list<'a>>) =
    match L with 
    | [] -> []
    | x::[] -> x
    | x::xs -> List.append x (concat xs)

// Toffoli depth computed the depth of a Toffoli network by maximally parallelizing all gates in the network. 
// Clifford gates are discarded in the depth count. (Sasha 2/19/17)
let ToffoliDepthStreamed (toffCirc:TofGate list) = 
    let NumToffoliQubits = 
        toffCirc |> List.map (fun x -> match x with | MyTOFF(a,b,c) -> [a;b;c] | MyCNOT(a,b) -> [a;b] | MyNOT(a) -> [a] ) 
                 |> List.max |> List.max 
    let state = [| for i in 0..NumToffoliQubits -> 0 |]

    let processGate maxDepth gate = 
        match gate with 
        | MyTOFF(a,b,c) -> 
            let m1 = max state.[a] state.[b]
            let m2 = max m1 state.[c]
            let m = m2+1
            state.[a] <- m
            state.[b] <- m
            state.[c] <- m
            max maxDepth m
        | _ -> maxDepth

    List.fold (fun depth gate -> 
            processGate depth gate) 0 toffCirc


let CircuitDepthStreamed (circuit:(list<int>*int) list) = 
    let numQubits = circuit |> List.map (fun x -> (fst x)) |> List.max |> List.max |> (fun x -> x+1)
    let state = Array.zeroCreate numQubits
    
    let scheduleGate maxDepth gate = 
        let m = ((fst gate) |> List.map (fun x -> state.[x]) |> List.max) + (snd gate)
        for i in (fst gate) do 
            state.[i] <- m
        max maxDepth m
        
    List.fold (fun depth gate -> scheduleGate depth gate) 0 circuit


// Toffoli depth computed the depth of a Toffoli network by maximally parallelizing all gates in the network. 
// Clifford gates are discarded in the depth count. Patterned after circ.fold() function in Liquid. 
let ToffoliDepthSlow (toffCirc:TofGate list) = 
    let getIndices ts =
        List.map (fun T -> 
                    match T with 
                    | MyNOT(a) -> Set.ofSeq [a]
                    | MyCNOT(a,b) -> Set.ofSeq [a; b]
                    | MyTOFF(a,b,c) -> Set.ofSeq [a; b; c]
                 ) ts
    let isToffoli t = 
        (Set.count t = 3)
                    
    let rec parallelizeCircuit depth cs =
        match cs with
        | []                           -> depth
        | current::cs          ->
            // Find all other gates that don't have anything in common
            let rec parallelizeGate (cs:Set<int> list) obstructions ps ss =
                match cs with
                | []           -> ps,(List.rev ss)
                | indexSet::cs  ->
                    let nextObstructions  = Set.union indexSet obstructions
                    if (Set.intersect obstructions indexSet = Set.empty) && (Set.count current = 3) then
                        parallelizeGate cs nextObstructions (indexSet::ps) ss
                    else
                        parallelizeGate cs nextObstructions ps ((indexSet)::ss)
            let ps, ss   = parallelizeGate cs current [] []
            if ps.Length = 0 then // cannot parallelize anything
                match Set.count current with 
                | 3 -> parallelizeCircuit (depth+1) ss // Toffoli depth goes up by 1
                | _ -> parallelizeCircuit depth ss     // gate is Clifford, i.e., depth stays the same 
                
            else                
                let ps      = current :: ps                
                match Set.count current with 
                | 3 -> parallelizeCircuit (depth+1) ss // Toffoli depth goes up by 1
                | _ -> match List.tryFind isToffoli ps with // if Clifford then it depends on which gates are in parallel list ps
                       | Some _ -> parallelizeCircuit (depth+1) ss  // there is a Toffoli gate that got parallelized
                       | None   -> parallelizeCircuit (depth) ss    // all gates were Clifford, i.e., depth stays the same
                           
                
    toffCirc |> getIndices |> parallelizeCircuit 0


let PrintToffoliNetworkMetrics (toffCirc:TofGate list) = 
    let qubits = toffCirc |> List.map (fun x -> match x with | MyTOFF(a,b,c) -> List.max [a;b;c] | MyCNOT(a,b) -> List.max [a;b] | MyNOT(a) -> List.max [a] ) 
                          |> List.max |> (fun x -> x + 1)
    let depth = ToffoliDepthStreamed toffCirc
    let size = toffCirc |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
    show "Toffoli network metrics: qubits=%A, size=%A, depth=%A\n" qubits size depth

let ReportToffoliNetworkMetrics (toffCirc:TofGate list) =     
    let qubits = toffCirc |> List.map (fun x -> match x with | MyTOFF(a,b,c) -> List.max [a;b;c] | MyCNOT(a,b) -> List.max [a;b] | MyNOT(a) -> List.max [a] ) 
                          |> List.max |> (fun x -> x + 1)
    let depth = ToffoliDepthStreamed toffCirc
    let size = toffCirc |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
    qubits, size, depth

///////////////////////////////////////////////////////////////////////////
//
// Quantum circuit gadgets
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Decomposing multiply controlled NOT gates with dirty ancillas
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

///////////////////////////////////////////////////////////////////////////
// Quantum counters based on binary shift registers
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
    

///////////////////////////////////////////////////////////////////////////
// Doubling and halfing
///////////////////////////////////////////////////////////////////////////

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


///////////////////////////////////////////////////////////////////////////
// SWAP and controlled SWAP
///////////////////////////////////////////////////////////////////////////

//let SWAP (xs:Qubits) (ys:Qubits) =
//    let n = xs.Length
//    // perform swap between xs and ys
//    for i in 0..(n-1) do
//        CNOT [xs.[i]; ys.[i]]
//        CNOT [ys.[i]; xs.[i]] 
//        CNOT [xs.[i]; ys.[i]]
//        
//let SWAPCompile (qs:Qubits) =
//    let n = int (qs.Length/2)
//    let xs = slice qs [0..(n-1)]
//    let ys = slice qs [n..(2*n-1)]    
//    SWAP xs ys 

let CSWAP (xs:Qubits) (ys:Qubits) (ctrls:Qubits) =
    let n = xs.Length
    // perform swap between x and y conditioned on ctrl
    for i in 0..(n-1) do
        BuildMultiplyControlledNOT [xs.[i]] [ys.[i]] [xs.[(i+1)%n]]
        BuildMultiplyControlledNOT (List.concat [[ys.[i]];ctrls]) [xs.[i]] [ys.[(i+1)%n]]
        BuildMultiplyControlledNOT [xs.[i]] [ys.[i]] [xs.[(i+1)%n]]

let CSWAPCompile (numctrls:int) (qs:Qubits) =
    let n = int (qs.Length-numctrls)/2
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let cs = slice qs [(2*n)..(2*n+numctrls-1)]
    CSWAP xs ys cs


///////////////////////////////////////////////////////////////////////////
// Efficient quantum state prep
///////////////////////////////////////////////////////////////////////////

let get_probabilities (ps:list<float>) (bit:int) =
    let num_probabilities = (1 <<< bit)
    let mask = num_probabilities - 1

    let get_probability p_index =
        let norm_list = List.map2 (fun i p -> match (i &&& mask) = p_index with
                                                | true -> p
                                                | false -> 0.) [0..ps.Length-1] ps
        let norm = List.reduce (+) norm_list
        let norm = match norm with
                    | 0. -> 1.
                    | _ -> norm
        let p_list = List.map2 (fun i p -> match (i >>> bit) &&& 1 = 0 with
                                            | true -> p / norm
                                            | false -> 0.) [0..ps.Length-1] norm_list
        List.reduce (+) p_list
    List.map get_probability [0..num_probabilities - 1]

let rec build_uniform_gate (angles:list<float>) (qs:Qubits) =
    if qs.Length <= 1 then
        if qs.Length = 1 then
            Rz (2. * angles.[0]) 1. "" [qs.[0]]
    else
        let even_angles = List.map (fun i -> angles.[2*i]) [0..angles.Length/2-1]
        let odd_angles = List.map (fun i -> angles.[2*i+1]) [0..angles.Length/2-1]
        let new_angles_plus = List.map2 (fun a1 a2 -> (a1 + a2) * 0.5) even_angles odd_angles
        build_uniform_gate new_angles_plus qs.[1..]
        CNOT [qs.[0]; qs.[qs.Length - 1]]
        let new_angles_minus = List.map2 (fun a1 a2 -> (a1 - a2) * 0.5) even_angles odd_angles
        build_uniform_gate new_angles_minus qs.[1..]
        CNOT [qs.[0]; qs.[qs.Length - 1]]

// prepares sum_i \sqrt(p_i) |i>
let StatePrep (p:list<float>) (qs:Qubits) =
    for i in 0..qs.Length-1 do
        let ps = get_probabilities p i
        let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
        Ybasis [qs.[i]]
        build_uniform_gate angles qs.[0..i]
        Ybasis' [qs.[i]]

// unprepares sum_i \sqrt(p_i) |i>
let StatePrep' (p:list<float>) (qs:Qubits) =
    for i in qs.Length-1..(-1)..0 do
        let ps = get_probabilities p i
        let angles = ps |> List.map (fun pr -> -Math.Acos(Math.Sqrt(pr)))
        Ybasis [qs.[i]]
        build_uniform_gate angles qs.[0..i]
        Ybasis' [qs.[i]]

[<LQD>]
let __TestLCUAmpExtraQubit () =
    let n = 3 // # bits
    let k = Ket(n)
    let qs = k.Qubits
    let p = [0.5; 0.5]

    let W (qs:Qubits) =
        Ry (2. * Math.Acos (Math.Sqrt 0.5)) 1. "" [qs.[1]]
        StatePrep p qs.[0..0]
        CNOT [qs.[0]; qs.[2]]
        X [qs.[0]]
        CZ [qs.[0]; qs.[2]]
        X [qs.[0]]
        StatePrep' p qs.[0..0]
    
    let W' (qs:Qubits) =
        Ry (-2. * Math.Acos (Math.Sqrt 0.5)) 1. "" [qs.[1]]
        StatePrep p qs.[0..0]
        CNOT [qs.[0]; qs.[2]]
        X [qs.[0]]
        CZ [qs.[0]; qs.[2]]
        X [qs.[0]]
        StatePrep' p qs.[0..0]

    let R (qs:Qubits) =
        X >< [qs.[0];qs.[1]]
        CZ [qs.[0];qs.[1]]
        X >< [qs.[0];qs.[1]]

    let ampl (qs:Qubits) =
        W qs
        for i in 0..0 do
            Gtheta (-Math.PI) 1. qs
            R qs
            W' qs
            R qs
            W qs
    let amplify  = Circuit.Compile ampl qs
    let mop (qs:Qubits) =
        M [qs.[0]]
        M [qs.[1]]

    let meas = Circuit.Compile mop qs

    amplify.Run qs
    show "Probability of success (00): %A" (k.Probs [qs.[0];qs.[1]]).[0]
    
    meas.Run qs
    let state = k.State qs.[qs.Length-1]
    let vec = (k.Probs [qs.[qs.Length-1]])
    show "State: %A Probabilities: %A" state vec
    let outcome = match [qs.[0].Bit;qs.[1].Bit] with
                    | [Bit.Zero;Bit.Zero] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!!"
    show "%A" outcome

[<LQD>]
let __TestLCUAmpPlusMinusID () =
    let n = 3 // # bits
    let k = Ket(n)
    let qs = k.Qubits
    let p = [0.5/Math.Sqrt 2.; 0.5/Math.Sqrt 2.; 0.5-0.5/Math.Sqrt 2.; 0.5-0.5/Math.Sqrt 2.]

    let W (qs:Qubits) =
        StatePrep p qs.[0..1]
        
        X [qs.[1]]
        CCNOT [qs.[0]; qs.[1]; qs.[2]]
        X [qs.[0]]
        Cgate CZ [qs.[0]; qs.[1]; qs.[2]]
        X [qs.[0]]
        X [qs.[1]]
        CZ [qs.[0]; qs.[1]] // distribute the remaining amplitudes equally among +/- ID contributions

        StatePrep' p qs.[0..1]

    let R (qs:Qubits) =
        X >< [qs.[0];qs.[1]]
        CZ [qs.[0];qs.[1]]
        X >< [qs.[0];qs.[1]]

    let ampl (qs:Qubits) =
        W qs
        for i in 0..0 do
            Gtheta (-Math.PI) 1. qs
            R qs
            W qs
            R qs
            W qs
        
    let amplify  = Circuit.Compile ampl qs
    let mop (qs:Qubits) =
        M [qs.[0]]
        M [qs.[1]]
    let meas = Circuit.Compile mop qs

    amplify.Run qs
    show "Probability of success (00): %A" (k.Probs [qs.[0];qs.[1]]).[0]
    
    meas.Run qs
    let state = k.State qs.[qs.Length-1]
    let vec = (k.Probs [qs.[qs.Length-1]])
    show "State: %A Probabilities: %A" state vec
    let outcome = match [qs.[0].Bit;qs.[1].Bit] with
                    | [Bit.Zero;Bit.Zero] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!!"
    show "%A" outcome

[<LQD>]
let __TestLCU () =
    let n = 6 // # bits
    let p = [0.25;0.1;0.1;0.05;0.1;0.025;0.125;0.25]
    let p = [0.24;0.20;0.3;0.2;0.05;0.;0.01;0.]
    let k = Ket(n)
    let qs = k.Qubits
    
    let W (qs:Qubits) =
        for i in 0..2 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i+3]]
            build_uniform_gate angles qs.[3..i+3]
            Ybasis' [qs.[i+3]]
        for i in 0..2 do
            CNOT [qs.[3+i]; qs.[i]]
        for i in 2..(-1)..0 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> -Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i+3]]
            build_uniform_gate angles qs.[3..i+3]
            Ybasis' [qs.[i+3]]

    let ops (qs:Qubits) =
        W qs
        

    let circ  = Circuit.Compile ops qs
    let mop (qs:Qubits) =
        M >< qs.[3..]
    let meas = Circuit.Compile mop qs
    circ.Run qs
    let norm2 = List.fold (fun s p -> s + p * p) 0. p
        
    show "Probability of success (000): %A == %A" (k.Probs qs.[3..]).[0] norm2
    meas.Run qs
    let vec = (k.Probs qs.[0..2])
    let one_norm = Array.fold (fun s p -> s + Math.Sqrt p) 0. vec
    show "%A" (List.map (fun i -> (Math.Sqrt vec.[i])/one_norm) [0..7])
    let measurement = List.map (fun i -> (k.Bit qs.[i+3]).v) [0..2]
    let status = match measurement with
                    | [0;0;0] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!! Please, try again."
    show "%A : %A" measurement status
    show "Done"

[<LQD>]
let __TestStatePrep () =
    let n = 4 // # bits
    //let p = [0.25;0.1;0.1;0.05;0.1;0.025;0.125;0.25]
    let p = [0.016393443;0.06557377;0.049180328;0.147540984;0.016393443;0.049180328;0.032786885;0.016393443;0.081967213;0.131147541;0.06557377;0.049180328;0.147540984;0.032786885;0.081967213;0.016393443]
    let k = Ket(n)
    let qs = k.Qubits
    
    let ops (qs:Qubits) =
        for i in 0..n-1 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i]]
            build_uniform_gate angles qs.[0..i]
            Ybasis' [qs.[i]]

    let circ  = Circuit.Compile ops qs
    circ.Run qs
    let vec = (k.Probs qs)
    for i in 0..vec.Length-1 do
        if Math.Abs (vec.[i] - p.[i]) > 1.e-9 then
            show "Failed!!!"
    show "%A" (vec)
    show "Done"
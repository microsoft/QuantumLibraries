// Quantum circuits for Montgomery arithmetic and elliptic curve point additions
// =============================================================================
// M. Roetteler, M. Naehrig, MSR, May 2016: Montgomery arithmetic, Shor for ECC
// M. Roetteler, Th. Haener, MSR, August 2016: Dirty qubits, Shor for factoring

// Todos: 
// - implement in place modular adder: DONE
// - implement out of place integer adder: DONE
// - implement out of place modular adder: DONE
// - introduce type Adder and type Multiplier 
// - introduce type Bitvector: DONE
// - add static methods for adder, using optional arguments: Cuccaro, Ripple, Takahasi
// - add further optional arguments: "mod 2^n", "controlled (by a list of qubits)", "inverse", "inplace"
// - implement the following transformation: take arbitrary circuit as input, then add 1 control

// - refactor inverse and multiplier to make it input, output, ancilla format
// - rewrite simulator to adapt it to binary ops: DONE
// - write function for circuit "fold" to compute the depth: DONE
// - implement mapper for general multipliers: DONE
// - implement entire affine point addition: DONE
// - fix problem with gate counts for Mersenne: DONE

// - reimplement MultiplyControlledNOT for cases where dirty ancillas are already present

module Microsoft.Research.Liquid.Montgomery

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions


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
    let qubits = toffCirc |> List.map (fun x -> match x with | MyTOFF(a,b,c) -> [a;b;c] | MyCNOT(a,b) -> [a;b] | MyNOT(a) -> [a] ) 
                          |> List.max |> List.max |> (fun x -> x + 1)
    let depth = ToffoliDepthStreamed toffCirc
    let size = toffCirc |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
    show "Toffoli network metrics: qubits=%A, size=%A, depth=%A\n" qubits size depth


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
//
// Quantum integer arithmetic
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Integer addition
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

// Case-by-case definitions of various integer adder types  

// Cuccaro adder (described e.g. in [Cuccaro et al, quant-ph/0410184])
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

// Have to implement this!!!
//XXX 
//let BuildTakahashiAdderInverseConstant (ms:bigint) (xs:Qubits) = 
//    let n = xs.Length 
//    
//    for i in 0..(n-1) do 
//        CNOT [ys.[i]; xs.[i]]
//    for i in (n-2)..(-1)..1 do 
//        CNOT [ys.[i]; ys.[i+1]]
//    for i in 1..(n-1) do 
//        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
//        CNOT [ys.[i]; xs.[i]]
//    for i in (n-1)..(-1)..0 do 
//        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
//    for i in 1..(n-1) do 
//        CNOT [ys.[i]; ys.[i+1]]
//    for i in 1..(n-1) do 
//        CNOT [ys.[i]; xs.[i]]
        
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

// Works for small size n=2 as well but then makes use of an ancilla qubit a to 
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

// Works for small sizes n=2 as well but then makes use of an ancilla qubit a
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

let BuildStrictlyLargerThanComparator (xs:Qubits) (ys:Qubits) = 
    let n = xs.Length 
    X >< xs // using the trick that x-y = (x'+y)', where ' denotes the one's complement
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
        // CCNOT [ys.[0]; ctrl.[0]; xs.[0]] taken out
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        // CCNOT [ys.[i]; ctrl.[0]; xs.[i]] taken out 
    MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CNOT [ys.[n-1]; ys.[n]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    X >< xs

let BuildCtrlStrictlyLargerThanComparator (xs:Qubits) (ys:Qubits) (ctrl:Qubits) = 
    let n = xs.Length 
    X >< xs // using the trick that x-y = (x'+y)', where ' denotes the one's complement
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
        // CCNOT [ys.[0]; ctrl.[0]; xs.[0]] taken out
    for i in (n-2)..(-1)..1 do 
        CNOT [ys.[i]; ys.[i+1]]
    for i in 1..(n-1) do 
        CCNOT [ys.[i-1]; xs.[i-1]; ys.[i]]
        // CCNOT [ys.[i]; ctrl.[0]; xs.[i]] taken out 
    MultiplyControlledNOT [ys.[n-1]; xs.[n-1]; ctrl.[0]; ys.[n]; ys.[n-2]] // needs one dirty ancilla taken to be ys.[n-2]
    for i in (n-2)..(-1)..0 do 
        CCNOT [ys.[i]; xs.[i]; ys.[i+1]]
    for i in 1..(n-2) do 
        CNOT [ys.[i]; ys.[i+1]]
    CCNOT [ys.[n-1]; ctrl.[0]; ys.[n]]
    for i in 1..(n-1) do 
        CNOT [ys.[i]; xs.[i]]
    X >< xs

// wrapper functions to unpack from a qubit register as required for the addition         

// CuccaroAdder is a wrapper function for the Cuccaro in-place addition circuit
let CuccaroAdder (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ys = slice qs [n..2*n-1]
    let z  = slice qs [2*n]
    let a  = slice qs [2*n+1]
    BuildCuccaroAdder xs ys z a

// RippleAdder is a wrapper function for the carry ripple out-of-place adder
let RippleAdder (qs:Qubits) = 
    let n = (qs.Length-1)/3
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let z  = slice qs [2*n]
    let cs = slice qs [2*n+1..3*n]
    BuildRippleAdder xs ys z cs

// TakahashiAdder is a wrapper function for the Takahashi in-place addition circuit
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
    
let StrictlyLargerThanComparator (qs:Qubits) = 
    let n = (qs.Length-1)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    BuildStrictlyLargerThanComparator xs ys

let CtrlStrictlyLargerThanComparator (qs:Qubits) = 
    let n = (qs.Length-2)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let ctrl = slice qs [(2*n+1)]
    BuildCtrlStrictlyLargerThanComparator xs ys ctrl

    
///////////////////////////////////////////////////////////////////////////
// Functions to run integer adders
///////////////////////////////////////////////////////////////////////////

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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constant incrementer with dirty qubits
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  * Addition of constants in 10nlog(n) using 1 and 8nlog(n) using 2 dirty ancilla qubits
//  * Incrementer (by Gidney) using n dirty ancillae in ~4n Toffolis
//  * Carry computation using n dirty ancillae in ~4n Toffolis
//  * Addition of a constant modulo N in ~twice the cost of addition of constant
//  * Modular multiplication using repeated modular addition and shift including uncomputation, i.e.
//          |x>|0> --> |(ax)mod N>|0>
//      in 32n^2 log(n) operations

#if FALSE
// Uses qubits xg of the form [x_0,...,x_n-1,garbage,...,garbage] starting at the LSB
// to add the constant c to x.
// This is an n^2 algorithm, see ConstantAdderInplace for an nlog(n) version
// using only 1 or 2 dirty qubits.
let ConstantAdderUsingGarbage (c:int) (xg:Qubits) = 
    let n = xg.Length/2
    let x = slice xg [0..(n-1)]
    let g = slice xg [n..(2*n-1)]
    for r in (n-1)..(-1)..0 do
        if r > 0 then
            CNOT [g.[r-1]; x.[r]]
        for i in (r-1)..(-1)..0 do
            if ((c >>> i)&&&1 = 1) then
                CNOT [x.[i]; g.[i]]
                if i > 0 then
                    CNOT [g.[i-1]; g.[i]]
            if i > 0 then
                CCNOT [g.[i-1]; x.[i]; g.[i]]
        for i in 0..(r-2) do
            CCNOT [g.[i]; x.[i+1]; g.[i+1]]
            if ((c >>> (i+1))&&&1 = 1) then
                CNOT [g.[i]; g.[i+1]]
        if r > 0 then
            CNOT [g.[r-1]; x.[r]]
        if (c >>> r)&&&1 = 1 then
            X [x.[r]]

        for i in (r-1)..(-1)..0 do
            if ((c >>> i)&&&1 = 1) then
                CNOT [x.[i]; g.[i]]
                if i > 0 then
                    CNOT [g.[i-1]; g.[i]]
            if i > 0 then
                CCNOT [g.[i-1]; x.[i]; g.[i]]
        for i in 0..(r-2) do
            CCNOT [g.[i]; x.[i+1]; g.[i+1]]
            if ((c >>> (i+1))&&&1 = 1) then
                CNOT [g.[i]; g.[i+1]]
#endif

// n-bit incrementer using n dirty bits
let DirtyQubitsIncrementer (xg:Qubits) = 
    let n = int (xg.Length/2)
    if n > 3 then
        let g = slice xg [n..(2*n-1)]
        let x = slice xg [0..(n-1)]
        let gx = List.concat [g;x]
        TakahashiModAdderInverse gx
        
        for i in 0..(n-1) do
            X [g.[i]]
        TakahashiModAdderInverse gx
        for i in 0..(n-1) do
            X [g.[i]]

    elif n = 2 then
        CNOT [xg.[0]; xg.[1]]
        X [xg.[0]]
    elif n = 3 then
        CCNOT [xg.[0]; xg.[1]; xg.[2]]
        CNOT [xg.[0]; xg.[1]]
        X [xg.[0]]
    elif n = 1 then
        X [xg.[0]]

// Controlled n-bit incrementer using n dirty qubits
// Last qubit is the control
let CtrlDirtyQubitsIncrementer (xgc:Qubits) = 
    let n = int ((xgc.Length-1)/2)
    let xg = slice xgc [0..(2*n-1)]
    let c = xgc.[2*n]
    if n > 2 then
        let g = slice xg [n..(2*n-1)]
        let x = slice xg [0..(n-1)]
        let gxc = List.concat [g;x;[c]]
        CtrlTakahashiModAdderInverse gxc
        
        for i in 0..(n-1) do
            X [g.[i]]
        CtrlTakahashiModAdderInverse gxc
        for i in 0..(n-1) do
            X [g.[i]]

    elif n = 2 then
        CCNOT [c; xg.[0]; xg.[1]]
        CNOT [c; xg.[0]]
    elif n = 1 then
        CNOT [c;xg.[0]]

// computes the carry of the computation x += c using dirty qubits g
// at round r (i.e. determines if the r-th bit of a would be flipped due to a carry
// propagation from less-significant bits)
let ComputeCarryUsingGarbage (c:bigint) (g:Qubits) (x:Qubits) (ctrls:Qubits) = 
    let r = x.Length - 1
    let one = bigint 1
    if r = 1 then
        if (c&&&one) = one then
            BuildMultiplyControlledNOT (List.concat [[x.[0]];ctrls]) [x.[1]] [g.[0]]
    else
        // make this controlled on a control-qubit register to have a conditional addition
        BuildMultiplyControlledNOT (List.concat [[g.[r-2]];ctrls]) [x.[r]] [x.[r-1]]
        for i in (r-1)..(-1)..2 do
            if ((c >>> i)&&&one = one) then
                CNOT [x.[i]; g.[i-1]]
                X [x.[i]]
            CCNOT [g.[i-2]; x.[i]; g.[i-1]]

        if ((c>>>1)&&&one) = one then
            CNOT [x.[1];g.[0]]
        if (c&&&one) = one then
            if ((c>>>1)&&&one) = one then
                X [x.[1]]
            CCNOT [x.[0];x.[1];g.[0]]

        for i in 1..(r-2) do
            CCNOT [g.[i-1]; x.[i+1]; g.[i]]

        BuildMultiplyControlledNOT (List.concat [[g.[r-2]];ctrls]) [x.[r]] [x.[r-1]]
        
        for i in (r-2)..(-1)..1 do
            CCNOT [g.[i-1]; x.[i+1]; g.[i]]

        if (c&&&one) = one then
            CCNOT [x.[0];x.[1];g.[0]]
            if ((c>>>1)&&&one) = one then
                X [x.[1]]
        if ((c>>>1)&&&one) = one then
            CNOT [x.[1];g.[0]]
            
        for i in 2..(r-1) do
            CCNOT [g.[i-2]; x.[i]; g.[i-1]]
            if ((c >>> i)&&&one = one) then
                X [x.[i]]
                CNOT [x.[i]; g.[i-1]]
        

// Adds the constant c to the x-register using 1 garbage qubit
// ctrls: Control qubits (if needed)

// Runs in 10 nlog(n) Toffoli gates for 1 garbage qubit
// and in 8 nlog(n) with 2 garbage qubits (chooses automatically)
let ConstantAdderInplace (c:bigint) (x:Qubits) (garbage:Qubits) (ctrls:Qubits) =
    let n = x.Length
    let one = bigint 1
    let g = garbage.[0]

    // Recursively splits the x-register into 2 parts and applies the incrementer/addition
    // decomposition. The dirty qubit (g) is used do hold the carry. Conditioned 
    // on g, the higherbits register is incremented. If two garbage qubits are available,
    // the increment operation is carried out on [[g];higherbits], followed by X [g]
    // which is equivalent to a controlled incrementer but costs less than 2 ctrltakahashi
    // adders.
    let rec CircGen (c:bigint) (x:Qubits) (frombit:int) = 
        let n = x.Length           
        // Handle cases with more than 2 qubits by recursive splitting & incrementing
        // If there were a carry when adding c to x when going from bit L to L+1.
        if n >= 2 then
            let L = x.Length-(int x.Length/2)
            let lowerbits = slice x [0..(L-1)]
            let higherbits = slice x [L..(x.Length-1)]
            
            let incrinput_lowbits = slice lowerbits [0..((higherbits.Length)-1)]
            let incrinput = List.concat [higherbits; incrinput_lowbits; [g]]
            // (conditionally) increment
            if garbage.Length >=2 then
                DirtyQubitsIncrementer (List.concat [[g];higherbits;incrinput_lowbits;[garbage.[1]]])
                X [g]
            else
                CtrlDirtyQubitsIncrementer incrinput

            // conditionally invert
            for i in 0..(higherbits.Length-1) do
                CNOT [g; higherbits.[i]]

            // compute carry
            let carryinputx = List.concat [lowerbits; [g]]

            ComputeCarryUsingGarbage (c >>> frombit) higherbits carryinputx ctrls
            
            // (conditionally) in/de-crement
            if garbage.Length >=2 then
                DirtyQubitsIncrementer (List.concat [[g];higherbits;incrinput_lowbits;[garbage.[1]]])
                X [g]
            else
                CtrlDirtyQubitsIncrementer incrinput

            // uncompute carry (i.e. repeat circuit from before)
            ComputeCarryUsingGarbage (c >>> frombit) higherbits carryinputx ctrls

            // conditionally invert            
            for i in 0..(higherbits.Length-1) do
                CNOT [g; higherbits.[i]]
            
            // generate circuit for lower half (recursive call)
            CircGen c lowerbits (frombit)
            // same for upper half
            CircGen c higherbits (L+frombit)

    // do recursive splitting and increment
    CircGen c x 0

    // perform addition
    for i in 0..(n-1) do
        // Finally: Do "addition" on 1 bit (carry has been taken care of already)
        if ((c>>>i) &&& one) = one then
            BuildMultiplyControlledNOT ctrls [x.[i]] [g]

// Wrapper for circuit.compile
let ConstantAdderInplaceCompile (numctrls:int) (c:bigint) (xgc:Qubits) =
    let ctrls = slice xgc [(xgc.Length-numctrls)..(xgc.Length-1)]
    let xg = slice xgc [0..(xgc.Length-numctrls-1)]
    let x = slice xg [0..(xg.Length-3)]
    let g = xg.[xg.Length-2]
    let g2 = xg.[xg.Length-1]
    ConstantAdderInplace c x [g;g2] ctrls

// Wrapper for circuit.compile
let ConstantSubtractorInplaceCompile (numctrls:int) (c:bigint) (xgc:Qubits) =
    let ctrls = slice xgc [(xgc.Length-numctrls)..(xgc.Length-1)]
    let xg = slice xgc [0..(xgc.Length-numctrls-1)]
    let x = slice xg [0..(xg.Length-3)]
    let g = xg.[xg.Length-2]
    let g2 = xg.[xg.Length-1]
    X >< x
    ConstantAdderInplace c x [g;g2] ctrls
    X >< x


// Wrapper for circuit.compile
let BuildCarryUsingGarbage (c:bigint) (numctrls:int) (qs:Qubits) =
    let ctrls = slice qs [(qs.Length-numctrls)..(qs.Length-1)]
    let n = int (((qs.Length-numctrls)+2)/2)
    let x = slice qs [(n-2)..(2*n-3)]
    let g = slice qs [0..(n-3)]

    ComputeCarryUsingGarbage c g x ctrls

// Adds the constant c to the register x modulo N
// xgc consists of a qubit in |0> (labelled z), the x register, n-1 dirty qubits, and control qubits (# control qubits = numctrls)
// It uses the Takahashi trick for modular reduction, i.e. computes the final carry and, conditioned on that, either
// adds a or adds a-N
let ConstantAdderInplaceModN (numctrls:int) (c:bigint) (N:bigint) (zxgc:Qubits) =
    let n = int ((zxgc.Length-numctrls)/2)
    let x = slice zxgc [1..n]
    let g = slice zxgc [(n+1)..(2*n-1)]
    let ctrls = slice zxgc [(2*n)..(2*n+numctrls-1)]
    let z = [zxgc.[0]]
    let one = bigint 1
    let xz = List.concat [x;z]

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage (N-c) g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

    ConstantAdderInplace c x [g.[0];g.[1]] z
    BuildMultiplyControlledNOT ctrls z [g.[0]]
    
    for i in 0..(n-1) do
        X [x.[i]]
    ConstantAdderInplace (N-c) x [g.[0];g.[1]] z
    ComputeCarryUsingGarbage c g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

// Subtracts the constant c from the register x modulo N
// xgc consists of a qubits in |0> (labelled z), the x register, n-1 dirty qubits, and control qubits (# control qubits = numctrls)
// It uses the Takahashi trick for modular reduction, i.e. computes the final carry and, conditioned on that, either
// adds a or adds a-N
let ConstantSubtractorInplaceModN (numctrls:int) (c:bigint) (N:bigint) (zxgc:Qubits) =
    let n = int ((zxgc.Length-numctrls)/2)
    let x = slice zxgc [1..n]
    let g = slice zxgc [(n+1)..(2*n-1)]
    let ctrls = slice zxgc [(2*n)..(2*n+numctrls-1)]
    let z = [zxgc.[0]]
    let one = bigint 1
    let xz = List.concat [x;z]

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage c g xz ctrls
    ConstantAdderInplace ((one<<<n)-(N-c)) x [g.[0];g.[1]] z
    for i in 0..(n-1) do
        X [x.[i]]

    BuildMultiplyControlledNOT ctrls z [g.[0]]
    ConstantAdderInplace ((one<<<n)-c) x [g.[0];g.[1]] z

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage (N-c) g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

///////////////////////////////////////////////////////////////////////////
// Integer multiplication
///////////////////////////////////////////////////////////////////////////

let BuildIntegerMultiplication (xs:Qubits) (ys:Qubits) (rs:Qubits) = 
    let n = xs.Length
    if not (xs.Length = ys.Length) then
        failwith "Register sizes need to be equal to multiply!"
    
    for i in 0..(n-1) do
        BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] ys.[0..(n-i-1)] [xs.[i]]

let IntegerMultiplication (qs:Qubits) =
    let n = qs.Length/3
    BuildIntegerMultiplication qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1]
    ()

let BuildCtrlIntegerMultiplication (xs:Qubits) (ys:Qubits) (rs:Qubits) (cs:Qubits) = 
    let n = xs.Length
    if not (xs.Length = ys.Length) then
        failwith "Register sizes need to be equal to multiply!"
    
    for i in 0..(n-1) do
        BuildMultiCtrlTakahashiModAdder rs.[i..(rs.Length-1)] ys.[0..(n-i-1)] ([xs.[i]] @ cs)

let CtrlIntegerMultiplication (qs:Qubits) =
    let n = (qs.Length-1)/3
    BuildCtrlIntegerMultiplication qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1] qs.[3*n..3*n]
    ()

let BuildFixedPointMultiplication (pointpos:int) (xs:Qubits) (ys:Qubits) (rs:Qubits) = 
    let n = xs.Length
    if not (xs.Length = ys.Length) then
        failwith "Register sizes need to be equal to multiply!"

    // Do the additions that require right-shifts
    for i in 0..(n-pointpos-1) do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = ys.[ystart..yend]
        
        let carry = rs.[numitems]
        if numitems > 0 then
            BuildCtrlTakahashiAdderSmall rs.[0..numitems-1] (List.concat [addend; [carry]]) [xs.[i]] [xs.[(i+1)%xs.Length]]

    // Do the addition that require left-shifts
    for i in 0..(pointpos-1) do
        let numitems = n-i
        let yend = numitems-1
        let addend = ys.[0..yend]

        if i = pointpos-1 then
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]
        else
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]

let FixedPointMultiplication (qs:Qubits) =
    let n = qs.Length/3
    BuildIntegerMultiplication qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1]
    ()


///////////////////////////////////////////////////////////////////////////
// Integer division with remainder
///////////////////////////////////////////////////////////////////////////

// Performs integer division, i.e. calculates x / y
// Takes x,y,0 to r,y,q where r denotes the remainder and q the quotient
let BuildIntegerDivision (xs:Qubits) (ys:Qubits) (rs:Qubits) =
    let n = xs.Length

    for i in xs.Length-1..(-1)..0 do
        let anc = slice rs [0..i-1]
        let subtr = slice xs [i..xs.Length-1]
        let subtraction_reg = List.concat [subtr;anc]
        BuildTakahashiAdderInverse subtraction_reg (List.concat [ys;[rs.[i]]])
        BuildCtrlTakahashiModAdder subtraction_reg ys [rs.[i]]
        X [rs.[i]]

let IntegerDivision (qs:Qubits) =
    let n = qs.Length/3
    BuildIntegerDivision qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1]
    ()

let BuildCtrlIntegerDivision (xs:Qubits) (ys:Qubits) (rs:Qubits) (cs:Qubits) =
    let n = xs.Length

    for i in xs.Length-1..(-1)..0 do
        let anc = slice rs [0..i-1]
        let subtr = slice xs [i..xs.Length-1]
        let subtraction_reg = List.concat [subtr;anc]
        BuildCtrlTakahashiAdderInverse subtraction_reg (List.concat [ys;[rs.[i]]]) cs
        BuildMultiCtrlTakahashiModAdder subtraction_reg ys ([rs.[i]] @ cs)
        X [rs.[i]]

let CtrlIntegerDivision (qs:Qubits) =
    let n = (qs.Length-1)/3
    BuildCtrlIntegerDivision qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1] qs.[3*n..3*n]
    ()


///////////////////////////////////////////////////////////////////////////
//
// Quantum modular arithmetic
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Modular adders
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

let BuildModularAdderConstantPrime (ms:bigint) (xs:Qubits) (ys:Qubits) (tmp:Qubits) =
    // Computes the sum x + y mod m, for constant m.
    let n = xs.Length
    // name: || xs | ys  | tmp ||
    // bits: || n  | n+1 | 11  ||
    // init: || x  | y 0 | gg  ||
    let ctrls = []
                    
    BuildTakahashiAdder xs ys
    X >< ys
    ConstantAdderInplace ms ys tmp ctrls
    X >< ys 
    ConstantAdderInplace ms (slice ys [0..(n-1)]) tmp ([ys.[n]] @ ctrls)
    BuildCtrlStrictlyLargerThanComparator xs ys ctrls
    
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

let BuildCtrlModularAdderConstantPrime (ms:bigint) (xs:Qubits) (ys:Qubits) (tmp:Qubits) (ctrls:Qubits) =
    // Computes the sum x + y mod m, for constant m, conditioned on the qubits in ctrls all being 1.
    let n = xs.Length
    // name: || xs | ys  | tmp ||
    // bits: || n  | n+1 | 11  ||
    // init: || x  | y 0 | gg  ||
    
    BuildCtrlTakahashiAdder xs ys ctrls
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms (xs @ [ys.[n]]) tmp ctrls
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms xs tmp ([ys.[n]] @ ctrls)
    BuildCtrlStrictlyLargerThanComparator xs ys ctrls
    X [ys.[n]]
      
let BuildCtrlNegator (xs:Qubits) (ms:Qubits) (c:Qubits) =
    let n = xs.Length
    
    BuildCtrlTakahashiModAdderInverse ms xs c
    BuildCtrlTakahashiModAdder xs ms c
    for i in 0..(n-1) do 
        CNOT  [xs.[i]; ms.[i]]
        CCNOT [c.[0]; ms.[i]; xs.[i]]
        CNOT  [xs.[i]; ms.[i]]

let ModularAdder (qs:Qubits) =
    let n = (qs.Length-2)/3
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

let ModularAdderConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-3)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let tmp = slice qs [(2*n+1)..(2*n+2)]
    BuildModularAdderConstantPrime ms xs ys tmp

let CtrlModularAdderConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-4)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)]
    let tmp = slice qs [(2*n+1)..(2*n+2)]
    let c = slice qs [(2*n+3)]
    BuildCtrlModularAdderConstantPrime ms xs ys tmp c

let ModADD (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile ModularAdder qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let CtrlModADD (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile CtrlModularAdder qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
let ModADDConstantPrime (verbose:bool) (ms:bigint) (qs:Qubits)  = 
    let a = Circuit.Compile (ModularAdderConstantPrime ms) qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let CtrlModADDConstantPrime (verbose:bool) (ms:bigint) (qs:Qubits)  = 
    let a = Circuit.Compile (CtrlModularAdderConstantPrime ms) qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let CtrlNeg (verbose:bool) (qs:Qubits) = 
    let a = Circuit.Compile CtrlNegator qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
///////////////////////////////////////////////////////////////////////////
// Functions to run modular adders
///////////////////////////////////////////////////////////////////////////

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
        | "ModADDConstantPrime" -> 
            // name: || xs | ys | tmp ||
            // bits: || n  | n  | 11  ||
            // init: || x  | y  | 00  ||
            let k = Ket(2*n+2)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore            
            qs
        | "CtrlModADD" -> 
            // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | tmp | modulus | tmp | ctrl ||
            // bits: || n  | n  | 11  | n       | 1   | 1    ||
            // init: || x  | y  | 00  | p       | 0   | c    ||
            let k = Ket(3*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore
            m |> BoolInt n |> PrepBool qs (2*n+1) 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs (3*n+3) 1 |> ignore
            qs
        | "CtrlModADDConstantPrime" -> 
            // name: || xs | ys | tmp || ctrl ||
            // bits: || n  | n  | 111 || 1    ||
            // init: || x  | y  | 000 || c    ||
            let k = Ket(2*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore            
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
        | "ModADDConstantPrime" -> 
            let ADD = ModADDConstantPrime verbose1 m
            let rslt = Array.create (2*n+2) 0
            if verbose1 then 
                show "Number of qubits = %A" (2*n+2)
            ADD, rslt
        | "CtrlModADD" ->
            let ADD = CtrlModADD verbose1
            let rslt = Array.create (3*n+4) 0   // used to store final result vector after measurement 
            if verbose1 then 
                show "Number of qubits = %A" (3*n+4)
            ADD, rslt
        | "CtrlModADDConstantPrime" -> 
            let ADD = CtrlModADDConstantPrime verbose1 m
            let rslt = Array.create (2*n+4) 0
            if verbose1 then 
                show "Number of qubits = %A" (2*n+4)
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

    // Use the following function to print basic circuit information (qubits, size, depth)
    PrintToffoliNetworkMetrics ModularAdderCircuit
    
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
        | "ModADD" | "CtrlModADD" | "ModADDConstantPrime" -> finalState.[n..(2*n-1)] 
        | "CtrlNeg" -> finalState.[0..(n-1)]
        | _ -> failwith "Unknown modular adder type"
    // printfn "Final result = %A" res
    let resInt = (IntBool res.Length res)
    // printfn "As a number mod p this is = %A" resInt
    resInt


///////////////////////////////////////////////////////////////////////////
// Functions to compile modular adders
///////////////////////////////////////////////////////////////////////////

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
// Modular doubling 
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
// Modular multiplication: Montgomery multiplication
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
// Modular multiplication: Montgomery squarer
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
// Modular multiplication: modulo a Mersenne number 2^n-1
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
// Modular squaring modulo a Mersenne number 2^n-1
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
// Functions to run modular multipliers
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

///////////////////////////////////////////////////////////////////////////
// Functions to compile modular multipliers
///////////////////////////////////////////////////////////////////////////

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
// Modular inversion: standard encoding, extended Euclidean algorithm
///////////////////////////////////////////////////////////////////////////

#if FALSE

let BuildModularInverseForward (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits) = 
    let EuclidRound (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (fs:Qubits) (quo:Qubits) (ps:Qubits) (ks:Qubits)  = 

        // step 1: check termination condition: is vs the all zero string
        X >< vs
        BuildMultiplyControlledNOT vs [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X >< vs
        // step 2: check counter state and flip flag
        BuildMultiplyControlledNOT ks [fs.[0]] [ss.[0]] // using first bit of ss as dirty ancilla 
        X fs
        //
        // NEED to add control!!!!!! (on fs.[0])
        //
        // step 3: perform the division with remainder us = vs * quo + rem
        BuildCtrlIntegerDivision us vs quo [fs.[0]]
        // step 4: swap us <-> vs and rs <-> ss
        SWAP us vs
        SWAP rs ss
        // step 5: update the b register depending on whether round number is even or odd 
        BuildCtrlIntegerMultiplication quo ss ps [fs.[0]]
        BuildCtrlTakahashiModAdder rs ps [fs.[0]]
        // step 5: increase counter if required
        X fs
        BuildControlShiftRegisterCounterSpecial ks fs 
        X fs

    let n = us.Length
    for i in 0..(n-1) do 
        EuclidRound us vs ss rs fs a ks 

#endif
    
///////////////////////////////////////////////////////////////////////////
// Modular inversion: time optimized inversion with 14n qubits
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
        | "Standard" -> 
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
// Modular inversion: space-optimized inversion with 10n qubits
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
    
///////////////////////////////////////////////////////////////////////////
// Functions to run modular inverters
///////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////
// Functions to compiler modular inverters
///////////////////////////////////////////////////////////////////////////

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
// Shor's algorithm for factoring and dlog
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Shor's algorithm: factoring using n^3 gates and MG encoding and 6n qubits
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
// Shor's algorithm: factoring in n^3 log n gates and 2n+2 qubits (Th. Haener)
///////////////////////////////////////////////////////////////////////////

// (Controlled) Modular multiplication of a quantum number x by a classical constant c
// Transforms |x>|y=0>|0>|ctrl> to |c*x mod N>|y=0>|0>|ctrl> if ctrl = 1...1
let ConstantMultiplierModN (c:bigint) (cinvModN:bigint) (N:bigint) (x:Qubits) (y:Qubits) (zero:Qubits) (ctrls:Qubits) =
    let n = x.Length
    let one = bigint 1
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        ConstantAdderInplaceModN (1+ctrls.Length) ((c<<<i)%(N)) N (List.concat [zero;y;freexs;ctrls;[x.[i]]])

    CSWAP x y ctrls

    // uncompute y
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        ConstantSubtractorInplaceModN (1+ctrls.Length) ((cinvModN<<<i)%(N)) N (List.concat [zero;y;freexs;ctrls;[x.[i]]])

// Same as the above, use with Circuit.Compile (is just a wrapper)
let ConstantMultiplierModNCompile (ccount:int) (a:bigint) (ainv:bigint) (N:bigint) (qs:Qubits) =
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]

    ConstantMultiplierModN a ainv N x y zero ctrls

// This does the same as the functions above but RUNS the circuit
// Use this for large bit-sizes, as the circuits get too large otherwise.
// It generates the circuit for each modular addition, executes it, counts the Toffolis, and then does the next one.
let ConstantMultiplierModNRun (ccount:int) (a:bigint) (ainv:bigint) (N:bigint) (qs:Qubits) (initialState:int[]) = 
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]
    let one = bigint 1

    let mutable currentState = initialState
    let mutable toffcount = 0
    let mutable depth = 0

    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        let circuit = Circuit.Compile (ConstantAdderInplaceModN (1+ctrls.Length) ((a<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) |> MyCircuitExport
        currentState <- MyCircuitSimulateFast circuit currentState
        toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
        //show "Compute: %A / %A" i n

    let circuit = Circuit.Compile (CSWAPCompile (ctrls.Length)) (List.concat [x; y; ctrls]) |> MyCircuitExport
    currentState <- MyCircuitSimulateFast circuit currentState
    toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    
    // uncompute y
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        let circuit = Circuit.Compile (ConstantSubtractorInplaceModN (1+ctrls.Length) ((ainv<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) |> MyCircuitExport
        currentState <- MyCircuitSimulateFast circuit currentState
        toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
        //show "Uncompute: %A / %A" i n

    show "N=%A TC=%A" n toffcount
    currentState.[0..(2*n+ccount)]


let ConstantMultiplierModNRunWholeCircuit (ccount:int) (a:bigint) (ainv:bigint) (n:int) (N:bigint) (qs:Qubits) = 
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]
    let one = bigint 1

    let totalCirc1 = 
        [ for i in 0..(n-1) do
            let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
            let circuit = Circuit.Compile (ConstantAdderInplaceModN (1+ctrls.Length) ((a<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) 
            yield circuit ]
        |> Seq    
                
    let totalCirc2 = 
        Circuit.Compile (CSWAPCompile (ctrls.Length)) (List.concat [x; y; ctrls])
    
    let totalCirc3 = 
        [ for i in 0..(n-1) do
            let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
            let circuit = Circuit.Compile (ConstantSubtractorInplaceModN (1+ctrls.Length) ((ainv<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) 
            yield circuit ]
        |> Seq 
        
    let totalCirc = [ totalCirc1; totalCirc2; totalCirc3 ] |> Seq 
    let totalToffCirc = MyCircuitExport totalCirc
           
    PrintToffoliNetworkMetrics totalToffCirc 
//    let compressed1 = totalCirc.Fold(false) // this is too slow with the old circ.fold even for n = 8
//    let compressed2 = totalCirc.Fold(true)  // this is too slow with the old circ.fold even for n = 8
//    let count0 = totalCirc.GateCount(false)
//    let count1 = totalCirc.GateCount(true)
//    let countSeq = compressed1.GateCount(false, fun x -> x.Name = "CCNOT")
//    let countPar = compressed1.GateCount(true, fun x -> x.Name = "CCNOT")
//
//    totalCirc.RenderHT("circ0")
//    compressed1.RenderHT("circ1")
//
//    show "w/o fold: Gate count serial %A" count0
//    show "w/o fold: Gate count parallel %A" count1
//    show "Gate count serial %A" countSeq
//    show "Gate count parallel %A" countPar
//    show "done"


///////////////////////////////////////////////////////////////////////////
//
// Quantum ECC point addition for affine Weierstrass curves
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Shor's algorithm: dlog over ECC curves
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
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
open CircBase           // basic quantum circuits
open Integer            // integer arithmetic

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
                    
    BuildTakahashiAdder xs ys 
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms (xs @ [ys.[n]]) tmp []
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms xs tmp ([ys.[n]])
    BuildStrictlyLargerThanComparator xs ys 
    X [ys.[n]] 
        
let BuildModularAdderConstantPrimeConstantNumber (ms:bigint) (ss:bigint) (xs:Qubits) (tmp:Qubits) (gs:Qubits) =
    // Computes the sum x + s mod m, for constant m and constant s.
    // name: || xs | tmp | gs  ||
    // bits: || n  | 1   | n-1 ||
    // init: || x  | 0   | g   ||
    ConstantAdderInplaceModN 0 ss ms (tmp @ xs @ gs)

let BuildCtrlModularAdderConstantPrimeConstantNumber (ms:bigint) (ss:bigint) (xs:Qubits) (tmp:Qubits) (gs:Qubits) (cs:Qubits) =
    // Computes the sum x + s mod m, for constant m and constant s.
    // name: || xs | tmp | gs  | cs ||
    // bits: || n  | 1   | n-1 | 1  ||
    // init: || x  | 0   | g   | c  ||
    ConstantAdderInplaceModN 1 ss ms (tmp @ xs @ gs @ cs)

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
    // name: || xs | ys  | tmp | ctrl ||
    // bits: || n  | n+1 | 11  | 1    ||
    // init: || x  | y 0 | gg  | c    ||
    
    BuildCtrlTakahashiAdder xs ys ctrls
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms (xs @ [ys.[n]]) tmp ctrls
    X >< (xs @ [ys.[n]])
    ConstantAdderInplace ms xs tmp ([ys.[n]] @ ctrls)
    BuildCtrlStrictlyLargerThanComparator xs ys ctrls
    BuildMultiplyControlledNOT ctrls [ys.[n]] [xs.[0]] // using xs.[0] as dirty ancilla
      
let BuildCtrlNegator (xs:Qubits) (ms:Qubits) (c:Qubits) =
    let n = xs.Length
    
    BuildCtrlTakahashiModAdderInverse ms xs c
    BuildCtrlTakahashiModAdder xs ms c
    for i in 0..(n-1) do 
        CNOT  [xs.[i]; ms.[i]]
        CCNOT [c.[0]; ms.[i]; xs.[i]]
        CNOT  [xs.[i]; ms.[i]]

let BuildCtrlNegatorConstantPrime (ms:bigint) (xs:Qubits) (tmp:Qubits) (c:Qubits) =
    let n = xs.Length
    let flip = function | 0 -> 1 | 1 -> 0
    let ms2 = BoolInt n ms |> List.map flip |> List.toArray |> fun x -> IntBool n x
    ConstantAdderInplace ms2 xs tmp c
    X >< xs 
    
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

let CtrlNegatorConstantPrime (ms:bigint) (qs:Qubits) =
    let n = qs.Length-3
    let xs = slice qs [0..(n-1)]
    let c = slice qs [n]
    let ts = slice qs [(n+1)..(n+2)]
    BuildCtrlNegatorConstantPrime ms xs ts c

let ModularAdderConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-3)/2
    let ys = slice qs ([0..(n-1)] @ [2*n])
    let xs = slice qs [n..(2*n-1)] 
    let tmp = slice qs [(2*n+1)..(2*n+2)]
    BuildModularAdderConstantPrime ms xs ys tmp

let ModularAdderConstantPrimeConstantNumber (ms:bigint) (ss:bigint) (qs:Qubits) =
    let n = (qs.Length)/2
    let xs = slice qs [0..(n-1)] 
    let tmp = slice qs [n]
    let gs = slice qs [(n+1)..(2*n-1)]
    BuildModularAdderConstantPrimeConstantNumber ms ss xs tmp gs

let CtrlModularAdderConstantPrimeConstantNumber (ms:bigint) (ss:bigint) (qs:Qubits) =
    let n = (qs.Length-1)/2
    let xs = slice qs [0..(n-1)] 
    let tmp = slice qs [n]
    let gs = slice qs [(n+1)..(2*n-1)]
    let cs = slice qs [2*n]
    BuildCtrlModularAdderConstantPrimeConstantNumber ms ss xs tmp gs cs

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

let ModADDConstantPrimeConstantNumber (verbose:bool) (ms:bigint) (ss:bigint) (qs:Qubits)  = 
    let a = Circuit.Compile (ModularAdderConstantPrimeConstantNumber ms ss) qs 
    let gates = MyCircuitExport a
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        PrintToffoliNetworkMetrics gates
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         

let CtrlModADDConstantPrimeConstantNumber (verbose:bool) (ms:bigint) (ss:bigint) (qs:Qubits)  = 
    let a = Circuit.Compile (CtrlModularAdderConstantPrimeConstantNumber ms ss) qs 
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

let CtrlNegConstantPrime (verbose:bool) (ms:bigint) (qs:Qubits)  = 
    let a = Circuit.Compile (CtrlNegatorConstantPrime ms) qs 
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
            let k = Ket(2*n+3)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt n |> PrepBool qs n 1 |> ignore            
            qs
        | "ModADDConstantPrimeConstantNumber" -> 
            // name: || xs | tmp | gs  ||
            // bits: || n  | 1   | n-1 ||
            // init: || x  | 00  | g   ||
            let k = Ket(2*n)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
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
        | "CtrlNegConstantPrime" -> 
            // name: || xs | ctrl | anc ||
            // bits: || n  | 1    | 2   ||
            // init: || x  | c    | 00  ||
            let k = Ket(n+3)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            ctrl |> BoolInt 1 |> PrepBool qs n 1 |> ignore
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
            let rslt = Array.create (2*n+3) 0
            if verbose1 then 
                show "Number of qubits = %A" (2*n+3)
            ADD, rslt
        | "ModADDConstantPrimeConstantNumber" -> 
            let ADD = ModADDConstantPrimeConstantNumber verbose1 m s2
            let rslt = Array.create (2*n) 0
            if verbose1 then 
                show "Number of qubits = %A" (2*n)
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
        | "CtrlNegConstantPrime" -> 
            let ADD = CtrlNegConstantPrime verbose1 m
            let rslt = Array.create (n+3) 0   // used to store final result vector after measurement 
            if verbose1 then 
                show "Number of qubits = %A" (n+3)
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
        | "CtrlNeg" | "CtrlNegConstantPrime" | "ModADDConstantPrimeConstantNumber" -> finalState.[0..(n-1)]
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
         | "ModADDConstantPrime" -> 
            // name: || xs | ys  | tmp ||
            // bits: || n  | n+1 | 11  ||
            // init: || x  | y   | 00  ||
            let k = Ket(2*n+3)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (2*n+3)
            qs
         | "ModADDConstantPrimeConstantNumber" -> 
            // name: || xs | tmp | gs  ||
            // bits: || n  | 1   | n-1 ||
            // init: || x  | 0   | g   ||
            let k = Ket(2*n)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (2*n)
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
         | "CtrlModADDConstantPrime" -> 
            // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | tmp || ctrl ||
            // bits: || n  | n  | 111 || 1    ||
            // init: || x  | y  | 000 || c    ||
            let k = Ket(2*n+4)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (2*n+4)
            qs                                                 
         | "CtrlModADDConstantPrimeConstantNumber" -> 
            // name: || xs | tmp | gs  | cs ||
            // bits: || n  | 1   | n-1 | 1  ||
            // init: || x  | 0   | g   | c  ||
            let k = Ket(2*n+1)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (2*n+1)
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
         | "CtrlNegConstantPrime" -> 
            // Prepare the initial state for controlled modular negator. The partitioning of the quantum register is as follows: 
            // name: || xs | ctrl | anc ||
            // bits: || n  | 1    | 2   ||
            // init: || x  | c    | 0   ||
            let k = Ket(n+3)
            let qs = k.Qubits
            if verbose2 then 
                show "Number of qubits = %A" (n+3)
            qs
     
         | _ -> failwith "Unknown modular adder."
    let qs = arrangeAdderInputs      // arranges inputs in pattern needed for the adder
    
    let arrangeAdderCircuit = 
        match name with 
        | "ModADD" -> 
            let ADD = ModADD verbose1
            let rslt = Array.create (3*n+2) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "ModADDConstantPrime" -> 
            let ADD = ModADDConstantPrime verbose1 m
            let rslt = Array.create (2*n+3) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "ModADDConstantPrimeConstantNumber" -> 
            let ADD = ModADDConstantPrimeConstantNumber verbose1 m s2 
            let rslt = Array.create (2*n) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlModADD" ->
            let ADD = CtrlModADD verbose1
            let rslt = Array.create (3*n+4) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlModADDConstantPrime" ->
            let ADD = CtrlModADDConstantPrime verbose1 m
            let rslt = Array.create (2*n+4) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlModADDConstantPrimeConstantNumber" -> 
            let ADD = CtrlModADDConstantPrimeConstantNumber verbose1 m s2 
            let rslt = Array.create (2*n+1) 0   // used to store final result vector after measurement 
            ADD, rslt       
        | "CtrlNeg" ->
            let ADD = CtrlNeg verbose1
            let rslt = Array.create (2*n+1) 0   // used to store final result vector after measurement 
            ADD, rslt
        | "CtrlNegConstantPrime" ->
            let ADD = CtrlNegConstantPrime verbose1 m
            let rslt = Array.create (n+3) 0   // used to store final result vector after measurement 
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

let BuildModularDoublerConstantPrime (ms:bigint) (xs:Qubits) (a:Qubits) (tmp:Qubits) =
    // Computes the double 2*x mod m, for constant m.
    let n = xs.Length
    // name: || xs | a | tmp ||
    // bits: || n  | 1 | 11  ||
    // init: || x  | 0 | gg  ||
    let ctrls = []
                    
    BuildBinaryDoubling (xs @ [a.[0]]) // Double input by bit shift.
    // Subtract m to see whether modular reduction is necessary.
    X >< (xs @ [a.[0]])
    ConstantAdderInplace ms (xs @ [a.[0]]) tmp ctrls 
    X >< (xs @ [a.[0]]) 
    // Conditional addition of m, if result of subtraction was negative.
    ConstantAdderInplace ms xs tmp ([a.[0]] @ ctrls) 
    CNOT [xs.[0]; a.[0]]  // Uncompute the control bit by checking whether the LSB is zero. 
    X [a.[0]]

let BuildModularHalver (xs:Qubits) (ms:Qubits) (a:Qubits) =
    let n = xs.Length
    X [a.[1]]
    CNOT [xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    BuildCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) [a.[1]]
    CNOT [a.[0]; a.[1]] // Copy out control bit.
    BuildTakahashiAdder xs (ms @ [a.[0]])
    BuildBinaryHalfing (xs @ [a.[0]]) // Halving of input by bit shift.

let BuildModularHalverConstantPrime (ms:bigint) (xs:Qubits) (a:Qubits) (tmp:Qubits) =
    // Computes the half x/2 mod m, for constant m.
    let n = xs.Length
    // name: || xs | a | tmp ||
    // bits: || n  | 1 | 11  ||
    // init: || x  | 0 | gg  ||
    let ctrls = []
      
    X [a.[0]]  
    CNOT [xs.[0]; a.[0]]
    X >< xs   
    ConstantAdderInplace ms xs tmp ([a.[0]] @ ctrls)
    X >< xs 
    ConstantAdderInplace ms (xs @ [a.[0]]) tmp ctrls                   
    BuildBinaryHalfing (xs @ [a.[0]]) // Halve by bit shift.
    
let BuildCtrlModularDoubler (xs:Qubits) (ms:Qubits) (cs:Qubits) (a:Qubits) =
    let n = xs.Length
    BuildCtrlBinaryDoubling (xs @ [a.[0]]) cs [a.[1]] // Halving of input by bit shift. Using a.[1] as dirty ancilla.
    BuildCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) cs
    CCNOT [cs.[0]; a.[0]; a.[1]] // Copy out control bit.
    BuildMultiCtrlTakahashiAdder xs (ms @ [a.[0]]) [cs.[0]; a.[1]]
    CCNOT [cs.[0]; xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    CNOT [cs.[0]; a.[1]]

let BuildCtrlModularDoublerConstantPrime (ms:bigint) (xs:Qubits) (a:Qubits) (tmp:Qubits) (ctrls:Qubits) =
    // Computes the double 2*x mod m, for constant m.
    let n = xs.Length
    // name: || xs | a | tmp ||
    // bits: || n  | 1 | 11  ||
    // init: || x  | 0 | gg  ||
                    
    BuildCtrlBinaryDoubling (xs @ [a.[0]]) ctrls tmp // Double input by bit shift.
    // Subtract m to see whether modular reduction is necessary.
    X >< (xs @ [a.[0]])
    ConstantAdderInplace ms (xs @ [a.[0]]) tmp ctrls 
    X >< (xs @ [a.[0]]) 
    // Conditional addition of m, if result of subtraction was negative.
    ConstantAdderInplace ms xs tmp ([a.[0]] @ ctrls) 
    BuildMultiplyControlledNOT (ctrls @ [xs.[0]]) [a.[0]] tmp // using tmp as dirty ancilla
    BuildMultiplyControlledNOT (ctrls) [a.[0]] tmp // using tmp as dirty ancilla

let BuildCtrlModularHalver (xs:Qubits) (ms:Qubits) (cs:Qubits) (a:Qubits) =
    let n = xs.Length
    CNOT [cs.[0]; a.[1]]
    CCNOT [cs.[0]; xs.[0]; a.[1]]  // Uncompute the control bit by checking whether the LSB is zero. 
    BuildMultiCtrlTakahashiAdderInverse xs (ms @ [a.[0]]) [cs.[0]; a.[1]]
    CCNOT [cs.[0]; a.[0]; a.[1]] // Copy out control bit.
    BuildCtrlTakahashiAdder xs (ms @ [a.[0]]) cs
    BuildCtrlBinaryHalfing (xs @ [a.[0]]) cs [a.[1]] // Halving of input by bit shift. Using a.[1] as dirty ancilla.

let BuildCtrlModularHalverConstantPrime (ms:bigint) (xs:Qubits) (a:Qubits) (tmp:Qubits) (ctrls:Qubits) =
    // Computes the half x/2 mod m, for constant m.
    let n = xs.Length
    // name: || xs | a | tmp ||
    // bits: || n  | 1 | 11  ||
    // init: || x  | 0 | gg  ||
      
    BuildMultiplyControlledNOT ctrls [a.[0]] tmp // using tmp as dirty ancilla
    BuildMultiplyControlledNOT (ctrls @ [xs.[0]]) [a.[0]] tmp // using tmp as dirty ancilla
    X >< xs   
    ConstantAdderInplace ms xs tmp ([a.[0]] @ ctrls)
    X >< xs 
    ConstantAdderInplace ms (xs @ [a.[0]]) tmp ctrls                   
    BuildCtrlBinaryHalfing (xs @ [a.[0]]) ctrls tmp // Conditionally halve by bit shift.

let ModularDoubler (qs:Qubits) =
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ms = slice qs [n..(2*n-1)]
    let a = slice qs [2*n..(2*n+1)]
    BuildModularDoubler xs ms a

let ModularDoublerConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-3)
    let xs = slice qs [0..n-1]
    let a = slice qs [n]
    let tmp = slice qs [(n+1)..(n+2)]
    BuildModularDoublerConstantPrime ms xs a tmp

let CtrlModularDoublerConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-4)
    let xs = slice qs [0..n-1]
    let a = slice qs [n]
    let tmp = slice qs [(n+1)..(n+2)]
    let c = slice qs [(n+3)]
    BuildCtrlModularDoublerConstantPrime ms xs a tmp c

let ModularHalver (qs:Qubits) =
    let n = (qs.Length-2)/2
    let xs = slice qs [0..n-1]
    let ms = slice qs [n..(2*n-1)]
    let a = slice qs [2*n..(2*n+1)]
    BuildModularHalver xs ms a

let ModularHalverConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-3)
    let xs = slice qs [0..n-1]
    let a = slice qs [n]
    let tmp = slice qs [(n+1)..(n+2)]
    BuildModularHalverConstantPrime ms xs a tmp

let CtrlModularHalverConstantPrime (ms:bigint) (qs:Qubits) =
    let n = (qs.Length-4)
    let xs = slice qs [0..n-1]
    let a = slice qs [n]
    let tmp = slice qs [(n+1)..(n+2)]
    let c = slice qs [(n+3)]
    BuildCtrlModularHalverConstantPrime ms xs a tmp c

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
// Modular multiplication: Montgomery multiplication (constant prime) 
///////////////////////////////////////////////////////////////////////////

let BuildMontgomeryMultiplierForwardConstantPrime (ms:bigint) (xs:Qubits) (ys:Qubits) (acc:Qubits) (mg:Qubits) = 
    let n = xs.Length
    
    let AdderRound (ms:bigint) (a:Qubits) (ys:Qubits) (acc:Qubits) (b:Qubits) = 
        let tmp = [ys.[0]; ys.[1]] // need 2 dirty ancillas
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ys @ [acc.[n+1]]) a // invariant: MSB of acc is equal to 0
        CNOT [acc.[0]; b.[0]] 
        //BuildCtrlTakahashiAdder (slice acc [0..n]) (ms @ [acc.[n+1]]) b // invariant: LSB of acc is equal to 0
        ConstantAdderInplace ms acc tmp b // invariant: LSB of acc is equal to 0
    
    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+1]) @ (slice acc [0..(i-1)]) // reindex the accumulator so we never have to rewire
        AdderRound ms [xs.[i]] ys acci [mg.[i]] 
   
let MontgomeryMultiplierForwardConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-4)/5
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n)]
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let mg = slice qs [(3*n+4)..(4*n+3)]
    BuildMontgomeryMultiplierForwardConstantPrime ms xs ys acc mg
           
let MontgomeryConditionalSubtractionConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-4)/5
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let tmp = [qs.[0]; qs.[1]] // need 2 dirty ancillas
    // check if final result in the accumulator overflows and if so, then subtract the modulus once
    //BuildTakahashiAdderInverse (slice acci [0..n]) (ms @ [acci.[n+1]]) 
    //BuildCtrlTakahashiModAdder (slice acci [0..n]) ms [acci.[n+1]]
    X >< acci
    ConstantAdderInplace ms acci tmp []
    X >< acci 
    ConstantAdderInplace ms (slice acci [0..(n-1)]) tmp [acci.[n]]
        
let CopyCircuitConstantPrime (qs:Qubits) = 
    let n, flag = 
        match ((qs.Length-5) % 5) with 
        | 0 -> (qs.Length-5)/5, true // this is a control copy circuit conditioned on the last qubit
        | _ -> (qs.Length-4)/5, false// this is a regular copy circuit
    let cs = slice qs [(qs.Length-1)..(qs.Length-1)]
    let acc = slice qs [(2*n+1)..(3*n+2)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let res = slice qs [(4*n+4)..(5*n+3)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) (cs:Qubits) = 
        for i in 0..(n-1) do 
           match flag with 
           | true -> CCNOT [cs.[0]; acc.[i]; res.[i]]        
           | _    ->  CNOT [acc.[i]; res.[i]]        
    CopyRegs acci res cs

let MontgomeryMultiplierFastConstantPrime (verbose:bool) (ms:bigint) (qs:Qubits) = 
    let a = Circuit.Compile (MontgomeryMultiplierForwardConstantPrime ms) qs |> MyCircuitExport
    let b = Circuit.Compile (MontgomeryConditionalSubtractionConstantPrime ms) qs |> MyCircuitExport
    let c = Circuit.Compile CopyCircuitConstantPrime qs |> MyCircuitExport
    let d = List.rev b // Run the final modular reduction in reverse
    let e = List.rev a // Run the sequence of Montgomery steps in reverse
    let gates = [a; b; c; d; e] |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    
    if verbose then 
        PrintToffoliNetworkMetrics gates
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
// Modular multiplication: Montgomery squarer (constant prime)
///////////////////////////////////////////////////////////////////////////

let BuildMontgomerySquarerForwardConstantPrime (ms:bigint) (xs:Qubits) (acc:Qubits) (mg:Qubits) = 
    let n = xs.Length-2    

    let AdderRound (ms:bigint) (a:Qubits) (ys:Qubits) (acc:Qubits) (b:Qubits) = 
        let tmp = [ys.[0]; ys.[1]] // need 2 dirty ancillas
        BuildCtrlTakahashiAdder (slice acc [0..n]) (ys @ [acc.[n+1]]) a // invariant: MSB of acc is equal to 0
        CNOT [acc.[0]; b.[0]] 
        ConstantAdderInplace ms acc tmp b // invariant: LSB of acc is equal to 0
    
    for i in 0..(n-1) do 
        let acci = (slice acc [i..n+1]) @ (slice acc [0..(i-1)]) // reindex the accumulator so we never have to rewire
        CNOT [xs.[i]; xs.[n+1]]
        AdderRound ms [xs.[n+1]] (slice xs [0..n]) acci [mg.[i]] 
        CNOT [xs.[i]; xs.[n+1]]

let MontgomerySquarerForwardConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-5)/4
    let xs = slice qs [0..(n+1)] // n bits of input x and 2 ancilla bits
    let acc = slice qs [(n+2)..(2*n+3)]
    let mg = slice qs [(2*n+5)..(3*n+4)]
    BuildMontgomerySquarerForwardConstantPrime ms xs acc mg
           
let MontgomerySquarerConditionalSubtractionConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-5)/4
    let acc = slice qs [(n+2)..(2*n+3)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let tmp = [qs.[0]; qs.[1]] // need 2 dirty ancillas
    // check if final result in the accumulator overflows and if so, then subtract the modulus once
    X >< acci
    ConstantAdderInplace ms acci tmp []
    X >< acci 
    ConstantAdderInplace ms (slice acci [0..(n-1)]) tmp [acci.[n]]
    
let CopySquarerCircuitConstantPrime (qs:Qubits) = 
    let n = (qs.Length-5)/4
    let acc = slice qs [(n+2)..(2*n+3)]
    let acci = (slice acc [n..(n+1)]) @ (slice acc [0..(n-1)]) // reindex the accumulator so we never have to rewire        
    let res = slice qs [(3*n+5)..(4*n+4)]
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (acc:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [acc.[i]; res.[i]]        
    CopyRegs acci res
 
let MontgomerySquarerConstantPrime (verbose:bool) (ms:bigint) (qs:Qubits) = 
    let a = Circuit.Compile (MontgomerySquarerForwardConstantPrime ms) qs
    let b = Circuit.Compile (MontgomerySquarerConditionalSubtractionConstantPrime ms) qs
    let c = Circuit.Compile (CopySquarerCircuitConstantPrime) qs
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

let BuildModDBLADDMultiplierConstantPrime (ms:bigint) (xs:Qubits) (ys:Qubits) (acc:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..1 do 
        BuildCtrlModularAdderConstantPrime ms (slice acc [0..(n-1)]) (ys @ slice acc [n]) tmp [xs.[i]]
        BuildModularDoublerConstantPrime ms (slice acc [0..(n-1)]) [acc.[n]] tmp
    BuildCtrlModularAdderConstantPrime ms (slice acc [0..(n-1)]) (ys @ [acc.[n]]) tmp [xs.[0]] 

let ModDblAddMultiplierConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-3)/3
    let xs = slice qs [0..(n-1)]
    let ys = slice qs [n..(2*n-1)]
    let acc = slice qs [(2*n)..(3*n)]
    let tmp = slice qs [(3*n+1)..(3*n+2)]
    BuildModDBLADDMultiplierConstantPrime ms xs ys acc tmp

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

let BuildModDBLADDSquarerConstantPrime (ms:bigint) (xs:Qubits) (acc:Qubits) (tmp:Qubits) = 
    let n = xs.Length
    
    for i in (n-1)..(-1)..1 do
        CNOT [xs.[i]; tmp.[2]]
        BuildCtrlModularAdderConstantPrime ms (slice acc [0..(n-1)]) (xs @ slice acc [n]) tmp [tmp.[2]]
        CNOT [xs.[i]; tmp.[2]]
        BuildModularDoublerConstantPrime ms (slice acc [0..(n-1)]) [acc.[n]] tmp
    CNOT [xs.[0]; tmp.[2]]
    BuildCtrlModularAdderConstantPrime ms (slice acc [0..(n-1)]) (xs @ [acc.[n]]) tmp [tmp.[2]]
    CNOT [xs.[0]; tmp.[2]]

let ModDblAddSquarerConstantPrime (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-5)/2
    let xs = slice qs [0..(n-1)]
    let acc = slice qs [(n)..(2*n)]
    let tmp = slice qs [(2*n+1)..(2*n+4)]
    BuildModDBLADDSquarerConstantPrime ms xs acc tmp

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
        | "MontgomeryConstantPrime" -> 
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | n  | 1       | n+2  | 1       | n         | n      ||
            // init: || x  | y  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(5*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
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
        | "MontgomerySquarerConstantPrime" -> 
            // Prepare the initial state for Montgomery squarer. The partitioning of the quantum register is as follows: 
            // name: || xs | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | 2       | n+2  | 1       | n         | n      ||
            // init: || x  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(4*n+5)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
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
        | "MontgomeryConstantPrime" -> 
            let MUL = MontgomeryMultiplierFastConstantPrime verbose1 m
            MUL
        | "MontgomerySquarer" -> 
            let MUL = MontgomerySquarer verbose1 
            MUL
        | "MontgomerySquarerConstantPrime" -> 
            let MUL = MontgomerySquarerConstantPrime verbose1 m 
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
        | "MontgomeryConstantPrime" -> 
            finalState.[(4*n+4)..(5*n+3)]
        | "MontgomerySquarer" -> 
            finalState.[(4*n+5)..(5*n+4)]
        | "MontgomerySquarerConstantPrime" -> 
            finalState.[(3*n+5)..(4*n+4)]
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


let RunMultiplierWholeCircuit (name:string) (n:int) (s1:bigint) (s2:bigint) (m:bigint) (verbosity:string) = 
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
        | "MontgomeryConstantPrime" -> 
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | n  | 1       | n+2  | 1       | n         | n      ||
            // init: || x  | y  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(5*n+4)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
            s2 |> BoolInt (n+1) |> PrepBool qs n 1 |> ignore
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
        | "MontgomerySquarerConstantPrime" -> 
            // Prepare the initial state for Montgomery squarer. The partitioning of the quantum register is as follows: 
            // name: || xs | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | 2       | n+2  | 1       | n         | n      ||
            // init: || x  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(4*n+5)
            let qs = k.Qubits
            s1 |> BoolInt n |> PrepBool qs 0 1 |> ignore
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
        | "MontgomeryConstantPrime" -> 
            let MUL = MontgomeryMultiplierFastConstantPrime verbose1 m
            MUL
        | "MontgomerySquarer" -> 
            let MUL = MontgomerySquarer verbose1 
            MUL
        | "MontgomerySquarerConstantPrime" -> 
            let MUL = MontgomerySquarerConstantPrime verbose1 m 
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
    
    let q, s, d = ReportToffoliNetworkMetrics MultiplierCircuit
    n, q, s, d

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
        | "MontgomeryConstantPrime" -> 
            // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
            // name: || xs | ys | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | n  | 1       | n+2  | 1       | n         | n      ||
            // init: || x  | y  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(5*n+4)
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
        | "MontgomerySquarerConstantPrime" -> 
            // Prepare the initial state for Montgomery squarer. The partitioning of the quantum register is as follows: 
            // name: || xs | ancilla | acc  | ancilla | mg-rounds | result || 
            // bits: || n  | 2       | n+2  | 1       | n         | n      ||
            // init: || x  | 0       | 0    | 0       | 0..0      | 0      || 
            let k = Ket(4*n+5)
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
        | "MontgomeryConstantPrime" -> 
            let MUL = MontgomeryMultiplierFastConstantPrime verbose1 m
            MUL
        | "CtrlMontgomery" -> 
            let MUL = CtrlMontgomeryMultiplier verbose1 
            MUL
        | "MontgomerySquarer" -> 
            let MUL = MontgomerySquarer verbose1 
            MUL
        | "MontgomerySquarerConstantPrime" -> 
            let MUL = MontgomerySquarerConstantPrime verbose1 m 
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
        PrintToffoliNetworkMetrics MultiplierCircuit
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
    
    
///////////////////////////////////////////////////////////////////////////
// Modular inversion: space-optimized inversion with 10n qubits
///////////////////////////////////////////////////////////////////////////

// BuildMontgomeryInverseForward implements the inverse mod p of numbers that are given in Montgomery foom x * 2^n mod p. The output is 
// again in Montogomery form. The implementation follows closely the Kalinski paper which in turn is a binary extended Euclidean algorithm. 
// The implementation uses a total of 14n+6+2*ceil(log(n)) qubits. The registers are: 
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
        let round = (slice mg [2*i..2*i+1]) // get 2 bits per round for match statement on (u, v)
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
    show "INVERTER: Time to concatenate: %f" watch.Elapsed.TotalSeconds
    // Compute the number of Toffoli gates in the circuit and write to console
    PrintToffoliNetworkMetrics gates
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
    show "INVERTER: Time to concatenate: %f" watch.Elapsed.TotalSeconds
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         


///////////////////////////////////////////////////////////////////////////
// Super efficient MG inverse (8n qubits; modulus still there)
///////////////////////////////////////////////////////////////////////////

// BuildSuperEfficientMontgomeryInverseForward implements the inverse mod p of numbers that are given in Montgomery foom x * 2^n mod p. The output is 
// again in Montogomery form. The implementation follows closely the Kalinski paper which in turn is a binary extended Euclidean algorithm. 
// The implementation uses a total of 8n+6+2*ceil(log(n)) qubits. The registers are: 
// us = intially equal to modulus p, vs = intitially equal to input x, ss = intially 1, rs = intiially 0, ms = modulus p, mg = 8n qubits 
// holding the results of the dispatch predicates for each round, fs = flag that determines whether we are in compute or in counter mode, 
// ks = a register of qubits to implement a counter from 0 up to worst case iteration bound of 2n, and a = 1 ancilla. 
// The method uses the fact that the rs register can be cleaned up in each round, eliminating 2n qubits from the Efficient MG Inverse circuit. 
// This method still uses the extra register to store the modulus. 

let BuildSuperEfficientMontgomeryInverseForward (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (mg:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits) = 
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
        // now, use the more efficient encoding: ue -> 01 (MSB/LSB), ve -> 10, ug -> 11, ul -> 00 which can 
        // be computed with 4 Toffoli gates and 1 ancilla which is returned clean after the case statement.
        CCNOT [fs.[0]; us.[0]; a.[2]] // if ue=1 then a.[2] round.[0] = 10
        X us 
        // step 3b: or else, check whether v is even
        X vs 
        X [a.[2]]
        BuildMultiplyControlledNOT [fs.[0]; a.[2]; vs.[0]] [round.[0]] [ss.[0]] // if ve=1 and ue = 0 then a.[2] round.[0] = 01  
        X [a.[2]]
        X vs
        // step 3c: or else, check whether u > v
        BuildCtrlTakahashiAdderInverse vs (us @ a) [fs.[0]] 
        CNOT [a.[2]; a.[1]]
        CNOT [round.[0]; a.[1]] // compute partity of a.[2] round.[0] into ancilla a.[1]
        X [a.[1]] // apply following gate only if ancilla is 0
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [a.[2]] [ss.[0]]
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [round.[0]] [ss.[0]] // if ue=ve=0 and ug=1 then a.[2] round.[0]= 11
        X [a.[1]] 
        CNOT [a.[2]; a.[1]]
        CNOT [round.[0]; a.[1]] // clean up ancilla a.[1] after each round. Note that parity is preseved in all cases        
        // step 3d: else, we have that u <= v
        // note: no computation needed in this case: ue=ve=ug=0 then a.[2] round.[0] = 00
        BuildCtrlTakahashiAdder vs (us @ a) [fs.[0]] // clean up the ancilla a.[0] after each round
        // step 4: dispatch 4 mutually exclusive cases
        // case 4a: u even: conditioning on pattern 10 
        X [round.[0]]
        BuildCtrlBinaryDoubling ss [a.[2]; round.[0]; fs.[0]] [a.[0]] // using a as dirty ancilla
        BuildCtrlBinaryHalfing  us [a.[2]; round.[0]; fs.[0]] [a.[0]]
        X [round.[0]]
        // case 4b: v even: conditioning on pattern 01 
        X [a.[2]]
        BuildCtrlBinaryDoubling rs [a.[2]; round.[0]; fs.[0]] [a.[0]]  
        BuildCtrlBinaryHalfing  vs [a.[2]; round.[0]; fs.[0]] [a.[0]]
        X [a.[2]]
        // case 4c: u odd, v odd, u > v: conditioning on pattern 11
        BuildMultiCtrlTakahashiModAdderInverse us vs [a.[2]; round.[0]; fs.[0]] // replace (u,v) -> (u-v,v)
        BuildCtrlBinaryHalfing  us [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace (u-v) -> (u-v)/2
        BuildMultiCtrlTakahashiModAdder rs ss [a.[2]; round.[0]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling ss [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace s -> 2s        
        // case 4d: u odd, v odd, u <= v: conditioning on pattern 00
        X [a.[2]]
        X [round.[0]]
        BuildMultiCtrlTakahashiModAdderInverse vs us [a.[2]; round.[0]; fs.[0]] // replace (u,v) -> (u,v-u)
        BuildCtrlBinaryHalfing  vs [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace (v-u) -> (v-u)/2
        BuildMultiCtrlTakahashiModAdder ss rs [a.[2]; round.[0]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling rs [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace r -> 2r      
        X [a.[2]]
        X [round.[0]]
        // step 5: clean up the round bit round.[0] which can be shown is equal to rs.[0]
        CNOT [rs.[0]; a.[2]]
        // step 6: increase counter if required
        X fs
        BuildControlShiftRegisterCounterSpecial ks fs 
        X fs
    
    let n = us.Length
    for i in 0..(2*n-1) do 
        let round = [mg.[i]] // get 1 bits per round for match statement on (u, v)
        EuclidRound us vs ss rs round fs a ks 
   
let InverterSuperEfficientModularReduction (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let rs = slice qs [(3*n+1)..(4*n+1)] // n+1 bits for the non-reduced register holding r
    // check if final result in r overflows and if so, then subtract the modulus once
    BuildTakahashiAdderInverse (slice rs [0..(n-1)]) (ms @ [rs.[n]]) 
    BuildCtrlTakahashiModAdder (slice rs [0..(n-1)]) ms [rs.[n]]

let InverterSuperEfficientFixAlmostInverseStandardToMonty (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let ks = slice qs [8*n+4..8*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [8*n+7+counterSize..8*n+6+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [7*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [8*n+4+counterSize..8*n+6+counterSize] // 3 clean ancillas 
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

let InverterSuperEfficientFixAlmostInverseMontyToMonty (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let ks = slice qs [8*n+4..8*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [8*n+7+counterSize..8*n+6+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [7*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [8*n+4+counterSize..8*n+6+counterSize] // 3 clean ancillas 
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

let InverterSuperEfficientFixAlmostInverseStandardToStandard (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let ks = slice qs [8*n+4..8*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [8*n+7+counterSize..8*n+6+2*counterSize] // register for counter in phase 2
    let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [7*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [8*n+4+counterSize..8*n+6+counterSize] // 3 clean ancillas 
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
         
let InverterSuperEfficientPrepCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    let initCounter =  (pown 2I counterSize) - 1I
    
    let ss = slice qs [(2*n)..(3*n)] // prepare initial state on s register: n+1 bits
    let fs = slice qs [(7*n+2)..(7*n+3)] // prepare flags for counter 1 and counter 2: 2 bit
    let cs1 = slice qs [(8*n+4)..(8*n+counterSize+3)]
    let cs2 = slice qs [(8*n+counterSize+7)..(8*n+2*counterSize+6)]
    let us = slice qs [0..(n-1)] // n bits for the register holding the us
    let ms = slice qs [4*n+2..5*n+1] // n bits for the register holding the modulus
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (ms:Qubits) (us:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [ms.[i]; us.[i]]        
    CopyRegs ms us 
    X ss; X >< fs; X >< cs1; X >< cs2 // set first bit in ss and all remaining qubits in fs, cs1, and cs2 to 1
    
let InverterSuperEfficientCopyCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let rs = slice qs [(3*n+1)..(4*n)] // n bits for the reduced register holding r
    let res = slice qs [7*n+4..8*n+3] // n bits for the final result
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (rs:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [rs.[i]; res.[i]]        
    CopyRegs rs res
 
let MontgomerySuperEfficientInverseForward (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/8
    let us = slice qs [0..(n-1)] // u (n bits)
    let vs = slice qs [n..(2*n-1)] // v (n bits)
    let ss = slice qs [(2*n)..(3*n)] // s (n+1 bits)
    let rs = slice qs [(3*n+1)..(4*n+1)] // r (n+1 bits); holds the almost MG inverse
    let mg = slice qs [(5*n+2)..(7*n+1)] // round status indicators
    let fs = slice qs [7*n+2] // flag needed for 1st while loop
    let a  = slice qs [8*n+4+counterSize..8*n+6+counterSize] // need 3 ancillas
    let ks = slice qs [8*n+4..8*n+3+counterSize] // shift register counter 
    BuildSuperEfficientMontgomeryInverseForward us vs ss rs mg fs a ks
 
let MontgomerySuperEfficientInverse (name:string) (counterSize:int) (verbose:bool) (qs:Qubits) =     
    let a = Circuit.Compile (InverterSuperEfficientPrepCircuit counterSize) qs
    let b = Circuit.Compile (MontgomerySuperEfficientInverseForward counterSize) qs
    let c = Circuit.Compile (InverterSuperEfficientModularReduction counterSize) qs
    let d = match name with 
            | "MontgomerySM" -> Circuit.Compile (InverterSuperEfficientFixAlmostInverseStandardToMonty counterSize) qs
            | "MontgomeryMM" -> Circuit.Compile (InverterSuperEfficientFixAlmostInverseMontyToMonty counterSize) qs
            | "MontgomerySS" -> Circuit.Compile (InverterSuperEfficientFixAlmostInverseStandardToStandard counterSize) qs
            | _ -> failwith "MG inverse: this shoud never happen"             
    let e = Circuit.Compile (InverterSuperEfficientCopyCircuit counterSize) qs
    let f = d.Reverse() // Run the correction to fix almost inverse in reverse
    let g = c.Reverse() // Run the final modular reduction in reverse
    let h = b.Reverse() // Run the sequence of Montgomery steps in reverse
    let i = a.Reverse() // Run the sequence of state preparation of the us in reverse
    let gates = 
        [a; b; c; d; e; f; g; h; i] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         


///////////////////////////////////////////////////////////////////////////
// Ultra efficient MG inverse (7n qubits; modulus eliminated)
///////////////////////////////////////////////////////////////////////////

// BuildUltraEfficientMontgomeryInverseForward implements the inverse mod p of numbers that are given in Montgomery foom x * 2^n mod p. The output is 
// again in Montogomery form. The implementation follows closely the Kalinski paper which in turn is a binary extended Euclidean algorithm. 
// The implementation uses a total of 8n+6+2*ceil(log(n)) qubits. The registers are: 
// us = intially equal to modulus p, vs = intitially equal to input x, ss = intially 1, rs = intiially 0, ms = modulus p, mg = 8n qubits 
// holding the results of the dispatch predicates for each round, fs = flag that determines whether we are in compute or in counter mode, 
// ks = a register of qubits to implement a counter from 0 up to worst case iteration bound of 2n, and a = 1 ancilla. 
// The method uses the fact that the rs register can be cleaned up in each round, eliminating 2n qubits from the Efficient MG Inverse circuit. 
// This method still uses the extra register to store the modulus. 

let BuildUltraEfficientMontgomeryInverseForward (us:Qubits) (vs:Qubits) (ss:Qubits) (rs:Qubits) (mg:Qubits) (fs:Qubits) (a:Qubits) (ks:Qubits) = 
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
        CCNOT [fs.[0]; us.[0]; a.[2]] // if ue=1 then a.[2] round.[0] = 10
        X us 
        // step 3b: or else, check whether v is even
        X vs 
        X [a.[2]]
        BuildMultiplyControlledNOT [fs.[0]; a.[2]; vs.[0]] [round.[0]] [ss.[0]] // if ve=1 and ue = 0 then a.[2] round.[0] = 01  
        X [a.[2]]
        X vs
        // step 3c: or else, check whether u > v
        BuildCtrlTakahashiAdderInverse vs (us @ a) [fs.[0]] 
        CNOT [a.[2]; a.[1]]
        CNOT [round.[0]; a.[1]] // compute partity of a.[2] round.[0] into ancilla a.[1]
        X [a.[1]] // apply following gate only if ancilla is 0
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [a.[2]] [ss.[0]]
        BuildMultiplyControlledNOT [fs.[0]; a.[0]; a.[1]] [round.[0]] [ss.[0]] // if ue=ve=0 and ug=1 then a.[2] round.[0]= 11
        X [a.[1]] 
        CNOT [a.[2]; a.[1]]
        CNOT [round.[0]; a.[1]] // clean up ancilla a.[1] after each round. Note that parity is preseved in all cases        
        // step 3d: else, we have that u <= v
        // note: no computation needed in this case: ue=ve=ug=0 then a.[2] round.[0] = 00
        BuildCtrlTakahashiAdder vs (us @ a) [fs.[0]] // clean up the ancilla a.[0] after each round
        // step 4: dispatch 4 mutually exclusive cases
        // case 4a: u even: conditioning on pattern 10 
        X [round.[0]]
        BuildCtrlBinaryDoubling ss [a.[2]; round.[0]; fs.[0]] [a.[0]] // using a as dirty ancilla
        BuildCtrlBinaryHalfing  us [a.[2]; round.[0]; fs.[0]] [a.[0]]
        X [round.[0]]
        // case 4b: v even: conditioning on pattern 01 
        X [a.[2]]
        BuildCtrlBinaryDoubling rs [a.[2]; round.[0]; fs.[0]] [a.[0]]  
        BuildCtrlBinaryHalfing  vs [a.[2]; round.[0]; fs.[0]] [a.[0]]
        X [a.[2]]
        // case 4c: u odd, v odd, u > v: conditioning on pattern 11
        BuildMultiCtrlTakahashiModAdderInverse us vs [a.[2]; round.[0]; fs.[0]] // replace (u,v) -> (u-v,v)
        BuildCtrlBinaryHalfing  us [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace (u-v) -> (u-v)/2
        BuildMultiCtrlTakahashiModAdder rs ss [a.[2]; round.[0]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling ss [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace s -> 2s        
        // case 4d: u odd, v odd, u <= v: conditioning on pattern 00
        X [a.[2]]
        X [round.[0]]
        BuildMultiCtrlTakahashiModAdderInverse vs us [a.[2]; round.[0]; fs.[0]] // replace (u,v) -> (u,v-u)
        BuildCtrlBinaryHalfing  vs [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace (v-u) -> (v-u)/2
        BuildMultiCtrlTakahashiModAdder ss rs [a.[2]; round.[0]; fs.[0]] // replace r -> r+s
        BuildCtrlBinaryDoubling rs [a.[2]; round.[0]; fs.[0]] [a.[0]] // replace r -> 2r      
        X [a.[2]]
        X [round.[0]]
        // step 5: clean up the round bit round.[0] which can be shown is equal to rs.[0] + 1
        CNOT [rs.[0]; a.[2]]
        // step 6: increase counter if required
        X fs
        BuildControlShiftRegisterCounterSpecial ks fs 
        X fs
    
    let n = us.Length
    for i in 0..(2*n-1) do 
        let round = [mg.[i]] // get 1 bits per round for match statement on (u, v)
        EuclidRound us vs ss rs round fs a ks 
   
let InverterUltraEfficientModularReduction (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    //let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let rs = slice qs [(3*n+1)..(4*n+1)] // n+1 bits for the non-reduced register holding r
    let ts = slice qs [0..1] // need 2 dirty ancillas
    // check if final result in r overflows and if so, then subtract the modulus once
    
    //BuildTakahashiAdderInverse (slice rs [0..(n-1)]) (ms @ [rs.[n]]) 
    //BuildCtrlTakahashiModAdder (slice rs [0..(n-1)]) ms [rs.[n]]
    ConstantSubtractorInplaceCompile 0 ms (rs @ ts)
    ConstantAdderInplace ms (slice rs [0..(n-1)]) ts [rs.[n]] 

let InverterUltraEfficientFixAlmostInverseStandardToMonty (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let ks = slice qs [7*n+4..7*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [7*n+7+counterSize..7*n+6+2*counterSize] // register for counter in phase 2
    //let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [6*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [7*n+4+counterSize..7*n+6+counterSize] // 3 clean ancillas 
    let ts = slice qs [0..1] // need 2 dirty ancillas
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times
    for i in 0..(n-1) do 
        // mapping -x^-1 2^k mod p to -x^-1 2^k 2^-n mod p 
        //BuildModularHalver rs ms a         
        BuildModularHalverConstantPrime ms rs a ts
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -x^-1 2^k 2^-n mod p to -x^-1 2^k 2^-n 2^(2n-k) = -x^-1 2^n mod p 
    for i in 0..(n-1) do 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        //BuildCtrlModularDoubler rs ms fs a 
        BuildCtrlModularDoublerConstantPrime ms rs a ts fs
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1 2^n)=x^-1 2^n mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    BuildCtrlNegatorConstantPrime ms rs ts [] 
    
let InverterUltraEfficientFixAlmostInverseMontyToMonty (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let ks = slice qs [7*n+4..7*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [7*n+7+counterSize..7*n+6+2*counterSize] // register for counter in phase 2
    //let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [6*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [7*n+4+counterSize..7*n+6+counterSize] // 3 clean ancillas 
    let ts = slice qs [0..1] // need 2 dirty ancillas
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times.
    for i in 0..(n-1) do 
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -(x^-1 2^-n) 2^k mod p to -x^-1 2^k 2^-n 2^(2n-k) = -x^-1 2^n mod p 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        
        //BuildCtrlModularDoubler rs ms fs a 
        BuildCtrlModularDoublerConstantPrime ms rs a ts fs        
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1 2^n)=x^-1 2^n mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    BuildCtrlNegatorConstantPrime ms rs ts [] 
    
let InverterUltraEfficientFixAlmostInverseStandardToStandard (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let ks = slice qs [7*n+4..7*n+3+counterSize] // register holding the counter value of while loop in phase 1 
    let ls = slice qs [7*n+7+counterSize..7*n+6+2*counterSize] // register for counter in phase 2
    //let ms = slice qs [(4*n+2)..(5*n+1)] // n bits for the modulus
    let fs = slice qs [6*n+3] // flag needed for 2nd while loop
    let rs = slice qs [(3*n+1)..(4*n)] // the reduced part of the r register
    let a  = slice qs [7*n+4+counterSize..7*n+6+counterSize] // 3 clean ancillas 
    let ts = slice qs [0..1] // need 2 dirty ancillas
    // note: the ks register has the value 2n-k at this point as the loop terminated in k steps, i.e., 
    // the shift register state in ks was incremented exactly 2n-k times
    for i in 0..(2*n-1) do 
         // mapping -x^-1 2^k mod p to -x^-1 2^k 2^-2n mod p 
         //BuildModularHalver rs ms a         
         BuildModularHalverConstantPrime ms rs a ts
    for i in 0..(n-1) do 
        // count down from 2n-k and in each step double the value of the r register, 
        // mapping -x^-1 2^k 2^-2n mod p to -x^-1 2^k 2^-2n 2^(2n-k) = -x^-1 mod p 
        BuildMultiplyControlledNOT ks [fs.[0]] [a.[0]] // using a as dirty ancilla 
        X fs
        BuildMultiplyControlledNOT ls [fs.[0]] [a.[0]] // using a as dirty ancilla 
        BuildControlShiftRegisterCounterSpecialInverse ks fs
        //BuildCtrlModularDoubler rs ms fs a 
        BuildCtrlModularDoublerConstantPrime ms rs a ts fs        
        X fs
        BuildControlShiftRegisterCounterSpecial ls fs
        X fs
    // finally, compute p-(-x^-1)=x^-1 mod p 
    // use one's complement identity (a-b) = (a'+b)', where ' denotes one's complement
    BuildCtrlNegatorConstantPrime ms rs ts [] 
         
let InverterUltraEfficientPrepCircuit (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1
    let initCounter =  (pown 2I counterSize) - 1I
   
    let ss = slice qs [(2*n)..(3*n)] // prepare initial state on s register: n+1 bits
    let fs = slice qs [(6*n+2)..(6*n+3)] // prepare flags for counter 1 and counter 2: 2 bit
    let cs1 = slice qs [(7*n+4)..(7*n+counterSize+3)]
    let cs2 = slice qs [(7*n+counterSize+7)..(7*n+2*counterSize+6)]
    let us = slice qs [0..(n-1)] // n bits for the register holding the us    
    // copy out the qubits holding the result using a cascade of CNOTs
    let a = BoolInt n ms 
    let ApplyX (x:int) (q:Qubits) = 
        match x with 
        | 1 -> X q
        | _ -> ()        
    us |> List.iteri (fun i q -> ApplyX a.[i] [q]) 
    X ss; X >< fs; X >< cs1; X >< cs2 // set first bit in ss and all remaining qubits in fs, cs1, and cs2 to 1
    
let InverterUltraEfficientCopyCircuit (counterSize:int) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let rs = slice qs [(3*n+1)..(4*n)] // n bits for the reduced register holding r
    let res = slice qs [6*n+4..7*n+3] // n bits for the final result
    // copy out the qubits holding the result using a cascade of CNOTs
    let CopyRegs (rs:Qubits) (res:Qubits) = 
        for i in 0..(n-1) do 
           CNOT [rs.[i]; res.[i]]        
    CopyRegs rs res
 
let MontgomeryUltraEfficientInverseForward (counterSize:int) (ms:bigint) (qs:Qubits) = 
    let n = (qs.Length-2*counterSize-7)/7
    let us = slice qs [0..(n-1)] // u (n bits)
    let vs = slice qs [n..(2*n-1)] // v (n bits)
    let ss = slice qs [(2*n)..(3*n)] // s (n+1 bits)
    let rs = slice qs [(3*n+1)..(4*n+1)] // r (n+1 bits); holds the almost MG inverse
    let mg = slice qs [(4*n+2)..(6*n+1)] // round status indicators
    let fs = slice qs [6*n+2] // flag needed for 1st while loop
    let a  = slice qs [7*n+4+counterSize..7*n+6+counterSize] // need 3 ancillas
    let ks = slice qs [7*n+4..7*n+3+counterSize] // shift register counter 
    BuildUltraEfficientMontgomeryInverseForward us vs ss rs mg fs a ks
 
let MontgomeryUltraEfficientInverseManual (name:string) (counterSize:int) (verbose:bool) (ms:bigint) (qs:Qubits) =     
    let watch = Diagnostics.Stopwatch() 
    
    watch.Reset()
    watch.Start()   
    let a = Circuit.Compile (InverterUltraEfficientPrepCircuit counterSize ms) qs
    watch.Stop()
    if verbose then show "INVERTER: Time for step a: %f" watch.Elapsed.TotalSeconds

    watch.Reset()
    watch.Start()   
    let b = Circuit.Compile (MontgomeryUltraEfficientInverseForward counterSize ms) qs
    watch.Stop()
    if verbose then show "INVERTER: Time for step b: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let c = Circuit.Compile (InverterUltraEfficientModularReduction counterSize ms) qs
    watch.Stop()
    if verbose then show "INVERTER: Time for step c: %f" watch.Elapsed.TotalSeconds
    
    watch.Reset()
    watch.Start()   
    let d = match name with 
            | "MontgomerySM" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseStandardToMonty counterSize ms) qs
            | "MontgomeryMM" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseMontyToMonty counterSize ms) qs
            | "MontgomerySS" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseStandardToStandard counterSize ms) qs
            | _ -> failwith "MG inverse: this shoud never happen"             
    watch.Stop()
    if verbose then show "INVERTER: Time for step d: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let e = Circuit.Compile (InverterUltraEfficientCopyCircuit counterSize) qs
    watch.Stop()
    if verbose then show "INVERTER: Time for step e: %f" watch.Elapsed.TotalSeconds
    
    let f = d.Reverse() // Run the correction to fix almost inverse in reverse
    watch.Stop()
    if verbose then show "INVERTER: Time for step f: %f" watch.Elapsed.TotalSeconds
    
    let g = c.Reverse() // Run the final modular reduction in reverse
    watch.Stop()
    if verbose then show "INVERTER: Time for step g: %f" watch.Elapsed.TotalSeconds
    
    let h = b.Reverse() // Run the sequence of Montgomery steps in reverse
    watch.Stop()
    if verbose then show "INVERTER: Time for step h: %f" watch.Elapsed.TotalSeconds
    
    let i = a.Reverse() // Run the sequence of state preparation of the us in reverse
    watch.Stop()
    if verbose then show "INVERTER: Time for step i: %f" watch.Elapsed.TotalSeconds
        
    watch.Reset()
    watch.Start()   
    let gates = 
        [a; b; c; d; e; f; g; h; i] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    watch.Stop()
    if verbose then show "INVERTER: Time to concatenate: %f" watch.Elapsed.TotalSeconds
    // Compute the number of Toffoli gates in the circuit and write to console    
    //PrintToffoliNetworkMetrics gates
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
    
        
let MontgomeryUltraEfficientInverse (name:string) (counterSize:int) (verbose:bool) (ms:bigint) (qs:Qubits) =       
    let a = Circuit.Compile (InverterUltraEfficientPrepCircuit counterSize ms) qs
    let b = Circuit.Compile (MontgomeryUltraEfficientInverseForward counterSize ms) qs
    let c = Circuit.Compile (InverterUltraEfficientModularReduction counterSize ms) qs
    let d = match name with 
            | "MontgomerySM" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseStandardToMonty counterSize ms) qs
            | "MontgomeryMM" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseMontyToMonty counterSize ms) qs
            | "MontgomerySS" -> Circuit.Compile (InverterUltraEfficientFixAlmostInverseStandardToStandard counterSize ms) qs
            | _ -> failwith "MG inverse: this shoud never happen"             
    let e = Circuit.Compile (InverterUltraEfficientCopyCircuit counterSize) qs
    let f = d.Reverse() // Run the correction to fix almost inverse in reverse
    let g = c.Reverse() // Run the final modular reduction in reverse
    let h = b.Reverse() // Run the sequence of Montgomery steps in reverse
    let i = a.Reverse() // Run the sequence of state preparation of the us in reverse
    let gates = 
        [a; b; c; d; e; f; g; h; i] 
        |> List.map (fun x -> MyCircuitExport x) 
        |> List.concat
    // Compute the number of Toffoli gates in the circuit and write to console
    if verbose then 
        show "Number of Toffoli gates = %A " (gates |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    gates                         
      
///////////////////////////////////////////////////////////////////////////
// Functions to run modular inverters
///////////////////////////////////////////////////////////////////////////

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
            0I |> BoolInt n |> PrepBool qs 0 1 |> ignore           // u register: n bits
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


let RunSuperEfficientInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    // this function saves one of the two bits that we need to compute in each round. 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    
    // Prepare the initial state. The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | modulus | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | n       | 2n        | 2     | n      | log(n)   | 3        | log(n)   ||
    // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 000      | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let k = Ket(8*n+7+2*counterSize)
            let qs = k.Qubits
            0I |> BoolInt n |> PrepBool qs 0 1 |> ignore            // u register: n bits
            x  |> BoolInt n |> PrepBool qs n 1 |> ignore            // v register: n bits
            0I |> BoolInt (n+1) |> PrepBool qs (2*n) 1 |> ignore    // s register: n+1 bits
            0I |> BoolInt (n+1) |> PrepBool qs (3*n+1) 1 |> ignore  // r register: n+1 bits
            m  |> BoolInt n |> PrepBool qs (4*n+2) 1 |> ignore      // m register: n bits (the modulus)
            0I |> BoolInt (2*n) |> PrepBool qs (5*n+2) 1 |> ignore  // information to unravel the compressed rounds: 2*n bits
            0I |> BoolInt 2 |> PrepBool qs (7*n+2) 1 |> ignore     // flags for counter 1 and counter 2: 2 bit
            0I |> BoolInt n |> PrepBool qs (7*n+4) 1 |> ignore  // final result: n bits
            0I |> BoolInt counterSize |> PrepBool qs (8*n+4) 1 |> ignore // counter: ceil(log(n))+1 bits
            0I |> BoolInt 3 |> PrepBool qs (8*n+4+counterSize) 1 |> ignore  // 3 ancillas for comparator, inplace doubler, and MultiplyCNOTs, and round qubit
            0I |> BoolInt counterSize |> PrepBool qs (8*n+7+counterSize) 1 |> ignore // counter: ceil(log(n))+1 bits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier
    
    let arrangeInverterCircuit = 
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let INV = MontgomerySuperEfficientInverse name counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (8*n+7+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    watch.Reset()
    watch.Start()   
    let INV, rslt = arrangeInverterCircuit
    let MontgomerySuperEfficientInverterCircuit = INV qs 
    watch.Stop()
    show "Time to compile %f" watch.Elapsed.TotalSeconds
    
    if verbose1 then 
        show "Number of qubits = %A" (8*n+7+2*counterSize)

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let qs2 = arrangeInverterInputs 
    let initialState = Array.zeroCreate (8*n+7+2*counterSize)
    List.iter (fun i -> M !!(qs2,i)) [0..qs2.Length-1]
    for i in 0..(qs2.Length-1) do 
        initialState.[i] <- qs2.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MontgomerySuperEfficientInverterCircuit initialState
    if verbose2 then 
        show "The initial state is %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[(7*n+4)..(8*n+3)] 
    show "Final result = %A" res
    let resInt = (IntBool res.Length res)
    show "As a number mod p this is = %A" resInt
    resInt


let RunUltraEfficientInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    // this function saves one of the two bits that we need to compute in each round. 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    
    // Prepare the initial state. The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | 2n        | 2     | n      | log(n)   | 3        | log(n)   ||
    // init: || p  | x  | 1   | 0   | 0..0      | 11    | 0      | 1..1     | 000      | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let k = Ket(7*n+7+2*counterSize)
            let qs = k.Qubits
            0I |> BoolInt n |> PrepBool qs 0 1 |> ignore            // u register: n bits
            x  |> BoolInt n |> PrepBool qs n 1 |> ignore            // v register: n bits
            0I |> BoolInt (n+1) |> PrepBool qs (2*n) 1 |> ignore    // s register: n+1 bits
            0I |> BoolInt (n+1) |> PrepBool qs (3*n+1) 1 |> ignore  // r register: n+1 bits            
            0I |> BoolInt (2*n) |> PrepBool qs (4*n+2) 1 |> ignore  // information to unravel the compressed rounds: 2*n bits
            0I |> BoolInt 2 |> PrepBool qs (6*n+2) 1 |> ignore     // flags for counter 1 and counter 2: 2 bit
            0I |> BoolInt n |> PrepBool qs (6*n+4) 1 |> ignore  // final result: n bits
            0I |> BoolInt counterSize |> PrepBool qs (7*n+4) 1 |> ignore // counter: ceil(log(n))+1 bits
            0I |> BoolInt 3 |> PrepBool qs (7*n+4+counterSize) 1 |> ignore  // 3 ancillas for comparator, inplace doubler, and MultiplyCNOTs, and round qubit
            0I |> BoolInt counterSize |> PrepBool qs (7*n+7+counterSize) 1 |> ignore // counter: ceil(log(n))+1 bits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier
    
    let arrangeInverterCircuit = 
        match name with 
        | "MontgomerySM" | "MontgomerySS" | "MontgomeryMM" -> 
            let INV = MontgomeryUltraEfficientInverse name counterSize verbose1 m // constructing the multiplier circuit
            let rslt = Array.create (7*n+7+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    watch.Reset()
    watch.Start()   
    let INV, rslt = arrangeInverterCircuit
    let MontgomeryUltraEfficientInverterCircuit = INV qs 
    watch.Stop()
    show "Time to compile %f" watch.Elapsed.TotalSeconds
    
    if verbose1 then 
        show "Number of qubits = %A" (7*n+7+2*counterSize)

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let qs2 = arrangeInverterInputs 
    let initialState = Array.zeroCreate (7*n+7+2*counterSize)
    List.iter (fun i -> M !!(qs2,i)) [0..qs2.Length-1]
    for i in 0..(qs2.Length-1) do 
        initialState.[i] <- qs2.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast MontgomeryUltraEfficientInverterCircuit initialState
    if verbose2 then 
        show "The initial state is %A" initialState
        show "And the final state is %A" finalState

    let res = finalState.[(6*n+4)..(7*n+3)] 
    show "Final result = %A" res
    let resInt = (IntBool res.Length res)
    show "As a number mod p this is = %A" resInt
    resInt


///////////////////////////////////////////////////////////////////////////
// Functions to compile modular inverters
///////////////////////////////////////////////////////////////////////////

let CompileEfficientModularInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    //let initCounter =  (pown 2I counterSize) - 1I

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
            let INV = MontgomeryEfficientInverse name counterSize verbose1 // constructing the multiplier circuit
            //let INV = MontgomeryEfficientInverseManualReverse name counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (10*n+6+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    let INV, rslt = arrangeInverterCircuit
    let MontgomeryEfficientInverterCircuit = INV qs 
    
    MontgomeryEfficientInverterCircuit


let CompileUltraEfficientModularInverter (name:string) (n:int) (m:bigint) (x:bigint) (verbosity:string) = 
    let verbose1 = (verbosity = "verbose") || (verbosity = "all") // level 1 = print circuit size and width
    let verbose2 = (verbosity = "all")                            // level 2 = print input and output states 
    let watch = Diagnostics.Stopwatch()

    let counterSize =  double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 
    //let initCounter =  (pown 2I counterSize) - 1I

    // Prepare the initial state. The partitioning of the quantum register is as follows: 
    // name: || us | vs | ss  | rs  | mg-rounds | flags | result | counter1 | ancillas | counter2 || 
    // bits: || n  | n  | n+1 | n+1 | 2n        | 2     | n      | log(n)   | 3        | log(n)   ||
    // init: || p  | x  | 1   | 0   | 0..0      | 11    | 0      | 1..1     | 000      | 1..1     || 

    let arrangeInverterInputs =   
        match name with 
        | "MontgomerySM" | "MontgomeryMM" | "MontgomerySS" -> 
            let k = Ket(7*n+7+2*counterSize)
            let qs = k.Qubits
            qs
        | _ -> failwith "Unknown inverter."
                      
    let qs = arrangeInverterInputs          // arranges inputs in pattern needed for the multiplier

    if verbose2 then 
        show "Number of qubits = %A" (7*n+7+2*counterSize)
    
    show "bitsize=%A" n

    let arrangeInverterCircuit = 
        match name with 
        | "MontgomerySM" | "MontgomeryMM" | "MontgomerySS" -> 
            let INV = MontgomeryUltraEfficientInverseManual name counterSize verbose1 m // constructing the multiplier circuit
            //let INV = MontgomeryEfficientInverseManualReverse name counterSize verbose1 // constructing the multiplier circuit
            let rslt = Array.create (7*n+7+2*counterSize) 0 // will be used to store final result 
            INV, rslt
        | _ -> failwith "This inverter is not implemented yet."

    let INV, rslt = arrangeInverterCircuit
    let MontgomeryUltraEfficientInverterCircuit = INV qs 
    
    MontgomeryUltraEfficientInverterCircuit


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
            if verbose2 then show "CSV: Time to construct MUL MM %f" watch.Elapsed.TotalSeconds            
            MULmon 
        | "DoubleAdd"  ->
            let MULdadd = CompileModularMultiplier "CtrlDoubleAdd" n 1I 1I N verbosity
            if verbose2 then printf "MUL DA done. "
            watch.Stop()
            if verbose2 then show "CSV: Time to construct MUL DA %f" watch.Elapsed.TotalSeconds               
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

    let watch = Diagnostics.Stopwatch()        
    
    watch.Start()
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
    watch.Stop()
    let t1 = watch.Elapsed.TotalSeconds
            
    let totalCirc = [ totalCirc1; totalCirc2; totalCirc3 ] |> Seq 
    let totalToffCirc = MyCircuitExport totalCirc
    
    let initialState = Array.zeroCreate (2*n+1+1)
    watch.Reset()
    watch.Start()
    let finalState = MyCircuitSimulateFast totalToffCirc initialState
    watch.Stop()
    let t2 = watch.Elapsed.TotalSeconds

    let qubits, size, depth = ReportToffoliNetworkMetrics totalToffCirc
    n, qubits, size, depth, t1, t2    
           
    //PrintToffoliNetworkMetrics totalToffCirc 
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


let ConstantMultiplierModNRunWholeCircuitFast (ccount:int) (a:bigint) (ainv:bigint) (n:int) (N:bigint) (qs:Qubits) = 
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]
    let one = bigint 1

    let watch = Diagnostics.Stopwatch()        
    watch.Start()
    
    let totalCirc1 = 
        [ for i in 0..(n-1) do
             printf "%i " i
             let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
             yield ((Circuit.Compile (ConstantAdderInplaceModN (1+ctrls.Length) ((a<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]])) |> MyCircuitExport)
        ] 
        |> List.concat
    
    let totalCirc2 = 
        (Circuit.Compile (CSWAPCompile (ctrls.Length)) (List.concat [x; y; ctrls])) |> MyCircuitExport
    
    let totalCirc3 = 
        [ for i in 0..(n-1) do
            let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
            yield ((Circuit.Compile (ConstantSubtractorInplaceModN (1+ctrls.Length) ((ainv<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]])) |> MyCircuitExport)
        ]
        |> List.concat
              
//    let totalCirc3 = 
//        List.rev totalCirc1
//    
    let totalToffCirc = 
        totalCirc1 @ totalCirc2 @ totalCirc3
        
    
    watch.Stop()
    let t1 = watch.Elapsed.TotalSeconds
                        
    let initialState = Array.zeroCreate (2*n+1+1)
    watch.Reset()
    watch.Start()
    let finalState = MyCircuitSimulateFast totalToffCirc initialState
    watch.Stop()
    let t2 = watch.Elapsed.TotalSeconds

    let qubits, size, depth = ReportToffoliNetworkMetrics totalToffCirc
    n, qubits, size, depth, t1, t2    

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
        if verbose2 then show "CSV: Time to remap construct INV: %f" watch.Elapsed.TotalSeconds            
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


///////////////////////////////////////////////////////////////////////////
// Shor's algorithm: dlog over ECC curves (constant prime)
///////////////////////////////////////////////////////////////////////////

let ec_add_affine_constant_prime (n:int) (p:bigint) (x2:bigint) (y2:bigint) verbosity = 
    let verbose1 = (verbosity = "verbose" || verbosity = "all")
    let verbose2 = (verbosity = "all")
    let m = double n |> log |> fun x -> x/(log 2.0) |> ceil |> int |> fun x -> x + 1 // counter size
    
    // Construction of the point addition formula for ECC in affine Weierstrass form, 
    // first point P=(X1,Y1), second point Q=(x2,y2), where Q is classically known 
    // and P is overwritten by P+Q. The entire point addition is controlled by a bit <ctrl>. 
    let qubits = [|0..(9*n+2*m+7)|] // total #qubits = 9n+2m+8
        // name: || ctrl | X1 | Y1 | t0 | lam | ancillas || 
        // bits: || 1    | n  | n  | n  | n   | 5n+2m+7 ||
        // maximum ancillas required: n + (n+1) + (n+1) + 2n + 2 + log(n) + 3 + log(n) = 5n + 2m + 7
    let ctrl = qubits.[0..0]    
    let X1   = qubits.[1..n]
    let Y1   = qubits.[(n+1)..(2*n)]
    let t0   = qubits.[(2*n+1)..(3*n)]
    let lam  = qubits.[(3*n+1)..(4*n)]
    let anc  = qubits.[(4*n+1)..(9*n+2*m+7)] 

    // Contruct all the arithmetic functions needed to implement the ECC point additions reversibly. 
    // First the circuits are constructed as Toffoli networks operating on qubits in a range [0..(R-1)], 
    // where R is the largest qubit index used by the gates. Then these circuits are embedded using a map. 
    
    let watch = Diagnostics.Stopwatch()        
    
    if verbose2 then printf "Creating circuits: "
    
    watch.Start()
    let NEG = CompileModularAdder "CtrlNegConstantPrime" n 1I 1I p 1I verbosity
    if verbose2 then printf "cNEG done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct cNEG %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let ADD = CompileModularAdder "ModADDConstantPrime" n 1I 1I p 1I verbosity
    if verbose2 then printf "ADD done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct ADD %f" watch.Elapsed.TotalSeconds            
    watch.Reset()
       
    watch.Start()
    let cADD = CompileModularAdder "CtrlModADDConstantPrime" n 1I 1I p 1I verbosity
    if verbose2 then printf "cADD done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct cADD %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let ADDc = CompileModularAdder "ModADDConstantPrimeConstantNumber" n 1I x2 p 1I verbosity
    if verbose2 then printf "ADDc done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct ADDc %f" watch.Elapsed.TotalSeconds            
    watch.Reset()
        
    watch.Start()
    let cADDc = CompileModularAdder "CtrlModADDConstantPrimeConstantNumber" n 1I x2 p 1I verbosity
    if verbose2 then printf "cADDc done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct cADDc %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let cADDc2 = CompileModularAdder "CtrlModADDConstantPrimeConstantNumber" n 1I y2 p 1I verbosity // need this for subtractor cSUBc with y2
    if verbose2 then printf "ADDc done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct ADDc %f" watch.Elapsed.TotalSeconds            
    watch.Reset()
    
    watch.Start()
    let SUB  = List.rev ADD
    if verbose2 then printf "SUB done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct SUB %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let cSUB = List.rev cADD
    if verbose2 then printf "cSUB done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct cSUB %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let SUBc  = List.rev ADDc
    if verbose2 then printf "SUBc done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct SUBc %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let cSUBc  = List.rev cADDc2
    if verbose2 then printf "cSUBc done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct cSUBc %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let MUL = CompileModularMultiplier "MontgomeryConstantPrime" n 1I 1I p verbosity
    if verbose2 then printf "MUL MM done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct MUL MM %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let SQU = CompileModularMultiplier "MontgomerySquarerConstantPrime" n 1I 1I p verbosity
    if verbose2 then printf "SQU done. "
    watch.Stop()
    if verbose2 then show "CSV: Time to construct SQU: %f" watch.Elapsed.TotalSeconds            
    watch.Reset()

    watch.Start()
    let INV = CompileUltraEfficientModularInverter "MontgomeryMM" n p 1I verbosity
    if verbose2 then printf "INV done.\n"
    watch.Stop()
    show "CSV: (bitsize=%A); Time to construct INV: %f" n watch.Elapsed.TotalSeconds            
    watch.Reset()

    // Mapping step: in the following register allotments, the following conventions are used: 
    // "-" stands for an input that is returned unchanged, "*" stands for an input that is overwritten and 
    // carries the output all other bits and registers are considered ancillas.

    let ctrl_neg_modp_prime_const (xs:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on n+3 qubits
        // old: 
        // name: ||*xs*|-modulus-|-ctrl-||
        // bits: || n  | n       | 1    ||
        // init: || x  | p       | c    ||   
        // new: 
        // name: ||*xs*|-ctrl-| anc ||
        // bits: || n  | 1    | 2   ||
        // init: || x  | c    | 0   ||
        let ds = [(n+1)..(n+2)] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[n] <- cs.[0]
        RemapToffoliCircuit NEG embed

    let add_modp_prime_const (xs:int []) (s:bigint) = // result is stored in xs
        // circuit operates nontrivially on 2n qubits
        // old: 
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp ||
        // bits: || n  | n  | 1   | n       | 1   ||
        // init: || x  | y  | 0   | p       | 0   ||
        // new: 
        // name: ||-xs-| tmp | gs  ||
        // bits: || n  | 1   | n-1 ||
        // init: || x  | 00  | g   ||
        let ds = concat [ [n]; [(n+1)..(2*n-1)] ] // pick up all ancillas       
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do             
            embed.[i] <- xs.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit ADDc embed

    let sub_modp_prime_const (xs:int []) (s:bigint) = // result is stored in xs
        // circuit operates nontrivially on 2n qubits
        // old: 
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp ||
        // bits: || n  | n  | 1   | n       | 1   ||
        // init: || x  | y  | 0   | p       | 0   ||
        // new: 
        // name: ||-xs-| tmp | gs  ||
        // bits: || n  | 1   | n-1 ||
        // init: || x  | 0   | g   ||
        let ds = concat [ [n]; [(n+1)..(2*n-1)] ] // pick up all ancillas       
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do             
            embed.[i] <- xs.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit SUBc embed
    
    let ctrl_add_modp_prime (xs:int []) (ys:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 2n+4 qubits
        // old: 
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        // new: 
        // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
        // name: ||-xs-|*ys*| tmp ||-ctrl-||
        // bits: || n  | n  | 111 || 1    ||
        // init: || x  | y  | 000 || c    ||
        let ds = concat [ [2*n..(2*n+2)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]            
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[2*n+3] <- cs.[0]
        RemapToffoliCircuit cADD embed
    
    let ctrl_add_modp_prime_const (xs:int []) (s:bigint) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 2n+1 qubits
        // old:   
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        // new: 
        // name: ||-xs-| tmp | gs  |-cs-||
        // bits: || n  | 1   | n-1 | 1  ||
        // init: || x  | 0   | g   | c  ||
        let ds = concat [ [n]; [(n+1)..(2*n-1)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]            
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[2*n] <- cs.[0]
        RemapToffoliCircuit cADDc embed
                
    let ctrl_sub_modp_prime_const (xs:int []) (s:bigint) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 2n+1 qubits
        // old:   
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        // new: 
        // name: ||-xs-| tmp | gs  |-cs-||
        // bits: || n  | 1   | n-1 | 1  ||
        // init: || x  | 0   | g   | c  ||
        let ds = concat [ [n]; [(n+1)..(2*n-1)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]            
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[2*n] <- cs.[0]
        RemapToffoliCircuit cSUBc embed

    let ctrl_sub_modp_prime (xs:int []) (ys:int []) (cs:int []) = // result is stored in xs
        // circuit operates nontrivially on 2n+4 qubits
        // old: 
        // name: ||-xs-|*ys*| tmp |-modulus-| tmp |-ctrl-||
        // bits: || n  | n  | 1   | n       | 2   | 1    ||
        // init: || x  | y  | 0   | p       | 00  | c    ||
        // new: 
        // Prepare the initial state for controlled modular adder. The partitioning of the quantum register is as follows: 
        // name: ||-xs-|*ys*| tmp ||-ctrl-||
        // bits: || n  | n  | 111 || 1    ||
        // init: || x  | y  | 000 || c    ||
        let ds = concat [ [2*n..(2*n+2)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- ys.[i]
            embed.[n+i] <- xs.[i]            
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        embed.[2*n+3] <- cs.[0]
        RemapToffoliCircuit cSUB embed

    let mul_modp_prime (xs:int []) (ys:int []) (res:int []) = // result is stored in res
        // circuit operates nontrivially on 5n+4 qubits
        // old:
        // name: ||-xs-|-ys-| ancilla | acc  |-modulus-| ancilla | mg-rounds |*result*|| 
        // bits: || n  | n  | 1       | n+2  | n       | 1       | n         | n      ||
        // init: || x  | y  | 0       | 0    | p       | 0       | 0..0      | 0      || 
        // new: 
        // Prepare the initial state for Montgomery. The partitioning of the quantum register is as follows: 
        // name: ||-xs-|-ys-| ancilla | acc  | ancilla | mg-rounds |*result*|| 
        // bits: || n  | n  | 1       | n+2  | 1       | n         | n      ||
        // init: || x  | y  | 0       | 0    | 0       | 0..0      | 0      || 
        
        let ds = concat [ [2*n]; [2*n..(3*n+2)]; [3*n+3]; [3*n+4..4*n+3] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i]
            embed.[n+i] <- ys.[i]
            embed.[4*n+4+i] <- res.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit MUL embed

    let squ_modp_prime (xs:int []) (res:int []) = // result is stored in res
        // circuit operates nontrivially on 4n+5 qubits
        // old:
        // name: ||-xs-| ancilla | acc  |-modulus-| ancilla | mg-rounds |*result*|| 
        // bits: || n  | 2       | n+2  | n       | 1       | n         | n      ||
        // init: || x  | 0       | 0    | p       | 0       | 0..0      | 0      ||         
        // new: 
        // name: ||-xs-| ancilla | acc  | ancilla | mg-rounds |*result*|| 
        // bits: || n  | 2       | n+2  | 1       | n         | n      ||
        // init: || x  | 0       | 0    | 0       | 0..0      | 0      || 
    
        let ds = concat [ [n..(n+1)]; [(n+2)..(2*n+3)]; [2*n+4]; [(2*n+5)..(3*n+4)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[i] <- xs.[i] 
            embed.[3*n+5+i] <- res.[i]             
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
        RemapToffoliCircuit SQU embed
        
    let inv_modp_prime (xs:int []) (res:int []) =  // result is stored in res
         // circuit operates nontrivially on 7n+2m+7 qubits
        // Prepare the initial state. The partitioning of the quantum register is as follows: 
        // old: 
        // name: || us |-vs-| ss  | rs  |-modulus-| mg-rounds | flags |*result*| counter1 | ancillas | counter2 || 
        // bits: || n  | n  | n+1 | n+1 | n       | 4n        | 2     | n      | log(n)   | 2        | log(n)   ||
        // init: || p  | x  | 1   | 0   | p       | 0..0      | 11    | 0      | 1..1     | 00       | 1..1     || 
        //        let ds = concat [ [0..(n-1)]; [(2*n)..(3*n)]; [(3*n+1)..(4*n+1)]; [(5*n+2)..(9*n+1)]; [(9*n+2)..(9*n+3)]; 
        //                          [(10*n+4)..(10*n+m+3)]; [(10*n+m+4)..(10*n+m+5)]; [(10*n+m+6)..(10*n+2*m+5)] ] // pick up all ancillas
        // new: 
        // name: || us |-vs-| ss  | rs  | mg-rounds | flags |*result*| counter1 | ancillas | counter2 || 
        // bits: || n  | n  | n+1 | n+1 | 2n        | 2     | n      | log(n)   | 3        | log(n)   ||
        // init: || p  | x  | 1   | 0   | 0..0      | 11    | 0      | 1..1     | 000      | 1..1     || 
        let ds = concat [ [0..(n-1)]; [(2*n)..(3*n)]; [(3*n+1)..(4*n+1)]; [(4*n+2)..(6*n+1)]; [(6*n+2)..(6*n+3)]; 
                          [(7*n+4)..(7*n+m+3)]; [(7*n+m+4)..(7*n+m+6)]; [(7*n+m+7)..(7*n+2*m+6)] ] // pick up all ancillas
        let embed = Dictionary<int,int>()
        for i in [0..(n-1)] do 
            embed.[n+i] <- xs.[i]            
            embed.[6*n+4+i] <- res.[i]
        for i in [0..(ds.Length-1)] do 
            embed.[ds.[i]] <- anc.[i]
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
        if verbose2 then show "CSV: Time to remap construct INV: %f" watch.Elapsed.TotalSeconds            
        watch.Reset()
        a

    let opList = [
        sub_modp_prime_const X1 x2;              // x1 <- x1 - x2
        ctrl_sub_modp_prime_const Y1 y2 ctrl;    // y1 <- y1 - y2    // CTRL 
        inv_modp_prime X1 t0;                    // t0 <- 1/(x1 - x2) 
        mul_modp_prime Y1 t0 lam;                // lam <- (y1 - y2)/(x1 - x2)
        mul_modp_prime lam X1 Y1;                // y1 <- 0
        inv_modp_prime X1 t0;                    // t0 <- 0 
        squ_modp_prime lam t0;                   // t0 <- ((y1 - y2)/(x1 - x2))^2 = lam^2
        ctrl_sub_modp_prime X1 t0 ctrl;          // x1 <- x1 - x2 - ((y1 - y2)/(x1 - x2))^2 // CTRL
        ctrl_add_modp_prime_const X1 x2 ctrl;    // CTRL
        ctrl_add_modp_prime_const X1 x2 ctrl;    // CTRL
        ctrl_add_modp_prime_const X1 x2 ctrl;    // x1 <- x1 - x2 - ((y1 - y2)/(x1 - x2))^2 + 3*x2 = -(x3 - x2) // CTRL
        squ_modp_prime lam t0;                   // t0 <- 0
        mul_modp_prime lam X1 Y1;                // y1 <- y3 + y2
        inv_modp_prime X1 t0;                    // t0 <- -1/(x3 - x2)
        mul_modp_prime t0 Y1 lam;                // lam <- 0
        inv_modp_prime X1 t0;                    // t0 <- 0        
        ctrl_neg_modp_prime_const X1 ctrl;       // x1 <- x3 - x2 // CTRL
        ctrl_sub_modp_prime_const Y1 y2 ctrl;    // y3 done // CTRL
        add_modp_prime_const X1 x2;              // x3 done
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
                printf "X1 = %A\n" X1
                printf "Y1 = %A\n" Y1
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
        | "AffineWeierstrassConstantPrime" -> 
            let init = Array.zeroCreate (9*n+2*m+8) // total #qubits = 9n+2m+8
            // name: || ctrl | X1 | Y1 | t0 | lam | ancillas || 
            // bits: || 1    | n  | n  | n  | n   | 4n+2m+7  ||
            // maximum ancillas required: n + (n+1) + (n+1) + 2n + 2 + log(n) + 3 + log(n) = 5n + 2m + 7
            match P, Q with 
            | AffinePoint(x1,y1), AffinePoint(x2,y2) ->
                let cs = BoolInt n ctrl;
                let X1 = BoolInt n x1; 
                let Y1 = BoolInt n y1;
                init.[0] <- cs.[0]
                for i in [0..(n-1)] do 
                    init.[i+1]   <- X1.[i]
                    init.[n+i+1] <- Y1.[i]                
            | _ -> failwith "wrong ECC point format; need P, Q both to be affine." 
            init
        | _ -> failwith "unknown ECC type."
   
    let ecc_circuit = // constructing the multiplier circuit
        match P, Q with 
            | AffinePoint(x1,y1), AffinePoint(x2,y2) ->
                match name with 
                | "AffineWeierstrass" -> ec_add_affine n p verbosity
                | "AffineWeierstrassConstantPrime" -> ec_add_affine_constant_prime n p x2 y2 verbosity
                | _ -> failwith "unknown ECC type"
            | _ -> failwith "unknown ECC point type"
        
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
    let mutable ecc_finalstate = 
        match name with 
        | "AffineWeierstrass" -> Array.zeroCreate (14*n+2*m+7) // total #qubits = 14n+2m+7
        | "AffineWeierstrassConstantPrime" -> Array.zeroCreate (9*n+2*m+8) // total #qubits = 8n+2m+8
        | _ -> failwith "Unknown ECC point addition method."
    for op in ecc_circuit do
        ecc_finalstate <-  MyCircuitSimulateFast op ecc_initialstate
        ecc_initialstate <- ecc_finalstate
    watch.Stop()
    show "CSV: Time to simulate with fast simulator %f" watch.Elapsed.TotalSeconds            
    let sim = watch.Elapsed.TotalSeconds
    
    let qubits, size, depth = ecc_circuit |> List.concat |> ReportToffoliNetworkMetrics

    if verbose2 && (name = "AffineWeierstrassConstantPrime") then 
        show "The initial state is   %A" ecc_initialstate
        show "And the final state is %A" ecc_finalstate                
        let x1dump = ecc_finalstate.[1..n] |> IntBool n
        let y1dump = ecc_finalstate.[(n+1)..(2*n)] |> IntBool n
        let t0dump = ecc_finalstate.[(2*n+1)..(3*n)] |> IntBool n
        let lamdump = ecc_finalstate.[(3*n+1)..(4*n)] |> IntBool n
        let ancdump = ecc_finalstate.[(4*n+1)..(9*n+2*m+7)] |> IntBool n
        show "x1=%A y1=%A t0=%A lam=%A anc=%A" x1dump y1dump t0dump lamdump ancdump

    if verbose2 && (name = "AffineWeierstrass") then 
        show "The initial state is   %A" ecc_initialstate
        show "And the final state is %A" ecc_finalstate                
        let x1dump = ecc_finalstate.[(n+1)..(2*n)] |> IntBool n
        let y1dump = ecc_finalstate.[(2*n+1)..(3*n)] |> IntBool n
        let t0dump = ecc_finalstate.[(5*n+1)..(6*n)] |> IntBool n
        let lamdump = ecc_finalstate.[(6*n+1)..(7*n)] |> IntBool n
        let ancdump = ecc_finalstate.[(7*n+1)..(14*n+2*m+6)] |> IntBool n
        show "x1=%A y1=%A t0=%A lam=%A anc=%A" x1dump y1dump t0dump lamdump ancdump

    let res = 
        match name with
        | "AffineWeierstrass" -> 
            AffinePoint( (IntBool n ecc_finalstate.[(n+1)..(2*n)]), (IntBool n ecc_finalstate.[(2*n+1)..(3*n)]))        
        | "AffineWeierstrassConstantPrime" -> 
            AffinePoint( (IntBool n ecc_finalstate.[1..n]), (IntBool n ecc_finalstate.[(n+1)..(2*n)]))
        | _ -> failwith "unknown ECC type."
    res, qubits, size, depth, sim


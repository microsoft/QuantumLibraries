module Microsoft.Research.Liquid.Integer

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open CircBase           // basic quantum circuits


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



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inplace work starts here
// Includes:
//  * Addition of constants in 10nlog(n) using 1 and 8nlog(n) using 2 dirty ancilla qubits
//  * Incrementer (by Gidney) using n dirty ancillae in ~4n Toffolis
//  * Carry computation using n dirty ancillae in ~4n Toffolis
//  * Addition of a constant modulo N in ~twice the cost of addition of constant
//  * Modular multiplication using repeated modular addition and shift including uncomputation, i.e.
//          |x>|0> --> |(ax)mod N>|0>
//      in 32n^2 log(n)
//
//  (+ various functions for testing)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Uses qubits xg of the form [x_0,...,x_n-1,garbage,...,garbage] starting at the LSB
// to add the constant c to x.
// This is an n^2 algorithm, see ConstantAdderInplace for an nlog(n) version
// using only 1 or 2 borrowed dirty qubits.
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

// n-bit incrementer using n borrowed bits
let BorrowedQubitsIncrementer (xg:Qubits) = 
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

// Controlled n-bit incrementer using n borrowed bits
// Last qubit is the control
let CtrlBorrowedQubitsIncrementer (xgc:Qubits) = 
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


// Controlled SWAP operation used by the constant modular multiplier
let CSWAP (x:Qubits) (y:Qubits) (ctrls:Qubits) =
    let n = x.Length
    // perform swap between x and y conditioned on ctrl
    for i in 0..(n-1) do
        BuildMultiplyControlledNOT [x.[i]] [y.[i]] [x.[(i+1)%n]]
        BuildMultiplyControlledNOT (List.concat [[y.[i]];ctrls]) [x.[i]] [y.[(i+1)%n]]
        BuildMultiplyControlledNOT [x.[i]] [y.[i]] [x.[(i+1)%n]]

let CSWAPCompile (numctrls:int) (xyc:Qubits) =
    let n = int ((xyc.Length-numctrls)/2)
    let x = slice xyc [0..(n-1)]
    let y = slice xyc [n..(2*n-1)]
    let c = slice xyc [(2*n)..(2*n+numctrls-1)]
    CSWAP x y c

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



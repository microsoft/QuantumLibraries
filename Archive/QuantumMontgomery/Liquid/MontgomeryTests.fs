module  Microsoft.Research.Liquid.MontgomeryTests

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util
open Operations

open Montgomery


///////////////////////////////////////////////////////////////////////////
//
// Tests: Some useful tools 
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
let __TestCircuitDepth () = 
    let gates = [([1;2;3],1); ([1;2],0); ([2;3;5],3); ([6;10;99;4],1); ([5..97],0); ([1;6;2],1); ([3;5;32],3)]
    printf "the depth is %A\n" (CircuitDepthStreamed gates)
    
[<LQD>]
let __TestTdepth () = 
    let k = Ket(10)
    let qs = k.Qubits

    let toffLiq1 (qs:Qubits) = 
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        X [qs.[3]]
        CCNOT [qs.[2]; qs.[5]; qs.[6]]
        CNOT [qs.[2]; qs.[4]]
        CCNOT [qs.[5]; qs.[4]; qs.[6]]
    
    let toffLiq2 (qs:Qubits) = 
        CNOT [qs.[1]; qs.[2]] 
        CNOT [qs.[2]; qs.[3]] 
        CNOT [qs.[3]; qs.[4]] 
        CNOT [qs.[4]; qs.[5]] 
        CCNOT [qs.[4]; qs.[5]; qs.[6]]
        CCNOT [qs.[7]; qs.[8]; qs.[9]]

    let toffLiq3 (qs:Qubits) = 
        CCNOT [qs.[1]; qs.[2]; qs.[5]]
        X [qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[4]]
        CNOT [qs.[1]; qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        CCNOT [qs.[5]; qs.[6]; qs.[7]]

    let toffLiq4 (qs:Qubits) = 
        // 1
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        X [qs.[3]]
        CCNOT [qs.[2]; qs.[5]; qs.[6]]
        CNOT [qs.[2]; qs.[4]]
        CCNOT [qs.[5]; qs.[4]; qs.[6]]
        // 3
        CCNOT [qs.[1]; qs.[2]; qs.[5]]
        X [qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[4]]
        CNOT [qs.[1]; qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        CCNOT [qs.[5]; qs.[6]; qs.[7]]
        // 2
        CNOT [qs.[1]; qs.[2]] 
        CNOT [qs.[2]; qs.[3]] 
        CNOT [qs.[3]; qs.[4]] 
        CNOT [qs.[4]; qs.[5]] 
        CCNOT [qs.[4]; qs.[5]; qs.[6]]
        CCNOT [qs.[7]; qs.[8]; qs.[9]]
        // 1
        CCNOT [qs.[1]; qs.[2]; qs.[5]]
        X [qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[4]]
        CNOT [qs.[1]; qs.[2]]
        CCNOT [qs.[1]; qs.[2]; qs.[3]]
        CCNOT [qs.[5]; qs.[6]; qs.[7]]

    let toffLiqCirc1 = Circuit.Compile toffLiq1 qs    
    let toffLiqCirc2 = Circuit.Compile toffLiq2 qs    
    let toffLiqCirc3 = Circuit.Compile toffLiq3 qs    
    let toffLiqCirc4 = Circuit.Compile toffLiq4 qs    

    let toffLiqCircFold1 = toffLiqCirc1.Fold(true)
    let toffLiqCircFold2 = toffLiqCirc2.Fold(true)
    let toffLiqCircFold3 = toffLiqCirc3.Fold(true)
    let toffLiqCircFold4 = toffLiqCirc4.Fold(true)

    let toffCirc1 = [ MyTOFF(1,2,3); MyNOT(3); MyTOFF(2,5,6); MyCNOT(2,4); MyTOFF(5,4,6) ]
    let toffCirc2 = [ MyCNOT(1,2); MyCNOT(2,3); MyCNOT(3,4); MyCNOT(4,5); MyTOFF(4,5,6); MyTOFF(7,8,9) ]
    let toffCirc3 = [ MyTOFF(1,2,5); MyNOT(2); MyTOFF(1,2,4); MyCNOT(1,2); MyTOFF(1,2,3); MyTOFF(5,6,7) ]
    let toffCirc4 = toffCirc1 @ toffCirc3 @ toffCirc2 @ toffCirc1 
    printfn "doing 1"
    let a1 = toffLiqCirc1.GateCount(true, fun x -> x.Name = "CCNOT")
    let b1 = toffLiqCircFold1.GateCount(true, fun x -> x.Name = "CCNOT")
    let c1 = ToffoliDepthStreamed toffCirc1
    
    printfn "\n"
    printfn "doing 2"
    let a2 = toffLiqCirc2.GateCount(true, fun x -> x.Name = "CCNOT")
    let b2 = toffLiqCircFold2.GateCount(true, fun x -> x.Name = "CCNOT")
    let c2 = ToffoliDepthStreamed toffCirc2
    
    printfn "\n"
    printfn "doing 3"
    let a3 = toffLiqCirc3.GateCount(true, fun x -> x.Name = "CCNOT")
    let b3 = toffLiqCircFold3.GateCount(true, fun x -> x.Name = "CCNOT")
    let c3 = ToffoliDepthStreamed toffCirc3
    
    printfn "\n"
    printfn "doing 4"
    let a4 = toffLiqCirc4.GateCount(true, fun x -> x.Name = "CCNOT")
    let b4 = toffLiqCircFold4.GateCount(true, fun x -> x.Name = "CCNOT")
    let c4 = ToffoliDepthStreamed toffCirc4

    printfn "The Toffoli (unfolded) depths is %A %A %A %A" a1 a2 a3 a4
    printfn "The Toffoli (folded) depths is %A %A %A %A" b1 b2 b3 b4
    printfn "The Toffoli (streamed) depth is %A %A %A %A" c1 c2 c3 c4
    
    toffLiqCirc2.RenderHT("circ0")
    toffLiqCircFold2.RenderHT("circ1")

    show "done"


///////////////////////////////////////////////////////////////////////////
//
// Tests: Quantum circuit gadgets
//
///////////////////////////////////////////////////////////////////////////

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
// Tests: Quantum integer arithmetic
//
///////////////////////////////////////////////////////////////////////////

[<LQD>]
let __TestIntegerMultiplication (a:int) (b:int) =
    let n = 6 // # bits
    let k = Ket(3*n)
    let qs = k.Qubits

    let circuit = Circuit.Compile IntegerMultiplication qs
                    |> MyCircuitExport

    let k = Ket(qs.Length)
    let qs = k.Qubits
    // prepare x
    bigint a |> BoolInt n |> PrepBool qs 0 1 |> ignore
    bigint b |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(2*n-1) do 
        initialState.[i] <- qs.[i].Bit.v
    for i in 2*n..(initialState.Length-1) do
        initialState.[i] <- 0

    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in (2*n)..(3*n-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                |> List.sum
    show "%A * %A = %A" a b res
    

[<LQD>]
let __TestFixedPointMultiplication (a:int) (b:int) =
    let n = 6 // # bits
    let k = Ket(3*n)
    let qs = k.Qubits

    let circuit = Circuit.Compile FixedPointMultiplication qs
                    |> MyCircuitExport
    
    let k = Ket(qs.Length)
    let qs = k.Qubits
    // prepare x
    bigint a |> BoolInt n |> PrepBool qs 0 1 |> ignore
    bigint b |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(2*n-1) do 
        initialState.[i] <- qs.[i].Bit.v
    for i in 2*n..(initialState.Length-1) do
        initialState.[i] <- 0

    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in (2*n)..(3*n-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                |> List.sum
    show "%A * %A = %A" a b res
 
 
[<LQD>]
let __TestIntegerDivision (a:int) (b:int) =
    let n = 10 // # bits
    let k = Ket(3*n)
    let qs = k.Qubits

    let circuit = Circuit.Compile IntegerDivision qs
                    |> MyCircuitExport
    
    let k = Ket(qs.Length)
    let qs = k.Qubits
    // prepare x
    bigint a |> BoolInt n |> PrepBool qs 0 1 |> ignore
    bigint b |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(2*n-1) do 
        initialState.[i] <- qs.[i].Bit.v
    for i in 2*n..(initialState.Length-1) do
        initialState.[i] <- 0

    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in (2*n)..(3*n-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                |> List.sum
    let rem = [for i in (0)..(n-1) do yield (bigint finalState.[i] <<< i)]
                |> List.sum
    show "%A / %A = %A, rem %A" a b res rem
      

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
    printf "final state = %A" finalState
    (finalState.[2*n] = 1)

[<LQD>]
let __TestComparatorAllInputs (n:int) (s:int) (c:int) = 
    for i in 0..(pown 2 n)-1 do 
        let res = TestComparator n s i c
        show "Is %i > %i? : %b" s i res 


[<LQD>]
let __RunInplaceAdderFixedLength (n:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true
    let a = 6I
    let x = 11
    let k = Ket(n+numctrls+2)
    let qs = k.Qubits
    let circuit = Circuit.Compile (ConstantAdderInplaceCompile numctrls a) qs |> MyCircuitExport
    // Use the following function to print basic circuit information (qubits, size, depth)
    PrintToffoliNetworkMetrics circuit

// Adds (classical) value a to quantum register containing x using the garbage register g which 
// consists of two qubits and is left unchanged.
// NOTE: There is also a version which only requires 1 dirty qubit at a cost of an increased Toff count:
// 1 dirty qubit: 10 nlog(n), 2 dirty qubits: 8 nlog(n)
[<LQD>]
let __RunInplaceAdder (nmax:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true

    for n in 0..1..nmax do
        //let a = bigint ((1L <<< n)-1L)
        let a = 6I
            
        //let x = 212311%(1<<<n)
        let x = 11

        let k = Ket(n+numctrls+2)
        let qs = k.Qubits
        let circuit = Circuit.Compile (ConstantAdderInplaceCompile numctrls a) qs
                    |> MyCircuitExport

        
        // Use the following function to print basic circuit information (qubits, size, depth)
        PrintToffoliNetworkMetrics circuit

        // tweak the circuit into a Toffoli network that can be simulated efficienly    
        if simulate = true then
            let initialState = Array.zeroCreate (n+2+numctrls)
            let xb = BoolInt n (bigint x)
            for i in 0..(xb.Length-1) do 
                initialState.[i] <- xb.[i]
            initialState.[n] <- (g>>>1)&&&1
            initialState.[n+1] <- g &&& 1
            for i in 0..(numctrls-1) do
                initialState.[n+2+i] <- (ctrlvalue >>> i )&&&1
        
            let finalState = MyCircuitSimulateFast circuit initialState
            printfn "the constant is %A" a
            printfn "the initial state is %A" initialState
            printfn "the final state is %A" finalState
            let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                        |> List.sum
            let gout = (finalState.[n]+2*finalState.[n+1])
            if ( (int ((bigint x)+a))&&&((1<<<n)-1) = res) = false then
                show "%A + %A != %A" x a res
            if ( (int (gout)) = g) = false then
                show "Garbage changed!"
        let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
        show "Toffoli-count@n=\t%A\t%A" n toffcount
    ()

// Test for inplace Adder
[<LQD>]
let __RunInplaceAdderTest (nmax:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true

    for n in 0..nmax do
        show "Testing for bit size n = %A" n
        for a in 0..((1<<<n)-1) do
            for x in 0..((1<<<n)-1) do
                let k = Ket(n+numctrls+2)
                let qs = k.Qubits
                let circuit = Circuit.Compile (ConstantAdderInplaceCompile numctrls (bigint a)) qs
                            |> MyCircuitExport

                // tweak the circuit into a Toffoli network that can be simulated efficienly    
                let initialState = Array.zeroCreate (n+2+numctrls)
                let xb = BoolInt n (bigint x)
                for i in 0..(xb.Length-1) do 
                    initialState.[i] <- xb.[i]
                initialState.[n] <- ((g>>>1)&&&1)
                initialState.[n+1] <- (g &&& 1)
                for i in 0..(numctrls-1) do
                    initialState.[n+2+i] <- ((ctrlvalue >>> i )&&&1)
        
                let finalState = MyCircuitSimulateFast circuit initialState
    
                let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                            |> List.sum
                let gout = (finalState.[n]+2*finalState.[n+1])
                let ctrlout = [for i in 0..1..(numctrls-1) do yield (finalState.[i+n+2] <<< (i))]
                            |> List.sum
                if ( (int ((bigint x)+(bigint a)))&&&((1<<<n)-1) = res) = false then
                    show "%A + %A != %A, ctrl = %A" x a res ctrlout
                if ( (int (gout)) = g) = false then
                    show "Garbage changed!"
    ()

let TestConstantAdderInplace (n:int) (s1:int) (s2:int) (c:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let b  = double (c+1) |> log |> fun x -> x/(log 2.0) |> ceil |> int 
    let k = Ket(n+2+b)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt (n+1) |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs (n+1) 1 |> ignore
    c  |> fun x -> (bigint x) |> BoolInt b |> PrepBool qs (n+2) 1 |> ignore
    
    //t ConstantAdderCircuit = Circuit.Compile (ConstantAdderInplaceCompile b (bigint s2)) qs
    let ConstantSubtractorCircuit = Circuit.Compile (ConstantSubtractorInplaceCompile b (bigint s2)) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+2+b)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ConstantSubtractorCircuit initialState
    printf "final state = %A" finalState
    IntBool (n+1) finalState.[0..n]

[<LQD>]
let __TestConstantAdderInplaceAllInputs (n:int) (s:int) (c:int) = 
    for i in 0..(pown 2 n)-1 do 
        let res = TestConstantAdderInplace (n+1) s i c
        //show "Is %i + %i = %A? : %b" s i res ((bigint s) + (bigint i) = res)
        show "Is %i - %i = %A? : %b" s i res ((bigint s) - (bigint i) = res)

let TestConstantPrimeAdderInplace (n:int) (p:int) (s1:int) (s2:int) (c:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let b  = double (c+1) |> log |> fun x -> x/(log 2.0) |> ceil |> int 
    let k = Ket(2*n+3+b)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    s2 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs n 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 3 |> PrepBool qs (2*n) 1 |> ignore
    c  |> fun x -> (bigint x) |> BoolInt b |> PrepBool qs (2*n+3) 1 |> ignore
    
    // ConstantAdderCircuit = Circuit.Compile (ConstantAdderInplaceCompile b (bigint s2)) qs
    let ModularAdderCircuit = Circuit.Compile (CtrlModularAdderConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (2*n+3+b)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ModularAdderCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[n..(2*n-1)]

[<LQD>]
let __TestModularAdderConstantPrimeAllInputs (n:int) (p:int) (s:int) (c:int) = 
    for i in 0..p-1 do 
        let res = RunModularAdder "ModADDConstantPrimeConstantNumber" (n:int) (bigint s) (bigint i) (bigint p) (bigint c) "verbose" 
        show "Is %i + %i = %A? : %b" s i res ((((bigint s) + (bigint i)) % (bigint p)) = res)
        
[<LQD>]
let __TestModularAdderConstantPrimeConstantNumberAllInputs (n:int) (p:int) (s:int) (c:int) = 
    for i in 0..p-1 do 
        let res = TestConstantPrimeAdderInplace n p s i c
        show "Is %i + %i = %A? : %b" s i res ((((bigint s) + (bigint i)) % (bigint p)) = res)

let TestConstantPrimeDoublerInplace (n:int) (p:int) (s1:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let k = Ket(n+3)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    
    // ConstantAdderCircuit = Circuit.Compile (ConstantAdderInplaceCompile b (bigint s2)) qs
    let ModularDoublerCircuit = Circuit.Compile (ModularDoublerConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+3)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ModularDoublerCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[0..(n-1)]

[<LQD>]
let __TestModularDoublerConstantPrimeAllInputs (n:int) (p:int) = 
    for i in 0..p-1 do 
        let res = TestConstantPrimeDoublerInplace n p i
        show "Is 2*%i = %A? : %b" i res ((((bigint i) + (bigint i)) % (bigint p)) = res)


let TestCtrlConstantPrimeDoublerInplace (n:int) (p:int) (s1:int) (c:int) = 
    // Prepare the initial state
    let k = Ket(n+4)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    c  |> fun x -> (bigint x) |> BoolInt c |> PrepBool qs (n+3) 1 |> ignore
    
    // ConstantAdderCircuit = Circuit.Compile (ConstantAdderInplaceCompile b (bigint s2)) qs
    let ModularDoublerCircuit = Circuit.Compile (CtrlModularDoublerConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+4)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ModularDoublerCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[0..(n-1)]

[<LQD>]
let __TestCtrlModularDoublerConstantPrimeAllInputs (n:int) (p:int) (c:int) = 
    for i in 0..p-1 do 
        let res = TestCtrlConstantPrimeDoublerInplace n p i c
        show "Is %i + %i*%i = %A? : %b" i c i res ((((bigint i) + (bigint c)*(bigint i)) % (bigint p)) = res)


let TestConstantPrimeHalverInplace (n:int) (p:int) (s1:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let k = Ket(n+3)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    
    let ModularHalverCircuit = Circuit.Compile (ModularHalverConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+3)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast ModularHalverCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[0..(n-1)]

[<LQD>]
let __TestModularHalverConstantPrimeAllInputs (n:int) (p:int) = 
    for i in 0..p-1 do 
        let res = TestConstantPrimeHalverInplace n p i
        show "Is %i/2 = %A? : %b" i res (((res + res) % (bigint p)) = (bigint i))

let TestCtrlConstantPrimeHalverInplace (n:int) (p:int) (s1:int) (c:int) = 
    // Prepare the initial state
    //let k = Ket(2*n+2)
    let k = Ket(n+4)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 1 |> PrepBool qs n 1 |> ignore
    c  |> fun x -> (bigint x) |> BoolInt c |> PrepBool qs (n+3) 1 |> ignore
    
    let ModularHalverCircuit = Circuit.Compile (CtrlModularHalverConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (n+4)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    let finalState = MyCircuitSimulateFast ModularHalverCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[0..(n-1)]

[<LQD>]
let __TestCtrlModularHalverConstantPrimeAllInputs (n:int) (p:int) (c:int) = 
    for i in 0..p-1 do 
        let res = TestCtrlConstantPrimeHalverInplace n p i c
        show "Is %i = %A + %i*%A? : %b" i res c res (((res + (bigint c)*res) % (bigint p)) = (bigint i))

let TestConstantPrimeDBLADDMultiplier (n:int) (p:int) (s1:int) (s2:int) = 
    // Prepare the initial state
    let k = Ket(3*n+3)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    s2 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs n 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 3 |> PrepBool qs (2*n) (n+1) |> ignore
    
    let ModularMultiplierCircuit = Circuit.Compile (ModDblAddMultiplierConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly  
    
    PrintToffoliNetworkMetrics ModularMultiplierCircuit
    //show "Number of Toffoli gates = %A " (ModularMultiplierCircuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
      
    let initialState = Array.zeroCreate (3*n+3)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ModularMultiplierCircuit initialState
    //printf "final state = %A" finalState
    IntBool n finalState.[(2*n)..(3*n-1)]

[<LQD>]
let __TestDBLADDMultiplierConstantPrimeAllInputs (n:int) (p:int) (s:int) = 
    for i in 0..p-1 do 
        let res = TestConstantPrimeDBLADDMultiplier n p s i
        show "Is %i * %i = %A? : %b" s i res ((((bigint s) * (bigint i)) % (bigint p)) = res)

let TestConstantPrimeDBLADDSquarer (n:int) (p:int) (s1:int) = 
    // Prepare the initial state
    let k = Ket(2*n+5)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt (n+5) |> PrepBool qs n (n+5) |> ignore
    
    let ModularSquarerCircuit = Circuit.Compile (ModDblAddSquarerConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (2*n+5)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    
    //let finalState = MyCircuitSimulateFast ConstantAdderCircuit initialState
    let finalState = MyCircuitSimulateFast ModularSquarerCircuit initialState
    printf "final state = %A" finalState
    IntBool n finalState.[(n)..(2*n-1)]

[<LQD>]
let __TestDBLADDSquarerConstantPrimeAllInputs (n:int) (p:int) = 
    for i in 0..p-1 do 
        let res = TestConstantPrimeDBLADDSquarer n p i
        show "Is %i^2 = %A? : %b" i res ((((bigint i) * (bigint i)) % (bigint p)) = res)


[<LQD>]
let __RunIncrementerDirtyQubits (x:int) (g:int) = 
    show "Testing garbage incrementer."
    let bits x = (int (System.Math.Log((float x), 2.)) )+1
    let n : int = max (bits x) (bits g)
    let n : int = max n 2
    show "Doing computation with n = %A" n
    let k = Ket(2*n+1)
    let qs = k.Qubits
    
    let circuit = Circuit.Compile CtrlDirtyQubitsIncrementer qs
                    |> MyCircuitExport
    // prepare x
    x |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    // prepare garbage register
    g |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    initialState.[2*n] <- 0
    
    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                |> List.sum
    let gout = [for i in n..1..(2*n-1) do yield (finalState.[i] <<< (i-n))]
                |> List.sum
    show "Result of calculation: %A + 1 = %A" x res
    show "Garbage in: %A, Garbage out: %A" g gout
    ()


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

let RunModularAdderConstantTest n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running modular adder (constant p) circuit test for %A + %A mod %A" s1 s2 p
        show "==================================================================="
    let result =  RunModularAdder "ModADDConstantPrime" n s1 s2 p 1I verbosity
    if (s1 + s2) % p = result then 
        if verbose then 
            show "Test passed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s1 s2 p = %A %A %A %A. Result: sum %A\n" n s1 s2 p result 
        0 // return integer code zero in case of failure

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

let RunMultiplierTest (name:string) n s1 s2 p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running %A circuit test for %A * %A mod %A" name s1 s2 p
        show "====================================================="
    let result =  RunMultiplier name n s1 s2 p verbosity
    if name = "Montgomery" || name = "MontgomeryConstantPrime" then 
        let decres =  RunMultiplier name n result 1I p verbosity
        let decode1 = RunMultiplier name n s1 1I p verbosity
        let decode2 = RunMultiplier name n s2 1I p verbosity
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

let RunMultiplierTestWholeCircuit (name:string) n s1 s2 p verbosity =  
    let b, q, s, d =  RunMultiplierWholeCircuit name n s1 s2 p verbosity
    b, q, s, d

let TestConstantPrimeDBLADDMultiplierWholeCircuit (n:int) (p:int) (s1:int) (s2:int) = 
    // Prepare the initial state
    let k = Ket(3*n+3)
    let qs = k.Qubits
    s1 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    s2 |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs n 1 |> ignore
    0  |> fun x -> (bigint x) |> BoolInt 3 |> PrepBool qs (2*n) (n+1) |> ignore
    
    let ModularMultiplierCircuit = Circuit.Compile (ModDblAddMultiplierConstantPrime (bigint p) ) qs
                                    |> MyCircuitExport
    // tweak the circuit into a Toffoli network that can be simulated efficienly  
    
    let q, s, d = ReportToffoliNetworkMetrics ModularMultiplierCircuit
    n, q, s, d


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

//
//let RunSquarerTest (name:string) n s1 p verbosity =  
//    let verbose = (verbosity = "verbose") || (verbosity = "all")
//    if verbose then 
//        show "Running %A circuit test for %A^2 mod %A" name s1 p
//        show "====================================================="
//    let result =  RunMultiplier name n s1 0I p verbosity
//    if (s1 * s1) % p = result then 
//        if verbose then 
//            show "Test passed: n s1 p = %A %A %A. Results: product %A\n" n s1 p result
//        1 // return integer code one in case of success
//    else 
//        if verbose then 
//            show "Test failed: n s1 p = %A %A %A\n" n s1 p 
//        0 // return integer code zero in case of failure
//        
let RunSquarerTest name n s p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running Montgomery squarer circuit test for %A ^2 mod %A" s p
        show "====================================================="
    let result =  RunMultiplier "MontgomerySquarer" n s s p verbosity
    let decres =  RunMultiplier "Montgomery" n result 1I p "default"
    let decode = RunMultiplier "Montgomery" n s 1I p "default"
    if (decode * decode) % p = decres then 
        if verbose then 
            show "Test passed: n s p = %A %A %A. Results: product decprod decs %A %A %A \n" n s p result decres decode
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s p = %A %A %A\n" n s p 
        0 // return integer code zero in case of failure

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
    RunSquarerTest "MontgomerySquarer" 4 12I 13I "all" |> ignore
    RunSquarerTest "Montgomery" 4 11I 13I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 4 7I 13I "verbose"   |> ignore
    RunSquarerTest "Montgomery" 5 13I 29I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 5 10I 17I "verbose"  |> ignore
    RunSquarerTest "Montgomery" 6 7I 63I "verbose"   |> ignore

[<LQD>]
let __RunSmallSquarerTestsConstantPrime () = 
    RunSquarerTest "MontgomerySquarerConstantPrime" 4 12I 13I "all" |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 4 11I 13I "verbose"  |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 4 7I 13I "verbose"   |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 5 13I 29I "verbose"  |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 5 10I 17I "verbose"  |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 6 7I 63I "verbose"   |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 10 542I 1021I "verbose"       |> ignore
    RunSquarerTest "MontgomeryConstantPrime" 20 500I 1000003I "verbose"    |> ignore
    
[<LQD>]
let __RunSmallMultiplierTests () = 
//    show "\n"
    RunMultiplierTest "Montgomery" 3 5I 2I 7I "all"   |> ignore
    RunMultiplierTest "Montgomery" 10 542I 7I 1021I "verbose"       |> ignore
    RunMultiplierTest "Montgomery" 20 500I 7I 1000003I "verbose"    |> ignore
    //RunMultiplierTest "Mersenne" 3 5I 2I 7I "all"   |> ignore
    //RunMultiplierTest "DoubleAdd" 3 5I 2I 7I "all"   |> ignore
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
let __RunSmallMultiplierTestsConstantPrime () = 
    RunMultiplierTest "MontgomeryConstantPrime" 3 5I 2I 7I "all"   |> ignore
    RunMultiplierTest "MontgomeryConstantPrime" 10 542I 7I 1021I "verbose"       |> ignore
    RunMultiplierTest "MontgomeryConstantPrime" 20 500I 7I 1000003I "verbose"    |> ignore
    
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

let RunSuperEfficientInverterTest n p s verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running super efficient Montgomery inverter circuit test for %A mod %A" s p
        show "=========================================================================="
    let resultSM =  RunSuperEfficientInverter "MontgomerySM" n p s verbosity
    let resultMM =  RunSuperEfficientInverter "MontgomeryMM" n p s verbosity
    let resultSS =  RunSuperEfficientInverter "MontgomerySS" n p s verbosity   
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

let RunUltraEfficientInverterTest n p s verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running super efficient Montgomery inverter circuit test for %A mod %A" s p
        show "=========================================================================="
    let resultSM =  RunUltraEfficientInverter "MontgomerySM" n p s verbosity
    let resultMM =  RunUltraEfficientInverter "MontgomeryMM" n p s verbosity
    let resultSS =  RunUltraEfficientInverter "MontgomerySS" n p s verbosity   
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
    RunEfficientInverterTest 4 13I 8I "verbose" |> ignore
    RunEfficientInverterTest 5 17I 10I "verbose" |> ignore
    RunEfficientInverterTest 20 1000003I 500I "verbose" |> ignore    


[<LQD>]
let __RunSmallSuperEfficientInverterTests () = 
    RunSuperEfficientInverterTest 4 11I 6I "all" |> ignore
    RunSuperEfficientInverterTest 4 13I 8I "verbose" |> ignore
    RunSuperEfficientInverterTest 5 17I 10I "verbose" |> ignore
    RunSuperEfficientInverterTest 20 1000003I 500I "verbose" |> ignore    

[<LQD>]
let __RunSmallUltraEfficientInverterTests () = 
    RunUltraEfficientInverterTest 4 11I 6I "all" |> ignore
    RunUltraEfficientInverterTest 4 13I 8I "verbose" |> ignore
    RunUltraEfficientInverterTest 5 17I 10I "verbose" |> ignore
    RunUltraEfficientInverterTest 20 1000003I 500I "verbose" |> ignore    

[<LQD>]
let __RunScaledUltraEfficientInverterTests () = 
    RunUltraEfficientInverterTest 4 11I 6I "all" |> ignore
    RunUltraEfficientInverterTest 8 13I 8I "all" |> ignore
    RunUltraEfficientInverterTest 16 17I 10I "all" |> ignore
    RunUltraEfficientInverterTest 32 17I 50I "all" |> ignore    


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

[<LQD>]
let __RunLargeSuperEfficientInverterTests () = 
    let p192 = 6277101735386680763835789423207666416083908700390324961279I
    let p224 = 26959946667150639794667015087019630673557916260026308143510066298881I
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I
    let p384 = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319I
    let p521 = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151I 
    
    let watch = Diagnostics.Stopwatch()
    watch.Reset()
    watch.Start()
    RunSuperEfficientInverterTest 192 p192 2I "verbose" |> ignore
    watch.Stop()
    show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds

///////////////////////////////////////////////////////////////////////////
//
// Tests: Quantum modular arithmetic
//
///////////////////////////////////////////////////////////////////////////

// Runs modular multiplication for 2^k-bit numbers up to 16kbit 
// Reports Toffoli count and result of the computation & uncomputation (y should be 0 afterwards)
[<LQD>]
let __RunModMultiplier (i:int)=
    let N16k = 1189731495357231765085759326628007130763444687096510237472674821233261358180483686904488595472612039915115437484839309258897667381308687426274524698341565006080871634366004897522143251619531446845952345709482135847036647464830984784714280967845614138476044338404886122905286855313236158695999885790106357018120815363320780964323712757164290613406875202417365323950267880089067517372270610835647545755780793431622213451903817859630690311343850657539360649645193283178291767658965405285113556134369793281725888015908414675289832538063419234888599898980623114025121674472051872439321323198402942705341366951274739014593816898288994445173400364617928377138074411345791848573595077170437644191743889644885377684738322240608239079061399475675334739784016491742621485229014847672335977897158397334226349734811441653077758250988926030894789604676153104257260141806823027588003441951455327701598071281589597169413965608439504983171255062282026626200048042149808200002060993433681237623857880627479727072877482838438705048034164633337013385405998040701908662387301605018188262573723766279240798931717708807901740265407930976419648877869604017517691938687988088008944251258826969688364194133945780157844364946052713655454906327187428531895100278695119323496808703630436193927592692344820812834297364478686862064169042458555136532055050508189891866846863799917647547291371573500701015197559097453040033031520683518216494195636696077748110598284901343611469214274121810495077979275556645164983850062051066517084647369464036640569339464837172183352956873912042640003611618789278195710052094562761306703551840330110645101995435167626688669627763820604342480357906415354212732946756073006907088870496125050068156659252761297664065498347492661798824062312210409274584565587264846417650160123175874034726261957289081466197651553830744424709698634753627770356227126145052549125229448040149114795681359875968512808575244271871455454084894986155020794806980939215658055319165641681105966454159951476908583129721503298816585142073061480888021769818338417129396878371459575846052583142928447249703698548125295775920936450022651427249949580708203966082847550921891152133321048011973883636577825533325988852156325439335021315312134081390451021255363707903495916963125924201167877190108935255914539488216897117943269373608639074472792751116715127106396425081353553137213552890539802602978645319795100976432939091924660228878912900654210118287298298707382159717184569540515403029173307292433820279730892035211089568317921472966412729818117888454622835234533481549572743196734571634129403963684506083367368909845783940774676061371524391659800174583099528253482739757522802678131797463892371605031905393218989649729585777477658297715183959931786937290802047985884915659668363557899564591609695996145630027327988546172359576120674838045740053640396189717133917288042702062394781105995775027475794051243596930032470884082214287907917283972886566217666099928929736659305895401188482643507350834137594417765254338122252480578715215887377442751706844188139483627322937448198918074963942373224651327246325315257389993715375779389735367127951944954922262846304494110166728898867798152486928838574921534377364723691418182036282103573740854602704670142191797776310269453779790001271690513762843505679799117285781202701919548686629332950416374405421584385147988271411879919984987621119528596275512991867785554210230568216704790576492500881265464921474895559459290185727309851736982251106571217022721811129816003904051940819258133219624777452655172398944653308802365249212729736452165680043648457346007468020497859675817059619455164828406637438597547809499476599841667405756636926594683316183551527318347414506716327282224104675912411296799047208452645924316981738234697893596112182781775518514959248984587435401557122123845271399302125402393337922749394455393324480804314744560744187982171228782269878349258607387176290025943618757718706454898247726013046950133059199119527299386145131246957734198972549900572364145426997855548761897904248475449664993993862428294060559878742860518378058716440614430083917881823389319994271944195984564613910857996369511365798097068122437154330466861277932794002897625255843879189558691957373003353663905696047447765667201581975169298500940792736093264986662256899062338573645314621899957140466596072241410164464110856818238166950911474536753400753329335262575148713324594079410166286488150567351475933291175562137155949959416056516881232905680654150013798790046699053524795019980299093858029924914578907051418730788232093062110578135615360678690695102511579902265905593484061796350040671437474498455768307744249390537909776724456222613234554195683081115022999893977896277499918028773373149365321042211436343920206922932005094285083739085355321384854749817569151104942327971701097722576153091798048279364363438424522407742110033167468367004558861384900599732150535646653656240721316539419706881676263475309286637885522950176418670198965931I
    let N8k = 1090748135619415929462984244733782862448264161996232692431832786189721331849119295216264234525201987223957291796157025273109870820177184063610979765077554799078906298842192989538609825228048205159696851613591638196771886542609324560121290553901886301017900252535799917200010079600026535836800905297805880952350501630195475653911005312364560014847426035293551245843928918752768696279344088055617515694349945406677825140814900616105920256438504578013326493565836047242407382442812245131517757519164899226365743722432277368075027627883045206501792761700945699168497257879683851737049996900961120515655050115561271491492515342105748966629547032786321505730828430221664970324396138635251626409516168005427623435996308921691446181187406395310665404885739434832877428167407495370993511868756359970390117021823616749458620969857006263612082706715408157066575137281027022310927564910276759160520878304632411049364568754920967322982459184763427383790272448438018526977764941072715611580434690827459339991961414242741410599117426060556483763756314527611362658628383368621157993638020878537675545336789915694234433955666315070087213535470255670312004130725495834508357439653828936077080978550578912967907352780054935621561090795845172954115955492451891216570613296084365855602561987825554164828983776100070639535058975031421325993557963092582246894596377877165664134176626579977819203711571822341535408979557604056707123186626577676416188241336594881408068968961998717359787691335649357729429863437279636352583621433720199268128212734394658767583155547333863583775463140939229166786852053411359347506292273580759659221329615535915374296921312179266278677936047679491864303619458615535439684608832386143593871684717434990189386528537845712625304664823679920329431996145195217962866020168617725228731562929836932081530970550958194455848937712504788617035029734883083588942942311417818515459029999559814744757117521882790786148595159691476349244130338981981174851473059507740147241521672189074910330952319180638318199504090691810250168776309204618406888213015285193048095145396995645340002985115979602262321310662883088594437748318499511605309954705034226005722298162852349163242408574091397080892968895451060657757201199454907101458809408532429510189097387093136564036047797365062857395373088494023754833694236076111826663219393138920039495403085305056656519919814219844108778484235361412730008772810730108656701948930183929294294950489931276241738768497131921385470975219731917721I
    let N4k = 1044388881413152506691752710716624382579964249047383780384233483283953907971557456848826811934997558340890106714439262837987573438185793607263236087851365277945956976543709998340361590134383718314428070011855946226376318839397712745672334684344586617496807908705803704071284048740118609114467977783598029006686938976881787785946905630190260940599579453432823469303026696443059025015972399867714215541693835559885291486318237914434496734087811872639496475100189041349008417061675093668333850551032972088269550769983616369411933015213796825837188091833656751221318492846368125550225998300412344784862595674492194616443716246933213029777899798818461970932131191348804759825573468593610446932516843547031871082533988701912296616359089685540536703762087577894961935770084452075260024141226911707933733334829279200075988660067382055419272677487938772582789085572793356060472402630868147530739655123992685056875065481628365716808346058511893983996604674482947443330804864016159772117403184411098879310084911985820576707955476414019490642324082779864275147623753677862525374215759379112292614276897779676360046934875731128797682090481735253923263326681166800529934324565723067016719014885465039728219454090779142376015916232389820650894246737I
    let N2k = 25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357I
    let N1k = 135066410865995223349603216278805969938881475605667027524485143851526510604859533833940287150571909441798207282164471551373680419703964191743046496589274256239341020864383202110372958725762358509643110564073501508187510676594629205563685529475213500852879416377328533906109750544334999811150056977236890927563I
    let N512 = 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006083527I
    let Ns = [251I;65521I;4294967291I;18446744073709551557I;340282366920938463463374607431768211297I;115792089237316195423570985008687907853269984665640564039457584007913129639747I;N512;N1k;N2k;N4k;N8k;N16k]
    let ctrlcount = 1 // one control qubit (phase estimation)

    let N = Ns.[i]
    let x = (N >>> 5) + 5I
    let a = (N >>> 4) + 17I
    let n = int (System.Math.Ceiling(bigint.Log(N,2.)))
    let exact = ((a*x) % N)
    
    show "Running n = %A bits" n
    
    let k = Ket(2*n+1+ctrlcount)
    let qs = k.Qubits
    let initialState = Array.zeroCreate (qs.Length)
    let xb = BoolInt n x
    for i in 0..(xb.Length-1) do 
        initialState.[i] <- xb.[i]
        initialState.[i+n] <- 0
    initialState.[2*n] <- 0
    for i in 0..(ctrlcount-1) do
        initialState.[2*n+1+i] <- 1 // control qubit is 1 --> do it.

    let finalState = ConstantMultiplierModNRun ctrlcount a (ModularInverse a N) N qs initialState

    let res = [for i in 0..1..(n-1) do yield ((bigint finalState.[i]) <<< i)]
                                |> List.sum
    let y = [for i in 0..1..(n-1) do yield ((bigint finalState.[i+n]) <<< (i))]
                                |> List.sum
    show "Result: %A * %A mod %A = %A, y = %A" a x N res y
    show "Result - Exact: %A" (res-exact)
          
[<LQD>]
let __RunModMultiplierSweep (n:int) =
    for i in 3..n do 
        let k = Ket(2*i+1+1)
        let qs = k.Qubits
        let N = (pown 2I i)-1I
        let a = 3I
        let ctrlcount = 1
       
        ConstantMultiplierModNRunWholeCircuit ctrlcount a (ModularInverse a N) i N qs |> ignore
        show "done"  

[<LQD>]
let __RunModMultiplierSweepBatch (i:int) (runs:int) =
    let k = Ket(2*i+1+1)
    let qs = k.Qubits
    let N = (pown 2I i)-1I
    let a = 3I
    let ctrlcount = 1
    let mutable ave1 = 0.0
    let mutable ave2 = 0.0
    let mutable btot = 0
    let mutable qtot = 0
    let mutable stot = 0
    let mutable dtot = 0
    for j in 1..runs do
        let b, q, s, d, sim1, sim2 = ConstantMultiplierModNRunWholeCircuit ctrlcount a (ModularInverse a N) i N qs
        ave1 <- ave1 + sim1
        ave2 <- ave2 + sim2
        btot <- (int b)
        qtot <- q
        stot <- s
        dtot <- d
    show "CSVfinal: bitsize, qubits, size, depth, average-sim, average-total, %A, %A, %A, %A, %A, %A" btot qtot stot dtot (ave1/(float runs)) (ave2/(float runs))
    show "done"
    
[<LQD>]
let __RunModularNegatorTests () = 
    RunModularNegatorTest 4 3I 11I "all" |> ignore
    RunModularNegatorTest 4 8I 11I "all" |> ignore
    RunModularNegatorTest 5 8I 29I "all" |> ignore
    RunModularNegatorTest 6 8I 41I "all" |> ignore
    
let RunModularNegatorConstantPrimeTest n s p verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running controlled modular negator constant prime circuit test for - %A mod %A" s p
        show "==============================================================================="
    let result =  RunModularAdder "CtrlNegConstantPrime" n s s p 1I verbosity
    if (s + result) % p = 0I then 
        if verbose then 
            show "Test passed: n s p = %A %A %A. Result: inverse sum %A %A\n" n s p result ((s + result) % p)
        1 // return integer code one in case of success
    else 
        if verbose then 
            show "Test failed: n s p = %A %A %A. Result: inverse sum %A %A\n" n s p result ((s + result) % p)
        0 // return integer code zero in case of failure


[<LQD>]
let __RunModularNegatorConstantPrimeTests () = 
    RunModularNegatorConstantPrimeTest 4 3I 11I "all" |> ignore
    RunModularNegatorConstantPrimeTest 4 8I 11I "all" |> ignore
    RunModularNegatorConstantPrimeTest 5 8I 29I "all" |> ignore
    RunModularNegatorConstantPrimeTest 6 8I 41I "all" |> ignore
    RunModularNegatorConstantPrimeTest 6 8I 41I "all" |> ignore
    let p192 = 6277101735386680763835789423207666416083908700390324961279I
    RunModularNegatorTest 192 77766543I p192 "all" |> ignore

[<LQD>]
let __RunModularAdderSweep (n:int) =
    for i in 3..n do
        RunModularAdderTest i 3I 5I ((pown 2I i)-1I) "all" |> ignore
    
[<LQD>]
let __RunModularAdderConstantTests () = 
    RunModularAdderConstantTest 4 3I 5I 11I "all" |> ignore

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

// Adds (classical) value a to quantum register containing x using the dirty qubit register g which 
// consists of two qubits and is left unchanged.
[<LQD>]
let __RunInplaceAdderModN (nmax:int) (N:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = false
    let output = false//true

    let nmin = int (System.Math.Ceiling (System.Math.Log(float N, 2.)))
    if simulate = true then
        for n in nmin..nmax do
            for x in 0..(N-1) do
                for a in 0..(N-1) do
                    let k = Ket(n*2+numctrls)
                    let qs = k.Qubits
                    let circuit = Circuit.Compile (ConstantAdderInplaceModN numctrls (bigint a) (bigint N)) qs |> MyCircuitExport      
    
                    let initialState = Array.zeroCreate (qs.Length)
                    let xb = BoolInt n (bigint x)
                    initialState.[0] <- 0
                    for i in 0..(xb.Length-1) do 
                        initialState.[1+i] <- xb.[i]
                        if (i < xb.Length-1) then
                            initialState.[1+i+n] <- (g>>>i)&&&1
                    for i in 0..(numctrls-1) do
                        initialState.[2*n+i] <- (ctrlvalue >>> i )&&&1
        
                    let finalState = MyCircuitSimulateFast circuit initialState
    
                    let res = [for i in 1..1..(n) do yield (finalState.[i] <<< (i-1))]
                                |> List.sum
                    let gout = [for i in (n+1)..1..(2*n-1) do yield (finalState.[i] <<< (i-1-n))]
                                |> List.sum
                    if (((x + int a)%(1<<<n)) = res) = false && (((x + (int a) - N)+(1<<<n))%(1<<<n) = res) = false  then
                        if not (ctrlvalue &&& ((1<<<numctrls)-1) = ((1<<<numctrls)-1)) then
                            assert(x=res)
                            if (x = res) = false then
                                show "Control fail!"
                        else
                            show "%A + %A != %A for n = %A" x a res n
                    elif output then
                        show "%A + %A = %A :-)" x a res
                    if ( (int (gout)) = g) = false then
                        show "garbage %A !!!! %A" g gout
                    if finalState.[0] = 1 then
                        show "z = %A for x = %A and c = %A, res = %A" finalState.[0] x a res
    else // determine toffoli counts
        let x = 0
        for ni in nmax..nmax do
        //for ni in 0..100..nmax do
            let n = max ni 3
            let N = (1<<<n)-1          
            let a = N

            let k = Ket(n*2+numctrls)
            let qs = k.Qubits
            let circuit = Circuit.Compile (ConstantAdderInplaceModN numctrls (bigint a) (bigint N)) qs |> MyCircuitExport
            PrintToffoliNetworkMetrics circuit
            let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
            show "Toffoli-count@n=\t%A\t%A" n toffcount
    ()

    
///////////////////////////////////////////////////////////////////////////
//
// Tests: Quantum ECC point addition for affine Weierstrass curves
//
///////////////////////////////////////////////////////////////////////////

let RunECCPointAdditionTest name (n:int) P Q R (p:bigint) (b:bigint) verbosity =  
    let verbose = (verbosity = "verbose") || (verbosity = "all")
    if verbose then 
        show "Running ECC point addition for affine Weierstrass curve y^2 = x^3-3x+%A over GF(%A)" b p 
        show "====================================================================================="
        show "Adding the points %A + %A = %A." P Q R
    // let result0 =  RunShorEllipticCurvePointAddition name n p P Q 0I verbosity
    let result1, _, _, _, _ =  RunShorEllipticCurvePointAddition name n p P Q 1I verbosity
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
    RunECCPointAdditionTest "AffineWeierstrass" 4 P Q R p b "all" |> ignore
    
//    let p = 1000003I
//    let b = 5I
//    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=5 over GF(1000003)
//    let P = AffinePoint(34236I,177368I) // MG representation corresponding to (563093 : 779719 : 1)
//    let Q = AffinePoint(481427I,185548I) // MG representation corresponding to (159317 : 146176 : 1)
//    let R = AffinePoint(39120I,600367I) // MG representation corresponding to (187760 : 407071 : 1)
//    RunECCPointAdditionTest "AffineWeierstrass" 20 P Q R p b "verbose" |> ignore

[<LQD>]
let __RunSmallECCPointAdditionTestsConstantPrime () = 
    let p = 11I 
    let b = 5I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=5 over GF(11)
    let P = AffinePoint(0I,2I) // MG representation corresponding to (0 : 7 : 1) 
    let Q = AffinePoint(5I,8I) // MG representation corresponding to (1 : 6 : 1) 
    let R = AffinePoint(0I,9I) // MG representation corresponding to (0 : 4 : 1) 
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 4 P Q R p b "verbose" |> ignore
    
    let p = 1000003I
    let b = 5I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=5 over GF(1000003)
    let P = AffinePoint(34236I,177368I) // MG representation corresponding to (563093 : 779719 : 1)
    let Q = AffinePoint(481427I,185548I) // MG representation corresponding to (159317 : 146176 : 1)
    let R = AffinePoint(39120I,600367I) // MG representation corresponding to (187760 : 407071 : 1)
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 20 P Q R p b "verbose" |> ignore


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
let __RunMediumECCPointAdditionTestsConstantPrime () = 
    let p = 1125899906842679I // 51 bit prime number
    let b = 111I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=111 over GF(112589990684267)
    let P = AffinePoint(451254016350526I,119715075055576I) // MG representation corresponding to (36839505554729 : 449271643872930 : 1)
    let Q = AffinePoint(102913185642772I,40875630210012I) // MG representation corresponding to (1002138888044907 : 490930180893078 : 1)
    let R = AffinePoint(939594529930895I,740730201703437I) // MG representation corresponding to (247344574010328 : 331036061128227 : 1)
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 51 P Q R p b "verbose" |> ignore

    let p = 1267650600228229401496703205653I // 101 bit prime number
    let b = 500I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b=500 over GF(1267650600228229401496703205653)
    let P = AffinePoint(1050973377032356303253151680466I,901012499748817281889417696786I) // MG representation corresponding to (332176921040233098041995516868 : 1211107952384860587565637691876 : 1)
    let Q = AffinePoint(1193382548062512834590307231000I,282766732732562641760008066570I) // MG representation corresponding to (784979104567596392202663529808 : 777470103510587425752922530425 : 1)
    let R = AffinePoint(52487571062535416645933108150I,549778904875717537058134673149I) // MG representation corresponding to (1025008991572534722840933398961 : 1218606481979730240903762949350 : 1)
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 101 P Q R p b "verbose" |> ignore
    
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

[<LQD>]
let __RunLargeECCPointAdditionTestsConstantPrime () = 
    let p256 = 115792089210356248762697446949407573530086143415290314195533631308867097853951I // 256 bit prime number
    //let p384 = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319I // 384 bit prime number
    let b = 41058363725152142129326129780047268409114441015993725554835256314039467401291I
    // we now consider the curve y^2 = x^3 +ax + b, where a=-3 and b as above over GF(p256)
    // MG representation corresponding to (66790887510078570501074792547226042371050454650016343615879109835932523443065 : 106377730366318381478632396320403129066886107237681426316144016021673156344273 : 1)
    let P = AffinePoint(71528713545620684423536968445225241465408169212990955778170977604457805197600I,103275792837057906611974512458293333482266004537501318958510683525869453045362I) 
    // MG representation corresponding to (91776155225397256558779214231279668060857253502872526815652263470579408313795 : 103169371820847966113293214366645798727644724219171959581593614101799275769904 : 1)
    let Q = AffinePoint(55060438472976370102375399605978382023362429340763975384115074388233082886021I,31367148643791062976004566319732149799763757143358924179768669057101515390805I) 
    // MG representation corresponding to (455438157303235004258818524917930849403739804230436973770526759698198641623 : 107400709691385187225095324036499254562775224598883627815051692242731812890812 : 1)    
    let R = AffinePoint(76441926148553224235917927366185924393481358162275007948038451322418981861496I,83037667562094125842330603399066966639949874188752994804195456661992373332472I)
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 256 P Q R p256 b "verbose" |> ignore

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
    RunECCPointAdditionTest "AffineWeierstrassConstantPrime" 521 P Q R p521 b "verbose" |> ignore

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
let __RunLargeFactoringTests () = // tested ok
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
let __RunVeryLargeFactoringTests () = // not yet tested
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
// Read ECC data file with tests
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
        let result1, _, _, _, _ =  RunShorEllipticCurvePointAddition "AffineWeierstrass" (int n) p pm qm 1I "none"
        match (result1 = rm) with
        | true -> show "Test passed." 
        | _ -> show "Test failed: p=%A n=%A a=%A n=%A p=%A q=%A r=%A" p n a b pm qm rm   
        //show "result0=%A result1=%A R=%A.\n" result0 result1 R        
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
        
    member this.RunConstantPrime() =
        let watch = Diagnostics.Stopwatch()
        watch.Start()
        show "Running ECC point addition (constant prime) for affine Weierstrass curve y^2 = x^3-%Ax+%A over GF(%A)" a b p 
        show "====================================================================================================="
        show "Adding the points %A + %A = %A." pm qm rm
        let result1, qubits, size, depth, sim =  RunShorEllipticCurvePointAddition "AffineWeierstrassConstantPrime" (int n) p pm qm 1I "none"
        match (result1 = rm) with
        | true -> show "Test passed." 
        | _ -> show "Test failed: p=%A n=%A a=%A n=%A p=%A q=%A r=%A" p n a b pm qm rm   
        //show "result0=%A result1=%A R=%A.\n" result0 result1 R        
        watch.Stop()
        show "Time to compile & simulate: %f" watch.Elapsed.TotalSeconds
        n, qubits, size, depth, sim, watch.Elapsed.TotalSeconds
    
[<LQD>]
let __RunMeasureECCPointAdditionBatch (path:string) =
    let lines = File.ReadAllLines path
   
    for data in lines do 
        let T = new eccinstance()
        T.Init(data)
        T.Run()                
    show "done"

[<LQD>]
let __RunMeasureECCPointAdditionConstantPrimeBatch (path:string) (n:int) =
    let lines = path |> File.ReadAllLines |> Array.take n
    let mutable ave1 = 0.0
    let mutable ave2 = 0.0
    let mutable btot = 0
    let mutable qtot = 0
    let mutable stot = 0
    let mutable dtot = 0
    
    for data in lines do 
        let T = new eccinstance()
        T.Init(data)
        let b, q, s, d, sim1, sim2 = T.RunConstantPrime()        
        ave1 <- ave1 + sim1
        ave2 <- ave2 + sim2
        btot <- (int b)
        qtot <- q
        stot <- s
        dtot <- d
    show "CSVfinal: bitsize, qubits, size, depth, average-sim, average-total, %A, %A, %A, %A, %A, %A" btot qtot stot dtot (ave1/(float n)) (ave2/(float n))
    show "done"

///////////////////////////////////////////////////////////////////////////
// Read RSA data file with tests
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

[<LQD>]
let __ReadQCCircuit (path:string) =
    let lines = File.ReadAllLines path
    for data in lines do 
        let T = new RSAinstance()
        T.Init(data)
        T.Run()
    show "done"
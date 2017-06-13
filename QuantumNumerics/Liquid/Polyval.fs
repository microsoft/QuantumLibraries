module Microsoft.Research.Liquid.Polyval

open System
open System.Collections.Generic
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

open Montgomery
open Inplace


// computes the label for the x register
// i.e., computes the interval label s.t. x \in I_{label}
//
// intervals are given as [a_i, b_i) with strictly increasing values
// and x >= a_0 is assumed. Furthermore, at least the highest bit of 
// x must be 0 (in order to check for overflow).
//
// n : number of bits to represent x
// pointpos: binary point position of x
// b is a list of right interval boundaries
// l are log(M) qubits, where M is the number of intervals
// x are n qubits representing the value x
// g are n dirty qubits
let ComputeLabelUsingGarbage (n:int) (pointpos:int) (a:list<double>) (l:Qubits) (x:Qubits) (g:Qubits) = 
    let one = bigint 1
    let ComputeCarry (c:bigint) (x:Qubits) =
        let r = x.Length - 1
        if r = 1 then
            if (c&&&one) = one then
                CNOT [x.[0];x.[1]]
        else
            CNOT [g.[r-2];x.[r]]
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

            CNOT [g.[r-2];x.[r]]

    let UncomputeCarry (c:bigint) (x:Qubits) =
        let r = x.Length - 1
        if r = 1 then
            if (c&&&one) = one then
                CNOT [x.[0];x.[1]]
        else
            CNOT [g.[r-2];x.[r]]
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
            CNOT [g.[r-2];x.[r]]
    
    for i in 0..n-1 do
        X [x.[i]]
    for m in 0..a.Length-1 do
        // compare x to the m-th value
        let c = bigint (a.[m] * double (one <<< (n - pointpos)))
        ComputeCarry (c) x
        for bit in 0..(l.Length-1) do
            if (((m^^^(m+1))>>>bit)&&&1) = 1 then
                CNOT [x.[x.Length-1]; l.[bit]]
        UncomputeCarry (c) x

    for i in 0..n-1 do
        X [x.[i]]

// qs = [LABEL (log M), ANC1 (1), X (n), RESULT (n), ANC2 (n*(d-1)+1)]
let ParallelPolyvalCompile (n:int) (pointpos:int) (b:list<double>) (coeffs:list<list<double>>) (qs:Qubits) = 
    let logM = int (Math.Ceiling (Math.Log((double b.Length), 2.)))
    ComputeLabelUsingGarbage n pointpos b qs.[0..logM-1] qs.[logM+1..logM+n] qs.[logM+n+1..logM+2*n]
    // last argument to parallel polyval are qubits: [log(M),1,n (x), n*(d-1), n (result), 1]
    let d = coeffs.[0].Length-1
    EvaluatePolynomialsInParallel pointpos n coeffs (List.concat [qs.[0..logM+n];qs.[logM+2*n+1..logM+2*n+n*(d-1)];qs.[logM+n+1..logM+2*n];[qs.[logM+2*n+n*(d-1)+1]]])

[<LQD>]
let __ParallelPolyval (x:double) =
    let n = 30
    let pointpos = 7

    let b = [0.1; 0.5; 1.]
    let coeffs = [[1.;2.0;3.;4.0];[1.;2.;3.;5.0];[1.;2.;3.;5.0]]
    let deg = coeffs.[0].Length-1
    let logM = 2
    
    let k = Ket(logM+1 + n + n * (deg) + 1)
    let qs = k.Qubits
    let circuit = Circuit.Compile (ParallelPolyvalCompile n pointpos b coeffs) qs
                    |> MyCircuitExport

    (bigint (x*double (1I<<<(n-pointpos)))) |> BoolInt n |> PrepBool qs (logM+1) 1 |> ignore

    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(initialState.Length-1) do
        initialState.[i] <- 0
    for i in 0..(n-1) do 
        initialState.[logM+i+1] <- qs.[logM+1+i].Bit.v

    //show "Calculating P(x) = %A" expres
    let finalState = MyCircuitSimulateFast circuit initialState
    let lbl = [for i in (0)..(logM) do yield (bigint finalState.[i] <<< (i))]
                |> List.sum
    let xres = [for i in 0..n-1 do yield (bigint finalState.[i+logM+1] <<< i)] |> List.sum
    let res = [for i in 0..n-1 do yield (bigint finalState.[i+logM+n+1] <<< i)] |> List.sum
    let rest = [for i in 0..(deg-1)*n do yield (bigint finalState.[i+logM+2*n+1] <<< i)] |> List.sum
    let conv (arg:bigint) = ((double arg)/(double (1I <<< (n-pointpos))))
    let mutable true_res = x**3.0 + 2.*x**2. + 3.*x + 4.
    if x >= 0.1 then
        true_res <- true_res + 1.
    printfn "%A\t%1.7f\t%1.7f\t%1.2f\t%A\t" lbl (conv xres) (conv res) true_res rest
module Microsoft.Research.Liquid.DataTypes

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open CircBase           // basic quantum circuits
open Integer            // integer arithmetic

//
// SORTING
//
//
let QBubbleSort (bits:int) (num:int) (qs:Qubits) =
    let mutable anc_index = 0
    for it in 0..num-1 do
        for k in 0..num-2-it do
            let anc = [qs.[bits * num + anc_index]]
            let a = qs.[k * bits..k * bits + bits - 1]
            let b = qs.[(k + 1) * bits .. (k + 1) * bits + bits - 1]

            BuildTakahashiAdderInverse a (List.concat [b; anc])
            BuildTakahashiModAdder a b
            CSWAP a b anc
            anc_index <- anc_index + 1
    ()

let QMergeSort (bits:int) (num:int) (qs:Qubits) =
    ()

[<LQD>]
let __RunBSort () =     
    let g = 3
    let n = 5
    let a = [1;5;3;8;3;2;5;1;9;4]

    let k = Ket(n*a.Length + a.Length * (a.Length - 1) / 2)
    let qs = k.Qubits
    let circuit = Circuit.Compile (QBubbleSort n a.Length) qs
                |> MyCircuitExport

    let initialState = Array.zeroCreate (qs.Length)
    for i in 0..a.Length-1 do
        for b in 0..n-1 do
            initialState.[i * n + b] <- (a.[i] >>> b) &&& 1
        
    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in 0..a.Length-1 do yield ([for b in 0..(n-1) do yield (finalState.[i * n + b] <<< b)]
                                                    |> List.sum)]
    show "%A" res
    ()

[<LQD>]
let __RunMSort () =     
    let g = 3
    let n = 5
    let a = [1;5;3;8;3;2;5;1]

    let k = Ket(n*a.Length + a.Length * (a.Length - 1) / 2)
    let qs = k.Qubits
    let circuit = Circuit.Compile (QBubbleSort n a.Length) qs
                |> MyCircuitExport

    let initialState = Array.zeroCreate (qs.Length)
    for i in 0..a.Length-1 do
        for b in 0..n-1 do
            initialState.[i * n + b] <- (a.[i] >>> b) &&& 1
        
    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in 0..a.Length-1 do yield ([for b in 0..(n-1) do yield (finalState.[i * n + b] <<< b)]
                                                    |> List.sum)]
    show "%A" res
    ()
// Sample user module to link into the LiquidX application

module Microsoft.Research.Liquid.UserSample

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions

open Util               // General utilities
open Operations         // Gate definitions
open Shor

// Optional extras:
//open Native             // Support for Native Interop
//open ARPACK             // Calls to ARPACK for sparse linear algebra
//open HamiltonianGates   // Extra gates for doing Hamiltonian simulations
//open HubbardGates       // Extra gates for doing Hubbard simulations
//open Shor               // Shor's algoritmm implementation
//open Tests              // All the built-in tests

[<LQD>]
let UserSample(cnt:int) =
    let k   = Ket(cnt)
    let qs  = k.Qubits
    let test (qs:Qubits) =
        X !!(qs,0)
        Y !!(qs,1)
        for i in 0..qs.Length-1 do LabelL (sprintf "\\ket{%d}" i) !!(qs,i)
        H >< qs
        for i in 1..2..qs.Length-1 do CNOT !!(qs,i-1,i)
        for i in 2..2..qs.Length-1 do CNOT !!(qs,i-1,i)
        QFT qs
        QFT' qs
    
    let circ    = Circuit.Compile test qs
    //let circ    = circ.Fold()
    //let gp      = GrowPars(20,2,true)
    //let circ    = circ.GrowGates(k,gp)
    //circ.Dump()
    circ.RenderHT "Test"
    circ.Run qs
  //  let circ = Circuit.Compile QFT qs
  //  circ.RenderHT "Test"
 //   let k   = k.Single()
 //   k.Dump(showInd)

/// <summary>
/// Take QFT; local copy 
/// </summary>
/// <param name="qs">Qubits to take QFT of</param>
let QFT (qs:Qubits) =
    let gate (qs:Qubits)    =
        let op (qs:Qubits) =
            let n   = qs.Length
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
/// HHL_BandMat yields the quantum circuit for Harrow, Hassidim, Lloyd algorithm for the special 
/// case of band matrices of the form [[1,1,0,...,0],[1,0,1...,0],[0,1,0,1,..],[0,..,1,-1].                
/// </summary>
(*[<LQD>]
let HHL_BandMat() =
  let A = 1
  A 
*)

/// <summary>
/// DCT yields circuit for Discrete Cosine Transform (type IV) on n qubits.
/// </summary>
/// Typs: 
(*[<LQD>]
let DCT(n:int) =
  let k = Ket(n)
  let qs = k.Qubits  
  let B = H * S
  B >< qs
*)

// Experiments with tail recursion

let fact x = 
    let rec fact x i = 
        match x with 
        | 0|1 -> i
        | _ -> fact (x-1) (x*i)
    fact x 1 

let sqrt2 =
    let rec calc i v =
        show "Level is %d and current value is %.5f%%" i v
        if i = 2 then v
        else calc (i+1) (0.5*v + 1./v)
    calc 1 (sqrt 2.0)

let logab a b = 
    (log a)/(log b)
    
//let myNewton x s = 
//    let p = (logab x 2.0) |> ceil |> int
//    let xinit = 1.0 
//    let rec bisect s i = 
//        match s with  
//        | 0|1 -> i
//        | _ -> 5*i
//       | _ -> bisect (s-1) (-x*(pown i 2)+2*i)
//    bisect s xinit
    
[<LQD>]
let mtest(u:int) = 
    let mvar = fact u
    let uvar = sqrt2
    let vvar = sqrt 2.0
    show "Factorial of %d is %d" u mvar
    show "Square root of 2 approx is %.30f" uvar
    show "Square root of 2 approx is %.30f" vvar
    show "Log[2](5) is %.30f" (logab 5.0 2.0)
//    show "myNewton of %.5f with %.1f iterations is %0.20f" 3.0 5.0 (myNewton 3.0 5)
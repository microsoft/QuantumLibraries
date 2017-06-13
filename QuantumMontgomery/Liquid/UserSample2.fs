// Sample user module to link into the LiquidX application

module Microsoft.Research.Liquid.UserSample

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions

open Util               // General utilities
open Operations         // Gate definitions

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
    let circ    = Circuit.Compile test qs
    let circ    = circ.Fold()
    //let gp      = GrowPars(20,2,true)
    //let circ    = circ.GrowGates(k,gp)
    circ.Dump()
    circ.RenderHT "Test"
    circ.Run qs
    let k   = k.Single()
    k.Dump(showInd)


/// <summary>
/// HHL_BandMat yields the quantum circuit for Harrow, Hassidim, Lloyd algorithm for the special 
/// case of band matrices of the form [[1,1,0,...,0],[1,0,1...,0],[0,1,0,1,..],[0,..,1,-1].                
/// </summary>
/// [<LQD>]
///let HHL_BandMat() =
///  let A = 1
///  A 

/// <summary>
/// DCT yields circuit for Discrete Cosine Transform (type IV) on n qubits.
/// </summary>
/// Typs: 
[<LQD>]
let M (qs:Qubits) = 
    let gate =
        Gate.Build("M",fun () -> 
            new Gate(
                Name    = "Meads",
                Help    = "Collapse State",
                Mat     = CSMat(2),
                Draw    = "\\meter",
                Op      = Measure
        ))
    gate.Run qs

[<LQD>]
let B (qs:Qubits) = 
    let gate = 
        Gate.Build("B",fun() -> 
            new Gate(
                Mat     = (
                    let vR  = 1.0/sqrt2
                    let vI  = 0.0
                    let wR  = 0.0
                    let wI  = 1.0/sqrt2
                    CSMat(2,[(0,0,vR,vI);(0,1,wR,wI);(1,0,vR,vI);(1,1,wR,wI)]))
            )
        )
    gate.Run qs 

[<LQD>]
let DCT(n:int) =
    let k = Ket(n)
    let qs = k.Qubits  
    let circ = Circuit.Compile B qs
    circ.Run qs 



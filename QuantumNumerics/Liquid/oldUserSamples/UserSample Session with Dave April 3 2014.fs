// Sample user module to link into the LiquidX application

module Microsoft.Research.Liquid.UserSample

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions

open Util               // General utilities
open Operations         // Gate definitions
//open Shor

// Optional extras:
//open Native             // Support for Native Interop
//open ARPACK             // Calls to ARPACK for sparse linear algebra
//open HamiltonianGates   // Extra gates for doing Hamiltonian simulations
//open HubbardGates       // Extra gates for doing Hubbard simulations
//open Shor               // Shor's algoritmm implementation
//open Tests              // All the built-in tests

//let fact x = 
//    let rec fact x i = 
//        match x with 
//        | 0|1 -> i
//        | _ -> fact (x-1) (x*i)
//    fact x 1 
//
//let sqrt2 =
//    let rec calc i v =
//        show "Level is %d and current value is %.5f%%" i v
//        if i = 2 then v
//        else calc (i+1) (0.5*v + 1./v)
//    calc 1 (sqrt 2.0)

// RUS circuit, April 2014
//
//[<LQD>]
//let R (k:int) (qs:Qubits) =
//    let gate (k:int) =
//        Gate.Build("R_" + k.ToString(),fun () ->
//            new Gate(
//                Name    = sprintf "R%d" k,
//                Help    = sprintf "2pi/2^k: %d" k,
//                Mat     = (
//                    let phi     = (2.0*Math.PI)/(pown 2.0 k)
//                    let phiR    = Math.Cos phi
//                    let phiI    = Math.Sin phi
//                    CSMat(2,[(0,0,1.,0.);(1,1,phiR,phiI)])),
//                Draw    = "\\gate{" + ("R" + (k.ToString())) + "}"
//               ))
//    (gate k).Run qs

let rdv =
    let rnd = Random(1234)
    fun () -> rnd.NextDouble()

[<LQD>]
let Ugate (qs:Qubits) =
    let gate = 
        Gate.Build("Ugate", fun() -> 
            new Gate(
                Qubits  = qs.Length, 
                Name    = "Ugate", 
                Help    = "random local SU(2) unitary",
                Mat     = (
                    let a   = rdv()
                    let br  = rdv()
                    let bi  = rdv()
                    let mat = CMat(2,[0,0,a,0.0;0,1,-br,bi;1,0,br,bi;1,1,a,0.0]) 
                    show "@@@DBG initial error = %.11g" (mat.UnitaryError())
                    let iter,err = mat.FixUnitary(1.0e-13,20)
                    show "@@@DBG final   error = %.11g iters = %d" err iter
                    CSMat(mat)
                    ),
                Draw    = "\\gate{U}"
                ))
    gate.Run qs

[<LQD>]
let Martin() =
    let ket         = Ket(4)
    show "Dumping ket:"
    ket.Dump()
    let cv      = ket.Single()
    show "Dumping vector:"
    cv.Dump()
    show "Dumping ket:"
    ket.Dump()
    let qs      = ket.Qubits
    let ops (qs:Qubits) =
        H !!(qs,2)
        CNOT !!(qs,2,3)
        X !!(qs,0)
        (Cgate Z) !!(qs,0,3)
        X !!(qs,0)
        Ugate !!(qs,1)
        (Adj Ugate) !!(qs,2)
        Ugate !!(qs,3)
        (Cgate SWAP) !!(qs,0,1,2)
        X !!(qs,0)
        (Cgate Z) !!(qs,0,3)
        X !!(qs,0)
        CNOT !!(qs,2,3)
        H !!(qs,2)
    let circ    = Circuit.Compile ops qs
    let circ    = circ.Fold()
    circ.RenderHT("Ugates")
    circ.Dump()

    circ.Run qs
    show "Dumping ket after run:"
    ket.Dump(showInd,0,true,true)
    show "Prob of measuring each qubit:"
    for q in qs do
        show "Q%d = %.12g" q.Id q.Prob1
    show "Done"


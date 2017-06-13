module Microsoft.Research.Liquid.UserSample

open System
open System.Collections.Generic
open System.IO
open System.Text
open System.Text.RegularExpressions
open System.Drawing
open System.Windows.Forms
open System.Windows.Forms.DataVisualization.Charting
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open Shor               // factoring using semi-classical qft
open Plots              // .net based plot lib written by Dave

open GaussSums
open HiddenShift
open Montgomery
open SignalProcessing

open AltOutput

//open Adders

///////////////////////////////////////////////////////////////////////////////
//
// Symplectic canonical form
//
///////////////////////////////////////////////////////////////////////////////

//type bit = { true | false }
//    member x.(+) (a:bit) (b:bit) = 
//        a ^^^ b
//    member x.(*) (a:bit) (b:bit) = 
//        a &&& b
//
//
//type SquareMat private(A : [] of [] of bit) = 
//    member x.LPUdecomp() =
//        let columnadd i j = 
//            for k in 1..n do 
//                A.[let n = Length(A)
//        let J = []
//        for j in 1..n do 
//            let a = A.[1..n,j]
//            let ind = Array.tryFind (fun x -> x <> 0) a 
//            if Some(ind) then
//                List.append J [[ind;j]]                 
//                for i in ind..n do 
//                    rowadd 
//type Symplectic private(A, B, C, D)
//    member x.

///////////////////////////////////////////////////////////////////////////////
//
// Depth-optimized Repeat-Until-Success implementations, April 2014
// Literature: [Bocharov, Roetteler, Svore, arxiv:1404.5320], Section App I.3
//
///////////////////////////////////////////////////////////////////////////////

// Random number generator
let rdv =
    let rnd = Random(1234)
    fun () -> rnd.NextDouble()

// Applies a single qubit rotation around the X axis with angle <phi>. 
let Xrotation (phi:float) (qs:Qubits) =
    let gate (phi:float) =
        new Gate(
            Qubits  = qs.Length, 
            Name    = sprintf "Xrotation around %0.4f" phi, 
            Help    = "predefined X rotation",
            Mat     = (             
                let a = Math.Cos phi
                let b = Math.Sin phi
                let mat = CMat(2,[0,0,a,0.0;0,1,-b,0.0;1,0,b,0.0;1,1,a,0.0]) 
                CSMat(mat)
                ),
            Draw    = "\\gate{U_x}"
            )
    (gate phi).Run qs

// Applies a single qubit rotation around the Z axis with angle <phi>. 
let Zrotation (phi:float) (qs:Qubits) =
    let gate (phi:float) = 
        new Gate(
            Qubits  = qs.Length, 
            Name    = sprintf "Zrotation around %0.4f" phi, 
            Help    = "predefined Z rotation",
            Mat     = (             
                let a = Math.Cos phi
                let b = Math.Sin phi
                let mat = CMat(2,[0,0,a,b;0,1,0.0,0.0;1,0,0.0,0.0;1,1,a,-b]) 
                CSMat(mat)
                ),
            Draw    = "\\gate{U_z}"
            )
    (gate phi).Run qs
      
// Create a single qubit SU(2) unitary from three given Euler angles 
let EulerAngleSU2Rotation (phi_x1:float) (phi_z:float) (phi_x2:float) (qs:Qubits) =     
    Zrotation phi_x1 qs
    Xrotation phi_z qs
    Zrotation phi_x2 qs
    
// Create a inverse of single qubit SU(2) unitary from three given Euler angles 
let EulerAngleSU2Rotation' (phi_x1:float) (phi_z:float) (phi_x2:float) (qs:Qubits) =     
    Zrotation -phi_x2 qs
    Xrotation -phi_z qs
    Zrotation -phi_x1 qs
    
[<LQD>]
// RUSCircuitTest() tests the functionality of a Jack-Of-Dagger design (see [BRS'14]) that 
// allows to implement axial rotations probabilistically.
let __RUSCircuitTest() = 
    let ket         = Ket(4)
    show "Dumping ket before run:"
    ket.Dump(showInd,0,false,true)
    let qs      = ket.Qubits
    
    let a   = rdv()*Math.PI
    let b   = rdv()*Math.PI
    let c   = rdv()*Math.PI
    let d   = rdv()*Math.PI
    let e   = rdv()*Math.PI
    
    let MyU = EulerAngleSU2Rotation a b -a     // the valid SU(2) elements for which circuit works
    let MyU' = EulerAngleSU2Rotation' a b -a   // the inverse of the previously defined SU(2) element
    let MyRand = EulerAngleSU2Rotation c d e    // circuit to rotate initial qubit; prepares random state 
    let MyRand' = EulerAngleSU2Rotation' c d e  // and the inverse to check if qubit is returned in 0
    
    // circuit to compute Jack-Of-Daggers matrices for SU(2) elements as above
    let ops (qs:Qubits) =
        X !!(qs,0) //this checks the "b=1" branch of the circuit. For the other branch, comment this line
        MyRand !!(qs,1)
        H !!(qs,2)
        X !!(qs,3)
        CNOT !!(qs,2,3)
        (Cgate Z) !!(qs,0,2)
        (Cgate SWAP) !!(qs,0,1,2)
        MyU !!(qs,1)
        MyU' !!(qs,2)
        MyU !!(qs,3)
        (Cgate SWAP) !!(qs,0,1,2)
        (Cgate Z) !!(qs,0,2)
        CNOT !!(qs,2,3)
        X !!(qs,3)
        H !!(qs,2)
        MyU !!(qs,1) //this checks the "b=1" branch of the circuit. For the other branch, use MyU' here
        MyRand' !!(qs,1)

    let circ  = Circuit.Compile ops qs
    // let circ    = circ.Fold()
    circ.RenderHT("Ugates")
    circ.Dump()

    // test the circuit and read out the state vector
    circ.Run qs
    show "Dumping ket after run:"
    ket.Dump(showInd,0,false,true)
    // we are looking for the vector |b>|0>|0>|0> here
    show "Prob of measuring each qubit:"
    for q in qs do
        show "Q%d = %.12g" q.Id q.Prob1
    circ.RenderHT "RUS"
    show "Done"

[<LQD>]
// Does random tests to check that the EPR singlet state |01>-|10> is an eigenstate 
// of a unitary of test form U \otimes U^dagger, where U is any unitary in SU(2). 
let __EPRprojectionTest() =
    let ket         = Ket(2)
    let qs      = ket.Qubits
    let cv      = ket.Single()
    show "Dumping vector:"
    cv.Dump()
    show "Dumping ket before run:"
    //ket.Dump()
    ket.Dump(showInd,0,false,true)
    for q in qs do
        show "Q%d = %.12g" q.Id q.Prob1
        
    show "Now constructing Euler rotations..."    
    let MyU = EulerAngleSU2Rotation 3.0 2.5 -3.0
    let MyU' = EulerAngleSU2Rotation' 3.0 2.5 -3.0
    
    show "Now compiling circuit U \otimes U^dagger..."
    let ops (qs:Qubits) =
        //X !!(qs,0)
        X !!(qs,1) 
        H !!(qs,0)
        CNOT !!(qs,0,1)
        MyU' !!(qs,1)
        MyU !!(qs,0)
        CNOT !!(qs,0,1)
        H !!(qs,0)
        X !!(qs,1) 

    let circ  = Circuit.Compile ops qs
    circ.RenderHT("Ugates")
    circ.Dump()
    
    show "Now applying circuit to EPR pair"
    circ.Run qs
    show "Dumping ket after run:"
    ket.Dump(showInd,0,false,true)
    show "Prob of measuring each qubit:"
    for q in qs do
        show "Q%d = %.12g" q.Id q.Prob1
    circ.RenderHT "RUS"
    show "Done"

///////////////////////////////////////////////////////////////////////////////
//
// Shors algorithm: Experiments with approximate QFTs 
// Literature: [Coppersmith quant-ph/0201067], [Roetteler, Beth, AAECC 19(3): 177-193 (2008)]
//
///////////////////////////////////////////////////////////////////////////////
   
/// <summary> Run all steps of Shor's algorithm with approximate QFT. Minor modification of Dave's code. </summary>
/// <param name="doCirc">Doing Circuits or direct execution</param>
/// <param name="rslt">Bit results from each step</param>
/// <param name="N">Number to factor</param>
/// <param name="l">Prune level</param>
/// <param name="a">Classical value on left</param>
/// <param name="qs">Qubits to operate on</param>
let MyShorRun doCirc (rslt:int[]) (N:int) (level:int) (a:int) (qs:Qubits) =
    let n,m,xs,bs,anc   = getGroups qs

    let k               = qs.Head.Ket

    // Need a private ket for the async circuit compile
    let k'              = Ket(k.Count,Zero,k.Parts)
    let qs'             = List.mapi (fun i (q:Qubit) -> k'.[q.Id]) qs

    let swComp          = Diagnostics.Stopwatch()
    let swWrap          = Diagnostics.Stopwatch()
    let swRun           = Diagnostics.Stopwatch()

    let getUa verbose bitPos =
        let aExp        = modExp a (1<<<bitPos) N
        if doCirc then
            if verbose then show "         - Compiling circuit"
            swComp.Start()
            let circUa      = CompileUa N aExp qs'
            swComp.Stop()
            if verbose then
                show "%8f = mins for compile" swComp.Elapsed.TotalMinutes
                show "%8d = cnt of gates" (circUa.GateCount()*n*2)
                let hits,misses = Gate.CacheStats()
                show "%8d = cache hits" hits
                show "%8d = cache misses" misses
                show "%8d = compiled memory (MB)" (let ps = procStats(true) in ps.wsetMB)
                show "         - Wrapping circuit pieces"
            swWrap.Start()
            let gp          = GrowPars(30,(if verbose then 2 else 0),false) 
            let circUa      = circUa.GrowGates(k,gp)
            swWrap.Stop()
            if verbose then
                show "%8f = mins for growing gates" swWrap.Elapsed.TotalMinutes
                show "%8d = cnt of gates" (circUa.GateCount()*n*2)
                show "%8d = grown memory (MB)" (let ps = procStats(true) in ps.wsetMB)
            let runUa() =
                swRun.Start()
                circUa.Run qs
                swRun.Stop()
            runUa
        else
            let runUa() =
                swRun.Start()
                Ua N aExp qs
                swRun.Stop()
            runUa

    // Have to start with m measured
    M m 
    let L       = 2*n
    let PI2     = 2.0 * Math.PI
    
    // Compute accumulated phase
    let computePhiK j =
        
        if j=1 then 0.0
        else
            let rec doPhi k rtn =
                if (k > min j level) then PI2 * rtn
                // modify phase only if k <= prune level  
                elif rslt.[j-k] = 1 then
                    doPhi (k+1) (rtn+1.0/(pown 2. k)) else doPhi (k+1) rtn
            doPhi 2 0.0
    
    for j in 1..L do
        let bitPos  = L-j
        let runUa   = getUa (j=1) bitPos
        
        Reset Zero m
        H m
        if j>1 then EShift (computePhiK j) m
        runUa()
        H m
        M m
        let r           = m.Head.Bit.v
        rslt.[j-1]     <- r
 
        let ps      = procStats(true)
        show "        Bit: %3d [MB:%5d m=%d]" bitPos ps.wsetMB r
 
    show "%8f = mins for running" swRun.Elapsed.TotalMinutes

/// <summary> Top level call to Shor factoring algorithm with approximate QFT. Minor modification of Dave's code. </summary>
/// <param name="N">Number to factor (0=run tests and dump diagrams on components)</param>
/// <param name="doCirc">Compile to circuit (and optimize) or run native gates</param>
/// <param name="chkAbort">Routine to throw exception if abort desired (for long runs)</param>
let MyShor (N:int) (level:int) (doCirc:bool) (chkAbort:unit->unit) = 
    // experiments with pruned QFT. Prune level defined by integer l. 
    let notLegal (msg:string) =
        let sb  = StringBuilder()
        show "Legal numbers include:"
        for bits in 4..16 do
            let rec tryNum cnt N =
                if N >= (1<<<(bits-1)) && cnt < 10 then
                    if N % 2 <> 0 && false = isPrime N && false = isPrimePower N then
                        sprintf " %5d" N |> sb.Append |> ignore
                        tryNum (cnt+1) (N-1)
                    else tryNum cnt (N-1)
            sb.Length  <- 0
            tryNum 0 ((1<<<bits)-1)
            show "%2d bits: %O" bits sb
        ShorResult.Fail(msg,N)

    if N = 0 then
        testShorGates(doCirc)

        let N           = 15
        let a           = coPrime N
        let n           = CVec.Bits N
        let qs          = genQubits N
        show "   Drawing circuits"
        let rslt        = Array.create (2*N) 0
        let circUa      = CompileUa N a qs
        let cnt         = circUa.GateCount()*2*n
        let hits,misses = Gate.CacheStats()
        show "    Gate cache hits=%d misses=%d" hits misses
        show "    Total primitive gates: %d" cnt
        let circUa      = circUa.Fold()
        showLog "======== Dump of circuit ===="
        circUa.Dump()
        showLog "============================="
        for i in 0..3 do
            let nam = sprintf "Shor_%1d" i
            match i with
            | 0 -> circUa.RenderHT(nam,i)
            | 1 -> circUa.RenderHT(nam,i)
            | 2 -> circUa.RenderHT(nam,i,50.0)
            | _ -> circUa.RenderHT(nam,i,25.0)
        ShorResult.Fail("Test Mode",N)
    elif N % 2 = 0 then notLegal "Num must be odd"
    elif isPrime N then notLegal "Num must not be prime"
    elif isPrimePower N then notLegal "Num must not be a prime power"
    else
        // Do one round
        let n               = CVec.Bits N
        show "======== Doing Shor Round ========="
        show "%8d = N = Number to factor" N
        let a           = coPrime N
        show "%8d = a = coPrime of N" a
        show "%8d = n = number of bits for N" n
        show "%8d = 2^n" (1<<<n)
        let qs          = genQubits N
        let k           = qs.Head.Ket
        show "%8d = total qubits" qs.Length
        show "%8d = starting memory (MB)" (let ps = procStats(true) in ps.wsetMB)

        let good,total,pcnt     = getRndSucc N a
        show "%7.2f%% = prob of random result (%d/%d)" pcnt good total
        show "%7.2f%% = prob of Shor (worst case)" (100.0/(Math.Log(float n,2.)))
        let rslt        = Array.create (2*n) 0
        k.ChkAbort     <- chkAbort
        let sw          = Diagnostics.Stopwatch()
        sw.Start()
        MyShorRun doCirc rslt N level a qs
        sw.Stop()
        let elapsed     = sw.Elapsed.TotalSeconds
        show "%8g = Total Elapsed time (seconds)" elapsed

        let L           = 2*n-1
        let m           = Array.mapi (fun i bit -> bit <<< i) rslt |> Array.sum

        show "%8d = Max Entangled" k.MaxEntangled
        let permG,permS,permN   = k.Perms
        show "%8d = Gates Permuted" permG
        show "%8d = State Permuted" permS
        show "%8d = None  Permuted" permN
        let rslt    = chkResult N a (uint64 m) elapsed true
        show "CSV N a m den f1 f2 good,%d,%d,%d,%d,%d,%d,%d" N a m rslt.den rslt.f1 rslt.f2 (if rslt.msg = "" then 1 else 0)
        rslt

[<LQD>]
// Compute the diophantine reconstuction for a range of input rational numbers 
// to test if the best approximations for the given precision are found. 
let __CFTests() =
    let twoN    = 1UL <<< 8
    let Tarr    = Array.init 30 (fun index -> 137.0 + (double index))
    for i in [0..4] do
        let Tarr3   = Tarr 
                      |> Array.map (fun x -> (rat (x/float twoN) ((1UL <<< i) * twoN))) 
                      |> Array.map (fun (x, y) -> printf "(%d, %d), " x y; x, y)
                      |> Array.map (fun (x, y) -> (float x)/(float y))
        printf "\n"
        Array.iter (fun x -> printf "%f, " x) Tarr3
        printf "\n"
    show "Done continued fraction tests"

[<LQD>]
// Compute phase corrections for the approximate QFT and test if they 
// match the given prune level. 
let __phaseTest() = 
    let rslt     = [| 1; 0; 1; 0; 1; 1; |] 
    let PI2      = 2.0 * Math.PI   
    let PrnLev   = 2 // prune level for QFT phases
 
    let computePhiK j prnLev =
        if j=1 then 0.0
        else
            let rec doPhi k rtn =
                if (k > min j (prnLev+1)) then PI2 * rtn 
                    // modify phase only if |j-k| <= prnKev 
                elif rslt.[j-k] = 1 then
                    doPhi (k+1) (rtn+1.0/(pown 2. k)) else doPhi (k+1) rtn
            doPhi 2 0.0 
    printf "measured bits: %A\n" rslt
    printf "Prune level = %i\n" PrnLev
    for i in [1..6] do
        printf "Total phase of bit %i is: %f\n" i (computePhiK i PrnLev)
    show "Done phase test"
    
[<LQD>]
// Execute <runs> many iterations of Shor's algorithm to factor the integer <N> where 
// the phases of the underlying QFT are pruned at <level>.
let __ShorApprox N (runs:int) (level:int) = // do n runs of Shor for N for given prune level
    for i in [1..runs] do
        let rslt = MyShor N level true (fun () -> ())
        let msg     = if rslt.msg = "" then "SUCCESS!!" else "FAILED: " + rslt.msg
        let str     = 
            sprintf "GOT:%5d=%4dx%4d co=%5d n,q=%2d,%2d mins=%.2f %s"
                rslt.N rslt.f1 rslt.f2 rslt.coPrime rslt.n rslt.qCnt 
                    (rslt.elapsed/60.0) msg
        show "%s" str
         

///////////////////////////////////////////////////////////////////////////
//
// Shor's algorithm for discrete logs over prime fields, October 2014
// Literature: [Proos, Zalka, arxiv:quant-ph/0301141], Section 4.1
//
///////////////////////////////////////////////////////////////////////////
[<LQD>]
/// <summary> Simulate Shor's algorithm for discrete logarithms modulo a prime </summary>
/// <param name="p">Prime</param>
/// <param name="u">Generator of cylic group</param>
/// <param name="v">Element to compute dlog of</param>
let __ShorDlogExact (p:int) (u:int) (v:int) =
    
    // Prepare the initial state; QFT of length p method
    let n       = CVec.Bits p // now 2^(n-1) <= p <= 2^n
    let ket     = Ket(2*n) // n qubits for measurement and n qubits for state 
    let qs      = ket.Qubits
    let cv      = ket.Single()
    let rslt    = Array.create (2*n) 0 // will be used to store final result
    X !!(qs,qs.Length-1) // prepare the initial state = neutral group element
             
    // Assemble the unitary gates needed for the circuit
    let GatePerm (a:int[]) (qs:Qubits) = 
        let gate (a:int[]) =  
            Gate(
                Qubits  = qs.Length,
                Name    = sprintf "Perm(%A)" a,
                Help    = sprintf "Permutation matrix corresponding to %A" a,
                Mat     = (
                    let np2     = double a.Length |> log |> fun x -> x/(log 2.0) |> ceil |> int |> pown 2 
                    // pick up next power of two np2 that is greater or equal qs.Length
                    let mat     = CSMat(np2,true)
                    for i in 0..a.Length-1 do   // permutation matrix corresponding to array a
                        mat.r(a.[i],i) <- 1.0
                    for i in a.Length..np2-1 do // fill up to size np2 with identity matrix
                        mat.r(i,i) <- 1.0
                    mat),
                Draw    = "\\gate{\pi_L}"
                )
        (gate a).Run qs
    
    // QFT for any length <a>. Naive implementation from the definition of the DFT matrix 
    let QFTArbLen (k:int) (qs:Qubits) = 
        let len = pown 2 qs.Length
        if k > len then 
            failwithf "Length %i of the QFT must be less or equal to number of qubits %i" k qs.Length   
        let gate (k:int) =  
            Gate(
                Qubits  = qs.Length,
                Name    = sprintf "QFTArbLen(%A)" k,
                Help    = sprintf "QFT of length %A, implemented directly" k,
                Mat     = (
                    // pick up next power of two np2 that is greater or equal qs.Length
                    let PI2k = 2.0 * Math.PI / (float k)
                    let mat     = CSMat(len,true)
                    for i in 0..k-1 do   // construct the DFT of size kxk directly
                        for j in 0..k-1 do 
                            mat.r(i,j) <- Math.Cos (-PI2k * (float i) * (float j))
                            mat.i(i,j) <- Math.Sin (-PI2k * (float i) * (float j))
                    for i in k..len-1 do // fill up to size np2 with identity matrix
                        mat.r(i,i) <- 1.0
                    mat),
                Draw    = "\\gate{DFT_k}"
                )
        (gate k).Run qs
    
    // implement the group shifts. NOTE: modify this part for different arithmetic, e.g. ECC
    let Ua (p:int) (a:int) = 
        Array.init p (fun i -> a*i % p) |> GatePerm

    // phase estimation for the first register
    List.iter (fun i -> M !!(qs,i); Reset Zero !!(qs,i)) [0..n-1] // initialize measurement qubits 
    QFTArbLen (p-1) !!(qs,[0..n-1]) // prepares input state on measurement register
    for j in 0..n-1 do 
        let aExp   = modExp u (1 <<< j) p // powers of generator u
        //printf "current bit = %i current modulus = %i\n" j aExp
        Cgate (fun qs -> Ua p aExp qs) !!(qs,[n-j-1] @ [n..(2*n-1)])
    QFTArbLen (p-1) !!(qs,[0..n-1])
    
    // measure the first register 
    List.iter (fun i -> M !!(qs,i)) [0..n-1] 
    for i in 0..n-1 do 
        rslt.[i] <- qs.[i].Bit.v
    
    // phase estimation for the second register
    List.iter (fun i -> Reset Zero !!(qs,i)) [0..n-1] // initialize measurement qubits
    QFTArbLen (p-1) !!(qs,[0..n-1]) // prepares input state on measurement register
    for j in 0..n-1 do 
        let aExp    = modExp v (1 <<< j) p // powers of generator v
        //printf "current bit = %i current modulus = %i\n" j aExp
        Cgate (fun qs -> Ua p aExp qs) !!(qs,[n-j-1] @ [n..(2*n-1)])
    QFTArbLen (p-1) !!(qs,[0..n-1])
    
    // measure the second register 
    List.iter (fun i -> M !!(qs,i)) [0..n-1] 
    for i in 0..n-1 do
        rslt.[n+i] <- qs.[i].Bit.v
    
    let ps      = procStats(true)
    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt
    
    // Check if result was correct dlog
    let extendedEuclid (a:int) (b:int) =
        let rec iterDivQuo r2 s2 t2 r1 s1 t1 = 
            if r1=0 then (r2, s2) // just return gcd and inverse of a 
            else
                let q = r2 / r1
                let r = r2 - q*r1
                let s = s2 - q*s1
                let t = t2 - q*t1
                iterDivQuo r1 s1 t1 r s t 
        iterDivQuo a 1 0 b 0 1
        
    let modulo m n = ((m % n) + n) % n // make sure result of mod n is in [0..n-1]
    
    let checkResult =  
        let m1 = List.map (fun i -> rslt.[n-1-i] <<< i)   [0..n-1] |> List.fold (+) 0
        let m2 = List.map (fun i -> rslt.[2*n-1-i] <<< i) [0..n-1] |> List.fold (+) 0 
        let inv = extendedEuclid m1 (p-1)
        if (fst inv) <> 1 then 
            show "FAIL! measured modulus is not coprime with (p-1)"
            show "RES: Shor u v p dlog %i %i %i %i %i" u v p 0 0
        else 
            let a = modulo ((snd inv) * m2) (p-1)
            show "SUCCESS! The dlog of %i with respect to %i modulo %i is %i" v u p a
            show "Check: %i = %i^%i mod %i: %b" v u a p (v=(modExp u a p))
            show "RES: Shor u v p dlog %i %i %i %i %i" u v p a 1    
    checkResult
    show "Shor dlog done"

[<LQD>]
let __ShorDlogRuns (runs:int) p u v = // do n runs of Shor dlog for N for given prune level
    for i in [1..runs] do
        __ShorDlogExact p u v 
    ()

    
/////////////////////////////////////////////////////////////////////////////
// 
// Implementations of Euclid for modular inversion
//
/////////////////////////////////////////////////////////////////////////////

// A simple implementation of Euclid vs extended binary Euclid
    
let rec gcd u v =
    match v with | 0 -> u | _ -> gcd v (u%v) 

let modinvRegular p x = 
    let rec xgcd u v r s = 
        printf "%A %A %A %A\n" u v r s
        match u, v, r, s with 
        | _, 0, _, _ -> r, u
        | _, _, _, _ -> let quo = u / v
                        let rem = u % v
                        xgcd v rem s (r-quo*s) 
    xgcd p x 0 1

let modinvBinary p x = 
    let rec xgcd u v b d =
        printf "%A %A %A %A\n" u v b d
        match u, v, b, d with 
        | 0, _, _, _                     -> d, v
        | u, _, b, _ when u%2=0 && b%2=0 -> xgcd (u >>> 1) v (b >>> 1) d
        | u, _, b, _ when u%2=0          -> xgcd (u >>> 1) v ((b-p) >>> 1) d 
        | _, v, _, d when v%2=0 && d%2=0 -> xgcd u (v >>> 1) b (d >>> 1)
        | _, v, _, d when v%2=0          -> xgcd u (v >>> 1) b ((d-p) >>> 1)
        | u, v, _, _ when u >= v         -> xgcd (u-v) v (b-d) d
        | _, _, _, _                     -> xgcd u (v-u) b (d-b)
    xgcd p x 0 1

let gcdcall () = 
    let x = 1009
    let y = 31
    let z = gcd x y
    printf "the gcd of %A and %A is %A\n" x y z

[<LQD>]
let __xgcdcall () = 
    let x = 1009
    for i in 800..802 do 
        let z, g = modinvBinary x i 
        let z1, g1 = modinvRegular x i 
        printf "Binary: the inverse of %A mod %A is %A and the gcd is %A\n" i x (((z%x)+x)%x) g
        printf "Regular: the inverse of %A mod %A is %A and the gcd is %A\n" i x (((z1%x)+x)%x) g1

let maxpower (x:int) = 
    let rec powerRed x k = 
        match x % 2 with 
        | 1 -> k 
        | _ -> powerRed (x >>> 1) (k+1) 
    powerRed x 0 

[<LQD>]
let __powerTest (x:int) = 
    printf "%A \t %A \n" x (maxpower x)
    show "done"

let predictorGreedy (r:int) (s:int) (u:int) (v:int) = 
    let mutable value = 0

    if r % 2 = 0 then 
        if ((r >>> 1) % 2 = 0) then 
            if (2*s > r) && (2*v >= u) then value <- 4 else value <- 2            
        else
            value <- 4
    if s % 2 = 0 then 
        if (s >>> 1) % 2 = 0 then 
            if (2*r > s) then 
                if (u % 2 = 1) && (2*v <= u) then value <- 1 else value <- 3
            else 
                value <- 1
        else 
            value <- 3
    value

let predictorSlow (r:int) (s:int) = 
    let mutable value = 0

    match r, s with 
    | r, _ when r % 2 = 0 -> 2
    | _, s when s % 2 = 0 -> if (r > s) then 3 else 1
    | r, s when r >= s -> 3
    | _, _ -> 4
    
let predictor (r:int) (s:int) (u:int) (v:int) (a:int) (b:int) (p:int) (round:int) = 
    let mutable value = 0

    if r % 2 = 0 then 
        if ((r >>> 1) % 2 = 0) then 
            if ((abs(r-s)) > s) then value <- 4 else value <- 2
        else value <- 4
    if s % 2 = 0 then 
        if (s >>> 1) % 2 = 0 then 
            if (2*r < s) then value <- 1
            else
                if ((pown 2 round) > p) then value <- 1 else value <- 3 //if (b < a) then value <- 1 else value <- 3
        else value <- 3
    value

let mginvGreedy p x = 
    let rec xmg u v r s k case =
        printf "%A \t %A \t %A \t %A \t %A \t %A \t rnd: %A \t fwd: %A\n " p x u v r s k case
        //if (case <> (predictor r s u)) && ((case = 1) || (case = 3)) then printf "Fail! %A %A\n" p x
        match u, v, r, s with 
        | _, 0, r, _                     -> r
        | u, _, _, _ when u%2=0          -> let pow = (maxpower u) 
                                            xmg (u >>> pow) v r (s <<< pow) (k+pow) 1
        | _, v, _, _ when v%2=0          -> let pow = (maxpower v) 
                                            xmg u (v >>> pow) (r <<< pow) s (k+pow) 2
        | u, v, _, d when u > v          -> xmg ((u-v) >>> 1) v (r+s) (s <<< 1) (k+1) 3
        | _, _, _, _                     -> xmg u ((v-u) >>> 1) (r <<< 1) (r+s) (k+1) 4
    xmg p x 0 1 0 0

let mginvComplete p x = 
    let rec xmg u v r s a b k case =
        printf "%A \t %A \t %A \t %A \t %A \t %A \t %A \t %A \t rnd: %A \t pred: %A \t case: %A \n" p x u v r s a b k (predictor r s u v a b p k) case 
        //printf "%A \t %A \t %A \t %A \t %A \t %A \t val1: %A \t val2: %A \t rnd: %A \t fwd: %A \n" p x u v r s (x*r+u*(pown 2 k)) (x*s-v*(pown 2 k)) k  case 
        //printf "%A \t %A \t %A \t %A \t %A \t %A \t rnd: %A \t fwd: %A pred: %A\n" p x u v r s k case (predictor r s u v)
        //if (case <> (predictor r s u v a b)) && ((case = 1) || (case = 3)) then printf "Fail! prime %A u %A v %A r %A s %A round %A case %A input %A\n" p u v r s k case x
        match u, v, r, s with 
        | _, 0, r, _                     -> r
        | u, _, _, _ when u%2=0          -> xmg (u >>> 1) v r (s <<< 1) a (b <<< 1) (k+1) 1
        | _, v, _, _ when v%2=0          -> xmg u (v >>> 1) (r <<< 1) s (a <<< 1) b (k+1) 2
        | u, v, _, d when u > v          -> xmg ((u-v) >>> 1) v (r+s) (s <<< 1) (a+b) (b <<< 1) (k+1) 3
        | _, _, _, _                     -> xmg u ((v-u) >>> 1) (r <<< 1) (r+s) (a <<< 1) (a+b) (k+1) 4
    xmg p x 0 1 1 0 0 0

let mginvFull p x = 
    let rec xmg u v r s k case =
        printf "%A \t %A \t %A \t %A \t %A \t %A \t rnd: %A \t fwd: %A \t pred: %A \n" p x u v r s k (predictor r s u v) case 
        //printf "%A \t %A \t %A \t %A \t %A \t %A \t val1: %A \t val2: %A \t rnd: %A \t fwd: %A \n" p x u v r s (x*r+u*(pown 2 k)) (x*s-v*(pown 2 k)) k  case 
        //printf "%A \t %A \t %A \t %A \t %A \t %A \t rnd: %A \t fwd: %A pred: %A\n" p x u v r s k case (predictor r s u v)
        //if (case <> (predictor r s u v)) && ((case = 1) || (case = 3)) then printf "Fail! %A %A\n" p x
        match u, v, r, s with 
        | _, 0, r, _                     -> r
        | u, _, _, _ when u%2=0          -> xmg (u >>> 1) v r (s <<< 1) (k+1) 1
        | _, v, _, _ when v%2=0          -> xmg u (v >>> 1) (r <<< 1) s (k+1) 2
        | u, v, _, d when u > v          -> xmg ((u-v) >>> 1) v (r+s) (s <<< 1) (k+1) 3
        | _, _, _, _                     -> xmg u ((v-u) >>> 1) (r <<< 1) (r+s) (k+1) 4
    xmg p x 0 1 0 0

let mginvSlow p x = 
    let rec xmg u v r s k case =
        match u, v, r, s with 
        | _, 0, r, _                     -> r
        | u, _, _, _ when u%2=0          -> xmg (u >>> 1) v r (s <<< 1) (k+1) 1
        | _, v, _, _ when v%2=0          -> xmg u (v >>> 1) (r <<< 1) s (k+1) 2
        | u, v, _, d when u > v          -> xmg (u-v) v (r+s) s (k+1) 3
        | _, _, _, _                     -> xmg u (v-u) r (r+s) (k+1) 4
    xmg p x 0 1 0 0

let mginv p x = 
    let rec xmg u v r s k case =
        match u, v, r, s with 
        | _, 0, r, _                     -> r
        | u, _, _, _ when u%2=0          -> xmg (u >>> 1) v r (s <<< 1) (k+1) 1
        | _, v, _, _ when v%2=0          -> xmg u (v >>> 1) (r <<< 1) s (k+1) 2
        | u, v, _, d when u > v          -> xmg ((u-v) >>> 1) v (r+s) (s <<< 1) (k+1) 3
        | _, _, _, _                     -> xmg u ((v-u) >>> 1) (r <<< 1) (r+s) (k+1) 4
    xmg p x 0 1 0 0

[<LQD>]
let __mgcall p x = 
    let z = mginv p x
    printf "the MG inverse of %A and %A is %A\n" p x z


[<LQD>]
let __mgcallLoop p = 
//    for p in 3..2..100000 do 
        for x in 1..(p-1) do 
            printf "Input: %A\n" x
            printf "=================================================================\n"
            mginvComplete p x |> ignore 
            printf "\n\n"
        //printf "the MG inverse of %A and %A is %A\n" p x z

/////////////////////////////////////////////////////////////////////////////
// 
// Implementations of Cordics for trigonometric functions and their inverses 
//
/////////////////////////////////////////////////////////////////////////////

/// Numerical experiments with Cordics for trig and inverse trig functions. 
/// In particular studying the behavior under truncated threshold decisions.

/// cordic (x:double) (nmax:int) (prec:int)
///     computes the arcsin(<x>) using the cordic method. The total number of rounds is bounded by <nmax>. 
///     Thresholds are computed with precision <prec> which also limits the overall achievable precision.
let cordic (x:double) (nmax:int) (prec:int) = 
    let myRound (x:double) (k:int) = 
        x |> (fun i -> (pown 10.0 k) * i) |> Math.Round |> (fun i -> i/(pown 10.0 k))
    
    let pi = Math.PI
    let tanPrecomp i = Math.Atan (pown 2.0 (-i))  

    let dynamicGain k = 
        [0..k] 
        |> Seq.map (fun i -> Math.Sqrt (1.0 + (pown 2.0 (-2*i))))
        |> Seq.fold (*) 1.0
            
    let gain = 
        [0..nmax] 
        |> Seq.map (fun i -> Math.Sqrt (1.0 + (pown 2.0 (-2*i))))
        |> Seq.fold (*) 1.0

    let thresholdCordic  = function z when (z < 0.0) -> -1.0 | _ -> 1.0
    //let thresholdArcSin  = function z when ((myRound z prec) < (myRound (gain * x)) prec) -> 1.0 | _ -> -1.0 
    let thresholdArcSin (i:int) = function z when ((myRound z prec) < (myRound ((dynamicGain (i+1))* x)) prec) -> 1.0 | _ -> -1.0 
    
    let rec cordicRec x y z i =
        match x, y, z, i with 
        | _, _, _, i when (i = nmax) -> (x, y, z)
        | x, y, z, _ -> 
            let d = thresholdArcSin i y
            let x' = x - y * d * (pown 2.0 (-i))
            let y' = y + x * d * (pown 2.0 (-i))
            let z' = z - d * (tanPrecomp i)
            cordicRec x' y' z' (i+1)
        | _, _, _, _ -> failwith "this should never happen"
        
    let _, _, arcsinRaw = (cordicRec 1.0 0.0 0.0 0)
    -arcsinRaw
    
let decay (x:double) (i:int) (prec:int) = 
    Math.Abs(Math.Sin(cordic x i prec) - x)

let plotDecayNormal (x:double) (n:int) (prec:int) = 
    let L = [| for i in 1..n do yield decay x i prec|]
    pltT.line(y=L,markerSize=10,color=Color.Blue,markerStyle=MarkerStyle.Star5).display()

let plotDecayLogscale (x:double) (n:int) (prec:int) = 
    let L = [| for i in 1..n do yield decay x i prec|] |> Array.map Math.Log10
    pltT.line(y=L,markerSize=10,color=Color.Blue,markerStyle=MarkerStyle.Star5).display()
        
[<LQD>]
let __plotDecay x n prec = 
    plotDecayLogscale x n prec
    
        
///////////////////////////////////////////////////////////////////////////
//
// Quantum cordic circuit to compute trig functions and their inverses
//
///////////////////////////////////////////////////////////////////////////


// big construction site here...
//
//let arcsinCircuit (x:double) (nmax:int) (prec:int) = 
//    let myRound (x:double) (k:int) = 
//        x |> (fun i -> (pown 10.0 k) * i) |> Math.Round |> (fun i -> i/(pown 10.0 k))
//    
//    let pi = Math.PI
//    let tanPrecomp i = Math.Atan (pown 2.0 (-i))  
//    let gain = 
//        [0..nmax] 
//        |> Seq.map (fun i -> Math.Sqrt (1.0 + (pown 2.0 (-2*i))))
//        |> Seq.fold (*) 1.0
//
//    let thresholdCordic  = function z when (z < 0.0) -> -1.0 | _ -> 1.0
//    let thresholdArcSin  = function z when ((myRound z prec) < (myRound (gain * x)) prec) -> 1.0 | _ -> -1.0 ;
//    
//    let rec cordicRec x y z i =
//        match x, y, z, i with 
//        | _, _, _, i when (i = nmax) -> (x, y, z)
//        | x, y, z, _ -> 
//            let d = thresholdArcSin y
//            let x' = x - y * d * (pown 2.0 (-i))
//            let y' = y + x * d * (pown 2.0 (-i))
//            let z' = z - d * (tanPrecomp i)
//            cordicRec x' y' z' (i+1)
//        | _, _, _, _ -> failwith "this should never happen"
//        
//    let _, _, arcsinRaw = (cordicRec 1.0 0.0 0.0 0)
//    -arcsinRaw
//    
//
//[<LQD>]
//let __TestCordic (name:string) (n:int) (x:int) = 
//    // Tools to manipulate bools
//    let BoolInt n x = 
//        [ for i in [0..(n-1)] do yield ((x >>> i) % 2) ] 
//    let IntBool n (x: int array) = 
//        List.fold (+) 0 <| List.map (fun i -> (1 <<< i) * x.[i]) [0..x.Length-1]
//    let PrepBool (qs:Qubits) offset step ls = 
//        List.mapi (fun i x -> if x <> 0 then X !!(qs,(offset+i*step))) ls 
//
//    // Prepare the initial state
//    let arrangeCordicInputs =   
//        match name with 
//        | "ArcSin" | "ArcCos" | "Sin" | "Cos" -> 
//            let k = Ket(n*n)
//            let qs = k.Qubits
//            x |> BoolInt n |> PrepBool qs 0 1 |> ignore
//            qs
//        | _ -> failwith "Unknown trig function."
//                      
//    let qs = arrangeCordicInputs          // arranges inputs in pattern needed for the multiplier
//    
//    let arrangeCordicCircuit = 
//        match name with 
//        | "ArcSin" -> 
//            let CIRC = arcsinCordic         // constructing the multiplier circuit
//            let rslt = Array.create (n*n) 0 // will be used to store final result 
//            CIRC, rslt
//        | _ -> failwith "This trig function is not implemented yet."
//
//    let CIRC, rslt = arrangeCordicCircuit
//    CIRC qs 
//    
//    // measure the register 
//    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1] 
//    for i in 0..(qs.Length-1) do 
//        rslt.[i] <- qs.[i].Bit.v
//    
//    let ps      = procStats(true)
//    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt
//
//    let checkResult =  
//        let sInt = rslt.[n..(2*n)] |> Array.rev |> IntBool n // have to reverse order      
//        show "RES: Cordic for %s" name
//        show ": n x rslt %i %i %i" n x sInt
//
//    checkResult
//    show "Check cordic done"

   
//[<LQD>]
//let __TestAdder (n:int) (s1:int) (s2:int) =
//    // Tools to manipulate bools
//    let BoolInt n x = 
//        [ for i in [0..(n-1)] do yield ((x >>> i) % 2) ] 
//    let IntBool n (x: int array) = 
//        List.fold (+) 0 <| List.map (fun i -> (1 <<< i) * x.[i]) [0..x.Length-1]
//    let makeBit = function 0 -> Zero | 1 -> One | _ -> Unknown 
//    let QubitsInt (a:int) =
//        let ket = Ket 0
//        List.iter (function 0 -> ignore (ket.Add Qubit.Zero) | _ -> ignore (ket.Add Qubit.One)) (BoolInt n a)
//        ket.Qubits
//    // Prepare the initial state; QFT of length p method
//    let l1 = BoolInt n s1 |> List.map makeBit
//    let l2 = BoolInt n s2 |> List.map makeBit
//    let S1 = List.map (fun x -> Ket(1,x)) l1 // n qubits for first number 
//    let S2 = List.map (fun x -> Ket(1,x)) l2 // n qubits for second number 
//    let xs = List.map (fun (x:Ket) -> x.Qubits) S1 |> List.fold (@) []
//    let ys = List.map (fun (x:Ket) -> x.Qubits) S2 |> List.fold (@) []   
//    
//    let rs = QubitsInt s1
//    let ss = QubitsInt s2
//    let ket    = Ket(2) // 1 qubit for final carry
//    let z      = ket.Qubits
//    //let qs = rs @ ss @ z                 // constructing the inputs x, y, and 1 additional bit for final carry
//    
//    let rec formatCuccaro (qs1:Qubits) (qs2:Qubits) (input:Ket) = 
//        //ignore (input.Add Qubit.Zero)
//        match qs1, qs2 with
//        | [], [] -> input.Qubits
//        | q::tailA, qs2 ->
//            ignore (input.Add q.State)
//            //ignore (input.Add Qubit.Zero)
//            formatCuccaro tailA qs2 input
//        | [], q::tailB ->
//            //ignore (input.Add Qubit.Zero)
//            ignore (input.Add q.State)
//            formatCuccaro [] tailB input
////        | qA::tailA, qB::tailB ->
////            ignore (input.Add qA.State)
////            ignore (input.Add qB.State)
////            formatCuccaro tailA tailB input
////    
//    
//    
////    let qs = formatCuccaro rs ss ket 
//    
//
//    let k = Ket(2*n+2)
//    let qs = k.Qubits
//        
//    for i in [0..(n-1)] do
//        if l1.[i] = One then X !!(qs,i)
//    for i in [0..(n-1)] do
//        if l2.[i] = One then X !!(qs,i)
//
//    for i in 0..qs.Length-1 do 
//        qs.[i].Dump()
//        show "qs[%2d] = %s" i (qs.[i].ToString ())
//
//    let ADD = CuccaroAdder               // constructing the adder circuit
//    let rslt    = Array.create (2*n+2) 0 // will be used to store final result
//    
//    let circ = Circuit.Compile ADD qs    // compiler the circuit for the given input
//    circ.Run qs                          // run the circuit
//    
//    for i in 0..qs.Length-1 do 
//        qs.[i].Dump()
//        show "qs[%2d] = %s" i (qs.[i].ToString ())
//
//    // measure the register 
//    List.iter (fun i -> M !!(qs,i)) [0..2*n+1] 
//    for i in 0..(qs.Length-1) do 
//        rslt.[i] <- qs.[i].Bit.v
//    
//    let ps      = procStats(true)
//    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt
//
//    let checkResult =  
//        let sInt = rslt.[2..(n+1)] |> Array.rev |> IntBool n // have to reverse order      
////        if sInt <> s then 
////            show "FAIL! measured bitvector is different from hidden shift"
//        show "RES: Adder n s1 s2 rslt %i %i %i %i" n s1 s2 sInt
////        else 
////            show "SUCCESS! The hidden shift is s=%i"   sInt
////            show "RES: HShift n s rslt %i %i %i 1" n s sInt
//    
//    checkResult
//    show "Check adder done"



///////////////////////////////////////////////////////////////////////////
//
// Some circuits
//
///////////////////////////////////////////////////////////////////////////

let BuildReciprocal (qs:Qubits) = 
    let tx = BuildMultiplyControlledNOT
    
    let in0 = qs.[0] 
    let in1 = qs.[1]
    let in2 = qs.[2] 
    let in3 = qs.[3]
    let in4 = qs.[4]
    let out0 = qs.[5]
    let out1 = qs.[6]
    let out2 = qs.[7]
    let out3 = qs.[8]
    let out4 = qs.[9]

    X [in0]; X [in1]; X [in2]; X [in3]
    tx  [in4; in0; in1; in2; in3] [out0] [out4]
    X [in0]; X [in1]; X [in2]; X [in3]
    X [in3]; X [in4]; 
    tx  [in1; in3; in4] [out0] [out4]
    X [in3]; X [in4]; 
    X [in4]    
    tx  [in4] [out0] [out4]
    X [in4]
    tx  [out0] [out1] [out4]
    X [in0]; X [in1]; X [in2]; X [in4]
    tx  [in3; in0; in1; in2; in4] [out0] [out4]
    X [in0]; X [in1]; X [in2]; X [in4]
    tx  [out0] [out1] [out4]
    tx  [out0] [out1] [out4]
    tx  [out0] [out2] [out4]
    X [in0]; X [in1]; X [in3]; X [in4]
    tx  [in0; in1; in3; in4] [out0] [out4]
    X [in0]; X [in1]; X [in3]; X [in4]
    tx  [out0] [out1] [out4]
    tx  [out0] [out2] [out4]
    tx  [out0] [out2] [out4]
    X [in2]; X [in3]; X [in4]; 
    tx  [in0; in1; in2; in3; in4] [out0] [out4]
    X [in2]; X [in3]; X [in4]; 
    tx  [out0] [out2] [out4]
    tx  [out0] [out4] [out1]
    X [in1]; X [in2]; X [in3]; X [in4]; 
    tx  [in1; in2; in3; in4] [out0] [out4]
    X [in1]; X [in2]; X [in3]; X [in4]; 
    tx  [out0] [out4] [out1]
    X [in3]; X [in4]; 
    tx  [in2; in3; in4] [out1] [out4]
    X [in3]; X [in4];
    X [in0]; X [in2]; X [in3]; X [in4]; 
    tx  [in0; in2; in3; in4] [out3] [out4]
    X [in0]; X [in2]; X [in3]; X [in4]; 

[<LQD>]
let __size5 (n:int) (x:int) = 
    let b = BoolInt n (bigint x)
    let initialState = Array.zeroCreate (2*n)
    for i in 0..(n-1) do 
        initialState.[i] <- b.[i]

    let k = Ket((2*n))
    let qs = k.Qubits
    let cc = Circuit.Compile BuildReciprocal qs
    let c1 = cc.Fold(true)
    c1.RenderHT("mycirc")
    let circ = Circuit.Compile BuildReciprocal qs |> MyCircuitExport
    let finalState = MyCircuitSimulateFast circ initialState
    show "%A" finalState


//////////////////////////////////////////////////////////////////////////
//
// Computing circuit Boolean expressions from Toffoli network
//
//////////////////////////////////////////////////////////////////////////

type TofGate =
  | MyNOT of int
  | MyCNOT of int * int
  | MyTOFF of int * int * int

type BoolExp =
  | BFalse
  | BVar of int
  | BNot of BoolExp
  | BAnd of BoolExp * BoolExp
  | BXor of BoolExp * BoolExp

let rec prettyPrintBexp bexp = 
  match bexp with
  | BFalse -> "false"
  | BVar i -> sprintf "%d" i
  | BNot x -> "~" + (prettyPrintBexp x)
  | BAnd (x, y) -> "(" + (prettyPrintBexp x) + " && " + (prettyPrintBexp y) + ")"
  | BXor (x, y) -> "(" + (prettyPrintBexp x) + " <> " + (prettyPrintBexp y) + ")"

let prettyPrintExps m = 
  let printExp k bexp = printf "%d -> %s\n" k (prettyPrintBexp bexp)
  Map.iter printExp m

let collectDep lst =
  let findDef k m = 
    match Map.tryFind k m with
    | None   -> Set.singleton k
    | Some s -> s
  let applyGate deps = function
    | MyNOT _ -> deps
    | MyCNOT (c, t) -> Map.add t (Set.union (findDef c deps) (findDef t deps)) deps
    | MyTOFF (c1, c2, t) ->
        let c1deps = findDef c1 deps
        let c2deps = findDef c2 deps
        let tdeps  = findDef t deps
        Map.add t (Set.unionMany [c1deps; c2deps; tdeps]) deps
  List.fold applyGate Map.empty lst

let collectExp lst =
  let findDef k m = 
    match Map.tryFind k m with
    | None   -> BVar k
    | Some s -> s
  let applyGate exps = function
    | MyNOT _ -> exps
    | MyCNOT (c, t) -> Map.add t (BXor (findDef c exps, findDef t exps)) exps
    | MyTOFF (c1, c2, t) ->
        let c1exp = findDef c1 exps
        let c2exp = findDef c2 exps
        let texp  = findDef t exps
        Map.add t (BXor (texp, BAnd (c1exp, c2exp))) exps
  List.fold applyGate Map.empty lst

let truncate i lst =
  let rec trunc acc i lst = 
    if i <= 0 then List.rev acc else 
      match lst with
      | [] -> List.rev acc
      | x::xs -> trunc (x::acc) (i-1) xs
  trunc [] i lst

let test = 
  printf "%A\n" (collectDep <| truncate 2 [MyNOT 1; MyTOFF (3, 2, 1); MyCNOT (1, 2)])
  prettyPrintExps (collectExp <| truncate 2 [MyNOT 1; MyTOFF (3, 2, 1); MyCNOT (1, 2)])


//////////////////////////////////////////////////////////////////////////
//
// Bug reports 
//
//////////////////////////////////////////////////////////////////////////

#if FALSE
[<LQD>]
let __FoldBug () = 
    let k = Ket(6)
    let qs = k.Qubits

    let toffnetwork (qs:Qubits) = 
        CNOT [qs.[0]; qs.[1]] 
        CCNOT [qs.[0]; qs.[1]; qs.[2]]
        CCNOT [qs.[3]; qs.[4]; qs.[5]]

    let circuit = Circuit.Compile toffnetwork qs
    let circuitfold = circuit.Fold(true)

    let depth = circuit.GateCount(true, fun x -> x.Name = "CCNOT")
    let depthfold = circuitfold.GateCount(true, fun x -> x.Name = "CCNOT")
    show "Depth (original circuit) is %6d" depth
    show "Depth (folded circuit) is %6d" 
#endif

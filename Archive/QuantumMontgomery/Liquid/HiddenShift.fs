module Microsoft.Research.Liquid.HiddenShift

open System
//open System.Collections.Generic
//open System.IO
//open System.Text
//open System.Text.RegularExpressions
//open System.Drawing
//open System.Windows.Forms
//open System.Windows.Forms.DataVisualization.Charting
//open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

open GaussSums

///////////////////////////////////////////////////////////////////////////
//
// Quantum algorithm for shifted bent functions
//
///////////////////////////////////////////////////////////////////////////

[<LQD>]
/// <summary> Simulate the Boolean hidden shift algorithm </summary>
/// <param name="n">Number of qubits</param>
/// <param name="s">Hidden shift, encoded as integer</param>
let __BooleanHiddenShift (n:int) (s:int) =
    
    // Tools to manipulate bools
    let BoolInt n x = 
        [ for i in [0..(n-1)] do yield ((x >>> i) % 2) ] |> List.toArray
    let IntBool n (x: int array) = 
        List.fold (+) 0 <| List.map (fun i -> (1 <<< i) * x.[i]) [0..x.Length-1]
    let Xor (a:int) (b:int) = 
        (a+b-2*a*b) % 2
    let PermBool n k i = 
        (i + k) % (pown 2 n)
            
    // Implementations of some known families of bent functions
    // Class 1: The inner product function (ip)
    let ip (a: int array) = 
        List.map (fun i -> a.[2*i]*a.[2*i+1]) [0..(a.Length)/2-1] 
        |> List.fold (Xor) 0
    // Class 2: Maiorana-McFarland bent functions (mmf)
    let mmf permEven permOdd (a: int array)  = 
        let len = (a.Length)/2
        let aEven = [| for i in 0..len-1 do yield a.[2*i] |]  
                    |> IntBool len |> permEven |> BoolInt len
        let aOdd  = [| for i in 0..len-1 do yield a.[2*i+1] |]
                    |> IntBool len |> permOdd  |> BoolInt len
        List.map (fun i -> aEven.[i] * aOdd.[i]) [0..len-1] 
        |> List.fold (Xor) 0  
    // Class 3: Quadratic bent functions (qbf) 
    let qbf (a: int array) = "to be implemented"
    // Class 4: Dillon's partial spreads (dps)
    let dps (a: int array) = "to be implemented"
    // Class 5: Niho bent functions (nbf)
    let nbf (a: int array) = "to be implemented"

    // ShiftedBentFunction computes the truth table of the Boolean function 
    // on n variables given by x -> f(x+s)
    let ShiftedBentFunction (n:int) (s:int) f = 
        // Boolean shift by a bitvector s       
        let shift (s: int array) (v: int array) = 
            List.zip (Array.toList s) (Array.toList v) 
            |> List.map (fun (x, y) -> Xor x y) 
            |> List.toArray
        // Construct the truth table of the shifted Boolean function
        let TruthTable = List.map (fun x -> BoolInt n x) [0..(pown 2 n)-1] |> List.toArray
        let shiftBool  = BoolInt n s
        Array.map ((shift shiftBool) >> f) TruthTable
             
    // Assemble the unitary gates needed for the circuit
    // GateDiag <list of 0, 1> <qubits> 
    let GateDiag (a:int[]) (qs:Qubits) = 
        let gate (a:int[]) = 
            let len = qs.Length-1 
            Gate(
                Qubits  = qs.Length,
                Name    = sprintf "Diag(%A)" a,
                Help    = sprintf "Diagonal matrix corresponding to %A" a,
                Mat     = (
                    let mat     = CSMat(a.Length ,true)
                    for i in 0..a.Length-1 do   // diagonal phase matrix corresponding to array a
                        mat.r(i,i) <- (pown (-1.0) a.[i])
                    mat),
                Draw    = sprintf "\\multigate{#%d}{%s}" len "U"
                //Draw    = sprintf "\\gate{U}\dwx[#%d]\\go[#%d]\\gate{U}" len len
                //Draw    = "\\gate{L}"
                )
        (gate a).Run qs
    
    // Some permutations used to define Maiorana-McFarland class
    let perm0 = PermBool (n/2) 0    // identity
    let p     = 5                   // specifies the permutation 
    let perm1 = PermBool (n/2) p    // permute cyclically by p
    let perm2 = PermBool (n/2) ((pown 2 (n/2))-p) // inverse of perm1

    // U1 implements the shifted bent function 
    let U1 (n:int) (s:int) = 
        #if IP
            ShiftedBentFunction n s ip |> GateDiag
        #else
            ShiftedBentFunction n s (mmf perm0 perm1)  |> GateDiag
        #endif

    // U0 implements the dual of the unshifted bent function 
    let U0 (n:int) = 
        #if IP
            ShiftedBentFunction n 0 ip |> GateDiag
        #else
            ShiftedBentFunction n 0 (mmf perm2 perm0)  |> GateDiag
        #endif

    // Prepare the initial state; QFT of length p method
    let ket     = Ket(n) // n qubits for state 
    let qs      = ket.Qubits
    let cv      = ket.Single()
    let rslt    = Array.create n 0 // will be used to store final result
    
    // Hidden shift quantum algorithm for bent functions
    let ops (qs:Qubits) =
        H ><   !!(qs,[0..qs.Length-1])   // apply Hadamard to all qubits
        show "Dumping after step 1: "; ket.Dump()
        U1 n s !!(qs,[0..qs.Length-1])   // apply the shifted bent function 
        show "Dumping after step 2: "; ket.Dump()
        H ><   !!(qs,[0..qs.Length-1])   // apply Hadamard to all qubits
        show "Dumping after step 3: "; ket.Dump()
        U0 n   !!(qs,[0..qs.Length-1])   // apply the unshifted bent function 
        show "Dumping after step 4: "; ket.Dump()
        H ><   !!(qs,[0..qs.Length-1])   // apply Hadamard to all qubits
        show "Dumping after step 5: "; ket.Dump()
    
    let circ  = Circuit.Compile ops qs
    let circ    = circ.Fold()
    circ.RenderHT("HiddenShiftGates")
    circ.Dump()
    circ.Run qs

    // measure the register 
    List.iter (fun i -> M !!(qs,i)) [0..n-1] 
    for i in 0..n-1 do 
        rslt.[i] <- qs.[i].Bit.v
    
    let ps      = procStats(true)
    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt

    let checkResult =  
        let sInt = rslt |> Array.rev |> IntBool n // have to reverse order 
        if sInt <> s then 
            show "FAIL! measured bitvector is different from hidden shift"
            show "RES: HShift n s rslt %i %i %i 0" n s sInt
        else 
            show "SUCCESS! The hidden shift is s=%i"   sInt
            show "RES: HShift n s rslt %i %i %i 1" n s sInt
    
    checkResult
           
    show "Hidden shift done"


///////////////////////////////////////////////////////////////////////////
//
// Quantum algorithm for abelian difference sets
//
///////////////////////////////////////////////////////////////////////////

[<LQD>]
/// <summary> Simulate abelian difference set algorithm for cyclic groups </summary>
let __AbelianDifferenceSet (p:int) =
    // Prepare the initial state; QFT of length p method
    let n       = CVec.Bits p // now 2^(n-1) <= p <= 2^n
    let ket     = Ket(n) // n qubits for measurement and n qubits for state 
    let qs      = ket.Qubits
    let rslt    = Array.create n 0  // will be used to store final result
   
    ket.Dump(showInd,0,false,true)
    
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
    
    // QFT inverse for any length <a>. Naive implementation from the definition of the DFT matrix 
    let QFTArbLenInverse (k:int) (qs:Qubits) = 
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
                            mat.r(i,j) <- Math.Cos (PI2k * (float i) * (float j))
                            mat.i(i,j) <- Math.Sin (PI2k * (float i) * (float j))
                    for i in k..len-1 do // fill up to size np2 with identity matrix
                        mat.r(i,i) <- 1.0
                    mat),
                Draw    = "\\gate{DFT_k}"
                )
        (gate k).Run qs
       
    // Quantum algorithm for shifted difference sets
        
    QFTArbLen p !!(qs,[0..qs.Length-1]) // prepares input state on qubits
    show "Dumping after step 1: ";     ket.Dump(showInd,0,false,true)
    
    Udiag16shift1 n !!(qs,[0..qs.Length-1])   // apply the characteristic function of difference set
    //Udiag32shift3 n !!(qs,[0..qs.Length-1])   // apply the characteristic function of difference set
    //Udiag16shift2 n !!(qs,[0..qs.Length-1])   // apply the characteristic function of difference set
    //Udiag8shift2 n !!(qs,[0..qs.Length-1])   // apply the characteristic function of difference set
    show "Dumping after step 2: ";     ket.Dump(showInd,0,false,true)
    
    QFTArbLenInverse p !!(qs,[0..qs.Length-1]) // prepares input state on qubits
    show "Dumping after step 3: ";     ket.Dump(showInd,0,false,true)
    ket.Dump(showInd,0,true,true) 
    // note: change this to ket.Dump(showInd,0,true,true) for MCC plot, i.e., probabilities
    
    Gdiag16 n !!(qs,[0..qs.Length-1])   // apply the characteristic function of difference set
    show "Dumping after step 4: ";     ket.Dump(showInd,0,false,true)
    
    QFTArbLen p !!(qs,[0..qs.Length-1]) // prepares input state on qubits
    show "Dumping after step 5: ";     ket.Dump(showInd,0,false,true)
    ket.Dump(showInd,0,true,true) 
    
    // measure the register 
    List.iter (fun i -> M !!(qs,i)) [0..n-1] 
    for i in 0..n-1 do 
        rslt.[i] <- qs.[i].Bit.v
    
    let ps      = procStats(true)
    show "Result: [mb:%5d m=%A]" ps.wsetMB rslt
       
    show "Difference set done."
    
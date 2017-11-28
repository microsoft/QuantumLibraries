// note that this sample is still work in progress: Build Action is set to None until 
// this harness for the Shor.qs Q# piece is completed. 
// FIXME: complete this sample

namespace Microsoft.Quantum.Samples.ShorSample

open System
open System.Windows
open System.Windows.Controls
open System.Text
open System.Collections.Generic

open Microsoft.Quantum.Simulation.Simulators
open Microsoft.Quantum.Canon
open Microsoft.Quantum.Simulation.Core

open Microsoft.FSharp.Core
open Microsoft.FSharp.Collections

module ShorSample =
    // We begin by defining several utility functions that will be useful to perform the 
    // classical side-calculations that are needed for Shor's algorithm. These include: 
    // precomputation of the constant a used for modular exponentiation, decomposition of 
    // said constant into its parts so that iterative phase estimation can be initialized. 
    // Furthermore, we precompute the inverse of the constant `a` modulo the constant `N`
    // that we want to factor. Finally, utlities are provided that reconstruct the period 
    // of `a` from the estimate of the inverse period returned by the quantum algorithm. 

    // The following function finds the smallest non-trivial factor of a number by brute 
    // force trial division. This is used only for testing for small moduli. The input 
    // `n` is the number to factor. 
    let smallestFactor (n : int) =
        let bound = int(Math.Floor(Math.Sqrt(float n)))
        let rec impl i =
            if i > bound then n
            elif n % i = 0 then i
            else impl (i+1)
        if n < 2 then 1 else impl 2

    // Building on this, we can define a simple primality check. 
    let isPrime (n : int) =
        n > 1 && n = smallestFactor n

    // The following function calculates the greatest common denominator (GCD) of two numbers.
    let rec gcd a b =
        if b = 0 then -1
        elif a%b = 0 then b
        else gcd b (a%b)

    // A coprime is a number that has gcd equal to one with respect to the given modulus. 
    // The following function is a simple heuristic that works for small odd numbers. 
    let coPrime n =
        if n = 221 then 163
        elif n = 493 then 113
        else 2   // Just use 2 for now

    // The following function calculates x^a mod n, using the Russian peasant method.
    // Here x is the base to be exponentiated, a is the exponent the base should be 
    // raised to and n is the modulus. 
    let modExp x a n =
        let rec doLoop a tmp value =
            if a > 0 then
                let value = if a &&& 1 = 1 then (value*tmp) % n else value
                let tmp = (tmp*tmp) % n
                let a = a >>> 1
                doLoop a tmp value
            else value
        doLoop a (x%n) 1

    // This function finds the denominator of the best rational approximation for a given float
    // with a specified maximum denominator. The inputs are the floating-point number f that we 
    // want to approximate, the maximum denominator maxDen to allow. The function returns the 
    // best rational approximation with respect to said denominator. 
    let denominator (f : float) (maxDen : uint64) =
        let eps = 0.5 / Math.Pow(double maxDen,2.0) 
        let rec doNxt (q0 : uint64) (q1 : uint64) (y : float) =
            let z = y - Math.Floor(y)
            if z <= eps then q1
            else 
                let y = 1.0 / z
                let q2 = (uint64 (Math.Floor y)) * q1 + q0
                if q2 > maxDen then q1
                else doNxt q1 q2 y
        doNxt 0UL 1UL f

    // Similarly, the following function computes the best rational approximation, and returns 
    // it as a (numerator, denominator) tuple.
    let rat (f : double) (maxDen : uint64) =
        let mxErr = 0.5 / (float maxDen) 
        let rec iter (first : bool) (d : uint64) (n : uint64) h0 h1 _h2 k0 k1 _k2 =
            let a   = if n = 0UL then 0UL else d / n
            if a = 0UL && not first then h1,k1
            else
                let x = d
                let d = n
                let n = x % n
                let doMore, doFinal, x = if k1 * a + k0 < maxDen then true, true, a
                                         else let x = (maxDen-k0) / k1
                                              if x*2UL >= a || k1 >= maxDen then false, true, x
                                              else false, false, a
                if doFinal 
                then
                    let h2 = x * h1 + h0
                    let h0 = h1
                    let h1 = h2
                    let k2 = x * k1 + k0
                    let k0 = k1
                    let k1 = k2
                    let err = abs (f - ((float h1) / (float k1)))
                    if doMore && err > mxErr then iter false d n h0 h1 h2 k0 k1 k2
                                             else h1, k1
                else h1, k1
        let rec findInt (f : double) (n : uint64) =
            if f = Math.Floor f then (uint64 f), n
                                else findInt (f * 2.0) (n <<< 1)
        let d, n = findInt f 1UL
        iter true d n 0UL 1UL 0UL 1UL 0UL 0UL
    
    // Computing the (classical) reciprocal of a, mod b; that is, (1/a) mod b. 
    let RecipModN a b =
        let rec doSub a b =
            if b = 0 then 1,0
            else
                let q,r = (a/b), (a%b)
                let s,t = doSub b r
                t, (s-q*t)
        let inv, _   = doSub a b
        if inv < 0 then inv+b else inv

    // This function generates all the qubits needed to run Shor's algorithm.
    // Here n is the number to be factored. The function returns a list containing 
    // all of the required qubits.
    // FIXME: this needs to be adapted from the Liquid code
    let genQubits n =
        let bits = CVec.Bits n
        let m = Op.AllocateQubits(1)
        let xs = Op.AllocateQubits(bits)
        let bs = Op.AllocateQubits(bits+1)
        let anc = Op.AllocateQubits(1)   
        let qs = List.concat [m;xs;bs;anc]
        show "Allocated %d qubits" qs.Length
        qs

    // Break the full list of qubits apart into logical registers. Here qs is the complete list of qubits, 
    // as returned from genQubits. The function returns a tuple compaining: n=group size, m=control qubit, 
    // xs=n x-values, bs=n+1 b-values, and an ancilla. Note that both m and the ancilla are returned as 
    // single-element lists rather than as individual qubits.
    // FIXME: this needs to be adapted from the Liquid code
    let getGroups (qs : Qubit list) =
        let n = (qs.Length - 3)/2
        let m = [qs.Head]
        let xs = grab 1 n qs
        let bs = grab (n+1) (n+1) qs
        let anc = grab (2*n+2) 1 qs
        n,m,xs,bs,anc

    // Reset the qubits in a quantum register to match the bit pattern for a classical integer value.
    // Here a is a classical integer to set the register to. The function returns a quantum register qs. 
    // FIXME: this needs to be adapted from the Liquid code
    let setBits a (qs : Qubits) =
        let rec doQ a (qs : Qubit list) =
            match qs with
            | [] -> ()
            | q :: qs' -> let want = if (a &&& 1) = 0 then 1 else -1
                          if want <> M q then X q
                          doQ (a>>>1) qs'
        doQ a qs

    // Measure a quantum register and returns the values as a classical integer.
    // The head of the list goes into the 1s bit, the second element into the 2s bit, etc.
    // Here qs is a quantum register to be measured. The function returns measured values of  
    // the bits in the register, combined into a classical integer.
    // FIXME: this needs to be adapted from the Liquid code
    let getBits (qs : Qubit list) =
        let rec doQ pos rslt (qs : Qubit list) =
            match qs with
            | [] -> rslt
            | q :: qs' ->
                let cur = M q
                let rslt = if cur = -1 then (rslt ||| (1<<<pos)) else rslt
                doQ (pos+1) rslt qs'
        doQ 0 0 qs

    // This function returns a string containing the binary representation of a classical integer.
    // Inputs are width which is the number of bits (characters) in the returned string and 
    // bits, the integer to return a binary representation of. Output is the binary string itself. 
    let toBits width bits =
        let ary = Array.create width '0'
        let rec doBits pos =
            if pos < 0 then System.String(ary)
            else
                if (1<<<pos) &&& bits <> 0 then ary.[width-(pos+1)] <- '1'
                doBits (pos-1)
        doBits (width-1)

    // The following function performs tests of the various operations used in Shor's algorithm. 
    // Input to the function is a number to test factor. This function does not actually perform 
    // the factorization but rather excercises all the operations that are required to factor. 
    // FIXME: this needs to be adapted from the Liquid code
    let ShorGateTests (testInput : int) =
        let qs = genQubits testInput
        let _,m,xs,bs,anc = getGroups qs
        let mutable passed = true

        let testQFT() =
            show "Testing: QFT"
            let circTst (qs : Qubit list) = 
                QFT qs
                QFT' qs
            let rec loop i ok =
                if i < 16 then
                    show "Running test %d" i
                    setBits i bs
                    circTst bs
                    let j = getBits bs
                    if i <> j then
                        show "Test QFT: %s => %s XXXX" (toBits 4 i) (toBits 4 j)
                        loop (i+1) false
                    else loop (i+1) ok
                else ok
            if loop 0 true then show "QFT passed" else show "QFT failed"; passed <- false

        let testAdd() =
            show "Testing: AddA and AddA'"
            let circTst a b = 
                setBits b bs
                QFT bs
                AddA a bs
                QFT' bs
                let i = getBits bs
                i
            let circTst2 a b = 
                setBits b bs
                QFT bs
                AddA a bs
                AddA' a bs
                QFT' bs
                getBits bs
            for a in 0..9 do
                for b in 7..15 do
                    let sum = circTst a b
                    show "Test AddA: %2d+ %2d=%2d or %s+ %s=%s%s"
                        b a sum (toBits 5 b) (toBits 4 a) (toBits 5 sum)
                        (if sum <> a+b then passed <- false; " <=== BAD" else "")
                    let sum = circTst2 a b
                    show "Test AddA':  %2d+-%2d=%2d or %s+-%s=%s%s"
                        b a sum (toBits 5 b) (toBits 4 a) (toBits 5 sum)
                        (if sum <> b then passed <- false; " <=== BAD" else "")

        let testCCAdd() =
            show "Testing :  CCAdd"
            let circTst a b c = 
                let qs = [xs.[0];xs.[1]] @ bs
                setBits b bs
                setBits c xs
                QFT bs
                CCAdd a qs
                QFT' bs
                getBits bs
            let circTst2 a b c = 
                let qs = [xs.[0];xs.[1]] @ bs
                setBits b bs
                setBits c xs
                QFT bs
                CCAdd a qs
                CCAdd' a qs
                QFT' bs
                getBits bs
            for c in 3..3 do
                let a = 6
                let b = 3
                let sum = circTst a b c
                let good = if c = 3 then a+b else b
                show "Test CCAdd: %2d+%2d(%d)=%2d (good=%2d) or %s+%s=%s (good=%s) %s"
                    b a c sum good (toBits 4 b) (toBits 4 a) (toBits 4 sum) (toBits 4 good)
                    (if sum <> good then passed <- false; " <=== BAD" else "")
                let sum = circTst2 a b c
                let good = b
                show "Test CCAdd':  %2d+%2d(%d)=%2d (good=%2d) or %s+%s=%s (good=%s) %s"
                    b a c sum good (toBits 4 b) (toBits 4 a) (toBits 4 sum) (toBits 4 good)
                    (if sum <> good then passed <- false; " <=== BAD" else "")
    
        let testAddModN() =
            show "Testing: AddModN"
            let circTst a b N = 
                setBits 3 xs
                setBits 0 anc
                setBits b bs
                let qs = [xs.[0];xs.[1]] @ bs @ anc
                QFT bs
                AddModN N a qs
                QFT' bs
                getBits bs
            for a in 9..15 do
                for b in 3..9 do
                    let sum = circTst a b testInput
                    let good = (a+b) % testInput
                    if sum <> good then
                        show "Test AddModN : %2d+ %2d=%2d (good=%2d) or %s+ %s=%s (good=%s)%s"
                            b a sum good (toBits 4 b) (toBits 4 a) (toBits 4 sum) (toBits 4 good)
                            (if sum <> good then passed <- false; " <=== BAD" else "")

        let testMulModN() =
            show "Testing: MulModN"
            let circTst c x b a N = 
                setBits c m
                setBits x xs
                setBits b bs
                setBits 0 anc
                let qs = [m.[0]] @ xs @ bs @ anc

                MulModN N a qs
                let gotC = getBits [qs.[0]]
                let gotX = getBits xs
                let gotB = getBits bs
                gotC,gotX,gotB

            for loop in 0..0 do
                let c = 1
                let x = 7
                let b = 7
                let a = 7
                let gotC,gotX,gotB = circTst c x b a testInput
                let good = if c=0 then b else (b+a*x) % testInput
                if gotB <> good then
                    show "Test MulModN: c : %1d b : %2d + a : %2d*x : %2d %% N : %2d=bOut : %2d bGood : =%2d xOut=%2d cOut=%d%s"
                        c b a x testInput gotB good gotX gotC
                        (if gotB <> good then passed <- false; " <=== BAD" else "")

        let testUa() =
            show "Testing: Ua"
            let circTst c x b a N = 
                setBits c m
                setBits x xs
                setBits b bs
                setBits 0 anc
                let qs = [m.[0]] @ xs @ bs @ anc

                Ua N a qs
                let gotC = getBits m
                let gotX = getBits xs
                let gotB = getBits bs
                gotC,gotX,gotB

            for loop in 0..0 do
                let c               = 1
                let x               = 7
                let b               = 0
                let a               = 7
                let gotC,gotX,gotB  = circTst c x b a testInput
                let good            = if c=0 then x else (a*x) % testInput
                if gotX <> good then
                    show "Test Ua :  c : %1d a : %2d*x : %2d %% N : %2d=xOut : %2d xGood : =%2d bOut=%2d cOut=%d%s"
                        c a x testInput gotX good gotB gotC
                        (if gotX <> good then passed <- false; " <=== BAD" else "")

        testQFT()
        testAdd()
        testCCAdd()
        testAddModN()
        testMulModN()
        testUa()
        show "   DONE WITH TESTS"
        show ""
        passed

    // The following function computes the overall accumulated phase during phase estimation.
    let computePhiK (rslt : int[]) j =
        if j=1 then 0.0
        else
            let rec doPhi k rtn =
                if k > j then 2.0 * Math.PI * rtn
                elif rslt.[j-k] = 1 then
                    doPhi (k+1) (rtn+1.0/(pown 2. k)) else doPhi (k+1) rtn
            doPhi 2 0.0
    
    // This is the main function of the sample to demonstrate Shor's algorithm and which runs all 
    // steps of Shor's algorithm to factor an integer N. Here rslt are the bit results from each 
    // step of the algorithm, a is a classical value of the element of which we want to determine 
    // the order, and qs is a quantum register that the algorithm operates on. 
    // FIXME: this needs to be adapted from the Liquid code
    let ShorImpl (rslt : int[]) (N : int) (a : int) (qs : Qubit list) =
        show "Running Shor on qubits %O" qs
        let n,m,xs,_bs,_anc = getGroups qs

        match Op.Model("CircSimModel") with
        | Some o -> let mm = o :?> CircSimModel
                    show "Found model %s" mm.Name
                    mm.Dump()
        | None -> ()

        // Make sure qubits start in |0>
        for q in qs do
            let v = M q
            if v <> 1 then X q
   
        X xs.Head       // Need to start in |1>
 
        let m = m.Head

        let L = 2*n

        for j in 1..L do
            let bitPos = L-j
            let aExp = modExp a (1<<<bitPos) N

            // Make sure m is in |0>
            let v   = M m
            if v <> 1 then X m

            H m
            let theta = computePhiK rslt j
            if j > 1 then EShift theta m
            Ua N aExp qs
            H m
            let v = M m
            let r = if v = 1 then 0 else 1
            rslt.[j-1] <- r
 
            show "        Bit :  %3d [m=%d]" bitPos r

    // Finally, this function is used to check if the measured result is a legal solution,
    // i.e., leading to a non-trivial factor of N. 
    // FIXME: this needs to be adapted from the Liquid code
    let chkResult N a (m : uint64) verbose =
        let n               = CVec.Bits N
        let twoN            = 1UL <<< 2*n

        // Return a failure
        let rtnFail (msg : string) = 
            if verbose then show "%s" msg
            (0,0)

        if m = 0UL then rtnFail"m=0"
        else
            let c       = (double m) / (double twoN)
            let p,den   = rat c twoN
            if verbose then
                show "%8d = m = quantum result" m
                show "%8g = c = %d/%d =~ %d/%d " c m twoN p den

            // Can't use odd period, might be able to use den*2
            let _p,den   =
                if den % 2UL = 1UL && 2UL * den < twoN then 
                    if verbose then show "        Odd denominator, expanding"
                    2UL*p,2UL*den 
                else p,den 

            let den         = int den
            if den % 2 = 1 then
                rtnFail "FAIL Odd period found."
            else
                let den2    = den/2
                let e       = modExp a den2 N
                let v1      = (e+1) % N
                let v2      = (e+N-1) % N
                let f1      = gcd N v1
                let f2      = gcd N v2
                let factor  = max f1 f2
                if verbose then
                    show "%8d = %d/2 = exponent" den2 den
                    show "%8d = %d^%d + 1 mod %d" v1 a den2 N
                    show "%8d = %d^%d - 1 mod %d" v2  a den2 N
                    show "%8d = factor = max(%d,%d)" factor f1 f2
                if factor = -1 then
                    rtnFail "Tried to calc n mod 0 for some n"
                elif factor = N || factor = 1 then
                    rtnFail "Found trivial factors"
                elif factor = 0 then
                    rtnFail "FAIL zero factor found"
                else
                    (factor, (N/factor))

    // As Shor's algorithm is a probabilistic algorithm, it is interesting to compare the probability of 
    // succeeding with the probability of just randomly guessing a solution. This is the purpose of the 
    // present function which computed this probability for a given modulus N and a given generator a. 
    let getRndSucc N a =
        let n           = CVec.Bits N
        let mBits       = 1UL <<< 2*n
        let samples     = if mBits < 1000000UL then mBits else 1000000UL
        let rnd         = Random()
        let getSamp (samp : uint64) =
            if samples = mBits then samp else (uint64 (rnd.Next())) % mBits
        let rec chkOne (m : uint64) good =
            if m >= samples then good,samples,100.0*(float good)/(float samples)
            elif fst (chkResult N a (getSamp m) false) <> 0 then 
                //show "                                     GOOD!"
                chkOne (m+1UL) (good+1)
            else
                //show "                                     <bad>" 
                chkOne (m+1UL) good
        chkOne 0UL 0

    // Harness to run Shor's algorithm: the function just requires a number N to be factored as 
    // input. Selection of a is done using a heuristic. Results are logged. 
    // FIXME: this needs to be adapted from the Liquid code
    let ShorAlgorithm (N : int) =
        if N % 2 = 0 then notLegal "Num must be odd" N
        elif isPrime N then notLegal "Num must not be prime" N
        elif isPrimePower N then notLegal "Num must not be a prime power" N
        else
            show "======== Running Shor's Algorithm ========="
            let n = CVec.Bits N
            let twoN = 1 <<< 2*n
            show "%8d = N = Number to factor" N
            let qsLen = 2*(CVec.Bits N)+3
            let a = coPrime N
            show "%8d = a = coPrime of N" a
            show "%8d = n = number of bits for N" n
            show "%8d = 2^2*n" twoN
            show "%8d = total qubits" qsLen
            let good,total,pcnt = getRndSucc N a
            show "%7.2f%% = prob of random result (%d/%d)" pcnt good total
            show "%7.2f%% = prob of Shor (worst case)" (100.0/(Math.Log(float n,2.)))
            let rslt = Array.create (2*n) 0
            let sw = Diagnostics.Stopwatch()
            sw.Start()
            show "         - Running simulation"
            let qs = genQubits N
            show "%d qubits returned from genQubits :  %O" qs.Length qs
            ShorImpl rslt N a qs
            sw.Stop()
            let elapsed = sw.Elapsed.TotalSeconds
            show "%8g = Elapsed time (seconds)" elapsed
            let m = Array.mapi (fun i bit -> bit <<< i) rslt |> Array.sum
            let r = chkResult N a (uint64 m) true
            if fst r > 0 then show "Found factors %d and %d" (fst r) (snd r) else show "Factoring failed"
            r

   
   // We begin by instantiating the simulator which we will use to run the Q# parts. 
    let qsim = new QuantumSimulator()
    // FIXME: this needs to be called from ShorImpl(..)

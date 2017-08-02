module Microsoft.Research.Liquid.Float

open System
open System.Collections.Generic
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

//open Montgomery
open CircBase           // basic quantum circuits
open Integer            // integer arithmetic


// Shifts the bits of x to the right by the value given in s
// The last M bits of x are/may be occupied, the first 2^(len(s))-1 must be empty.
// Those are used to catch the overflow.
let ShiftRight (M:int) (x:Qubits) (s:Qubits) =
    for m in 0..s.Length-1 do
        let first_bit = max (x.Length - M - ((1 <<< m) - 1)) (1 <<< m)

        for i in first_bit..(x.Length-1) do
            CSWAP [x.[i]] [x.[i - (1 <<< m)]] [s.[m]]

// Does the same but copies the sign along while shifting
let ShiftRightCopySign (M:int) (xS:Qubits) (x:Qubits) (s:Qubits) =
    for m in 0..s.Length-1 do
        let first_bit = max (x.Length - M - ((1 <<< m) - 1)) (1 <<< m)

        for i in first_bit..(x.Length-1) do
            CSWAP [x.[i]] [x.[i - (1 <<< m)]] [s.[m]]

        for i in (x.Length-(1 <<< m))..x.Length-1 do
            CCNOT [s.[m]; xS.[0]; x.[i]]

let ShiftRightCompile (M:int) (qs:Qubits) = 
    ShiftRight M qs.[0..2*M-1] qs.[2*M..qs.Length-1]


// Shifts the bits of x to the left by the value given in s
// The first M bits of x are/may be occupied, the last 2^(len(s))-1 must be empty.
// Those are used to catch the underflow.
let ShiftLeft (M:int) (x:Qubits) (s:Qubits) =
    for m in 0..s.Length-1 do
        let last_bit = M - 1 + (1 <<< m) - 1

        for i in last_bit..(-1)..0 do
            CSWAP [x.[i]] [x.[i+(1 <<< m)]] [s.[m]]

// Shifts the bits of x to the left by the value given in s.
// This version can only be used for re-normalization; the shift register value s has to be s.t.
// no bits of x are set that are higher-order than the value of s.
// Those are used to catch the underflow.
let ShiftLeftShort (x:Qubits) (s:Qubits) =
    let M = x.Length
    for m in 0..s.Length-1 do
        let last_bit = (1<<<m)

        for i in (M-1)..(-1)..last_bit do
            CSWAP [x.[i]] [x.[i-(1 <<< m)]] [s.[m]]

let ShiftLeftCompile (M:int) (qs:Qubits) = 
    ShiftLeft M qs.[0..2*M-1] qs.[2*M..qs.Length-1]

// Ancilla will be 0 if the first 1 is found (and its position is in pos, consisting of log(x.Length) qubits)
// This result can then be used to have a normalized zero: If there's no 1, set exponent to 0.
let PositionOfFirstOne (x:Qubits) (pos:Qubits) (anc:Qubit) =
    X [anc] // set ancilla to 1 ( = first 1 hasn't been found yet)
    CNOT [x.[x.Length-1];anc] // first 1 is in first position, i.e., pos=0 is already correct.

    let bit_is_set_in N k = (N >>> k)&&&1I // true if bit k is set in N

    for i in 1..x.Length-1 do
        for p in 0..pos.Length-1 do
            if ((i>>>p)&&&1) = 1 then
               CCNOT [x.[x.Length-1-i];anc;pos.[p]] // initialize position if current bit is set && position hasn't been found yet

        // flip ancilla if this was indeed the first 1 (i.e., if pos = i)
        let ctrls = [for b in 0..pos.Length-1 do if ((i >>> b)&&&1 = 1) then yield pos.[b]]
        BuildMultiplyControlledNOT ctrls [anc] [x.[0]]

let PositionOfFirstOneCompile (n:int) (p:int) (x:Qubits) =
    PositionOfFirstOne x.[0..n-1] x.[n..n+p-1] x.[n+p]

// Takes the signed mantissa (1 sign bit + mantissa) and converts it to two's complement representation
// using M clean ancilla qubits in state 0, where M = size of the mantissa = SM.Length-1
// SM.[0] is the sign bit and SM.[1:] the mantissa itself.
let SignedMantissa2Complement (SM:Qubits) (anc:Qubits) =
    for i in 0..SM.Length-2 do
        CNOT [SM.[SM.Length-1]; SM.[i]]

    CNOT [SM.[SM.Length-1]; anc.[0]]
    BuildTakahashiModAdder SM.[0..SM.Length-2] anc
    CNOT [SM.[SM.Length-1]; anc.[0]]


let FloatAddition (xS:Qubits) (xM:Qubits) (xE:Qubits) (yS:Qubits) (yM:Qubits) (yE:Qubits) (rS:Qubits) (rM:Qubits) (rE:Qubits) (anc:Qubits) =
    let logM = int (Math.Ceiling (Math.Log((double xM.Length), 2.)))

    
    let overflow_mantissa = anc.[0..xM.Length-1]
    let position_first_one = anc.[xM.Length..xM.Length+xE.Length-1]
    let greater_than = anc.[xM.Length+xE.Length]
    let add_anc = anc.[xM.Length+xE.Length+1..xM.Length+xE.Length+5]
    let zero = anc.[xM.Length+xE.Length+6]
    let add_anc_copy = anc.[xM.Length+xE.Length+7]
    let anc_exp_comparison = anc.[xM.Length+xE.Length+8]

    let P = position_first_one // position of the first one after subtracting and re-shifting to the right if necessary --> also positive

    // PREPARE:
    // sort the two numbers by their exponent s.t. the first number >= second
    CNOT [yE.[yE.Length-1]; greater_than]
    CNOT [xE.[yE.Length-1]; greater_than]
    BuildTakahashiAdderInverse xE (List.concat [yE;[greater_than]])
    BuildTakahashiModAdder xE yE
    CSWAP xE yE [greater_than]
    CSWAP xM yM [greater_than]
    CSWAP xS yS [greater_than]

    // conditioned on the sign-bit, invert both mantissas w.r.t. two's complement (and add the implicit leading 1)
    SignedMantissa2Complement (List.concat [xM; xS]) rM
    SignedMantissa2Complement (List.concat [yM; yS]) rM

    // SHIFT:
    // calculate the positive difference in the exponent \delta E (call the log(M) sized sub-string D)
    BuildTakahashiModAdderInverse xE yE
    let D = xE.[0..logM-1]

    // Compare the difference with M for later
    // initialize a temporary register to M
    let tmp_reg = rM.[0..xE.Length-1]
    for i in 0..logM-1 do
        if ((xM.Length >>> i)&&&1) = 1 then
            X [tmp_reg.[i]]
    // subtract M, catching the overflow into an ancilla qubit
    BuildTakahashiAdderInverse xE (List.concat [tmp_reg;[anc_exp_comparison]])
    BuildTakahashiModAdder xE tmp_reg
    // reset to 0
    for i in 0..logM-1 do
        if ((xM.Length >>> i)&&&1) = 1 then
            X [tmp_reg.[i]]

    // for all m in log(M), perform a swap of distance 2^m conditioned on the m-th bit of D, shifting the second mantissa to the right by 2^m
    // use M clean additional qubits to catch the overflow and also copy the sign bit into the first m bits of the shifted mantissa,
    // conditioned on the m-th bit of D
    ShiftRightCopySign yM.Length yS (List.concat [overflow_mantissa;yM]) D

    // ADD:
    // Carry out an in-place addition, adding the first mantissa to the (first M bits of) the second (now shifted) mantissa.
    // Use 2 clean ancilla qubits: 1 for overlow and 1 for sign.
    let yMplus2anc = (List.concat [yM;add_anc.[0..1]])
    let xMplus2anc = (List.concat [xM;add_anc.[2..3]])
    // initialize pseudo-sign-bits
    CNOT [yS.[0];add_anc.[0]]
    CNOT [yS.[0];add_anc.[1]]
    CNOT [xS.[0];add_anc.[2]]
    CNOT [xS.[0];add_anc.[3]]
    BuildTakahashiModAdder yMplus2anc xMplus2anc

    // RENORMALIZE:
    // Conditioned on the new sign bit, invert the rest of the mantissa w.r.t. 2's complement.
    SignedMantissa2Complement yMplus2anc (List.concat [rM;[add_anc.[4]]])
    // Conditioned on the overflow qubit, shift the entire mantissa by one
    CNOT [add_anc.[0];add_anc_copy]
    ShiftRight (yM.Length+1) (List.concat [[zero];yMplus2anc.[0..yM.Length]]) [add_anc_copy]
    // ... and increment the (old) exponent by 1
    CNOT [add_anc_copy; tmp_reg.[0]]
    BuildTakahashiModAdder xE tmp_reg
    CNOT [add_anc_copy; tmp_reg.[0]]

    
    // Determine the position of the first 1 of the result and write it into a log(M)-sized register P
    PositionOfFirstOne yM P anc.[0] // anc will be 1 if x is all-zero
    // for all m in log(M), perform a swap of distance 2^m conditioned on the m-th bit of P, but this time to the left (instead of right)
    ShiftLeftShort yM P
    
    // copy out the old exponent (it will be the new one but may require modification), if the first 1 was found successfully
    // (i.e. the result is not all-zero)
    BuildTakahashiModAdder xE yE // uncompute the difference
    X [anc.[0]] // anc is now 1 <=> result is non-zero
    for i in 0..rE.Length-1 do
        CCNOT [anc.[0]; xE.[i]; rE.[i]]
    X [anc.[0]] // anc is now 1 <=> result is zero

    
    // copy out the new mantissa if the difference in exponents D < M; else copy out the first mantissa
    // Also, copy the new sign bit and ignore the leading 1 (don't copy).
    for i in 0..rM.Length-1 do
        CCNOT [yM.[i]; anc_exp_comparison; rM.[i]]
    CCNOT [add_anc.[1]; anc_exp_comparison; rS.[0]]
    X [anc_exp_comparison]
    for i in 0..rM.Length-1 do
        CCNOT [xM.[i]; anc_exp_comparison; rM.[i]]
    CCNOT [xS.[0]; anc_exp_comparison; rS.[0]]
    X [anc_exp_comparison]

    
    
    // The new exponent can be found by adding P if D<M.
    BuildCtrlTakahashiModAdderInverse rE P [anc_exp_comparison]
    // uncompute all intermediate results
    // undo sorting
    
    ()

(*
Mutliplies mantissa m1 by mantissa m2 into a 2*n-sized result register r, where n is the number of bits in m1 and m2.
*)
let MultiplyMantissas (m1:Qubits) (m2:Qubits) (r:Qubits) =
    for i in 0..m1.Length-1 do
        BuildCtrlTakahashiAdder  r.[i+1..i+m1.Length] (List.concat [m2;[r.[i+m1.Length+1]]]) [m1.[i]]

let FloatMultiplication (x_S:Qubits) (x_M:Qubits) (x_E:Qubits) (y_S:Qubits) (y_M:Qubits) (y_E:Qubits) (r_S:Qubits) (r_M:Qubits) (r_E:Qubits) (anc:Qubits) =
    let M = y_M.Length
    let extra_M = anc.[0..M-2]
    let L_overflow = anc.[M-1] // qubit to catch overflow when multiplying
    let R_overflow = anc.[M] // qubit to catch underflow when right-shifting for re-normalization
    let L_overflow_cpy = anc.[M+1]
    let tmp_mant = anc.[M+2..2*M+1]
    let tmp_exp = anc.[2*M+2..2*M+1+x_E.Length]
    let tmp_exp_anc = anc.[2*M+2+x_E.Length]

    let extended_res = (List.concat [[R_overflow];extra_M;tmp_mant;[L_overflow]]) // multiply into larger register (2n-1 bits +1 bit for overflow)
    MultiplyMantissas x_M y_M extended_res // multiply the mantissas

    // conditioned on the overflow-bit (copy it first!), shift everything by 1
    CNOT [L_overflow; L_overflow_cpy] // copy out overflow bit
    ShiftRight (M+1) extended_res [L_overflow_cpy] // conditional shift by 1
    CNOT [L_overflow_cpy; tmp_exp.[0]] // "increment" exponent by 1 (if there was overflow and we shifted already)
    BuildCtrlTakahashiModAdder tmp_exp x_E [tmp_mant.[tmp_mant.Length-1]] // if the result is non-zero: add exponent of the first value
    BuildCtrlTakahashiModAdder tmp_exp y_E [tmp_mant.[tmp_mant.Length-1]] // and add the second exponent as well

    // finally: take care of the sign-bit
    CNOT [x_S.[0]; r_S.[0]]
    CNOT [y_S.[0]; r_S.[0]]

    // do NOT copy out mantissa & exponent if either
    // 1) both exponents were negative, resulting exponent is positive (--> result should be 0)
    // or
    // 2) both exponents were positive, resulting exponent is negative (--> result should be inf)

    
    CNOT [x_E.[x_E.Length-1]; tmp_exp.[tmp_exp.Length-1]]
    CNOT [x_E.[x_E.Length-1]; y_E.[y_E.Length-1]]
    X [y_E.[y_E.Length-1]]
    X [tmp_exp_anc]
    CCNOT [y_E.[y_E.Length-1];tmp_exp.[tmp_exp.Length-1];tmp_exp_anc]

    CNOT [x_E.[x_E.Length-1]; tmp_exp.[tmp_exp.Length-1]]
    for i in 0..M-1 do
        CCNOT [tmp_exp_anc;tmp_mant.[i];r_M.[i]]
    for i in 0..tmp_exp.Length-1 do
        CCNOT [tmp_exp_anc;tmp_exp.[i];r_E.[i]]
    CNOT [x_E.[x_E.Length-1]; tmp_exp.[tmp_exp.Length-1]]

    CCNOT [y_E.[y_E.Length-1];tmp_exp.[tmp_exp.Length-1];tmp_exp_anc]
    X [y_E.[y_E.Length-1]]
    X [tmp_exp_anc]
    CNOT [x_E.[x_E.Length-1]; tmp_exp.[tmp_exp.Length-1]]
    CNOT [x_E.[x_E.Length-1]; y_E.[y_E.Length-1]]
    

let FloatMultiplicationCompile (M:int) (E:int) (qs:Qubits) =
    FloatMultiplication [qs.[0]] qs.[1..M] qs.[M+1..M+E] [qs.[M+E+1]] qs.[M+E+2..2*M+E+1] qs.[2*M+E+2..2*M+2*E+1] [qs.[2*(M+E)+2]] qs.[2*M+2*E+3..3*M+2*E+2] qs.[3*M+2*E+3..3*M+3*E+2] qs.[3*(M+E)+3..qs.Length-1]


[<LQD>]
let __FloatMultiplication (x_input:double) (exp_x:int) (y_input:double) (exp_y:int) =
    let M = 53
    let E = 11
    
    let k = Ket(4*(M+E+1)-1-E (* 3 floats (1sign+mantissa+exp) and 1mantissa to perform exact mul, including overflow *)
                + 1 (* copy out overflow from multiplication *)
                + 1 (* catch underflow when right-shifting for renormalization*)
                + M + E + 1 (* output + 1 qubit to catch over/underflow when adding exponents*))
    let qs = k.Qubits
    let circuit = Circuit.Compile (FloatMultiplicationCompile M E) qs
                    |> MyCircuitExport
    let toffcount = (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
    printfn "\nFloating-point multiplication:\n%A-bit mantissa, %A-bit exponent.\n\nToffoli count: %A, T gates: %A, #Qubits: %A\n" (M-1) E toffcount (7*toffcount) qs.Length

    (bigint exp_x) |> BoolInt E |> PrepBool qs 0 1 |> ignore

    let initialState = Array.zeroCreate (qs.Length)
    let x = Math.Abs x_input
    let y = Math.Abs y_input
    if x_input < 0. then
        initialState.[0] <- 1
    if y_input < 0. then
        initialState.[E+M+1] <- 1

    for i in 0..M-1 do
        initialState.[i+1] <- int (((bigint (0.5+x*double (1I<<<(M-1)))) >>> i)&&&1I)
    for i in 0..E-1 do
        initialState.[i+M+1] <- int ((bigint (exp_x) >>> i)&&&1I)
    for i in 0..M-1 do
        initialState.[i+E+M+2] <- int (((bigint (0.5+y*double (1I<<<(M-1)))) >>> i)&&&1I)
    for i in 0..E-1 do
        initialState.[i+E+2*M+2] <- int ((bigint (exp_y) >>> i)&&&1I)
    
    let finalState = MyCircuitSimulateFast circuit initialState
    let xm = [for i in 0..M-1 do yield (bigint finalState.[i+1] <<< i)] |> List.sum
    let xe = [for i in 0..E-1 do yield (bigint finalState.[i+M+1] <<< i)] |> List.sum
    let ym = [for i in 0..M-1 do yield (bigint finalState.[i+M+E+2] <<< i)] |> List.sum
    let ye = [for i in 0..E-1 do yield (bigint finalState.[i+2*M+E+2] <<< i)] |> List.sum
    let zm = [for i in 0..M-1 do yield (bigint finalState.[i+2*(M+E)+3] <<< i)] |> List.sum
    let ze = [for i in 0..E-1 do yield (bigint finalState.[i+M+2*(M+E)+3] <<< i)] |> List.sum
    let anc = [for i in 0..(qs.Length-1-4*(M+E)) do yield (bigint finalState.[i+4*(M+E)] <<< i)] |> List.sum
    
    printfn "register values: a: %A * 2**%A, b': %A * 2**%A, result: %A * 2**%A, anc:%A" (double xm/double (1I <<< (M-1))) xe (double ym/double (1I <<< (M-1))) ye (double zm/double (1I <<< (M-1))) ze anc
    let mutable exp_r = double (1I <<< (int32 ze))
    if ze >= (1I <<< (E-1)) then
        exp_r <- 1./(double (1I <<< ((1 <<< E)-(int32 ze))))
    if finalState.[2*(M+E+1)] = 1 then
        exp_r <- -1. * exp_r
    let actual_result = (x_input * Math.Exp((double exp_x)*Math.Log(2.)))*(y_input * Math.Exp((double exp_y)*Math.Log(2.)))
    printfn "\nPerformed multiplication: %A * %A = %A (should be %A)" (x_input * Math.Exp((double exp_x)*Math.Log(2.))) (y_input * Math.Exp((double exp_y)*Math.Log(2.))) ((double zm/double (1I <<< (M-1))) * exp_r) actual_result
    printfn "abs(result - exact): %A" (Math.Abs (((double zm/double (1I <<< (M-1))) * exp_r) - actual_result))



let FloatAdditionCompile (M:int) (E:int) (qs:Qubits) =
    FloatAddition [qs.[0]] qs.[1..M] qs.[M+1..M+E] [qs.[M+E+1]] qs.[M+E+2..2*M+E+1] qs.[2*M+E+2..2*M+2*E+1] [qs.[2*(M+E)+2]] qs.[2*M+2*E+3..3*M+2*E+2] qs.[3*M+2*E+3..3*M+3*E+2] qs.[3*(M+E)+3..qs.Length-1]


[<LQD>]
let __FloatAddition (x_input:double) (exp_x:int) (y_input:double) (exp_y:int) =
    let M = 24
    let E = 8
    
    let k = Ket(4*(M+E+1)-1 (* 3 floats (1sign+mantissa+exp) and 1mantissa, 1exp extra to catch overflows *)
                + 1 (* swap? *) 
                + 4 (* 2 addition ancilla each *)
                + 1 (* for two's complement inversion *)
                + 1 (* for overflow shift *)
                + 1 (* to copy the overflow indicator *)
                + 1 (* ancilla for exponent comparison *))
    let qs = k.Qubits
    let circuit = Circuit.Compile (FloatAdditionCompile M E) qs
                    |> MyCircuitExport
    let toffcount = (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)

    let logM = (float (Math.Floor (Math.Log((double (M-1)), 2.))))
    let estimate = 2.*(double (M-1))*(logM-2.),+3.5*(double (M-1))*logM + 19.*(double E) + 11.*(double (M-1)) - 5.
    printfn "\nFloating-point addition:\n%A-bit mantissa, %A-bit exponent.\n\nToffoli count: %A, T gates: %A, #Qubits: %A\nEstimate: %A Toffoli gates." (M-1) E toffcount (7*toffcount) qs.Length estimate

    (bigint exp_x) |> BoolInt E |> PrepBool qs 0 1 |> ignore

    let initialState = Array.zeroCreate (qs.Length)
    let x = Math.Abs x_input
    let y = Math.Abs y_input
    if x_input < 0. then
        initialState.[0] <- 1
    if y_input < 0. then
        initialState.[E+M+1] <- 1

    for i in 0..M-1 do
        initialState.[i+1] <- int (((bigint (x*double (1I<<<(M-1)))) >>> i)&&&1I)
    for i in 0..E-1 do
        initialState.[i+M+1] <- int ((bigint (exp_x) >>> i)&&&1I)
    for i in 0..M-1 do
        initialState.[i+E+M+2] <- int (((bigint (y*double (1I<<<(M-1)))) >>> i)&&&1I)
    for i in 0..E-1 do
        initialState.[i+E+2*M+2] <- int ((bigint (exp_y) >>> i)&&&1I)
    
    let finalState = MyCircuitSimulateFast circuit initialState
    let xm = [for i in 0..M-1 do yield (bigint finalState.[i+1] <<< i)] |> List.sum
    let xe = [for i in 0..E-1 do yield (bigint finalState.[i+M+1] <<< i)] |> List.sum
    let ym = [for i in 0..M-1 do yield (bigint finalState.[i+M+E+2] <<< i)] |> List.sum
    let ye = [for i in 0..E-1 do yield (bigint finalState.[i+2*M+E+2] <<< i)] |> List.sum
    let zm = [for i in 0..M-1 do yield (bigint finalState.[i+2*(M+E)+3] <<< i)] |> List.sum
    let ze = [for i in 0..E-1 do yield (bigint finalState.[i+M+2*(M+E)+3] <<< i)] |> List.sum
    let anc = [for i in 0..(qs.Length-1-4*(M+E)) do yield (bigint finalState.[i+4*(M+E)] <<< i)] |> List.sum
    
    printfn "register values: a: %A * 2**%A, b': %A * 2**%A, result: %A * 2**%A, anc:%A" (double xm/double (1I <<< (M-1))) xe (double ym/double (1I <<< (M-1))) ye (double zm/double (1I <<< (M-1))) ze anc
    let mutable exp_r = double (1 <<< (int32 ze))
    if ze >= (1I <<< (E-1)) then
        exp_r <- 1./(double (1I <<< (int32 (1I <<< E)-(int32 ze))))
    if finalState.[2*(M+E+1)] = 1 then
        exp_r <- -1. * exp_r
    let actual_result = (x_input * Math.Exp((double exp_x)*Math.Log(2.)))+(y_input * Math.Exp((double exp_y)*Math.Log(2.)))
    printfn "\nPerformed addition: %A + %A = %A (should be %A)" (x_input * Math.Exp((double exp_x)*Math.Log(2.))) (y_input * Math.Exp((double exp_y)*Math.Log(2.))) ((double zm/double (1I <<< (M-1))) * exp_r) actual_result
    printfn "abs(result - exact): %A" (Math.Abs (((double zm/double (1I <<< (M-1))) * exp_r) - actual_result))


[<LQD>]
let __ShiftCircuit (x:int) (s:int) =
    let M = 8
    let E = 3
    
    let k = Ket(2*M+E)
    let qs = k.Qubits
    let circuit = Circuit.Compile ((*ShiftRightCompile*)ShiftLeftCompile M) qs
                    |> MyCircuitExport

    let initialState = Array.zeroCreate (qs.Length)
    for i in 0..M-1 do
        initialState.[i+M] <- int ((x >>> i)&&&1)
    for i in 0..E-1 do
        initialState.[i+2*M] <- int ((s >>> i)&&&1)

    let finalState = MyCircuitSimulateFast circuit initialState
    let res_m = [for i in 0..2*M-1 do yield (bigint finalState.[i] <<< i)] |> List.sum
    
    printfn "%A\t" ((double res_m)/double (1 <<< M))

[<LQD>]
let __PositionCircuit (x:int) =
    let M = 8
    let E = 3
    
    let k = Ket(M+E+1)
    let qs = k.Qubits
    let circuit = Circuit.Compile (PositionOfFirstOneCompile M E) qs
                    |> MyCircuitExport

    let initialState = Array.zeroCreate (qs.Length)
    for i in 0..M-1 do
        initialState.[i] <- int ((x >>> i)&&&1)

    let finalState = MyCircuitSimulateFast circuit initialState
    let res_m = [for i in 0..E-1 do yield (bigint finalState.[i+M] <<< i)] |> List.sum
    
    printfn "%A\t" res_m
    
module Microsoft.Research.Liquid.Polynomials

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

open CircBase           // basic quantum circuits
open Integer            // integer arithmetic


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Polynomials start here
// Includes:
//  * FMA (fused-multiply-add) using fixed-point 
//  * Multiply (by constant / quantum number) using fixed-point 
//  * EvaluatePolynomial using fixed-point 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Computes a*x into rs using n bits (for input and output)
// NOTE: x must be positive!!! Use the controlled version to handle these cases separately and absorb the sign into the constant to multiply by.
// a: constant to multiply by
// pointpos: point position to use for a, i.e. #digits to the left of the point (for example 01.1001 has pointpos=2)
// xs: x register (quantum register to multiply by a)
// rs: result register (must be 0 initially and will hold a*x in fixedpoint notation, same point position as x)
let FixedPointMultiplyConst (a:double) (pointpos:int) (xs:Qubits) (rs:Qubits) = 
    let n = xs.Length
    let a = bigint (a*double (1I<<<(n-pointpos)))

    // Do the additions that require right-shifts
    for i in 0..(n-pointpos-1) do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        if numitems > 0 then
            let carry = rs.[numitems]
            if (a>>>i)&&&1I = 1I then
                BuildTakahashiAdder rs.[0..numitems-1] (List.concat [addend; [carry]])

    // Do the addition that require left-shifts
    for i in 0..(pointpos-1) do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]

        if i = pointpos-1 then
            if (a>>>(i+(n-pointpos))&&&1I)=1I then
                BuildTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend
        elif (a>>>(i+(n-pointpos)))&&&1I=1I then
            BuildTakahashiModAdder rs.[i..(rs.Length-1)] addend
// Uncomputes the above
let FixedPointMultiplyConstInverse (a:double) (pointpos:int) (xs:Qubits) (rs:Qubits) = 
    let n = xs.Length
    let a = bigint (a*double (1I<<<(n-pointpos)))

    // Do the addition that require left-shifts
    for i in (pointpos-1)..(-1)..0 do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]

        if i = pointpos-1 then
            if (a>>>(i+(n-pointpos))&&&1I)=1I then
                BuildTakahashiModAdder rs.[i..(rs.Length-1)] addend
        elif (a>>>(i+(n-pointpos)))&&&1I=1I then
            BuildTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend

    // Do the additions that require right-shifts
    for i in (n-pointpos-1)..(-1)..0 do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        if numitems > 0 then
            let carry = rs.[numitems]
            if (a>>>i)&&&1I = 1I then
                BuildTakahashiAdderInverse rs.[0..numitems-1] (List.concat [addend; [carry]])

// Controlled version of the multiplier by a constant (see above)
let FixedPointMultiplyConstCtrl (a:double) (pointpos:int) (xs:Qubits) (rs:Qubits) (ctrl:Qubits) = 
    let n = xs.Length
    let a = bigint (a*double (1I<<<(n-pointpos)))

    // Do the additions that require right-shifts
    for i in 0..(n-pointpos-1) do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        let carry = rs.[numitems]
        if numitems > 0 && (a>>>i)&&&1I = 1I then
            BuildCtrlTakahashiAdder rs.[0..numitems-1] (List.concat [addend; [carry]]) ctrl

    // Do the addition that require left-shifts
    for i in 0..(pointpos-1) do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]

        if i = pointpos-1 then
            if (a>>>(i+(n-pointpos))&&&1I)=1I then
                BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend ctrl
        elif (a>>>(i+(n-pointpos)))&&&1I=1I then
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend ctrl
// Uncomputes the above
let FixedPointMultiplyConstCtrlInverse (a:double) (pointpos:int) (xs:Qubits) (rs:Qubits) (ctrl:Qubits) = 
    let n = xs.Length
    let a = bigint (a*double (1I<<<(n-pointpos)))

    // Do the addition that require left-shifts
    for i in (pointpos-1)..(-1)..0 do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]

        if i = pointpos-1 then
            if (a>>>(i+(n-pointpos))&&&1I)=1I then
                BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend ctrl
        elif (a>>>(i+(n-pointpos)))&&&1I=1I then
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend ctrl
    // Do the additions that require right-shifts
    for i in (n-pointpos-1)..(-1)..0 do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        
        let carry = rs.[numitems]
        if numitems > 0 && (a>>>i)&&&1I = 1I then
            BuildCtrlTakahashiAdderInverse rs.[0..numitems-1] (List.concat [addend; [carry]]) ctrl

// Squares xs into rs using an ancilla qubit a
let FixedPointSquare (pointpos:int) (xs:Qubits) (a:Qubits) (rs:Qubits) (ctrls:Qubits) = 
    let n = xs.Length

    // Do the additions that require right-shifts
    for i in 0..(n-pointpos-1) do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        BuildMultiplyControlledNOT (List.concat [[xs.[i]]; ctrls]) [a.[0]] [rs.[0]]
        let carry = rs.[numitems]
        if numitems > 0 then
            BuildCtrlTakahashiAdderSmall rs.[0..numitems-1] (List.concat [addend; [carry]]) a [xs.[(i+1)%xs.Length]]
        BuildMultiplyControlledNOT (List.concat [[xs.[i]]; ctrls]) [a.[0]] [rs.[0]]

    // Do the addition that require left-shifts
    for i in 0..(pointpos-1) do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]
        BuildMultiplyControlledNOT (List.concat [[xs.[i+(n-pointpos)]]; ctrls]) [a.[0]] [rs.[0]]
        if i = pointpos-1 then
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend a
        else
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend a
        BuildMultiplyControlledNOT (List.concat [[xs.[i+(n-pointpos)]]; ctrls]) [a.[0]] [rs.[0]]

// Squares xs into rs using an ancilla qubit a
let FixedPointSquareInverse (pointpos:int) (xs:Qubits) (a:Qubits) (rs:Qubits) (ctrls:Qubits) = 
    let n = xs.Length

    // Do the addition that require left-shifts
    for i in (pointpos-1)..(-1)..0 do
        let numitems = n-i
        let yend = numitems-1
        let addend = xs.[0..yend]
        BuildMultiplyControlledNOT (List.concat [[xs.[i+(n-pointpos)]]; ctrls]) [a.[0]] [rs.[0]]
        if i = pointpos-1 then
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend a
        else
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend a
        BuildMultiplyControlledNOT (List.concat [[xs.[i+(n-pointpos)]]; ctrls]) [a.[0]] [rs.[0]]

    // Do the additions that require right-shifts
    for i in (n-pointpos-1)..(-1)..0 do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = xs.[ystart..yend]
        BuildMultiplyControlledNOT (List.concat [[xs.[i]]; ctrls]) [a.[0]] [rs.[0]]
        
        let carry = rs.[numitems]
        if numitems > 0 then
            BuildCtrlTakahashiAdderInverseSmall rs.[0..numitems-1] (List.concat [addend; [carry]]) a [xs.[(i+1)%xs.Length]]
        BuildMultiplyControlledNOT (List.concat [[xs.[i]]; ctrls]) [a.[0]] [rs.[0]]

    

// Computes x*y into rs using n bits (for input and output)
// NOTE: y must be positive!!! Take care of this by copying the MSB / pseudo-sign bit, inverting y conditioned on it, multiplying, 
// inverting the result conditionally on the sign bit, and then undoing all sign-changes, including setting the sign bit. 
// This costs the same as inverting both operands (modulo 2-s complement), which one would have to do anyway.

// pointpos: point position of x = #digits to the left of the point (for example 01.1001 has pointpos=2)
// xs, ys: 2 n-qubit registers to multiply
// rs: result register (must be 0 initially and will hold x*y in fixedpoint notation, same point position as y)
let FixedPointMultiply (pointpos:int) (xs:Qubits) (ys:Qubits) (rs:Qubits) = 
    let n = xs.Length
    if not (xs.Length = ys.Length) then
        failwith "Register sizes need to be equal to multiply!"

    // Do the additions that require right-shifts
    for i in 0..(n-pointpos-1) do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = ys.[ystart..yend]
        
        let carry = rs.[numitems]
        if numitems > 0 then
            BuildCtrlTakahashiAdderSmall rs.[0..numitems-1] (List.concat [addend; [carry]]) [xs.[i]] [xs.[(i+1)%xs.Length]]

    // Do the addition that require left-shifts
    for i in 0..(pointpos-1) do
        let numitems = n-i
        let yend = numitems-1
        let addend = ys.[0..yend]

        if i = pointpos-1 then
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]
        else
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]

let FixedPointMultiplyInverse (pointpos:int) (xs:Qubits) (ys:Qubits) (rs:Qubits) = 
    let n = xs.Length
    if not (xs.Length = ys.Length) then
        failwith "Register sizes need to be equal to multiply!"

    // Do the addition that require left-shifts
    for i in (pointpos-1)..(-1)..0 do
        let numitems = n-i
        let yend = numitems-1
        let addend = ys.[0..yend]

        if i = pointpos-1 then
            BuildCtrlTakahashiModAdder rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]
        else
            BuildCtrlTakahashiModAdderInverse rs.[i..(rs.Length-1)] addend [xs.[i+(n-pointpos)]]

    // Do the additions that require right-shifts
    for i in (n-pointpos-1)..(-1)..0 do
        let numitems = min n (pointpos+i)
        let ystart = n-numitems
        let yend = ystart+numitems-1
        let addend = ys.[ystart..yend]
        
        let carry = rs.[numitems]
        if numitems > 0 then
            BuildCtrlTakahashiAdderInverseSmall rs.[0..numitems-1] (List.concat [addend; [carry]]) [xs.[i]] [xs.[(i+1)%xs.Length]]


// Adds a fixed-point constant to a register
// Depending on the size of the ancilla register it either uses the inplace adder (nlogn) (for 1 or 2 ancillae)
// or Takahashi (O(n)) when n clean ancillae are given.
let AddFixedPointConstant (ppos:int) (a:double) (xs:Qubits) (anc:Qubits) (ctrls:Qubits) =
    let n = xs.Length
    if anc.Length = n then
        let intA = bigint (a*double (1I<<<(n-ppos)))
        for i in 0..(n-1) do
            if (intA>>>i)&&&1I = 1I then
                BuildMultiplyControlledNOT ctrls [anc.[i]] [xs.[0]]
        BuildTakahashiModAdder xs anc
        for i in 0..(n-1) do
            if (intA>>>i)&&&1I = 1I then
                BuildMultiplyControlledNOT ctrls [anc.[i]] [xs.[0]]
    elif anc.Length = 2 || anc.Length = 1 then
        ConstantAdderInplace (bigint (a*double (1I<<<(n-ppos)))) xs anc ctrls
    else
        failwith "AddFixedPointConstant: Provide either n clean (0) ancillae or 1 or 2 dirty ancillae to do inplace addition."

// Computes x*y+a into rs using n bits (for input and output)
// NOTE: Use the one below if ancillae are available!!! (it's faster)
// pointpos[x/y]: point position of x/y = #digits to the left of the point (for example 01.1001 has pointpos=2)
// xs, ys: 2 n-qubit registers to multiply
// rs: result register (must be 0 initially and will hold x*y+a in fixedpoint notation, same point position as y)
// subtractctrl: Either empty list [] or list containing one qubit which indicates whether to subtract or add the constant
let FusedMultiplyAddConstant (pointposx:int) (pointposy:int) (a:double) (xs:Qubits) (ys:Qubits) (rs:Qubits) (subtractctrl:Qubits) = 
    let n = xs.Length
    FixedPointMultiply pointposx xs ys rs
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    if (abs a) > 1.e-14 then
        AddFixedPointConstant pointposy a rs [ys.[0];xs.[0]] []
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
// Uncomputes the above
let FusedMultiplyAddConstantInverse (pointposx:int) (pointposy:int) (a:double) (xs:Qubits) (ys:Qubits) (rs:Qubits) (subtractctrl:Qubits) = 
    let n = xs.Length
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    if (abs a) > 1.e-14 then
        AddFixedPointConstant pointposy -a rs [ys.[0];xs.[0]] []
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    
    FixedPointMultiplyInverse pointposx xs ys rs

// Computes x*y+a into rs using n bits (for input and output) and n ancilla qubits to be faster than the implementation above
// pointpos[x/y]: point position of x / y = #digits to the left of the point (for example 01.1001 has pointpos=2)
// xs, ys: 2 n-qubit registers to multiply
// rs: result register (must be 0 initially and will hold a*x in fixedpoint notation, same point position as y)
// anc: n-qubit ancilla-register, all-zero
// subtractctrl: Either empty list [] or list containing one qubit which indicates whether to subtract or add the constant
let FastFusedMultiplyAddConstant (pointposx:int) (pointposy:int) (a:double) (xs:Qubits) (ys:Qubits) (rs:Qubits) (anc:Qubits) (subtractctrl:Qubits) = 
    let n = xs.Length
    FixedPointMultiply pointposx xs ys rs
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    AddFixedPointConstant pointposy a rs anc []
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
// Uncomputes the above
let FastFusedMultiplyAddConstantInverse (pointposx:int) (pointposy:int) (a:double) (xs:Qubits) (ys:Qubits) (rs:Qubits) (anc:Qubits) (subtractctrl:Qubits) = 
    let n = xs.Length
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    AddFixedPointConstant pointposy -a rs anc []
    if subtractctrl.Length > 0 then
        for i in 0..(n-1) do
            CNOT [subtractctrl.[0];rs.[i]]
    FixedPointMultiplyInverse pointposx xs ys rs


// Evaluates the polynomial 
// c_0 x^d + c_1 x^{d-1} + ... + c_d
// with coefficients c=coeffs
// 
// input:
// ppos - point position to use for fixedpoint representation
// coeffs - coefficient list of doubles
// space - n*(d+2)+1 qubits where the first n hold x and the last n+1 hold the result and a 0-ancilla (used as an explicit sign-bit)
let EvaluatePolynomial (ppos:int) (n:int) (coeffs:list<double>) (space:Qubits) =
    let deg = (coeffs.Length-1)
    
    if deg < 1 then
        failwith "Degree < 1? Why would you call this function?"

    let x = space.[0..(n-1)]
    let r = space.[n..(2*n-1)]

    // invert x if needed:
    // copy sign bit
    let signbit = space.[space.Length-1]
    CNOT [x.[n-1];signbit]
    // invert x by first bitflipping and then adding 1
    for i in 0..(n-1) do
        CNOT [signbit; x.[i]]
    CNOT [signbit;r.[0]]
    BuildTakahashiModAdder x r
    CNOT [signbit;r.[0]]

    // do polynomial evaluation on x>0
    if deg % 2 = 0 then
        FixedPointMultiplyConst coeffs.[0] ppos x r 
    else
        FixedPointMultiplyConstCtrl -coeffs.[0] ppos x r [signbit]
        X [signbit]
        FixedPointMultiplyConstCtrl coeffs.[0] ppos x r [signbit]
        X [signbit]

    if deg = 1 then // no work qubits left --> use inplace O(nlogn) addition
        AddFixedPointConstant ppos coeffs.[1] r [x.[0];x.[1]] []
    else
        let rnext = space.[2*n..3*n-1]
        if (deg-1) % 2 = 1 then
            for i in 0..(n-1) do
                CNOT [signbit; rnext.[i]]
            AddFixedPointConstant ppos coeffs.[1] r rnext []
            for i in 0..(n-1) do
                CNOT [signbit; rnext.[i]]
        else
            AddFixedPointConstant ppos coeffs.[1] r rnext []

    for i in 1..(deg-1) do
        let y = space.[i*n..(i+1)*n-1]
        let r = space.[(i+1)*n..(i+2)*n-1]
        if i < deg-1 then // we have ancilla qubits available --> do fast version
            let rnext = space.[(i+2)*n..(i+3)*n-1]
            if (deg-i-1) % 2 = 0 then
                FastFusedMultiplyAddConstant ppos ppos coeffs.[i+1] y x r rnext []
            else
                FastFusedMultiplyAddConstant ppos ppos coeffs.[i+1] y x r rnext [signbit]
        else
            if (deg-i-1) % 2 = 0 then
                FusedMultiplyAddConstant ppos ppos coeffs.[i+1] y x r []
            else
                FusedMultiplyAddConstant ppos ppos coeffs.[i+1] y x r [signbit]

    // uncompute intermediate results
    for i in (deg-2)..(-1)..1 do
        let y = space.[i*n..(i+1)*n-1]
        let r = space.[(i+1)*n..(i+2)*n-1]
        if i < deg-2 then // we have ancilla qubits available --> do fast version
            let rnext = space.[(i+2)*n..(i+3)*n-1]
            if (deg-i-1) % 2 = 0 then
                FastFusedMultiplyAddConstantInverse ppos ppos coeffs.[i+1] y x r rnext []
            else
                FastFusedMultiplyAddConstantInverse ppos ppos coeffs.[i+1] y x r rnext [signbit]
        else
            if (deg-i-1) % 2 = 0 then
                FusedMultiplyAddConstantInverse ppos ppos coeffs.[i+1] y x r []
            else
                FusedMultiplyAddConstantInverse ppos ppos coeffs.[i+1] y x r [signbit]
    

    if deg > 1 then
        if deg = 2 then // no work qubits left --> use inplace O(nlogn) addition
            for i in 0..(n-1) do
                CNOT [signbit; r.[i]]
            AddFixedPointConstant ppos -coeffs.[1] r [x.[0];x.[1]] []
            for i in 0..(n-1) do
                CNOT [signbit; r.[i]]
        else
            let rnext = space.[2*n..3*n-1]
            if (deg-1) % 2 = 1 then
                for i in 0..(n-1) do
                    CNOT [signbit; rnext.[i]]
                AddFixedPointConstant ppos -coeffs.[1] r rnext []
                for i in 0..(n-1) do
                    CNOT [signbit; rnext.[i]]
            else
                AddFixedPointConstant ppos -coeffs.[1] r rnext []
        if deg % 2 = 0 then
            FixedPointMultiplyConstInverse coeffs.[0] ppos x r 
        else
            X [signbit]
            FixedPointMultiplyConstCtrlInverse coeffs.[0] ppos x r [signbit]
            X [signbit]
            FixedPointMultiplyConstCtrlInverse -coeffs.[0] ppos x r [signbit]

    let r = space.[n..2*n-1]
    CNOT [signbit;r.[0]]
    BuildTakahashiModAdderInverse x r
    CNOT [signbit;r.[0]]
    
    for i in 0..(n-1) do
        CNOT [signbit; x.[i]]
    CNOT [x.[n-1];signbit]
    ()


[<LQD>]
let __RunPolynomial (n:int) (pointpos:int) (x:double) (*(a0:double) (a1:double) (a2:double)*) =
    //let a = [a0;a1;a2]
    let a = [0.0401228098468; 0.012316781492; 0.0324722511722; 0.0444170327186; 0.0750123462156; 0.166666371194; 1.00000000204; ]
    let k = Ket(a.Length*n+1)
    let qs = k.Qubits
    let pointposition = pointpos

    let circuit = Circuit.Compile (EvaluatePolynomial pointposition n a) qs
                    |> MyCircuitExport
    // prepare x
    (bigint (x*double (1I<<<(n-pointposition)))) |> fun x -> (x) |> BoolInt n |> PrepBool qs 0 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (a.Length*n+1)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(n-1) do 
        initialState.[i] <- qs.[i].Bit.v
    for i in n..(initialState.Length-1) do
        initialState.[i] <- 0

    let f = double (1I<<<(n-pointposition))
    let mutable expres = 0.
    for i in 0..(a.Length-1) do
        expres <- x*expres+a.[i]

    show "Calculating P(x) = %A" expres
    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in ((a.Length-1)*n)..(initialState.Length-2) do yield (bigint finalState.[i] <<< (i-(a.Length-1)*n))]
                |> List.sum
    let rest = [for i in n..((a.Length-1)*n-1) do yield (bigint finalState.[i] <<< (i-n))]
                |> List.sum

    let res = ((double res)/f)
    show "Result: %A" res
    //show "sin(x) = %A" (sin x)
    show "\nDifference : %A" (abs (expres-res))
    show "\nUncomp = %A" rest
    //show "Discrep: %A, Expected: %A" (abs (res-expres)) ((double n)*f)
    ()


// TODO: Optimize as follows:
//       Don't add the constant of each polynomial, but only set the extra register to those values
//       then, run the adder ONCE, prior to uncomputing the constants again.
//
// Evaluates M polynomials in parallel
// c_00 x^d + c_01 x^{d-1} + ... + c_0d
// c_10 x^d + c_11 x^{d-1} + ... + c_1d
// ...
// c_m0 x^d + c_m1 x^{d-1} + ... + c_md
// with coefficients c=coeffs 
// 
// input:
// ppos - point position to use for fixedpoint representation
// coeffs - list of coefficient lists: coeffs[1][0] gives the coefficient of highest power for the second polynomial
// space - log(M)+1+n*(d+1)+1 qubits where the first ceil(log(M)) hold the label (i.e. which polynomial to evaluate),
//            then 1 ancilla qubit to avoid multiply-controlled NOTS.
//          The next n hold x and the last n+1 hold the result and a 0-ancilla (used as an explicit sign-bit)
// space = [log(M),1,n (x), n*(d-1), n (result), 1]
let EvaluatePolynomialsInParallel (ppos:int) (n:int) (coeffs:list<list<double>>) (space:Qubits) =
    let M = coeffs.Length // # polynomials to evaluate
    let deg = (coeffs.[0].Length-1)
    
    if deg < 1 then
        failwith "Degree < 1? Why would you call this function?"

    let offset = int (Math.Ceiling (Math.Log((double M), 2.)))

    let label = space.[0..offset-1]
    let ctrl = space.[offset]
    let space = space.[offset+1..space.Length-1]
    let x = space.[0..(n-1)]
    let r = space.[n..(2*n-1)]

    let ResetLabel () = for i in 0..label.Length-1 do
                        if ((M-1)>>>i)&&&1 = 0 then
                            X [label.[i]]
    // invert x if needed:
    // copy sign bit
    let signbit = space.[space.Length-1]
    CNOT [x.[n-1];signbit]
    // invert x by first bitflipping and then adding 1
    for i in 0..(n-1) do
        CNOT [signbit; x.[i]]
    CNOT [signbit;r.[ppos-1]]
    BuildTakahashiModAdder x r
    CNOT [signbit;r.[ppos-1]]

    // do polynomial evaluation on x>0
    // do first multiplication separately
    for m in 0..M-1 do
        for i in 0..label.Length-1 do
            if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                X [label.[i]]
        if deg % 2 = 0 then
            BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
            FixedPointMultiplyConstCtrl coeffs.[m].[0] ppos x r [ctrl]
            BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
        else
            BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
            FixedPointMultiplyConstCtrl -coeffs.[m].[0] ppos x r [ctrl]
            BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
            X [signbit]
            BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
            FixedPointMultiplyConstCtrl coeffs.[m].[0] ppos x r [ctrl]
            BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
            X [signbit]
    // reset label register (ceil(log2(M)) >= M and not =)
    ResetLabel ()
    // Do first addition separately
    for m in 0..M-1 do
        for i in 0..label.Length-1 do
            if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                X [label.[i]]
        BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
        if deg = 1 then // no work qubits left --> use inplace O(nlogn) addition
            AddFixedPointConstant ppos coeffs.[m].[1] r [x.[0];x.[1]] [ctrl]
        else
            let rnext = space.[2*n..3*n-1]
            AddFixedPointConstant ppos coeffs.[m].[1] r rnext [ctrl]
        BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
    // reset label register (ceil(log2(M)) >= M and not =)
    ResetLabel ()
    
    for p in 1..(deg-1) do
        let y = space.[p*n..(p+1)*n-1]
        let r = space.[(p+1)*n..(p+2)*n-1]
        FixedPointMultiply ppos y x r
        if deg-p % 2 = 1 then
            for i in 0..n-1 do
                CNOT [signbit; r.[i]]
        if p >= deg-1 then
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                AddFixedPointConstant ppos c r [x.[0];x.[1]] [ctrl]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            ResetLabel ()
        else
            let rnext = space.[(p+2)*n..(p+3)*n-1]
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                let intC = bigint (c*double (1I<<<(n-ppos)))
                for i in 0..(n-1) do
                    if (intC>>>i)&&&1I = 1I then
                        CNOT [ctrl;rnext.[i]]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            ResetLabel ()
            BuildTakahashiModAdder r rnext
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                let intC = bigint (c*double (1I<<<(n-ppos)))
                for i in 0..(n-1) do
                    if (intC>>>i)&&&1I = 1I then
                        CNOT [ctrl;rnext.[i]]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            ResetLabel ()
        if deg-p % 2 = 1 then
            for i in 0..n-1 do
                CNOT [signbit; r.[i]]
    // uncompute intermediate results
    
    for p in (deg-2)..(-1)..1 do
        let y = space.[p*n..(p+1)*n-1]
        let r = space.[(p+1)*n..(p+2)*n-1]
        if deg-p % 2 = 1 then
            for i in 0..n-1 do
                CNOT [signbit; r.[i]]
        if p >= deg-2 then
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                AddFixedPointConstant ppos -c r [x.[0];x.[1]] [ctrl]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            ResetLabel ()
        else
            let rnext = space.[(p+2)*n..(p+3)*n-1]
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                let intC = bigint (c*double (1I<<<(n-ppos)))
                for i in 0..(n-1) do
                    if (intC>>>i)&&&1I = 1I then
                        CNOT [ctrl;rnext.[i]]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            BuildTakahashiModAdderInverse r rnext
            ResetLabel ()
            for m in 0..M-1 do
                for i in 0..label.Length-1 do
                    if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                        X [label.[i]]
                let c = coeffs.[m].[p+1]
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
                
                let intC = bigint (c*double (1I<<<(n-ppos)))
                for i in 0..(n-1) do
                    if (intC>>>i)&&&1I = 1I then
                        CNOT [ctrl;rnext.[i]]
                
                BuildMultiplyControlledNOT label [ctrl] [space.[0]]
            ResetLabel ()
        if deg-p % 2 = 1 then
            for i in 0..n-1 do
                CNOT [signbit; r.[i]]
        FixedPointMultiplyInverse ppos y x r
        

    let r = space.[n..(2*n-1)]
    if deg > 1 then
        // Undo first addition separately
        for m in 0..M-1 do
            for i in 0..label.Length-1 do
                if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                    X [label.[i]]
            let rnext = space.[2*n..3*n-1]
            BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
            AddFixedPointConstant ppos -coeffs.[m].[1] r rnext [ctrl]
            BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
        // reset label register (ceil(log2(M)) >= M and not =)
        ResetLabel ()
        for m in 0..M-1 do
            for i in 0..label.Length-1 do
                if (((m-1)^^^(m))>>>i)&&&1 = 1 || m = 0 then
                    X [label.[i]]
            if deg % 2 = 0 then
                BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
                FixedPointMultiplyConstCtrlInverse coeffs.[m].[0] ppos x r [ctrl]
                BuildMultiplyControlledNOT label [ctrl] [space.[2*n]]
            else
                BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
                FixedPointMultiplyConstCtrlInverse -coeffs.[m].[0] ppos x r [ctrl]
                BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
                X [signbit]
                BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
                FixedPointMultiplyConstCtrlInverse coeffs.[m].[0] ppos x r [ctrl]
                BuildMultiplyControlledNOT (List.concat [label;[signbit]]) [ctrl] [space.[2*n]]
                X [signbit]
        // reset label register (ceil(log2(M)) >= M and not =)
        ResetLabel ()
         
    CNOT [signbit;r.[ppos-1]]
    BuildTakahashiModAdderInverse x r
    CNOT [signbit;r.[ppos-1]]
    for i in 0..(n-1) do
        CNOT [signbit; x.[i]]
    CNOT [x.[n-1];signbit]
    ()

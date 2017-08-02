module Microsoft.Research.Liquid.Numerics

open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

open CircBase           // basic quantum circuits
open Integer            // integer arithmetic
open Polynomials            // polynomial evaluation


// WARNING: OLD! (But optimized as opposed to the new version (below))
// Performs initial Newton step for inverse square root, including finding 2^(log2(y)/2)
// Requires pointposition (fixed point representation), y = input, x0 = output (initial guess)
// and n+n/2+1 extraspace which should be available anyway, since we'll do newton afterwards
let ApproximateInverseSqrt (pointpos:int) (ys:Qubits) (x0:Qubits) (space:Qubits) =
    let n = ys.Length
    // do approximate log2 and first Newton iteration in one step
    let flag = space.[0]
    let tmplog2space = slice space [1..n/2]
    let extraspace = slice space [n/2+1..n+n/2]
    X [flag] // set flag to true, i.e. we haven't encountered the first 1 yet
    let yhalf = slice ys [1..(n-1)] // y/2 <-> right-shift
    let constant_left = 1.442 
    let constant_right = 1.73253
    // handle even #bits to the left of the binary point
    if pointpos%2 = 0 then
        CNOT [ys.[n-1]; flag]
        X [flag]
        let curlog2half = pointpos/2
        let idx = n-pointpos-curlog2half
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_left*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [flag; x0.[idx-k]]
            k <- k + 1
        let yadder = slice yhalf [(3*curlog2half)..(n-2)]
        let yadder = List.concat [yadder;extraspace.[0..3*curlog2half-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdderInverse x0.[0..(n-2)] yadder [flag]
        X [flag]
    // do digits to the left of the binary point in groups of two, since log2/2 is identical for those digits
    for i in (2-(pointpos%2))..2..(pointpos-2) do 
        let c1 = ys.[n-i]
        let c2 = ys.[n-i-1]
        let curlog2half = (pointpos-i)/2 // current log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[n/2-curlog2half]
        // check if (at least) one of the two is 1 --> -log2(y)/2 = curlog2half
        CCNOT [flag;c2;anc]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        // update flag (did we encounter the first 1 just now?)
        CNOT [anc;flag]
        // initialize to 1.5 shifted by log2(y)/2 (later, init to 1.7325 or 1.442, depending on whether log is positive or negative)
        let idx = n-pointpos-curlog2half
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_left*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [anc; x0.[idx-k]]
            k <- k + 1
        let yadder = slice yhalf [(3*curlog2half)..(n-2)]
        let yadder = List.concat [yadder;extraspace.[0..3*curlog2half-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdderInverse x0.[0..(n-2)] yadder [anc]
    // do digits to the right of the binary point
    for i in (n-1-pointpos)..(-2)..0 do 
        let c1 = ys.[i+1]
        let c2 = ys.[i]
        let curlog2half = (n-1-pointpos-i)/2 // current -log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[curlog2half]
        // check if (at least) one of the two is 1 --> -log2(y)/2 = curlog2half
        CCNOT [flag;c2;anc]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        // update flag (did we encounter the first 1 just now?)
        CNOT [anc;flag]
        // initialize to 1.5 shifted by log2(y)/2 (later, init to 1.7325 or 1.442, depending on whether log is positive or negative)
        let idx = n-pointpos+curlog2half
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_right*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [anc; x0.[idx-k]]
            k <- k + 1
        let yadder = List.concat [(slice yhalf [0..(max (n-2-3*curlog2half) 0)]);[x0.[n-1]]]
        BuildCtrlTakahashiAdderInverseSmall x0.[(min (3*curlog2half) (n-2))..(n-2)] yadder [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
    // do odd left overs (if there are any) and handle all-zero. Both imply that flag = 1 still
    // odd left overs:
    let lim = min pointpos ((n-pointpos+1)/2)
    if (n-1-pointpos)%2 = 1 then
        // all-zero
        X [ys.[0]]
        for i in 0..(lim-1) do
            CCNOT [flag; ys.[0]; x0.[n-pointpos+i]]
        X [ys.[0]]
        // non-zero
        let curlog2half = (n-1-pointpos)/2
        let idx = n-pointpos+curlog2half
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_right*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CCNOT [flag; ys.[0]; x0.[idx-k]]
            k <- k + 1
        let yadder = List.concat [(slice yhalf [0..(max (n-2-3*curlog2half) 0)]);[x0.[n-1]]]
        if yadder.Length > 1 then
            BuildMultiCtrlTakahashiAdderInverse x0.[(min (3*curlog2half) (n-2))..(n-2)] yadder [flag; ys.[0]]
    else // handle all-zero if there are no left-overs
        for i in 0..(lim-1) do
            CNOT [flag; x0.[n-pointpos+i]]

    // clear temporary space
    for i in ((n-1-pointpos)%2)..2..(n-1-pointpos) do 
        let c1 = ys.[i+1]
        let c2 = ys.[i]
        let curlog2half = (n-1-pointpos-i)/2 // current -log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[curlog2half]
        CNOT [anc;flag]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CCNOT [flag;c2;anc]
    for i in (pointpos-2)..(-2)..1 do 
        let c1 = ys.[n-i]
        let c2 = ys.[n-i-1]
        let curlog2half = (pointpos-i)/2 // current log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[n/2-curlog2half]
        CNOT [anc;flag]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CCNOT [flag;c2;anc]
    if pointpos%2 = 0 then
        CNOT [ys.[n-1]; flag]
    X [flag]

// WARNING: OLD! (But optimized as opposed to the new version (below))
// Performs the inverse of the above (used to uncompute)
// Requires pointposition (fixed point representation), y = input, x0 = output (initial guess)
// and n+n/2 extraspace which should be available anyway, since we'll do newton afterwards
let ApproximateInverseSqrtInverse (pointpos:int) (ys:Qubits) (x0:Qubits) (space:Qubits) =
    let n = ys.Length
    // do approximate log2 and first Newton iteration in one step
    let flag = space.[0]
    let tmplog2space = slice space [1..n/2]
    let extraspace = slice space [n/2+1..n+n/2]
    let yhalf = slice ys [1..(n-1)] // y/2 <-> right-shift
    let constant_left = 1.442
    let constant_right = 1.73253

    X [flag] // set flag to true, i.e. we haven't encountered the first 1 yet
    // undo clear temporary space
    if pointpos%2 = 0 then
        CNOT [ys.[n-1]; flag]
    for i in (2-(pointpos%2))..2..(pointpos-2) do 
        let c1 = ys.[n-i]
        let c2 = ys.[n-i-1]
        let curlog2half = (pointpos-i)/2 // current log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[n/2-curlog2half]
        CCNOT [flag;c2;anc]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CNOT [anc;flag]

    for i in (n-1-pointpos)..(-2)..((n-1-pointpos)%2) do 
        let c1 = ys.[i+1]
        let c2 = ys.[i]
        let curlog2half = (n-1-pointpos-i)/2 // current -log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[curlog2half]
        CCNOT [flag;c2;anc]
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CNOT [anc;flag]
        
    // undo odd left overs (if there are any) and handle all-zero. Both imply that flag = 1 still
    // undo odd left overs:
    let lim = min pointpos ((n-pointpos+1)/2)
    if (n-1-pointpos)%2 = 1 then
        // non-zero
        let curlog2half = (n-1-pointpos)/2
        let idx = n-pointpos+curlog2half
        let yadder = List.concat [(slice yhalf [0..(n-2-3*curlog2half)]);[x0.[n-1]]]
        if yadder.Length > 1 then
            BuildMultiCtrlTakahashiAdder x0.[3*curlog2half..(n-2)] yadder [flag; ys.[0]]
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_right*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CCNOT [flag; ys.[0]; x0.[idx-k]]
            k <- k + 1
        // all-zero
        X [ys.[0]]
        for i in 0..(lim-1) do
            CCNOT [flag; ys.[0]; x0.[n-pointpos+i]]
        X [ys.[0]]
    else // handle all-zero if there are no left-overs
        for i in 0..(lim-1) do
            CNOT [flag; x0.[n-pointpos+i]]

    // undo digits to the right of the binary point
    for i in ((n-1-pointpos)%2)..2..(n-1-pointpos) do 
        let c1 = ys.[i+1]
        let c2 = ys.[i]
        let curlog2half = (n-1-pointpos-i)/2 // current -log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[curlog2half]
        let yadder = List.concat [(slice yhalf [0..(n-2-3*curlog2half)]);[x0.[n-1]]]
        BuildCtrlTakahashiAdderSmall x0.[3*curlog2half..(n-2)] yadder [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        // uninitialize to 1.5 shifted by log2(y)/2 (later, init to 1.7325 or 1.442, depending on whether log is positive or negative)
        let idx = n-pointpos+curlog2half
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_right*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [anc; x0.[idx-k]]
            k <- k + 1
        // update flag (did we encounter the first 1 just now?)
        CNOT [anc;flag]
        // check if (at least) one of the two is 1 --> -log2(y)/2 = curlog2half
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CCNOT [flag;c2;anc]

    // undo digits to the left of the binary point in groups of two, since log2/2 is identical for those digits
    for i in (pointpos-2)..(-2)..(2-(pointpos%2)) do 
        let c1 = ys.[n-i]
        let c2 = ys.[n-i-1]
        let curlog2half = (pointpos-i)/2 // current log2(y)/2 if we encounter the first 1 now
        let anc = tmplog2space.[n/2-curlog2half]
        let yadder = slice yhalf [(3*curlog2half)..(n-2)]
        let yadder = List.concat [yadder;extraspace.[0..3*curlog2half-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdder x0.[0..(n-2)] yadder [anc]
        let idx = n-pointpos-curlog2half
        let mutable k = 0
        // initialize to 1.5 shifted by log2(y)/2 (later, init to 1.7325 or 1.442, depending on whether log is positive or negative)
        while idx-k >= 0 do
            if (bigint (constant_left*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [anc; x0.[idx-k]]
            k <- k + 1
        // update flag (did we encounter the first 1 just now?)
        CNOT [anc;flag]
        // check if (at least) one of the two is 1 --> -log2(y)/2 = curlog2half
        X [c2]
        BuildMultiplyControlledNOT [flag;c1;c2] [anc] [tmplog2space.[(curlog2half+1)%(n/2)]]
        X [c2]
        CCNOT [flag;c2;anc]

    // handle even #bits to the left of the binary point
    if pointpos%2 = 0 then
        X [flag]
        let curlog2half = pointpos/2
        let idx = n-pointpos-curlog2half
        let yadder = slice yhalf [(3*curlog2half)..(n-2)]
        let yadder = List.concat [yadder;extraspace.[0..3*curlog2half-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdder x0.[0..(n-2)] yadder [flag]
        let mutable k = 0
        while idx-k >= 0 do
            if (bigint (constant_left*double (1I<<<k)))&&&1I = 1I && idx-k < n then
                CNOT [flag; x0.[idx-k]]
            k <- k + 1
        CNOT [ys.[n-1]; flag]
        X [flag]

    X [flag]



// Performs initial Newton step for inverse square root, including finding 2^(log2(y)/2)
// Requires pointposition (fixed point representation), y = input, x0 = output (initial guess)
// and 2n+1 extraspace which should be available anyway, since we'll do newton afterwards

// TODO: Optimize such that it only uses n+n/2 ancilla qubits, see old version above
let ApproximateInverseSqrt2 (pointposy:int) (ys:Qubits) (pointposx:int) (x0:Qubits) (space:Qubits) = 
    let n = ys.Length
    let flag = space.[space.Length-1]
    let log2half = slice space [0..n-1]
    let paddingspace = slice space [n..2*n-1]

    let pposdelta = pointposx - pointposy
    let constant_left = 1.613
    let constant_middle = 1.5
    let constant_right = 1.62

    let get_constant_bit (i:int) (bitnum:int) = 
        let mutable ret = 0.
        if (pointposy - i) <= -1 then
            ret <- constant_left
        elif (pointposy - i) < 1 then
            ret <- constant_middle
        else
            ret <- constant_right
        bigint (ret*double (1I<<<bitnum))&&&1I

    // set flag to 1, indicating that we haven't found the first 1 yet:
    X [flag]

    // determine (integer rounded) exponent E of x = 2^-log2(y)/2
    // and calculate first iterate x0 = x * ( 1.5 - y*x^2 / 2 )
    // 1. Use CNOTs to initialize x0 to 1.5 * 2^E (actually it's not 1.5 but a constant which depends on E)
    // 2. Subtract x^3 / 2 * y = 2^{3E-1} * y
    for i in 0..n-1 do
        let idx = n - (i+1) // start from very left and move right
        let neglog2half = -int (floor (double (pointposy - i)/2.)) // rounded -log2y / 2 if current digits are the first 1
        
        let first_one = log2half.[i] // to store the intermediate result (i.e. 2^round(-log2(y)/2))
        // if flag is zero, we are already done, so don't do anything, else: check if current digit is 1
        CCNOT [flag; ys.[idx]; first_one]
        CNOT [first_one; flag]

        let shifted_pposx = n - pointposx + neglog2half // index in x0 where 1.5*x = 1.5*2^neglog2half goes
        let mutable k = 0
        while shifted_pposx - k >= 0 do // set to 1.5 if first 1. Actually, we use a different constant --> converges better
            if shifted_pposx - k < n && (get_constant_bit i k) = 1I then
                CNOT [first_one; x0.[shifted_pposx - k]]
            k <- k + 1

        // now, subtract x * (y * x^2/2) = y * 2^{-3*log2y/2 - 1}
        let mutable y_takhadder = []
        let addition_shift = -(3*neglog2half-1-pposdelta)
        if addition_shift <= 0 then
            y_takhadder <- List.concat [paddingspace.[0..(min (-addition_shift-1) (n-2))];ys.[0..(n-2+addition_shift)];[x0.[n-1]]]
        else
            y_takhadder <- List.concat [ys.[addition_shift..(n-2)];paddingspace.[0..addition_shift-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdderInverse x0.[0..(n-2)] y_takhadder [first_one]

    // flag is still 1 --> all-zero number goes to all-ones
    let lim = min (n-2) (n-pointposx+(n-pointposy)/2)
    for i in (n-pointposx)..lim do
        CNOT [flag;x0.[i]]
        
    // uncompute temporary space & flag
    for i in n-1..(-1)..0 do
        let idx = n - (i+1) // start from very right and move left
        let neglog2half = -int (floor (double (pointposy - i)/2.)) // rounded -log2y / 2 if current digits are the first 1
        
        let first_one = log2half.[i] // result to reset
        CNOT [first_one; flag]
        CCNOT [flag; ys.[idx]; first_one]
    X [flag]

// Inverse of the above
// TODO: Optimize both such that they only use n+n/2 ancilla qubits, see old versions above
let ApproximateInverseSqrt2Inverse (pointposy:int) (ys:Qubits) (pointposx:int) (x0:Qubits) (space:Qubits) = 
    let n = ys.Length
    let flag = space.[space.Length-1]
    let log2half = slice space [0..n-1]
    let paddingspace = slice space [n..2*n-1]

    let pposdelta = pointposx - pointposy
    let constant_left = 1.613
    let constant_middle = 1.5
    let constant_right = 1.62

    let get_constant_bit (i:int) (bitnum:int) = 
        let mutable ret = 0.
        if (pointposy - i) <= -1 then
            ret <- constant_left
        elif (pointposy - i) < 1 then
            ret <- constant_middle
        else
            ret <- constant_right
        bigint (ret*double (1I<<<bitnum))&&&1I

    // set flag to 1, indicating that we haven't found the first 1 yet:
    X [flag]

    // recompute temporary space & flag
    for i in 0..n-1 do
        let idx = n - (i+1) // start from very right and move left
        let neglog2half = -int (floor (double (pointposy - i)/2.)) // rounded -log2y / 2 if current digits are the first 1

        let first_one = log2half.[i] // result to reset
        
        CCNOT [flag; ys.[idx]; first_one]
        CNOT [first_one; flag]

    // flag is still 1 --> all-zero number goes to all-ones
    let lim = min (n-2) (n-pointposx+(n-pointposy)/2)
    for i in (n-pointposx)..lim do
        CNOT [flag;x0.[i]]

    for i in n-1..(-1)..0 do
        let idx = n - (i+1) // start from very left and move right
        let neglog2half = -int (floor (double (pointposy - i)/2.)) // rounded -log2y / 2 if current digits are the first 1

        let first_one = log2half.[i] // to store the intermediate result (i.e. 2^round(-log2(y)/2))
        // add x * (y * x^2/2) = y * 2^{-3*log2y/2 - 1}
        let mutable y_takhadder = []
        let addition_shift = -(3*neglog2half-1-pposdelta)
        if addition_shift <= 0 then
            y_takhadder <- List.concat [paddingspace.[0..(min (-addition_shift-1) (n-2))];ys.[0..(n-2+addition_shift)];[x0.[n-1]]]
        else
            y_takhadder <- List.concat [ys.[addition_shift..(n-2)];paddingspace.[0..addition_shift-1];[x0.[n-1]]]
        BuildCtrlTakahashiAdder x0.[0..(n-2)] y_takhadder [first_one]
        let shifted_pposx = n - pointposx + neglog2half // index in x0 where 1.5*x = 1.5*2^neglog2half goes
        let mutable k = 0
        while shifted_pposx - k >= 0 do // set to 1.5 if first 1. Actually, we use a different constant --> converges better
            if shifted_pposx - k < n && (get_constant_bit i k) = 1I then
                CNOT [first_one; x0.[shifted_pposx - k]]
            k <- k + 1
        CNOT [first_one; flag]
        CCNOT [flag; ys.[idx]; first_one]
    
    X [flag]

// Fast Inverse Square Root / Square root (depending on flag)
// Takes |y>|0.....0> to |y>|1/sqrt(y)>|0....0> or |y>|sqrt(y)>|0...0>
// Uses an approximate 2^(log2(y)/2) and first inline newton iteration as an initial guess
// Then performs the usual Newton iteration to improve the error
// Requires k + 5 registers for k Newton iterations (for temporary computation & uncompute) plus 1 qubit (for squaring)
// Result will be in space.[0..n-1]
let FastInverseSqrt (pointposy:int) (pointposx:int) (numit:int) (ys:Qubits) (space:Qubits) (sqrt:bool) =
    let n = ys.Length
    let x0 = slice space [0..n-1] // space for our initial guess

    // Get initial guess: Performs first Newton iteration on x0=2^log2(y)
    // Using bit-shifts and additions only
    ApproximateInverseSqrt2 pointposy ys pointposx x0 space.[n..4*n-1]
    for i in 1..numit do
        let last_x = space.[(i-1)*n..i*n-1]
        let next_x = space.[(i*n)..(i+1)*n-1]
        let tmpspace = space.[((i+1)*n)..((i+2)*n-1)]
        let tmpspace2 = space.[((i+3)*n)..((i+4)*n-1)]

        let xsquared = space.[((i+2)*n)..((i+3)*n-1)]

        FixedPointSquare pointposx last_x [space.[space.Length-1]] xsquared []

        // y*x^2/2
        FixedPointMultiply (pointposy-1) ys xsquared tmpspace

        // set to 1.5
        X [tmpspace2.[n-pointposx]]
        X [tmpspace2.[n-pointposx-1]]
        BuildTakahashiModAdderInverse tmpspace2 tmpspace

        // multiply by x (into new register)
        FixedPointMultiply (pointposx) tmpspace2 last_x next_x

        BuildTakahashiModAdder tmpspace2 tmpspace
        X [tmpspace2.[n-pointposx]]
        X [tmpspace2.[n-pointposx-1]]

        // uncompute tmpspace
        FixedPointMultiplyInverse (pointposy-1) ys xsquared tmpspace

        // uncompute x^2
        FixedPointSquareInverse pointposx last_x [space.[space.Length-1]] xsquared []

    if numit > 0 then // copy out result and uncompute
        let it = numit
        let res = slice space [(it+4)*n..(it+5)*n-1]
        let lastx = slice space [it*n..(it+1)*n-1]
        if sqrt = false then
            for i in 0..n-1 do
                CNOT [lastx.[i];res.[i]]
        else
            FixedPointMultiply pointposx lastx ys res // to calculate sqrt(y), simply multiply by y
        for i in it..(-1)..1 do
            let last_x = space.[(i-1)*n..i*n-1]
            let next_x = space.[(i*n)..(i+1)*n-1]
            let tmpspace = space.[((i+1)*n)..((i+2)*n-1)]
            let xsquared = space.[((i+2)*n)..((i+3)*n-1)]
            let tmpspace2 = space.[((i+3)*n)..((i+4)*n-1)]

            FixedPointSquare pointposx last_x [space.[space.Length-1]] xsquared []

            // y*x^2/2
            FixedPointMultiply (pointposy-1) ys xsquared tmpspace

            // set to 1.5
            X [tmpspace2.[n-pointposx]]
            X [tmpspace2.[n-pointposx-1]]
            BuildTakahashiModAdderInverse tmpspace2 tmpspace

            // multiply by x (into new register)
            FixedPointMultiplyInverse (pointposx) tmpspace2 last_x next_x

            BuildTakahashiModAdder tmpspace2 tmpspace
            X [tmpspace2.[n-pointposx]]
            X [tmpspace2.[n-pointposx-1]]

            // uncompute tmpspace
            FixedPointMultiplyInverse (pointposy-1) ys xsquared tmpspace

            // uncompute x^2
            FixedPointSquareInverse pointposx last_x [space.[space.Length-1]] xsquared []
        ApproximateInverseSqrt2Inverse pointposy ys pointposx x0 space.[n..4*n-1]
        CSWAP res x0 []
    
// pointposy: point position for input
// pointposx: point position to use for iterates
// reciprocal: set to true for 1/sqrt(x)
// NOTE: For sqrt(x), setting the iterate point position helps save qubits, since 1/sqrt(x) will be multiplied
// by the original input anyway --> one can calculate n MSBs of 1/sqrt(x) instead of having a balanced fixed-point
// number, i.e. pointpos ~ n/2
let FastSqrtFunctionsCompile (reciprocal:bool) (pointposy:int) (pointposx:int) (numit:int) (n:int) (space:Qubits) = 
    FastInverseSqrt pointposy pointposx numit space.[0..n-1] space.[n..space.Length-1] (not reciprocal)
    ()

(*
    Tests for 1/sqrt(x) and sqrt(x)
    N: #numbers to test
    numit: #newton iterations
    n: #bits per number
    pointpos: fixed point position to use
    do_sqrt: 0/1, whether to do sqrt. if 0, it tests 1/sqrt(x) instead

*)
[<LQD>]
let __TestFastInverseSqrt (N:int) (numit:int) (n:int) (pointpos:int) (do_sqrt:bool) =
    let k = Ket((5+numit+1)*n+1)
    let qs = k.Qubits
    let pointposition = pointpos
    let mutable pointposition2 = pointpos
    if do_sqrt then
        pointposition2 <- pointpos + 5
    let circuit = Circuit.Compile (FastSqrtFunctionsCompile (not do_sqrt) pointposition (pointposition2) numit n) qs
                    |> MyCircuitExport
    for i in 1..N do
        let x = 5.*(double (i))/double N
        let k = Ket(qs.Length)
        let qs = k.Qubits
        // prepare x
        (bigint (x*double (1I<<<(n-pointposition)))) |> BoolInt n |> PrepBool qs 0 1 |> ignore

        // tweak the circuit into a Toffoli network that can be simulated efficienly    
        let initialState = Array.zeroCreate (qs.Length)
        List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
        for i in 0..(n-1) do 
            initialState.[i] <- qs.[i].Bit.v
        for i in n..(initialState.Length-1) do
            initialState.[i] <- 0


        //show "Calculating P(x) = %A" expres
        let finalState = MyCircuitSimulateFast circuit initialState
    
        let res = [for i in (n)..(2*n-1) do yield (bigint finalState.[i] <<< (i-n))]
                    |> List.sum
        let rest = [for i in (2*n)..(finalState.Length-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                    |> List.sum
        let pointposx = pointposition2
        let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
        if do_sqrt = false then
            let res = ((double res)/(double (1I<<<(n-pointposx))))
            printfn "%1.15f\t%1.15f\t%A\t %A Toff:%A" x (abs (res-1./sqrt x)) rest res toffcount
        else
            let res = ((double res)/(double (1I<<<(n-pointpos))))
            printfn "%1.15f\t%1.15f\t%A\t %A Toff:%A" x (abs (res-sqrt x)) rest res toffcount
            
    (*
    show "Result: %A" res
    show "1/sqrt(x) = %A" (1./sqrt x)
    show "\nDifference : %A" (abs (res-(1./sqrt x)))
    show "\nRest : %A" rest
    *)
    //show "Discrep: %A, Expected: %A" (abs (res-expres)) ((double n)*f)
    ()

// Calculates arcsin(x)
// Input: 
//      pointpos: pointposition of fixedpoint x
//      x:        n qubits, number in [-1,1]
// Output: y, n qubits
// Temporary space: space
let ArcSine (pointpos:int) (numit:int) (x:Qubits) (y:Qubits) (space:Qubits) =
    let a = space.[space.Length-1]
    let signbit = space.[space.Length-2]
    let n = x.Length
    
    let itcount = numit // 1/sqrt(x) for arcsin

    if pointpos > 1 then // conditionally invert x
        CNOT [x.[x.Length-1]; signbit]
        for i in 0..n-1 do
            CNOT [signbit; x.[i]]
        CNOT [signbit; space.[0]]
        BuildTakahashiModAdder x.[0..n-2] space.[0..n-2]
        CNOT [signbit; space.[0]]

    // compute condition: 'x < 0.5' --> a
    X [x.[n-pointpos]]
    X [x.[n-pointpos-1]]
    CCNOT [x.[n-pointpos]; x.[n-pointpos-1]; a] // <=> both bits were 0 (prior to the two X-gates) <=> x < .5
    X [x.[n-pointpos]]
    X [x.[n-pointpos-1]]

    // compute z = .5 - (x>>1)
    let z = space.[0..n-1]
    X [z.[n-pointpos-1]] // = 0.5
    BuildTakahashiAdderInverse z.[0..n-2] (List.concat [x.[1..n-1]; [z.[n-1]]])


    // conditional copy of .5(1-x) into new register. conditional square x
    let polyinput = space.[n..2*n-1]
    X [a]
    for i in 0..n-1 do
        CCNOT [z.[i];a;polyinput.[i]]
    X [a]
    FixedPointSquare pointpos x [space.[2*n]] polyinput [a]

    
    // run polyval to get arcsin(z) (including uncompute)
    // needs n*(d+1)+1 qubits in space, first n hold x^2, last n+1 consist of result (first) and then a zero ancilla
    let arcsin_coeffs3 = [0.0683039213468; 0.0700620936908; 0.167031047683; 0.999992844795; ]
    let arcsin_coeffs4 = [0.0533211503668; 0.0381372893277; 0.0757617508325; 0.166631017736; 1.00000046342; ]
    let arcsin_coeffs5 = [0.0450602470503; 0.0221860463684; 0.0459698365257; 0.0748982216286; 0.166669973683; 0.999999969443; ]
    let arcsin_coeffs6 = [0.0401228098468; 0.012316781492; 0.0324722511722; 0.0444170327186; 0.0750123462156; 0.166666371194; 1.00000000204; ]
    let arcsin_coeffs7 = [0.0370749359021; 0.00522825693564; 0.0254620155695; 0.0299484109714; 0.0446768005991; 0.0749986030971; 0.166666692354; 0.999999999863; ]
    let arcsin_coeffs8 = [0.0353759984124; -0.000518617087625; 0.0216788087785; 0.0216420468609; 0.0304536885976; 0.044638873558; 0.0750001138673; 0.166666665299; 1.0; ]
    let mutable arcsin_coeffs = []
    if numit <= 3 then
        arcsin_coeffs <- arcsin_coeffs3
    elif numit <=4 then
        arcsin_coeffs <- arcsin_coeffs6
    else
        arcsin_coeffs <- arcsin_coeffs8
    let polyspace = List.concat [polyinput;space.[2*n..(1+arcsin_coeffs.Length)*n]]
    EvaluatePolynomial pointpos n arcsin_coeffs polyspace

    
    // uncompute polyinput
    FixedPointSquareInverse pointpos x [space.[2*n]] polyinput [a]
    X [a]
    for i in 0..n-1 do
        CCNOT [z.[i];a;polyinput.[i]]
    X [a]

    let emptyspace = List.concat [polyinput;space.[2*n..arcsin_coeffs.Length*n-1];space.[(arcsin_coeffs.Length+1)*n..(7+itcount)*n]] // n 4n x*n => x = itcount
    let arcsin = space.[arcsin_coeffs.Length*n..(arcsin_coeffs.Length+1)*n-1]
    
    // run sqrt(z)
    let sqrtspace = emptyspace.[0..(5+itcount)*n]
    FastInverseSqrt pointpos (pointpos+11) itcount z sqrtspace true

    // conditional copy of arcsin(z) into new register // OR // init new register to pi/2 and then subtract (arcsin(z)<<1)
    let emptyspace = emptyspace.[n..(itcount+5)*n]
    let sqrt = sqrtspace.[0..n-1]    
    let tmpout = emptyspace.[0..n-1]
    let output = y
    CSWAP x sqrt [a] // make choice using swaps --> we now have to multiply by and add either sqrt((1-x)/2) or x
    FixedPointMultiply pointpos sqrt arcsin tmpout // multiply polynomial by sqrt or x, depending on value of a
    //BuildTakahashiModAdder tmpout sqrt // add sqrt(...) or x
    CSWAP x sqrt [a]

    for i in 0..n-1 do
        CCNOT [a; tmpout.[i]; output.[i]]
    X [a]
    // set to pi/2 and subtract tmpout
    let intpiovertwo = bigint (1.57079632679*double (1I<<<(n-pointpos)))
    for i in 0..n-1 do
        if (intpiovertwo>>>i)&&&1I = 1I then
            CNOT [a; output.[i]]
    BuildCtrlTakahashiModAdderInverse output (List.concat [[emptyspace.[n]];tmpout.[0..n-2]]) [a]
    X [a]
    // invert the result conditioned on signbit
    // and undo inversion of x
    if pointpos > 1 then // conditionally invert x
        let emptyspace = emptyspace.[n..2*n-1]
        
        for i in 0..n-1 do
            CNOT [signbit; output.[i]]
            CNOT [signbit; x.[i]]
        CNOT [signbit; emptyspace.[0]]
        BuildTakahashiModAdder output emptyspace.[0..n-1]
        BuildTakahashiModAdder x.[0..n-1] emptyspace.[0..n-1]
        CNOT [signbit; emptyspace.[0]]

        CNOT [x.[n-1];signbit]
    

let ArcSineCompile (pointpos:int) (n:int) (numit:int) (qs:Qubits) =
    ArcSine pointpos numit qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..qs.Length-1]

[<LQD>]
let __TestArcSine (N:int) (numit:int) (n:int) (pointpos:int) (xstart:double) =
    //let n = 40
    //let pointpos = 15
    //let numit = 4 // # newton iterations for 1/sqrt(x)
    let k = Ket((9+numit)*n+3)
    let qs = k.Qubits
    let pointposition = pointpos

    let circuit = Circuit.Compile (ArcSineCompile pointposition n numit) qs
                    |> MyCircuitExport
    for i in 1..N do
        let x = xstart+(double i)/double N
        let k = Ket(qs.Length)
        let qs = k.Qubits
        // prepare x
        (bigint (x*double (1I<<<(n-pointposition)))) |> BoolInt n |> PrepBool qs 0 1 |> ignore

        // tweak the circuit into a Toffoli network that can be simulated efficienly    
        let initialState = Array.zeroCreate (qs.Length)
        List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
        for i in 0..(n-1) do 
            initialState.[i] <- qs.[i].Bit.v
        for i in n..(initialState.Length-1) do
            initialState.[i] <- 0


        //show "Calculating P(x) = %A" expres
        let finalState = MyCircuitSimulateFast circuit initialState
    
        let res = [for i in (n)..(2*n-1) do yield (bigint finalState.[i] <<< (i-n))]
                    |> List.sum
        let rest = [for i in (2*n)..(finalState.Length-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                    |> List.sum
        let res = ((double res)/(double (1I<<<(n-pointpos))))
        let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
        printfn "%1.15f\t%1.15f\t%A\t%A" x res (abs (res-asin(x))) toffcount


let Tanh (pointpos:int) (n:int) (deg:int) (x:Qubits) (y:Qubits) (space:Qubits) =
    let coeffs = [  
                    // deg=4
                    [
                    [0.00836605407729; -0.0433001854593; 0.129220873061; -0.332658282587; 0.999968474111; ];
                    [0.0017880958143; -0.0189714446127; 0.0934908201968; -0.307901330138; 0.993184739898; ];
                    [0.000377713402567; -0.00683523569103; 0.0535129186003; -0.248073286098; 0.958862596326; ];
                    [8.18608713022e-05; -0.00230062776775; 0.0271461608152; -0.179084801726; 0.89029684561; ];
                    [1.97004975078e-05; -0.000800827841966; 0.0134744060414; -0.123253816885; 0.804091464055; ];
                    [3.95000617653e-06; -0.000238110857742; 0.00587278260311; -0.077225073229; 0.698679945669; ];
                    [7.35026439544e-07; -6.52562921174e-05; 0.00236309835924; -0.0453106525641; 0.588968014312; ];
                    [1.37218178417e-07; -1.78307523866e-05; 0.000942311097591; -0.0262540976327; 0.492394339598; ];
                    ];
                    // deg=5
                    [
                    [-0.00283958763106; 0.0161145078516; -0.0509731825118; 0.132535390551; -0.333240063209; 0.999996855378; ];
                    [-0.000439464194481; 0.00537566547671; -0.030472522365; 0.111568948124; -0.321821365938; 0.997382643292; ];
                    [-6.79322708042e-05; 0.00144106041681; -0.0134281696714; 0.0737460188351; -0.278800647088; 0.977333964545; ];
                    [-1.08292315323e-05; 0.000359584137351; -0.00513247129476; 0.0414953931822; -0.215213694722; 0.926455932594; ];
                    [-1.95543362225e-06; 9.42815969379e-05; -0.00193403082968; 0.0220480896915; -0.155553055175; 0.852560429083; ];
                    [-2.82512911715e-07; 2.02458643595e-05; -0.00061183916508; 0.0101320559475; -0.101345274572; 0.752975409793; ];
                    [-3.71134142181e-08; 3.91935607024e-06; -0.000174098499311; 0.00421559217064; -0.0610099744809; 0.641965057256; ];
                    [-4.85194440293e-09; 7.49902564197e-07; -4.86261036593e-05; 0.00171239627242; -0.0358343102012; 0.539826458471; ];
                    ];
                    ]
                    
    let square_ancilla = space.[space.Length-1]
    // Intervals: (0,1),(1,1.5),(1.5,2),(2,2.5),(2.5,3),(3,3.75),(3.75,4.5),(4.5,5.5)

    FixedPointSquare pointpos x [square_ancilla] space.[4..n+3] [] // x^2 into first register of polynomial eval space
    // requires n*(d+1)+1 space: first n are x, last n+1 are result and 1 ancilla
    if space.Length-2 < n*(deg-1)+3 then
        failwith "Gaussian needs more space! Requires n*deg+3 qubits (besides the input and output (n each))."
    
    // initialize label:
    // 1.5 > x >= 1
    for i in (n-pointpos+1)..n-1 do
        X [x.[i]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 2 > x >= 1.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    // 2.5 > x >= 2
    for i in (n-pointpos-1)..(n-pointpos+1) do
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 3 > x >= 2.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    // 3.5 > x >= 3 or 3.75 > x >= 3.5
    X [x.[n-pointpos]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos-1]]
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[0]] [square_ancilla]
    
    // 4.5 > x >= 3.75
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    for i in (n-pointpos-1)..(n-pointpos+2) do // >= 4, < 4.5
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]

    // > 4.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos]]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos+1]]
    for i in (n-pointpos+3)..n-1 do
        X [x.[i]]
        
    EvaluatePolynomialsInParallel pointpos n coeffs.[deg-4] (List.concat [space.[0..3];space.[4..space.Length-2]; y; [square_ancilla]])
    
    for i in (n-pointpos+1)..n-1 do
        X [x.[i]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 2 > x >= 1.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    // 2.5 > x >= 2
    for i in (n-pointpos-1)..(n-pointpos+1) do
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 3 > x >= 2.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    // 3.5 > x >= 3 or 3.75 > x >= 3.5
    X [x.[n-pointpos]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos-1]]
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[0]] [square_ancilla]
    
    // 4.5 > x >= 3.75
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    for i in (n-pointpos-1)..(n-pointpos+2) do // >= 4, < 4.5
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]

    // > 4.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos]]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[1]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[0]] [square_ancilla]
    X [x.[n-pointpos+1]]
    for i in (n-pointpos+3)..n-1 do
        X [x.[i]]

    FixedPointSquareInverse pointpos x [square_ancilla] space.[4..n+3] [] // undo x^2 into first register of polynomial eval space

    FixedPointMultiply pointpos x y space.[0..n-1]
    CSWAP y space.[0..n-1] []

let TanhCompile (pointpos:int) (n:int) (deg:int) (qs:Qubits) =
    Tanh pointpos n deg qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..qs.Length-1]

let Gaussian (pointpos:int) (n:int) (deg:int) (x:Qubits) (y:Qubits) (space:Qubits) =
    let coeffs = [  
                    // deg=4
                    [
                    [0.0256435188419; -0.153964704725; 0.495703837424; -0.999483216799; 0.999990010674; ];
                    [0.008393327869; -0.0881842716417; 0.395338064465; -0.927159957139; 0.979482306457; ];
                    [0.00295372568759; -0.0432063601294; 0.25436104732; -0.728380662789; 0.872999192856; ];
                    [0.00123522568151; -0.0223929705091; 0.159440520264; -0.535176954357; 0.724907819034; ];
                    [0.000266607674502; -0.00653738958273; 0.0614585812752; -0.264298605099; 0.442308192305; ];
                    [2.26830754491e-05; -0.000783258225071; 0.0102461400599; -0.060343179529; 0.13544763756; ];
                    [3.04237123363e-07; -1.64801739133e-05; 0.000333719872925; -0.00299633043656; 0.0100752205571; ];
                    ];
                    // deg=5
                    [
                    [-0.00511499867207; 0.0383771404112; -0.165046205896; 0.49963285863; -0.999969353879; 0.999999584612; ];
                    [-0.00167170582074; 0.0219488103723; -0.131335237322; 0.462673086865; -0.97862587932; 0.994889848637; ];
                    [-0.000589703294161; 0.0107816785273; -0.0846492354622; 0.363741071845; -0.872295513351; 0.948515272573; ];
                    [-0.000246466514465; 0.00558464614435; -0.0530271425589; 0.267083652076; -0.723875308925; 0.85692784939; ];
                    [-5.2614331278e-05; 0.0016120965901; -0.020217212838; 0.130573555437; -0.437809607269; 0.61545859969; ];
                    [-4.4478088487e-06; 0.000191909540176; -0.00334817392589; 0.0296037738021; -0.133086731209; 0.244337121657; ];
                    [-5.42975342326e-08; 3.67197872567e-06; -9.91975542974e-05; 0.00133917978519; -0.00904395193377; 0.0244743120838; ];
                    ];
                    // deg=6
                    [
                    [0.000850906851876; -0.00766002729018; 0.0412297404529; -0.166518502838; 0.499975874825; -0.999998521177; 0.999999985183; ];
                    [0.000277808112028; -0.00437645371661; 0.0327583132698; -0.154018232118; 0.489017696333; -0.994676121135; 0.998896407364; ];
                    [9.81635503588e-05; -0.00215360860852; 0.0211388606782; -0.121145364401; 0.435909063948; -0.94822454535; 0.981721916152; ];
                    [4.10108299005e-05; -0.00111506122428; 0.0132363477274; -0.0889133202485; 0.361586661604; -0.856366703068; 0.934185678037; ];
                    [8.68793974356e-06; -0.000319378663661; 0.00500852948112; -0.0431672867668; 0.217372245033; -0.612019397686; 0.760418021443; ];
                    [7.31134372592e-07; -3.78482551648e-05; 0.000825599400538; -0.00973916504114; 0.0657398761848; -0.241695927482; 0.379897286383; ];
                    [8.33664517646e-09; -6.76088188369e-07; 2.28418695083e-05; -0.000411847382875; 0.00418385541378; -0.0227335185041; 0.0516958073412; ];
                    ];
                    ]
                    
    let square_ancilla = space.[space.Length-1]

    FixedPointSquare pointpos x [square_ancilla] space.[4..n+3] [] // x^2 into first register of polynomial eval space
    // requires n*(d+1)+1 space: first n are x, last n+1 are result and 1 ancilla
    if space.Length-2 < n*(deg-1)+3 then
        failwith "Gaussian needs more space! Requires n*deg+3 qubits (besides the input and output (n each))."
    
    // initialize label:
    // 1.5 > x >= 1
    for i in (n-pointpos+1)..n-1 do
        X [x.[i]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 1.75 > x >= 1.5
    X [x.[n-pointpos-1]]
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    // 2 > x >= 1.75
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[0]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    // 2.5 > x >= 2
    for i in (n-pointpos-1)..(n-pointpos+1) do
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    // 3 > x >= 2.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 4 > x >= 3
    X [x.[n-pointpos]]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[1]] [square_ancilla]
    for i in (n-pointpos+2)..n-1 do
        X [x.[i]]

    EvaluatePolynomialsInParallel pointpos n coeffs.[deg-4] (List.concat [space.[0..3];space.[4..space.Length-2]; y; [square_ancilla]])

    // initialize label:
    // 1.5 > x >= 1
    for i in (n-pointpos+1)..n-1 do
        X [x.[i]]
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 1.75 > x >= 1.5
    X [x.[n-pointpos-1]]
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    // 2 > x >= 1.75
    X [x.[n-pointpos-2]]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[0]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-2..n-1] [space.[1]] [square_ancilla]
    // 2.5 > x >= 2
    for i in (n-pointpos-1)..(n-pointpos+1) do
        X [x.[i]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    // 3 > x >= 2.5
    X [x.[n-pointpos-1]]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos-1..n-1] [space.[0]] [square_ancilla]
    // 4 > x >= 3
    X [x.[n-pointpos]]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[2]] [square_ancilla]
    BuildMultiplyControlledNOT x.[n-pointpos..n-1] [space.[1]] [square_ancilla]
    for i in (n-pointpos+2)..n-1 do
        X [x.[i]]

    FixedPointSquareInverse pointpos x [square_ancilla] space.[4..n+3] [] // undo x^2 into first register of polynomial eval space


let GaussianCompile (pointpos:int) (n:int) (deg:int) (qs:Qubits) =
    Gaussian pointpos n deg qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..qs.Length-1]

let Sine (pointpos:int) (n:int) (deg:int) (x:Qubits) (y:Qubits) (space:Qubits) =
    let coeffs = [-2.40335666236e-08; 2.75336986568e-06; -0.000198410305395; 0.00833333232747; -0.166666666531; 0.999999999997; ];
    let deg = coeffs.Length-1
    let square_ancilla = space.[space.Length-1]

    FixedPointSquare pointpos x [square_ancilla] space.[0..n-1] [] // x^2 into first register of polynomial eval space
    // requires n*(d+2)+1 space: first n are x, last n+1 are result and 1 ancilla    
    EvaluatePolynomial pointpos n coeffs (List.concat [space.[0..n*(deg)-1]; y; [square_ancilla]])

    FixedPointSquareInverse pointpos x [square_ancilla] space.[0..n-1] [] // undo x^2 into first register of polynomial eval space
    FixedPointMultiply pointpos x y space.[0..n-1]
    CSWAP y space.[0..n-1] []

let SineCompile (pointpos:int) (n:int) (deg:int) (qs:Qubits) =
    Sine pointpos n deg qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..qs.Length-1]


// Tests for Gaussian and Tanh
[<LQD>]
let __TestGaussianOrTanh (N:int) (deg:int) (n:int) (xstart:double) =
    //let n = 40
    //let pointpos = 15
    //let numit = 4 // # newton iterations for 1/sqrt(x)
    let k = Ket((2+deg)*n+1+4)
    let qs = k.Qubits
    let pointposition = 10

    let Tanh = false//true
    let Sin = true//true

    let mutable func = (GaussianCompile pointposition n deg)
    let mutable TestIntervalSize = 4.
    if Tanh then
        TestIntervalSize <- 5.5
        func <- (TanhCompile pointposition n deg)
    if Sin then
        TestIntervalSize <- Math.PI/2.
        func <- (SineCompile pointposition n deg)

    let circuit = Circuit.Compile func qs
                    |> MyCircuitExport
    for i in 1..N do
        let x = xstart+TestIntervalSize*(double (i-1))/double N
        let k = Ket(qs.Length)
        let qs = k.Qubits
        // prepare x
        (bigint (x*double (1I<<<(n-pointposition)))) |> BoolInt n |> PrepBool qs 0 1 |> ignore

        // tweak the circuit into a Toffoli network that can be simulated efficienly    
        let initialState = Array.zeroCreate (qs.Length)
        List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
        for i in 0..(n-1) do 
            initialState.[i] <- qs.[i].Bit.v
        for i in n..(initialState.Length-1) do
            initialState.[i] <- 0


        //show "Calculating P(x) = %A" expres
        let finalState = MyCircuitSimulateFast circuit initialState
        let x = [for i in (0)..(n-1) do yield (bigint finalState.[i] <<< (i))]
                    |> List.sum
        let res = [for i in (n)..(2*n-1) do yield (bigint finalState.[i] <<< (i-n))]
                    |> List.sum
        let rest = [for i in (2*n)..(finalState.Length-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                    |> List.sum
        let mutable res = ((double res)/(double (1I<<<(n-pointposition))))
        let mutable x = ((double x)/(double (1I<<<(n-pointposition))))
        if res > double (1I<<<(pointposition-1)) then
            res <- double (1I<<<(pointposition))-res
        let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
        if Tanh = false && Sin = false then
            printfn "%1.15f\t%1.15f\t%A\t%A\t%A" x res (abs (res-exp(-x*x))) rest toffcount
        elif Sin = false then
            printfn "%1.15f\t%1.15f\t%A\t%A\t%A" x res (abs (res-tanh(x))) rest toffcount
        else
            printfn "%1.15f\t%1.15f\t%A\t%A\t%A" x res (abs (res-sin(x))) rest toffcount



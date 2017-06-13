module Microsoft.Research.Liquid.Inplace

open System
open System.Collections.Generic
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions

open Montgomery



// Performs integer division, i.e. calculates x / y
// Takes x,y,0 to r,y,q
// where r denotes the remainder and q the quotient
let IntegerDivision (xs:Qubits) (ys:Qubits) (rs:Qubits) =
    let n = xs.Length

    for i in xs.Length-1..(-1)..0 do
        let anc = slice rs [0..i-1]
        let subtr = slice xs [i..xs.Length-1]
        let subtraction_reg = List.concat [subtr;anc]
        BuildTakahashiAdderInverse subtraction_reg (List.concat [ys;[rs.[i]]])
        BuildCtrlTakahashiModAdder subtraction_reg ys [rs.[i]]
        X [rs.[i]]

let IntegerDivisionCompile (qs:Qubits) =
    let n = qs.Length/3
    IntegerDivision qs.[0..n-1] qs.[n..2*n-1] qs.[2*n..3*n-1]
    ()

[<LQD>]
let __TestIntegerDivision (a:int) (b:int) =
    let n = 10 // # bits
    let k = Ket(3*n)
    let qs = k.Qubits

    let circuit = Circuit.Compile IntegerDivisionCompile qs
                    |> MyCircuitExport
    
    let k = Ket(qs.Length)
    let qs = k.Qubits
    // prepare x
    bigint a |> BoolInt n |> PrepBool qs 0 1 |> ignore
    bigint b |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(2*n-1) do 
        initialState.[i] <- qs.[i].Bit.v
    for i in 2*n..(initialState.Length-1) do
        initialState.[i] <- 0

    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in (2*n)..(3*n-1) do yield (bigint finalState.[i] <<< (i-2*n))]
                |> List.sum
    let rem = [for i in (0)..(n-1) do yield (bigint finalState.[i] <<< i)]
                |> List.sum
    printfn "%A / %A = %A, rem %A" a b res rem



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inplace work starts here
// Includes:
//  * Addition of constants in 10nlog(n) using 1 and 8nlog(n) using 2 dirty ancilla qubits
//  * Incrementer (by Gidney) using n dirty ancillae in ~4n Toffolis
//  * Carry computation using n dirty ancillae in ~4n Toffolis
//  * Addition of a constant modulo N in ~twice the cost of addition of constant
//  * Modular multiplication using repeated modular addition and shift including uncomputation, i.e.
//          |x>|0> --> |(ax)mod N>|0>
//      in 32n^2 log(n)
//
//  (+ various functions for testing)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Uses qubits xg of the form [x_0,...,x_n-1,garbage,...,garbage] starting at the LSB
// to add the constant c to x.
// This is an n^2 algorithm, see ConstantAdderInplace for an nlog(n) version
// using only 1 or 2 borrowed dirty qubits.
let ConstantAdderUsingGarbage (c:int) (xg:Qubits) = 
    let n = xg.Length/2
    let x = slice xg [0..(n-1)]
    let g = slice xg [n..(2*n-1)]
    for r in (n-1)..(-1)..0 do
        if r > 0 then
            CNOT [g.[r-1]; x.[r]]
        for i in (r-1)..(-1)..0 do
            if ((c >>> i)&&&1 = 1) then
                CNOT [x.[i]; g.[i]]
                if i > 0 then
                    CNOT [g.[i-1]; g.[i]]
            if i > 0 then
                CCNOT [g.[i-1]; x.[i]; g.[i]]
        for i in 0..(r-2) do
            CCNOT [g.[i]; x.[i+1]; g.[i+1]]
            if ((c >>> (i+1))&&&1 = 1) then
                CNOT [g.[i]; g.[i+1]]
        if r > 0 then
            CNOT [g.[r-1]; x.[r]]
        if (c >>> r)&&&1 = 1 then
            X [x.[r]]

        for i in (r-1)..(-1)..0 do
            if ((c >>> i)&&&1 = 1) then
                CNOT [x.[i]; g.[i]]
                if i > 0 then
                    CNOT [g.[i-1]; g.[i]]
            if i > 0 then
                CCNOT [g.[i-1]; x.[i]; g.[i]]
        for i in 0..(r-2) do
            CCNOT [g.[i]; x.[i+1]; g.[i+1]]
            if ((c >>> (i+1))&&&1 = 1) then
                CNOT [g.[i]; g.[i+1]]

// n-bit incrementer using n borrowed bits
let BorrowedQubitsIncrementer (xg:Qubits) = 
    let n = int (xg.Length/2)
    if n > 3 then
        let g = slice xg [n..(2*n-1)]
        let x = slice xg [0..(n-1)]
        let gx = List.concat [g;x]
        TakahashiModAdderInverse gx
        
        for i in 0..(n-1) do
            X [g.[i]]
        TakahashiModAdderInverse gx
        for i in 0..(n-1) do
            X [g.[i]]

    elif n = 2 then
        CNOT [xg.[0]; xg.[1]]
        X [xg.[0]]
    elif n = 3 then
        CCNOT [xg.[0]; xg.[1]; xg.[2]]
        CNOT [xg.[0]; xg.[1]]
        X [xg.[0]]
    elif n = 1 then
        X [xg.[0]]

// Controlled n-bit incrementer using n borrowed bits
// Last qubit is the control
let CtrlBorrowedQubitsIncrementer (xgc:Qubits) = 
    let n = int ((xgc.Length-1)/2)
    let xg = slice xgc [0..(2*n-1)]
    let c = xgc.[2*n]
    if n > 2 then
        let g = slice xg [n..(2*n-1)]
        let x = slice xg [0..(n-1)]
        let gxc = List.concat [g;x;[c]]
        CtrlTakahashiModAdderInverse gxc
        
        for i in 0..(n-1) do
            X [g.[i]]
        CtrlTakahashiModAdderInverse gxc
        for i in 0..(n-1) do
            X [g.[i]]

    elif n = 2 then
        CCNOT [c; xg.[0]; xg.[1]]
        CNOT [c; xg.[0]]
    elif n = 1 then
        CNOT [c;xg.[0]]

// computes the carry of the computation x += c using dirty qubits g
// at round r (i.e. determines if the r-th bit of a would be flipped due to a carry
// propagation from less-significant bits)
let ComputeCarryUsingGarbage (c:bigint) (g:Qubits) (x:Qubits) (ctrls:Qubits) = 
    let r = x.Length - 1
    let one = bigint 1
    if r = 1 then
        if (c&&&one) = one then
            BuildMultiplyControlledNOT (List.concat [[x.[0]];ctrls]) [x.[1]] [g.[0]]
    else
        // make this controlled on a control-qubit register to have a conditional addition
        BuildMultiplyControlledNOT (List.concat [[g.[r-2]];ctrls]) [x.[r]] [x.[r-1]]
        for i in (r-1)..(-1)..2 do
            if ((c >>> i)&&&one = one) then
                CNOT [x.[i]; g.[i-1]]
                X [x.[i]]
            CCNOT [g.[i-2]; x.[i]; g.[i-1]]

        if ((c>>>1)&&&one) = one then
            CNOT [x.[1];g.[0]]
        if (c&&&one) = one then
            if ((c>>>1)&&&one) = one then
                X [x.[1]]
            CCNOT [x.[0];x.[1];g.[0]]

        for i in 1..(r-2) do
            CCNOT [g.[i-1]; x.[i+1]; g.[i]]

        BuildMultiplyControlledNOT (List.concat [[g.[r-2]];ctrls]) [x.[r]] [x.[r-1]]
        
        for i in (r-2)..(-1)..1 do
            CCNOT [g.[i-1]; x.[i+1]; g.[i]]

        if (c&&&one) = one then
            CCNOT [x.[0];x.[1];g.[0]]
            if ((c>>>1)&&&one) = one then
                X [x.[1]]
        if ((c>>>1)&&&one) = one then
            CNOT [x.[1];g.[0]]
            
        for i in 2..(r-1) do
            CCNOT [g.[i-2]; x.[i]; g.[i-1]]
            if ((c >>> i)&&&one = one) then
                X [x.[i]]
                CNOT [x.[i]; g.[i-1]]
        

// Adds the constant c to the x-register using 1 garbage qubit
// ctrls: Control qubits (if needed)

// Runs in 10 nlog(n) Toffoli gates for 1 garbage qubit
// and in 8 nlog(n) with 2 garbage qubits (chooses automatically)
let ConstantAdderInplace (c:bigint) (x:Qubits) (garbage:Qubits) (ctrls:Qubits) =
    let n = x.Length
    let one = bigint 1
    let g = garbage.[0]

    // Recursively splits the x-register into 2 parts and applies the incrementer/addition
    // decomposition. The dirty qubit (g) is used do hold the carry. Conditioned 
    // on g, the higherbits register is incremented. If two garbage qubits are available,
    // the increment operation is carried out on [[g];higherbits], followed by X [g]
    // which is equivalent to a controlled incrementer but costs less than 2 ctrltakahashi
    // adders.
    let rec CircGen (c:bigint) (x:Qubits) (frombit:int) = 
        let n = x.Length           
        // Handle cases with more than 2 qubits by recursive splitting & incrementing
        // If there were a carry when adding c to x when going from bit L to L+1.
        if n >= 2 then
            let L = x.Length-(int x.Length/2)
            let lowerbits = slice x [0..(L-1)]
            let higherbits = slice x [L..(x.Length-1)]
            
            let incrinput_lowbits = slice lowerbits [0..((higherbits.Length)-1)]
            let incrinput = List.concat [higherbits; incrinput_lowbits; [g]]
            // (conditionally) increment
            if garbage.Length >=2 then
                BorrowedQubitsIncrementer (List.concat [[g];higherbits;incrinput_lowbits;[garbage.[1]]])
                X [g]
            else
                CtrlBorrowedQubitsIncrementer incrinput

            // conditionally invert
            for i in 0..(higherbits.Length-1) do
                CNOT [g; higherbits.[i]]

            // compute carry
            let carryinputx = List.concat [lowerbits; [g]]

            ComputeCarryUsingGarbage (c >>> frombit) higherbits carryinputx ctrls
            
            // (conditionally) in/de-crement
            if garbage.Length >=2 then
                BorrowedQubitsIncrementer (List.concat [[g];higherbits;incrinput_lowbits;[garbage.[1]]])
                X [g]
            else
                CtrlBorrowedQubitsIncrementer incrinput

            // uncompute carry (i.e. repeat circuit from before)
            ComputeCarryUsingGarbage (c >>> frombit) higherbits carryinputx ctrls

            // conditionally invert            
            for i in 0..(higherbits.Length-1) do
                CNOT [g; higherbits.[i]]
            
            // generate circuit for lower half (recursive call)
            CircGen c lowerbits (frombit)
            // same for upper half
            CircGen c higherbits (L+frombit)

    // do recursive splitting and increment
    CircGen c x 0

    // perform addition
    for i in 0..(n-1) do
        // Finally: Do "addition" on 1 bit (carry has been taken care of already)
        if ((c>>>i) &&& one) = one then
            BuildMultiplyControlledNOT ctrls [x.[i]] [g]

// Wrapper for circuit.compile
let ConstantAdderInplaceCompile (numctrls:int) (c:bigint) (xgc:Qubits) =
    let ctrls = slice xgc [(xgc.Length-numctrls)..(xgc.Length-1)]
    let xg = slice xgc [0..(xgc.Length-numctrls-1)]
    let x = slice xg [0..(xg.Length-3)]
    let g = xg.[xg.Length-2]
    let g2 = xg.[xg.Length-1]
    ConstantAdderInplace c x [g;g2] ctrls

// Wrapper for circuit.compile
let BuildCarryUsingGarbage (c:bigint) (numctrls:int) (qs:Qubits) =
    let ctrls = slice qs [(qs.Length-numctrls)..(qs.Length-1)]
    let n = int (((qs.Length-numctrls)+2)/2)
    let x = slice qs [(n-2)..(2*n-3)]
    let g = slice qs [0..(n-3)]

    ComputeCarryUsingGarbage c g x ctrls

// Adds the constant c to the register x modulo N
// xgc consists of a qubit in |0> (labelled z), the x register, n-1 dirty (borrowed) qubits, and control qubits (# control qubits = numctrls)
// It uses the Takahashi trick for modular reduction, i.e. computes the final carry and, conditioned on that, either
// adds a or adds a-N
let ConstantAdderInplaceModN (numctrls:int) (c:bigint) (N:bigint) (zxgc:Qubits) =
    let n = int ((zxgc.Length-numctrls)/2)
    let x = slice zxgc [1..n]
    let g = slice zxgc [(n+1)..(2*n-1)]
    let ctrls = slice zxgc [(2*n)..(2*n+numctrls-1)]
    let z = [zxgc.[0]]
    let one = bigint 1
    let xz = List.concat [x;z]

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage (N-c) g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

    ConstantAdderInplace c x [g.[0];g.[1]] z
    BuildMultiplyControlledNOT ctrls z [g.[0]]
    
    for i in 0..(n-1) do
        X [x.[i]]
    ConstantAdderInplace (N-c) x [g.[0];g.[1]] z
    ComputeCarryUsingGarbage c g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

// Subtracts the constant c from the register x modulo N
// xgc consists of a qubits in |0> (labelled z), the x register, n-1 dirty (borrowed) qubits, and control qubits (# control qubits = numctrls)
// It uses the Takahashi trick for modular reduction, i.e. computes the final carry and, conditioned on that, either
// adds a or adds a-N
let ConstantSubtractorInplaceModN (numctrls:int) (c:bigint) (N:bigint) (zxgc:Qubits) =
    let n = int ((zxgc.Length-numctrls)/2)
    let x = slice zxgc [1..n]
    let g = slice zxgc [(n+1)..(2*n-1)]
    let ctrls = slice zxgc [(2*n)..(2*n+numctrls-1)]
    let z = [zxgc.[0]]
    let one = bigint 1
    let xz = List.concat [x;z]

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage c g xz ctrls
    ConstantAdderInplace ((one<<<n)-(N-c)) x [g.[0];g.[1]] z
    for i in 0..(n-1) do
        X [x.[i]]

    BuildMultiplyControlledNOT ctrls z [g.[0]]
    ConstantAdderInplace ((one<<<n)-c) x [g.[0];g.[1]] z

    for i in 0..(n-1) do
        X [x.[i]]
    ComputeCarryUsingGarbage (N-c) g xz ctrls
    for i in 0..(n-1) do
        X [x.[i]]

// Controlled SWAP operation used by the constant modular multiplier
let CSWAP (x:Qubits) (y:Qubits) (ctrls:Qubits) =
    let n = x.Length
    // perform swap between x and y conditioned on ctrl
    for i in 0..(n-1) do
        BuildMultiplyControlledNOT [x.[i]] [y.[i]] [x.[(i+1)%n]]
        BuildMultiplyControlledNOT (List.concat [[y.[i]];ctrls]) [x.[i]] [y.[(i+1)%n]]
        BuildMultiplyControlledNOT [x.[i]] [y.[i]] [x.[(i+1)%n]]

let CSWAPCompile (numctrls:int) (xyc:Qubits) =
    let n = int ((xyc.Length-numctrls)/2)
    let x = slice xyc [0..(n-1)]
    let y = slice xyc [n..(2*n-1)]
    let c = slice xyc [(2*n)..(2*n+numctrls-1)]
    CSWAP x y c

// (Controlled) Modular multiplication of a quantum number x by a classical constant c
// Transforms |x>|y=0>|0>|ctrl> to |c*x mod N>|y=0>|0>|ctrl> if ctrl = 1...1
let ConstantMultiplierModN (c:bigint) (cinvModN:bigint) (N:bigint) (x:Qubits) (y:Qubits) (zero:Qubits) (ctrls:Qubits) =
    let n = x.Length
    let one = bigint 1
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        ConstantAdderInplaceModN (1+ctrls.Length) ((c<<<i)%(N)) N (List.concat [zero;y;freexs;ctrls;[x.[i]]])

    CSWAP x y ctrls

    // uncompute y
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        ConstantSubtractorInplaceModN (1+ctrls.Length) ((cinvModN<<<i)%(N)) N (List.concat [zero;y;freexs;ctrls;[x.[i]]])

// Same as the above, use with Circuit.Compile (is just a wrapper)
let ConstantMultiplierModNCompile (ccount:int) (a:bigint) (ainv:bigint) (N:bigint) (qs:Qubits) =
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]

    ConstantMultiplierModN a ainv N x y zero ctrls

// This does the same as the functions above but RUNS the circuit
// Use this for large bit-sizes, as the circuits get too large otherwise.
// It generates the circuit for each modular addition, executes it, counts the Toffolis, and then does the next one.
let ConstantMultiplierModNRun (ccount:int) (a:bigint) (ainv:bigint) (N:bigint) (qs:Qubits) (initialState:int[]) = 
    let n = int ((qs.Length-1-ccount)/2)
    let ctrls = slice qs [(2*n+1)..(2*n+ccount)]
    let x = slice qs [0..(n-1)]
    let y = slice qs [n..(2*n-1)]
    let zero = [qs.[2*n]]
    let one = bigint 1

    let mutable currentState = initialState
    let mutable toffcount = 0

    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        let circuit = Circuit.Compile (ConstantAdderInplaceModN (1+ctrls.Length) ((a<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) |> MyCircuitExport
        currentState <- MyCircuitSimulateFast circuit currentState
        toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
        //show "Compute: %A / %A" i n

    let circuit = Circuit.Compile (CSWAPCompile (ctrls.Length)) (List.concat [x; y; ctrls]) |> MyCircuitExport
    currentState <- MyCircuitSimulateFast circuit currentState
    toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)

    // uncompute y
    for i in 0..(n-1) do
        let freexs = slice x (List.concat [[0..(i-1)];[(i+1)..(n-1)]])
        let circuit = Circuit.Compile (ConstantSubtractorInplaceModN (1+ctrls.Length) ((ainv<<<i)%(N)) N) (List.concat [zero;y;freexs;ctrls;[x.[i]]]) |> MyCircuitExport
        currentState <- MyCircuitSimulateFast circuit currentState
        toffcount <- toffcount + (circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length)
        //show "Uncompute: %A / %A" i n

    show "N=%A TC=%A" n toffcount
    currentState.[0..(2*n+ccount)]

// (Classical algorithm) determines a^-1 mod N 
// This is used for the modular multiplication to uncompute x after having computed y=ax mod N
let ModularInverse (a:bigint) (n:bigint) =
    let mutable t = bigint 0
    let mutable r = n
    let mutable newt = bigint 1
    let mutable newr = a

    while (newr = bigint 0) = false do
        let q = r / newr
        let oldt = newt
        newt <- (t - q*newt)
        t <- oldt
        let oldr = newr
        newr <- (r - q*newr)
        r <- oldr
    let t = (t+n)%n
    t

// Runs modular multiplication for 2^k-bit numbers up to 16kbit 
// Reports Toffoli count and result of the computation & uncomputation (y should be 0 afterwards)
[<LQD>]
let __RunModMultiplier (i:int)=
    let N16k = 1189731495357231765085759326628007130763444687096510237472674821233261358180483686904488595472612039915115437484839309258897667381308687426274524698341565006080871634366004897522143251619531446845952345709482135847036647464830984784714280967845614138476044338404886122905286855313236158695999885790106357018120815363320780964323712757164290613406875202417365323950267880089067517372270610835647545755780793431622213451903817859630690311343850657539360649645193283178291767658965405285113556134369793281725888015908414675289832538063419234888599898980623114025121674472051872439321323198402942705341366951274739014593816898288994445173400364617928377138074411345791848573595077170437644191743889644885377684738322240608239079061399475675334739784016491742621485229014847672335977897158397334226349734811441653077758250988926030894789604676153104257260141806823027588003441951455327701598071281589597169413965608439504983171255062282026626200048042149808200002060993433681237623857880627479727072877482838438705048034164633337013385405998040701908662387301605018188262573723766279240798931717708807901740265407930976419648877869604017517691938687988088008944251258826969688364194133945780157844364946052713655454906327187428531895100278695119323496808703630436193927592692344820812834297364478686862064169042458555136532055050508189891866846863799917647547291371573500701015197559097453040033031520683518216494195636696077748110598284901343611469214274121810495077979275556645164983850062051066517084647369464036640569339464837172183352956873912042640003611618789278195710052094562761306703551840330110645101995435167626688669627763820604342480357906415354212732946756073006907088870496125050068156659252761297664065498347492661798824062312210409274584565587264846417650160123175874034726261957289081466197651553830744424709698634753627770356227126145052549125229448040149114795681359875968512808575244271871455454084894986155020794806980939215658055319165641681105966454159951476908583129721503298816585142073061480888021769818338417129396878371459575846052583142928447249703698548125295775920936450022651427249949580708203966082847550921891152133321048011973883636577825533325988852156325439335021315312134081390451021255363707903495916963125924201167877190108935255914539488216897117943269373608639074472792751116715127106396425081353553137213552890539802602978645319795100976432939091924660228878912900654210118287298298707382159717184569540515403029173307292433820279730892035211089568317921472966412729818117888454622835234533481549572743196734571634129403963684506083367368909845783940774676061371524391659800174583099528253482739757522802678131797463892371605031905393218989649729585777477658297715183959931786937290802047985884915659668363557899564591609695996145630027327988546172359576120674838045740053640396189717133917288042702062394781105995775027475794051243596930032470884082214287907917283972886566217666099928929736659305895401188482643507350834137594417765254338122252480578715215887377442751706844188139483627322937448198918074963942373224651327246325315257389993715375779389735367127951944954922262846304494110166728898867798152486928838574921534377364723691418182036282103573740854602704670142191797776310269453779790001271690513762843505679799117285781202701919548686629332950416374405421584385147988271411879919984987621119528596275512991867785554210230568216704790576492500881265464921474895559459290185727309851736982251106571217022721811129816003904051940819258133219624777452655172398944653308802365249212729736452165680043648457346007468020497859675817059619455164828406637438597547809499476599841667405756636926594683316183551527318347414506716327282224104675912411296799047208452645924316981738234697893596112182781775518514959248984587435401557122123845271399302125402393337922749394455393324480804314744560744187982171228782269878349258607387176290025943618757718706454898247726013046950133059199119527299386145131246957734198972549900572364145426997855548761897904248475449664993993862428294060559878742860518378058716440614430083917881823389319994271944195984564613910857996369511365798097068122437154330466861277932794002897625255843879189558691957373003353663905696047447765667201581975169298500940792736093264986662256899062338573645314621899957140466596072241410164464110856818238166950911474536753400753329335262575148713324594079410166286488150567351475933291175562137155949959416056516881232905680654150013798790046699053524795019980299093858029924914578907051418730788232093062110578135615360678690695102511579902265905593484061796350040671437474498455768307744249390537909776724456222613234554195683081115022999893977896277499918028773373149365321042211436343920206922932005094285083739085355321384854749817569151104942327971701097722576153091798048279364363438424522407742110033167468367004558861384900599732150535646653656240721316539419706881676263475309286637885522950176418670198965931I
    let N8k = 1090748135619415929462984244733782862448264161996232692431832786189721331849119295216264234525201987223957291796157025273109870820177184063610979765077554799078906298842192989538609825228048205159696851613591638196771886542609324560121290553901886301017900252535799917200010079600026535836800905297805880952350501630195475653911005312364560014847426035293551245843928918752768696279344088055617515694349945406677825140814900616105920256438504578013326493565836047242407382442812245131517757519164899226365743722432277368075027627883045206501792761700945699168497257879683851737049996900961120515655050115561271491492515342105748966629547032786321505730828430221664970324396138635251626409516168005427623435996308921691446181187406395310665404885739434832877428167407495370993511868756359970390117021823616749458620969857006263612082706715408157066575137281027022310927564910276759160520878304632411049364568754920967322982459184763427383790272448438018526977764941072715611580434690827459339991961414242741410599117426060556483763756314527611362658628383368621157993638020878537675545336789915694234433955666315070087213535470255670312004130725495834508357439653828936077080978550578912967907352780054935621561090795845172954115955492451891216570613296084365855602561987825554164828983776100070639535058975031421325993557963092582246894596377877165664134176626579977819203711571822341535408979557604056707123186626577676416188241336594881408068968961998717359787691335649357729429863437279636352583621433720199268128212734394658767583155547333863583775463140939229166786852053411359347506292273580759659221329615535915374296921312179266278677936047679491864303619458615535439684608832386143593871684717434990189386528537845712625304664823679920329431996145195217962866020168617725228731562929836932081530970550958194455848937712504788617035029734883083588942942311417818515459029999559814744757117521882790786148595159691476349244130338981981174851473059507740147241521672189074910330952319180638318199504090691810250168776309204618406888213015285193048095145396995645340002985115979602262321310662883088594437748318499511605309954705034226005722298162852349163242408574091397080892968895451060657757201199454907101458809408532429510189097387093136564036047797365062857395373088494023754833694236076111826663219393138920039495403085305056656519919814219844108778484235361412730008772810730108656701948930183929294294950489931276241738768497131921385470975219731917721I
    let N4k = 1044388881413152506691752710716624382579964249047383780384233483283953907971557456848826811934997558340890106714439262837987573438185793607263236087851365277945956976543709998340361590134383718314428070011855946226376318839397712745672334684344586617496807908705803704071284048740118609114467977783598029006686938976881787785946905630190260940599579453432823469303026696443059025015972399867714215541693835559885291486318237914434496734087811872639496475100189041349008417061675093668333850551032972088269550769983616369411933015213796825837188091833656751221318492846368125550225998300412344784862595674492194616443716246933213029777899798818461970932131191348804759825573468593610446932516843547031871082533988701912296616359089685540536703762087577894961935770084452075260024141226911707933733334829279200075988660067382055419272677487938772582789085572793356060472402630868147530739655123992685056875065481628365716808346058511893983996604674482947443330804864016159772117403184411098879310084911985820576707955476414019490642324082779864275147623753677862525374215759379112292614276897779676360046934875731128797682090481735253923263326681166800529934324565723067016719014885465039728219454090779142376015916232389820650894246737I
    let N2k = 25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357I
    let N1k = 135066410865995223349603216278805969938881475605667027524485143851526510604859533833940287150571909441798207282164471551373680419703964191743046496589274256239341020864383202110372958725762358509643110564073501508187510676594629205563685529475213500852879416377328533906109750544334999811150056977236890927563I
    let N512 = 13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006083527I
    let Ns = [251I;65521I;4294967291I;18446744073709551557I;340282366920938463463374607431768211297I;115792089237316195423570985008687907853269984665640564039457584007913129639747I;N512;N1k;N2k;N4k;N8k;N16k]
    let ctrlcount = 1 // one control qubit (phase estimation)

    let N = Ns.[i]
    let x = (N >>> 5) + 5I
    let a = (N >>> 4) + 17I
    let n = int (System.Math.Ceiling(bigint.Log(N,2.)))
    let exact = ((a*x) % N)
    
    show "Running n = %A bits" n
    
    let k = Ket(2*n+1+ctrlcount)
    let qs = k.Qubits
    let initialState = Array.zeroCreate (qs.Length)
    let xb = BoolInt n x
    for i in 0..(xb.Length-1) do 
        initialState.[i] <- xb.[i]
        initialState.[i+n] <- 0
    initialState.[2*n] <- 0
    for i in 0..(ctrlcount-1) do
        initialState.[2*n+1+i] <- 1 // control qubit is 1 --> do it.

    let finalState = ConstantMultiplierModNRun ctrlcount a (ModularInverse a N) N qs initialState

    let res = [for i in 0..1..(n-1) do yield ((bigint finalState.[i]) <<< i)]
                                |> List.sum
    let y = [for i in 0..1..(n-1) do yield ((bigint finalState.[i+n]) <<< (i))]
                                |> List.sum
    show "Result: %A * %A mod %A = %A, y = %A" a x N res y
    show "Result - Exact: %A" (res-exact)
        
// Adds (classical) value a to quantum register containing x using the dirty qubit register g which 
// consists of two qubits and is left unchanged.
[<LQD>]
let __RunInplaceAdderModN (nmax:int) (N:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true
    let output = false//true

    let nmin = int (System.Math.Ceiling (System.Math.Log(float N, 2.)))
    if simulate = true then
        for n in nmin..nmax do
            for x in 0..(N-1) do
                for a in 0..(N-1) do
                    let k = Ket(n*2+numctrls)
                    let qs = k.Qubits
                    let circuit = Circuit.Compile (ConstantAdderInplaceModN numctrls (bigint a) (bigint N)) qs |> MyCircuitExport

                    let initialState = Array.zeroCreate (qs.Length)
                    let xb = BoolInt n (bigint x)
                    initialState.[0] <- 0
                    for i in 0..(xb.Length-1) do 
                        initialState.[1+i] <- xb.[i]
                        if (i < xb.Length-1) then
                            initialState.[1+i+n] <- (g>>>i)&&&1
                    for i in 0..(numctrls-1) do
                        initialState.[2*n+i] <- (ctrlvalue >>> i )&&&1
        
                    let finalState = MyCircuitSimulateFast circuit initialState
    
                    let res = [for i in 1..1..(n) do yield (finalState.[i] <<< (i-1))]
                                |> List.sum
                    let gout = [for i in (n+1)..1..(2*n-1) do yield (finalState.[i] <<< (i-1-n))]
                                |> List.sum
                    if (((x + int a)%(1<<<n)) = res) = false && (((x + (int a) - N)+(1<<<n))%(1<<<n) = res) = false  then
                        if not (ctrlvalue &&& ((1<<<numctrls)-1) = ((1<<<numctrls)-1)) then
                            assert(x=res)
                            if (x = res) = false then
                                show "Control fail!"
                        else
                            show "%A + %A != %A for n = %A" x a res n
                    elif output then
                        show "%A + %A = %A :-)" x a res
                    if ( (int (gout)) = g) = false then
                        show "garbage %A !!!! %A" g gout
                    if finalState.[0] = 1 then
                        show "z = %A for x = %A and c = %A, res = %A" finalState.[0] x a res
    else // determine toffoli counts
        let x = 0
        for ni in 0..100..nmax do
            let n = max ni 3
            let N = (1<<<n)-1          
            let a = N

            let k = Ket(n*2+numctrls)
            let qs = k.Qubits
            let circuit = Circuit.Compile (ConstantAdderInplaceModN numctrls (bigint a) (bigint N)) qs |> MyCircuitExport
            let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
            show "Toffoli-count@n=\t%A\t%A" n toffcount
    ()

// Adds (classical) value a to quantum register containing x using the garbage register g which 
// consists of two qubits and is left unchanged.
// NOTE: There is also a version which only requires 1 dirty qubit at a cost of an increased Toff count:
// 1 dirty qubit: 10 nlog(n), 2 dirty qubits: 8 nlog(n)
[<LQD>]
let __RunInplaceAdder (nmax:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true

    for n in 0..1..nmax do
        let a = bigint ((1L <<< n)-1L)
            
        let x = 212311%(1<<<n)

        let k = Ket(n+numctrls+2)
        let qs = k.Qubits
        let circuit = Circuit.Compile (ConstantAdderInplaceCompile numctrls a) qs
                    |> MyCircuitExport

        // tweak the circuit into a Toffoli network that can be simulated efficienly    
        if simulate = true then
            let initialState = Array.zeroCreate (n+2+numctrls)
            let xb = BoolInt n (bigint x)
            for i in 0..(xb.Length-1) do 
                initialState.[i] <- xb.[i]
            initialState.[n] <- (g>>>1)&&&1
            initialState.[n+1] <- g &&& 1
            for i in 0..(numctrls-1) do
                initialState.[n+2+i] <- (ctrlvalue >>> i )&&&1
        
            let finalState = MyCircuitSimulateFast circuit initialState
    
            let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                        |> List.sum
            let gout = (finalState.[n]+2*finalState.[n+1])
            if ( (int ((bigint x)+a))&&&((1<<<n)-1) = res) = false then
                show "%A + %A != %A" x a res
            if ( (int (gout)) = g) = false then
                show "Garbage changed!"
        let toffcount = circuit |> List.filter (fun y -> match y with | MyTOFF(a,b,c) -> true | _ -> false) |> List.length
        show "Toffoli-count@n=\t%A\t%A" n toffcount
    ()

// Test for inplace Adder
[<LQD>]
let __RunInplaceAdderTest (nmax:int) (numctrls:int) (ctrlvalue:int) =     
    let g = 3
    let simulate = true

    for n in 0..nmax do
        show "Testing for bit size n = %A" n
        for a in 0..((1<<<n)-1) do
            for x in 0..((1<<<n)-1) do
                let k = Ket(n+numctrls+2)
                let qs = k.Qubits
                let circuit = Circuit.Compile (ConstantAdderInplaceCompile numctrls (bigint a)) qs
                            |> MyCircuitExport

                // tweak the circuit into a Toffoli network that can be simulated efficienly    
                let initialState = Array.zeroCreate (n+2+numctrls)
                let xb = BoolInt n (bigint x)
                for i in 0..(xb.Length-1) do 
                    initialState.[i] <- xb.[i]
                initialState.[n] <- ((g>>>1)&&&1)
                initialState.[n+1] <- (g &&& 1)
                for i in 0..(numctrls-1) do
                    initialState.[n+2+i] <- ((ctrlvalue >>> i )&&&1)
        
                let finalState = MyCircuitSimulateFast circuit initialState
    
                let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                            |> List.sum
                let gout = (finalState.[n]+2*finalState.[n+1])
                let ctrlout = [for i in 0..1..(numctrls-1) do yield (finalState.[i+n+2] <<< (i))]
                            |> List.sum
                if ( (int ((bigint x)+(bigint a)))&&&((1<<<n)-1) = res) = false then
                    show "%A + %A != %A, ctrl = %A" x a res ctrlout
                if ( (int (gout)) = g) = false then
                    show "Garbage changed!"
    ()

[<LQD>]
let __RunBorrowedIncrementer (x:int) (g:int) = 
    show "Testing garbage incrementer."
    let bits x = (int (System.Math.Log((float x), 2.)) )+1
    let n : int = max (bits x) (bits g)
    let n : int = max n 2
    show "Doing computation with n = %A" n
    let k = Ket(2*n+1)
    let qs = k.Qubits
    
    let circuit = Circuit.Compile CtrlBorrowedQubitsIncrementer qs
                    |> MyCircuitExport
    // prepare x
    x |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs 0 1 |> ignore
    // prepare garbage register
    g |> fun x -> (bigint x) |> BoolInt n |> PrepBool qs n 1 |> ignore

    // tweak the circuit into a Toffoli network that can be simulated efficienly    
    let initialState = Array.zeroCreate (qs.Length)
    List.iter (fun i -> M !!(qs,i)) [0..qs.Length-1]
    for i in 0..(qs.Length-1) do 
        initialState.[i] <- qs.[i].Bit.v
    initialState.[2*n] <- 0
    
    let finalState = MyCircuitSimulateFast circuit initialState
    
    let res = [for i in 0..1..(n-1) do yield (finalState.[i] <<< (i))]
                |> List.sum
    let gout = [for i in n..1..(2*n-1) do yield (finalState.[i] <<< (i-n))]
                |> List.sum
    show "Result of calculation: %A + 1 = %A" x res
    show "Garbage in: %A, Garbage out: %A" g gout
    ()


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
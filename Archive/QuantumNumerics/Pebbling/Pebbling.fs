module Pebbling

open System
open System.Collections.Generic

let flipCount (width : int ref) = 
    fun n ->
        width := !width + 1

let flipCountWH (width : int ref) (height : int ref) = 
    fun n ->
        //printfn "%d" n
        if n + 1 > !height then
            height := n + 1
        width := !width + 1

let flipCon depth  =
    let tiles : bool array = Array.zeroCreate depth
    fun n ->
        tiles.[n] <- not tiles.[n]
        for i in 0..(tiles.Length-1) do
            if tiles.[i] then
                printf "%-3d" i
            else
                printf "   "
        printf "\n"

let doFlips flip f =
        for i in f do
            flip i

let _simplePebble _start finish _pebbles =
    let rec algo start pebbles =
        let mutable flips = []
        let pebbleTriangle frist last = 
            let mutable flips = []
            for i in frist..(last-1) do
                flips <- i :: flips
            List.rev flips
        match pebbles with
        | 0 -> []
        | pebbles when pebbles >= finish - start -> 
            pebbleTriangle start finish
        | pebbles when 2*pebbles - 1  > finish - start ->
            flips <-pebbleTriangle start (finish - (pebbles-1))
            flips <- flips @ (flips |> List.rev |> List.tail) 
            flips <- flips @ (algo (finish - (pebbles-1)) (pebbles-1))
            flips
        | pebbles ->
            flips <-pebbleTriangle start (start+pebbles)
            flips <- flips @ (flips |> List.rev |> List.tail) 
            flips <- flips @ (algo (start+pebbles) (pebbles-1))
            flips 
    let flips = algo _start _pebbles
    (flips @  (flips |> List.rev |> List.tail ))

let simplePebble start finish pebbles flip  = doFlips  flip (_simplePebble start finish pebbles)

let pebList levels = 
    let rec _pebList n la lb lc =
        if n > 0 then
            let lower = la @ lb @ List.map (fun (a,b) -> (a,b-1)) lc 
            let upper = List.map (fun (a,b) -> (a + pown 2 (levels-n+1) ,b)) (la @ lb @ lc) 
            _pebList (n-1) lower upper (List.rev lower)
        else 
            (la @ lb @ lc)
    _pebList levels [(0,0)] [(1,0)] [(0,0)]


let _recursivePebble levels pebbles =
    let pList = pebList levels
    let mutable flips = []
    for (a,b) in pList do
        flips <- flips @ _simplePebble (pebbles*a) (pebbles*(a+1)) (pebbles+b)
    flips

let flipBit  n (inp : uint64) : uint64 = (1UL <<< n) ^^^ inp
let checkBit n (inp : uint64) : bool = ((1UL <<< n) &&& inp) <> 0UL

//Count bits by clearing the LSB until number = 0
let countBits (inp : uint64) =
    let mutable x = inp
    let mutable res = 0
    while x <> 0UL do
        x <- x &&& (x-1UL)
        res <- res + 1
    res


let dynamicPebble (pebbles : int) (n : int ) flips =
    let evalMoves xs = 
        let mutable  pos : uint64 = 0UL
        for i in xs do
            pos <- flipBit i pos
        pos
    
    let posMap = new Dictionary<uint64, int>(HashIdentity.Structural)
    let bestMoves = ref ([],System.Int32.MaxValue)
    let rec algo (prevMoves : (int list)*int*uint64) =  
        let (moves,length,pos) = prevMoves 
        if length < snd !bestMoves then 
            let posInDict = posMap.TryGetValue(pos)
            if not (fst posInDict) || length < snd posInDict then  
                if checkBit (n-1) pos then
                    bestMoves := moves,length
                    //printf "Current best:%A\n" prevMoves
                else    
                    posMap.[pos] <- length
                    let pebsUsed = countBits pos 
                    for i in (n-2) .. -1 .. 0 do
                        if (checkBit i pos) && (pebsUsed < pebbles ||checkBit (i+1) pos) then
                            algo ((i+1) :: moves , length + 1 , flipBit (i+1) pos)
                    if (checkBit 0 pos) || (pebsUsed < pebbles) then
                        algo ( 0 :: moves , length + 1 , flipBit 0 pos)
        ()
    posMap.[0UL] <- 0
    algo ([0],1,1UL)
    let strategy = ((List.rev (fst !bestMoves)) @ (fst !bestMoves))
    doFlips flips strategy
    strategy

let dynamicPebble2 (pebbles : int) (depth : int )  =
    let funcMap = new Dictionary<int*int, int*int list>(HashIdentity.Structural)
    let rec F n s =   
        let funcInDict = funcMap.TryGetValue((n,s))
        if fst funcInDict then
            snd funcInDict
        else 
            match (n,s) with
            | (1,s) when s > 0 -> 
                funcMap.[(1,s)] <- (1,[0])
                (1,[0])
            | (n,1) when n > 1 ->
                funcMap.[(n,1)] <- (System.Int32.MaxValue/3,[])
                (System.Int32.MaxValue/3,[])
            | (n,s) -> 
                let mutable min = (0,System.Int32.MaxValue/3)
                for m in 1..(n-1) do
                    let (f0,_) =  F m s
                    let (f1,_) =  F m (s-1)
                    let (f2,_) = F (n-m) (s-1)
                    let sumf = f0 + f1 + f2
                    if sumf < snd min then
                        min <- (m,sumf)
                let min_m = fst min
                if snd min >= (System.Int32.MaxValue/3-1) then
                    funcMap.[(n,s)] <- (System.Int32.MaxValue/3,[])
                    (System.Int32.MaxValue/3,[])
                else 
                    let f0 =  F min_m s
                    let f1 =  F min_m (s-1)
                    let f2 = F (n-min_m) (s-1)
                    let sumf = fst f0 + fst f1 + fst f2
                    let value = (sumf, snd f0 @ List.map (fun x -> x + min_m) (snd f2) @ List.rev (snd f1) )
                    funcMap.[(n,s)] <- value
                    value
    for s in 2..pebbles do
        for n in 3..depth do
            ignore (F n s)
    let mutable ret = []
    for i in funcMap do 
        ret <- (i.Key,i.Value) ::ret
    ret <- List.filter (fun (_,(a,_)) -> a < System.Int32.MaxValue/3) ret
    ret <- List.sortBy (fun (_,(a,_)) -> a) ret
    List.map (fun (a,(_,c)) -> (a,c)) ret
            
                     
let recursivePebble levels pebbles flip = doFlips flip (_recursivePebble levels pebbles )


let printPeblist l =
    for (a,b) in l do
        for i in 0..a do
            printf "  "
        printf "%-+2d\n" b
        
let runPebble (depth:int) (pebbles:int) = 
    //makeDP 4 16
    //let depth = 300
    //let pebbles = 25
    
    //let width = ref 25
    let strats =  (dynamicPebble2 pebbles depth)
    //let chartedPebbles = [8;10;50;100;250]
    //let chartedPebbles = [9]
    let chartedPebbles = [pebbles]
    let mutable charts = []
    let mutable dwPlot = [1,2]

    for i in chartedPebbles do
        let mutable dwPlot = []
        let  currStrats = List.filter (fun ((_,s),_) -> s = i) strats
        for ((n,s),strat) in currStrats do
            printf "steps: %A pebbles: %A\n" n s
            printf "number of moves: %A\n" strat.Length
            printf "moves: %A\n" strat
            dwPlot <- (n,List.length strat) :: dwPlot   
        dwPlot <- List.rev dwPlot
   

[<EntryPoint>]
let main argv = 
    runPebble (int argv.[0]) (int argv.[1])
    0 // return an integer exit code

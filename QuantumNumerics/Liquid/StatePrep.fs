module Microsoft.Research.Liquid.StatePrep

open System
open System.Collections.Generic
open System.Numerics

open Util               // general utilities
open Operations         // gate definitions
open HamiltonianGates


let get_probabilities (ps:list<float>) (bit:int) =
    let num_probabilities = (1 <<< bit)
    let mask = num_probabilities - 1

    let get_probability p_index =
        let norm_list = List.map2 (fun i p -> match (i &&& mask) = p_index with
                                                | true -> p
                                                | false -> 0.) [0..ps.Length-1] ps
        let norm = List.reduce (+) norm_list
        let norm = match norm with
                    | 0. -> 1.
                    | _ -> norm
        let p_list = List.map2 (fun i p -> match (i >>> bit) &&& 1 = 0 with
                                            | true -> p / norm
                                            | false -> 0.) [0..ps.Length-1] norm_list
        List.reduce (+) p_list
    List.map get_probability [0..num_probabilities - 1]

let rec build_uniform_gate (angles:list<float>) (qs:Qubits) =
    if qs.Length <= 1 then
        if qs.Length = 1 then
            Rz (2. * angles.[0]) 1. "" [qs.[0]]
    else
        let even_angles = List.map (fun i -> angles.[2*i]) [0..angles.Length/2-1]
        let odd_angles = List.map (fun i -> angles.[2*i+1]) [0..angles.Length/2-1]
        let new_angles_plus = List.map2 (fun a1 a2 -> (a1 + a2) * 0.5) even_angles odd_angles
        build_uniform_gate new_angles_plus qs.[1..]
        CNOT [qs.[0]; qs.[qs.Length - 1]]
        let new_angles_minus = List.map2 (fun a1 a2 -> (a1 - a2) * 0.5) even_angles odd_angles
        build_uniform_gate new_angles_minus qs.[1..]
        CNOT [qs.[0]; qs.[qs.Length - 1]]

// prepares sum_i \sqrt(p_i) |i>
let StatePrep (p:list<float>) (qs:Qubits) =
    for i in 0..qs.Length-1 do
        let ps = get_probabilities p i
        let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
        Ybasis [qs.[i]]
        build_uniform_gate angles qs.[0..i]
        Ybasis' [qs.[i]]

// unprepares sum_i \sqrt(p_i) |i>
let StatePrep' (p:list<float>) (qs:Qubits) =
    for i in qs.Length-1..(-1)..0 do
        let ps = get_probabilities p i
        let angles = ps |> List.map (fun pr -> -Math.Acos(Math.Sqrt(pr)))
        Ybasis [qs.[i]]
        build_uniform_gate angles qs.[0..i]
        Ybasis' [qs.[i]]

[<LQD>]
let __TestLCUAmpExtraQubit () =
    let n = 3 // # bits
    let k = Ket(n)
    let qs = k.Qubits
    let p = [0.5; 0.5]

    let W (qs:Qubits) =
        Ry (2. * Math.Acos (Math.Sqrt 0.5)) 1. "" [qs.[1]]
        StatePrep p qs.[0..0]
        CNOT [qs.[0]; qs.[2]]
        X [qs.[0]]
        CZ [qs.[0]; qs.[2]]
        X [qs.[0]]
        StatePrep' p qs.[0..0]
    
    let W' (qs:Qubits) =
        Ry (-2. * Math.Acos (Math.Sqrt 0.5)) 1. "" [qs.[1]]
        StatePrep p qs.[0..0]
        CNOT [qs.[0]; qs.[2]]
        X [qs.[0]]
        CZ [qs.[0]; qs.[2]]
        X [qs.[0]]
        StatePrep' p qs.[0..0]

    let R (qs:Qubits) =
        X >< [qs.[0];qs.[1]]
        CZ [qs.[0];qs.[1]]
        X >< [qs.[0];qs.[1]]

    let ampl (qs:Qubits) =
        W qs
        for i in 0..0 do
            Gtheta (-Math.PI) 1. qs
            R qs
            W' qs
            R qs
            W qs
    let amplify  = Circuit.Compile ampl qs
    let mop (qs:Qubits) =
        M [qs.[0]]
        M [qs.[1]]

    let meas = Circuit.Compile mop qs

    amplify.Run qs
    show "Probability of success (00): %A" (k.Probs [qs.[0];qs.[1]]).[0]
    
    meas.Run qs
    let state = k.State qs.[qs.Length-1]
    let vec = (k.Probs [qs.[qs.Length-1]])
    show "State: %A Probabilities: %A" state vec
    let outcome = match [qs.[0].Bit;qs.[1].Bit] with
                    | [Bit.Zero;Bit.Zero] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!!"
    show "%A" outcome

[<LQD>]
let __TestLCUAmpPlusMinusID () =
    let n = 3 // # bits
    let k = Ket(n)
    let qs = k.Qubits
    let p = [0.5/Math.Sqrt 2.; 0.5/Math.Sqrt 2.; 0.5-0.5/Math.Sqrt 2.; 0.5-0.5/Math.Sqrt 2.]

    let W (qs:Qubits) =
        StatePrep p qs.[0..1]
        
        X [qs.[1]]
        CCNOT [qs.[0]; qs.[1]; qs.[2]]
        X [qs.[0]]
        Cgate CZ [qs.[0]; qs.[1]; qs.[2]]
        X [qs.[0]]
        X [qs.[1]]
        CZ [qs.[0]; qs.[1]] // distribute the remaining amplitudes equally among +/- ID contributions

        StatePrep' p qs.[0..1]

    let R (qs:Qubits) =
        X >< [qs.[0];qs.[1]]
        CZ [qs.[0];qs.[1]]
        X >< [qs.[0];qs.[1]]

    let ampl (qs:Qubits) =
        W qs
        for i in 0..0 do
            Gtheta (-Math.PI) 1. qs
            R qs
            W qs
            R qs
            W qs
        
    let amplify  = Circuit.Compile ampl qs
    let mop (qs:Qubits) =
        M [qs.[0]]
        M [qs.[1]]
    let meas = Circuit.Compile mop qs

    amplify.Run qs
    show "Probability of success (00): %A" (k.Probs [qs.[0];qs.[1]]).[0]
    
    meas.Run qs
    let state = k.State qs.[qs.Length-1]
    let vec = (k.Probs [qs.[qs.Length-1]])
    show "State: %A Probabilities: %A" state vec
    let outcome = match [qs.[0].Bit;qs.[1].Bit] with
                    | [Bit.Zero;Bit.Zero] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!!"
    show "%A" outcome

[<LQD>]
let __TestLCU () =
    let n = 6 // # bits
    let p = [0.25;0.1;0.1;0.05;0.1;0.025;0.125;0.25]
    let p = [0.24;0.20;0.3;0.2;0.05;0.;0.01;0.]
    let k = Ket(n)
    let qs = k.Qubits
    
    let W (qs:Qubits) =
        for i in 0..2 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i+3]]
            build_uniform_gate angles qs.[3..i+3]
            Ybasis' [qs.[i+3]]
        for i in 0..2 do
            CNOT [qs.[3+i]; qs.[i]]
        for i in 2..(-1)..0 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> -Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i+3]]
            build_uniform_gate angles qs.[3..i+3]
            Ybasis' [qs.[i+3]]

    let ops (qs:Qubits) =
        W qs
        

    let circ  = Circuit.Compile ops qs
    let mop (qs:Qubits) =
        M >< qs.[3..]
    let meas = Circuit.Compile mop qs
    circ.Run qs
    let norm2 = List.fold (fun s p -> s + p * p) 0. p
        
    show "Probability of success (000): %A == %A" (k.Probs qs.[3..]).[0] norm2
    meas.Run qs
    let vec = (k.Probs qs.[0..2])
    let one_norm = Array.fold (fun s p -> s + Math.Sqrt p) 0. vec
    show "%A" (List.map (fun i -> (Math.Sqrt vec.[i])/one_norm) [0..7])
    let measurement = List.map (fun i -> (k.Bit qs.[i+3]).v) [0..2]
    let status = match measurement with
                    | [0;0;0] -> "SUCCESS!!!"
                    | _ -> "FAILURE!!! Please, try again."
    show "%A : %A" measurement status
    show "Done"

[<LQD>]
let __TestStatePrep () =
    let n = 4 // # bits
    //let p = [0.25;0.1;0.1;0.05;0.1;0.025;0.125;0.25]
    let p = [0.016393443;0.06557377;0.049180328;0.147540984;0.016393443;0.049180328;0.032786885;0.016393443;0.081967213;0.131147541;0.06557377;0.049180328;0.147540984;0.032786885;0.081967213;0.016393443]
    let k = Ket(n)
    let qs = k.Qubits
    
    let ops (qs:Qubits) =
        for i in 0..n-1 do
            let ps = get_probabilities p i
            let angles = ps |> List.map (fun pr -> Math.Acos(Math.Sqrt(pr)))
            Ybasis [qs.[i]]
            build_uniform_gate angles qs.[0..i]
            Ybasis' [qs.[i]]

    let circ  = Circuit.Compile ops qs
    circ.Run qs
    let vec = (k.Probs qs)
    for i in 0..vec.Length-1 do
        if Math.Abs (vec.[i] - p.[i]) > 1.e-9 then
            show "Failed!!!"
    show "%A" (vec)
    show "Done"
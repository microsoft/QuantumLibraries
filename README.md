Library and samples in LIQUi|>, Solid, and Qb
=============================================

#Qb libs and samples
This repo contains Qb libraries and samples. The following are targets for December 2017: 

- Libraries
  - Phase Estimation library
  - Amplitude Estimation library
  - Basic Arithmetic library  

- Scenarios/Samples:
  -Shor’s algorithm
    - rotation based
    - calls Phase Estimation and Basic Arithmetic library

  - Hamiltonian simulation
    - Trotter based.  
              
  - “Nielsen and Chuang” tutorial samples
    - Teleport
    - Superdense coding
    - QFT (small code snippet)
    - Grover’s algorithm (for some hard coded oracle function f)
    - Deutsch-Jozsa algorithm (for some hard coded oracle function f)
    - Hidden shift algorithm (for some hard coded oracle function f)
    - Basic phase estimation (to lead up to Shor and Ham. Sim.)
    - Basic error correction
    - RUS circuit example
    - Circuit identities
      - Toffoli via Cliff + T
      - Phase Toffoli via Cliff + T
      - Multi-target CNOTs
      - Measuring an operator 
              
#LIQUi|> libs  
This repo also contains LIQUi|> quantum circuit generator libraries and samples for several intended target uses that 
extend the libraries already built-in to LIQi|>. These libraries are the following:
- Basic circuits (LIQUi|>): 
  - basic tools for Toffoli gates (simulator, depth calculation, metrics)
  - multiply controlled gates
  - tools for LFSR based counters
- Integers (LIQUi|>): 
  - integer addition, subtraction, multiplication, division
- Quantum Montgomery library (LIQUi|>): 
  - modular addition, subtraction, multiplication, division
  - elliptic curve point addition
- Numerics library (LIQUi|>): 
  - Implementation of various useful functions such as 1/x, 1/sqrt(x), arcsin
- Polynomials (LIQUi|>): 
  - basic fixed point arithmetic, polynomial evaluation
- Floats (LIQUi|>):  
  - floating point arithmetic, basic operations of addition and multiplication
- Signal processing (LIQUi|>): 
  - Quantum Fourier transforms
  - Other signal transforms (DCTs, DSTs, wavelets)
- Data types (LIQUi|>): 
  - circuits to implement more advanced quantum data types such as sets
- Patterns: (to be implemented)
  - basic algorithmic patterns, including phase estimation (PE), amplitude amplification (AA), amplitude estimation (AE),
    oblivious amplitude estimation (OAA), quantum rejection sampling (QRS), linear combination of unitaries (LCU)

#Build and Test
To build the VS project libs, go into the subfolder for each lib, open the .sln file and build the solution. Ideally, you want to build this with Target F# runtim (at least) F# 4.0 (FSharp.Core, 4.4.0.0) and Target Framework 4.6. 

The solutions come with several tests which should be easily identifiable as they start with two underscores, for instance function __RunSmallAdderTests() in QuantumMontgomeryTests.fs which can be run in the command line as

C:\QuantumMontgomery\Liquid\bin\Release>liquid __RunSmallAdderTests()


Library development for LIQUi|>, Solid, and Qb
==============================================

#Introduction 
This repo contains quantum circuit generator libraries for several intended target uses. Currently supported are: 
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
- Patterns: 
  - basic algorithmic patterns, including phase estimation (PE), amplitude amplification (AA), amplitude estimation (AE)
    oblivious amplitude estimation (OAA), quantum rejection sampling (QRS), linear combination of unitaries (LCU)

#Build and Test
To build the VS project libs, go into the subfolder for each lib, open the .sln file and build the solution. Ideally, you want to build this with Target F# runtim (at least) F# 4.0 (FSharp.Core, 4.4.0.0) and Target Framework 4.6. 

The solutions come with several tests which should be easily identifiable as they start with two underscores, for instance function __RunSmallAdderTests() in QuantumMontgomeryTests.fs which can be run in the command line as

C:\QuantumMontgomery\Liquid\bin\Release>liquid __RunSmallAdderTests()


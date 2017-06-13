Library development for LIQUi|> and Solid
=========================================

#Introduction 
This repo contains quantum circuit generator libraries for several intended target uses. Currently supported are: 
- Quantum Montgomery library (LIQUi|>): 
  - integer addition, subtraction, multiplication, division
  - modular addition, subtraction, multiplication, division
  - elliptic curve point addition
- Quantum numerics library (LIQUi|>): 
  - fixed point arithmetic, implementation of various useful functions such as 1/x, 1/sqrt(x), arcsin, polynomials. 
  - floating point arithmetic: basic operations of addition and multiplication

#Build and Test
To build the VS project libs, go into the subfolder for each lib, open the .sln file and build the solution. Ideally, you want to build this with Target F# runtim (at least) F# 4.0 (FSharp.Core, 4.4.0.0) and Target Framework 4.6. 

The solutions come with several tests which should be easily identifiable as they start with two underscores, for instance function __RunSmallAdderTests() in QuantumMontgomeryTests.fs which can be run in the command line as

C:\QuantumMontgomery\Liquid\bin\Release>liquid __RunSmallAdderTests()


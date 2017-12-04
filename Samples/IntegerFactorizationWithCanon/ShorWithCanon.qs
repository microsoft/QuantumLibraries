// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace IntegerFactorizationWithCanon
{
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Extensions.Math;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Canon;

    /// # Summary 
    /// Uses Shor's algorithm to factor a `number`
    ///
    /// # Input 
    /// ## number
    /// An integer to be factored
    /// ## useRobustPhaseEstimation
    /// If set to true, we use Microsoft.Quantum.Canon.RobustPhaseEstimation and 
    /// Microsoft.Quantum.Canon.QuantumPhaseEstimation
    ///
    /// # Output 
    /// Pair of numbers p > 1 and q > 1 such that p⋅q = `number`
    operation ShorWithCanon ( number : Int, useRobustPhaseEstimation : Bool ) : (Int,Int) {
        body
        {
            //  first check the most trivial case 
            if (  number % 2 == 0 ) {
                Message("An even number has been passed; 2 is the factor.");
                return (number / 2, 2);
            }
            // Trying to guess a co-prime to number
            // First get a random integer in the interval [1,number-1]
            let coprimeCandidate = RandomInt(number - 2) + 1;
            if( IsCoprime(coprimeCandidate, number) ) {

                Message($"Estimating period of {coprimeCandidate}");

                let period = EstimatePeriod(coprimeCandidate, number, useRobustPhaseEstimation);
                if( period % 2 == 0 ) {
                    let halfPower = ExpMod(coprimeCandidate,period/2,number);
                    if( halfPower != number - 1 ) {
                        let factor = MaxI(GCD(halfPower-1, number), GCD(halfPower+1, number));
                        return (factor,number/factor);
                    }
                    else {
                        fail "Residue xᵃ = -1 (mod N) where a is a period.";
                    }
                }
                else {
                    fail "Period is odd.";
                }
                return (0,0);                
            }
            else { // In this case we guessed a divisor by accident
                let gcd = GCD(number,coprimeCandidate);
                Message($"We have guessed a divisor of {number} to be {gcd} by accident.");
                return ( gcd, number / gcd );
            }
        }
    }

    /// # Summary 
    /// Interprets `target` as encoding unsigned little-endian integer k 
    /// and performs transformation |k⟩ ↦ |gᵖ⋅k mod N ⟩ where 
    /// p is `power`, g is `generator` and N is `modulus`.
    /// 
    /// # Input 
    /// ## generator 
    /// The unsigned integer multiplicative order ( period )
    /// of which is being estimated. Must be co-prime to `modulus`.
    /// ## modulus
    /// The modulus which defines the residue ring Z mod `modulus` 
    /// in which the multiplicative order of `generator` is being estimated.
    /// ## power 
    /// Power of `generator` by which `target` is multiplied.
    /// ## target 
    /// Register interpreted as LittleEndian which is multiplied by 
    /// given power of the generator. The multiplication is performed modulo 
    /// `modulus`.
    operation OrderFindingOracle(
        generator : Int, modulus : Int, power : Int , target : Qubit[] ) : () {
        body {
            AssertBoolEqual(
                IsCoprime(generator,modulus), true,
                "`generator` and `modulus` must be co-prime" );

            ModularMultiplyByConstantLE(
                ExpMod(generator,power,modulus),
                modulus,
                LittleEndian(target)
                );
        }
        adjoint auto
        controlled auto
        adjoint controlled auto
    }

    /// # Summary 
    /// Finds a multiplicative order of the generator 
    /// in the residue ring Z mod `modulus`.
    /// 
    /// # Input 
    /// ## generator 
    /// The unsigned integer multiplicative order ( period )
    /// of which is being estimated. Must be co-prime to `modulus`.
    /// ## modulus
    /// The modulus which defines the residue ring Z mod `modulus` 
    /// in which the multiplicative order of `generator` is being estimated.
    /// ## useRobustPhaseEstimation
    /// If set to true, we use Microsoft.Quantum.Canon.RobustPhaseEstimation and 
    /// Microsoft.Quantum.Canon.QuantumPhaseEstimation
    /// 
    /// # Output 
    /// The period ( multiplicative order ) of the generator mod `modulus`
    operation EstimatePeriod( generator : Int, modulus : Int, useRobustPhaseEstimation : Bool ) : Int {
        body{
            AssertBoolEqual(
                IsCoprime(generator,modulus), true,
                "`generator` and `modulus` must be co-prime" );

            mutable result = 1;
            let bitsize = BitSize( modulus );
            let bitsPrecision = 2*bitsize + 1;

            repeat {    
                mutable dyadicFractionNum = 0;
                using( eignestateRegister = Qubit[bitsize]  ) {
                    let eignestateRegisterLE = LittleEndian(eignestateRegister);

                    // initialize eignestateRegister to 1 which is a superposition of 
                    // the eigenstates we are estimating the phases of
                    InPlaceXorLE(1,eignestateRegisterLE);

                    // Find the numerator of a dyadic fraction that approximates 
                    // s/r where r is the multiplicative order ( period ) of g
                    if( useRobustPhaseEstimation )
                    {
                        let phase = RobustPhaseEstimation(
                            bitsPrecision, 
                            DiscreteOracle(OrderFindingOracle(generator,modulus,_,_)),
                            eignestateRegisterLE
                            );
                        set dyadicFractionNum = 
                            Round( phase * ToDouble(2 ^ bitsPrecision  ) / 2.0 / PI() ) ;
                    }
                    else {
                        // when using QuantumPhaseEstimation we will need extra `bitsPrecision`

                        using ( dyadicFractionNumerator = Qubit[bitsPrecision] ) {
                            let dyadicFractionNumeratorBE = BigEndian(dyadicFractionNumerator);

                            QuantumPhaseEstimation(
                                DiscreteOracle( OrderFindingOracle(generator,modulus,_,_) ),
                                eignestateRegisterLE,
                                dyadicFractionNumeratorBE);

                            set dyadicFractionNum = MeasureIntegerBE(dyadicFractionNumeratorBE);
                        }
                    }
                    ResetAll(eignestateRegister);
                }

                if ( dyadicFractionNum == 0 ) {
                    fail "We measured 0 for the numerator";
                }

                Message($"Estimated eigenvalue is {dyadicFractionNum}/2^{bitsPrecision}.");

                let (numerator,period) = 
                    ContinuedFractionConvergent(
                        Fraction(dyadicFractionNum, 2^(bitsPrecision)), 
                        modulus);
                let (numeratorAbs,periodAbs) = (AbsI(numerator), AbsI(period));

                Message($"Estimated period is {periodAbs}, we have projected on eigenstate " +
                        $"marked by {numeratorAbs}.");

                set result = periodAbs * result / GCD(result,periodAbs);
            }
            until( ExpMod(generator,result,modulus) == 1 ) 
            fixup {
                Message($"It looks like the period has divisors and we have " + 
                        $"found only a divisor of the period. Trying again ...");
            }

            return result;
        }
    }
}

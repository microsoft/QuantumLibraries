# Q# API Design Meeting / 24 July 2020
Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @msoeken.

## Agenda

- New `Microsoft.Quantum.Random` namespace

## Discussion

Proposal: https://github.com/microsoft/QuantumLibraries/issues/304

- This proposal requires changes to the `QSharpCore` project in the qsharp-runtime repo, and thus should be done as soon as possible to help manage dependency.
    - How difficult would it be to adapt simulation runtime and translators? Should be fine, as runtime changes localized to `SimulatorBase`.
    - **Resolution:** Create issue on qsharp-runtime to track needed modifications to QSharpCore and `SimulatorBase`.

- Impact on / from Q# language discussions:
    - Possible language discussion area: type classes.
        - E.g.: Cannot meaningfully compare equality on callables, so even equality requires type-class–like concept.
        - Q# API design should consider what reasonable type classes may look like.

- API for target machines to expose random seed setting.
    - For Q# standalone implementation, should we propagate the random seed settings via project properties?
    - Inappropriate to expose at Q# level, since not all targets need expose a PRNG at all; could be QRNG-backed.
    - The proposed API intentionally does not distinguish between an engine that produces randomness, and distributions based on that source of randomness.
    - **Resolution:** Proposal unmodified, as random seeds are simulator-specific and shouldn't be exposed at Q# level.
    
- No distribution over complex numbers in current proposal.
    - Should we introduce `ComplexDistribution`?
    - Ensure coverage for all reasonable arithmetic data types (e.g.: `BigDiscreteDistribution`, eventually `Distribution<BigInt>`).
      Current arithmetic data types in Q#:
        - `Fraction`, `BigFraction`
        - `Complex`, `ComplexPolar`
        - `Int`, `BigInt`
        - `Double`
    - Don't try to be minimal, but aim for coverage. Minimal only in sense that isn't redundant.
    - **Resolution:** Proposal modified to improve coverage of arithmetic data types.

- Introducing new arithmetic UDTs comes at a large cost: need to support in all future arithmetic APIs. Opportunity to improve existing support to better cover `BigInt` (https://github.com/microsoft/QuantumLibraries/issues/95).
    - **Proposal**: Adopt as API design principle that functionality should cover all relevant built-in and standard library types.

**Consensus**:
- The proposed API changes are reasonable and meet with Q# design principles.

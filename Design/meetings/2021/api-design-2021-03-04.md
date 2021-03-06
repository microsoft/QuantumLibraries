# Q# API Design Meeting / 4 March 2021

Attendees (in order by username): @bettinaheim, @cgranade, @guenp, @msoeken.

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/408
- https://github.com/microsoft/QuantumLibraries/issues/413
- https://github.com/microsoft/QuantumLibraries/issues/424
- https://github.com/microsoft/QuantumLibraries/issues/418
- https://github.com/microsoft/QuantumLibraries/issues/405
- Preliminary discussion of https://github.com/microsoft/QuantumLibraries/issues/419

## Discussion

### Library support for multidimensional arrays (microsoft/qsharp-language#49)

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/408

**Discussion**:

- Previous consensus: Defer until next review, read in detail and discuss async until then.
- Good to take inspiration from other high-level multidimensional array APIs like NumPy.
  - Some differences reflecting different terminology in Q# from Python (e.g.: "size" vs "shape").
  - **Action item**: Rename "reshape" → "resize" for consistency with QEP 3.
- Not specific to this proposal, but lack of overloading can make names confusing.
- Selection of different functions and operations makes sense.
- Nice to provide guidance on how users can write their own.
- Perhaps reasonable to schedule out-of-band review to do deeper dive.
- Good to have consensus by approximately July.
- Add application-driven examples (e.g.: BSM).

**Consensus**: Defer until out-of-band review, revise in the meantime.

<!--  -->

### Additions to Microsoft.Quantum.Math for floating point support

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/413

**Discussion**:

- Tense of `FiniteFact` may not match `EqualityFact`; perhaps `FinitenessFact`?
  - "Finite" is a synonym of "finiteness," such that a slight inconsistency in tense may be reasonable here for the sake of brevity.

**Consensus**: Approve as-is.

<!--  -->

### Autoemulation attribute

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/424

**Discussion**:

- "Emulate with" is a verb phrase.
- "Emulatable with" may make clear that the attribute isn't doing the substitution, but is providing metadata to the relevant rewrite step.
- "Emulation" may not be the right root either; one way to use it, but functionality is more general: target-specific substitution. Simulators are one specific kind of target.
  - Alternatively, can narrow definition of attribute; e.g.: "simulator specific implementation."
  - Only supported on C# runtime for now, but keep more general in case other runtimes support it eventually?
- Microsoft.Quantum.Targeting as namespace instead of Microsoft.Quantum.Core?
- **Action item**: Modify to `Microsoft.Quantum.Targeting.SubstitutableOnTarget("Alternative", "TargetName")`.

**Consensus**: Approved with modification above.

<!--  -->

### Deprecate Microsoft.Quantum.Environment operations

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/418

**Discussion**:

- Long-term goal would be to do qubit allocation statically, would be very difficult to do dynamically with this API.
- May not be usecases for this feature; need other ways of exposing width–depth tradeoffs to static analysis.
- Dynamic qubit management would require JIT anyway.

**Consensus**: Approved as-is.

<!--  -->

### DoubleAsInt

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/405

**Discussion**:

- Previous consensus: Introduce `DoubleAsInt`, effectively as an alias for `Truncate`.
- Revisiting due to additional feedback from @bettinaheim.
- **Action item**: Duplicate table of different truncate/round/ceil/floor outputs in each different function.

**Consensus**: Resolve with additional API documentation to clarify.

<!--  -->


### Formatting functions

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/419

**Administrative note**: Proposal still in draft, but request to consider deprecating previous functions immediately.

**Consensus**: Deprecate formatting functions immediately.
 
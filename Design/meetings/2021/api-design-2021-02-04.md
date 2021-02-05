# Q# API Design Meeting / 4 February 2021

Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @guenp, @msoeken.

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/405
- https://github.com/microsoft/QuantumLibraries/issues/406
- https://github.com/microsoft/QuantumLibraries/issues/407
- https://github.com/microsoft/QuantumLibraries/issues/408
- https://github.com/microsoft/QuantumLibraries/issues/409

## Discussion

### DoubleAsInt

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/405

**Discussion**:

- In C#, does `(int)x` for `double x` truncate `x`?
- Does the omission of `DoubleAsInt` force people to be more intentional about choosing how they want to convert?
- From language side, would be happy to do lossless conversions automatically. Functions are already used to mark conversions that may fail or that may lose information.
- Can consider `MaybeDoubleAsInt` as alternative, but that introduces other problems such as comparing equality on doubles.
- Can consider code action and API documentation to suggest that `Truncate` is more explicit than `DoubleAsInt`.

**Consensus**: Introduce `DoubleAsInt`, effectively as an alias for `Truncate`.

<!--  -->

### CurriedOp

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/406

**Consensus**: Approve as proposed.

<!--  -->

### QIR Attributes

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/407

**Discussion**:

- Proposed attributes won't be used directly by most Q# users. `@TargetInstruction` may be used directly by targeting package authors, while `@Inline` would be an advanced feature.
- Should we consider the opposite, introducing `@DontInline` to suppress inlining instead?
- Plan is for compiler to have good heuristics for determining when to inline and not.
- Part of more general plan to offer "hints" to the compiler in the form of attributes.
- Some users like to express finer-grained control over compilation process; heuristic attributes can be a good way of exposing that.
- Can only ever be hints, since inlining will always be disabled for function values (e.g.: passing a function to another function or partially applying a function).
- Only used in QIR emission, as documented in API documentation comment.

**Consensus**: Approve as proposed.

<!--  -->

### Library support for multidimensional arrays (microsoft/qsharp-language#49)

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/408

**Discussion**:

- May be good to defer detailed discussion until language feature approved.
- Discussing now can help iterate on proposal, and to help identify if any gaps are missing that would require modifying corresponding language proposal.
- Good to take inspiration from other high-level multidimensional array APIs like NumPy. E.g.:
  - Add `Zeros` and `Ones` functions
  - Add random sampling operations

**Consensus**: Defer until next review, read in detail and discuss async until then

<!--  -->

### Support for single-qubit Clifford group

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/409

**Discussion**:

- `ApplyPauli` → `ApplyP` could be confusing to users. Consider keeping `ApplyPauli` as an alias.
- `ApplyP` may be more consistent if thinking of "Pauli" as a type, but "Pauli" is likely to be read as a domain-specific term instead of as the type. As a result, `ApplyPauli` may not be inconsistent.
- `ApplyP` may not be accessible to new developers who aren't used to thinking of `Pauli` of a type in its own right.
- Aliases should be discouraged in general; creates opportunity for ambiguity and confusion. This is an exceptional case, however, and may contravene general push against aliases and/or type suffix notation.
- **Action item**: Modify proposal to introduce `ApplyP` without deprecating `ApplyPauli`.

**Consensus**: Approved with modification noted above.

<!--  -->

# Q# API Design Meeting / 11 May 2021

Attendees (in order by username): @bettinaheim, @cgranade, @msoeken.

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/448
- https://github.com/microsoft/QuantumLibraries/issues/377
- https://github.com/microsoft/QuantumLibraries/issues/453
- https://github.com/microsoft/QuantumLibraries/issues/405

## Discussion

### New special functions in Microsoft.Quantum.Math

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/448

**Discussion**:

- Rename `FactorialD` to `ApproximateFactorial`?
- Change input type for `ApproximateFactorial` to `Int`.
- Type suffixes are good to keep for future expansion and for consistency.

**Consensus**: Approve with above change to `FactorialD`.

<!--  -->

### Consistency between ApplyIf operations for Result and Bool variants #377

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/377

**Discussion**:

- Cannot make this change without breaking code.
- Proposed ordering for `ApplyIf` consistent with `ApplyIfElse`, more closely follows Q# style guide.
- Worth the breaking change to offer a more consistent experience to developers.

**Consensus**: Approve as written, but defer implementation until next breaking changes needed anyway.

<!--  -->

### Attributes for denoting API stability levels

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/453

**Discussion**:

- `@Stable()` implies two purposes: to mark stability, and to mark when that item became stable.
- What should we do with existing (unlabeled) functions/operations/types?
  - Leave unlabeled, but start auditing in run-up to Q# 1.0 and marking things we're good with being in 1.0 as `@Stable()`. Goal to have entire stdlib covered with `@Stable()` by 1.0.
- Modify `@Unstable()` to require link to GitHub issue for tracking discussions.
  - Good to enforce at compile time that link goes to tracking discussion — helps avoid things staying as unstable without either promoting to stability or removing.
  - Could even imagine IRewriteStep that checks that links in `@Unstable()` resolve to open GitHub issues.
- Should we add attributes to API before compiler support, or gate on https://github.com/microsoft/qsharp-compiler/issues/998 having completed?
  - Good to at least be able to document `@Unstable` and `@Stable`.
  - https://github.com/microsoft/qsharp-compiler/issues/998 could be good for external contribution ("help wanted" issue tag).
- **Action item**: Double check that "unstable" is broadly used in this connotation.
- `Microsoft.Quantum.Documentation` is a reasonable namespace.

**Consensus**: Approved with modification to `@Unstable()`, but defer implementation until at least documentation generation support available on the compiler.

<!--  -->

### DoubleAsInt should be added to balance out IntAsDouble

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/405

**Discussion**:

- **Previous consensus:** Approve as is.
- Revisiting due to comment from @bettinaheim.
- Having both `XAsY` and `YAsX` directions more important than slight redundancy.

**Consensus**: Previous consensus holds, approve as written.

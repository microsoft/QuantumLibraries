# Q# API Design Meeting / 30 August 2021

Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @guenp, @msoeken.

## Agenda

- Deep dive discussion of https://github.com/microsoft/QuantumLibraries/issues/408

## Discussion

### Library support for multidimensional arrays

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/408

**Discussion**:

- Can 2, 3, etc. suffixes be removed and/or consolidated?
  - Expecting to be able to use one function `Concatenated` for each different dimensionality.
  - Consider exposing `*Unsafe` functions to allow consolidating functions of different dimensionality?
    - Even if unsafe functions were exposed, users would still have to manually specify type parameters in many cases. As a result, exposing unsafe functions may not necessarily be easier to use.
  - Future Q# language proposals could help (e.g.: bounded polymorphism).
- Should single-dimensional case have `1` as suffix?
  - Expect `1` to be specific to one-dimensional case, no suffix to be polymorphic over types differing in dimension.
  - Changing proposal to add `1` suffix would allow introducing polymorphic variants later when language features support.
  - Unsuffixed versions used commonly in Microsoft.Quantum.Arrays (e.g.: `ElementAt` already exists without `1` suffix).
- Rank 8 limitation?
  - Arbitrary rank limitations for multidimensional arrays are common in other languages.
  - Defer considering extending on user feedback indicating need for dimension > 8.
- Relationship to default values, e.g.: in `DiagonalArray𝑁`?
  - Default values in general aren't well-defined for all Q# types; need a feature like bounded polymorphism to implement.
  - Could also consider Rust-style, e.g.: `fn eye<'T: Zero + One>(n : usize) -> Array2<T> { ... }`.
- Rename `Length` → `Size`?
  - Not ideal to have `Length` as special case of `Size𝑁`, rather than `Size` to match patterns elsewhere in the proposal.
  - Changing `Length` → `Size` is possible to do without breaking changes, but would be a much larger change.
  - Defer any decision on `Length` vs `Size` to future proposals.
- Consider additional feedback from expected users before implementing final design.
  - E.g.: "How do you think about matrices? How do you write them on paper? How do you write them in {Python, Julia, ...}?"

**Consensus**: Approved as-is.

# Q# API Design Meeting / 10 September 2020

Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @guenp, @msoeken.

## Agenda

- Improvements to `Microsoft.Quantum.Array` namespace
- Removing `Microsoft.Quantum.Simulation.QuantumProcessor.Extensions` namespace

## Discussion

### Microsoft.Quantum.Arrays

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/313

**Discussion**:

- Zip is very commonly used in existing  Q# code, is it worth being consistent in this case at the cost of causing significant breakage?
  - Consistency is important for reaching Q# 1.0.
  - Can we resolve using code actions as aliases instead of deprecation stubs?
  - No conventions and nomenclature are truly universal across languages.
  - Helped by IntelliSense features as well, need matching investment in IQ#.
  - Modify unknown identifier error message to provide actionable guidance ("did you mean...?")
  - Aliases can be helpful, but are maintenance nightmare.
  - **Action item**: discuss aliases, error messages, and/or other possible onboarding paths further.
- Likely missing: function `((Int -> 'T), Range) -> 'T[]`.
  - There as `MappedOverRange`, but should we have a better name?
  - Perhaps "over" is confusing?
  - Inconsistent w/ use of suffixes for input types.
  - **Action item**: Propose better name here, addressing "over" and type suffix conventions.
- `CumulativeFolded`: Should it return initial state or not?
  - `np.cumsum` excludes initial state
  - **Action item**: Examine other languages here.
- `ColumnAt`: No matching `RowAt`, because it's not directly needed. Should we add anyway?
  - `RowAt` is distinct from `ElementAt` in that it requires two levels of indices, whereas `ElementAt` does not.
  - **Action item**: Modify proposal to add `RowAt`.
- `Count`: simpler for what is currently a common usecase for `Filtered`.
  - Would be simplified by https://github.com/microsoft/qsharp-compiler/issues/557; defer introducing specializations until result of that discussion.
- `EmptyArray`: doesn't solve problem that `[]` can't be used as an empty array of a given type, while `[3]` works as an array literal.
  - **Action item**: File as issue on qsharp-compiler?
- `Interleaved`: Should it support more than two inputs?
  - If array-valued input, can use Transposed and Flattened together.
  - Possibly change to take array-valued input, but then implementation can check length of input and do a more efficient thing if possible.
  - **Action item**: Modify proposal to take array of inputs.
- `Unique`: doesn't work as written, would be slow if implemented.
  - **Action item**: Modify proposal to defer.
- `Zipped`/`Unzipped`: Not elegant to have `Unzipped3`, `Unzipped4`, ...
  - Also solvable by https://github.com/microsoft/qsharp-compiler/issues/557, but solution would be unwieldy.
- Missing: function to sort arrays.
  - **Action item**: Modify proposal to add `Sorted`.
- `ApplyTo{Head,Tail,Most,Rest}`
  - ApplyToHead is redundant with ApplyToFirstQubit, but worth keeping for consistency with Tail, Most, Rest.
- `ArrayAsHeadAndRest` / `ArrayAsMostAndTail`:
  - **Action item**: Rename to `HeadAndRest` / `MostAndTail`, move to `Microsoft.Quantum.Arrays`.
- Missing fact to check that a 2D array is rectangular / square.
  - `RectangularArrayFact` / `SquareArrayFact`?
  - **Action item**: Modify proposal to add these two facts.

**Other action items**:

- Comment on user feedback with ElementAt.
- Rename `folder` input on Microsoft.Quantum.Arrays.Folded to `fn` for consistency with `CumulativeFolded`.
- Discuss grouping removal of deprecation stubs, rather than splitting removals over multiple releases. Would 1.0 be a reasonable timeframe?

**Consensus**: Proceed as per above action items.

### Remove Microsoft.Quantum.Simulation.QuantumProcessor.Extensions

**Proposal**: The functions and operations in the Microsoft.Quantum.Simulation.QuantumProcessor.Extensions namespace are completely redundant with existing Q# standard library callables, and exist to help manage dependency ordering. It is proposed to move duplicated callables to the qsharp-runtime repo, and deprecate the copies in Microsoft.Quantum.Simulation.QuantumProcessor.Extensions.

**Discussion**:

- Does this affect the rewrite steps?
  - Yes, but we can port those rewrite steps, and should do so. It will take minimal effort to do so.
- Impact on samples/katas?
  - Minimal; only one sample uses at all.
- Observation: calling deprecated operations from deprecated operations still results in warnings (e.g.: ResultStack).
  - Proposal: don't raise deprecation warnings from within deprecated callables.
  - Action item: make feature request on qsharp-compiler.

**Consensus**: Proceed as proposed.

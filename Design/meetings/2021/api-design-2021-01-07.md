# Q# API Design Meeting / 7 January 2021

Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @guenp, @msoeken.

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/141
- https://github.com/microsoft/QuantumLibraries/issues/337
- https://github.com/microsoft/QuantumLibraries/pull/391
- https://github.com/microsoft/QuantumLibraries/pull/392

## Discussion

### ApplyWith2

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/141

**Discussion**:

- Haven't used existing `ApplyWith` since introduction of `within`/`apply`.
- Name `ApplyWith2` may be confusing.
- Canon namespace may not be the right place; too general for something useful but niche.
- Introducing into stdlib can force users to learn new operation names when reading code written against the library.
- **Action item**: Reply in thread with examples using `Delay`/`Delayed` and using `within`/`apply`.

**Consensus**: Useful, but niche and existing workarounds can be used.

<!--  -->

### Enhanced numerics

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/337

**Discussion**:

- Modifying (e.g.) `PrepareQI` to `PrepareQInt` could be helpful since "quantum integer" is a new enough concept users may be surprised.
  - Counterpoint: `QI` is more consistent with existing style guide and conventions.
  - Suffixes may be dropped if bounded polymorphism is adopted (https://github.com/microsoft/qsharp-compiler/issues/557). On the other hand, bounded polymorphism is a very large feature and is not planned for inclusion in Q# 1.0 (two years at _earliest_).
  - Spelling out full type names can improve searchability in docs.
  - For consistency, could require modifying other existing Q# APIs (e.g.: `AbsD` → `AbsDouble`).
  - **Action item**: Follow up with docs team to see if searchability can be improved with "tags" or search keywords.
  - **Action item**: Add detail to proposal to highlight changes and additions to documentation.

**Consensus**: Requires further discussion, continue review in February meeting.

<!--  -->

### ApplyUnitary

**Proposal**: https://github.com/microsoft/QuantumLibraries/pull/391

**Discussion**:

- API may need to change if/when multidimensional arrays are introduced (https://github.com/microsoft/qsharp-language/issues/39).
- M.Q.Synthesis right namespace.
- LittleEndian may not be best type for second input.
  - Resolved: needed as indices of matrix refer to integer representation of quantum states of second input.

**Consensus**: Approve API change as proposed.

<!--  -->

### Improved LocalUnivariateMinimum

**Proposal**: https://github.com/microsoft/QuantumLibraries/pull/392

**Discussion**:

- **Action item**: Double-check that no sample or kata code is broken. Also check MLADS notebooks to inform decision.
- Formally a breaking change, but relatively minor and mitigated by existing best practices.
- **Action item**: Capture breaking change in release notes and version number.

**Consensus**: Approve API change as proposed, ensuring that Jan release version number is bumped accordingly.

<!--  -->

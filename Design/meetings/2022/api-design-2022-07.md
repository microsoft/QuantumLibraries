# Q# API Design Discussions / July 2022

Reviewers (in order by username): @tcNickolas, @msoeken, @cgranade

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/598
- https://github.com/microsoft/QuantumLibraries/issues/601
- https://github.com/microsoft/QuantumLibraries/issues/602
- https://github.com/microsoft/QuantumLibraries/issues/607

## Discussion

### Fixed point truncation

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/598

**Reviews**:
* @tcNickolas, *approve*, left a comment
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**: Approved

---

### Table lookup

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/601

**Reviews**:
* @tcNickolas, *approve*, left comments
* @msoeken, *approve*, left comments
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**: Postponed, requires more discussion

---

### Windowed unitary application

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/602

**Reviews**:
* @tcNickolas, *approve*
* @msoeken, *approve*, left comments
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**: Approved as `ApplyToEachWindow` without `argumentTransform`, without `target` argument, with type parameter for input array, and with all 4 functor variants.

---

### New operation: Apply arithmetic function via table lookup

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/607

**Reviews**:
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.
* @msoeken, *approve*, left comments

**Consensus**: Approved after comments

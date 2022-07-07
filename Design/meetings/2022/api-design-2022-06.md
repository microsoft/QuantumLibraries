# Q# API Design Discussions / June 2022

Reviewers (in order by username): @tcNickolas, @cgranade

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/593
- https://github.com/microsoft/QuantumLibraries/issues/594
- https://github.com/microsoft/QuantumLibraries/issues/595

## Discussion

### Comparison operations

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/593

**Reviews**:
- @cgranade, *comment*: Does the `action` input to `ApplyControlledOnLessThanFxP` need to be adjointable? Should there be a separate `Adj` variant? We may also want to describe quickly how this proposal works with the outstanding proposal for numerics refactoring all-up (https://github.com/microsoft/QuantumLibraries/issues/337).
* @tcNickolas, *approve*, left a comment
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

---

### Functions for smallest and largest representable fixed point

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/594

**Reviews**: @cgranade, *approve*, same comment as above (https://github.com/microsoft/QuantumLibraries/issues/337).
* @tcNickolas, *approve*, left a comment
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

---

### Conversion functions for signed integers

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/595

**Reviews**:
* @cgranade: *comment*. We should coordinate this with the Q# language discussion around bit vectors in general. @swernli and @bettinaheim in particular may want to provide input here.
* @tcNickolas, *approve*, left a comment
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

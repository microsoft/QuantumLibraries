# Q# API Design Discussions / May 2022

Reviewers (in order by username): cgranade

## Agenda

- https://github.com/microsoft/QuantumLibraries/issues/579
- https://github.com/microsoft/qsharp-runtime/issues/1005
- https://github.com/microsoft/QuantumLibraries/issues/423

## Discussion

### Move Exp into libraries instead of intrinsics

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/579

**Reviews**:
- **Approved with suggestions.** I think this makes a lot of sense; Exp is extremely useful, but it's much higher-level than any other intrinsic, making it an odd fit for the Microsoft.Quantum.Intrinsics namespace. If I understand correctly, we did so due to the full-state simulator having some nice shortcuts for simulating `Exp`, but I don't think that implementation detail justifies including `Exp` in the basic set of intrinsics.
 I would suggest, however, that if we move `Exp` and make a breaking change that way, we also take the opportunity to give it a more clear name that reflects its action a bit better. In particular, `Exp` simulates evolution under a Pauli Hamiltonian such that `Microsoft.Quantum.Simulation.EvolveUnderPauli` may make sense; alternatively, it can be thought of as applying a joint rotation, perhaps justifying `Microsoft.Quantum.Canon.ApplyJointRotation`?

> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

---

### Add MeasureEachZ operation to runtime implementation

**Proposal**: https://github.com/microsoft/qsharp-runtime/issues/1005

**Reviews**:
- **Approved.** I think this is a pretty straightforward one, and works great for replacing the old `MultiM` with something more maintainable and easier to use.
> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

---

### Simple and consistent arithmetic API

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/423

**Reviews**:

- **Comment.**
  - I'm a bit concerned that names like `GreaterThan` read like functions (noun or adjective) rather than as operations.
  - How does it look to chain carry-out qubits to the `carryIn` of subsequent calls? The asymmetry that carry inputs are separate arguments while carry outputs are the MSB of LE registers feels like it could make that difficult to chain such that an example may be helpful.
  - I'd personally recommend staging this together with https://github.com/microsoft/QuantumLibraries/issues/337 to minimize the number of breaking changes and to make sure that both APIs are as consistent with each other as possible. #337 in particular feels like it could be helpful in making the input data types consistent not only across arithmetic APIs but other places where arithmetic data is used in general-purpose APIs.

> Please add a bullet point including your alias, your review result (*approve*, *reject*, *comment*), and a comment (optional when result is *approve*).  Alternatively, add a line to the PR discussion incl. a reference to this issue.

**Consensus**:

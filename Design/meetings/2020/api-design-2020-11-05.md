# Q# API Design Meeting / 9 October 2020 and 5 November 2020

Attendees (in order by username): @bettinaheim, @cgranade, @efratshabtai, @guenp, @msoeken.

## Agenda

- Improved representation of numeric data in quantum registers (https://github.com/microsoft/QuantumLibraries/issues/337)
- Improved `Microsoft.Quantum.Preparation` namespace (https://github.com/microsoft/QuantumLibraries/issues/344)

## Discussion

### Numeric data

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/337

**Discussion**:

- Modify `QInt` → `QSignedInt`? `ControlledOnInt` could be confusing otherwise.
- Add https://github.com/microsoft/qsharp-language/pull/41 as related issue.
- What are the disadvantages to making parallel Microsoft.Quantum.Numerics namespace? Mainly about explaining change to user.
- How much work to remove `BigEndian`? Probably not much, since most BE support already removed.
- **Action item**: Add detail to proposal about existing split between numerics package and arithmetic namespace.
- **Action item**: Explain more about what existing code is broken by this proposal.
- Fold into Q# standard library instead of Microsoft.Quantum.Numerics package? Probably not, lots of changes to tutorials, docs, etc. needed to support that.
  - On the other hand, if we provide language-level support for quantized plus, etc. operators, may be odd to have in a separate package. What are right operators, types?
- **Action item**: Add examples to the proposal from existing sample (before / after).
- **Action item**: Expand proposal to use full signatures (or at least more concrete examples) for new operations.
- **Action item**: Add to conceptual overview feedback to the end that users shouldn't need to know internal representation of integer data to use integer UDTs.

**Consensus**: Proposal seems reasonable, and should proceed after completing outstanding action items. Continue revising and follow up in next API review meeting.

### Microsoft.Quantum.Preparation

**Proposal**: https://github.com/microsoft/QuantumLibraries/issues/344

**Discussion**:

- **Action item**: Add suggestion for anonymous DUs (e.g.: `ComplexPolar[] | Double[]`) to DU proposal at https://github.com/microsoft/qsharp-compiler/issues/557.
    - See https://github.com/microsoft/qsharp-compiler/issues/406#issuecomment-723195125 for anonymous DU suggestion.
- **Action item**: Clarify in proposal that `PrepareQubit` is renamed to `PrepareSingleQubitPositivePauliEigenstate` rather than removed entirely.
- Is `PrepareSingleQubitPositivePauliEigenstate` too long? Can we shorten to `PreparePauliEigenstate`? Single-qubit is implied by there being a single eigenstate in the first place (e.g.: 𝑋𝑋 has a positive eigenspace of dimension 2, not 1).
- **Action item**: Add references to mixed state preparation algorithms to proposal to help provide context.
    - https://arxiv.org/pdf/1805.03662.pdf?page=15
- Clarify what's meant by "purified mixed state?" Refers to open systems theory, ρ = 𝑝ᵢ |ψᵢ⟩⟨ψᵢ| maps to |φ⟩ = √𝑝ᵢ |ψᵢ⟩ ⊗ |𝑖⟩.
- Generalize sign in mixed state preparation to allow more general data to be added to preparation; sign is then a particular example.
- **Action item**: Clarify that in the context of this proposal, "mixed state" is implicitly diagonal in the computational basis.

**Consensus**: Approved, modulo remaining action items. No further review required.

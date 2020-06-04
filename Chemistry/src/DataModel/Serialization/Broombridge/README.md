# Broombridge Serialization Design

The C# files in this folder implement serialization and deserialization for Broombridge files, using the YamlDotNet package.
In particular, the data model used for loading from and saving to YAML is implemented by the `BroombridgeDataStructuresv0.1.cs` and `BroombridgeDataStructuresv0.2.cs` files.
These files define internal classes for use with YamlDotNet, and that then get converted into instances of the public `ElectronicStructureProblem` classes.
For compatability with earlier versions of Broombridge, these files include full data models for historical Broombridge versions; where newer versions have not modified data structures, data models for later versions may have properties whose types are associated with earlier versions.
For instance, the `basis_set` property in Broombridge 0.2 is the same as it was in 0.1, such that the 0.2 `ProblemDescription` class has a `BasisSet` property whose type is `V0_1.BasisSet`.

This separation prevents users from depending on the exact details of how electronic structure problems are represented in each version of Broombridge, focusing instead on the _logical_ structure of problems once deserialized.

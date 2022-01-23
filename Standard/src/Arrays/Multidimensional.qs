// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Logical;
    open Microsoft.Quantum.Math;

    /// # Summary
    /// Returns the transpose of a matrix represented as an array
    /// of arrays.
    ///
    /// # Description
    /// Input as an $r \times c$ matrix with $r$ rows and $c$ columns.  The matrix
    /// is row-based, i.e., `matrix[i][j]` accesses the element at row $i$ and column $j$.
    ///
    /// This function returns the $c \times r$ matrix that is the transpose of the
    /// input matrix.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `matrix`.
    ///
    /// # Input
    /// ## matrix
    /// Row-based $r \times c$ matrix
    ///
    /// # Output
    /// Transposed $c \times r$ matrix
    ///
    /// # Example
    /// ```qsharp
    /// // same as [[1, 4], [2, 5], [3, 6]]
    /// let transposed = Transposed([[1, 2, 3], [4, 5, 6]]);
    /// ```
    function Transposed<'T>(matrix : 'T[][]) : 'T[][] {
        let numRows = Length(matrix);
        Fact(numRows > 0, "Matrix must have at least 1 row");
        let numColumns = Length(Head(matrix));
        Fact(numColumns > 0, "Matrix must have at least 1 column");
        RectangularArrayFact(matrix, "Matrix is not a rectangular array");

        return Mapped(ColumnAtUnchecked(_, matrix), SequenceI(0, numColumns - 1));
    }

    /// # Summary
    /// Returns the at the given index of an array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## index
    /// Index of element
    /// ## array
    /// The array being indexed.
    ///
    /// # Remark
    /// This function is more general than `LookupFunction`, since
    /// it can also be used for partial application on a fixed index.
    /// Note that the type parameter must explicitly be provided in
    /// this case as it cannot be deduced automatically.
    ///
    /// # Example
    /// Get the third number in four famous integer sequences. (note
    /// that the 0 index corresponds to the _first_ value of the sequence.)
    /// ```qsharp
    /// let lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76];
    /// let prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    /// let fibonacci = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34];
    /// let catalan = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862];
    /// let famous2 = Mapped(ElementAt<Int>(2, _), [lucas, prime, fibonacci, catalan]);
    /// // same as: famous2 = [3, 5, 1, 2]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.LookupFunction
    /// - Microsoft.Quantum.Arrays.ElementsAt
    function ElementAt<'T>(index : Int, array : 'T[]) : 'T {
        Fact(index >= 0 and index < Length(array), "Index is out of bound");
        return array[index];
    }

    /// # Summary
    /// Returns the array's elements at a given range
    /// of indices.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## range
    /// Range of array indexes
    /// ## array
    /// Array
    ///
    /// # Example
    /// Get the odd indexes in famous integer sequences. (note
    /// that the 0 index corresponds to the _first_ value of the sequence.)
    /// ```qsharp
    /// let lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76];
    /// let prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    /// let fibonacci = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34];
    /// let catalan = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862];
    /// let famousOdd = Mapped(ElementsAt<Int>(0..2..9, _), [lucas, prime, fibonacci, catalan]);
    /// // same as: famousOdd = [[2, 3, 7, 18, 47], [2, 5, 11, 17, 23], [0, 1, 3, 8, 21], [1, 2, 14, 132, 1430]]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.ElementAt
    /// - Microsoft.Quantum.Arrays.LookupFunction
    function ElementsAt<'T>(range : Range, array : 'T[]) : 'T[] {
        return array[range];
    }

    /// # Summary
    /// Extracts a column from a matrix.
    ///
    /// # Description
    /// This function extracts a column in a matrix in row-wise order.
    /// Extracting a row corresponds to element access of the first index
    /// and therefore requires no further treatment.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `matrix`.
    ///
    /// # Input
    /// ## column
    /// Column of the matrix
    /// ## matrix
    /// 2-dimensional matrix in row-wise order
    ///
    /// # Example
    /// ```qsharp
    /// let matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    /// let column = ColumnAt(0, matrix);
    /// // same as: column = [1, 4, 7]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Transposed
    /// - Microsoft.Quantum.Arrays.Diagonal
    function ColumnAt<'T>(column : Int, matrix : 'T[][]) : 'T[] {
        RectangularArrayFact(matrix, "Matrix is not a rectangular array");
        return ColumnAtUnchecked(column, matrix);
    }

    /// # Summary
    /// This function does not check for matrix shape
    ///
    /// # Description
    /// This function can be used in other multidimensional functions,
    /// which already check the input matrix for a valid rectangular shape.
    internal function ColumnAtUnchecked<'T>(column : Int, matrix : 'T[][]) : 'T[] {
        return Mapped(
                Compose(
                    ElementAt(column, _),
                    LookupFunction(matrix)
                ), RangeAsIntArray(IndexRange(matrix)));
    }

    /// # Summary
    /// Returns an array of diagonal elements of a 2-dimensional array
    ///
    /// # Description
    /// If the 2-dimensional array has not a square shape, the diagonal over
    /// the minimum over the number of rows and columns will be returned.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `matrix`.
    ///
    /// # Input
    /// ## matrix
    /// 2-dimensional matrix in row-wise order
    ///
    /// # Example
    /// ```qsharp
    /// let matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    /// let diagonal = Diagonal(matrix);
    /// // same as: column = [1, 5, 9]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Transposed
    function Diagonal<'T>(matrix : 'T[][]) : 'T[] {
        RectangularArrayFact(matrix, "Matrix is not a rectangular array");

        let numRows = Length(matrix);
        let numColumns = numRows == 0 ? 0 | Length(Head(matrix));

        return MappedOverRange(ElementAtDiagonal(_, matrix), 0..(MinI(numRows, numColumns) - 1));
    }

    internal function ElementAtDiagonal<'T>(index : Int, matrix : 'T[][]) : 'T {
        return matrix[index][index];
    }

    /// # Summary
    /// Represents a condition that a 2-dimensional array has a rectangular shape
    ///
    /// # Description
    /// This function asserts that each row in an array has the same length.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## array
    /// A 2-dimensional array of elements
    /// ## message
    /// A message to be printed if the array is not a rectangular array
    ///
    /// # Example
    /// ```qsharp
    /// RectangularArrayFact([[1, 2], [3, 4]], "Array is not rectangular");       // okay
    /// RectangularArrayFact([[1, 2, 3], [4, 5, 6]], "Array is not rectangular"); // okay
    /// RectangularArrayFact([[1, 2], [3, 4, 5]], "Array is not rectangular");    // will fail
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.SquareArrayFact
    function RectangularArrayFact<'T>(array : 'T[][], message : String) : Unit {
        if (Length(array) == 0) {
            return ();
        } else {
            let numColumns = Length(Head(array));
            for i in IndexRange(Rest(array)) {
                if Length(array[i+1]) != numColumns {
                    fail message;
                }
            }
            // qsharp-compiler Issue #964: QIR generation fails when passing a generic callable as a
            // parameter with an inherited type specifier. https://github.com/microsoft/qsharp-compiler/issues/964
            // if (Any(Compose(NotEqualI(numColumns, _), Length<'T>), Rest(array))) {
            //     fail message;
            // }
        }
    }

    /// # Summary
    /// Represents a condition that a 2-dimensional array has a square shape
    ///
    /// # Description
    /// This function asserts that each row in an array has
    /// as many elements as there are rows (elements) in the array.
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of each element of `array`.
    ///
    /// # Input
    /// ## array
    /// A 2-dimensional array of elements
    /// ## message
    /// A message to be printed if the array is not a square array
    ///
    /// # Example
    /// ```qsharp
    /// SquareArrayFact([[1, 2], [3, 4]], "Array is not a square");       // okay
    /// SquareArrayFact([[1, 2, 3], [4, 5, 6]], "Array is not a square"); // will fail
    /// SquareArrayFact([[1, 2], [3, 4, 5]], "Array is not a square");    // will fail
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.RectangularArrayFact
    function SquareArrayFact<'T>(array : 'T[][], message : String) : Unit {
        if (Length(array) == 0) {
            return ();
        } else {
            let numColumns = Length(array);
            for i in IndexRange(Rest(array)) {
                if Length(array[i+1]) != numColumns {
                    fail message;
                }
            }
            // qsharp-compiler Issue #964: QIR generation fails when passing a generic callable as a
            // parameter with an inherited type specifier. https://github.com/microsoft/qsharp-compiler/issues/964
            // if (Any(Compose(NotEqualI(numColumns, _), Length<'T>), array)) {
            //     fail message;
            // }
        }
    }
}

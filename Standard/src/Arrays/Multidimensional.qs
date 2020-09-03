// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Arrays {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Diagnostics;

    /// # Summary
    /// Returns the transpose of a matrix
    ///
    /// # Description
    /// Input as an $r \times c$ matrix with $r$ rows and $$ columns.  The matrix
    /// is row-based, i.e., `matrix[i][j]` access the element at row $i$ and column $j$.
    ///
    /// This function returns the $c \times r$ matrix that is the transpose of the
    /// input matrix.
    ///
    /// # Input
    /// ## matrix
    /// Row-based $r \times c$ matrix
    ///
    /// # Output
    /// Transposed $c \times r$ matrix
    ///
    /// # Example
    /// ```Q#
    /// // same as [[1, 4], [2, 5], [3, 6]]
    /// let transposed = Transposed([[1, 2, 3], [4, 5, 6]]);
    /// ```
    function Transposed<'T>(matrix : 'T[][]) : 'T[][] {
        let numRows = Length(matrix);
        Fact(numRows > 0, "Matrix must have at least 1 row");
        let numColumns = Length(Head(matrix));
        Fact(numColumns > 0, "Matrix must have at least 1 column");

        return Mapped(ColumnVectorAt(_, matrix), SequenceI(0, numColumns - 1));
    }

    /// # Summary
    /// Returns the array's element at given index
    ///
    /// # Input
    /// ## index
    /// Index of element
    /// ## array
    /// Array
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
    /// ```Q#
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
        Fact(index < Length(array), "Index is out of bound");
        return array[index];
    }

    /// # Summary
    /// Returns the array's elements at a given range
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
    /// ```Q#
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
    /// Extracts a column in a matrix
    ///
    /// # Dimension
    /// This function extracts a column in a matrix in row-wise order.
    /// Extracting a row corrsponds to element access of the first index
    /// and therefore requires no further treatment.
    ///
    /// # Input
    /// ## column
    /// Column of the matrix
    /// ## matrix
    /// 2-dimensional matrix in row-wise order
    ///
    /// # Example
    /// ```Q#
    /// let matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
    /// let column = ColumnVectorAt(0, matrix);
    /// // same as: column = [1, 4, 7]
    /// ```
    ///
    /// # See Also
    /// - Microsoft.Quantum.Arrays.Transposed
    function ColumnVectorAt<'T>(column : Int, matrix : 'T[][]) : 'T[] {
        return Mapped(
                Compose(
                    ElementAt<'T>(column, _),
                    LookupFunction(matrix)
                ), RangeAsIntArray(IndexRange(matrix)));
    }
}

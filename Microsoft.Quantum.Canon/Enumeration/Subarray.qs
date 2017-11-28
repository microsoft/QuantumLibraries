namespace Microsoft.Quantum.Canon {
    /// # Summary 
    /// The `Subarrary` function takes an array and a list of locations and 
    /// produces a new array formed from the elements of the original 
    /// array that match the given locations.
    ///
    /// # Remarks
    /// The function is defined for generic types, i.e., whenever we have 
    /// an array `'T[]` and a list of locations `Int[]` defining the subarray. 
    ///
    /// # Type Parameters
    /// ## 'T
    /// The type of `array` elements.
    ///
    /// # Input
    /// ## indices
    /// A list of integers that is used to define the subarray. 
    /// ## array
    /// An array of elements over `'T`.
    ///
    /// # Output
    /// An array `'T[]` of elements whose indices correspond to the subarray.
    ///
    /// # See also
    /// - @"Microsoft.Quantum.Canon.Permute": The functions `Permute` and `Subarray` are
    ///	  identical. (TODO: fix this redundancy?)
    function Subarray<'T>(indices : Int[], array : 'T[]) : 'T[] {
        let nSliced = Length(indices);
        mutable sliced = new 'T[nSliced];
        for( idx in 0..nSliced - 1 ) {
            set sliced[idx] = array[indices[idx]];
        }
        return sliced;
    }
}
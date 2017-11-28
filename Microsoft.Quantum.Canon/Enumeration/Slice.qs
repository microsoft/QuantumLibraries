namespace Microsoft.Quantum.Canon {
    /// # Summary 
    /// The `Slice` function takes an array and a list of locations and 
    /// produces a new array formed from the elements of the original 
    /// array that match the given locations.
    ///
    /// # Remark
    /// The function is defined for generic types, i.e., whenever we have 
    /// an array `'T[]` and a list of locations `Int[]` we can slice. 
    /// 
    /// # Input
    /// ## indices
    /// A list of integers that is used to define the slice. 
    /// ## array
    /// An array of elements over `'T`.
    /// # Output
    /// An array `'T[]` of elements whose indices correspond to the slice. 
    function Slice<'T>(indices : Int[], array : 'T[]) : 'T[] {
        let nSliced = Length(indices);
        mutable sliced = new 'T[nSliced];
        for( idx in 0..nSliced - 1 ) {
            set sliced[idx] = array[indices[idx]];
        }
        return sliced;
	}
}
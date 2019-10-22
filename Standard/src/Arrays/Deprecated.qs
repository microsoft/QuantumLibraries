// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Reversed" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Reversed")
    function Reverse<'T>(array : 'T[]) : 'T[] {
        return Reversed(array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Filtered" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Filtered")
    function Filter<'T> (predicate : ('T -> Bool), array : 'T[]) : 'T[] {
        return Filtered(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.All" instead.
    @Deprecated("Microsoft.Quantum.Arrays.All")
    function ForAll<'T> (predicate : ('T -> Bool), array : 'T[]) : Bool {
        return All(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Any" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Any")
    function ForAny<'T> (predicate : ('T -> Bool), array : 'T[]) : Bool {
        return Any(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Mapped" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Mapped")
    function Map<'T, 'U> (mapper : ('T -> 'U), array : 'T[]) : 'U[] {
        return Mapped(mapper, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.MappedByIndex" instead.
    @Deprecated("Microsoft.Quantum.Arrays.MappedByIndex")
    function MapIndex<'T, 'U> (mapper : ((Int, 'T) -> 'U), array : 'T[]) : 'U[] {
        return MappedByIndex(mapper, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Padded" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Padded")
    function Pad<'T> (nElementsTotal : Int, defaultElement : 'T, inputArray : 'T[]) : 'T[] {
        return Padded(nElementsTotal, defaultElement, inputArray);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Partitioned" instead.
    @Deprecated("Microsoft.Quantum.Arrays.Partitioned")
    function SplitArray<'T>(nElements: Int[], arr: 'T[]) : 'T[][] {
        return Partitioned(nElements, arr);
    }

}

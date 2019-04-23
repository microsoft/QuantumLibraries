// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Arrays;

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Reversed" instead.
    function Reverse<'T>(array : 'T[]) : 'T[] {
        _Renamed("Microsoft.Quantum.Canon.Reverse", "Microsoft.Quantum.Arrays.Reversed");
        return Reversed(array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Filtered" instead.
    function Filter<'T> (predicate : ('T -> Bool), array : 'T[]) : 'T[] {
        _Renamed("Microsoft.Quantum.Canon.Filter", "Microsoft.Quantum.Arrays.Filtered");
        return Filtered(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.All" instead.
    function ForAll<'T> (predicate : ('T -> Bool), array : 'T[]) : Bool {
        _Renamed("Microsoft.Quantum.Canon.ForAll", "Microsoft.Quantum.Arrays.All");
        return All(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Any" instead.
    function ForAny<'T> (predicate : ('T -> Bool), array : 'T[]) : Bool {
        _Renamed("Microsoft.Quantum.Canon.ForAny", "Microsoft.Quantum.Arrays.Any");
        return Any(predicate, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Mapped" instead.
    function Map<'T, 'U> (mapper : ('T -> 'U), array : 'T[]) : 'U[] {
        _Renamed("Microsoft.Quantum.Canon.Map", "Microsoft.Quantum.Arrays.Mapped");
        return Mapped(mapper, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.MappedByIndex" instead.
    function MapIndex<'T, 'U> (mapper : ((Int, 'T) -> 'U), array : 'T[]) : 'U[] {
        _Renamed("Microsoft.Quantum.Canon.MapIndex", "Microsoft.Quantum.Arrays.MappedByIndex");
        return MappedByIndex(mapper, array);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Padded" instead.
    function Pad<'T> (nElementsTotal : Int, defaultElement : 'T, inputArray : 'T[]) : 'T[] {
        _Renamed("Microsoft.Quantum.Canon.Pad", "Microsoft.Quantum.Arrays.Padded");
        return Padded(nElementsTotal, defaultElement, inputArray);
    }

    /// # Deprecated
    /// Please use @"Microsoft.Quantum.Arrays.Partitioned" instead.
    function SplitArray<'T>(nElements: Int[], arr: 'T[]) : 'T[][] {
        _Renamed("Microsoft.Quantum.Canon.SplitArray", "Microsoft.Quantum.Arrays.Partitioned");
        return Partitioned(nElements, arr);
    }

}

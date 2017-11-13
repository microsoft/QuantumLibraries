// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {

    /// # Summary
    /// Iterates a variable, say `arr`, through a Cartesian product
    /// [ 0, bounds[0]-1 ] × [ 0, bounds[1]-1 ] × [ 0, bounds[Length(bounds)-1]-1 ]
    /// and calls op(arr) for every element of the cartesian product
    operation IterateThroughCartesianProduct( bounds : Int[], op : ((Int[]) => ()) ) : ()
    {
        body
        {
            mutable arr = new Int[Length(bounds)];
            mutable finished = false;
            repeat
            {
                if( !finished )
                {
                    op(arr);
                }
            }
            until(finished)
            fixup
            {
                //computes the next element in the Cartesian product
                set arr[0] = arr[0] + 1;
                for( i in 0 .. Length(arr) -2 )
                {
                    if( arr[i] == bounds[i] )
                    {
                        set arr[i+1] = arr[i+1] + 1;
                        set arr[i] = 0;
                    }
                }
                if( arr[Length(arr) -1] == bounds[Length(arr) -1] )
                {
                    set finished = true;
                }
            }
        }
    }

    /// # Summary
    /// Creates an array of given length with all elments equal to given value.
    // FIXME: Make generic, move to Arrays.qb.
    function MakeConstArray( length : Int, value : Int ) : Int[]
    {
        mutable arr = new Int[length];
        for( i in 0 .. length - 1 )
        {
            set arr[i] = value;
        }
        return arr;
    }

    /// # Summary
    /// Iterates a variable, say arr, through cartesian product 
    /// [ 0, bound - 1 ] × [ 0, bound - 1 ] × [ 0, bound - 1 ]
    /// and calls op(arr) for every element of the cartesian product
    operation IterateThroughCartesianPower( power : Int, bound : Int, op : ((Int[]) => ()) ) : ()
    {
        body
        {
            IterateThroughCartesianProduct(MakeConstArray(power, bound), op);
        }
    }

}

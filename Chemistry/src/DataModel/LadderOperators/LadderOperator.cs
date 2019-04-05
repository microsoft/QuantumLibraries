// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using Microsoft.Quantum.Simulation.Core;

using System;
using System.Linq;
using System.Collections.Generic;
using Microsoft.Quantum.Chemistry;
using System.Numerics;

namespace Microsoft.Quantum.Chemistry
{

    public struct LadderOperator
    {
        public enum Type
        {
            u = 0, d = 1, identity
        }

        
        public Type type;
        public int index;
        

        public static (Type, int) AntiCommutation(LadderOperator x, LadderOperator y)
        {
            // {a_x, a_y^\dag} = \delta_{xy}
            // {a_x, a_y} = 1
            // {a_x^\dag, a_y^\dag} = 1
            if (x.index == y.index)
            {
                if(x.type != y.type)
                {
                    if(x.type == Type.d){
                        return (Type.identity, 1);
                    }
                    else
                    {
                        return (Type.identity, -1);
                    }
                }
                else
                {
                    return (Type.identity, 0);
                }
            }
            return (Type.identity, 0);
        }

        public LadderOperator(Type setType, int setIndex)
        {
            type = setType;
            index = setIndex;
        }

        public LadderOperator((Type, int) set) : this(set.Item1, set.Item2) { }

        #region Equality Testing
        public override bool Equals(object obj)
        {
            return (obj is LadderOperator x) ? Equals(x) : false;
        }

        public bool Equals(LadderOperator x)
        {
            return (type == x.type) && (index == x.index);
        }

        public override int GetHashCode()
        {
            return index * Enum.GetNames(typeof(Type)).Length + (int) type;
        }

        public static bool operator ==(LadderOperator x, LadderOperator y)
        {
            return x.Equals(y);
        }

        public static bool operator !=(LadderOperator x, LadderOperator y)
        {
            return !(x.Equals(y));
        }
#endregion
    }
    

}




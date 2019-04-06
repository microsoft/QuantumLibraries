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
        /// <summary>
        /// Enum for raising or lowering operator.
        /// </summary>
        public enum Type
        {
            u = 0, d = 1, identity
        }

        /// <summary>
        /// Type specifying raising or lowering operator.
        /// </summary>
        public Type type;

        /// <summary>
        /// System index operator acts on.
        /// </summary>
        public int index;
        
        /// <summary>
        /// Anti-commutation of ladder operators {x,y}.
        /// </summary>
        /// <param name="x">Left ladder operator.</param>
        /// <param name="y">Right ladder operator.</param>
        /// <returns>Result of {x,y}.</returns>
        [Obsolete]
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

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="setType">Set raising or lowering operator.</param>
        /// <param name="setIndex">Set system index.</param>

        public LadderOperator(Type setType, int setIndex)
        {
            type = setType;
            index = setIndex;
        }

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="set">Tuple where first item sets the operator type,
        /// and the second item indexes the system.</param>
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

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() {
            string op = type.ToString();
            return $"{index}{op}";
        }
    }
    

}




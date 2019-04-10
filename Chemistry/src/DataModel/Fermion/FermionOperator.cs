// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;
using Microsoft.Quantum.Chemistry.LadderOperators;

namespace Microsoft.Quantum.Chemistry.Fermion
{

        /*
    /// <summary>
    /// Data strcture for raising and lowering operators.
    /// </summary>
    public class FermionOperator : 
        LadderOperator<int>, 
        ILadderOperator<int>
    {

        #region Constructors
        public FermionOperator() : base() { }

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="setType">Set raising or lowering operator.</param>
        /// <param name="setIndex">Set system index.</param>
        public FermionOperator(RaisingLowering setType, int setIndex) : base(setType, setIndex) { }

        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="set">Tuple where first item sets the operator type,
        /// and the second item indexes the system.</param>
        public FermionOperator((RaisingLowering, int) set) : this(set.Item1, set.Item2) { }
        
        
        //* There are no current use cases where an implicit ladder operator is useful.
        //* All ladder operators invocations so far occur inside lists, and implicit
        //* conversion in a list violates type safety.
        //* 
        /// <summary>
        /// Implicit operator for creating a Ladder operator.
        /// </summary>
        /// <param name="setOperator">TUple where the first parameter
        /// is the raising or lowering index, and the second parameter
        /// is the position index of the ladder operator.</param>
        public static implicit operator FermionOperator((RaisingLowering, int) setOperator)
        {
            return new FermionOperator(setOperator.Item1, setOperator.Item2);
        }

        /// <summary>
        /// Implicit operator for creating a Ladder operator.
        /// </summary>
        /// <param name="setOperator">TUple where the first parameter
        /// is the raising or lowering index, and the second parameter
        /// is the position index of the ladder operator.</param>
        public static implicit operator (RaisingLowering, int)(FermionOperator fermionOperator)
        {
            return (fermionOperator.GetRaisingLowering(), fermionOperator.GetIndex());
        }
        
        #endregion
        */

   // }


}




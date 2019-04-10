// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System;

namespace Microsoft.Quantum.Chemistry.LadderOperators
{


    /// <summary>
    /// Data strcture for raising and lowering operators.
    /// </summary>
    public class LadderOperator<TIndex> : 
        IEquatable<LadderOperator<TIndex>>
        where TIndex : IEquatable<TIndex>
    {
        /// <summary>
        /// LadderType specifying raising or lowering operator.
        /// </summary>
        public RaisingLowering Type { get; set; }

        public RaisingLowering GetRaisingLowering() => Type;
            

        /// <summary>
        /// System index operator acts on.
        /// </summary>
        public TIndex Index { get; set; }

        #region Constructors
        public LadderOperator() { }

        public TIndex GetIndex() => Index;


        /// <summary>
        /// Constructor for ladder operator.
        /// </summary>
        /// <param name="setType">Set raising or lowering operator.</param>
        /// <param name="setIndex">Set system index.</param>
        public LadderOperator(RaisingLowering setType, TIndex setIndex)
        {
            Type = setType;
            Index = setIndex;
        }

        /// <summary>
        /// Implicit constructor for ladder operator.
        /// </summary>
        /// <param name="set">Tuple where first item sets the operator type,
        /// and the second item indexes the system.</param>
        public LadderOperator((RaisingLowering, TIndex) set) : this(set.Item1, set.Item2) { }

         #endregion

        #region Equality Testing
        public override bool Equals(object obj) => (obj is LadderOperator<TIndex> x) ? Equals(x) : false;

        public bool Equals(LadderOperator<TIndex> x) => (Type == x.Type) && (Index.Equals(x.Index));

        public override int GetHashCode() => Index.GetHashCode() * Enum.GetNames(typeof(RaisingLowering)).Length + (int) Type;

        public static bool operator ==(LadderOperator<TIndex> x, LadderOperator<TIndex> y) => x.Equals(y);

        public static bool operator !=(LadderOperator<TIndex> x, LadderOperator<TIndex> y) => !(x.Equals(y));
        #endregion

        /// <summary>
        /// Returns a human-readable description this object.
        /// </summary>
        public override string ToString() {
            string op = Type.ToString();
            return $"{Index.ToString()}{op}";
        }
    }
    

}




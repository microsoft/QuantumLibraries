// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

using System.Collections.Generic;
using System.Linq;

#nullable enable

namespace Microsoft.Quantum.Chemistry
{
    using System;

    public static partial class Extensions
    {

        #region Extension methods
        /// <summary>
        /// Computes the L_p norm of coefficients of all terms in a Hamiltonian.
        /// </summary>
        /// <param name="power">Selects type of norm.</param>
        /// <returns>L_p norm of Hamiltonian coefficients.</returns>
        public static double Norm(this IEnumerable<double> x, double power = 1.0)
        {
            double sum = x.Select(y => Math.Pow(Math.Abs(y), power)).Sum();
            return Math.Pow(sum, 1.0 / power);
        }

        /// <summary>
        /// Construct Double that implements the ITermValue interface. 
        /// </summary>
        /// <param name="x">Input double.</param>
        /// <returns>Double representing the input double.</returns>
        public static DoubleCoeff ToDoubleCoeff(this double x)
        {
            return new DoubleCoeff(x);
        }

        /// <summary>
        /// Computes `x^y` for an integer base `x` and exponent `y`
        /// </summary>
        /// <param name="x">Base.</param>
        /// <param name="exponent">Exponent.</param>
        /// <returns>An integer `x^y`</returns>
        public static int Pow(this int x, int exponent)
        {
            int output = 1;
            for (int i = 0; i < exponent; i++)
            {
                output *= x;
            }
            return output;
        }

        /// <summary>
        /// String representation of elements of an <see cref="IEnumerable{T}"/> of <see cref="int"/>.
        /// </summary>
        /// <param name="ints"><see cref="IEnumerable{T}"/> of <see cref="int"/>.</param>
        /// <returns>String representation of input elements.</returns>
        public static string ToCustomString(this IEnumerable<int> ints)
        {
            return "[" + string.Join(", ", ints) + "]";
        }


        internal static TValue GetValueOrDefault<TKey, TValue>(
                this Dictionary<TKey, TValue> dictionary, TKey key,
                TValue defaultValue = default
        ) =>
            dictionary.TryGetValue(key, out var value)
            ? value
            : defaultValue;

        #endregion


        #region Array Comparers
        /// <summary>
        /// Compares two equal-length sequences of integers by perform an element-wise
        /// comparison, starting from the first element.
        /// </summary>
        /// <param name="xArr">First sequence of elements.</param>
        /// <param name="yArr">Second sequence of elements.</param>
        /// <param name="comparer">Comparer used to define ordering of elements.</param>
        /// <returns>
        /// Returns <c>1</c> if <paramref name="xArr"/> is greater than <paramref name="yArr"/>.
        /// Returns <c>-1</c> if <paramref name="xArr"/> is less than <paramref name="yArr"/>.
        /// Returns <c>0</c> if <paramref name="xArr"/> is equal to <paramref name="yArr"/>.
        /// </returns>
        /// <example>
        /// CompareIntArray(new Int64[] {5}, new Int64[] {7}) == -1;
        /// CompareIntArray(new Int64[] {5,7}, new Int64[] {5,6}) == 1;
        /// CompareIntArray(new Int64[] {2,1,3}, new Int64[] {2,2,3}) == -1;
        /// </example>
        public static int CompareArray<TElement>(IEnumerable<TElement> xArr, IEnumerable<TElement> yArr, IComparer<TElement> comparer = null)
        {
            var useComparer = comparer == null ? Comparer<TElement>.Default : comparer;
            foreach (var item in xArr.Zip(yArr, (x, y) => (x, y)))
            {
                if (useComparer.Compare(item.y, item.x) > 0)
                {
                    return -1;
                }
                else if (useComparer.Compare(item.y, item.x) < 0)
                {
                    return 1;
                }
            }
            return 0;
        }


        /// <summary>
        /// Checks whether a sequence e.g. of integers {x,y,z,...} is sorted in
        /// ascending order x <= y <= z <= ...
        /// </summary>
        /// <param name="enumerable">Sequence of elements to be checked.</param>
        /// <param name="comparer">Comparer used to define ordering of elements.</param>
        /// <returns>
        /// <c>true</c> if the sequenceis sorted in ascending order, and <c>false</c> otherwise.
        /// </returns>
        public static bool IsInAscendingOrder<TElement>(this IEnumerable<TElement> enumerable, IComparer<TElement> comparer = null)
        {
            // return true for empty list.
            if (!enumerable.Any())
            {
                return true;
            }
            var curr = enumerable.First();
            foreach (var next in enumerable.Skip(1))
            {
                var useComparer = comparer == null ? Comparer<TElement>.Default : comparer;
                if (useComparer.Compare(next, curr) < 0)
                {
                    return false;
                }
                curr = next;
            }
            return true;
        }

        public class ArrayEqualityComparer<TElement> : IEqualityComparer<IEnumerable<TElement>>
        {
            IComparer<TElement> useComparer;
            public ArrayEqualityComparer(IComparer<TElement> comparer = null)
            {
                var useComparer = comparer == null ? Comparer<TElement>.Default : comparer;
            }

            public bool Equals(IEnumerable<TElement> x, IEnumerable<TElement> y)
            {
                return CompareArray(x, y, useComparer) == 0;
            }
            public int GetHashCode(IEnumerable<TElement> x)
            {
                int h = 19;
                foreach (var i in x)
                {
                    h = h * 31 + (i.GetHashCode());
                }
                return h;
            }
        }

        public class ArrayLexicographicComparer<TElement> : IComparer<IEnumerable<TElement>>
        {
            IComparer<TElement> useComparer;
            public ArrayLexicographicComparer(IComparer<TElement> comparer = null)
            {
                var useComparer = comparer == null ? Comparer<TElement>.Default : comparer;
            }
            public int Compare(IEnumerable<TElement> x, IEnumerable<TElement> y)
            {
                return CompareArray(x, y);
            }
        }


        #endregion


        #region Map
        /// <summary>
        ///      Given a value of an enumeration type, and an action for each
        ///      possible value of that enumeration type, performs the action
        ///      corresponding to the given value.
        /// </summary>
        public static void Map<E>
            (this E @enum, params (E, Action)[] actions)
            where E : struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    action();
                    return;
                }
            }
        }

        /// <summary>
        ///      Given a value of an enumeration type, and an function for each
        ///      possible value of that enumeration type, returns the value
        //       returned by the function corresponding to the given value.
        /// </summary>
        public static T Map<E, T>
            (this E @enum, params (E, Func<T>)[] actions)
            where E : struct, IConvertible
        {
            foreach (var (value, action) in actions)
            {
                if (@enum.Equals(value))
                {
                    return action();
                }
            }

            throw new ArgumentException($"Expected {@enum} to be a member of {@enum.GetType()}.");
        }

        #endregion

        /// <summary>
        /// Clones the values of an array, assuming that each element is a value type.
        /// </summary>
        /// <returns>A deep copy of this object.</returns>
        /// <typeparam name="T">Type of array to clone.</typeparam>
        /// <param name="array">Input array to clone.</param>
        /// <returns>Clone of the input array.</returns>
        public static TValue[] Clone<TValue>(this TValue[] array) => array.Select(el => el).ToArray();
        
        /// <summary>
        ///       Searches base types of a given type to find the type that immediately derives from
        ///       <see href="System.Object" />.
        /// </summary>
        internal static Type GetBasestType(this Type t) =>
            (t.BaseType == typeof(object) || t == typeof(object)) ? t : GetBasestType(t.BaseType);

        internal static int[] ToZeroBasedIndices(this IEnumerable<int> indices) =>
            indices.Select(idx => idx - 1).ToArray();

        internal static int[] ToOneBasedIndices(this IEnumerable<int> indices) =>
            indices.Select(idx => idx + 1).ToArray();

        internal static IEnumerable<TOutput> SelectMaybe<TInput, TOutput>(
                this IEnumerable<TInput> source, Func<TInput, TOutput?> selector
        ) 
        where TOutput : class =>
            source.Select(selector).Where(output => output != null).Select(item => item!);



    }

}

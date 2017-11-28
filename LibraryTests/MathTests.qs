// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Math;

	function NativeFnsAreCallableTest() : () {
		let arg = PI() / 2.0;
		AssertAlmostEqual(Sin(arg), 1.0);
		AssertAlmostEqual(Cos(arg), 0.0);

		let arcArg = 1.0;
		AssertAlmostEqual(ArcCos(arcArg), 0.0);
		AssertAlmostEqual(ArcSin(arcArg), arg);
	}

    function RealModTest() : () {
        AssertAlmostEqual(RealMod(5.5 * PI(), 2.0 * PI(), 0.0), 1.5 * PI());
        AssertAlmostEqual(RealMod(0.5 * PI(), 2.0 * PI(), -PI() / 2.0), 0.5 * PI());
    }

}

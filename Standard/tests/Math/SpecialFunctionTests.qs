// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

namespace Microsoft.Quantum.Canon {
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;

    @Test("QuantumSimulator")
    function LogFactorialDIsCorrect() : Unit {
        NearEqualityFactD(LogFactorialD(2000), 13206.524350513806756309398150028795369763937803753991516096227076);
        NearEqualityFactD(LogFactorialD(4321), 31856.241848248713008413364405861276737565269737465573546318165371);
    }

    // TODO: enable these tests when DoubleFactorialL goes public.
    // @Test("QuantumSimulator")
    // function DoubleFactorialLIsCorrect() : Unit {
    //     EqualityFactL(DoubleFactorialL(120), 9593444981835986954891939947669322185182489942608389896364094195294295395488811817369600000000000000L, "120!! was incorrect.");
    //     EqualityFactL(DoubleFactorialL(71), 2395415678676082004163677716234578672981800778515625L, "71!! was incorrect.");
    // }

    @Test("QuantumSimulator")
    function LogGammaDIsCorrect() : Unit {
        NearEqualityFactD(LogGammaD(3.14), 0.8261387047770284764910316081071750904374173003906590856150708884);
        NearEqualityFactD(LogGammaD(0.782), 0.1698067219140444477017554445169821679365849471502663300627878103);
        NearEqualityFactD(LogGammaD(1234.567), 7551.0278099842760398085493506933061185258592164059260052791257174);
    }

    @Test("QuantumSimulator")
    function FactorialIIsCorrect() : Unit {
        EqualityFactI(FactorialI(10), 3628800, "10! was incorrect.");
        EqualityFactI(FactorialI(6), 720, "6! was incorrect.");
        EqualityFactI(FactorialI(7), 5040, "7! was incorrect.");
    }

    function ApproximateFactorialCase(param : Int, expected : Double, multTolerance : Double) : Unit {
        let actual = ApproximateFactorial(param);
        EqualityWithinToleranceFact(actual / expected, 1.0, multTolerance);
    }

    @Test("QuantumSimulator")
    function ApproximateFactorialIsCorrect() : Unit {
        EqualityWithinToleranceFact(ApproximateFactorial(0), 1.0, 5e-2);
        EqualityWithinToleranceFact(ApproximateFactorial(1), 1.0, 5e-2);
        EqualityWithinToleranceFact(ApproximateFactorial(6), 720.0, 5e-2);
        EqualityWithinToleranceFact(ApproximateFactorial(7), 5040.0, 5e-2);
        EqualityWithinToleranceFact(ApproximateFactorial(10), 3628800.0, 5e-2);

        // For bigger cases, check multiplicative tolerance.
        ApproximateFactorialCase(12, 479001600.0, 0.001);
        ApproximateFactorialCase(72, 6.12344584e103, 0.001);
    }

    @Test("QuantumSimulator")
    function BinomIsCorrect() : Unit {
        EqualityFactI(Binom(31, 7), 2629575, "(31 7) was incorrect.");
        EqualityFactI(Binom(23, 9), 817190, "(23 9) was incorrect.");
        EqualityFactI(Binom(13, 5), 1287, "(13 5) was incorrect.");
        EqualityFactI(Binom(4, 0), 1, "(4 0) was incorrect.");
        EqualityFactI(Binom(4, 4), 1, "(4 4) was incorrect.");
    }

    @Test("QuantumSimulator")
    function FactorialLIsCorrectForPositiveInputs() : Unit {
        EqualityFactL(FactorialL(17), 355687428096000L, "17! was incorrect.");
        EqualityFactL(FactorialL(170), 7257415615307998967396728211129263114716991681296451376543577798900561843401706157852350749242617459511490991237838520776666022565442753025328900773207510902400430280058295603966612599658257104398558294257568966313439612262571094946806711205568880457193340212661452800000000000000000000000000000000000000000L, "170! was incorrect.");
        EqualityFactL(FactorialL(1700), 2998353205558418094701131115664613459850171438409583070169314496531140242247368919694363792466983858311481231384824780959878108306827711261613844950739492636665721700772411752144098241167488305729445936574215494519458901409182830930052432152888534372498719031822059189856486665368978846365916661190594744982844200865289620103875510740805972724042495289521016631330406109507024067707915350735127122360241028056095854097413440581325296450863778401319147282178998472716588603052487921303646737473677644631686179918257946503312668550160579326391433406171476262399591087772652113924101508920242485412020900808485404292120984052033241617306255312533207484832128428195892676872105110140425556455696529279622032604660056325281754536030524623417498267078586947004378177409019923967194915696628992157697278869086126927032833011711955548381226527524812138447954742569945776938385311777767433547406086141485401850467176039463883379028287886858613513479174940817288766757735829602002033350789525582106192132760910347342888017787564899411470636433664584027199788416867503748374613661323468015256415459848591528917807037879523611928200042173943044300636215362905927029758745268523183917902736812992812116021462911914158053109991970492273574063100978164267137822925925750887236251141544140977729807763929738072205138418727384401736398798293021947953818496159536819197920323874192790804534811682139886380738866046298387808445963982740817094044347446194416861486682286669874036824997312223696043267667524355036285105373699358134183372816096130534848906626018354319480185574360027946122190404332385803574864371070904464931872060834497981688550024104868739797741032409568429873434135618867622987851822414549096269190230820579964250910803704537565022307977457848649684204073383618515750592073448011954237678777827190916855790693831751140055032069717791825641145765238359129184568155465082398021387458757620934724321522795479371345125245515630995764560860862177231767696370263571179590907427725334652049718714691434415245678677487645441930874298204422913067827541850501885062300561585783037593795556535806675207057903884759500318075370777578492947381353837948967949286707525641972039074429323958151648800781566742997106427516689225841502764672105477104043382243019765946354325876422058090805188201790869918378645609051595809794403616441083516892525721992542946052592586463652255822964226148533904561064448625365949567932434140075846561831396505495737084626283083536723968347249473678956686493406788493952887051424049943753686860623346024699879773896264258525016468569296697745118183433124633920817140914110013999920362817995623302805517872734433183900413959556688890727800342538355059327508437096455153979928231506284255305592872121260922162226246308649832639732513126384972278355422436537394980375435960229288914481792794160724330203530033050984956555476368473697189674393388904870226316034762849494002436552890015086230398037141133837322469982044591734120914430399009404735953541265745148783679471071238882265006972882149508542806168239526629026007353045918384421414627048521770655466471170838327136349890231683971406138839537662225670750380294396864538729837113953538834662025033712532940785279960514669647190161970629158868087745892654708607537079992778103346600668423025477309980577879390401282485764463853318226596394201946144769908271691100310876682192884598748598309347481447633210989448429433214875408830118136815583374910637973420215675516805056868951875737630547161110315193297758192866984855044493168306811234085739956708142097794702335796035000426311590898319244138497431000093470301766845294541173428415111902926621489060722506907902037324215357470916109367344160884681872626809010762744056000299634611504744022217013138236225513177049543691804649273976160443897737847494498185929082013187196360218603243891156037661848297391553314018511623058189762002868068452658110796554902617412351624054876664743693168108605718383937799529227544669830826890532627267670454954296385596937136498550423995897097114268502182819280701871121326174020025039571810614711340873607591102437762508961805451573035509215030399220712619771588619461005461971756343142397328143334540870848854046912720665114562295663970751651201719187953580622932267380854578541830375207903653678094544651903436363109591764494064424751789517265887982423080062723879125692857698815574016000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000L, "1700! was incorrect.");
    }

    @Test("QuantumSimulator")
    function HalfIntegerBinomIsCorrect() : Unit {
        NearEqualityFactD(HalfIntegerBinom(7), 0.01611328125);
        NearEqualityFactD(HalfIntegerBinom(31), 0.001654486661428540561491473681599018163979053497314453125);
    }

}

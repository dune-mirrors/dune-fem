#ifndef DUNE_FEM_MODIFIEDNEWTONCOTES_IMPLEMENTATION_HH
#define DUNE_FEM_MODIFIEDNEWTONCOTES_IMPLEMENTATION_HH

#include "gausspoints.hh"

namespace Dune
{

  namespace Fem
  {

    inline ModifiedNewtonCotes :: ModifiedNewtonCotes ()
      : QuadPtsBase( MAXP )
    {
      int m = 0;
      O[m] = -1;

      m = 1;
      G[m][0] = 0.5;
      W[m][0] = 1.0;
      O[m] = 1;

      m = 2;
      G[m][0] = 0.25;
      G[m][1] = 0.75;
      W[m][0] = 0.5;
      W[m][1] = 0.5;
      O[m] = 2;

      m = 3;
      G[m][0] = 0.166666666666666666666666666666666666666666666666666666666666666666;
      G[m][1] = 0.5;
      G[m][2] = 0.833333333333333333333333333333333333333333333333333333333333333334;
      //W[m][0] = 0.3749999999999997783023397701641302596498899930650837633629968379;
      //W[m][1] = 0.2500000000000004392319841173274051234505041792370948406966572771;
      //W[m][2] = 0.3749999999999997783023397701641302596498899930650837633629968379;
      W[m][0] = 0.375;
      W[m][1] = 0.25;
      W[m][2] = 0.375;
      O[m] = 3;

      m = 4;
      G[m][0] = 0.125;
      G[m][1] = 0.375;
      G[m][2] = 0.625;
      G[m][3] = 0.875;
      //W[m][0] = 0.2708333333333333115704388165333897852564625667803215258049242423;
      //W[m][1] = 0.2291666666666666875884831345081582901801801089084509647253787880;
      //W[m][2] = 0.2291666666666666875884831345081582901801801089084509647253787880;
      //W[m][3] = 0.2708333333333333115704388165333897852564625667803215258049242423;
      W[m][0] = 0.2708333333333333333333333333333333333333333333333333333333333334;
      W[m][1] = 0.2291666666666666666666666666666666666666666666666666666666666666;
      W[m][2] = 0.2291666666666666666666666666666666666666666666666666666666666666;
      W[m][3] = 0.2708333333333333333333333333333333333333333333333333333333333334;
      O[m] = 4;

      m = 5;
      G[m][0] = 0.1;
      G[m][1] = 0.3;
      G[m][2] = 0.5;
      G[m][3] = 0.7;
      G[m][4] = 0.9;
      W[m][0] = 0.2387152777777788056375025915626439460696655368777231098503898678;
      W[m][1] = 0.08680555555555111847922462128482839513603995502075648650891286466;
      W[m][2] = 0.3489583333333401433260932639018047916720023338236316409380052477;
      W[m][3] = 0.08680555555555111847922462128482839513603995502075648650891286466;
      W[m][4] = 0.2387152777777788056375025915626439460696655368777231098503898678;
      O[m] = 5;

      m = 6;
      G[m][0] = 0.08333333333333333;
      G[m][1] = 0.25;
      G[m][2] = 0.41666666666666663;
      G[m][3] = 0.5833333333333334;
      G[m][4] = 0.75;
      G[m][5] = 0.9166666666666666;
      W[m][0] = 0.1929687500000000013847540667623024057526432086974825455366484308;
      W[m][1] = 0.1085937500000000797992249010895169245009405066912896530066483272;
      W[m][2] = 0.1984374999999999192022829688934024943142733024164833667399857226;
      W[m][3] = 0.1984374999999999192022829688934024943142733024164833667399857226;
      W[m][4] = 0.1085937500000000797992249010895169245009405066912896530066483272;
      W[m][5] = 0.1929687500000000013847540667623024057526432086974825455366484308;
      O[m] = 6;

      m = 7;
      G[m][0] = 0.07142857142857142;
      G[m][1] = 0.2142857142857143;
      G[m][2] = 0.3571428571428572;
      G[m][3] = 0.5;
      G[m][4] = 0.6428571428571429;
      G[m][5] = 0.7857142857142858;
      G[m][6] = 0.9285714285714286;
      W[m][0] = 0.1790002893518445971649372821051280941351976494768261140343260599;
      W[m][1] = 0.006380208333377521015320139693810201727906493521771217358988132379;
      W[m][2] = 0.4051432291665542717164287055762303838130447443895170311477506895;
      W[m][3] = -0.1810474537035528045415947699275381783042663091407185772232933536;
      W[m][4] = 0.4051432291665542717164287055762303838130447443895170311477506895;
      W[m][5] = 0.006380208333377521015320139693810201727906493521771217358988132379;
      W[m][6] = 0.1790002893518445971649372821051280941351976494768261140343260599;
      O[m] = 7;

      m = 8;
      G[m][0] = 0.0625;
      G[m][1] = 0.1875;
      G[m][2] = 0.3125;
      G[m][3] = 0.4375;
      G[m][4] = 0.5625;
      G[m][5] = 0.6875;
      G[m][6] = 0.8125;
      G[m][7] = 0.9375;
      W[m][0] = 0.1527503926917994333110139110397758603228589968808048027591302014;
      W[m][1] = 0.03685567542989155415424477213309793456504446866282408262526113435;
      W[m][2] = 0.2437639508928620613458099888235322048854978318946847990746808363;
      W[m][3] = 0.06662998098544694983085463715945480772427634642016934288730688334;
      W[m][4] = 0.06662998098544694983085463715945480772427634642016934288730688334;
      W[m][5] = 0.2437639508928620613458099888235322048854978318946847990746808363;
      W[m][6] = 0.03685567542989155415424477213309793456504446866282408262526113435;
      W[m][7] = 0.1527503926917994333110139110397758603228589968808048027591302014;
      O[m] = 8;

      m = 9;
      G[m][0] = 0.05555555555555555;
      G[m][1] = 0.16666666666666666;
      G[m][2] = 0.2777777777777778;
      G[m][3] = 0.38888888888888884;
      G[m][4] = 0.5;
      G[m][5] = 0.6111111111111112;
      G[m][6] = 0.7222222222222222;
      G[m][7] = 0.8333333333333333;
      G[m][8] = 0.9444444444444444;
      W[m][0] = 0.1451278250558263581123826914500552046686621766752494114420174508;
      W[m][1] = -0.04548130580375798865456098187886034148580562861680329758851389815;
      W[m][2] = 0.5062688337060229530505184317051100327509695056244393372180093760;
      W[m][3] = -0.5627887834834912483498188998354131185300678863530942438743892934;
      W[m][4] = 0.9137468610507998475455952640389667849117671704543959668383070226;
      W[m][5] = -0.5627887834834912483498188998354131185300678863530942438743892934;
      W[m][6] = 0.5062688337060229530505184317051100327509695056244393372180093760;
      W[m][7] = -0.04548130580375798865456098187886034148580562861680329758851389815;
      W[m][8] = 0.1451278250558263581123826914500552046686621766752494114420174508;
      O[m] = 9;

      m = 10;
      G[m][0] = 0.05;
      G[m][1] = 0.15;
      G[m][2] = 0.25;
      G[m][3] = 0.35;
      G[m][4] = 0.44999999999999996;
      G[m][5] = 0.5499999999999999;
      G[m][6] = 0.65;
      G[m][7] = 0.75;
      G[m][8] = 0.85;
      G[m][9] = 0.95;
      W[m][0] = 0.1278639428409557892417565762028831157751716134235817613410369588;
      W[m][1] = -0.01326074886155737724154929435336486457217401837489202425911590484;
      W[m][2] = 0.3302022405823394108301385935919416367634035498558943128345843991;
      W[m][3] = -0.1688483236538980141900724423136292583376997100092413468376685321;
      W[m][4] = 0.2240428890921601965079773208525276407403592555489023238112887705;
      W[m][5] = 0.2240428890921601965079773208525276407403592555489023238112887705;
      W[m][6] = -0.1688483236538980141900724423136292583376997100092413468376685321;
      W[m][7] = 0.3302022405823394108301385935919416367634035498558943128345843991;
      W[m][8] = -0.01326074886155737724154929435336486457217401837489202425911590484;
      W[m][9] = 0.1278639428409557892417565762028831157751716134235817613410369588;
      O[m] = 10;

    }
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GAUSSPOINTS_IMPLEMENTATION_HH

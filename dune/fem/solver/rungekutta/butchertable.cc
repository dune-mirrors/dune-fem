#include <config.h>

#include <dune/fem/solver/rungekutta/butchertable.hh>

namespace DuneODE
{

  // implicit34ButcherTable
  // ----------------------

  // R. Alexander:  Diagonally Implicit Runge-Kutta Methods for Stiff ODEs (1977)
  static const double dirk34_alpha = 2.*std::cos(M_PI/18.)/std::sqrt(3.);
  static const double dirk34_alpha2 = dirk34_alpha * dirk34_alpha;
  static const double DIRK34_A[] =
    {(1.+dirk34_alpha)*0.5, 0., 0.,
     -0.5*dirk34_alpha, (1.+dirk34_alpha)*0.5, 0.,
     1+dirk34_alpha, -(1+2*dirk34_alpha), (1.+dirk34_alpha)*0.5};
  static const double DIRK34_b[] =
    {1./(6.*dirk34_alpha2), 1.-1./(3.*dirk34_alpha2),
    1./(6.*dirk34_alpha2)};
  static const double DIRK34_c[] =
    {(1.+dirk34_alpha)*0.5, 0.5, (1.-dirk34_alpha)*0.5};

  SimpleButcherTable< double > implicit34ButcherTable ()
  {
    return SimpleButcherTable< double >( 3, 4, DIRK34_A, DIRK34_b, DIRK34_c );
  }



  // implicit3ButcherTable
  // ---------------------

  static const double delta_dirk = 1.0/2.0 + sqrt(3.0)/6.0;
  static const double DIRK3_A[] =
    {delta_dirk, 0.0,
     1.0-2.0*delta_dirk, delta_dirk
    };
  static const double DIRK3_b[] =
    {(0.5-delta_dirk)/(1.0-2.0*delta_dirk), (0.5-delta_dirk)/(1.0-2.0*delta_dirk)};
  static const double DIRK3_c[] =
    {delta_dirk, 1.0-delta_dirk};

  SimpleButcherTable< double > implicit3ButcherTable ()
  {
    return SimpleButcherTable< double >( 2, 3, DIRK3_A, DIRK3_b, DIRK3_c );
  }



  // implicitEulerButcherTable
  // -------------------------

  static const double ImplicitEuler_A[] = {1.0};
  static const double ImplicitEuler_b[] = {1.0};
  static const double ImplicitEuler_c[] = {1.0};

  SimpleButcherTable< double > implicitEulerButcherTable ()
  {
    return SimpleButcherTable< double >( 1, 1, ImplicitEuler_A, ImplicitEuler_b, ImplicitEuler_c );
  }



  // gauss2ButcherTable
  // ------------------

  //class Gauss2 (Crank-Nicholson)
  static const double Gauss2_A[] = {0.5};
  static const double Gauss2_b[] = {1.0};
  static const double Gauss2_c[] = {0.5};

  SimpleButcherTable< double > gauss2ButcherTable ()
  {
    return SimpleButcherTable< double >( 1, 2, Gauss2_A, Gauss2_b, Gauss2_c );
  }



  // semiImplicitEulerButcherTable
  // -----------------------------

  // semi implicit Euler, 1 stage, 1st order
  // implicit part
  static const double SIEuler_A[] = {1.0};
  static const double SIEuler_c[] = {1.0};
  // explicit part
  static const double SIEuler_Aex[] = {0.0};
  static const double SIEuler_cex[] = {0.0};
  // common part
  static const double SIEuler_b[] = {1.0};

  SimpleButcherTable< double > semiImplicitEulerButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 1, 1, expl ? SIEuler_Aex : SIEuler_A, SIEuler_b, expl ? SIEuler_cex : SIEuler_c );
  }



  // semiImplicit23ButcherTable
  // --------------------------

  // semi implicit 3 stages, 2nd order
  // implicit part
  static const double SIRK23_A[] =
    {0.5, 0.0, 0.0,
     -1.0, 0.5, 0.0,
     0.25, 0.25, 0.5
    };
  static const double SIRK23_c[] =
    {0.5, -0.5, 1.0};
  // explicit part
  static const double SIRK23_Aex[] =
    {0.0, 0.0, 0.0,
     1.0, 0.0, 0.0,
     0.0, 0.5, 0.0
    };
  static const double SIRK23_cex[] =
    {0.0, 1.0, 0.5};
  // common part
  static const double SIRK23_b[] =
    {0.25, 0.25, 0.5};

  SimpleButcherTable< double > semiImplicit23ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 3, 2, expl ? SIRK23_Aex : SIRK23_A, SIRK23_b, expl ? SIRK23_cex : SIRK23_c );
  }



  // semiImplicit33ButcherTable
  // --------------------------

  // SIRK33, 3 stages, 3rd order
  // YZ33 from Dennis Diss
  // implicit part
  static const double SIRK33_A[] =
    {3.0/4.0, 0.0, 0.0,
     5589.0/6524.0, 75.0/233.0, 0.0,
     7691.0/26096.0, -26335.0/78288.0, 65.0/168.0
    };
  static const double SIRK33_c[] =
    {3.0/4.0,
     5589.0/6524.0 + 75.0/233.0,
     7691.0/26096.0 - 26335.0/78288.0 + 65.0/168.0
    };
  // explicit part
  static const double SIRK33_Aex[] =
    {0.0, 0.0, 0.0,
     8.0/7.0, 0.0, 0.0,
     71.0/252.0, 7.0/36.0, 0.0
    };
  static const double SIRK33_cex[] =
    {0.0, 8.0/7.0, 71.0/252.0 + 7.0/36.0};
  // common part
  static const double SIRK33_b[] =
    {1.0/8.0, 1.0/8.0, 3.0/4.0};

  SimpleButcherTable< double > semiImplicit33ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 3, 3, expl ? SIRK33_Aex : SIRK33_A, SIRK33_b, expl ? SIRK33_cex : SIRK33_c );
  }



  // semiImplicitSSP222ButcherTable
  // ------------------------------

  // IMEX_SSP222, 2 stages, 2nd order
  // implicit part
  static const double delta = 1.0 - 1.0/sqrt(2.0);
  static const double IMEX_SSP222_A[] =
    {delta, 0.0,
     1.0-2.0*delta, delta
    };
  static const double IMEX_SSP222_c[] =
    {delta, 1.0-delta};
  // explicit part
  static const double IMEX_SSP222_Aex[] =
    {0.0, 0.0,
     1.0, 0.0
    };
  static const double IMEX_SSP222_cex[] =
    {0.0, 1.0};
  // common part
  static const double IMEX_SSP222_b[] =
    {0.5, 0.5};

  SimpleButcherTable< double > semiImplicitSSP222ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 2, 2, expl ? IMEX_SSP222_Aex : IMEX_SSP222_A, IMEX_SSP222_b, expl ? IMEX_SSP222_cex : IMEX_SSP222_c );
  }


#if 0 // these are not working

  // semiImplicitARK34ButcherTable
  // -----------------------------

  // ARK34, 4 stages, 3rd order
  // Kennedy and Carpenter, Appl. Num. Math. 44, 2003
  // explicit part
  static const double ARK34_Aex[] =
    {
      0.0,  0.0,  0.0,  0.0,
      1767732205903.0/20278366411180.0,  0.0,  0.0,  0.0,
      5535828885825.0/10492691773637.0,  788022342437.0/10882634858940.0,  0.0, 0.0,
      6485989280629.0/16251701735622.0,  -4246266847089.0/9704473918619.0, 10755448449292.0/10357097424841.0, 0.0,
    };

  // implicit part
  static const double ARK34_A[] =
    {
      0.0,  0.0,  0.0,  0.0,
      1767732205903.0/4055673282236.0,  1767732205903.0/4055673282236.0,  0.0,  0.0,
      2746238789719.0/10658868560708.0, -640167445237.0/6845629431997.0,  1767732205903.0/4055673282236.0, 0.0,
      1471266399579.0/7840856788654.0,  -4482444167858.0/7529755066697.0, 11266239266428.0/11593286722821.0, 1767732205903.0/4055673282236.0
    };

  // common part
  static const double ARK34_b[] =
    { 1471266399579.0/7840856788654.0, -4482444167858.0/7529755066697.0, 11266239266428.0/11593286722821.0, 1767732205903.0/4055673282236.0 };
    //{ 2756255671327.0/12835298489170.0, -10771552573575.0/22201958757719.0, 9247589265047.0/10645013368117.0, 2193209047091.0/5459859503100.0 };
  static const double ARK34_c[] =
    {0.0, 1767732205903.0/2027836641118.0 , 3.0/5.0, 1.0};

  SimpleButcherTable< double > semiImplicitARK34ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 4, 3, expl ? ARK34_Aex : ARK34_A, ARK34_b, ARK34_c );
  }



  // semiImplicitARK46ButcherTable
  // -----------------------------

  // ARK46, 6 stages, 4rd order
  // Kennedy and Carpenter, Appl. Num. Math. 44, 2003
  // explicit part
  static const double ARK46_Aex[] =
    {
      0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
      0.5,  0.0, 0.0, 0.0, 0.0, 0.0,
      13861.0/62500.0, 6889.0/62500.0, 0.0, 0.0, 0.0, 0.0,
      -116923316275.0/2393684061468.0, -2731218467317.0/15368042101831.0, 9408046702089.0/11113171139209.0, 0.0, 0.0, 0.0,
      -451086348788.0/2902428689909.0, -2682348792572.0/7519795681897.0, 12662868775082.0/11960479115383.0, 3355817975965.0/11060851509271.0, 0.0, 0.0,
      647845179188.0/3216320057751.0, 73281519250.0/8382639484533.0, 552539513391.0/3454668386233.0, 3354512671639.0/8306763924573.0, 4040.0/17871.0, 0.0
    };
  static const double ARK46_cex[] = {
      0.0, 0.5, 13861.0/62500.0 + 6889.0/62500.0,
      -116923316275.0/2393684061468.0 -2731218467317.0/15368042101831.0 + 9408046702089.0/11113171139209.0,
      -451086348788.0/2902428689909.0 -2682348792572.0/7519795681897.0 + 12662868775082.0/11960479115383.0 + 3355817975965.0/11060851509271.0,
      647845179188.0/3216320057751.0 + 73281519250.0/8382639484533.0 + 552539513391.0/3454668386233.0 + 3354512671639.0/8306763924573.0 + 4040.0/17871.0
  };

  // implicit part
  static const double ARK46_A[] =
    {
      0.0,  0.0,  0.0, 0.0, 0.0, 0.0,
      0.25, 0.25, 0.0, 0.0, 0.0, 0.0,
      8611.0/62500.0, -1743.0/31250.0, 0.25 , 0.0, 0.0, 0.0,
      5012029.0/34652500.0, -654441.0/2922500.0, 174375.0/388108.0, 0.25, 0.0, 0.0,
      15267082809.0/155376265600.0, -71443401.0/120774400.0, 730878875.0/902184768.0, 2285395.0/8070912.0, 0.25, 0.0,
      82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 0.25
    };

  // common part
  static const double ARK46_b[] =
    { 82889.0/524892.0, 0.0, 15625.0/83664.0, 69875.0/102672.0, -2260.0/8211.0, 0.25 };
  static const double ARK46_c[] =
    {0.0, 0.5 , 8611.0/62500.0 -1743.0/31250.0 + 0.25,
     5012029.0/34652500.0 -654441.0/2922500.0 + 174375.0/388108.0 + 0.25,
     15267082809.0/155376265600.0 -71443401.0/120774400.0 + 730878875.0/902184768.0 + 2285395.0/8070912.0 + 0.25,
     82889.0/524892.0 + 15625.0/83664.0 + 69875.0/102672.0  -2260.0/8211.0 + 0.25
    };

  SimpleButcherTable< double > semiImplicitARK46ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 6, 4, expl ? ARK46_Aex : ARK46_A, ARK46_b, expl ? ARK46_cex : ARK46_c );
  }

#endif

  // semiImplicitIERK45ButcherTable
  // ------------------------------

  // IERK45 5 stages, 4rd order
  // Lindblad et. al. ECCOMAS CFD 2006
  // explicit part
  static const double IERK45_Aex[] =
    {
      0.0,  0.0,  0.0,  0.0,  0.0,
      0.39098372452428, 0.0,  0.0,  0.0,  0.0,
      1.09436646160460, 0.33181504274704, 0.0, 0.0, 0.0,
      0.14631668003312, 0.69488738277516, 0.46893381306619, 0.0, 0.0,
      -1.33389883143642, 2.90509214801204, -1.06511748457024, 0.27210900509137, 0.0
    };
  static const double IERK45_cex[] =
    { 0.0, 0.39098372452428, IERK45_Aex[ 10 ] + IERK45_Aex[ 11 ],
      IERK45_Aex[ 15 ] + IERK45_Aex[ 16 ] + IERK45_Aex[ 17 ],
      IERK45_Aex[ 20 ] + IERK45_Aex[ 21 ] + IERK45_Aex[ 22 ] + IERK45_Aex[ 23 ]
    };

  // implicit part
  static const double IERK45_A[] =
    {
      0.25, 0.0,  0.0,  0.0,  0.0,
      0.34114705729739, 0.25, 0.0,  0.0,  0.0,
      0.80458720789763, -0.07095262154540, 0.25,  0.0,  0.0,
      -0.52932607329103, 1.15137638494253, -0.80248263237803, 0.25, 0.0,
      0.11933093090075, 0.55125531344927, -0.1216872844994, 0.20110104014943, 0.25
    };
  static const double IERK45_c[] =
    {
      0.25, 0.34114705729739 + 0.25, 0.80458720789763  -0.07095262154540 + 0.25,
      -0.52932607329103 + 1.15137638494253 -0.80248263237803+ 0.25,
      0.11933093090075 + 0.55125531344927  -0.1216872844994 + 0.20110104014943 + 0.25
    };


  // common part
  static const double IERK45_b[] =
    { IERK45_A[20] , IERK45_A[21], IERK45_A[22], IERK45_A[23] , IERK45_A[24] };


  SimpleButcherTable< double > semiImplicitIERK45ButcherTable ( bool expl )
  {
    return SimpleButcherTable< double >( 5, 4, expl ? IERK45_Aex : IERK45_A, IERK45_b, expl ? IERK45_cex : IERK45_c );
  }

  // row2ButcherTable
  // ---------------------
  static const double ROW2_A[] =
    {0.5};
  static const double ROW2_b[] =
    {1.};
  static const double ROW2_c[] =
    {0.};
  static const double ROW2_B[] =
    {0.0};
  ROWSimpleButcherTable< double > row2ButcherTable ()
  {
    return ROWSimpleButcherTable< double >( 1, 2, ROW2_A, ROW2_b, ROW2_c, ROW2_B );
  }

  // row3ButcherTable
  // ---------------------
  static const double delta_row = 1.0/2.0 + sqrt(3.0)/6.0;
  static const double ROW3_A[] =
    {delta_row, 0.0,
     -4./3.*delta_row, delta_row
    };
  static const double ROW3_b[] =
    {1./4.,3./4.};
  static const double ROW3_c[] =
    {0,2./3.};
  static const double ROW3_B[] =
    {0.0, 0.0,
     2./3., 0.0
    };

  ROWSimpleButcherTable< double > row3ButcherTable ()
  {
    return ROWSimpleButcherTable< double >( 2, 3, ROW3_A, ROW3_b, ROW3_c, ROW3_B );
  }

#if 0
  static const double delta_row = 7.8867513459481287e-01;
  static const double ROW3_A[] =
    {delta_row, 0.0, 0.0,
     -1.5773502691896257e+00, delta_row, 0.0,
     -6.7075317547305480e-01, -1.7075317547305482e-01, delta_row
    };
  static const double ROW3_b[] =
    {1.0566243270259355e-01, 4.9038105676657971e-02, 8.4529946162074843e-01};
  static const double ROW3_c[] =
    {0,0};
  static const double ROW3_B[] =
    {0.0, 0.0, 0.0,
     1.5773502691896257e+00, 0.0, 0.0,
     5.0000000000000000e-01, 0.0, 0.0 
    };

  ROWSimpleButcherTable< double > row3ButcherTable ()
  {
    return ROWSimpleButcherTable< double >( 3, 3, ROW3_A, ROW3_b, ROW3_c, ROW3_B );
  }
#endif


} // namespace DuneODE

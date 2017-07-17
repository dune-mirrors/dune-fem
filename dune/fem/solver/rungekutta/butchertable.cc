#include <config.h>

#include <dune/fem/solver/rungekutta/butchertable.hh>

namespace DuneODE
{
  //////////////////////////////////////////////////////
  //
  //  Explicit Butcher Tables
  //
  //////////////////////////////////////////////////////

  // explicitEulerButcherTable
  // -------------------------

  SimpleButcherTable< double > explicitEulerButcherTable ()
  {
    static const double A[] = {0.0};
    static const double b[] = {1.0};
    static const double c[] = {1.0};

    return SimpleButcherTable< double >( 1, 1, A, b, c );
  }


  // TVD2ButcherTable (Heun)
  // -----------------------

  SimpleButcherTable< double > tvd2ButcherTable ()
  {
    static const double A[] = {0.0, 0.0,
                               1.0, 0.0};
    static const double b[] = {0.5, 0.5};
    static const double c[] = {0.0, 1.0};

    return SimpleButcherTable< double >( 2, 2, A, b, c );
  }


  // TVD3ButcherTable
  // ----------------

  SimpleButcherTable< double > tvd3ButcherTable ()
  {
    static const double A[] = {0.0,  0.0,  0.0,
                               1.0,  0.0,  0.0,
                               0.25, 0.25, 0.0};
    static const double b[] = {1./6., 1./6., 2./3.};
    static const double c[] = {0.0,   1.0,   0.5};

    return SimpleButcherTable< double >( 3, 3, A, b, c );
  }


  // RK4bButcherTable
  // ----------------

  SimpleButcherTable< double > rk4ButcherTable ()
  {
    //class ExplicitRK4b
    static const double RK4_A[] =
      {0.0, 0.0, 0.0, 0.0,
       0.5, 0.0, 0.0, 0.0,
       0.0, 0.5, 0.0, 0.0,
       0.0, 0.0, 1.0, 0.0
      };
    static const double RK4_b[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    static const double RK4_c[] = {0.0, 0.5, 0.5, 1.0};

    return SimpleButcherTable< double >( 4, 4, RK4_A, RK4_b, RK4_c );
  }

  // Explicit6ButcherTable
  // ----------------

  SimpleButcherTable< double > expl6ButcherTable ()
  {
    //class ExplicitButcher6
    static const double Butcher6_A[] =
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       2.0/9.0, 4.0/9.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       7.0/36.0, 2.0/9.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,
       -35.0/144.0, -55.0/36.0, 35.0/48.0, 15.0/8.0, 0.0, 0.0, 0.0,
       -1.0/360.0, -11.0/36.0, -1.0/8.0, 1.0/2.0, 1.0/10.0, 0.0, 0.0,
       -41.0/260.0, 22.0/13.0, 43.0/156.0, -118.0/39.0, 32.0/195.0, 80.0/39.0, 0.0
      };
    static const double Butcher6_b[] =
      {13.0/200.0, 0.0, 11.0/40.0, 11.0/40.0, 4.0/25.0, 4.0/25.0, 13.0/200.0};
    static const double Butcher6_c[] =
      {0.0, 0.5, 2.0/3.0, 1.0/3.0, 5.0/6.0, 1.0/6.0, 1.0};

    return SimpleButcherTable< double >( 7, 6, Butcher6_A, Butcher6_b, Butcher6_c );
  }




  //////////////////////////////////////////////////////
  //
  //  Implicit Butcher Tables
  //
  //////////////////////////////////////////////////////

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

} // namespace DuneODE

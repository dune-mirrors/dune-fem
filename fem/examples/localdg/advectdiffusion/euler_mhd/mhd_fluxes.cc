#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <string>

// #include "mhd_eqns.hh"

namespace Mhd {
  /***************************************************************************
     new structures
  ****************************************************************************/

  //  typedef MhdSolver::VEC1D VEC;
  //  typedef MhdSolver::MAT MAT;

  /***************************************************************************
     global variables (initialized in initmhd8.c, mhd8.c)
  ****************************************************************************/

  const int ewave=MhdSolver::ewave;
  //EOS #define gamma MhdSolver::GAMMA
  // #define maxtime MhdSolver::maxtime
  // #define initialize MhdSolver::initialize
  // #define test_problem MhdSolver::test_problem
#define _MHD8

  /***************************************************************************
     global control variables
  ****************************************************************************/

#define docheck MhdSolver::docheck
#define tstep_auto MhdSolver::tstep_auto
#define use_hll_fix MhdSolver::use_hll_fix
#define verbose_mode MhdSolver::verbose_mode
#define use_point_fix MhdSolver::use_point_fix 
#define use_entropy_fix MhdSolver::use_entropy_fix
#define use_velocity_fix MhdSolver::use_velocity_fix
#define check_rmv_errors MhdSolver::check_rmv_errors
#define use_reference_fix MhdSolver::use_reference_fix
#define check_2D_conservation MhdSolver::check_2D_conservation
#define hll_use_roe_mean_values MhdSolver::hll_use_roe_mean_values
#define visc MhdSolver::visc
#define algorithm MhdSolver::algorithm
   
  /***************************************************************************
     global counters 
  ****************************************************************************/

#define undo MhdSolver::undo
#define hllcount MhdSolver::hllcount
#define point_fixes MhdSolver::point_fixes
#define entropy_fixes MhdSolver::entropy_fixes
#define velocity_fixes MhdSolver::velocity_fixes
#define reference_fixes MhdSolver::reference_fixes
#define conservation_errors_2D MhdSolver::conservation_errors_2D
#define max_check_rmv_error MhdSolver::max_check_rmv_error
#define max_conservation_error_2D MhdSolver::max_conservation_error_2D

  /***************************************************************************
     unused global variables 
  ****************************************************************************/

#define filenametime "/dev/null"
#define filenamegamma "/dev/null"

  /***************************************************************************
     functions
  ****************************************************************************/

  /* vecmat.c */

  double vmult(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2);
  void vadd(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2, MhdSolver::VEC1D *pret);
  void vsub(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2, MhdSolver::VEC1D *pret);
  void vcopy(const MhdSolver::VEC1D pfrom, MhdSolver::VEC1D *pto);
  void svmult(const double ps, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret);
  void svdiv(const double ps, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret);
  void madd(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret);
  void msub(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret);
  void mcopy(const MhdSolver::MAT pfrom, MhdSolver::MAT * pto);
  void smmult(const double ps, const MhdSolver::MAT pm, MhdSolver::MAT *pret);
  void smdiv(const double ps, const MhdSolver::MAT pm, MhdSolver::MAT *pret);
  void mvmult(const MhdSolver::MAT pmat, const MhdSolver::VEC1D pvec, MhdSolver::VEC1D *pret);
  void vmmult(const MhdSolver::VEC1D pvec, const MhdSolver::MAT pmat, MhdSolver::VEC1D *pret);
  void mmult(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret);

  /* init*.c (e.g. initmhd.c, initeuler.c) */

  void choose_test_problem(int pselect);

  /* bct.c */

  void flux_bct(const MhdSolver::VEC1D pul,
                const MhdSolver::VEC1D pur,
                MhdSolver::VEC1D *pret,
                double gamma);
  void initialize_bct(void);
  void statistics_bct(void);

  /* dw.c */

  void flux_dw(const MhdSolver::VEC1D pul,
               const MhdSolver::VEC1D pur,
               MhdSolver::VEC1D *pret,
               double gamma);
  void initialize_dw(void);
  void statistics_dw(void);

  /* hlle.c */

  void flux_hlle(const MhdSolver::VEC1D pul,
                 const MhdSolver::VEC1D pur,
                 MhdSolver::VEC1D *pret,
                 double gamma);
  void initialize_hlle(void);
  void statistics_hlle(void);

  /* hllemg.c */

  void flux_hllemg(const MhdSolver::VEC1D pul,
                   const MhdSolver::VEC1D pur,
                   MhdSolver::VEC1D *pret,
                   double gamma);
  void initialize_hllemg(void);
  void statistics_hllemg(void);

  /* hlleml.c */

  void flux_hlleml(const MhdSolver::VEC1D pul,
                   const MhdSolver::VEC1D pur,
                   MhdSolver::VEC1D *pret,
                   double gamma);
  void initialize_hlleml(void);
  void statistics_hlleml(void);

  /* roe.c */

  void flux_roe(const MhdSolver::VEC1D pul,
                const MhdSolver::VEC1D pur,
                MhdSolver::VEC1D *pret,
                double gamma);
  void initialize_roe(void);
  void statistics_roe(void);
  /***************************************************************************
     global definitions
  ****************************************************************************/

  const double eps=1e-12;
  const double del=5e-12;

  const int dim=8;

  /* pi */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


  /***************************************************************************
     simple vector routines
  ****************************************************************************/

  double vmult(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2)
  {
    int li;
    double lerg = 0;
  
    for (li=0;li<dim;li++) lerg += pu1[li] * pu2[li];

    return lerg;
  }

  void vadd(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2, MhdSolver::VEC1D *pret)
  {
    int li;
  
    for (li=0;li<dim;li++) (*pret)[li] = pu1[li] + pu2[li];
  }

  void vsub(const MhdSolver::VEC1D pu1, const MhdSolver::VEC1D pu2, MhdSolver::VEC1D *pret)
  {
    int li;
  
    for (li=0;li<dim;li++) (*pret)[li] = pu1[li] - pu2[li];
  }

  void vcopy(const MhdSolver::VEC1D pfrom, MhdSolver::VEC1D *pto)
  {
    int li;

    for (li=0;li<dim;li++) (*pto)[li] = pfrom[li];
  }

  void svmult(const double ps, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret)
  {
    int li;
  
    for (li=0;li<dim;li++) (*pret)[li] = ps * pu[li];
  }

  void svdiv(const double ps, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret)
  {
    int li;

    for (li=0;li<dim;li++) (*pret)[li] = pu[li] / ps;
  }

  /***************************************************************************
     simple matrix routines
  ****************************************************************************/

  void madd(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret)
  {
    int li,lj;
  
    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        (*pret)[li][lj] = pmat1[li][lj] + pmat2[li][lj];
  }

  void msub(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret)
  {
    int li,lj;
  
    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        (*pret)[li][lj] = pmat1[li][lj] - pmat2[li][lj];
  }

  void mcopy(const MhdSolver::MAT pfrom, MhdSolver::MAT *pto)
  {
    int li,lj;

    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        (*pto)[li][lj] = pfrom[li][lj];
  }

  void smmult(const double ps, const MhdSolver::MAT pm, MhdSolver::MAT *pret)
  {
    int li,lj;

    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        (*pret)[li][lj] = pm[li][lj] * ps;
  }

  void smdiv(const double ps, const MhdSolver::MAT pm, MhdSolver::MAT *pret)
  {
    int li,lj;

    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        (*pret)[li][lj] = pm[li][lj] / ps;
  }

  /***************************************************************************
     matrix - vector routines
  ****************************************************************************/

  void mvmult(const MhdSolver::MAT pmat, const MhdSolver::VEC1D pvec, MhdSolver::VEC1D *pret)
  {
    int li,lj;

    for (li=0;li<dim;li++)
      {
        (*pret)[li] = 0.0;
        for (lj=0;lj<dim;lj++) (*pret)[li] += pmat[li][lj] * pvec[lj];
      }
  }

  void vmmult(const MhdSolver::VEC1D pvec, const MhdSolver::MAT pmat, MhdSolver::VEC1D *pret)
  {
    int li,lj;

    for (lj=0;lj<dim;lj++)
      {
        (*pret)[lj] = 0.0;
        for (li=0;li<dim;li++) (*pret)[lj] += pvec[li] * pmat[li][lj];
      }
  }

  /***************************************************************************
     matrix - matrix routines
  ****************************************************************************/

  void mmult(const MhdSolver::MAT pmat1, const MhdSolver::MAT pmat2, MhdSolver::MAT *pret)
  {
    int li,lj,lk;

    for (li=0;li<dim;li++)
      for (lj=0;lj<dim;lj++)
        {
          (*pret)[li][lj] = 0.0;
          for (lk=0;lk<dim;lk++) (*pret)[li][lj] += pmat1[li][lk] * pmat2[lk][lj];
        }
  }


  /***************************************************************************
     switch functions for numerical fluxes
  ****************************************************************************/

  double phi(double px)
  {
    double lret;

    if (px < eps)
      {
        lret = 0.0;
      }
    else if (px < del)
      {
        lret = (px - eps) / (del - eps);
      }
    else
      {
        lret = 1.0;
      }

    return lret;
  }

  double psi(double px)
  {
    double lret;

    if (px < -del)
      {
        lret = 1.0;
      }
    else if (px < -eps)
      {
        lret = 0.5 * (px + del) / (eps - del) + 1.0; 
      }
    else if (px < eps)
      {
        lret = 0.5;
      }
    else if (px < del)
      {
        lret = 0.5 * (px - del) / (eps - del);
      }
    else
      {
        lret = 0.0;
      }

    return lret;
  }

  /***************************************************************************
     functions depending on underlying equations                           
  ****************************************************************************/

  /* flux function */

  void f(const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,lBx,lBy,lBz,lrhoE,lp;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    lux = lrhoux / lrho;
    luy = lrhouy / lrho;
    luz = lrhouz / lrho;
    lp  = (gamma - 1.0) 
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));

    (*pret)[0] = lrhoux;
    (*pret)[1] =   lrhoux * lux
      + lp
      + ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI );
    (*pret)[2] = lrhoux * luy - ( lBx * lBy ) / ( 4.0 * M_PI );
    (*pret)[3] = lrhoux * luz - ( lBx * lBz ) / ( 4.0 * M_PI );
    (*pret)[4] = 0.0;
    (*pret)[5] = lux * lBy - luy * lBx;
    (*pret)[6] = lux * lBz - luz * lBx;
    (*pret)[7] = (lrhoE + lp + ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI )) * lux
      - (lBx / ( 4.0 * M_PI )) * ( luy * lBy + luz * lBz );
  }

  /* Jacobian of the flux function */

  void Df(const MhdSolver::VEC1D pu, MhdSolver::MAT *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,lBx,lBy,lBz,lrhoE,lp,lpifac;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    lux = lrhoux / lrho;
    luy = lrhouy / lrho;
    luz = lrhouz / lrho;
    lp  = (gamma - 1.0) 
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));

    lpifac = 1.0 / (4.0 * M_PI);

    (*pret)[0][0] = 0.0;
    (*pret)[0][1] = 1.0;
    (*pret)[0][2] = 0.0;
    (*pret)[0][3] = 0.0;
    (*pret)[0][4] = 0.0;
    (*pret)[0][5] = 0.0;
    (*pret)[0][6] = 0.0;
    (*pret)[0][7] = 0.0;

    (*pret)[1][0] =   0.5 * (gamma - 3.0) * lux * lux
      + 0.5 * (gamma - 1.0) * (luy * luy + luz * luz);
    (*pret)[1][1] = (3.0 - gamma) * lux;
    (*pret)[1][2] = (1.0 - gamma) * luy;
    (*pret)[1][3] = (1.0 - gamma) * luz;
    (*pret)[1][4] = 0.0;
    (*pret)[1][5] = (2.0 - gamma) * lpifac * lBy;
    (*pret)[1][6] = (2.0 - gamma) * lpifac * lBz;
    (*pret)[1][7] = gamma - 1.0;

    (*pret)[2][0] = - lux * luy;
    (*pret)[2][1] = luy;
    (*pret)[2][2] = lux;
    (*pret)[2][3] = 0.0;
    (*pret)[2][4] = 0.0;
    (*pret)[2][5] = - lpifac * lBx;
    (*pret)[2][6] = 0.0;
    (*pret)[2][7] = 0.0;

    (*pret)[3][0] = - lux * luz;
    (*pret)[3][1] = luz;
    (*pret)[3][2] = 0.0;
    (*pret)[3][3] = lux;
    (*pret)[3][4] = 0.0;
    (*pret)[3][5] = 0.0;
    (*pret)[3][6] = - lpifac * lBx;
    (*pret)[3][7] = 0.0;

    (*pret)[4][0] = 0.0;
    (*pret)[4][1] = 0.0;
    (*pret)[4][2] = 0.0;
    (*pret)[4][3] = 0.0;
    (*pret)[4][4] = 0.0;
    (*pret)[4][5] = 0.0;
    (*pret)[4][6] = 0.0;
    (*pret)[4][7] = 0.0;

    (*pret)[5][0] = (luy * lBx - lux * lBy) / lrho;
    (*pret)[5][1] = lBy / lrho;
    (*pret)[5][2] = - lBx / lrho;
    (*pret)[5][3] = 0.0;
    (*pret)[5][4] = 0.0;
    (*pret)[5][5] = lux;
    (*pret)[5][6] = 0.0;
    (*pret)[5][7] = 0.0;

    (*pret)[6][0] = (luz * lBx - lux * lBz) / lrho;
    (*pret)[6][1] = lBz / lrho;
    (*pret)[6][2] = 0.0;
    (*pret)[6][3] = - lBx / lrho;
    (*pret)[6][4] = 0.0;
    (*pret)[6][5] = 0.0;
    (*pret)[6][6] = lux;
    (*pret)[6][7] = 0.0;

    (*pret)[7][0] =   lux * (   0.5 * (gamma - 1.0) 
                                * (lux * lux + luy * luy + luz * luz) 
                                - (   lrhoE + lp 
                                      + 0.5 * lpifac * (lBy * lBy + lBz * lBz)) / lrho)
      + lpifac * lBx * (luy * lBy + luz * lBz) / lrho;
    (*pret)[7][1] =   (lrhoE + lp + 0.5 * lpifac * (lBy * lBy + lBz * lBz)) / lrho
      + (1.0 - gamma) * lux * lux;
    (*pret)[7][2] = (1.0 - gamma) * lux * luy - lpifac * lBx * lBy / lrho;
    (*pret)[7][3] = (1.0 - gamma) * lux * luz - lpifac * lBx * lBz / lrho;
    (*pret)[7][4] = 0.0;
    (*pret)[7][5] = lpifac * ((2.0 - gamma) * lux * lBy - luy * lBx);
    (*pret)[7][6] = lpifac * ((2.0 - gamma) * lux * lBz - luz * lBx);
    (*pret)[7][7] = gamma * lux;
  }

  /* calculate single eigenvalue */

  void lambda(const int pidx, const MhdSolver::VEC1D pu, double *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,lBx,lBy,lBz;
    double lrhoE,lp,lcs2,lvax2,lva2,lrad1,lrad2;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    lux = lrhoux / lrho;
    luy = lrhouy / lrho;
    luz = lrhouz / lrho;

    switch (pidx)
      {
      case 0 :
        lp    = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2  = gamma * lp / lrho;
        lva2  = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );

        lrad1 = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad1) <= eps) lrad1 = 0.0;
        lrad2 = 0.5 * ( ( lva2 + lcs2 ) + sqrt( lrad1 ) );
        if (fabs(lrad2) <= eps) lrad2 = 0.0;
        (*pret) = lux - sqrt ( lrad2 );
        break;
      case 1 :
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        (*pret) = lux - sqrt(lvax2);
        break;
      case 2 :
        lp    = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2  = gamma * lp / lrho;
        lva2  = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );

        lrad1 = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad1) <= eps) lrad1 = 0.0;
        lrad2 = 0.5 * ( ( lva2 + lcs2 ) - sqrt( lrad1 ) );
        if (fabs(lrad2) <= eps) lrad2 = 0.0;
        (*pret) = lux - sqrt ( lrad2 );
        break;
      case 3 :
        (*pret) = lux;
        break;
      case 4 :
        (*pret) = 0.0;
        break;
      case 5 :
        lp    = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2  = gamma * lp / lrho;
        lva2  = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );

        lrad1 = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad1) <= eps) lrad1 = 0.0;
        lrad2 = 0.5 * ( ( lva2 + lcs2 ) - sqrt( lrad1 ) );
        if (fabs(lrad2) <= eps) lrad2 = 0.0;
        (*pret) = lux + sqrt ( lrad2 );
        break;
      case 6 :
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        (*pret) = lux + sqrt(lvax2);
        break;
      case 7 :
        lp    = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2  = gamma * lp / lrho;
        lva2  = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );

        lrad1 = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad1) <= eps) lrad1 = 0.0;
        lrad2 = 0.5 * ( ( lva2 + lcs2 ) + sqrt( lrad1 ) );
        if (fabs(lrad2) <= eps) lrad2 = 0.0;
        (*pret) = lux + sqrt ( lrad2 );
        break;
      default :
        printf("Invalid argument \"pidx == %d\" in function lambda!\n\n",pidx);
        exit(17);
      }
  }

  /* calculate vector of eigenvalues */

  void lambda_vec(const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,lBx,lBy,lBz,lrhoE,lp;
    double lrad1,lrad2,lvf,lvax,lvs,lcs2,lvax2,lva2;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    lux = lrhoux / lrho;
    luy = lrhouy / lrho;
    luz = lrhouz / lrho;

    lp    = (gamma - 1.0)
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
    lvax2 = lBx * lBx / ( 4.0 * M_PI * lrho );
    lcs2  = gamma * lp / lrho;
    lva2  = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );

    lrad1 = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
    if (fabs(lrad1) <= eps) lrad1 = 0.0;
    lrad2 = 0.5 * ( ( lva2 + lcs2 ) + sqrt( lrad1 ) );
    if (fabs(lrad2) <= eps) lrad2 = 0.0;
    lvf   = sqrt( lrad2 );
    lrad2 = 0.5 * ( ( lva2 + lcs2 ) - sqrt( lrad1 ) );
    if (fabs(lrad2) <= eps) lrad2 = 0.0;
    lvs   = sqrt( lrad2 );
    if (fabs( lva2 ) <= eps) lvax2 = 0.0;
    lvax  = sqrt( lvax2 );

    assert(lvf > 0.0);

    (*pret)[0] = lux - lvf;
    (*pret)[1] = lux - lvax;
    (*pret)[2] = lux - lvs;
    (*pret)[3] = lux;
    (*pret)[4] = 0.0;
    (*pret)[5] = lux + lvs;
    (*pret)[6] = lux + lvax;
    (*pret)[7] = lux + lvf;
  }

  /* calculate single right eigenvector */

  void r(const int pidx, const int psign, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp,lsgnBx;
    double lcs2,lcs,lvax2,lvax,lva2,lvf2,lvf,lvs2,lvs;
    double lalphas,lalphaf,lbetay,lbetaz;
    double lsqrtBy2pBz2,lsqrtvf2mvs2,lsqrtvfs,lrad;
    double lmultvs,lmultvax,lmultvf;
    int li;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    switch (pidx)
      {
      case 0 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);

        lmultvf  = lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                  + lalphas * lalphas * (lvf2 + lvax2) );

        (*pret)[0] = -lrho * lalphaf;
        (*pret)[1] = -lrhoux * lalphaf + lrho * lalphaf * lvf;
        (*pret)[2] = -lrhouy * lalphaf - lrho * lalphas * lbetay * lvax * lsgnBx;
        (*pret)[3] = -lrhouz * lalphaf - lrho * lalphas * lbetaz * lvax * lsgnBx;
        (*pret)[4] = 0.0;
        (*pret)[5] = -lalphas * lbetay * lvf * sqrt( 4.0 * M_PI * lrho );
        (*pret)[6] = -lalphas * lbetaz * lvf * sqrt( 4.0 * M_PI * lrho );
        (*pret)[7] = - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
          * lalphaf
          + lrhoux * lalphaf * lvf
          - lrhouy * lalphas * lbetay * lvax * lsgnBx
          - lrhouz * lalphas * lbetaz * lvax * lsgnBx
          - lBy * lalphas * lbetay * lvf * sqrt( lrho / (4.0 * M_PI) )
          - lBz * lalphas * lbetaz * lvf * sqrt( lrho / (4.0 * M_PI) )
          - (lalphaf * gamma * lp) / (gamma - 1.0);

        for (li=0;li<8;li++) (*pret)[li] *= lmultvf;

        break;
      case 1 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);

        lmultvax = lvf / sqrt(2.0);

        (*pret)[0] = 0.0;
        (*pret)[1] = 0.0;
        (*pret)[2] = -lrho * lbetaz;
        (*pret)[3] =  lrho * lbetay;
        (*pret)[4] = 0.0;
        (*pret)[5] = -lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetaz;
        (*pret)[6] =  lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetay;
        (*pret)[7] = - lrhouy * lbetaz
          + lrhouz * lbetay;

        for (li=0;li<8;li++) (*pret)[li] *= lmultvax;

        break;
      case 2 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert( lvf > 0.0 );
        if (fabs(lvs2) <= eps) lvs = 0.0; else lvs  = sqrt( lvs2 );
        if (fabs(lcs2) <= eps) lcs = 0.0; else lcs  = sqrt( lcs2 );

        lmultvs  = lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                  + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

        (*pret)[0] = -lrho * lalphas;
        (*pret)[1] = -lrhoux * lalphas + lrho * lalphas * lvs;
        (*pret)[2] = -lrhouy * lalphas + lrho * lalphaf * lbetay * lcs * lsgnBx;
        (*pret)[3] = -lrhouz * lalphas + lrho * lalphaf * lbetaz * lcs * lsgnBx;
        (*pret)[4] = 0.0;
        (*pret)[5] =  lalphaf * lbetay * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
        (*pret)[6] =  lalphaf * lbetaz * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
        (*pret)[7] = - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
          * lalphas
          + lrhoux * lalphas * lvs
          + lrhouy * lalphaf * lbetay * lcs * lsgnBx
          + lrhouz * lalphaf * lbetaz * lcs * lsgnBx
          + lBy * lalphaf * lbetay * lcs2 
          * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
          + lBz * lalphaf * lbetaz * lcs2
          * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
          - (lalphas * gamma * lp) / (gamma - 1.0);

        for (li=0;li<8;li++) (*pret)[li] *= lmultvs;

        break;
      case 3 :
        (*pret)[0] = -lrho;
        (*pret)[1] = -lrhoux;
        (*pret)[2] = -lrhouy;
        (*pret)[3] = -lrhouz;
        (*pret)[4] = 0.0;
        (*pret)[5] = 0.0;
        (*pret)[6] = 0.0;
        (*pret)[7] = -0.5 * (   lrhoux * lrhoux + lrhouy * lrhouy 
                                + lrhouz * lrhouz ) / lrho;
        break;
      case 4:
        (*pret)[0] = 0.0;
        (*pret)[1] = 0.0;
        (*pret)[2] = 0.0;
        (*pret)[3] = 0.0;
        (*pret)[4] = 1.0;
        (*pret)[5] = 0.0;
        (*pret)[6] = 0.0;
        (*pret)[7] = 0.0;
      case 5 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lvs2) <= eps) lvs = 0.0; else lvs  = sqrt( lvs2 );
        if (fabs(lcs2) <= eps) lcs = 0.0; else lcs  = sqrt( lcs2 );

        lmultvs  = lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                  + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

        (*pret)[0] =  lrho * lalphas;
        (*pret)[1] =  lrhoux * lalphas + lrho * lalphas * lvs;
        (*pret)[2] =  lrhouy * lalphas + lrho * lalphaf * lbetay * lcs * lsgnBx;
        (*pret)[3] =  lrhouz * lalphas + lrho * lalphaf * lbetaz * lcs * lsgnBx;
        (*pret)[4] =  0.0;
        (*pret)[5] = -lalphaf * lbetay * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
        (*pret)[6] = -lalphaf * lbetaz * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
        (*pret)[7] =   0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
          * lalphas
          + lrhoux * lalphas * lvs
          + lrhouy * lalphaf * lbetay * lcs * lsgnBx
          + lrhouz * lalphaf * lbetaz * lcs * lsgnBx
          - lBy * lalphaf * lbetay * lcs2 
          * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
          - lBz * lalphaf * lbetaz * lcs2
          * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
          + (lalphas * gamma * lp) / (gamma - 1.0);

        for (li=0;li<8;li++) (*pret)[li] *= lmultvs;

        break;
      case 6 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);

        lmultvax = lvf / sqrt(2.0);

        (*pret)[0] = 0.0;
        (*pret)[1] = 0.0;
        (*pret)[2] =  lrho * lbetaz;
        (*pret)[3] = -lrho * lbetay;
        (*pret)[4] = 0.0;
        (*pret)[5] = -lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetaz;
        (*pret)[6] =  lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetay;
        (*pret)[7] =   lrhouy * lbetaz
          - lrhouz * lbetay;

        for (li=0;li<8;li++) (*pret)[li] *= lmultvax;

        break;
      case 7 :
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);

        lmultvf  = lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                  + lalphas * lalphas * (lvf2 + lvax2) );

        (*pret)[0] =  lrho * lalphaf;
        (*pret)[1] =  lrhoux * lalphaf + lrho * lalphaf * lvf;
        (*pret)[2] =  lrhouy * lalphaf - lrho * lalphas * lbetay * lvax * lsgnBx;
        (*pret)[3] =  lrhouz * lalphaf - lrho * lalphas * lbetaz * lvax * lsgnBx;
        (*pret)[4] =  0.0;
        (*pret)[5] =  lalphas * lbetay * lvf * sqrt( 4.0 * M_PI * lrho );
        (*pret)[6] =  lalphas * lbetaz * lvf * sqrt( 4.0 * M_PI * lrho );
        (*pret)[7] =   0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
          * lalphaf
          + lrhoux * lalphaf * lvf
          - lrhouy * lalphas * lbetay * lvax * lsgnBx
          - lrhouz * lalphas * lbetaz * lvax * lsgnBx
          + lBy * lalphas * lbetay * lvf * sqrt( lrho / (4.0 * M_PI) )
          + lBz * lalphas * lbetaz * lvf * sqrt( lrho / (4.0 * M_PI) )
          + (lalphaf * gamma * lp) / (gamma - 1.0);

        for (li=0;li<8;li++) (*pret)[li] *= lmultvf;

        break;
      default :
        printf("Invalid argument \"pidx == %d\" in function r!\n\n",pidx);
        exit(17);
      }
  }

  /* calculate vector of right eigenvectors */

  void r_vec(const int psign, const MhdSolver::VEC1D pu, MhdSolver::VEC1D pret[dim], double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp,lsgnBx;
    double lcs2,lcs,lvax2,lvax,lva2,lvf2,lvf,lvs2,lvs;
    double lalphas,lalphaf,lbetay,lbetaz;
    double lsqrtBy2pBz2,lsqrtvf2mvs2,lsqrtvfs,lrad;
    double lmultvs,lmultvax,lmultvf;
    int li;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
    lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
    if (lsqrtBy2pBz2 > eps)
      {
        lbetay = lBy / lsqrtBy2pBz2;
        lbetaz = lBz / lsqrtBy2pBz2; 
      }
    else lbetay = lbetaz = 1.0 / sqrt(2.0);
    lp   = (gamma - 1.0) 
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
    lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
    lcs2     = gamma * lp / lrho;
    lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
    lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
    if (fabs(lrad) <= eps) lrad = 0.0;
    lsqrtvfs = sqrt( lrad );
    lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
    lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
    lrad     = lvf2 - lvs2;
    if (lrad <= eps) lrad = 0.0;
    lsqrtvf2mvs2 = sqrt( lrad );
    if (lsqrtvf2mvs2 > eps)
      {
        lrad    = lvf2 - lvax2;
        if (lrad > eps)
          lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
        else lalphaf = 0.0;
        lrad    = lvf2 - lcs2;
        if (lrad > eps)
          lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
        else lalphas = 0.0;
      }
    else if (psign > 0)
      {
        lalphas = 1.0;
        lalphaf = 0.0;
      }
    else
      {
        lalphas = 0.0;
        lalphaf = 1.0;
      }
    lvf  = sqrt( lvf2 );
    assert(lvf > 0.0);
    if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
    if (fabs(lvs2)  <= eps) lvs  = 0.0; else lvs  = sqrt( lvs2 );
    if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );

    lmultvs  = lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                              + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
    lmultvax = lvf / sqrt(2.0);
    lmultvf  = lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                              + lalphas * lalphas * (lvf2 + lvax2) );

    pret[0][0] = -lrho * lalphaf;
    pret[0][1] = -lrhoux * lalphaf + lrho * lalphaf * lvf;
    pret[0][2] = -lrhouy * lalphaf - lrho * lalphas * lbetay * lvax * lsgnBx;
    pret[0][3] = -lrhouz * lalphaf - lrho * lalphas * lbetaz * lvax * lsgnBx;
    pret[0][4] = 0.0;
    pret[0][5] = -lalphas * lbetay * lvf * sqrt( 4.0 * M_PI * lrho );
    pret[0][6] = -lalphas * lbetaz * lvf * sqrt( 4.0 * M_PI * lrho );
    pret[0][7] = - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
      * lalphaf
      + lrhoux * lalphaf * lvf
      - lrhouy * lalphas * lbetay * lvax * lsgnBx
      - lrhouz * lalphas * lbetaz * lvax * lsgnBx
      - lBy * lalphas * lbetay * lvf * sqrt( lrho / (4.0 * M_PI) )
      - lBz * lalphas * lbetaz * lvf * sqrt( lrho / (4.0 * M_PI) )
      - (lalphaf * gamma * lp) / (gamma - 1.0);

    pret[1][0] = 0.0;
    pret[1][1] = 0.0;
    pret[1][2] = -lrho * lbetaz;
    pret[1][3] =  lrho * lbetay;
    pret[1][4] = 0.0;
    pret[1][5] = -lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetaz;
    pret[1][6] =  lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetay;
    pret[1][7] = - lrhouy * lbetaz
      + lrhouz * lbetay;

    pret[2][0] = -lrho * lalphas;
    pret[2][1] = -lrhoux * lalphas + lrho * lalphas * lvs;
    pret[2][2] = -lrhouy * lalphas + lrho * lalphaf * lbetay * lcs * lsgnBx;
    pret[2][3] = -lrhouz * lalphas + lrho * lalphaf * lbetaz * lcs * lsgnBx;
    pret[2][4] = 0.0;
    pret[2][5] =  lalphaf * lbetay * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
    pret[2][6] =  lalphaf * lbetaz * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
    pret[2][7] = - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
      * lalphas
      + lrhoux * lalphas * lvs
      + lrhouy * lalphaf * lbetay * lcs * lsgnBx
      + lrhouz * lalphaf * lbetaz * lcs * lsgnBx
      + lBy * lalphaf * lbetay * lcs2 
      * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
      + lBz * lalphaf * lbetaz * lcs2
      * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
      - (lalphas * gamma * lp) / (gamma - 1.0);

    pret[3][0] = -lrho;
    pret[3][1] = -lrhoux;
    pret[3][2] = -lrhouy;
    pret[3][3] = -lrhouz;
    pret[3][4] = 0.0;
    pret[3][5] = 0.0;
    pret[3][6] = 0.0;
    pret[3][7] = -0.5 * (   lrhoux * lrhoux + lrhouy * lrhouy 
                            + lrhouz * lrhouz ) / lrho;

    pret[4][0] = 0.0;
    pret[4][1] = 0.0;
    pret[4][2] = 0.0;
    pret[4][3] = 0.0;
    pret[4][4] = 1.0;
    pret[4][5] = 0.0;
    pret[4][6] = 0.0;
    pret[4][7] = 0.0;

    pret[5][0] =  lrho * lalphas;
    pret[5][1] =  lrhoux * lalphas + lrho * lalphas * lvs;
    pret[5][2] =  lrhouy * lalphas + lrho * lalphaf * lbetay * lcs * lsgnBx;
    pret[5][3] =  lrhouz * lalphas + lrho * lalphaf * lbetaz * lcs * lsgnBx;
    pret[5][4] = 0.0;
    pret[5][5] = -lalphaf * lbetay * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
    pret[5][6] = -lalphaf * lbetaz * lcs2 * sqrt( lrho * 4.0 * M_PI ) / lvf;
    pret[5][7] =   0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
      * lalphas
      + lrhoux * lalphas * lvs
      + lrhouy * lalphaf * lbetay * lcs * lsgnBx
      + lrhouz * lalphaf * lbetaz * lcs * lsgnBx
      - lBy * lalphaf * lbetay * lcs2 
      * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
      - lBz * lalphaf * lbetaz * lcs2
      * sqrt ( lrho / ( 4.0 * M_PI ) ) / lvf
      + (lalphas * gamma * lp) / (gamma - 1.0);

    pret[6][0] = 0.0;
    pret[6][1] = 0.0;
    pret[6][2] =  lrho * lbetaz;
    pret[6][3] = -lrho * lbetay;
    pret[6][4] = 0.0;
    pret[6][5] = -lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetaz;
    pret[6][6] =  lsgnBx * sqrt( 4.0 * M_PI * lrho ) * lbetay;
    pret[6][7] =   lrhouy * lbetaz
      - lrhouz * lbetay;

    pret[7][0] =  lrho * lalphaf;
    pret[7][1] =  lrhoux * lalphaf + lrho * lalphaf * lvf;
    pret[7][2] =  lrhouy * lalphaf - lrho * lalphas * lbetay * lvax * lsgnBx;
    pret[7][3] =  lrhouz * lalphaf - lrho * lalphas * lbetaz * lvax * lsgnBx;
    pret[7][4] =  0.0;
    pret[7][5] =  lalphas * lbetay * lvf * sqrt( 4.0 * M_PI * lrho );
    pret[7][6] =  lalphas * lbetaz * lvf * sqrt( 4.0 * M_PI * lrho );
    pret[7][7] =   0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
      * lalphaf
      + lrhoux * lalphaf * lvf
      - lrhouy * lalphas * lbetay * lvax * lsgnBx
      - lrhouz * lalphas * lbetaz * lvax * lsgnBx
      + lBy * lalphas * lbetay * lvf * sqrt( lrho / (4.0 * M_PI) )
      + lBz * lalphas * lbetaz * lvf * sqrt( lrho / (4.0 * M_PI) )
      + (lalphaf * gamma * lp) / (gamma - 1.0);

    for (li=0;li<8;li++)
      {
        pret[0][li] *= lmultvf;
        pret[1][li] *= lmultvax;
        pret[2][li] *= lmultvs;
        pret[5][li] *= lmultvs;
        pret[6][li] *= lmultvax;
        pret[7][li] *= lmultvf;
      }
  }

  /* calculate single left eigenvector */

  void l(const int pidx, const int psign, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp,lsgnBx;
    double lcs2,lcs,lvax2,lvax,lva2,lvf2,lvf,lvs2,lvs;
    double lalphas,lalphaf,lbetay,lbetaz;
    double lsqrtBy2pBz2,lsqrtvf2mvs2,lsqrtvfs,lrad;
    double lmultvs,lmultvax,lmultvf,lsqrt4pirho;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    switch (pidx)
      {
      case 0:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvf  = 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                    + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

        (*pret)[0] = lmultvf  * ( - lalphaf * lvf * lux / lrho
                                  + lalphas * lbetay * lvax * lsgnBx * luy / lrho
                                  + lalphas * lbetaz * lvax * lsgnBx * luz / lrho
                                  + lalphaf / lrho * 0.5 * (1.0 - gamma) 
                                  * (lux * lux + luy * luy + luz * luz));
        (*pret)[1] = lmultvf  * (   lalphaf * lvf / lrho
                                    + lalphaf / lrho * (gamma - 1.0) * lux );
        (*pret)[2] = lmultvf  * ( - lalphas * lbetay * lvax * lsgnBx / lrho
                                  + lalphaf / lrho * (gamma - 1.0) * luy );
        (*pret)[3] = lmultvf  * ( - lalphas * lbetaz * lvax * lsgnBx / lrho
                                  + lalphaf / lrho * (gamma - 1.0) * luz );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvf  * ( - lalphas * lbetay * lvf / lsqrt4pirho
                                  + lalphaf / lrho * (gamma - 1.0) * lBy
                                  / (4.0 * M_PI) );
        (*pret)[6] = lmultvf  * ( - lalphas * lbetaz * lvf / lsqrt4pirho
                                  + lalphaf / lrho * (gamma - 1.0) * lBz
                                  / (4.0 * M_PI) );
        (*pret)[7] = lmultvf  * (   lalphaf / lrho * (1.0 - gamma) );
        break;
      case 1:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvax = 1.0 / ( sqrt(2.0) * lvf );

        (*pret)[0] = lmultvax * (   lbetaz * luy / lrho
                                    - lbetay * luz / lrho );
        (*pret)[1] = 0.0;
        (*pret)[2] = lmultvax * ( - lbetaz / lrho );
        (*pret)[3] = lmultvax * (   lbetay / lrho );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvax * ( - lsgnBx * lbetaz / lsqrt4pirho );
        (*pret)[6] = lmultvax * (   lsgnBx * lbetay / lsqrt4pirho );
        (*pret)[7] = 0.0;
        break;
      case 2:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lvs2)  <= eps) lvs  = 0.0; else lvs  = sqrt( lvs2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvs  = 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                 + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

        (*pret)[0] = lmultvs  * ( - lalphas * lvs * lux / lrho
                                  - lalphaf * lbetay * lcs * lsgnBx * luy / lrho
                                  - lalphaf * lbetaz * lcs * lsgnBx * luz / lrho
                                  + lalphas / lrho * (1.0 - gamma) * 0.5
                                  * (lux * lux + luy * luy + luz * luz));
        (*pret)[1] = lmultvs  * (   lalphas * lvs / lrho
                                    + lalphas / lrho * (gamma - 1.0) * lux );
        (*pret)[2] = lmultvs  * (   lalphaf * lbetay * lcs * lsgnBx / lrho
                                    + lalphas / lrho * (gamma - 1.0) * luy );
        (*pret)[3] = lmultvs  * (   lalphaf * lbetaz * lcs * lsgnBx / lrho
                                    + lalphas / lrho * (gamma - 1.0) * luz );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvs  * (   lalphaf * lbetay * lcs2 / (lvf * lsqrt4pirho)
                                    + lalphas / lrho * (gamma - 1.0) * lBy
                                    / (4.0 * M_PI) );
        (*pret)[6] = lmultvs  * (   lalphaf * lbetaz * lcs2 / (lvf * lsqrt4pirho)
                                    + lalphas / lrho * (gamma - 1.0) * lBz
                                    / (4.0 * M_PI) );
        (*pret)[7] = lmultvs  * (   lalphas / lrho * (1.0 - gamma) );
        break;
      case 3:
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));

        (*pret)[0] = - 1.0 / lrho
          + (gamma - 1.0) / (gamma * lp) * 0.5 
          * (lux * lux + luy * luy + luz * luz);
        (*pret)[1] =   (1.0 - gamma) / (gamma * lp) * lux;
        (*pret)[2] =   (1.0 - gamma) / (gamma * lp) * luy;
        (*pret)[3] =   (1.0 - gamma) / (gamma * lp) * luz;
        (*pret)[4] =   0.0;
        (*pret)[5] =   (1.0 - gamma) / (gamma * lp) * lBy / (4.0 * M_PI);
        (*pret)[6] =   (1.0 - gamma) / (gamma * lp) * lBz / (4.0 * M_PI);
        (*pret)[7] =   (gamma - 1.0) / (gamma * lp);
        break;
      case 4:
        (*pret)[0] = 0.0;
        (*pret)[1] = 0.0;
        (*pret)[2] = 0.0;
        (*pret)[3] = 0.0;
        (*pret)[4] = 1.0;
        (*pret)[5] = 0.0;
        (*pret)[6] = 0.0;
        (*pret)[7] = 0.0;
        break;
      case 5:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lvs2)  <= eps) lvs  = 0.0; else lvs  = sqrt( lvs2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvs  = 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                 + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

        (*pret)[0] = lmultvs  * ( - lalphas * lvs * lux / lrho
                                  - lalphaf * lbetay * lcs * lsgnBx * luy / lrho
                                  - lalphaf * lbetaz * lcs * lsgnBx * luz / lrho
                                  + lalphas / lrho * (gamma - 1.0) * 0.5
                                  * (lux * lux + luy * luy + luz * luz));
        (*pret)[1] = lmultvs  * (   lalphas * lvs / lrho
                                    + lalphas / lrho * (1.0 - gamma) * lux );
        (*pret)[2] = lmultvs  * (   lalphaf * lbetay * lcs * lsgnBx / lrho
                                    + lalphas / lrho * (1.0 - gamma) * luy );
        (*pret)[3] = lmultvs  * (   lalphaf * lbetaz * lcs * lsgnBx / lrho
                                    + lalphas / lrho * (1.0 - gamma) * luz );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvs  * ( - lalphaf * lbetay * lcs2 / (lvf * lsqrt4pirho)
                                  + lalphas / lrho * (1.0 - gamma) * lBy
                                  / (4.0 * M_PI) );
        (*pret)[6] = lmultvs  * ( - lalphaf * lbetaz * lcs2 / (lvf * lsqrt4pirho)
                                  + lalphas / lrho * (1.0 - gamma) * lBz
                                  / (4.0 * M_PI) );
        (*pret)[7] = lmultvs  * (   lalphas / lrho * (gamma - 1.0) );
        break;
      case 6:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvax = 1.0 / ( sqrt(2.0) * lvf );

        (*pret)[0] = lmultvax * ( - lbetaz * luy / lrho
                                  + lbetay * luz / lrho );
        (*pret)[1] = 0.0;
        (*pret)[2] = lmultvax * (   lbetaz / lrho );
        (*pret)[3] = lmultvax * ( - lbetay / lrho );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvax * ( - lsgnBx * lbetaz / lsqrt4pirho );
        (*pret)[6] = lmultvax * (   lsgnBx * lbetay / lsqrt4pirho );
        (*pret)[7] = 0.0;
        break;
      case 7:
        if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
        lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
        if (lsqrtBy2pBz2 > eps)
          {
            lbetay = lBy / lsqrtBy2pBz2;
            lbetaz = lBz / lsqrtBy2pBz2; 
          }
        else lbetay = lbetaz = 1.0 / sqrt(2.0);
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        lrad     = lvf2 - lvs2;
        if (lrad <= eps) lrad = 0.0;
        lsqrtvf2mvs2 = sqrt( lrad );
        if (lsqrtvf2mvs2 > eps)
          {
            lrad    = lvf2 - lvax2;
            if (lrad > eps)
              lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphaf = 0.0;
            lrad    = lvf2 - lcs2;
            if (lrad > eps)
              lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
            else lalphas = 0.0;
          }
        else if (psign > 0)
          {
            lalphas = 1.0;
            lalphaf = 0.0;
          }
        else
          {
            lalphas = 0.0;
            lalphaf = 1.0;
          }
        lvf  = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

        lmultvf  = 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                    + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

        (*pret)[0] = lmultvf  * ( - lalphaf * lvf * lux / lrho
                                  + lalphas * lbetay * lvax * lsgnBx * luy / lrho
                                  + lalphas * lbetaz * lvax * lsgnBx * luz / lrho
                                  + lalphaf / lrho * 0.5 * (gamma - 1.0) 
                                  * (lux * lux + luy * luy + luz * luz));
        (*pret)[1] = lmultvf  * (   lalphaf * lvf / lrho
                                    + lalphaf / lrho * (1.0 - gamma) * lux );
        (*pret)[2] = lmultvf  * ( - lalphas * lbetay * lvax * lsgnBx / lrho
                                  + lalphaf / lrho * (1.0 - gamma) * luy );
        (*pret)[3] = lmultvf  * ( - lalphas * lbetaz * lvax * lsgnBx / lrho
                                  + lalphaf / lrho * (1.0 - gamma) * luz );
        (*pret)[4] = 0.0;
        (*pret)[5] = lmultvf  * (   lalphas * lbetay * lvf / lsqrt4pirho
                                    + lalphaf / lrho * (1.0 - gamma) * lBy
                                    / (4.0 * M_PI) );
        (*pret)[6] = lmultvf  * (   lalphas * lbetaz * lvf / lsqrt4pirho
                                    + lalphaf / lrho * (1.0 - gamma) * lBz
                                    / (4.0 * M_PI) );
        (*pret)[7] = lmultvf  * (   lalphaf / lrho * (gamma - 1.0) );
        break;
      default:
        abort();
      }
  }

  /* calculate matrix of left eigenvectors */

  void l_mat(const int psign, const MhdSolver::VEC1D pu, MhdSolver::MAT *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp,lsgnBx;
    double lcs2,lcs,lvax2,lvax,lva2,lvf2,lvf,lvs2,lvs;
    double lalphas,lalphaf,lbetay,lbetaz;
    double lsqrtBy2pBz2,lsqrtvf2mvs2,lsqrtvfs,lrad;
    double lmultvs,lmultvax,lmultvf,lsqrt4pirho;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    if (lBx >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;
    lsqrtBy2pBz2 = sqrt( lBy * lBy + lBz * lBz );
    if (lsqrtBy2pBz2 > eps)
      {
        lbetay = lBy / lsqrtBy2pBz2;
        lbetaz = lBz / lsqrtBy2pBz2; 
      }
    else lbetay = lbetaz = 1.0 / sqrt(2.0);
    lp   = (gamma - 1.0) 
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
    lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
    lcs2     = gamma * lp / lrho;
    lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
    lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
    if (fabs(lrad) <= eps) lrad = 0.0;
    lsqrtvfs = sqrt( lrad );
    lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
    lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
    lrad     = lvf2 - lvs2;
    if (lrad <= eps) lrad = 0.0;
    lsqrtvf2mvs2 = sqrt( lrad );
    if (lsqrtvf2mvs2 > eps)
      {
        lrad    = lvf2 - lvax2;
        if (lrad > eps)
          lalphaf = sqrt( lrad ) / lsqrtvf2mvs2;
        else lalphaf = 0.0;
        lrad    = lvf2 - lcs2;
        if (lrad > eps)
          lalphas = sqrt( lrad ) / lsqrtvf2mvs2;
        else lalphas = 0.0;
      }
    else if (psign > 0)
      {
        lalphas = 1.0;
        lalphaf = 0.0;
      }
    else
      {
        lalphas = 0.0;
        lalphaf = 1.0;
      }
    lvf  = sqrt( lvf2 );
    assert(lvf > 0.0);
    if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
    if (fabs(lvs2)  <= eps) lvs  = 0.0; else lvs  = sqrt( lvs2 );
    if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
    lsqrt4pirho = sqrt( 4.0 * M_PI * lrho );

    lmultvs  = 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                             + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
    lmultvax = 1.0 / ( sqrt(2.0) * lvf );
    lmultvf  = 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

    (*pret)[0][0] = lmultvf  * ( - lalphaf * lvf * lux / lrho
                                 + lalphas * lbetay * lvax * lsgnBx * luy / lrho
                                 + lalphas * lbetaz * lvax * lsgnBx * luz / lrho
                                 + lalphaf / lrho * 0.5 * (1.0 - gamma) 
                                 * (lux * lux + luy * luy + luz * luz));
    (*pret)[0][1] = lmultvf  * (   lalphaf * lvf / lrho
                                   + lalphaf / lrho * (gamma - 1.0) * lux );
    (*pret)[0][2] = lmultvf  * ( - lalphas * lbetay * lvax * lsgnBx / lrho
                                 + lalphaf / lrho * (gamma - 1.0) * luy );
    (*pret)[0][3] = lmultvf  * ( - lalphas * lbetaz * lvax * lsgnBx / lrho
                                 + lalphaf / lrho * (gamma - 1.0) * luz );
    (*pret)[0][4] = 0.0;
    (*pret)[0][5] = lmultvf  * ( - lalphas * lbetay * lvf / lsqrt4pirho
                                 + lalphaf / lrho * (gamma - 1.0) * lBy
                                 / (4.0 * M_PI) );
    (*pret)[0][6] = lmultvf  * ( - lalphas * lbetaz * lvf / lsqrt4pirho
                                 + lalphaf / lrho * (gamma - 1.0) * lBz
                                 / (4.0 * M_PI) );
    (*pret)[0][7] = lmultvf  * (   lalphaf / lrho * (1.0 - gamma) );

    (*pret)[1][0] = lmultvax * (   lbetaz * luy / lrho
                                   - lbetay * luz / lrho );
    (*pret)[1][1] = 0.0;
    (*pret)[1][2] = lmultvax * ( - lbetaz / lrho );
    (*pret)[1][3] = lmultvax * (   lbetay / lrho );
    (*pret)[1][4] = 0.0;
    (*pret)[1][5] = lmultvax * ( - lsgnBx * lbetaz / lsqrt4pirho );
    (*pret)[1][6] = lmultvax * (   lsgnBx * lbetay / lsqrt4pirho );
    (*pret)[1][7] = 0.0;

    (*pret)[2][0] = lmultvs  * ( - lalphas * lvs * lux / lrho
                                 - lalphaf * lbetay * lcs * lsgnBx * luy / lrho
                                 - lalphaf * lbetaz * lcs * lsgnBx * luz / lrho
                                 + lalphas / lrho * (1.0 - gamma) * 0.5
                                 * (lux * lux + luy * luy + luz * luz));
    (*pret)[2][1] = lmultvs  * (   lalphas * lvs / lrho
                                   + lalphas / lrho * (gamma - 1.0) * lux );
    (*pret)[2][2] = lmultvs  * (   lalphaf * lbetay * lcs * lsgnBx / lrho
                                   + lalphas / lrho * (gamma - 1.0) * luy );
    (*pret)[2][3] = lmultvs  * (   lalphaf * lbetaz * lcs * lsgnBx / lrho
                                   + lalphas / lrho * (gamma - 1.0) * luz );
    (*pret)[2][4] = 0.0;
    (*pret)[2][5] = lmultvs  * (   lalphaf * lbetay * lcs2 / (lvf * lsqrt4pirho)
                                   + lalphas / lrho * (gamma - 1.0) * lBy
                                   / (4.0 * M_PI) );
    (*pret)[2][6] = lmultvs  * (   lalphaf * lbetaz * lcs2 / (lvf * lsqrt4pirho)
                                   + lalphas / lrho * (gamma - 1.0) * lBz
                                   / (4.0 * M_PI) );
    (*pret)[2][7] = lmultvs  * (   lalphas / lrho * (1.0 - gamma) );

    (*pret)[3][0] = - 1.0 / lrho
      + (gamma - 1.0) / (gamma * lp) * 0.5 
      * (lux * lux + luy * luy + luz * luz);
    (*pret)[3][1] =   (1.0 - gamma) / (gamma * lp) * lux;
    (*pret)[3][2] =   (1.0 - gamma) / (gamma * lp) * luy;
    (*pret)[3][3] =   (1.0 - gamma) / (gamma * lp) * luz;
    (*pret)[3][4] =   0.0;
    (*pret)[3][5] =   (1.0 - gamma) / (gamma * lp) * lBy / (4.0 * M_PI);
    (*pret)[3][6] =   (1.0 - gamma) / (gamma * lp) * lBz / (4.0 * M_PI);
    (*pret)[3][7] =   (gamma - 1.0) / (gamma * lp);

    (*pret)[4][0] = 0.0;
    (*pret)[4][1] = 0.0;
    (*pret)[4][2] = 0.0;
    (*pret)[4][3] = 0.0;
    (*pret)[4][4] = 1.0;
    (*pret)[4][5] = 0.0;
    (*pret)[4][6] = 0.0;
    (*pret)[4][7] = 0.0;

    (*pret)[5][0] = lmultvs  * ( - lalphas * lvs * lux / lrho
                                 - lalphaf * lbetay * lcs * lsgnBx * luy / lrho
                                 - lalphaf * lbetaz * lcs * lsgnBx * luz / lrho
                                 + lalphas / lrho * (gamma - 1.0) * 0.5
                                 * (lux * lux + luy * luy + luz * luz));
    (*pret)[5][1] = lmultvs  * (   lalphas * lvs / lrho
                                   + lalphas / lrho * (1.0 - gamma) * lux );
    (*pret)[5][2] = lmultvs  * (   lalphaf * lbetay * lcs * lsgnBx / lrho
                                   + lalphas / lrho * (1.0 - gamma) * luy );
    (*pret)[5][3] = lmultvs  * (   lalphaf * lbetaz * lcs * lsgnBx / lrho
                                   + lalphas / lrho * (1.0 - gamma) * luz );
    (*pret)[5][4] = 0.0;
    (*pret)[5][5] = lmultvs  * ( - lalphaf * lbetay * lcs2 / (lvf * lsqrt4pirho)
                                 + lalphas / lrho * (1.0 - gamma) * lBy
                                 / (4.0 * M_PI) );
    (*pret)[5][6] = lmultvs  * ( - lalphaf * lbetaz * lcs2 / (lvf * lsqrt4pirho)
                                 + lalphas / lrho * (1.0 - gamma) * lBz
                                 / (4.0 * M_PI) );
    (*pret)[5][7] = lmultvs  * (   lalphas / lrho * (gamma - 1.0) );

    (*pret)[6][0] = lmultvax * ( - lbetaz * luy / lrho
                                 + lbetay * luz / lrho );
    (*pret)[6][1] = 0.0;
    (*pret)[6][2] = lmultvax * (   lbetaz / lrho );
    (*pret)[6][3] = lmultvax * ( - lbetay / lrho );
    (*pret)[6][4] = 0.0;
    (*pret)[6][5] = lmultvax * ( - lsgnBx * lbetaz / lsqrt4pirho );
    (*pret)[6][6] = lmultvax * (   lsgnBx * lbetay / lsqrt4pirho );
    (*pret)[6][7] = 0.0;

    (*pret)[7][0] = lmultvf  * ( - lalphaf * lvf * lux / lrho
                                 + lalphas * lbetay * lvax * lsgnBx * luy / lrho
                                 + lalphas * lbetaz * lvax * lsgnBx * luz / lrho
                                 + lalphaf / lrho * 0.5 * (gamma - 1.0) 
                                 * (lux * lux + luy * luy + luz * luz));
    (*pret)[7][1] = lmultvf  * (   lalphaf * lvf / lrho
                                   + lalphaf / lrho * (1.0 - gamma) * lux );
    (*pret)[7][2] = lmultvf  * ( - lalphas * lbetay * lvax * lsgnBx / lrho
                                 + lalphaf / lrho * (1.0 - gamma) * luy );
    (*pret)[7][3] = lmultvf  * ( - lalphas * lbetaz * lvax * lsgnBx / lrho
                                 + lalphaf / lrho * (1.0 - gamma) * luz );
    (*pret)[7][4] = 0.0;
    (*pret)[7][5] = lmultvf  * (   lalphas * lbetay * lvf / lsqrt4pirho
                                   + lalphaf / lrho * (1.0 - gamma) * lBy
                                   / (4.0 * M_PI) );
    (*pret)[7][6] = lmultvf  * (   lalphas * lbetaz * lvf / lsqrt4pirho
                                   + lalphaf / lrho * (1.0 - gamma) * lBz
                                   / (4.0 * M_PI) );
    (*pret)[7][7] = lmultvf  * (   lalphaf / lrho * (gamma - 1.0) );
  }

  /* calculate gradient of a single eigenvalue */

  void grad_lambda(int pidx, int psign, const MhdSolver::VEC1D pu, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp;
    double lcs2,lcs,lvax2,lvax,lva2,lva,lvf2,lvf,lvs2,lvs=0.0,lsqrtvfs;
    double lmgva2,lmgvax2,lmgcs2,lrad;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    (*pret)[0] = -lux * ltau;
    (*pret)[1] = ltau;
    (*pret)[2] = 0.0;
    (*pret)[3] = 0.0;
    (*pret)[4] = 0.0;
    (*pret)[5] = 0.0;
    (*pret)[6] = 0.0;
    (*pret)[7] = 0.0;

    switch (pidx)
      {
      case 0 :
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        if (fabs(lva2)  <= eps) lva  = 0.0; else lva  = sqrt( lva2 );
        lvf      = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lrad) > eps)
          {
            lmgva2   = 1.0 + ( lva2 + lcs2 ) / lsqrtvfs;
            lmgvax2  = -2.0 * lcs2 / lsqrtvfs;
            lmgcs2   = 1.0 + ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

            (*pret)[0] -= (  lmgva2  * (- lva2  * ltau) 
                             + lmgvax2 * (- lvax2 * ltau)
                             + lmgcs2  * ltau 
                             * ((- lcs2 
                                 + gamma * (gamma - 1.0) * 0.5
                                 * (  lux * lux
                                      + luy * luy
                                      + luz * luz))))
              * 0.25 / lvf;
            (*pret)[1] -= lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
              * 0.25 / lvf ;
            (*pret)[2] -= lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
              * 0.25 / lvf;
            (*pret)[3] -= lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
              * 0.25 / lvf;
            (*pret)[4] -= 0.0;
            (*pret)[5] -= (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                             + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                         / (4.0 * M_PI)))
              * 0.25 / lvf;
            (*pret)[6] -= (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                             + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                         / (4.0 * M_PI)))
              * 0.25 / lvf;
            (*pret)[7] -= lmgcs2 * (gamma * (gamma - 1.0) * ltau)
              * 0.25 / lvf;
          }
        else if (psign > 0)
          {
            (*pret)[0] -= -0.5 * lvax * ltau;
            (*pret)[1] -= 0.0;
            (*pret)[2] -= 0.0;
            (*pret)[3] -= 0.0;
            (*pret)[4] -= 0.0;
            (*pret)[5] -= 0.0;
            (*pret)[6] -= 0.0;
            (*pret)[7] -= 0.0;
          }
        else
          {
            (*pret)[0] -= ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvf;
            (*pret)[1] -= (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvf;
            (*pret)[2] -= (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvf;
            (*pret)[3] -= (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvf;
            (*pret)[4] -= 0.0;
            (*pret)[5] -= (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvf;
            (*pret)[6] -= (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvf;
            (*pret)[7] -= (gamma * (gamma - 1.0) * ltau) * 0.25 / lvf;
          }
        break; 
      case 1 :
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        (*pret)[0] -= -0.5 * lvax * ltau;
        (*pret)[1] -= 0.0;
        (*pret)[2] -= 0.0;
        (*pret)[3] -= 0.0;
        (*pret)[4] -= 0.0;
        (*pret)[5] -= 0.0;
        (*pret)[6] -= 0.0;
        (*pret)[7] -= 0.0;
        break; 
      case 2 :
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        if (fabs(lva2)  <= eps) lva  = 0.0; else lva  = sqrt( lva2 );
        lvf      = sqrt( lvf2 );
        assert(lvf > 0.0);

        if (fabs(lrad) > eps)
          {
            if (fabs(lvs2) > eps)
              {
                lvs      = sqrt( lvs2 );
                lmgva2   = 1.0 - ( lva2 + lcs2 ) / lsqrtvfs;
                lmgvax2  = 2.0 * lcs2 / lsqrtvfs;
                lmgcs2   = 1.0 - ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

                (*pret)[0] -= (  lmgva2  * (- lva2  * ltau) 
                                 + lmgvax2 * (- lvax2 * ltau)
                                 + lmgcs2  * ltau 
                                 * ((- lcs2 
                                     + gamma * (gamma - 1.0) * 0.5
                                     * (  lux * lux
                                          + luy * luy
                                          + luz * luz))))
                  * 0.25 / lvs;
                (*pret)[1] -= lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
                  * 0.25 / lvs ;
                (*pret)[2] -= lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
                  * 0.25 / lvs;
                (*pret)[3] -= lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
                  * 0.25 / lvs;
                (*pret)[4] -= 0.0;
                (*pret)[5] -= (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                                 + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                             / (4.0 * M_PI)))
                  * 0.25 / lvs;
                (*pret)[6] -= (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                                 + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                             / (4.0 * M_PI)))
                  * 0.25 / lvs;
                (*pret)[7] -= lmgcs2 * (gamma * (gamma - 1.0) * ltau)
                  * 0.25 / lvs;
              }
          }
        else if (psign > 0)
          {
            (*pret)[0] -= ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvs;
            (*pret)[1] -= (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvs;
            (*pret)[2] -= (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvs;
            (*pret)[3] -= (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvs;
            (*pret)[4] -= 0.0;
            (*pret)[5] -= (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            (*pret)[6] -= (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            (*pret)[7] -= (gamma * (gamma - 1.0) * ltau) * 0.25 / lvs;
          }
        else
          {
            (*pret)[0] -= -0.5 * lvax * ltau;
            (*pret)[1] -= 0.0;
            (*pret)[2] -= 0.0;
            (*pret)[3] -= 0.0;
            (*pret)[4] -= 0.0;
            (*pret)[5] -= 0.0;
            (*pret)[6] -= 0.0;
            (*pret)[7] -= 0.0;
          }
        break; 
      case 3:
        break;
      case 4:
        (*pret)[0] = 0.0;
        (*pret)[1] = 0.0;
        break;
      case 5 :
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        if (fabs(lva2)  <= eps) lva  = 0.0; else lva  = sqrt( lva2 );
        lvf      = sqrt( lvf2 );
        assert(lvf > 0.0);
   
        if (fabs(lrad) > eps)
          {
            if (fabs(lvs2) > eps)
              {
                lvs      = sqrt( lvs2 );
                lmgva2   = 1.0 - ( lva2 + lcs2 ) / lsqrtvfs;
                lmgvax2  = 2.0 * lcs2 / lsqrtvfs;
                lmgcs2   = 1.0 - ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

                (*pret)[0] += (  lmgva2  * (- lva2  * ltau) 
                                 + lmgvax2 * (- lvax2 * ltau)
                                 + lmgcs2  * ltau 
                                 * ((- lcs2 
                                     + gamma * (gamma - 1.0) * 0.5
                                     * (  lux * lux
                                          + luy * luy
                                          + luz * luz))))
                  * 0.25 / lvs;
                (*pret)[1] += lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
                  * 0.25 / lvs ;
                (*pret)[2] += lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
                  * 0.25 / lvs;
                (*pret)[3] += lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
                  * 0.25 / lvs;
                (*pret)[4] += 0.0;
                (*pret)[5] += (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                                 + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                             / (4.0 * M_PI)))
                  * 0.25 / lvs;
                (*pret)[6] += (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                                 + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                             / (4.0 * M_PI)))
                  * 0.25 / lvs;
                (*pret)[7] += lmgcs2 * (gamma * (gamma - 1.0) * ltau)
                  * 0.25 / lvs;
              }
          }
        else if (psign > 0)
          {
            (*pret)[0] += ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvs;
            (*pret)[1] += (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvs;
            (*pret)[2] += (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvs;
            (*pret)[3] += (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvs;
            (*pret)[4] += 0.0;
            (*pret)[5] += (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            (*pret)[6] += (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            (*pret)[7] += (gamma * (gamma - 1.0) * ltau) * 0.25 / lvs;
          }
        else
          {
            (*pret)[0] += -0.5 * lvax * ltau;
            (*pret)[1] += 0.0;
            (*pret)[2] += 0.0;
            (*pret)[3] += 0.0;
            (*pret)[4] += 0.0;
            (*pret)[5] += 0.0;
            (*pret)[6] += 0.0;
            (*pret)[7] += 0.0;
          }
        break; 
      case 6 :
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        (*pret)[0] += -0.5 * lvax * ltau;
        (*pret)[1] += 0.0;
        (*pret)[2] += 0.0;
        (*pret)[3] += 0.0;
        (*pret)[4] += 0.0;
        (*pret)[5] += 0.0;
        (*pret)[6] += 0.0;
        (*pret)[7] += 0.0;
        break; 
      case 7 :
        lp   = (gamma - 1.0) 
          * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
             - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
        lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
        lcs2     = gamma * lp / lrho;
        lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
        lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
        if (fabs(lrad) <= eps) lrad = 0.0;
        lsqrtvfs = sqrt( lrad );
        lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
        if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
        if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
        if (fabs(lva2)  <= eps) lva  = 0.0; else lva  = sqrt( lva2 );
        lvf      = sqrt( lvf2 );
        assert(lvf > 0.0);
        if (fabs(lrad) > eps)
          {
            lmgva2   = 1.0 + ( lva2 + lcs2 ) / lsqrtvfs;
            lmgvax2  = -2.0 * lcs2 / lsqrtvfs;
            lmgcs2   = 1.0 + ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

            (*pret)[0] += (  lmgva2  * (- lva2  * ltau) 
                             + lmgvax2 * (- lvax2 * ltau)
                             + lmgcs2  * ltau 
                             * ((- lcs2 
                                 + gamma * (gamma - 1.0) * 0.5
                                 * (  lux * lux
                                      + luy * luy
                                      + luz * luz))))
              * 0.25 / lvf;
            (*pret)[1] += lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
              * 0.25 / lvf ;
            (*pret)[2] += lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
              * 0.25 / lvf;
            (*pret)[3] += lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
              * 0.25 / lvf;
            (*pret)[4] += 0.0;
            (*pret)[5] += (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                             + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                         / (4.0 * M_PI)))
              * 0.25 / lvf;
            (*pret)[6] += (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                             + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                         / (4.0 * M_PI)))
              * 0.25 / lvf;
            (*pret)[7] += lmgcs2 * (gamma * (gamma - 1.0) * ltau)
              * 0.25 / lvf;
          }
        else if (psign > 0)
          {
            (*pret)[0] += -0.5 * lvax * ltau;
            (*pret)[1] += 0.0;
            (*pret)[2] += 0.0;
            (*pret)[3] += 0.0;
            (*pret)[4] += 0.0;
            (*pret)[5] += 0.0;
            (*pret)[6] += 0.0;
            (*pret)[7] += 0.0;
          }
        else
          {
            (*pret)[0] += ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvf;
            (*pret)[1] += (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvf;
            (*pret)[2] += (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvf;
            (*pret)[3] += (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvf;
            (*pret)[4] += 0.0;
            (*pret)[5] += (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvf;
            (*pret)[6] += (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvf;
            (*pret)[7] += (gamma * (gamma - 1.0) * ltau) * 0.25 / lvf;
          }
        break; 
      default :
        printf("Invalid argument \"pidx == %d\" in function grad_lambda!\n\n",
               pidx);
        exit(17);
      }
  }

  /* calculate gradients of all eigenvalues */

  void grad_lambda_vec(const int psign, const MhdSolver::VEC1D pu, MhdSolver::VEC1D pret[dim], double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lux,luy,luz,ltau,lBx,lBy,lBz,lrhoE,lp;
    double lcs2,lcs,lvax2,lvax,lva2,lva,lvf2,lvf,lvs2,lvs,lsqrtvfs;
    double lmgva2,lmgvax2,lmgcs2,lrad;
    int li;

    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];

    assert(lrho > 0.0);

    ltau   = 1.0 / lrho;
    lux    = ltau * lrhoux;
    luy    = ltau * lrhouy;
    luz    = ltau * lrhouz;

    for (li=0;li<dim;li++)
      {
        pret[li][0] = -lux * ltau;
        pret[li][1] = ltau;
        pret[li][2] = 0.0;
        pret[li][3] = 0.0;
        pret[li][4] = 0.0;
        pret[li][5] = 0.0;
        pret[li][6] = 0.0;
        pret[li][7] = 0.0;
      }

    pret[4][0] = 0.0;
    pret[4][1] = 0.0;

    lp   = (gamma - 1.0) 
      * (lrhoE - 0.5 * ( lrhoux * lux + lrhouy * luy + lrhouz * luz )
         - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
    lvax2    = lBx * lBx / ( 4.0 * M_PI * lrho );
    lcs2     = gamma * lp / lrho;
    lva2     = ( lBx * lBx + lBy * lBy + lBz * lBz ) / ( 4.0 * M_PI * lrho );
    lrad     = (lva2 + lcs2) * (lva2 + lcs2) - 4.0 * lvax2 * lcs2;
    if (fabs(lrad) <= eps) lrad = 0.0;
    lsqrtvfs = sqrt( lrad );
    lvf2     = 0.5 * ( (lva2 + lcs2) + lsqrtvfs );
    lvs2     = 0.5 * ( (lva2 + lcs2) - lsqrtvfs );
    if (fabs(lvax2) <= eps) lvax = 0.0; else lvax = sqrt( lvax2 );
    if (fabs(lcs2)  <= eps) lcs  = 0.0; else lcs  = sqrt( lcs2 );
    if (fabs(lva2)  <= eps) lva  = 0.0; else lva  = sqrt( lva2 );
    if (fabs(lvs2)  <= eps) lvs  = 0.0; else lvs  = sqrt( lvs2 );
    lvf      = sqrt( lvf2 );

    assert(lvf > 0.0);

    if (fabs(lrad) > eps)
      {
        lmgva2   = 1.0 + ( lva2 + lcs2 ) / lsqrtvfs;
        lmgvax2  = -2.0 * lcs2 / lsqrtvfs;
        lmgcs2   = 1.0 + ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

        pret[0][0] -= (  lmgva2  * (- lva2  * ltau) 
                         + lmgvax2 * (- lvax2 * ltau)
                         + lmgcs2  * ltau 
                         * ((- lcs2 
                             + gamma * (gamma - 1.0) * 0.5
                             * (  lux * lux
                                  + luy * luy
                                  + luz * luz))))
          * 0.25 / lvf;
        pret[0][1] -=  lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
          * 0.25 / lvf ;
        pret[0][2] -=  lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
          * 0.25 / lvf;
        pret[0][3] -=  lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
          * 0.25 / lvf;
        pret[0][4] -= 0.0;
        pret[0][5] -= (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvf;
        pret[0][6] -= (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvf;
        pret[0][7] -= lmgcs2 * (gamma * (gamma - 1.0) * ltau)
          * 0.25 / lvf;
      }
    else if (psign > 0)
      {
        pret[0][0] -=-0.5 * lvax * ltau;
        pret[0][1] -= 0.0;
        pret[0][2] -= 0.0;
        pret[0][3] -= 0.0;
        pret[0][4] -= 0.0;
        pret[0][5] -= 0.0;
        pret[0][6] -= 0.0; 
        pret[0][7] -= 0.0;
      }
    else
      {
        pret[0][0] -= ltau * (- lcs2 
                              + gamma * (gamma - 1.0) * 0.5
                              * (  lux * lux
                                   + luy * luy
                                   + luz * luz))
          * 0.25 / lvf;
        pret[0][1] -= (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvf;
        pret[0][2] -= (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvf;
        pret[0][3] -= (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvf;
        pret[0][4] -= 0.0;
        pret[0][5] -= (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
          * 0.25 / lvf;
        pret[0][6] -= (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
          * 0.25 / lvf;
        pret[0][7] -= (gamma * (gamma - 1.0) * ltau) * 0.25 / lvf;
      }

    pret[1][0] -=-0.5 * lvax * ltau;
    pret[1][1] -= 0.0;
    pret[1][2] -= 0.0;
    pret[1][3] -= 0.0;
    pret[1][4] -= 0.0;
    pret[1][5] -= 0.0;
    pret[1][6] -= 0.0;
    pret[1][7] -= 0.0;
  
    if ((fabs(lrad) > eps) && (fabs(lvs2) > eps))
      {
        lmgva2   = 1.0 - ( lva2 + lcs2 ) / lsqrtvfs;
        lmgvax2  = 2.0 * lcs2 / lsqrtvfs;
        lmgcs2   = 1.0 - ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

        pret[2][0] -= (  lmgva2  * (- lva2  * ltau) 
                         + lmgvax2 * (- lvax2 * ltau)
                         + lmgcs2  * ltau 
                         * ((- lcs2 
                             + gamma * (gamma - 1.0) * 0.5
                             * (  lux * lux
                                  + luy * luy
                                  + luz * luz))))
          * 0.25 / lvs;
        pret[2][1] -= lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
          * 0.25 / lvs ;
        pret[2][2] -= lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
          * 0.25 / lvs;
        pret[2][3] -= lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
          * 0.25 / lvs;
        pret[2][4] -= 0.0;
        pret[2][5] -= (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvs;
        pret[2][6] -= (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvs;
        pret[2][7] -= lmgcs2 * (gamma * (gamma - 1.0) * ltau)
          * 0.25 / lvs;

        lmgva2   = 1.0 - ( lva2 + lcs2 ) / lsqrtvfs;
        lmgvax2  = 2.0 * lcs2 / lsqrtvfs;
        lmgcs2   = 1.0 - ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

        pret[5][0] += (  lmgva2  * (- lva2  * ltau) 
                         + lmgvax2 * (- lvax2 * ltau)
                         + lmgcs2  * ltau 
                         * ((- lcs2 
                             + gamma * (gamma - 1.0) * 0.5
                             * (  lux * lux
                                  + luy * luy
                                  + luz * luz))))
          * 0.25 / lvs;
        pret[5][1] += lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
          * 0.25 / lvs ;
        pret[5][2] += lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
          * 0.25 / lvs;
        pret[5][3] += lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
          * 0.25 / lvs;
        pret[5][4] += 0.0;
        pret[5][5] += (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvs;
        pret[5][6] += (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvs;
        pret[5][7] += lmgcs2 * (gamma * (gamma - 1.0) * ltau)
          * 0.25 / lvs;
      }
    else if (fabs(lrad) <= eps)
      {
        if (psign > 0)
          {
            pret[2][0] -= ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvs;
            pret[2][1] -= (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvs;
            pret[2][2] -= (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvs;
            pret[2][3] -= (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvs;
            pret[2][4] -= 0.0;
            pret[2][5] -= (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            pret[2][6] -= (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            pret[2][7] -= (gamma * (gamma - 1.0) * ltau) * 0.25 / lvs;

            pret[5][0] += ltau * (- lcs2 
                                  + gamma * (gamma - 1.0) * 0.5
                                  * (  lux * lux
                                       + luy * luy
                                       + luz * luz))
              * 0.25 / lvs;
            pret[5][1] += (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvs;
            pret[5][2] += (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvs;
            pret[5][3] += (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvs;
            pret[5][4] += 0.0;
            pret[5][5] += (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            pret[5][6] += (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
              * 0.25 / lvs;
            pret[5][7] += (gamma * (gamma - 1.0) * ltau) * 0.25 / lvs;
          }
        else
          {
            pret[2][0] -= -0.5 * lvax * ltau;
            pret[2][1] -= 0.0;
            pret[2][2] -= 0.0;
            pret[2][3] -= 0.0;
            pret[2][4] -= 0.0;
            pret[2][5] -= 0.0;
            pret[2][6] -= 0.0;
            pret[2][7] -= 0.0;

            pret[5][0] += -0.5 * lvax * ltau;
            pret[5][1] += 0.0;
            pret[5][2] += 0.0;
            pret[5][3] += 0.0;
            pret[5][4] += 0.0;
            pret[5][5] += 0.0;
            pret[5][6] += 0.0;
            pret[5][7] += 0.0;
          }
      }

    pret[6][0] += -0.5 * lvax * ltau;
    pret[6][1] += 0.0;
    pret[6][2] += 0.0;
    pret[6][3] += 0.0;
    pret[6][4] += 0.0;
    pret[6][5] += 0.0;
    pret[6][6] += 0.0;
    pret[6][7] += 0.0;

    if (fabs(lrad) > eps)
      {
        lmgva2   = 1.0 + ( lva2 + lcs2 ) / lsqrtvfs;
        lmgvax2  = -2.0 * lcs2 / lsqrtvfs;
        lmgcs2   = 1.0 + ( lva2 + lcs2 - 2.0 * lvax2 ) / lsqrtvfs;

        pret[7][0] += (  lmgva2  * (- lva2  * ltau) 
                         + lmgvax2 * (- lvax2 * ltau)
                         + lmgcs2  * ltau 
                         * ((- lcs2 
                             + gamma * (gamma - 1.0) * 0.5
                             * (  lux * lux
                                  + luy * luy
                                  + luz * luz))))
          * 0.25 / lvf;
        pret[7][1] += lmgcs2 * (- gamma * (gamma - 1.0) * lux * ltau)
          * 0.25 / lvf ;
        pret[7][2] += lmgcs2 * (- gamma * (gamma - 1.0) * luy * ltau)
          * 0.25 / lvf;
        pret[7][3] += lmgcs2 * (- gamma * (gamma - 1.0) * luz * ltau)
          * 0.25 / lvf;
        pret[7][4] += 0.0;
        pret[7][5] += (  lmgva2 * lBy * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBy * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvf;
        pret[7][6] += (  lmgva2 * lBz * ltau / (2.0 * M_PI)
                         + lmgcs2 * (- gamma * (gamma - 1.0) * lBz * ltau
                                     / (4.0 * M_PI)))
          * 0.25 / lvf;
        pret[7][7] += lmgcs2 * (gamma * (gamma - 1.0) * ltau)
          * 0.25 / lvf;
      }
    else if (psign > 0)
      {
        pret[7][0] +=-0.5 * lvax * ltau;
        pret[7][1] += 0.0;
        pret[7][2] += 0.0;
        pret[7][3] += 0.0;
        pret[7][4] += 0.0;
        pret[7][5] += 0.0;
        pret[7][6] += 0.0; 
        pret[7][7] += 0.0;
      }
    else
      {
        pret[7][0] += ltau * (- lcs2 
                              + gamma * (gamma - 1.0) * 0.5
                              * (  lux * lux
                                   + luy * luy
                                   + luz * luz))
          * 0.25 / lvf;
        pret[7][1] += (- gamma * (gamma - 1.0) * lux * ltau) * 0.25 / lvf;
        pret[7][2] += (- gamma * (gamma - 1.0) * luy * ltau) * 0.25 / lvf;
        pret[7][3] += (- gamma * (gamma - 1.0) * luz * ltau) * 0.25 / lvf;
        pret[7][4] += 0.0;
        pret[7][5] += (- gamma * (gamma - 1.0) * lBy * ltau / (4.0 * M_PI))
          * 0.25 / lvf;
        pret[7][6] += (- gamma * (gamma - 1.0) * lBz * ltau / (4.0 * M_PI))
          * 0.25 / lvf;
        pret[7][7] += (gamma * (gamma - 1.0) * ltau) * 0.25 / lvf;
      }
  }

  /* Roe mean values (Brio - Wu; requires gamma = 2) */

  void roe_mean_values(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrhol,lrhouxl,lrhouyl,lrhouzl,luxl,luyl,luzl,lBxl,lByl,lBzl,lrhoEl;
    double lrhor,lrhouxr,lrhouyr,lrhouzr,luxr,luyr,luzr,lBxr,lByr,lBzr,lrhoEr;
    double lPl,lHl,lPr,lHr,lmvrho,lmvux,lmvuy,lmvuz,lmvBx,lmvBy,lmvBz,lmvH;
    double lsqrtrhol,lsqrtrhor,lmult;

    lrhol   = pul[0];
    lrhouxl = pul[1];
    lrhouyl = pul[2];
    lrhouzl = pul[3];
    lBxl    = pul[4];
    lByl    = pul[5];
    lBzl    = pul[6];
    lrhoEl  = pul[7];

    assert(lrhol > 0.0);

    luxl = lrhouxl / lrhol;
    luyl = lrhouyl / lrhol;
    luzl = lrhouzl / lrhol;
    lPl  = (gamma - 1.0) 
      * (lrhoEl - 0.5 * ( lrhouxl * luxl + lrhouyl * luyl + lrhouzl * luzl ))
      + (2.0 - gamma) * ( lByl * lByl + lBzl * lBzl ) / ( 8.0 * M_PI );
    lHl  = (lrhoEl + lPl) / lrhol;

    lrhor   = pur[0];
    lrhouxr = pur[1];
    lrhouyr = pur[2];
    lrhouzr = pur[3];
    lBxr    = pur[4];
    lByr    = pur[5];
    lBzr    = pur[6];
    lrhoEr  = pur[7];

    assert(lrhor > 0.0);

    luxr = lrhouxr / lrhor;
    luyr = lrhouyr / lrhor;
    luzr = lrhouzr / lrhor;
    lPr  = (gamma - 1.0) 
      * (lrhoEr - 0.5 * ( lrhouxr * luxr + lrhouyr * luyr + lrhouzr * luzr ))
      + (2.0 - gamma) * ( lByr * lByr + lBzr * lBzr ) / ( 8.0 * M_PI );
    lHr  = (lrhoEr + lPr) / lrhor;

    /* calculate Roe mean values */

    lsqrtrhol = sqrt(lrhol);
    lsqrtrhor = sqrt(lrhor);
    lmult     = 1.0 / (lsqrtrhol + lsqrtrhor);

    lmvrho = sqrt(lrhol * lrhor);
    lmvux  = (lsqrtrhol*luxl + lsqrtrhor*luxr) * lmult;
    lmvuy  = (lsqrtrhol*luyl + lsqrtrhor*luyr) * lmult;
    lmvuz  = (lsqrtrhol*luzl + lsqrtrhor*luzr) * lmult;
    lmvBx  = (lBxl/lsqrtrhol + lBxr/lsqrtrhor) * lmult * lmvrho;
    lmvBy  = (lByl/lsqrtrhol + lByr/lsqrtrhor) * lmult * lmvrho;
    lmvBz  = (lBzl/lsqrtrhol + lBzr/lsqrtrhor) * lmult * lmvrho;
    lmvH   = (lsqrtrhol*lHl  + lsqrtrhor*lHr)  * lmult;

    (*pret)[0] = lmvrho;
    (*pret)[1] = lmvux * lmvrho;
    (*pret)[2] = lmvuy * lmvrho;
    (*pret)[3] = lmvuz * lmvrho;
    (*pret)[4] = lmvBx;
    (*pret)[5] = lmvBy;
    (*pret)[6] = lmvBz;
    (*pret)[7] = (  lmvrho * ( lmvH + (gamma - 1.0) * 0.5
                               * (  lmvux*lmvux 
                                    + lmvuy*lmvuy + lmvuz*lmvuz ))
                    - (2.0 - gamma) * ( lmvBy * lmvBy + lmvBz * lmvBz )
                    / ( 8.0 * M_PI ) ) / gamma;
  }

  /* Roe mean values (Aslan) */

  void roe_mean_values_as(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    double lrhol,lrhouxl,lrhouyl,lrhouzl,luxl,luyl,luzl,lBxl,lByl,lBzl,lrhoEl;
    double lrhor,lrhouxr,lrhouyr,lrhouzr,luxr,luyr,luzr,lBxr,lByr,lBzr,lrhoEr;
    double lPl,lPr,lmvrho,lmvux,lmvuy,lmvuz,lmvBx,lmvBy,lmvBz,lmvP;
    double lsqrtrhol,lsqrtrhor,lmult;

    lrhol   = pul[0];
    lrhouxl = pul[1];
    lrhouyl = pul[2];
    lrhouzl = pul[3];
    lBxl    = pul[4];
    lByl    = pul[5];
    lBzl    = pul[6];
    lrhoEl  = pul[7];

    assert(lrhol > 0.0);

    luxl = lrhouxl / lrhol;
    luyl = lrhouyl / lrhol;
    luzl = lrhouzl / lrhol;
    lPl  = (gamma - 1.0) 
      * (lrhoEl - 0.5 * ( lrhouxl * luxl + lrhouyl * luyl + lrhouzl * luzl ))
      + (2.0 - gamma) * ( lByl * lByl + lBzl * lBzl ) / ( 8.0 * M_PI );

    lrhor   = pur[0];
    lrhouxr = pur[1];
    lrhouyr = pur[2];
    lrhouzr = pur[3];
    lBxr    = pur[4];
    lByr    = pur[5];
    lBzr    = pur[6];
    lrhoEr  = pur[7];

    assert(lrhor > 0.0);

    luxr = lrhouxr / lrhor;
    luyr = lrhouyr / lrhor;
    luzr = lrhouzr / lrhor;
    lPr  = (gamma - 1.0) 
      * (lrhoEr - 0.5 * ( lrhouxr * luxr + lrhouyr * luyr + lrhouzr * luzr ))
      + (2.0 - gamma) * ( lByr * lByr + lBzr * lBzr ) / ( 8.0 * M_PI );

    /* calculate Roe mean values */

    lsqrtrhol = sqrt(lrhol);
    lsqrtrhor = sqrt(lrhor);
    lmult     = 1.0 / (lsqrtrhol + lsqrtrhor);

    lmvrho = sqrt(lrhol * lrhor);
    lmvux  = (lsqrtrhol*luxl + lsqrtrhor*luxr) * lmult;
    lmvuy  = (lsqrtrhol*luyl + lsqrtrhor*luyr) * lmult;
    lmvuz  = (lsqrtrhol*luzl + lsqrtrhor*luzr) * lmult;
    lmvBx  = (lBxl/lsqrtrhol + lBxr/lsqrtrhor) * lmult * lmvrho;
    lmvBy  = (lByl/lsqrtrhol + lByr/lsqrtrhor) * lmult * lmvrho;
    lmvBz  = (lBzl/lsqrtrhol + lBzr/lsqrtrhor) * lmult * lmvrho;
    lmvP   = (lPl/lsqrtrhol  + lPr/lsqrtrhor)  * lmult * lmvrho;

    (*pret)[0] = lmvrho;
    (*pret)[1] = lmvux * lmvrho;
    (*pret)[2] = lmvuy * lmvrho;
    (*pret)[3] = lmvuz * lmvrho;
    (*pret)[4] = lmvBx;
    (*pret)[5] = lmvBy;
    (*pret)[6] = lmvBz;
    (*pret)[7] =    lmvP / (gamma - 1.0)
      + 0.5 * lmvrho * (   lmvux * lmvux
                           + lmvuy * lmvuy
                           + lmvuz * lmvuz )
      + ( lmvBy * lmvBy + lmvBz * lmvBz ) / (8.0 * M_PI);
  }

  /***************************************************************************
     output routines
  ****************************************************************************/


  /* calculate minimum */

  static double min(const double p1, const double p2)
  {
    double lret;

    if (p1 <= p2) lret = p1; else lret = p2;

    return lret;
  }

  /* calculate maximum */

  static double max(const double p1, const double p2)
  {
    double lret;

    if (p1 >= p2) lret = p1; else lret = p2;

    return lret;
  }

  /* timestep: use local cubes instead of a global one */

  double timestep(const MhdSolver::VEC1D *pu, const double pcfl,
                  const double pxstep, const int pnpoints, double gamma)
  {
    double lrhol,lrhouxl,lrhouyl,lrhouzl,lBxl,lByl,lBzl,lrhoEl,luxl,luyl,luzl,lpl;
    double lrhor,lrhouxr,lrhouyr,lrhouzr,lBxr,lByr,lBzr,lrhoEr,luxr,luyr,luzr,lpr;
    double lBx2l,lBy2l,lBz2l,lBx2r,lBy2r,lBz2r,lmax,lval;
    
    int li;

    lrhor   = pu[0][0];
    lrhouxr = pu[0][1];
    lrhouyr = pu[0][2];
    lrhouzr = pu[0][3];
    lBxr    = pu[0][4];
    lByr    = pu[0][5];
    lBzr    = pu[0][6];
    lrhoEr  = pu[0][7];

    assert(lrhor > 0.0);

    luxr = lrhouxr / lrhor;
    luyr = lrhouyr / lrhor;
    luzr = lrhouzr / lrhor;
    lpr  = (gamma - 1.0) 
      * (lrhoEr - 0.5 * ( lrhouxr * luxr + lrhouyr * luyr + lrhouzr * luzr )
         - ( lByr * lByr + lBzr * lBzr ) / ( 8.0 * M_PI ));

    assert(lpr > 0.0);

    lBx2r = lBxr * lBxr;
    lBy2r = lByr * lByr;
    lBz2r = lBzr * lBzr;

    /* timestep <= pcfl * pxstep */

    lmax = 1.0;

    for (li=1;li<pnpoints;li++)
      {
        lrhol   = lrhor;
        lrhouxl = lrhouxr;
        lrhouyl = lrhouyr;
        lrhouzl = lrhouzr;
        lBxl    = lBxr;
        lByl    = lByr;
        lBzl    = lBzr;
        lrhoEl  = lrhoEr;
        luxl    = luxr;
        luyl    = luyr;
        luzl    = luzr;
        lpl     = lpr;

        lBx2l   = lBx2r;
        lBy2l   = lBy2r;
        lBz2l   = lBz2r;

        lrhor   = pu[li][0];
        lrhouxr = pu[li][1];
        lrhouyr = pu[li][2];
        lrhouzr = pu[li][3];
        lBxr    = pu[li][4];
        lByr    = pu[li][5];
        lBzr    = pu[li][6];
        lrhoEr  = pu[li][7];

        assert(lrhor > 0.0);

        luxr = lrhouxr / lrhor;
        luyr = lrhouyr / lrhor;
        luzr = lrhouzr / lrhor;
        lpr  = (gamma - 1.0) 
          * (lrhoEr - 0.5 * (lrhouxr * luxr + lrhouyr * luyr + lrhouzr * luzr)
             - ( lByr * lByr + lBzr * lBzr ) / ( 8.0 * M_PI ));

        assert(lpr > 0.0);

        lBx2r = lBxr * lBxr;
        lBy2r = lByr * lByr;
        lBz2r = lBzr * lBzr;

        lval =   max(fabs(luxl),fabs(luxr)) 
          + sqrt((  (  max(lBx2l,lBx2r) + max(lBy2l,lBy2r)
                       + max(lBz2l,lBz2r)) / (4.0 * M_PI)
                    + gamma * max(lpl,lpr) ) / min(lrhol,lrhor));
        if (lval > lmax) lmax = lval;
      }  

    return pcfl * pxstep / lmax;
  }


  /***************************************************************************
     initialize HLLE-Solver
  ****************************************************************************/

  void initialize_hlle(void)
  {
    hll_use_roe_mean_values = 0;
  }

  /***************************************************************************
     functions implementing the HLLE-Solver
  ****************************************************************************/

  /* calculate min(parg,0) */

  static double minuspart(double parg)
  {
    double lret;

    if (parg >= 0.0) lret = 0.0; else lret = parg;

    return lret;
  }

  /* calculate max(parg,0) */

  static double pluspart(double parg)
  {
    double lret;

    if (parg <= 0.0) lret = 0.0; else lret = parg;

    return lret; 
  }


  /* numerical flux for HLLE-Solver */

  void flux_hlle(const MhdSolver::VEC1D pul,const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D lful,lfur,ldiff,lmv,lhvec;
    double lmult, lmultful, lmultfur, lmultdiff;
    double llambdamv_min, llambdal_min, llambdamv_max, llambdar_max, lbl, lbr;

    if (hll_use_roe_mean_values)
      {
        roe_mean_values(pul,pur,&lmv,gamma);
      }
    else
      {
        vadd(pul,pur,&lmv);
        svmult(0.5,lmv,&lmv);
      }

    /* Einfeldt, Munz, Roe, Sj"ogreen version */

    f(pul,&lful,gamma);
    f(pur,&lfur,gamma);
    vsub(pur,pul,&ldiff);

    lambda(0,      lmv,&llambdamv_min,gamma);
    lambda(dim - 1,lmv,&llambdamv_max,gamma);
    lambda(0,      pul,&llambdal_min,gamma);
    lambda(dim - 1,pur,&llambdar_max,gamma);

    lbl = min(llambdamv_min,llambdal_min);
    lbr = max(llambdamv_max,llambdar_max);

    lmult = pluspart(lbr) - minuspart(lbl);

    assert(lmult > 0.0);

    lmultful  =   pluspart(lbr)  / lmult;
    lmultfur  = - minuspart(lbl) / lmult;
    lmultdiff =   pluspart(lbr) * minuspart(lbl) / lmult;

    svmult(lmultful,lful,pret);
    svmult(lmultfur,lfur,&lhvec);
    vadd(*pret,lhvec,pret);
    svmult(lmultdiff,ldiff,&lhvec);
    vadd(*pret,lhvec,pret);   
  }


  /* print statistics */

  void statistics_hlle(void)
  {
    printf("\n\nundo          : %d\n\n",undo);
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      }
  } 

#ifdef _MHD8
#define _MHD8DC
#endif


  /***************************************************************************
     initialize DW-Solver
  ****************************************************************************/

  void initialize_dw(void)
  {
    visc = 0.05;
    use_hll_fix = 1;
    hll_use_roe_mean_values = 1;
  }

  /***************************************************************************
     functions implementing the Solver of Dai and Woodward
  ****************************************************************************/

  /* numerical flux for Solver of Dai and Woodward */

  void flux_dw(const MhdSolver::VEC1D pul,const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D lr[dim],lu[dim],luref,lrhs,llambdam,lmv,lhvec1;
    MhdSolver::MAT lT,lTinv;
    double lmult,llambda_max_mv,llambda_ewave_mv,lvmax_mv;
    double llambda_ewave_l,llambda_ewave_r,lphi,lpsi,lch = 0.0;
    int li,lj;

    /* calculate mean value */

    roe_mean_values(pul,pur,&lmv,gamma);

#ifdef _MHD8DC

    lmv[4] = (pul[4] + pur[4]) * 0.5;

#endif

    /* calculate upwind states */

    lambda_vec(lmv,&llambdam,gamma);
    for (li=0;li<dim;li++)
      {
        lpsi = psi(llambdam[li]);
        for (lj=0;lj<dim;lj++)
          lu[li][lj] = (1 - lpsi) * pul[lj] + lpsi * pur[lj];

#ifdef _MHD8DC

        lu[li][4] = lmv[4];

#endif
      }

    /* calculate parameter for hll-fix */

    if ((use_hll_fix) && (dim > 1))
      {
        lambda(dim - 1,lmv,&llambda_max_mv,gamma);
        lambda(ewave,  lmv,&llambda_ewave_mv,gamma);

        lvmax_mv = llambda_max_mv - llambda_ewave_mv;

        lambda(ewave,  pul,&llambda_ewave_l,gamma);
        lambda(ewave,  pur,&llambda_ewave_r,gamma);

        if (lvmax_mv < llambda_ewave_r - llambda_ewave_l)
          lch = llambda_ewave_r - llambda_ewave_l - lvmax_mv;
      }

    /* calculate flux */

    lphi = phi(lch);
  
    if (lphi >= 1.0 - eps)
      {
        flux_hlle(pul,pur,pret,gamma);
      }
    else
      {

        /* approximate left eigenvalues and solve linear system */

        r_vec(0,lmv,lr,gamma);
        l_mat(0,lmv,&lTinv,gamma);
        for (li=0;li<dim;li++)
          {
            for (lj=0;lj<dim;lj++)
              {
                lT[li][lj] = lr[lj][li];
              }
          }

        for (li=0;li<dim;li++)
          {
            lrhs[li] = 0.0;
            for (lj=0;lj<dim;lj++) lrhs[li] += lTinv[li][lj] * lu[li][lj];
          }

        mvmult(lT,lrhs,&luref);

        /* calculate flux */

        f(luref,pret,gamma);

        /* artificial viscosity */

        if (fabs(visc) > eps)
          {
            vsub(pul,pur,&lhvec1);
            lmult = visc * sqrt(vmult(lhvec1,lhvec1));
            if (lmult > 1.0) lmult = 1.0;
            svmult(lmult,lhvec1,&lhvec1);
            vadd(*pret,lhvec1,pret);
          }

        if (lphi > eps)
          {
            flux_hlle(pul,pur,&lhvec1,gamma);
            svmult(1.0 - lphi,*pret,pret);
            svmult(lphi,lhvec1,&lhvec1);
            vadd(*pret,lhvec1,pret);
          }
      }

    /* statistics and output */

    if (lch > eps)
      {
        if (verbose_mode) fprintf(stderr,"flux_dw: hll-fix used!\n");
        hllcount++;
      }
  }


  /* print statistics */

  void statistics_dw(void)
  {
    printf("\n\nundo          : %d\n\n",undo);
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      }

    if (use_hll_fix)
      printf("hllcount : %d\n\n",hllcount);
  }

  /***************************************************************************
     initialize Roe-Solver
  ****************************************************************************/

  void initialize_roe(void)
  {
    docheck = 0;
    use_hll_fix = 1;
    hll_use_roe_mean_values = 1;
  }

  /***************************************************************************
     functions implementing the Roe-Solver
  ****************************************************************************/

  /* check Roe matrix ( F(U_l) - F(U_r) ?=? A(U_l,U_r) (U_l - U_r) ) */

  double check_roe_matrix(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, const MhdSolver::MAT pmat, double gamma)
  {
    MhdSolver::VEC1D ldiff_u,ldiff_fu,lful,lfur,lh;

    f(pul,&lful,gamma);
    f(pur,&lfur, gamma);
    vsub(lful,lfur,&ldiff_fu);
    vsub(pul,pur,&ldiff_u);
    mvmult(pmat,ldiff_u,&lh);
    vsub(ldiff_fu,lh,&lh);

    return(sqrt(vmult(lh,lh)));
  }

  /* calc parameter for entropy-fix */

  static double check_entropy_fix(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, double gamma)
  {
    int li;
    MhdSolver::VEC1D llambdal,llambdar;
    double lret = 0;

    lambda_vec(pul,&llambdal,gamma);
    lambda_vec(pur,&llambdar,gamma);
    for (li=0;li<dim;li++)
      {
        if ((llambdal[li] < 0.0) && (llambdar[li] > 0.0))
          lret += llambdar[li] - llambdal[li];
      }

    return lret;
  }

  /* numerical flux for Roe-Solver */

  void flux_roe(const MhdSolver::VEC1D pul,const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D lful,lfur,ldiff,lr[dim],lmv,lalpha,llambda,lhvec1,lhvec2;
    MhdSolver::MAT lroemat,linv;
    double lcheck,llambda_max_mv,llambda_ewave_mv,lvmax_mv;
    double llambda_ewave_l,llambda_ewave_r,lphi,lce = 0.0,lch = 0.0;
    int li;

    roe_mean_values(pul,pur,&lmv,gamma);

    /* calculate parameter for hll-fix */

    if ((use_hll_fix) && (dim > 1))
      {
        lambda(dim - 1,lmv,&llambda_max_mv,gamma);
        lambda(ewave,  lmv,&llambda_ewave_mv,gamma);

        lvmax_mv = llambda_max_mv - llambda_ewave_mv;

        lambda(ewave,  pul,&llambda_ewave_l,gamma);
        lambda(ewave,  pur,&llambda_ewave_r,gamma);

        if (lvmax_mv < llambda_ewave_r - llambda_ewave_l) 
          lch = llambda_ewave_r - llambda_ewave_l - lvmax_mv;
      }

    /* calculate roe-flux */

    if (docheck)
      {
        Df(lmv,&lroemat,gamma);
        lcheck = check_roe_matrix(pul,pur,lroemat,gamma);
        if (lcheck > eps)
          {
            printf("   check_roe_matrix: lcheck = %le\n",lcheck);
            check_rmv_errors++;
            if (lcheck > max_check_rmv_error) max_check_rmv_error = lcheck;
          }
      }
    vsub(pur,pul,&ldiff);

    r_vec(0,lmv,lr,gamma);
    l_mat(0,lmv,&linv,gamma);

    lambda_vec(lmv,&llambda,gamma);

    mvmult(linv,ldiff,&lalpha);

    if (use_entropy_fix)
      lce = check_entropy_fix(pul, pur,gamma);

    lphi = phi(lce + lch);

    if (lphi >= 1.0 - eps)
      flux_hlle(pul,pur,pret,gamma);
    else
      {
        svmult(fabs(llambda[0])*lalpha[0],lr[0],&lhvec2);
        for (li=1;li<dim;li++)
          {
            svmult(fabs(llambda[li])*lalpha[li],lr[li],&lhvec1);
            vadd(lhvec1,lhvec2,&lhvec2);
          }
        svmult(0.5,lhvec2,&lhvec2);

        f(pul,&lful,gamma);
        f(pur,&lfur,gamma);
        vadd(lful,lfur,&lhvec1);
        svmult(0.5,lhvec1,&lhvec1);

        vsub(lhvec1,lhvec2,pret);

        if (lphi > eps)
          {
            flux_hlle(pul,pur,&lhvec1,gamma);
            svmult(1.0 - lphi,*pret,pret);
            svmult(lphi,lhvec1,&lhvec1);
            vadd(*pret,lhvec1,pret);
          }
      }

    /* statistics and output */

    if (lce > 0.5 * eps)
      {
        if (verbose_mode) fprintf(stderr,"flux_roe: hll-entropy-fix used!\n");
        entropy_fixes++;
      }
    if (lch > 0.5 * eps)
      {
        if (verbose_mode) fprintf(stderr,"flux_roe: hll-fix used!\n");
        hllcount++;
      }

  }


  /* print statistics */

  void statistics_roe(void)
  {
    printf("\n\nundo          : %d\n\n",undo);
    if (use_hll_fix)
      printf("hll_fixes     : %d\n\n",hllcount);
    if (use_entropy_fix)
      printf("entropy_fixes : %d\n\n",entropy_fixes);
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      }
    if (docheck)
      {
        printf("check_rmv_errors          : %d\n",check_rmv_errors);
        printf("max_check_rmv_error       : %le\n\n",max_check_rmv_error);
      }
  }

  /***************************************************************************
     initialize BCT-Solver
  ****************************************************************************/

  void initialize_bct(void)
  {
    use_hll_fix = 1;
    use_point_fix = 1;
    use_velocity_fix = 1;
    use_reference_fix = 1;
    hll_use_roe_mean_values = 0;
  }

  /***************************************************************************
     functions implementing the BCT-algorithm
  ****************************************************************************/
  typedef                            /* type of evaluation in eval_plmm_int: */
  enum {e_normal, e_colltrans}     /*   e_normal   : normal integration    */
  EVALUATION_MODE;                 /*   e_colltrans: special treatment as  */
                                   /*                "collapsed transonic" */
                                   /*                wave.                 */

  typedef int COLLINF[dim-1];        /* "collapse information"               */

  typedef struct evalinf             /* struct for approximated eigenvalue   */
  {
    unsigned short int parts;
    double x[4],y[4];
    EVALUATION_MODE emode;
    double nu;
  }EVALINF;

  typedef MhdSolver::VEC1D EVEC_APPROX[dim];      /* array of approximated eigenvectors   */

  typedef EVALINF EVAL_APPROX[dim];  /* array of approximated eigenvalues    */


  /* calculate sign */

  static double sgn(double px)
  {
    double lret = 1.0;

    if (px < 0.0) lret = -1.0;

    return lret;
  }

  /* choose reference state */
  
  void ref_state(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur,
                 double *ppsi, double *psigmabar, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D ldiff;
    double llambdal,llambdar;
    int li;

    assert((ewave >= 0) && (ewave <= dim));

    lambda(ewave,pul,&llambdal,gamma);
    lambda(ewave,pur,&llambdar,gamma);

    *psigmabar = 0.5 * (llambdal + llambdar);

    if (*psigmabar >= 0.0)
      for(li=0;li<dim;li++) (*pret)[li] = pul[li];
    else
      for(li=0;li<dim;li++) (*pret)[li] = pur[li];

    *ppsi = psi(*psigmabar);

    if ((*ppsi > eps) && (*ppsi < 1.0 - eps))
      {
        vsub(pul,pur,&ldiff);
        if (sqrt(vmult(ldiff,ldiff)) <= eps) *ppsi = psi(sgn(*psigmabar));  
      }
  }

  /* approximate phase-space-solution */

  void approx_phasespacesol(const int psign, const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, 
                            MhdSolver::VEC1D *palpha, EVEC_APPROX *pevec, double gamma)
  {
    MhdSolver::VEC1D ldiff, lavg;
    MhdSolver::MAT linv;
    int li;

    vsub(pur,pul,&ldiff);
    vadd(pul,pur,&lavg);
    svmult(0.5,lavg,&lavg);

    r_vec(psign,lavg,*pevec,gamma);
    l_mat(psign,lavg,&linv,gamma);

    mvmult(linv,ldiff,palpha);

    for (li=0;li<dim;li++)
      if ((*palpha)[li] < 0.0)
        {
          (*palpha)[li] *= -1.0;
          svmult(-1.0,(*pevec)[li],&((*pevec)[li]));
        } 
  }

  /* evaluate a cubic polynom */

  double eval_cubic_pol(const double pbeta[4], const double px)
  {
    return (  pbeta[3] * px * px * px 
              + pbeta[2] * px * px
              + pbeta[1] * px
              + pbeta[0] );
  }

  /* approximate eigenvalues along approximated phase-space-solution */

  void approx_eigenvalues(const int psign, const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur,
                          const MhdSolver::VEC1D palpha, const EVEC_APPROX pevec,
                          EVAL_APPROX *peval, double gamma)
  {
    int li,lj,lidx,ln[2];
    double la,lb,lc,ld,lalpha,lbeta[4],lrad,lx[2],lh;
    MhdSolver::VEC1D llambdal,llambdar,lgradlambdal[dim],lgradlambdar[dim];

    lambda_vec(pul,&llambdal,gamma);
    lambda_vec(pur,&llambdar,gamma);

    grad_lambda_vec(psign,pul,lgradlambdal,gamma);
    grad_lambda_vec(psign,pur,lgradlambdar,gamma);

    for (li=0;li<dim;li++)
      {

        if ((*peval)[li].emode == e_normal)
          {
            lalpha = palpha[li];

            if (fabs(lalpha) < eps) (*peval)[li].parts = 0;
            else 
              {
                la = llambdal[li];
                lc = llambdar[li];

                lb = vmult(lgradlambdal[li],pevec[li]);
                ld = vmult(lgradlambdar[li],pevec[li]);

      
                /* calculate cubic interpolation of lambda[li] */

                assert(lalpha != 0);

                lbeta[0] = la;
                lbeta[1] = lb;
                lbeta[2] = 3.0 * lc / (lalpha * lalpha)
                  - ld / lalpha
                  - 2.0 * lbeta[1] / lalpha
                  - 3.0 * lbeta[0] / (lalpha * lalpha);
                lbeta[3] = ( ld / lalpha - 2.0 * lbeta[2] - lbeta[1] / lalpha)
                  / (3.0 * lalpha);

                /* calculate piecewise linear interpolation */

                lrad = 4.0 * lbeta[2] * lbeta[2] - 12.0 * lbeta[1] * lbeta[3];
                if (fabs(lrad) <= eps) lrad = 0.0;
                if (lrad < 0.0) ln[0] = ln[1] = 0;
                else
                  {
                    if (fabs(lbeta[3]) > eps)
                      {
                        lx[0] = (-2.0 * lbeta[2] - sqrt(lrad)) / (6.0 * lbeta[3]);
                        lx[1] = (-2.0 * lbeta[2] + sqrt(lrad)) / (6.0 * lbeta[3]);
                        if (lx[0] > lx[1])
                          {
                            lh    = lx[0];
                            lx[0] = lx[1];
                            lx[1] = lh;
                          }
                        for (lj=0;lj<2;lj++)
                          if ((lx[lj] > 0.0) && (lx[lj] < lalpha)) ln[lj] = 1;
                          else ln[lj] = 0;
                        if ( fabs(lx[0]-lx[1]) <= eps ) ln[1] = 0;
                      }
                    else if (fabs(lbeta[2]) > eps)
                      {
                        lx[0] = - 0.5 * lbeta[1] / lbeta[2];
                        if ((lx[0] > 0.0) && (lx[0] < lalpha)) ln[0] = 1; else ln[0] = 0;
                        ln[1] = 0;
                      }
                    else ln[0] = ln[1] = 0;
                  }
                (*peval)[li].parts = (short int)(1 + ln[0] + ln[1]);

                (*peval)[li].x[0] = 0.0;
                (*peval)[li].y[0] = eval_cubic_pol(lbeta,0.0);
                lidx = 1;    
                for (lj=0;lj<2;lj++)
                  if (ln[lj])
                    {
                      (*peval)[li].x[lidx] = lx[lj];
                      (*peval)[li].y[lidx] = eval_cubic_pol(lbeta,lx[lj]);
                      lidx++;
                    }
                (*peval)[li].x[lidx] = lalpha;
                (*peval)[li].y[lidx] = eval_cubic_pol(lbeta,lalpha);
              }
          }
      }  
  }

  /* find zero for line segment [(plx,ply),(prx,pry)] with 0 <= plx < prx 
     and ((ply * pry <= eps) or (|ply| <= eps) or (|pry| <= eps))          */

  double find_zero(const double plx, const double ply,
                   const double prx, const double pry)
  {
    double lzero,lbeta[2];

    assert((0.0 <= plx) && (plx < prx));
    assert((ply * pry <= eps) || (fabs(ply) <= eps) || (fabs(pry) <= eps));

    if (fabs(ply) <= eps) lzero = plx;
    else if (fabs(pry) <= eps) lzero = prx;
    else if (ply * pry <= -eps)
      {
        if (plx < eps)
          {
            lbeta[0] = ply;
            lbeta[1] = (pry - ply) / prx;
          }
        else
          {
            lbeta[0] = (plx * pry - ply * prx) / (plx - prx);
            lbeta[1] = (ply - lbeta[0]) / plx;
          }
        lzero = - lbeta[0] / lbeta[1];
      }
    else if (fabs(ply) < fabs(pry)) lzero = plx; else lzero = prx; 

    return lzero;
  }

  /* evaluate piecewise linear min/max-integral */

  double eval_plmm_int(const double psigmabar, EVALINF plambdabar)
  {
    double lret = 0.0,lint,lleftx,lrightx,llefty,lrighty;
    int li;

    /* - max(f,0) = min(-f,0) */

    if (psigmabar < 0.0)
      for (li=0;li<=plambdabar.parts;li++) plambdabar.y[li] *= -1.0;

    for (li=0;li<plambdabar.parts;li++)
      {
        /* integrate line segment [(lleftx,llefty),(lrightx,lrighty)] */
        if (plambdabar.y[li] <= 0.0)
          {
            lleftx  = plambdabar.x[li];
            llefty  = plambdabar.y[li];
            if (plambdabar.y[li+1] <= 0.0)
              {
                lrightx = plambdabar.x[li+1];
                lrighty = plambdabar.y[li+1];
              }
            else
              {
                lrightx = find_zero(plambdabar.x[li],plambdabar.y[li],
                                    plambdabar.x[li+1],plambdabar.y[li+1]);
                lrighty = 0.0;
              }
          }
        else
          {
            lrightx = plambdabar.x[li+1];
            llefty  = 0.0;
            if (plambdabar.y[li+1] > 0.0)
              {
                lleftx  = lrightx;
                lrighty = 0.0;
              }
            else
              {
                lleftx  = find_zero(plambdabar.x[li],plambdabar.y[li],
                                    plambdabar.x[li+1],plambdabar.y[li+1]);
                lrighty = plambdabar.y[li+1];
              }
          }
        lint = 0.5 * (lrightx - lleftx) * (llefty + lrighty);

        lret += lint;
      }

    return lret;
  }

  /* get sign of vax2 - cs2 */

  int sign_vax2mcs2(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, double ppoint, double gamma)
  {
    double lrho,lrhoux,lrhouy,lrhouz,lrhoE,lp,lBx,lvax2,lcs2,ldiff;
    int lret,lincr;

    assert((dim == 7) || (dim == 8));

    if (dim == 8) lincr = 1;
    else lincr = 0;

    assert(fabs(pul[4+lincr]) <= 10.0 * eps);
    assert(fabs(pul[5+lincr]) <= 10.0 * eps);
    assert(fabs(pur[4+lincr]) <= 10.0 * eps);
    assert(fabs(pur[5+lincr]) <= 10.0 * eps);

    lrho   = (1.0 - ppoint) * pul[0] + ppoint * pur[0];
    lrhoux = (1.0 - ppoint) * pul[1] + ppoint * pur[1];
    lrhouy = (1.0 - ppoint) * pul[2] + ppoint * pur[2];
    lrhouz = (1.0 - ppoint) * pul[3] + ppoint * pur[3];
    lrhoE  = (1.0 - ppoint) * pul[6+lincr] + ppoint * pur[6+lincr];

    assert(lrho > 0.0);

    lp = (gamma - 1.0) * (lrhoE - 0.5 * (   lrhoux * lrhoux 
                                            + lrhouy * lrhouy
                                            + lrhouz * lrhouz ) / lrho );

    assert(lp > 0.0);
  
    lvax2 = -17.42;
#ifdef _MHD
    lvax2 = Bx * Bx / (4.0 * M_PI * lrho);
#endif
#ifdef _MHD8
    lBx   = 0.5 * (pul[4] + pur[4]);
    lvax2 = lBx * lBx / (4.0 * M_PI * lrho);
#endif
    assert(lvax2 > 0.0);

    lcs2  = gamma * lp / lrho;

    ldiff = lvax2 - lcs2;

    if (ldiff < -eps)
      lret = -1;
    else if (ldiff <= eps)
      lret = 0;
    else 
      lret = 1;

    return lret;
  }

  /* get data for velocity conservation fix */

  void velocity_fix(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur,
                    int *pnr_of_points, int *psign, double *ppoint, double gamma)
  {
    double la,lb,lc,ld,lrad,lsqrt,lpoint[2],lnr_of_points = 0;
    double lrhol,lrhouxl,lrhouyl,lrhouzl,lrhoEl;
    double lrhor,lrhouxr,lrhouyr,lrhouzr,lrhoEr;
    double lrhoElrhol,lrhoElrhor,lrhoErrhol,lrhoErrhor;
    double lrhoul2,lrhour2,lrhoulrhour,lh,lBx;

    int li,lj,lincr;

    assert((dim == 7) || (dim == 8));

    if (dim == 8) lincr = 1;
    else lincr = 0;

    lrhol   = pul[0];
    lrhouxl = pul[1];
    lrhouyl = pul[2];
    lrhouzl = pul[3];
    lrhoEl  = pul[6+lincr];

    lrhor   = pur[0];
    lrhouxr = pur[1];
    lrhouyr = pur[2];
    lrhouzr = pur[3];
    lrhoEr  = pur[6+lincr];

    assert(lrhol > 0.0);
    assert(lrhor > 0.0);

    lrhoElrhol  = lrhoEl * lrhol;
    lrhoElrhor  = lrhoEl * lrhor;
    lrhoErrhol  = lrhoEr * lrhol;
    lrhoErrhor  = lrhoEr * lrhor;
    lrhoul2     = lrhouxl * lrhouxl + lrhouyl * lrhouyl + lrhouzl * lrhouzl;
    lrhour2     = lrhouxr * lrhouxr + lrhouyr * lrhouyr + lrhouzr * lrhouzr;
    lrhoulrhour = lrhouxl * lrhouxr + lrhouyl * lrhouyr + lrhouzl * lrhouzr;

#ifdef _MHD
    ld = Bx * Bx / (4.0 * M_PI * gamma * (gamma - 1.0));
#endif
#ifdef _MHD8
    lBx= 0.5 * (pul[4] + pur[4]);
    ld = lBx * lBx / (4.0 * M_PI * gamma * (gamma - 1.0));
#endif

    la =   lrhoElrhol - lrhoElrhor - lrhoErrhol + lrhoErrhor
      - 0.5 * ( lrhoul2 + lrhour2 ) + lrhoulrhour;
    lb = - 2.0 * lrhoElrhol + lrhoElrhor + lrhoErrhol
      + lrhoul2 - lrhoulrhour + ld * (lrhol - lrhor);
    lc =   lrhoElrhol - 0.5 * lrhoul2 - ld * lrhol;

    if (fabs(la) <= eps)
      {
        if (fabs(lb) <= eps) lnr_of_points = 0;
        else
          {
            lnr_of_points = 1;
            lpoint[0] = - lc / lb;
          }
      }
    else
      {
        lrad = lb * lb - 4.0 * la * lc;

        if (lrad < -eps)
          {
            lnr_of_points = 0;
          }
        else if (lrad <= eps)
          {
            lnr_of_points = 1;
            lpoint[0] = - lb * 0.5 / la;
          }
        else
          {
            lnr_of_points = 2;
            lsqrt = sqrt(lrad);
            lpoint[0] = ( - lb - lsqrt ) * 0.5 / la;
            lpoint[1] = ( - lb + lsqrt ) * 0.5 / la;
            if (lpoint[0] > lpoint[1])
              {
                lh        = lpoint[0];
                lpoint[0] = lpoint[1];
                lpoint[1] = lh;
              }
          }
      }

    lj = 0;
    *pnr_of_points = 0;
    for (li=0;li<lnr_of_points;li++)
      {
        if ((lpoint[li] > 0.0) && (lpoint[li] < 1.0))
          {
            assert(0 == sign_vax2mcs2(pul,pur,lpoint[li],gamma));
            ppoint[lj] = lpoint[li];
            lj++;
            (*pnr_of_points)++;
          }
      }
    psign[0] = sign_vax2mcs2(pul,pur,0.0,gamma);
    psign[1] = sign_vax2mcs2(pul,pur,1.0,gamma);
  }

  /* calc parameter for entropy-fix */


  /* numerical flux for BCT-algorithm */

  void flux_bct(const MhdSolver::VEC1D pul, const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D lalpha, luref, lhvec, ldvec, lul, lur, lmv;
    EVEC_APPROX levec;
    EVAL_APPROX leval;
    COLLINF lcollapse;
    double lsigmabar,lpoint[2],lint,lvmax_mv,lphi,lpsi,lce = 0.0,lch = 0.0;
    double llambda_max_mv, llambda_ewave_mv, llambda_ewave_l, llambda_ewave_r;
    double lBx, lBy, lBz;
    int li,lad,lrfix_request = 0,lpfix_request = 0,lpsasteps = 1, lzcount, lpfpoint=-1;
    int lvfix_request = 0,lnr_of_points,lsign[2],lincr,lcsign = 0;

    assert((dim == 7) || (dim == 8));
    if (dim == 8)
      {
        lincr = 1;
        lBx = 0.5 * (pul[4] + pur[4]);
      }
    else lincr = 0;

    /* calculate parameter for hll-fix */

    if ((use_hll_fix) && (dim > 1))
      {
        vadd(pul,pur,&lmv);
        svmult(0.5,lmv,&lmv);

        lambda(dim - 1,lmv,&llambda_max_mv,gamma);
        lambda(ewave,  lmv,&llambda_ewave_mv,gamma);

        lvmax_mv = llambda_max_mv - llambda_ewave_mv;

        lambda(ewave,  pul,&llambda_ewave_l,gamma);
        lambda(ewave,  pur,&llambda_ewave_r,gamma);

        if (lvmax_mv < llambda_ewave_r - llambda_ewave_l)
          lch = llambda_ewave_r - llambda_ewave_l - lvmax_mv;
      }
  
    /* initialization */

    for (li=0;li<dim;li++) (*pret)[li] = 0.0;

    ref_state(pul,pur,&lpsi,&lsigmabar,&luref,gamma);
    if ((lpsi > eps) && (lpsi < 1.0 - eps) && (use_reference_fix))
      {
        f(pul,pret,gamma);
        f(pur,&lhvec,gamma);
        svmult(1.0 - lpsi,*pret,pret);
        svmult(lpsi,lhvec,&lhvec);
        vadd(*pret,lhvec,pret);
        if (verbose_mode) fprintf(stderr,"flux_bct: using reference fix!\n");
        reference_fixes++;
        lrfix_request = 1;
      }
    else f(luref,pret,gamma);

    /* (adaptive) approximation of phase-space solution */

    vsub(pur,pul,&ldvec);

    if ((use_point_fix) || (use_velocity_fix))
      {
        lzcount  = 0;
        lpfpoint = -1;
        lBy = pul[4+lincr];
        lBz = pul[5+lincr];
        if (sqrt(lBy * lBy + lBz * lBz) <= eps)
          {
            lzcount++;
            lpfpoint = 0;
          }
        lBy += pur[4+lincr];
        lBz += pur[5+lincr];
        lBy *= 0.5;
        lBz *= 0.5;
        if (sqrt(lBy * lBy + lBz * lBz) <= eps)
          {
            lzcount++;
            lpfpoint = 1;
          }
        lBy = pur[4+lincr];
        lBz = pur[5+lincr];
        if (sqrt(lBy * lBy + lBz * lBz) <= eps)
          {
            lzcount++;
            lpfpoint = 2;
          }
        if (lzcount > 0)
          {
            if (lzcount == 1)
              {
                if (use_point_fix)
                  {
                    lpfix_request = 1;
                    if (verbose_mode) fprintf(stderr,"flux_bct: using point fix!\n");
                    point_fixes++;
                    if (lpfpoint == 1)
                      {
                        lpsasteps = 2;
                        lpoint[0] = 0.5;
                      }
                  }
              }
#ifdef _MHD8
            else if ((use_velocity_fix) && (lBx > eps))
#endif
#ifdef _MHD
              else if ((use_velocity_fix) && (Bx > eps))
#endif
                {
                  lvfix_request = 1;
                  velocity_fix(pul,pur,&lnr_of_points,lsign,lpoint,gamma);
                  if (lnr_of_points > 0)
                    {
                      lpsasteps = lnr_of_points + 1;
                      velocity_fixes++;
                      if (verbose_mode) fprintf(stderr,"flux_bct: using velocity fix!\n");
                    }
                }
          }
      }

    for (lad=0;lad<lpsasteps;lad++)
      {
        if (lad == 0) vcopy(pul,&lul);
        else
          {
            svmult(lpoint[lad-1],ldvec,&lhvec);
            vadd(pul,lhvec,&lul);
          }
        if (lad == lpsasteps-1) vcopy(pur,&lur);
        else
          {
            svmult(lpoint[lad],ldvec,&lhvec);
            vadd(pul,lhvec,&lur);
          }

        if (lpfix_request)
          {
            switch (lpfpoint)
              {
              case 0 : assert(lad == 0);
                if (fabs(lur[4+lincr]) > 0.0) lul[4+lincr] = 1.5 * eps;
                if (lur[4+lincr] < 0.0) lul[4+lincr] *= -1.0;
                if (fabs(lur[5+lincr]) > 0.0) lul[5+lincr] = 1.5 * eps;
                break;
              case 1 : switch (lad)
                {
                case 0 : if (fabs(lul[4+lincr]) > 0.0) lur[4+lincr] = 1.5 * eps;
                  if (lul[4+lincr] < 0.0) lur[4+lincr] *= -1.0;
                  if (fabs(lul[5+lincr]) > 0.0) lur[5+lincr] = 1.5 * eps;
                  break;
                case 1 : if (fabs(lur[4+lincr]) > 0.0) lul[4+lincr] = 1.5 * eps;
                  if (lur[4+lincr] < 0.0) lul[4+lincr] *= -1.0;
                  if (fabs(lur[5+lincr]) > 0.0) lul[5+lincr] = 1.5 * eps;
                  break;
                default: {assert(0 > 1); abort();} //uwe
                }
                break;
              case 2 : assert(lad == 0);
                if (fabs(lul[4+lincr]) > 0.0) lur[4+lincr] = 1.5 * eps;
                if (lul[4+lincr] < 0.0) lur[4+lincr] *= -1.0;
                if (fabs(lul[5+lincr]) > 0.0) lur[5+lincr] = 1.5 * eps;
                break;
              default: {assert(0 > 1); abort();}
              }
          }
        else if (lvfix_request)
          {
            switch (lnr_of_points)
              {
              case 0:
                if (0 != lsign[0]) lcsign = lsign[0];
                else if (0 != lsign[1]) lcsign = lsign[1];
                else lcsign = sign_vax2mcs2(lul,lur,0.5,gamma);
                break;
              case 1:
                if (lad == 0)
                  {
                    if (0 != lsign[0]) lcsign = lsign[0];
                    else lcsign = sign_vax2mcs2(lul,lur,0.5,gamma);
                  }
                else
                  {
                    assert(lad == 1);
                    if (0 != lsign[1]) lcsign = lsign[1];
                    else lcsign = sign_vax2mcs2(lul,lur,0.5,gamma);
                  }
                break;
              case 2:
                if (lad == 0) lcsign = lsign[0];
                else if (lad == 1) lcsign = -lsign[0];
                else
                  {
                    assert(lad == 2);
                    lcsign = lsign[1];
                  }
                break;
              default:
                exit(17);
              }
          }
    
        approx_phasespacesol(lcsign,lul,lur,&lalpha,&levec,gamma);
        for (li=0;li<dim;li++)
          {
            lcollapse[li] = 0;
            leval[li].emode = e_normal;
          }
        approx_eigenvalues(lcsign,lul,lur,lalpha,levec,&leval,gamma);

        for (li=0;li<dim;li++)
          {
            if (lrfix_request)
              {
                lint = eval_plmm_int( 1.0,leval[li]);
                svmult(lint * (1.0 - lpsi),levec[li],&lhvec);
                vadd(*pret,lhvec,pret);
                lint = eval_plmm_int(-1.0,leval[li]);
                svmult(lint * lpsi,levec[li],&lhvec);
                vadd(*pret,lhvec,pret);
              }
            else
              {
                lint = eval_plmm_int(lsigmabar,leval[li]);
                svmult(lint,levec[li],&lhvec);
                vadd(*pret,lhvec,pret);
              }
          }
      }

    if (use_entropy_fix)
      lce = check_entropy_fix(pul,pur,gamma);

    lphi = phi(lce + lch);

    if (lphi > eps)
      {
        flux_hlle(pul,pur,&lhvec,gamma);
        svmult(1.0 - lphi,*pret,pret);
        svmult(lphi,lhvec,&lhvec);
        vadd(*pret,lhvec,pret);
      }

    /* statistics and output */

    if (lch > 0.5 * eps)
      {
        if (verbose_mode) fprintf(stderr,"flux_bct: hll-fix used!\n");
        hllcount++;
      }
    if (lce > 0.5 * eps)
      {
        if (verbose_mode) fprintf(stderr,"flux_bct: hll-entropy-fix used!\n");
        entropy_fixes++;
      }
  }


  /* print statistics */

  void statistics_bct(void)
  {
    printf("\n\nundo : %d\n\n",undo);
    if (use_entropy_fix)  printf("entropy_fixes  : %d\n",entropy_fixes);
    if (use_point_fix) printf("point_fixes    : %d\n",point_fixes);
    if (use_velocity_fix) printf("velocity_fixes : %d\n",velocity_fixes);
    if (use_reference_fix) printf("reference_fixes: %d\n",reference_fixes);
    printf("\n");
    if (use_hll_fix)
      {
        printf("hllcount : %d\n\n",hllcount);
      }
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      } 
  }


  /***************************************************************************
     initialize HLLEMG-Solver
  ****************************************************************************/

  void initialize_hllemg(void)
  {
    hll_use_roe_mean_values = 0;
  }

  /***************************************************************************
     functions implementing the HLLEMG-Solver
  ****************************************************************************/


  /* numerical flux for HLLEMG-Solver */

  void flux_hllemg(const MhdSolver::VEC1D pul,const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::VEC1D lful,lfur,ldiff,lmv,lhvec,lhvec1,lhvec2,ll_l,ll_r,lr_l,lr_r,llambdamv;
    double lmult,lmultful,lmultfur,lmultdiff,lphi,lbl,lbr;
    double lalpha_l,lalpha_r,lvel,ldelta,llambdal_min,llambdar_max;

    assert(dim >= 3);

    if (hll_use_roe_mean_values)
      {
        roe_mean_values(pul,pur,&lmv,gamma);
      }
    else
      {
        vadd(pul,pur,&lmv);
        svmult(0.5,lmv,&lmv);
      }

    /* Einfeldt, Munz, Roe, Sj"ogreen version of mean values */

    f(pul,&lful,gamma);
    f(pur,&lfur,gamma);
    vsub(pur,pul,&ldiff);

    lambda_vec(lmv,&llambdamv,gamma);
    lambda(0,      pul,&llambdal_min,gamma);
    lambda(dim - 1,pur,&llambdar_max,gamma);

    l(0,      0,lmv,&ll_l,gamma);
    l(dim - 1,0,lmv,&ll_r,gamma);
    r(0,      0,lmv,&lr_l,gamma);
    r(dim - 1,0,lmv,&lr_r,gamma);

    lbl = min(llambdamv[0],llambdal_min);
    lbr = max(llambdamv[dim-1],llambdar_max);

    lmult = pluspart(lbr) - minuspart(lbl);

    assert(lmult > 0.0);

    lmultful  =   pluspart(lbr)  / lmult;
    lmultfur  = - minuspart(lbl) / lmult;
    lmultdiff =   pluspart(lbr) * minuspart(lbl) / lmult;

    lphi      = phi(llambdamv[1] - llambdamv[0]);

    lalpha_l  = vmult(ll_l,ldiff);
    lalpha_r  = vmult(ll_r,ldiff);
    lvel      = llambdamv[ewave] - llambdamv[0];
    ldelta    = lphi * lvel / (lvel + 0.5 * fabs(lbl + lbr));

    svmult(lmultful,lful,pret);
    svmult(lmultfur,lfur,&lhvec);
    vadd(*pret,lhvec,pret);
    svmult(1.0 - ldelta,ldiff,&lhvec);
    svmult(lalpha_l,lr_l,&lhvec1);
    svmult(lalpha_r,lr_r,&lhvec2);
    vadd(lhvec1,lhvec2,&lhvec1);
    svmult(ldelta,lhvec1,&lhvec1);
    vadd(lhvec,lhvec1,&lhvec);
    svmult(lmultdiff,lhvec,&lhvec);
    vadd(*pret,lhvec,pret);   
  }


  /* print statistics */

  void statistics_hllemg(void)
  {
    printf("\n\nundo          : %d\n\n",undo);
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      }
  } 

  /***************************************************************************
     initialize HLLEML-Solver
  ****************************************************************************/

  void initialize_hlleml(void)
  {
    hll_use_roe_mean_values = 0;
  }

  /***************************************************************************
     functions implementing the HLLEML-Solver
  ****************************************************************************/


  /* numerical flux for HLLEML-Solver */

  void flux_hlleml(const MhdSolver::VEC1D pul,const MhdSolver::VEC1D pur, MhdSolver::VEC1D *pret, double gamma)
  {
    MhdSolver::MAT lmat;
    MhdSolver::VEC1D lful,lfur,ldiff,lmv,lhvec,llambdamv,lvel,lalpha,ldelta,lr[dim];
    double lmult,lmultful,lmultfur,lmultdiff,lphi1,lphi2 = 0.0;
    double lbl,lbr,llambdal_min,llambdar_max,ldiv;
    int li,lj;

    assert(dim>=3);

    if (hll_use_roe_mean_values)
      {
        roe_mean_values(pul,pur,&lmv,gamma);
      }
    else
      {
        vadd(pul,pur,&lmv);
        svmult(0.5,lmv,&lmv);
      }

    f(pul,&lful,gamma);
    f(pur,&lfur,gamma);
    vsub(pur,pul,&ldiff);

    lambda_vec(lmv,&llambdamv,gamma);
    lambda(0,pul,&llambdal_min,gamma);
    lambda(dim-1,pur,&llambdar_max,gamma);

    l_mat(0,lmv,&lmat,gamma);
    r_vec(0,lmv,lr,gamma);

    lbl = min(llambdamv[0],llambdal_min);
    lbr = max(llambdamv[dim-1],llambdar_max);

    lmult = pluspart(lbr) - minuspart(lbl);

    assert(lmult > 0.0);

    lmultful  =   pluspart(lbr)  / lmult;
    lmultfur  = - minuspart(lbl) / lmult;
    lmultdiff =   pluspart(lbr) * minuspart(lbl) / lmult;

    lphi1 = phi(llambdamv[1] - llambdamv[0]);

#ifdef _MHD8
    lphi2 = phi(lmv[4] * lmv[4] + lmv[5] * lmv[5] + lmv[6] * lmv[6]);
#endif
#ifdef _MHD
    lphi2 = phi(Bx * Bx + lmv[4] * lmv[4] + lmv[5] * lmv[5]);
#endif

    for (li=1;li<dim-1;li++)
      {
        lalpha[li] = 0.0;

        for (lj=0;lj<dim;lj++)
          lalpha[li] += lmat[li][lj] * ldiff[lj];
      }

    for (li=0;li<ewave;li++)
      lvel[li] = llambdamv[ewave] - llambdamv[li];

    ldiv = lvel[0] + (1.0 - lphi2) * fabs(llambdamv[ewave]);
    for (li=1;li<ewave;li++)
      ldelta[li]  = lphi1 * lvel[0] / (ldiv + lvel[li]);
    ldelta[ewave] = lphi1 * lvel[0] / (lvel[0] + fabs(llambdamv[ewave]));

    svmult(lmultful,lful,pret);
    svmult(lmultfur,lfur,&lhvec);
    vadd(*pret,lhvec,pret);
    svmult(lmultdiff,ldiff,&lhvec);
    vadd(*pret,lhvec,pret);

    for (li=1;li<ewave;li++)
      {
        svmult(lmultdiff*ldelta[li]*lalpha[li],lr[li],&lhvec);
        vsub(*pret,lhvec,pret);
        svmult(lmultdiff*ldelta[li]*lalpha[dim-1-li],lr[dim-1-li],&lhvec);
        vsub(*pret,lhvec,pret);
      }
    svmult(lmultdiff*ldelta[ewave]*lalpha[ewave],lr[ewave],&lhvec);
    vsub(*pret,lhvec,pret);
  }


  /* print statistics */

  void statistics_hlleml(void)
  {
    printf("\n\nundo          : %d\n\n",undo);
    if (check_2D_conservation)
      {
        printf("conservation_errors_2D    : %d\n",conservation_errors_2D);
        printf("max_conservation_error_2D : %le\n\n",max_conservation_error_2D);
      }
  } 

} // end namespace Mhd 

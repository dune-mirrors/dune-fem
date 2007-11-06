static volatile char RCSId_eqns_mhd_cc [] = "$Id$";

// #include "mhd_eqns.hh"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

namespace Mhd {

  double MhdSolver::GAMMA_ID      = 1.4; //testweise us, statt: 1.4;
  double MhdSolver::R_ID          = 8.314e7;  //testweise us, statt 8.31;
  eos_cons_t *MhdSolver::EOS_CONS = (eos_cons_t *)0;
  eos_prim_t *MhdSolver::EOS_PRIM = (eos_prim_t *)0;
  MhdSolver::mflux_t MhdSolver::FLUXCHOICE = MhdSolver::mf_dw;
  const double MhdSolver::s2 = 1.41421356237309510000 ;

  const int     MhdSolver::ewave = 3;
  double  MhdSolver::maxtime;
  void  (*MhdSolver::initialize)(const double px,
                                 MhdSolver::VEC1D *pu);

  char   MhdSolver::algorithm[20];
  int    MhdSolver::docheck;
  int    MhdSolver::tstep_auto;
  int    MhdSolver::use_hll_fix;
  int    MhdSolver::verbose_mode;
  int    MhdSolver::use_point_fix;
  int    MhdSolver::use_entropy_fix;
  int    MhdSolver::use_velocity_fix;
  int    MhdSolver::check_rmv_errors;
  int    MhdSolver::use_reference_fix;
  int    MhdSolver::check_2D_conservation;
  int    MhdSolver::hll_use_roe_mean_values;
  double MhdSolver::visc;

  int    MhdSolver::undo;
  int    MhdSolver::hllcount;
  int    MhdSolver::point_fixes;
  int    MhdSolver::entropy_fixes;
  int    MhdSolver::velocity_fixes;
  int    MhdSolver::reference_fixes;
  int    MhdSolver::conservation_errors_2D;
  double MhdSolver::max_check_rmv_error;
  double MhdSolver::max_conservation_error_2D;

  const double MhdSolver::MhdEPS = 1e-12;

  void (*MhdSolver::ipflux_1d)(const MhdSolver::VEC1D,
                               const MhdSolver::VEC1D,
			       MhdSolver::VEC1D *,
                               double);
  void (*MhdSolver::init_ipflux)(void);
  flux_t *MhdSolver::rgflux_1d = (flux_t *)0;


  /*******************************************************************************
 ######  #       #    #  #    #
 #       #       #    #   #  #
 #####   #       #    #    ##
 #       #       #    #    ##
 #       #       #    #   #  #
 #       ######   ####   #    #
  *******************************************************************************/
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

  void MhdSolver::init(flux1d_t fluxchoice)
  {
    switch (fluxchoice)
      {
      case MhdSolver::mf_bct:
        ipflux_1d   = flux_bct;
        init_ipflux = initialize_bct;
        break;
      case MhdSolver::mf_dw:
        ipflux_1d   = flux_dw;
        init_ipflux = initialize_dw;
        break;
      case MhdSolver::mf_hlle:
        ipflux_1d   = flux_hlle;
        init_ipflux = initialize_hlle;
        break;
      case MhdSolver::mf_hllemg:
        ipflux_1d   = flux_hllemg;
        init_ipflux = initialize_hllemg;
        break;
      case MhdSolver::mf_hlleml:
        ipflux_1d   = flux_hlleml;
        init_ipflux = initialize_hlleml;
        break;
      case MhdSolver::mf_roe:
        ipflux_1d   = flux_roe;
        init_ipflux = initialize_roe;
        break;
      case MhdSolver::mf_rgdw:
        rgflux_1d   = new flux_dw_t;
        break;
      case MhdSolver::mf_rgdwc:
        rgflux_1d   = new flux_dwc_t;
        break;
      case MhdSolver::mf_rghlle:
        rgflux_1d   = new flux_hlle_t;
        break;
      case MhdSolver::mf_rghllem:
        rgflux_1d   = new flux_hllem_t;
        break;
      case MhdSolver::mf_rghllemc:
        rgflux_1d   = new flux_hllemc_t;
        break;
      case MhdSolver::mf_rglf:
        rgflux_1d   = new flux_lf_t;
        break;
      case MhdSolver::mf_rgvfroe:
        rgflux_1d   = new flux_vfroe_t;
        break;
      default:
        cerr << "ERROR (MhdSolver::init()): unknown 1d-flux!"
             << endl;
        abort();
      }

    if (fluxchoice < MhdSolver::first_rgflux)
      init_ipflux();
  }

  /******************************************************************************
   ******************************************************************************
   ***                                                                        ***
   ***                         CLASS IMPLEMENTATIONS                          ***
   ***                                                                        ***
   ******************************************************************************
   ******************************************************************************/

  /*******************************************************************************
   class Cons_mhd
  *******************************************************************************/

  double Cons_mhd::tau() const
  {
    return (1.0 / rho());
  }

  double Cons_mhd::eps() const
  {
    return (e() - 0.5 * u2() - B2() / (8.0 * M_PI * rho()));
  }

  void Cons_mhd::gradp(vec_t &pret) const
  {
    double ltau = tau();

    double ldpdtau = eos->dpdtau(*this);
    double ldpdeps = eos->dpdeps(*this);

    double ldepsdrho   =   ltau * (0.5 * u2() - eps());
    double ldepsdrhoux = - ltau * ux();
    double ldepsdrhouy = - ltau * uy();
    double ldepsdrhouz = - ltau * uz();
    double ldepsdBx    = - Bx() * ltau / (4.0 * M_PI);
    double ldepsdBy    = - By() * ltau / (4.0 * M_PI);
    double ldepsdBz    = - Bz() * ltau / (4.0 * M_PI);
    double ldepsdrhoe  =   ltau;

    pret[0] = - ltau * ltau * ldpdtau + ldepsdrho * ldpdeps;
    pret[1] = ldepsdrhoux * ldpdeps;
    pret[2] = ldepsdrhouy * ldpdeps;
    pret[3] = ldepsdrhouz * ldpdeps;
    pret[4] = ldepsdBx    * ldpdeps;
    pret[5] = ldepsdBy    * ldpdeps;
    pret[6] = ldepsdBz    * ldpdeps;
    pret[7] = ldepsdrhoe  * ldpdeps;
  }

  void Cons_mhd::dWdU(mat_t &pret) const
  {
    vec_t lgradp;
    double ltau = tau();

    gradp(lgradp);

    pret(0,0) = - ltau * ltau;
    pret(0,1) = 0.0;
    pret(0,2) = 0.0;
    pret(0,3) = 0.0;
    pret(0,4) = 0.0;
    pret(0,5) = 0.0;
    pret(0,6) = 0.0;
    pret(0,7) = 0.0;

    pret(1,0) = - ux() * ltau;
    pret(1,1) = ltau;
    pret(1,2) = 0.0;
    pret(1,3) = 0.0;
    pret(1,4) = 0.0;
    pret(1,5) = 0.0;
    pret(1,6) = 0.0;
    pret(1,7) = 0.0;

    pret(2,0) = - uy() * ltau;
    pret(2,1) = 0.0;
    pret(2,2) = ltau;
    pret(2,3) = 0.0;
    pret(2,4) = 0.0;
    pret(2,5) = 0.0;
    pret(2,6) = 0.0;
    pret(2,7) = 0.0;

    pret(3,0) = - uz() * ltau;
    pret(3,1) = 0.0;
    pret(3,2) = 0.0;
    pret(3,3) = ltau;
    pret(3,4) = 0.0;
    pret(3,5) = 0.0;
    pret(3,6) = 0.0;
    pret(3,7) = 0.0;

    pret(4,0) = 0.0;
    pret(4,1) = 0.0;
    pret(4,2) = 0.0;
    pret(4,3) = 0.0;
    pret(4,4) = 1.0;
    pret(4,5) = 0.0;
    pret(4,6) = 0.0;
    pret(4,7) = 0.0;

    pret(5,0) = 0.0;
    pret(5,1) = 0.0;
    pret(5,2) = 0.0;
    pret(5,3) = 0.0;
    pret(5,4) = 0.0;
    pret(5,5) = 1.0;
    pret(5,6) = 0.0;
    pret(5,7) = 0.0;

    pret(6,0) = 0.0;
    pret(6,1) = 0.0;
    pret(6,2) = 0.0;
    pret(6,3) = 0.0;
    pret(6,4) = 0.0;
    pret(6,5) = 0.0;
    pret(6,6) = 1.0;
    pret(6,7) = 0.0;

    pret(7,0) = lgradp[0];
    pret(7,1) = lgradp[1];
    pret(7,2) = lgradp[2];
    pret(7,3) = lgradp[3];
    pret(7,4) = lgradp[4];
    pret(7,5) = lgradp[5];
    pret(7,6) = lgradp[6];
    pret(7,7) = lgradp[7];
  }

  /*******************************************************************************
   class Prim_rp_mhd
  *******************************************************************************/

  double Prim_rp_mhd::tau() const
  {
    return (1.0 / rho());
  }

  double Prim_rp_mhd::p() const
  {
    assert(val[7] > 0.0);
    return val[7];
  }

  /*******************************************************************************
   class Prim_tp_mhd
  *******************************************************************************/

  double Prim_tp_mhd::tau() const
  {
    assert(val[0] > 0.0);
    return val[0];
  }

  double Prim_tp_mhd::p() const
  {
    assert(val[7] > 0.0);
    return val[7];
  }

  void Prim_tp_mhd::gradcs2(vec_t &pret) const
  {
    pret[0] = eos->dcs2dtau(*this);
    pret[1] = 0.0;
    pret[2] = 0.0;
    pret[3] = 0.0;
    pret[4] = 0.0;
    pret[5] = 0.0;
    pret[6] = 0.0;
    pret[7] = eos->dcs2dp(*this);
  }

  void Prim_tp_mhd::gradeps(vec_t &pret) const
  {
    pret[0] = eos->depsdtau(*this);
    pret[1] = 0.0;
    pret[2] = 0.0;
    pret[3] = 0.0;
    pret[4] = 0.0;
    pret[5] = 0.0;
    pret[6] = 0.0;
    pret[7] = eos->depsdp(*this);
  }

  void Prim_tp_mhd::dUdW(mat_t &pret) const
  {
    vec_t lgradeps,lgradrhoe;
    double lrho = rho();
    double ltau = tau();

    gradeps(lgradeps);

    lgradrhoe[0] = (-(eps() + 0.5 * u2()) / ltau + lgradeps[0]) / ltau;
    lgradrhoe[1] = (lgradeps[1] + ux()) / ltau;
    lgradrhoe[2] = (lgradeps[2] + uy()) / ltau;
    lgradrhoe[3] = (lgradeps[3] + uz()) / ltau;
    lgradrhoe[4] = lgradeps[4] / ltau + Bx() / (4.0 * M_PI);
    lgradrhoe[5] = lgradeps[5] / ltau + By() / (4.0 * M_PI);
    lgradrhoe[6] = lgradeps[6] / ltau + Bz() / (4.0 * M_PI);
    lgradrhoe[7] = lgradeps[7] / ltau;

    pret(0,0) = - lrho * lrho;
    pret(0,1) = 0.0;
    pret(0,2) = 0.0;
    pret(0,3) = 0.0;
    pret(0,4) = 0.0;
    pret(0,5) = 0.0;
    pret(0,6) = 0.0;
    pret(0,7) = 0.0;

    pret(1,0) = - ux() * lrho * lrho;
    pret(1,1) = lrho;
    pret(1,2) = 0.0;
    pret(1,3) = 0.0;
    pret(1,4) = 0.0;
    pret(1,5) = 0.0;
    pret(1,6) = 0.0;
    pret(1,7) = 0.0;

    pret(2,0) = - uy() * lrho * lrho;
    pret(2,1) = 0.0;
    pret(2,2) = lrho;
    pret(2,3) = 0.0;
    pret(2,4) = 0.0;
    pret(2,5) = 0.0;
    pret(2,6) = 0.0;
    pret(2,7) = 0.0;

    pret(3,0) = - uz() * lrho * lrho;
    pret(3,1) = 0.0;
    pret(3,2) = 0.0;
    pret(3,3) = lrho;
    pret(3,4) = 0.0;
    pret(3,5) = 0.0;
    pret(3,6) = 0.0;
    pret(3,7) = 0.0;

    pret(4,0) = 0.0;
    pret(4,1) = 0.0;
    pret(4,2) = 0.0;
    pret(4,3) = 0.0;
    pret(4,4) = 1.0;
    pret(4,5) = 0.0;
    pret(4,6) = 0.0;
    pret(4,7) = 0.0;

    pret(5,0) = 0.0;
    pret(5,1) = 0.0;
    pret(5,2) = 0.0;
    pret(5,3) = 0.0;
    pret(5,4) = 0.0;
    pret(5,5) = 1.0;
    pret(5,6) = 0.0;
    pret(5,7) = 0.0;

    pret(6,0) = 0.0;
    pret(6,1) = 0.0;
    pret(6,2) = 0.0;
    pret(6,3) = 0.0;
    pret(6,4) = 0.0;
    pret(6,5) = 0.0;
    pret(6,6) = 1.0;
    pret(6,7) = 0.0;

    pret(7,0) = lgradrhoe[0];
    pret(7,1) = lgradrhoe[1];
    pret(7,2) = lgradrhoe[2];
    pret(7,3) = lgradrhoe[3];
    pret(7,4) = lgradrhoe[4];
    pret(7,5) = lgradrhoe[5];
    pret(7,6) = lgradrhoe[6];
    pret(7,7) = lgradrhoe[7];
  }

  /*******************************************************************************
   class Mhd
  *******************************************************************************/

  int Mhd::nonlinear(int pidx) const
  {
    int lret = 1;

    if ((pidx == 1) || (pidx == 3) || (pidx == 4) || (pidx == 6))
      lret = 0;

    return lret;
  }

  void Mhd::f(const cons_t &pcons, cons_t &pret) const
  {
    pret[0] = pcons.rhoux();
    pret[1] =   pcons.rhoux() * pcons.ux() + pcons.p()
      + pcons.B2() / (8.0 * M_PI)
      - pcons.Bx() * pcons.Bx() / (4.0 * M_PI);
    pret[2] = pcons.rhoux() * pcons.uy() - pcons.Bx() * pcons.By() / (4.0 * M_PI);
    pret[3] = pcons.rhoux() * pcons.uz() - pcons.Bx() * pcons.Bz() / (4.0 * M_PI);
    pret[4] = 0.0;
    pret[5] = pcons.ux() * pcons.By() - pcons.uy() * pcons.Bx();
    pret[6] = pcons.ux() * pcons.Bz() - pcons.uz() * pcons.Bx();
    pret[7] =   (pcons.rhoe() + pcons.p() + pcons.B2() / (8.0 * M_PI)) * pcons.ux()
      - pcons.Bx() * pcons.Bu() / (4.0 * M_PI);

  }

  void Mhd::f(const prim_t &pprim, cons_t &pret) const
  {
    pret[0] = pprim.rhoux();
    pret[1] =   pprim.rhoux() * pprim.ux() + pprim.p()
      + pprim.B2() / (8.0 * M_PI)
      - pprim.Bx() * pprim.Bx() / (4.0 * M_PI);
    pret[2] = pprim.rhoux() * pprim.uy() - pprim.Bx() * pprim.By() / (4.0 * M_PI);
    pret[3] = pprim.rhoux() * pprim.uz() - pprim.Bx() * pprim.Bz() / (4.0 * M_PI);
    pret[4] = 0.0;
    pret[5] = pprim.ux() * pprim.By() - pprim.uy() * pprim.Bx();
    pret[6] = pprim.ux() * pprim.Bz() - pprim.uz() * pprim.Bx();
    pret[7] =   (pprim.rhoe() + pprim.p() + pprim.B2() / (8.0 * M_PI)) * pprim.ux()
      - pprim.Bx() * pprim.Bu() / (4.0 * M_PI);
  }

  double Mhd::lambda(int pidx, const prim_t &pprim) const
  {
    double lret = pprim.ux();

    switch (pidx)
      {
      case 0:
        lret -= pprim.vf();
        break;
      case 1:
        lret -= pprim.vax();
        break;
      case 2:
        lret -= pprim.vs();
        break;
      case 3:
        break;
      case 4:
        lret = 0.0;
        break;
      case 5:
        lret += pprim.vs();
        break;
      case 6:
        lret += pprim.vax();
        break;
      case 7:
        lret += pprim.vf();
        break;
      default:
        {assert(0>1); abort(); }
      }

    return lret;
  }

  void Mhd::lambda_vec(const prim_t &pprim, vec_t &pret) const
  {
    double lux  = pprim.ux();
    double lvs  = pprim.vs();
    double lvax = pprim.vax();
    double lvf  = pprim.vf();

    pret[0] = lux - lvf;
    pret[1] = lux - lvax;
    pret[2] = lux - lvs;
    pret[3] = lux;
    pret[4] = 0.0;
    pret[5] = lux + lvs;
    pret[6] = lux + lvax;
    pret[7] = lux + lvf;
  }

  void Mhd::r(int pidx, const prim_t &pprim, vec_t &pret) const
  {
    switch (pidx)
      {
      case 0:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lvf   = pprim.vf();
          double lvax  = pprim.vax();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   lalphaf * ltau;
          pret[1] =   lalphaf * lvf;
          pret[2] = - lalphas * lbetay * lvax * lsgnBx;
          pret[3] = - lalphas * lbetaz * lvax * lsgnBx;
          pret[4] =   0.0;
          pret[5] = - lalphas * lbetay * lvf * sqrt(4.0 * M_PI / ltau);
          pret[6] = - lalphas * lbetaz * lvf * sqrt(4.0 * M_PI / ltau);
          pret[7] = - lalphaf * pprim.gamma() * pprim.p(); 

          pret   *= lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                   + lalphas * lalphas * (lvf2 + lvax2) );

          break;
        }
      case 1:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lsqrtBy2pBz2,lbetay,lbetaz,lsgnBx;

          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0; 
          pret[1] =   0.0;
          pret[2] = - lbetaz;
          pret[3] =   lbetay;
          pret[4] =   0.0;
          pret[5] = - lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetaz;
          pret[6] =   lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetay;
          pret[7] =   0.0;

          pret   *= pprim.vf() / sqrt(2.0);

          break;
        }
      case 2:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lcs   = pprim.cs();
          double lvs   = pprim.vs();
          double lvf   = pprim.vf();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   lalphas * ltau; 
          pret[1] =   lalphas * lvs;
          pret[2] =   lalphaf * lbetay * lcs * lsgnBx;
          pret[3] =   lalphaf * lbetaz * lcs * lsgnBx;
          pret[4] =   0.0;
          pret[5] =   lalphaf * lbetay * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
          pret[6] =   lalphaf * lbetaz * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
          pret[7] = - lalphas * pprim.gamma() * pprim.p();

          pret   *= lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                   + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
          break;
        }
      case 3:
        {
          pret[0] = pprim.tau(); 
          pret[1] = 0.0;
          pret[2] = 0.0;
          pret[3] = 0.0;
          pret[4] = 0.0;
          pret[5] = 0.0;
          pret[6] = 0.0;
          pret[7] = 0.0;
          break;
        }
      case 4:
        {
          pret[0] = 0.0; 
          pret[1] = 0.0;
          pret[2] = 0.0;
          pret[3] = 1.0;
          pret[4] = 0.0;
          pret[5] = 0.0;
          pret[6] = 0.0;
          pret[7] = 0.0;
          break;
        }
      case 5:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lcs   = pprim.cs();
          double lvs   = pprim.vs();
          double lvf   = pprim.vf();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] = - lalphas * ltau; 
          pret[1] =   lalphas * lvs;
          pret[2] =   lalphaf * lbetay * lcs * lsgnBx;
          pret[3] =   lalphaf * lbetaz * lcs * lsgnBx;
          pret[4] =   0.0;
          pret[5] = - lalphaf * lbetay * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
          pret[6] = - lalphaf * lbetaz * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
          pret[7] =   lalphas * pprim.gamma() * pprim.p();

          pret   *= lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                   + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
          break;
        }
      case 6:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lsqrtBy2pBz2,lbetay,lbetaz,lsgnBx;

          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0; 
          pret[1] =   0.0;
          pret[2] =   lbetaz;
          pret[3] = - lbetay;
          pret[4] =   0.0;
          pret[5] = - lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetaz;
          pret[6] =   lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetay;
          pret[7] =   0.0;

          pret   *= pprim.vf() / sqrt(2.0);

          break;
        }
      case 7:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lvf   = pprim.vf();
          double lvax  = pprim.vax();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] = - lalphaf * ltau;
          pret[1] =   lalphaf * lvf;
          pret[2] = - lalphas * lbetay * lvax * lsgnBx;
          pret[3] = - lalphas * lbetaz * lvax * lsgnBx;
          pret[4] =   0.0;
          pret[5] =   lalphas * lbetay * lvf * sqrt(4.0 * M_PI / ltau);
          pret[6] =   lalphas * lbetaz * lvf * sqrt(4.0 * M_PI / ltau);
          pret[7] =   lalphaf * pprim.gamma() * pprim.p(); 

          pret   *= lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                   + lalphas * lalphas * (lvf2 + lvax2) );

          break;
        }
      default:
        {assert(0>1); abort();}
      }
  }

  void Mhd::r_mat(const prim_t &pprim, mat_t &pret) const
  {
    double lgp   = pprim.gamma() * pprim.p();
    double ltau  = pprim.tau();
    double lBy   = pprim.By();
    double lBz   = pprim.Bz();
    double lcs   = pprim.cs();
    double lvs   = pprim.vs();
    double lvax  = pprim.vax();
    double lvf   = pprim.vf();
    double lcs2  = pprim.cs2();
    double lvs2  = pprim.vs2();
    double lvax2 = pprim.vax2();
    double lvf2  = pprim.vf2();
    double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
    double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;
    double lmultvs,lmultvax,lmultvf;

    lrad = lvf2 - lvs2;
    if (lrad <= MhdEPS)
      lrad = 0.0;
    lsqrtvf2mvs2 = sqrt( lrad );
    if (lsqrtvf2mvs2 > MhdEPS)
      {
        lrad = lvf2 - lvax2;
        if (lrad > MhdEPS)
          lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
        else lalphaf = 0.0;
        lrad = lvf2 - lcs2;
        if (lrad > MhdEPS)
          lalphas = sqrt(lrad) / lsqrtvf2mvs2;
        else lalphas = 0.0;
      }
    else
      {
        lalphas = 0.0;
        lalphaf = 1.0;
      }
    lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
    if (lsqrtBy2pBz2 > MhdEPS)
      {
        lbetay = lBy / lsqrtBy2pBz2;
        lbetaz = lBz / lsqrtBy2pBz2; 
      }
    else lbetay = lbetaz = 1.0 / sqrt(2.0);
    if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

    lmultvs  = lvf2 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                              + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
    lmultvax = lvf / sqrt(2.0);
    lmultvf  = lvf  / sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                              + lalphas * lalphas * (lvf2 + lvax2) );

    pret(0,0) =   lalphaf * ltau;
    pret(1,0) =   lalphaf * lvf;
    pret(2,0) = - lalphas * lbetay * lvax * lsgnBx;
    pret(3,0) = - lalphas * lbetaz * lvax * lsgnBx;
    pret(4,0) =   0.0;
    pret(5,0) = - lalphas * lbetay * lvf * sqrt(4.0 * M_PI / ltau);
    pret(6,0) = - lalphas * lbetaz * lvf * sqrt(4.0 * M_PI / ltau);
    pret(7,0) = - lalphaf * lgp;

    pret(0,1) =   0.0; 
    pret(1,1) =   0.0;
    pret(2,1) = - lbetaz;
    pret(3,1) =   lbetay;
    pret(4,1) =   0.0;
    pret(5,1) = - lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetaz;
    pret(6,1) =   lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetay;
    pret(7,1) =  0.0;

    pret(0,2) =   lalphas * ltau; 
    pret(1,2) =   lalphas * lvs;
    pret(2,2) =   lalphaf * lbetay * lcs * lsgnBx;
    pret(3,2) =   lalphaf * lbetaz * lcs * lsgnBx;
    pret(4,2) =   0.0;
    pret(5,2) =   lalphaf * lbetay * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
    pret(6,2) =   lalphaf * lbetaz * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
    pret(7,2) = - lalphas * lgp;

    pret(0,3) = ltau; 
    pret(1,3) = 0.0;
    pret(2,3) = 0.0;
    pret(3,3) = 0.0;
    pret(4,3) = 0.0;
    pret(5,3) = 0.0;
    pret(6,3) = 0.0;
    pret(7,3) = 0.0;

    pret(0,4) = 0.0;
    pret(1,4) = 0.0;
    pret(2,4) = 0.0;
    pret(3,4) = 0.0;
    pret(4,4) = 1.0;
    pret(5,4) = 0.0;
    pret(6,4) = 0.0;
    pret(7,4) = 0.0;

    pret(0,5) = - lalphas * ltau; 
    pret(1,5) =   lalphas * lvs;
    pret(2,5) =   lalphaf * lbetay * lcs * lsgnBx;
    pret(3,5) =   lalphaf * lbetaz * lcs * lsgnBx;
    pret(4,5) =   0.0;
    pret(5,5) = - lalphaf * lbetay * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
    pret(6,5) = - lalphaf * lbetaz * lcs2 / lvf * sqrt(4.0 * M_PI / ltau);
    pret(7,5) =   lalphas * lgp;

    pret(0,6) =   0.0; 
    pret(1,6) =   0.0;
    pret(2,6) =   lbetaz;
    pret(3,6) = - lbetay;
    pret(4,6) =   0.0;
    pret(5,6) = - lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetaz;
    pret(6,6) =   lsgnBx * sqrt(4.0 * M_PI / ltau) * lbetay;
    pret(7,6) =   0.0;

    pret(0,7) = - lalphaf * ltau;
    pret(1,7) =   lalphaf * lvf;
    pret(2,7) = - lalphas * lbetay * lvax * lsgnBx;
    pret(3,7) = - lalphas * lbetaz * lvax * lsgnBx;
    pret(4,7) =   0.0;
    pret(5,7) =   lalphas * lbetay * lvf * sqrt(4.0 * M_PI / ltau);
    pret(6,7) =   lalphas * lbetaz * lvf * sqrt(4.0 * M_PI / ltau);
    pret(7,7) =   lalphaf * lgp; 

    for (int i=0;i<8;i++)
      {
        pret(i,0) *= lmultvf;
        pret(i,1) *= lmultvax;
        pret(i,2) *= lmultvs;
        pret(i,5) *= lmultvs;
        pret(i,6) *= lmultvax;
        pret(i,7) *= lmultvf;
      }
  }

  void Mhd::l(int pidx, const prim_t &pprim, vec_t &pret) const
  {
    switch (pidx)
      {
      case 0:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lvf   = pprim.vf();
          double lvax  = pprim.vax();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0; 
          pret[1] =   lalphaf * lvf;
          pret[2] = - lalphas * lbetay * lvax * lsgnBx;
          pret[3] = - lalphas * lbetaz * lvax * lsgnBx;
          pret[4] =   0.0;
          pret[5] = - lalphas * lbetay * lvf * sqrt(ltau / (4.0 * M_PI));
          pret[6] = - lalphas * lbetaz * lvf * sqrt(ltau / (4.0 * M_PI));
          pret[7] = - lalphaf * ltau;

          pret   *= 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                     + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

          break;
        }
      case 1:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lsqrtBy2pBz2,lbetay,lbetaz,lsgnBx;

          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0;
          pret[1] =   0.0;
          pret[2] = - lbetaz;
          pret[3] =   lbetay;
          pret[4] =   0.0;
          pret[5] = - lsgnBx * lbetaz * sqrt(ltau / (4.0 * M_PI));
          pret[6] =   lsgnBx * lbetay * sqrt(ltau / (4.0 * M_PI));
          pret[7] =   0.0;

          pret   *= 1.0 / (sqrt(2.0) * pprim.vf());

          break;
        }
      case 2:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lcs   = pprim.cs();
          double lvs   = pprim.vs();
          double lvf   = pprim.vf();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0;
          pret[1] =   lalphas * lvs;
          pret[2] =   lalphaf * lbetay * lcs * lsgnBx;
          pret[3] =   lalphaf * lbetaz * lcs * lsgnBx;
          pret[4] =   0.0;
          pret[5] =   lalphaf * lbetay * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
          pret[6] =   lalphaf * lbetaz * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
          pret[7] = - lalphas * ltau;

          pret   *= 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                  + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

          break;
        }
      case 3:
        {
          pret[0] =   1.0 / pprim.tau(); 
          pret[1] =   0.0;
          pret[2] =   0.0;
          pret[3] =   0.0;
          pret[4] =   0.0;
          pret[5] =   0.0;
          pret[6] =   0.0;
          pret[7] =   1.0 / (pprim.gamma() * pprim.p());
          break;
        }
      case 4:
        {
          pret[0] =   0.0;
          pret[1] =   0.0;
          pret[2] =   0.0;
          pret[3] =   0.0;
          pret[4] =   1.0;
          pret[5] =   0.0;
          pret[6] =   0.0;
          pret[7] =   0.0;
          break;
        }
      case 5:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lcs   = pprim.cs();
          double lvs   = pprim.vs();
          double lvf   = pprim.vf();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0;
          pret[1] =   lalphas * lvs;
          pret[2] =   lalphaf * lbetay * lcs * lsgnBx;
          pret[3] =   lalphaf * lbetaz * lcs * lsgnBx;
          pret[4] =   0.0;
          pret[5] = - lalphaf * lbetay * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
          pret[6] = - lalphaf * lbetaz * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
          pret[7] =   lalphas * ltau;

          pret   *= 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                                  + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );

          break;
        }
      case 6:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lsqrtBy2pBz2,lbetay,lbetaz,lsgnBx;

          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0;
          pret[1] =   0.0;
          pret[2] =   lbetaz;
          pret[3] = - lbetay;
          pret[4] =   0.0;
          pret[5] = - lsgnBx * lbetaz * sqrt(ltau / (4.0 * M_PI));
          pret[6] =   lsgnBx * lbetay * sqrt(ltau / (4.0 * M_PI));
          pret[7] =   0.0;

          pret   *= 1.0 / (sqrt(2.0) * pprim.vf());

          break;
        }
      case 7:
        {
          double ltau  = pprim.tau();
          double lBy   = pprim.By();
          double lBz   = pprim.Bz();
          double lvf   = pprim.vf();
          double lvax  = pprim.vax();
          double lcs2  = pprim.cs2();
          double lvs2  = pprim.vs2();
          double lvax2 = pprim.vax2();
          double lvf2  = pprim.vf2();
          double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
          double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;

          lrad = lvf2 - lvs2;
          if (lrad <= MhdEPS)
            lrad = 0.0;
          lsqrtvf2mvs2 = sqrt( lrad );
          if (lsqrtvf2mvs2 > MhdEPS)
            {
              lrad = lvf2 - lvax2;
              if (lrad > MhdEPS)
                lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphaf = 0.0;
              lrad = lvf2 - lcs2;
              if (lrad > MhdEPS)
                lalphas = sqrt(lrad) / lsqrtvf2mvs2;
              else lalphas = 0.0;
            }
          else
            {
              lalphas = 0.0;
              lalphaf = 1.0;
            }
          lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
          if (lsqrtBy2pBz2 > MhdEPS)
            {
              lbetay = lBy / lsqrtBy2pBz2;
              lbetaz = lBz / lsqrtBy2pBz2; 
            }
          else lbetay = lbetaz = 1.0 / sqrt(2.0);
          if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

          pret[0] =   0.0; 
          pret[1] =   lalphaf * lvf;
          pret[2] = - lalphas * lbetay * lvax * lsgnBx;
          pret[3] = - lalphas * lbetaz * lvax * lsgnBx;
          pret[4] =   0.0;
          pret[5] =   lalphas * lbetay * lvf * sqrt(ltau / (4.0 * M_PI));
          pret[6] =   lalphas * lbetaz * lvf * sqrt(ltau / (4.0 * M_PI));
          pret[7] =   lalphaf * ltau;

          pret   *= 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                     + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

          break;
        }
      default:
        {assert(0>1); abort();}
      }
  }

  void Mhd::l_mat(const prim_t &pprim, mat_t &pret) const
  {
    double ltau  = pprim.tau();
    double lBy   = pprim.By();
    double lBz   = pprim.Bz();
    double lcs   = pprim.cs();
    double lvs   = pprim.vs();
    double lvax  = pprim.vax();
    double lvf   = pprim.vf();
    double lcs2  = pprim.cs2();
    double lvs2  = pprim.vs2();
    double lvax2 = pprim.vax2();
    double lvf2  = pprim.vf2();
    double lrad,lsqrtvf2mvs2,lsqrtBy2pBz2;
    double lalphaf,lalphas,lbetay,lbetaz,lsgnBx;
    double lmultvs,lmultvax,lmultvf;

    lrad = lvf2 - lvs2;
    if (lrad <= MhdEPS)
      lrad = 0.0;
    lsqrtvf2mvs2 = sqrt( lrad );
    if (lsqrtvf2mvs2 > MhdEPS)
      {
        lrad = lvf2 - lvax2;
        if (lrad > MhdEPS)
          lalphaf = sqrt(lrad) / lsqrtvf2mvs2;
        else lalphaf = 0.0;
        lrad = lvf2 - lcs2;
        if (lrad > MhdEPS)
          lalphas = sqrt(lrad) / lsqrtvf2mvs2;
        else lalphas = 0.0;
      }
    else
      {
        lalphas = 0.0;
        lalphaf = 1.0;
      }
    lsqrtBy2pBz2 = sqrt(lBy * lBy + lBz * lBz);
    if (lsqrtBy2pBz2 > MhdEPS)
      {
        lbetay = lBy / lsqrtBy2pBz2;
        lbetaz = lBz / lsqrtBy2pBz2; 
      }
    else lbetay = lbetaz = 1.0 / sqrt(2.0);
    if (pprim.Bx() >= 0.0) lsgnBx = 1.0; else lsgnBx = -1.0;

    lmultvs  = 1.0 / sqrt(   lalphaf * lalphaf * lcs2 * (lvf2 + lcs2)
                             + lalphas * lalphas * lvf2 * (lvs2 + lcs2) );
    lmultvax = 1.0 / ( sqrt(2.0) * lvf );
    lmultvf  = 1.0  / ( sqrt(   lalphaf * lalphaf * (lvf2 + lcs2)
                                + lalphas * lalphas * (lvf2 + lvax2) ) * lvf );

    pret(0,0) =   0.0; 
    pret(0,1) =   lalphaf * lvf;
    pret(0,2) = - lalphas * lbetay * lvax * lsgnBx;
    pret(0,3) = - lalphas * lbetaz * lvax * lsgnBx;
    pret(0,4) =   0.0;
    pret(0,5) = - lalphas * lbetay * lvf * sqrt(ltau / (4.0 * M_PI));
    pret(0,6) = - lalphas * lbetaz * lvf * sqrt(ltau / (4.0 * M_PI));
    pret(0,7) = - lalphaf * ltau;

    pret(1,0) =   0.0;
    pret(1,1) =   0.0;
    pret(1,2) = - lbetaz;
    pret(1,3) =   lbetay;
    pret(1,4) =   0.0;
    pret(1,5) = - lsgnBx * lbetaz * sqrt(ltau / (4.0 * M_PI));
    pret(1,6) =   lsgnBx * lbetay * sqrt(ltau / (4.0 * M_PI));
    pret(1,7) =   0.0;

    pret(2,0) =   0.0;
    pret(2,1) =   lalphas * lvs;
    pret(2,2) =   lalphaf * lbetay * lcs * lsgnBx;
    pret(2,3) =   lalphaf * lbetaz * lcs * lsgnBx;
    pret(2,4) =   0.0;
    pret(2,5) =   lalphaf * lbetay * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
    pret(2,6) =   lalphaf * lbetaz * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
    pret(2,7) = - lalphas * ltau;

    pret(3,0) =   1.0 / pprim.tau(); 
    pret(3,1) =   0.0;
    pret(3,2) =   0.0;
    pret(3,3) =   0.0;
    pret(3,4) =   0.0;
    pret(3,5) =   0.0;
    pret(3,6) =   0.0;
    pret(3,7) =   1.0 / (pprim.gamma() * pprim.p());

    pret(4,0) =   0.0;
    pret(4,1) =   0.0;
    pret(4,2) =   0.0;
    pret(4,3) =   0.0;
    pret(4,4) =   1.0;
    pret(4,5) =   0.0;
    pret(4,6) =   0.0;
    pret(4,7) =   0.0;

    pret(5,0) =   0.0;
    pret(5,1) =   lalphas * lvs;
    pret(5,2) =   lalphaf * lbetay * lcs * lsgnBx;
    pret(5,3) =   lalphaf * lbetaz * lcs * lsgnBx;
    pret(5,4) =   0.0;
    pret(5,5) = - lalphaf * lbetay * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
    pret(5,6) = - lalphaf * lbetaz * lcs2 / lvf * sqrt(ltau / (4.0 * M_PI));
    pret(5,7) =   lalphas * ltau;

    pret(6,0) =   0.0;
    pret(6,1) =   0.0;
    pret(6,2) =   lbetaz;
    pret(6,3) = - lbetay;
    pret(6,4) =   0.0;
    pret(6,5) = - lsgnBx * lbetaz * sqrt(ltau / (4.0 * M_PI));
    pret(6,6) =   lsgnBx * lbetay * sqrt(ltau / (4.0 * M_PI));
    pret(6,7) =   0.0;

    pret(7,0) =   0.0; 
    pret(7,1) =   lalphaf * lvf;
    pret(7,2) = - lalphas * lbetay * lvax * lsgnBx;
    pret(7,3) = - lalphas * lbetaz * lvax * lsgnBx;
    pret(7,4) =   0.0;
    pret(7,5) =   lalphas * lbetay * lvf * sqrt(ltau / (4.0 * M_PI));
    pret(7,6) =   lalphas * lbetaz * lvf * sqrt(ltau / (4.0 * M_PI));
    pret(7,7) =   lalphaf * ltau;

    for (int i=0;i<8;i++)
      {
        pret(0,i) *= lmultvf;
        pret(1,i) *= lmultvax;
        pret(2,i) *= lmultvs;
        pret(5,i) *= lmultvs;
        pret(6,i) *= lmultvax;
        pret(7,i) *= lmultvf;
      }
  }

  double Mhd::dt_local(double phl, double phr,
                       const cons_t &pconsl, const cons_t &pconsr) const
  {
    const double lmaxstep = 100.0;

    double lrho   = min(pconsl.rho(),pconsr.rho());
    double lux    = max(fabs(pconsl.ux()),fabs(pconsr.ux()));
    double lBx    = max(fabs(pconsl.Bx()),fabs(pconsr.Bx()));
    double lBy    = max(fabs(pconsl.By()),fabs(pconsr.By()));
    double lBz    = max(fabs(pconsl.Bz()),fabs(pconsr.Bz()));
    double lp     = max(pconsl.p(),pconsr.p());
    double lgamma = max(pconsl.gamma(),pconsr.gamma());  
    double lret,lval,lrad;

    assert(phl > 0.0);
    assert(phr > 0.0);

    lrad = ((lBx * lBx + lBy * lBy + lBz * lBz) / (4.0 * M_PI) + lgamma * lp) / lrho;

    if (lrad < MhdEPS)
      lval = lux;
    else
      lval = lux + sqrt(lrad);
  
    if (fabs(lval) < MhdEPS)
      {
        lret = lmaxstep;
      }
    else
      {
        lret = min(phl,phr) / lval;

        if (lret > lmaxstep)
          lret = lmaxstep;
      }

    return lret;
  }

  double Mhd::indicator(bool prel, int pindselect,
                        const cons_t &pconsl, const cons_t &pconsr) const
  {
    double ret;

    switch (pindselect)
      {
      case mi_rho:
        {
          double rhol = pconsl.rho();
          double rhor = pconsr.rho();
          double rel  = ((prel) ? (0.5*(rhol+rhor)) : 1.0);

          ret = (rhor - rhol) / rel;
        }
        break;
      case mi_rhoux:
        {
          double rhouxl = pconsl.rhoux();
          double rhouxr = pconsr.rhoux();
          double rel    = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.rhoux())
                        + fabs(mean.rhouy())
                        + fabs(mean.rhouz()));
            }

          ret = (rhouxr - rhouxl) / rel;
        }
        break;
      case mi_rhouy:
        {
          double rhouyl = pconsl.rhouy();
          double rhouyr = pconsr.rhouy();
          double rel    = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.rhoux())
                        + fabs(mean.rhouy())
                        + fabs(mean.rhouz()));
            }

          ret = (rhouyr - rhouyl) / rel;
        }
        break;
      case mi_rhouz:
        {
          double rhouzl = pconsl.rhouz();
          double rhouzr = pconsr.rhouz();
          double rel    = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.rhoux())
                        + fabs(mean.rhouy())
                        + fabs(mean.rhouz()));
            }

          ret = (rhouzr - rhouzl) / rel;
        }
        break;
      case mi_bx:
        {
          double Bxl = pconsl.Bx();
          double Bxr = pconsr.Bx();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.Bx())
                        + fabs(mean.By())
                        + fabs(mean.Bz()));
            }

          ret = (Bxr - Bxl) / rel;
        }
        break;
      case mi_by:
        {
          double Byl = pconsl.By();
          double Byr = pconsr.By();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.Bx())
                        + fabs(mean.By())
                        + fabs(mean.Bz()));
            }

          ret = (Byr - Byl) / rel;
        }
        break;
      case mi_bz:
        {
          double Bzl = pconsl.Bz();
          double Bzr = pconsr.Bz();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.Bx())
                        + fabs(mean.By())
                        + fabs(mean.Bz()));
            }

          ret = (Bzr - Bzl) / rel;
        }
        break;
      case mi_rhoe:
        {
          double rhoel = pconsl.rhoe();
          double rhoer = pconsr.rhoe();
          double rel   = ((prel) ? (0.5*(rhoel+rhoer)) : 1.0);

          ret = (rhoer - rhoel) / rel;
        }
        break;
      case mi_ux:
        {
          double uxl = pconsl.ux();
          double uxr = pconsr.ux();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.ux())
                        + fabs(mean.uy())
                        + fabs(mean.uz()));
            }

          ret = (uxr - uxl) / rel;
        }
        break;
      case mi_uy:
        {
          double uyl = pconsl.uy();
          double uyr = pconsr.uy();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.ux())
                        + fabs(mean.uy())
                        + fabs(mean.uz()));
            }

          ret = (uyr - uyl) / rel;
        }
        break;
      case mi_uz:
        {
          double uzl = pconsl.uz();
          double uzr = pconsr.uz();
          double rel = 1.0;

          if (prel)
            {
              cons_t mean((pconsl+pconsr)*0.5);
              rel = max(MhdEPS,
                        fabs(mean.ux())
                        + fabs(mean.uy())
                        + fabs(mean.uz()));
            }

          ret = (uzr - uzl) / rel;
        }
        break;
      case mi_e:
        {
          double el  = pconsl.e();
          double er  = pconsr.e();
          double rel = ((prel) ? (0.5*(el+er)) : 1.0);

          ret = (er - el) / rel;
        }
        break;
      case mi_p:
        {
          double pl  = pconsl.p();
          double pr  = pconsr.p();
          double rel = ((prel) ? (0.5*(pl+pr)) :1.0);

          ret = (pr - pl) / rel;
        }
        break;
      default:
        cerr << "ERROR (Mhd::indicator()): illegal choice of component!" << endl;
        abort();
      }

    return ret;
  }

  void Mhd::Df(const cons_t &pvar, mat_t &pret) const
  {
    vec_t lgradp;
    double lp    = pvar.p();
    double lux   = pvar.ux();
    double luy   = pvar.uy();
    double luz   = pvar.uz();
    double lBx   = pvar.Bx();
    double lBy   = pvar.By();
    double lBz   = pvar.Bz();
    double lB2   = pvar.B2();
    double ltau  = pvar.tau();
    double lrhoe = pvar.rhoe();

    pvar.gradp(lgradp);

    pret(0,0) = 0.0;
    pret(0,1) = 1.0;
    pret(0,2) = 0.0;
    pret(0,3) = 0.0;
    pret(0,4) = 0.0;
    pret(0,5) = 0.0;
    pret(0,6) = 0.0;
    pret(0,7) = 0.0;

    pret(1,0) = - lux * lux + lgradp[0];
    pret(1,1) =   2.0 * lux + lgradp[1];
    pret(1,2) =   lgradp[2];
    pret(1,3) =   lgradp[3];
    pret(1,4) =   0.0;
    pret(1,5) =   lBy / (4.0 * M_PI) + lgradp[5];
    pret(1,6) =   lBz / (4.0 * M_PI) + lgradp[6];
    pret(1,7) =   lgradp[7];

    pret(2,0) = - lux * luy;
    pret(2,1) =   luy;
    pret(2,2) =   lux;
    pret(2,3) =   0.0;
    pret(2,4) =   0.0;
    pret(2,5) = - lBx / (4.0 * M_PI);
    pret(2,6) =   0.0;
    pret(2,7) =   0.0;

    pret(3,0) = - lux * luz;
    pret(3,1) =   luz;
    pret(3,2) =   0.0;
    pret(3,3) =   lux;
    pret(3,4) =   0.0;
    pret(3,5) =   0.0;
    pret(3,6) = - lBx / (4.0 * M_PI);
    pret(3,7) =   0.0;

    pret(4,0) =   0.0;
    pret(4,1) =   0.0;
    pret(4,2) =   0.0;
    pret(4,3) =   0.0;
    pret(4,4) =   0.0;
    pret(4,5) =   0.0;
    pret(4,6) =   0.0;
    pret(4,7) =   0.0;

    pret(5,0) =   (luy * lBx - lux * lBy) * ltau;
    pret(5,1) =   lBy * ltau;
    pret(5,2) = - lBx * ltau;
    pret(5,3) =   0.0;
    pret(5,4) =   0.0;
    pret(5,5) =   lux;
    pret(5,6) =   0.0;
    pret(5,7) =   0.0;

    pret(6,0) =   (luz * lBx - lux * lBz) * ltau;
    pret(6,1) =   lBz * ltau;
    pret(6,2) =   0.0;
    pret(6,3) = - lBx * ltau;
    pret(6,4) =   0.0;
    pret(6,5) =   0.0;
    pret(6,6) =   lux;
    pret(6,7) =   0.0;

    pret(7,0) = - lux * ltau * (lrhoe + lp + lB2 / (8.0 * M_PI))
      + lux * lgradp[0]
      + lBx * pvar.Bu() * ltau / (4.0 * M_PI);
    pret(7,1) =   ltau * (lrhoe + lp + lB2 / (8.0 * M_PI))
      + lux * lgradp[1]
      - lBx * lBx * ltau / (4.0 * M_PI);
    pret(7,2) = - lBx * lBy * ltau / (4.0 * M_PI)
      + lux * lgradp[2];
    pret(7,3) = - lBx * lBz * ltau / (4.0 * M_PI)
      + lux * lgradp[3];
    pret(7,4) =   0.0;
    pret(7,5) =   lux * (lBy / (4.0 * M_PI) + lgradp[5])
      - lBx * luy / (4.0 * M_PI);
    pret(7,6) =   lux * (lBz / (4.0 * M_PI) + lgradp[6])
      - lBx * luz / (4.0 * M_PI);
    pret(7,7) =   lux + lux * lgradp[7];
  }

  void Mhd::grad_lambda_vec(const cons_t &pvar, mat_t &pret) const
  {
    mat_t ldWdU_t;
    vec_t lgradux,lgradcs2prim,lgradcs2,lgradvax2,lgradva2,lgradvf2,lgradvs2;
    prim_t lprim(pvar);
    double ltau   = pvar.tau();
    double lux    = pvar.ux();
    double lBx    = pvar.Bx();
    double lBy    = pvar.By();
    double lBz    = pvar.Bz();
    double lB2    = pvar.B2();
    double lcs2   = pvar.cs2();
    double lva2   = pvar.va2();
    double lvax2  = pvar.vax2();
    double lvf    = pvar.vf();
    double lvs    = pvar.vs();
    double lvax   = pvar.vax();
    double lvsqrt = pvar.vsqrt();
    int i;

    pvar.dWdU_t(ldWdU_t);
    lprim.gradcs2(lgradcs2prim);

    lgradcs2 = ldWdU_t * lgradcs2prim;

    lgradux[0]   = - lux * ltau;
    lgradux[1]   =   ltau;
    lgradux[2]   =   0.0;
    lgradux[3]   =   0.0;
    lgradux[4]   =   0.0;
    lgradux[5]   =   0.0;
    lgradux[6]   =   0.0;
    lgradux[7]   =   0.0;

    lgradvax2[0] = - lBx * lBx * ltau * ltau / (4.0 * M_PI);
    lgradvax2[1] =   0.0;
    lgradvax2[2] =   0.0;
    lgradvax2[3] =   0.0;
    lgradvax2[4] =   lBx * ltau / (2.0 * M_PI);
    lgradvax2[5] =   0.0;
    lgradvax2[6] =   0.0;
    lgradvax2[7] =   0.0;

    lgradva2[0]  = - lB2 * ltau * ltau / (4.0 * M_PI);
    lgradva2[1]  =   0.0;
    lgradva2[2]  =   0.0;
    lgradva2[3]  =   0.0;
    lgradva2[4]  =   lBx * ltau / (2.0 * M_PI);
    lgradva2[5]  =   lBy * ltau / (2.0 * M_PI);
    lgradva2[6]  =   lBz * ltau / (2.0 * M_PI);
    lgradva2[7]  =   0.0;

    for (i=0;i<8;i++)
      {
        lgradvf2[i] = (   lgradva2[i]  * (1.0 + (lva2 + lcs2) / lvsqrt)
                          + lgradcs2[i]  * (1.0 + (lva2 + lcs2 - 2.0 * lvax2) / lvsqrt)
                          - lgradvax2[i] * 2.0 * lcs2 / lvsqrt)
          * 0.5;
        lgradvs2[i] = (   lgradva2[i]  * (1.0 - (lva2 + lcs2) / lvsqrt)
                          + lgradcs2[i]  * (1.0 - (lva2 + lcs2 - 2.0 * lvax2) / lvsqrt)
                          + lgradvax2[i] * 2.0 * lcs2 / lvsqrt)
          * 0.5;
      }

    for (i=0;i<8;i++)
      {
        pret(0,i) = lgradux[i] - lgradvf2[i]  / (2.0 * lvf);
        pret(1,i) = lgradux[i] - lgradvax2[i] / (2.0 * lvax);
        pret(2,i) = lgradux[i] - lgradvs2[i]  / (2.0 * lvs);
        pret(3,i) = lgradux[i];
        pret(4,i) = 0.0;
        pret(5,i) = lgradux[i] + lgradvs2[i]  / (2.0 * lvs);
        pret(6,i) = lgradux[i] + lgradvax2[i] / (2.0 * lvax);
        pret(7,i) = lgradux[i] + lgradvf2[i]  / (2.0 * lvf);
      }
  }

  /****************************************************************************
 ######   ####    ####
 #       #    #  #
 #####   #    #   ####
 #       #    #       #
 #       #    #  #    #
 ######   ####    ####
  *****************************************************************************/

  /******************************************************************************
   ******************************************************************************
   ***                                                                        ***
   ***                         CLASS IMPLEMENTATIONS                          ***
   ***                                                                        ***
   ******************************************************************************
   ******************************************************************************/

  /*******************************************************************************
   class Eos_tabular_tau
  *******************************************************************************/

  Eos_tabular_tau :: Eos_tabular_tau(eos_cons_t *peos_c,eos_prim_t *peos_p,
                                     double ptaumin,double ptaumax,
                                     double pepsmin,double pepsmax,
                                     int pdimtau,int pdimeps) 
    : taumin(ptaumin), taumax(ptaumax), epsmin(pepsmin), epsmax(pepsmax), 
      dimtau(pdimtau), dimeps(pdimeps)
  {
    assert(peos_c && peos_p);
    assert(taumin<taumax && epsmin<epsmax);
    assert(dimtau>0 && dimeps>0);
    int t,e;
    double tau,eps;
    htau=(taumax-taumin)/((double)dimtau-1.0);
    heps=(epsmax-epsmin)/((double)dimeps-1.0);
    _ptable=new double*[dimtau];
    for (t=0,tau=taumin;t<dimtau;t++,tau+=htau)
      {
        _ptable[t]=new double[dimeps];
        for (e=0,eps=epsmin;e<dimeps;e++,eps+=heps)
          {
            cons_t u(peos_c,peos_p); 
            u.rho()=1/tau;
            u.rhoe()=eps/tau;
            _ptable[t][e]=peos_c->p(u);
          }
      }
  }

  inline void Eos_tabular_tau :: findintable(double tau,double eps,
                                             int &t,int &e,
                                             double &psw,double &pse,double &pnw,double &pne,
                                             double &lambdat,double &lambdae)
    //        6  |:   1   :|  5
    //           |:.......:|
    //        ---|---------|---
    //        '''|:        |:''
    //        4  |:   0    |: 2 
    //        ...|:........|:..
    //        ---|---------|---
    //           |:       :|
    //        7  |:   3   :|  8
  {
    t=(int)((tau-taumin)/htau);
    e=(int)((eps-epsmin)/heps);
    if (t<0 || t>dimtau-1 || e<0 || e>dimeps-1)
      cerr << "Ausserhalb der Tabelle! " << endl 
           << tau << " " << eps << " " << t << " " << e << endl;
    if (t<0) // Bereiche 7,3,8
      {
        t=0;
        if (e<0) // Bereich 7
          e=0;
        else if (e>=dimeps-1) // Bereich 8
          e=dimeps-2;
      }
    else if (t>=dimtau-1) // Bereiche 6,1,5
      {
        t=dimtau-2;
        if (e<0) // Bereich 6
          e=0;
        else if (e>=dimeps-1) // Bereich 5
          e=dimeps-2;
      }
    else if (e<0) // Bereich 4
      e=0;
    else if (e>=dimeps-1) // Bereich 2
      e=dimeps-2;
    lambdat=(tau-(double)t*htau-taumin)/htau;
    lambdae=(eps-(double)e*heps-epsmin)/heps;
    psw=ptable(t,e),
      pse=ptable(t+1,e),
      pnw=ptable(t,e+1),
      pne=ptable(t+1,e+1);
  }

  double Eos_tabular_tau :: p(double tau,double eps)
  {
    double p=-1;
    double psw,pse,pnw,pne,lambdat,lambdae;
    int t,e;
    findintable(tau,eps,t,e,psw,pse,pnw,pne,lambdat,lambdae);
    p=(pne+psw-pse-pnw)*lambdat*lambdae+
      (pse-psw)*lambdat+ 
      (pnw-psw)*lambdae+
      psw;
    return p;
  }

  double Eos_tabular_tau :: eps(double tau,double pres)
  {
    double eps=0.0;
    int t=(int)((tau-taumin)/htau),e=0;
    if (t<0)
      {
        t=0;
        while (e<dimeps-2 && p(tau,e*heps+epsmin)<pres)
          e++;
      }
    else if (t>=dimtau-1)
      {
        t=dimtau-2;
        while (e<dimeps-2 && p(tau,e*heps+epsmin)<pres)
          e++;
      }
    else
      {
        while (e<dimeps-2 && ptable(t,e)<pres)
          e++;
      }
    double lambdat=(tau-(double)t*htau-taumin)/htau,
      lambdae;
    double psw=ptable(t,e),
      pse=ptable(t+1,e),
      pnw=ptable(t,e+1),
      pne=ptable(t+1,e+1);
    lambdae=(pres-(pse-psw)*lambdat-psw)/((pne+psw-pse-pnw)*lambdat+(pnw-psw));
    eps=lambdae*heps+(double)e*heps+epsmin;
    assert(eps>0.0);
    return eps;
  }

  double Eos_tabular_tau :: dp_deps(double tau,double eps)
  {
    double dp;
    double psw,pse,pnw,pne,lambdat,lambdae;
    int t,e;
    findintable(tau,eps,t,e,psw,pse,pnw,pne,lambdat,lambdae);
    dp=(pne+psw-pse-pnw)*lambdat+
      (pnw-psw);
    assert(dp>0.0);
    return dp/heps;
  }

  double Eos_tabular_tau :: dp_dtau(double tau,double eps)
  {
    double dp;
    double psw,pse,pnw,pne,lambdat,lambdae;
    int t,e;
    findintable(tau,eps,t,e,psw,pse,pnw,pne,lambdat,lambdae);
    dp=(pne+psw-pse-pnw)*lambdae+
      (pse-psw);
    return dp/htau;
  }

  double Eos_tabular_tau :: cs2(double tau,double eps)
  {
    //double lc2=tau*tau*(-dp_dtau(tau,eps)+p(tau,eps)*dp_deps(tau,eps));
    double dpe,dpt,p;
    double psw,pse,pnw,pne,lambdat,lambdae;
    int t,e;
    findintable(tau,eps,t,e,psw,pse,pnw,pne,lambdat,lambdae);
    p=(pne+psw-pse-pnw)*lambdat*lambdae+
      (pse-psw)*lambdat+ 
      (pnw-psw)*lambdae+
      psw;
    dpe=(pne+psw-pse-pnw)*lambdat+
      (pnw-psw);
    assert(dpe>0.0);
    dpe/=heps;
    dpt=(pne+psw-pse-pnw)*lambdae+
      (pse-psw);
    dpt/=htau;
    double lc2=tau*tau*(-dpt+p*dpe);
    assert(lc2>0);
    return lc2;
  }

  /*******************************************************************************
   class Eos_tabular_rho
  *******************************************************************************/

  Eos_tabular_rho :: Eos_tabular_rho(eos_cons_t *peos_c,eos_prim_t *peos_p,
                                     double ptaumin,double ptaumax,
                                     double pepsmin,double pepsmax,
                                     int pdimtau,int pdimeps) 
    : rhomin(1.0/ptaumax), rhomax(1.0/ptaumin), 
      epsmin(pepsmin), epsmax(pepsmax), 
      dimrho(pdimtau), dimeps(pdimeps)
  {
    assert(peos_c && peos_p);
    assert(rhomin<rhomax && epsmin<epsmax);
    assert(dimrho>0 && dimeps>0);
    int r,e;
    double rho,eps;
    hrho=(rhomax-rhomin)/((double)dimrho-1.0);
    heps=(epsmax-epsmin)/((double)dimeps-1.0);
    _ptable=new double*[dimrho];
    for (r=0,rho=rhomin;r<dimrho;r++,rho+=hrho)
      {
        _ptable[r]=new double[dimeps];
        for (e=0,eps=epsmin;e<dimeps;e++,eps+=heps)
          {
            cons_t u(peos_c,peos_p); 
            u.rho()=rho;
            u.rhoe()=eps*rho;
            _ptable[r][e]=peos_c->p(u);
          }
      }
  }

  inline void Eos_tabular_rho :: findintable(double rho,double eps,
                                             int &r,int &e,
                                             double &psw,double &pse,double &pnw,double &pne,
                                             double &lambdar,double &lambdae)
    //        6  |:   1   :|  5
    //           |:.......:|
    //        ---|---------|---
    //        '''|:        |:''
    //        4  |:   0    |: 2 
    //        ...|:........|:..
    //        ---|---------|---
    //           |:       :|
    //        7  |:   3   :|  8
  {
    r=(int)((rho-rhomin)/hrho);
    e=(int)((eps-epsmin)/heps);
    if (r<0 || r>dimrho-1 || e<0 || e>dimeps-1)
      cerr << "Ausserhalb der Tabelle! " << endl 
           << rho << " " << eps << " " << r << " " << e << endl;
    if (r<0) // Bereiche 7,3,8
      {
        r=0;
        if (e<0) // Bereich 7
          e=0;
        else if (e>=dimeps-1) // Bereich 8
          e=dimeps-2;
      }
    else if (r>=dimrho-1) // Bereiche 6,1,5
      {
        r=dimrho-2;
        if (e<0) // Bereich 6
          e=0;
        else if (e>=dimeps-1) // Bereich 5
          e=dimeps-2;
      }
    else if (e<0) // Bereich 4
      e=0;
    else if (e>=dimeps-1) // Bereich 2
      e=dimeps-2;
    lambdar=(rho-(double)r*hrho-rhomin)/hrho;
    lambdae=(eps-(double)e*heps-epsmin)/heps;
    psw=ptable(r,e),
      pse=ptable(r+1,e),
      pnw=ptable(r,e+1),
      pne=ptable(r+1,e+1);
  }

  double Eos_tabular_rho :: p(double tau,double eps)
  {
    double p=-1;
    double psw,pse,pnw,pne,lambdar,lambdae;
    int r,e;
    findintable(1.0/tau,eps,r,e,psw,pse,pnw,pne,lambdar,lambdae);
    p=(pne+psw-pse-pnw)*lambdar*lambdae+
      (pse-psw)*lambdar+ 
      (pnw-psw)*lambdae+
      psw;
    return p;
  }

  double Eos_tabular_rho :: eps(double tau,double pres)
  {
    double eps=0.0;
    int r=(int)((1.0/tau-rhomin)/hrho),e=0;
    if (r<0)
      {
        r=0;
        while (e<dimeps-2 && p(tau,e*heps+epsmin)<pres)
          e++;
      }
    else if (r>=dimrho-1)
      {
        r=dimrho-2;
        while (e<dimeps-2 && p(tau,e*heps+epsmin)<pres)
          e++;
      }
    else
      {
        while (e<dimeps-2 && ptable(r,e)<pres)
          e++;
      }
    double lambdar=(1.0/tau-(double)r*hrho-rhomin)/hrho,
      lambdae;
    double psw=ptable(r,e),
      pse=ptable(r+1,e),
      pnw=ptable(r,e+1),
      pne=ptable(r+1,e+1);
    lambdae=(pres-(pse-psw)*lambdar-psw)/((pne+psw-pse-pnw)*lambdar+(pnw-psw));
    eps=lambdae*heps+(double)e*heps+epsmin;
    assert(eps>0.0);
    return eps;
  }

  double Eos_tabular_rho :: dp_deps(double tau,double eps)
  {
    double dp;
    double psw,pse,pnw,pne,lambdar,lambdae;
    int r,e;
    findintable(1.0/tau,eps,r,e,psw,pse,pnw,pne,lambdar,lambdae);
    dp=(pne+psw-pse-pnw)*lambdar+
      (pnw-psw);
    assert(dp>0.0);
    return dp/heps;
  }

  double Eos_tabular_rho :: dp_dtau(double tau,double eps)
  {
    double dp;
    double psw,pse,pnw,pne,lambdar,lambdae;
    int r,e;
    findintable(1.0/tau,eps,r,e,psw,pse,pnw,pne,lambdar,lambdae);
    dp=(pne+psw-pse-pnw)*lambdae+
      (pse-psw);
    return dp/hrho*(-1.0*tau*tau);
  }

  double Eos_tabular_rho :: cs2(double tau,double eps)
  {
    //double lc2=tau*tau*(-dp_dtau(tau,eps)+p(tau,eps)*dp_deps(tau,eps));
    double dpe,dpt,p;
    double psw,pse,pnw,pne,lambdar,lambdae;
    int r,e;
    findintable(1.0/tau,eps,r,e,psw,pse,pnw,pne,lambdar,lambdae);
    p=(pne+psw-pse-pnw)*lambdar*lambdae+
      (pse-psw)*lambdar+ 
      (pnw-psw)*lambdae+
      psw;
    dpe=(pne+psw-pse-pnw)*lambdar+
      (pnw-psw);
    assert(dpe>0.0);
    dpe/=heps;
    dpt=(pne+psw-pse-pnw)*lambdae+
      (pse-psw);
    dpt/=hrho*(-1.0*tau*tau);
    double lc2=tau*tau*(-dpt+p*dpe);
    assert(lc2>0);
    return lc2;
  }

  /******************************************************************/
  /******************************************************************/
  void MhdSolver::init_eos(const Eosmode& peosmode)
  {
    Eosmode::meos_t eos       = peosmode.eos;
    Eosmode::meos_t func_eos  = peosmode.func_eos;
    Eosmode::meos_t adtab_eos = peosmode.adtab_eos;

    if (EOS_PRIM != NULL)
      {
        cerr << "MhdSolver::init_eos(): WARNING" << endl
             << "                       (reinitialization of EOS_PRIM)"
             << endl;
        delete EOS_PRIM;
        EOS_PRIM = NULL;
      }

    if (EOS_CONS != NULL)
      {
        cerr << "MhdSolver::init_eos(): WARNING" << endl
             << "                       (reinitialization of EOS_CONS)"
             << endl;
        delete EOS_CONS;
        EOS_PRIM = NULL;
      }

    assert(EOS_CONS == NULL);
    assert(EOS_PRIM == NULL);

    GAMMA_ID = -1.0;
    R_ID     = -1.0;
    switch (eos)
      {
      case Eosmode::me_default:
        cerr << "MhdSolver::init_eos(): illegal value \"me_default\"!" << endl;
        abort();
        break;
      case Eosmode::me_ideal:
        EOS_CONS  = new eos_cons_ideal_t(peosmode.gamma_id,peosmode.R_id);
        EOS_PRIM  = new eos_prim_ideal_t(peosmode.gamma_id,peosmode.R_id);
        GAMMA_ID  = peosmode.gamma_id;
        R_ID      = peosmode.R_id;
        break;
      case Eosmode::me_waals:
        EOS_CONS  = new eos_cons_waals_t;
        EOS_PRIM  = new eos_prim_waals_t;
        break;
      case Eosmode::me_osborne:
        EOS_CONS  = new eos_cons_osborne_t;
        EOS_PRIM  = new eos_prim_osborne_t;
        break;
      case Eosmode::me_tmv:
        EOS_CONS  = new eos_cons_tmv_t;
        EOS_PRIM  = new eos_prim_tmv_t;
        break;
      case Eosmode::me_nconv:
        EOS_CONS  = new eos_cons_nconv_t;
        EOS_PRIM  = new eos_prim_nconv_t;
        break;
      case Eosmode::me_file:
        EOS_CONS  = new eos_cons_file_t(peosmode.constabfile);
        EOS_PRIM  = new eos_prim_file_t(peosmode.primtabfile);
        break;
      case Eosmode::me_func:
        {
          eos_cons_t *func_eos_cons;
          eos_prim_t *func_eos_prim;
          Init_eos *init_cons,*init_prim;

          switch (func_eos)
            {
            case Eosmode::me_ideal:
              func_eos_cons = new eos_cons_ideal_t(peosmode.func_idealgamma,
                                                   peosmode.func_idealR);
              func_eos_prim = new eos_prim_ideal_t(peosmode.func_idealgamma,
                                                   peosmode.func_idealR);
              break;
            case Eosmode::me_waals:
              func_eos_cons = new eos_cons_waals_t();
              func_eos_prim = new eos_prim_waals_t();
              break;
            case Eosmode::me_osborne:
              func_eos_cons = new eos_cons_osborne_t();
              func_eos_prim = new eos_prim_osborne_t();
              break;
            case Eosmode::me_tmv:
              func_eos_cons = new eos_cons_tmv_t();
              func_eos_prim = new eos_prim_tmv_t();
              break;
            case Eosmode::me_nconv:
              func_eos_cons = new eos_cons_nconv_t();
              func_eos_prim = new eos_prim_nconv_t();
              break;
            case Eosmode::me_sun:
              func_eos_cons = (eos_cons_t*)0;
              func_eos_prim = (eos_prim_t*)0;
              break;
            default:
              {cout << "mhd.cc Z.405\n"<<flush; abort();} //uwe
            }

          if (func_eos == Eosmode::me_sun)
            {
              init_cons = new Sun_init_conseos;
              init_prim = new Sun_init_primeos;
            }
          else
            {
              assert(func_eos_cons);
              assert(func_eos_prim);
              init_cons = new wrap_eos_init_conseos_t(func_eos_cons,func_eos_prim);
              init_prim = new wrap_eos_init_primeos_t(func_eos_prim,func_eos_cons);
            }

          EOS_CONS = new eos_cons_func_t(init_cons);
          EOS_PRIM = new eos_prim_func_t(init_prim);
        }
        break;
      case Eosmode::me_sun:
        adtab_eos = Eosmode::me_sun;
        eos       = Eosmode::me_adtab;
      case Eosmode::me_adtab:
        {
          eos_cons_t *adtab_eos_cons;
          eos_prim_t *adtab_eos_prim;
          Init_eos *init_cons,*init_prim;
          double minfirst,maxfirst;

          switch (adtab_eos)
            {
            case Eosmode::me_ideal:
              adtab_eos_cons = new eos_cons_ideal_t(peosmode.adtab_idealgamma,
                                                    peosmode.adtab_idealR);
              adtab_eos_prim = new eos_prim_ideal_t(peosmode.adtab_idealgamma,
                                                    peosmode.adtab_idealR);
              break;
            case Eosmode::me_waals:
              adtab_eos_cons = new eos_cons_waals_t();
              adtab_eos_prim = new eos_prim_waals_t();
              break;
            case Eosmode::me_osborne:
              adtab_eos_cons = new eos_cons_osborne_t();
              adtab_eos_prim = new eos_prim_osborne_t();
              break;
            case Eosmode::me_tmv:
              adtab_eos_cons = new eos_cons_tmv_t();
              adtab_eos_prim = new eos_prim_tmv_t();
              break;
            case Eosmode::me_nconv:
              adtab_eos_cons = new eos_cons_nconv_t();
              adtab_eos_prim = new eos_prim_nconv_t();
              break;
            case Eosmode::me_sun:
              adtab_eos_cons = (eos_cons_t*)0;
              adtab_eos_prim = (eos_prim_t*)0;
              break;
            default:
              {cout << "mhd.cc Z.464\n"<<flush; abort();} //uwe
            }
          if (adtab_eos == Eosmode::me_sun)
            {
              init_cons = new Sun_init_conseos;
              init_prim = new Sun_init_primeos;
            }
          else
            {
              assert(adtab_eos_cons);
              assert(adtab_eos_prim);
              init_cons = new wrap_eos_init_conseos_t(adtab_eos_cons,adtab_eos_prim);
              init_prim = new wrap_eos_init_primeos_t(adtab_eos_prim,adtab_eos_cons);
            }

          minfirst = (peosmode.adtab_constau
                      ? peosmode.adtab_mintau
                      : 1.0/peosmode.adtab_maxtau);
          maxfirst = (peosmode.adtab_constau
                      ? peosmode.adtab_maxtau
                      : 1.0/peosmode.adtab_mintau);

          EOS_CONS
            = new eos_cons_adtab_t(peosmode.adtab_constau,
                                   peosmode.adtab_minlevel,peosmode.adtab_maxlevel,
                                   peosmode.adtab_dim, peosmode.adtab_dim,
                                   peosmode.adtab_addim, peosmode.adtab_addim,
                                   minfirst, maxfirst,
                                   peosmode.adtab_mineps, peosmode.adtab_maxeps,
                                   peosmode.adtab_tol, *init_cons);

          minfirst = (peosmode.adtab_primtau
                      ? peosmode.adtab_mintau
                      : 1.0/peosmode.adtab_maxtau);
          maxfirst = (peosmode.adtab_primtau
                      ? peosmode.adtab_maxtau
                      : 1.0/peosmode.adtab_mintau);

          EOS_PRIM
            = new eos_prim_adtab_t(peosmode.adtab_primtau,
                                   peosmode.adtab_minlevel, peosmode.adtab_maxlevel,
                                   peosmode.adtab_dim, peosmode.adtab_dim,
                                   peosmode.adtab_addim, peosmode.adtab_addim,
                                   minfirst, maxfirst,
                                   peosmode.adtab_minp, peosmode.adtab_maxp,
                                   peosmode.adtab_tol, *init_prim);
        }
        break;
      default:
        cerr << "MhdSolver::init_eos(): unknown EOS value \""
             << eos << "\"!" << endl;
        abort();
        break;
      }
    if (peosmode.tabmode)
      {
        Eos_tabular *eos_tab;

        switch (peosmode.tabmode)
          {
          case Eosmode::mt_tau:
            eos_tab = new Eos_tabular_tau(EOS_CONS,EOS_PRIM,
                                          peosmode.firstmin,
                                          peosmode.firstmax,
                                          peosmode.epsmin,
                                          peosmode.epsmax,
                                          peosmode.firstdim,
                                          peosmode.epsdim);
            break;
          case Eosmode::mt_rho:
            eos_tab = new Eos_tabular_rho(EOS_CONS,EOS_PRIM,
                                          peosmode.firstmin,
                                          peosmode.firstmax,
                                          peosmode.epsmin,
                                          peosmode.epsmax,
                                          peosmode.firstdim,
                                          peosmode.epsdim);
            break;
          default:
            cerr << "MhdSolver::init_eos(): unknown tabularization mode \""
                 << peosmode.tabmode << "\"!" << endl;
            abort();
          }

        delete EOS_PRIM;
        delete EOS_CONS;
        EOS_CONS = new Eostab_cons(*eos_tab);
        EOS_PRIM = new Eostab_prim(*eos_tab);
      }
  }
  /******************************************************************/
  /******************************************************************/
  /******************************************************************/

  /******************************************************************************
   ******************************************************************************
   ***                                                                        ***
   ***                         CLASS IMPLEMENTATIONS                          ***
   ***                                                                        ***
   ******************************************************************************
   ******************************************************************************/

  /*******************************************************************************
   class Eos_table
  *******************************************************************************/

  void Eos_table::tableinfo(ostream &pout) const
  {
    int level = ((entry_type == et_full) ? get_level() : -1);

    if (level >= 0)
      {
        pout << min1 << " " << min2 << " " << level << endl;
        pout << max1 << " " << min2 << " " << level << endl;
        pout << max1 << " " << max2 << " " << level << endl;
        pout << min1 << " " << max2 << " " << level << endl;
        pout << min1 << " " << min2 << " " << level << endl;
        pout << endl << endl;
      }
  }

  void Eos_table::tableinfo(int plevelsel, ostream &pout) const
  {
    int i,j;
    int show = ((plevelsel < 0) || (plevelsel == get_level()) ? 1 : 0);

    for (i=0;i<dim1;i++)
      for (j=0;j<dim2;j++)
        if (eos_entry[i*dim2+j]->get_table())
          eos_entry[i*dim2+j]->get_table()->tableinfo(plevelsel,pout);

    if (show)
      tableinfo(pout);
  }

  Eos_table::Eos_table(int pfirst_is_tau, int plevel,
                       int pminlevel, int pmaxlevel,
                       int pdim1, int pdim2,
                       int paddim1, int paddim2,
                       double pmin1, double pmax1,
                       double pmin2, double pmax2,
                       double ptol, Init_eos &ptable_init)
    : additional_values(ptable_init.additional_values()),
      first_is_tau(pfirst_is_tau),
      level(plevel), minlevel(pminlevel), maxlevel(pmaxlevel),
      dim1(pdim1), dim2(pdim2), addim1(paddim1), addim2(paddim2),
      min1(pmin1), max1(pmax1), min2(pmin2), max2(pmax2),
      tol(ptol), table_init(ptable_init)
  {
    int i,j;

#ifdef _MHD2D
    mutex = new pthread_mutex_t;
    pthread_mutex_init(mutex,(pthread_mutexattr_t *)0);
    assert(mutex);
#endif

    assert(plevel >= 0);
    assert(maxlevel >= minlevel);
    assert(minlevel >= 0);
    assert(dim1 > 1);
    assert(dim2 > 1);
    assert(addim1 > 1);
    assert(addim2 > 2);
    assert(min1 < max1);
    assert(min2 < max2);
    assert(tol >= 0.0);

    step1 = (max1 - min1) / ((double)dim1 - 1.0);
    step2 = (max2 - min2) / ((double)dim2 - 1.0);

    assert(step1 > 0.0);
    assert(step2 > 0.0);

    eos_entry = new Eos_thin_entry* [dim1*dim2];

    if (refine())
      {
        entry_type = et_thin;

        for (i=0;i<dim1;i++)
          for (j=0;j<dim2;j++)
            eos_entry[i*dim2+j] = new Eos_thin_entry();
      }
    else
      {
        entry_type = et_full;

        switch (additional_values)
          {
          case 0:
            for (i=0;i<dim1;i++)
              {
                for (j=0;j<dim2;j++)
                  {
                    double val[3],dummy0,dummy1;
                    double first  = (first_is_tau ?        min1 + (double)i * step1
                                     : 1.0 / (min1 + (double)i * step1));
                    double second = min2 + (double)j * step2;

                    table_init(first,second,val[0],val[1],val[2],dummy0,dummy1);

                    eos_entry[i*dim2+j] = new Eos_full_entry<3>(val);
                  }
              }
            break;
          case 1:
            for (i=0;i<dim1;i++)
              {
                for (j=0;j<dim2;j++)
                  {
                    double val[4],dummy;
                    double first  = (first_is_tau ?        min1 + (double)i * step1
                                     : 1.0 / (min1 + (double)i * step1));
                    double second = min2 + (double)j * step2;

                    table_init(first,second,val[0],val[1],val[2],val[3],dummy);

                    eos_entry[i*dim2+j] = new Eos_full_entry<4>(val);
                  }
              }
            break;
          case 2:
            for (i=0;i<dim1;i++)
              {
                for (j=0;j<dim2;j++)
                  {
                    double val[5];
                    double first  = (first_is_tau ?        min1 + (double)i * step1
                                     : 1.0 / (min1 + (double)i * step1));
                    double second = min2 + (double)j * step2;

                    table_init(first,second,val[0],val[1],val[2],val[3],val[4]);

                    eos_entry[i*dim2+j] = new Eos_full_entry<5>(val);
                  }
              }
            break;
          default:
            cerr << "ERROR (Eos_table::Eos_table): illegal value of"
                 << " \"additional_values\"!" << endl;
            abort();
          }
      }
  }

  Eos_table::~Eos_table()
  {
    for (int i=dim1-1;i>=0;i--)
      for (int j=dim2-1;j>=0;j--)
        delete eos_entry[i*dim2+j];
    delete [] eos_entry;

#ifdef _MHD2D
    pthread_mutex_destroy(mutex);
    delete mutex;
#endif
  }

  int Eos_table::get_table_value(int pcomp,
                                 double pval1,
                                 double pval2,
                                 double &pret) const
  {
    double lam1,lam2,val_sw,val_se,val_nw,val_ne;
    int idx1 = (int)((pval1 - min1) / step1);
    int idx2 = (int)((pval2 - min2) / step2);
    int range_error = 0;

    // possibly hit regions (table covers region 0):
    //
    //        6  |:   1   :|  5
    //           |:.......:|
    //        ---|---------|---
    //        '''|:''''''':|:''
    //        4  |:   0   :|: 2 
    //        ...|:........|:..
    //        ---|---------|---
    //           |:''''''':|
    //        7  |:   3   :|  8
    //

    if (   (idx1 < 0) || (idx1 > dim1-1)
           || (idx2 < 0) || (idx2 > dim2-1))
      {
        // value not in table range

        if (idx1 < 0)                // ranges 7,3,8
          {
            idx1 = 0;
            if (idx2 < 0)                 // range 7
              idx2 = 0;
            else if (idx2 >= dim2-1)      // range 8
              idx2 = dim2-2;
          }
        else if (idx1 >= dim1-1)     // ranges 6,1,5
          {
            idx1 = dim1-2;
            if (idx2 < 0)                 // range 6
              idx2 = 0;
            else if (idx2 >= dim2-1)      // range 5
              idx2 = dim2-2;
          }
        else if (idx2 < 0)              // range 4
          idx2 = 0;
        else if (idx2 >= dim2-1)     // range 2
          idx2 = dim2-2;

        range_error = get_level()+1;
      }

    if (range_error)
      {
        if (pcomp < 3 + additional_values)
          {
            double first = (first_is_tau ? pval1 : 1.0/pval1);
            double val[5];

            table_init(first,pval2,val[0],val[1],val[2],val[3],val[4]);

            pret = val[pcomp];
          }
        else
          {
            cerr << "ERROR (Eos_table::get_table_value): "
                 << "Unrecoverable range error!" << endl;
            abort();
          }
      }
    else if (entry_type == et_full)
      {
        lam1 = (pval1 - (double)idx1 * step1 - min1) / step1;
        lam2 = (pval2 - (double)idx2 * step2 - min2) / step2;

        val_sw = table(pcomp,idx1  ,idx2  );
        val_se = table(pcomp,idx1+1,idx2  );
        val_nw = table(pcomp,idx1  ,idx2+1);
        val_ne = table(pcomp,idx1+1,idx2+1);

        pret =   (val_ne + val_sw - val_se - val_nw) * lam1 * lam2
          + (val_se - val_sw) * lam1
          + (val_nw - val_sw) * lam2
          + val_sw;
      }
    else
      {
        if (!eos_entry[idx1*dim2+idx2]->get_table())
          {
#ifdef _MHD2D
            pthread_mutex_lock(mutex);
            if (!eos_entry[idx1*dim2+idx2]->get_table())
#endif
              eos_entry[idx1*dim2+idx2]->attach_table(
                                                      new Eos_table(first_is_tau,
                                                                    get_level()+1,
                                                                    minlevel,maxlevel,
                                                                    addim1,addim2,addim1,addim2,
                                                                    min1 + (double)(idx1  ) * step1,
                                                                    min1 + (double)(idx1+1) * step1,
                                                                    min2 + (double)(idx2  ) * step2,
                                                                    min2 + (double)(idx2+1) * step2,
                                                                    tol,table_init));
#ifdef _MHD2D
            pthread_mutex_unlock(mutex);
#endif
          }
        range_error
          = eos_entry[idx1*dim2+idx2]->get_table()
          ->get_table_value(pcomp,pval1,pval2,pret);
      }

    return range_error;  
  }

  int Eos_table::get_table_derivative(derivative_type_t pdtype, int pcomp,
                                      double ptau, double psecond,
                                      double &pret) const
  {
    double lam1,lam2,val_sw,val_se,val_nw,val_ne;
    double val1 = (first_is_tau ? ptau : 1.0/ptau);
    double val2 = psecond;
    int idx1 = (int)((val1 - min1) / step1);
    int idx2 = (int)((val2 - min2) / step2);
    int range_error = 0;

    if (   (idx1 < 0) || (idx1 > dim1-1)
           || (idx2 < 0) || (idx2 > dim2-1))
      {
        range_error = get_level()+1;
      }

    if (range_error)
      {
        cerr << "ERROR (Eos_table::get_table_derivative()): "
             << "cannot compute derivatives outside of table!" << endl;
        abort();
      }
    else if (entry_type == et_full)
      {
        lam1 = (val1 - (double)idx1 * step1 - min1) / step1;
        lam2 = (val2 - (double)idx2 * step2 - min2) / step2;

        val_sw = table(pcomp,idx1  ,idx2  );
        val_se = table(pcomp,idx1+1,idx2  );
        val_nw = table(pcomp,idx1  ,idx2+1);
        val_ne = table(pcomp,idx1+1,idx2+1);

        switch (pdtype)
          {
          case dt_tau:
            pret  =   (  (val_ne + val_sw - val_se - val_nw) * lam2
                         + (val_se - val_sw))
              / step1;
            pret /= (first_is_tau ? 1.0 : (- ptau * ptau));
            break;
          case dt_second:
            pret  =   (  (val_ne + val_sw - val_se - val_nw) * lam1
                         + (val_nw - val_sw))
              / step2;
            break;
          default:
            cerr << "ERROR (Eos_table::get_table_derivative()): "
                 << "unknown type of derivative!" << endl;
            abort();
          }
      }
    else
      {
        if (!eos_entry[idx1*dim2+idx2]->get_table())
          {
#ifdef _MHD2D
            pthread_mutex_lock(mutex);
            if (!eos_entry[idx1*dim2+idx2]->get_table())
#endif
              eos_entry[idx1*dim2+idx2]->attach_table(
                                                      new Eos_table(first_is_tau,
                                                                    get_level()+1,
                                                                    minlevel,maxlevel,
                                                                    addim1,addim2,addim1,addim2,
                                                                    min1 + (double)(idx1  ) * step1,
                                                                    min1 + (double)(idx1+1) * step1,
                                                                    min2 + (double)(idx2  ) * step2,
                                                                    min2 + (double)(idx2+1) * step2,
                                                                    tol,table_init));
#ifdef _MHD2D
            pthread_mutex_unlock(mutex);
#endif
          }
        range_error
          = eos_entry[idx1*dim2+idx2]->get_table()
          ->get_table_derivative(pdtype,pcomp,ptau,psecond,pret);
      }

    return range_error;
  }

  void Eos_table::tableinfo(int plevelsel, const char *pfname) const
  {
    cerr << "Eos_table::tableinfo(): creating file \"" << pfname << "\"" << endl;

    ofstream out(pfname,ios::out|ios::trunc);
    tableinfo(plevelsel,out);
    out.close();
  }


  /******************************************************************/
  /******************************************************************/
  /******************************************************************/
  /******************************************************************/
  /******************************************************************/



  /******************************************************************************
   ******************************************************************************
   ***                                                                        ***
   ***                         CLASS IMPLEMENTATIONS                          ***
   ***                                                                        ***
   ******************************************************************************
   ******************************************************************************/

  /*******************************************************************************
   class Sun_eos
  *******************************************************************************/

  float Sun_eos::zbrent(int flag, float x1, float x2,
                        float tol, float f1, float f2) const
  {
    int iter;
    float a=x1,b=x2,c=x2,min1,min2,fa,fb;
    float d=0.0, e=0.0;
	
    if(flag)
      {
        fa=ener_b(a,f1,f2);
        fb=ener_b(b,f1,f2);
      }
    else
      {
        fa=ener_a(a,f1,f2);
        fb=ener_a(b,f1,f2); 
      }  

    float fc,q,r,s,tol1,xm;
    float p =0.0;

    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
      cerr<<"Root must be bracketed in zbrent (" 
          << flag << ": fa = " << fa << ", fb = " << fb << ")" <<endl;
	
    fc=fb;
    for (iter=1;iter<=ITMAX;iter++)
      {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
          {
            c=a;
            fc=fa;
            e=d=b-a;
          }
        if (fabs(fc) < fabs(fb))
          {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
          }

        tol1=2.0*TOL*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) 
          return b;

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
          {	
            s=fb/fa;
            if (a == c)
              {
                p=2.0*xm*s;
                q=1.0-s;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
              }
            if (p > 0.0)
              q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2))
              {
                e=d;
                d=p/q;
              }
            else
              {
                d=xm;
                e=d;
              }
          }
        else
          {
            d=xm;
            e=d;
          }
        a=b;
        fa=fb;
        if (fabs(d) > tol1)
          b += d;
        else
          b += SIGN(tol1,xm);

        if(flag)
          fb=ener_b(b,f1,f2);
        else
          fb=ener_a(b,f1,f2);  
      }

    return 0.0;
  }

  //*****************************************************************************

  void Sun_eos::eos_a(float ptau, float peps,
                      double *pcs2,double *pp,double *pT,double *pdpdeps) const
  {
    //     *** For given values of the internal energy (eps)and specific volume
    //     *** (tau=1/rho), this subroutine determines the internal energy (eps)
    //     *** the velocity of sound squared (cq), the partial derivatives
    //     *** of p w.r.t. to tau (dptau) and eps (dpeps), and the "gamma"
    //     *** gam1 = cq/(p*tau) for a pure, ionizing hydrogen gas.
    //     *** ALL UNITS ARE CGS!!! ***
    //         Version 1.0  --  M. Schuessler  --  March 2000
    //         Status : 21.3.2000

    double boltem,arg,arg3,phi,nh,ne,x,xmx,arg5,f1,f2,nad,mu,gam1,dptau;

    //     Determine first the temperature (temp) through the pressure
    //     equation and the Saha equation. This must be done by iteration/
    //     since the equation is transcendental (Brent method, Numerical
    //     Recipes ch. 9.3)

    *pT=zbrent(0,1.e-8,0.5e8,1.e-8,peps,ptau);
 
    //     Next is calculation of the internal energy
    
    boltem = boltz*(*pT);
    arg = chi/boltem; 
    phi = cii*(pow(*pT,1.5))*exp(-arg);
    nh = 1./(ptau*amass*hmu);
    ne = -.5*phi + sqrt(.25*(pow(phi,2)) + nh*phi);

    //    Now we can determine the pressure...

    *pp = (ne+nh)*boltem;

    //    Now we can determine the velocity of sound
    //    (eq. 36 of my memo) and gam1
    //    and the molecular weight mu with n_alpha=n_h 
 
    x = ne/nh;
    mu=mu1/(1.+x);

    xmx = (1.-x)*x;
    arg5 = arg + 2.5;
    arg3 = arg + 1.5;

    *pcs2 = (*pp)*ptau*(5. + xmx*(pow(arg5,2)))/(3. + xmx*(arg3*arg5 - arg));
    gam1  = (*pcs2)/((*pp)*ptau);

    //     Now we need the adiabatic gradient, nabla_ad (nad)

    f1 = (1.+x)*(1. + .5*x*(1.-x)*arg5);
    f2 = 2.5*(1.+x) + x*(1.-x)*(pow(arg5,2)) - .5*x*pow(((1.-x)*arg5),2);

    nad = f1/f2;
      
    *pdpdeps = (gam1)*nad/ptau;
      
    dptau    = -(*pcs2)/(pow(ptau,2)) + (*pp)*(*pdpdeps);
  }

  //*****************************************************************************

  void Sun_eos::eos_b(float ptau, float pp,
                      double *pcs2, double *peps, double *pT) const
  {
    //     *** For given values of the pressure (p) und specific volume
    //     *** (tau=1/rho), this subroutine determines the internal energy (eps)
    //     *** the velocity of sound squared (cq), the partial derivatives
    //     *** of p w.r.t. to tau (dptau) and eps (dpeps), and the "gamma"
    //     *** gam1 = cq/(p*tau) for a pure, ionizing hydrogen gas.
    //     *** ALL UNITS ARE CGS!!! ***
    //         Version 1.0  --  M. Schuessler  --  March 2000
    //         Status : 21.3.2000

    double boltem,arg,arg3,phi,nh,ne,x,xmx,arg5,f1,f2,nad,mu,gam1,dpeps,dptau;

    //  Determine first the temperature (temp) through the pressure
    //  equation and the Saha equation. This must be done by iteration/
    //  since the equation is transcendental (Brent method, Numerical
    //  Recipes ch. 9.3)

    *pT=zbrent(1,1.e-8,0.5e8,1.e-8,pp,ptau);
 
    //  Next is calculation of the internal energy

    boltem = boltz*(*pT);
    arg = chi/boltem; 
    phi = cii*(pow(*pT,1.5))*exp(-arg);
    nh = 1./(ptau*amass*hmu);
    ne = -.5*phi + sqrt(.25*(pow(phi,2)) + nh*phi);

    *peps = ptau*(1.5*(ne+nh)*boltem + ne*chi);

    //  Now we can determine the velocity of sound
    //  (eq. 36 of my memo) and gam1
    //  and the molecular weight mu with n_alpha=n_h 
 
    x = ne/nh;
    mu=mu1/(1.+x);

    xmx = (1.-x)*x;
    arg5 = arg + 2.5;
    arg3 = arg + 1.5;

    *pcs2 = pp*ptau*(5. + xmx*(pow(arg5,2)))/(3. + xmx*(arg3*arg5 - arg));
    gam1  = (*pcs2)/(pp*ptau);

    //  Now we need the adiabatic gradient, nabla_ad (nad)

    f1 = (1.+x)*(1. + .5*x*(1.-x)*arg5);
    f2 = 2.5*(1.+x) + x*(1.-x)*(pow(arg5,2)) - .5*x*pow(((1.-x)*arg5),2);

    nad = f1/f2;
      
    dpeps = gam1*nad/ptau;
      
    dptau = -(*pcs2)/(pow(ptau,2)) + pp*(dpeps); 
  }

  //*****************************************************************************

  float Sun_eos::ener_a(float temp, float feps, float ftau) const
  {
    //  *** This function is used to determine the temperature (temp)
    //  *** for given internal energy (feps) and specific volume (ftau)
    //  *** by root finding through Numerical Recipes routine zbrent

    //  *** ALL UNITS ARE CGS!!! ***

    //      Version 1.0  --  M. Schuessler  --  March 2000
    //      Status : 20.3.2000

    double boltem = boltz*temp;
    double arg = chi/boltem; 
    double phi = cii*(pow(temp,(float)1.5))*exp(-arg);
    double nh = 1./(hmu*amass*ftau);
    double ne = -.5*phi + sqrt(.25*(pow(phi,2)) + nh*phi);

    float ener = ftau*(1.5*(ne+nh)*boltem + ne*chi) - feps;

    return ener;
  }

  //*****************************************************************************

  float Sun_eos::ener_b(float temp, float fp, float ftau) const
  {
    //  *** This function is used to determine the temperature (temp)
    //  *** for given pressure (fp) and specific volume (ftau)
    //  *** by root finding through Numerical Recipes routine zbrent

    //  *** ALL UNITS ARE CGS!!! ***

    //      Version 1.0  --  M. Schuessler  --  March 2000
    //      Status : 20.3.2000

    double boltem = boltz*temp;
    double arg = chi/boltem; 
    double phi = cii*(pow(temp,(float)1.5))*exp(-arg);
    double nh = 1./(hmu*amass*ftau);
    double ne = -.5*phi + sqrt(.25*(pow(phi,2)) + nh*phi);

    float ener=(ne+nh)*boltem -fp;

    return ener;
  }

  /*******************************************************************************
   class Sun_init_conseos
  *******************************************************************************/

  void Sun_init_conseos::operator()(double ptau, double peps,
                                    double &pcs2, double &pp,
                                    double &pT, double &pdpdeps,
                                    double &pdummy) const
  {
    const double Rgas = 8.314e7;
    double rho0_sun=2.53700e-6;
    double T0_sun=16404.;
    double mu0_sun=1.18596;
    double p0_sun=rho0_sun*T0_sun*Rgas/mu0_sun;

    ptau *= rho0_sun;
    peps /= eps0;

    double ltau = ptau / rho0;
    double leps = peps * eps0;
  
    eos_a(ltau,leps,&pcs2,&pp,&pT,&pdpdeps);
  
    pcs2    /= cs20;
    pp      /= p0;
    pT      /= T0;
    pdpdeps *= eps0/p0;
    pdummy   = -1.0;
  
    // DEBUGGING

    double lp = pp * p0;
    double ltcs2,lteps,ltT;

    eos_b(ltau,lp,&ltcs2,&lteps,&ltT);

    ltcs2 /= cs20;
    lteps /= eps0;
    ltT   /= T0;

    if (   (fabs(pcs2 - ltcs2)/(0.5*(pcs2+ltcs2)) > 5.0e-7)
           || (fabs(peps - lteps)/(0.5*(peps+lteps)) > 5.0e-7)
           || (fabs(pT - ltT)    /(0.5*(pT  +ltT  )) > 5.0e-7))
      {
        cerr << "Sun_init_conseos::operator(): " << endl;
        cerr << "   (tau,eps) ==> (cs2,p  ,T): (" << ptau << "," << peps << ") ==> ("
             << pcs2 << "," << pp << "," << pT << ")" << endl;
        cerr << "   (tau,p  ) ==> (cs2,eps,T): (" << ptau << "," << pp   << ") ==> ("
             << ltcs2 << "," << lteps << "," << ltT << ")" << endl;
      }
  
    pp *= p0_sun;
    pT *= T0_sun;

  }

  /*******************************************************************************
   class Sun_init_primeos
  *******************************************************************************/

  void Sun_init_primeos::operator()(double ptau, double pp,
                                    double &pcs2, double &peps,
                                    double &pT, double &pdummy0,
                                    double &pdummy1) const
  {
    const double Rgas = 8.314e7;
    double rho0_sun=2.53700e-6;
    double T0_sun=16404.;
    double mu0_sun=1.18596;
    double p0_sun=rho0_sun*T0_sun*Rgas/mu0_sun;

    ptau *= rho0_sun;
    pp /= p0_sun;

    double ltau = ptau / rho0;
    double lp   = pp * p0;
    eos_b(ltau,lp,&pcs2,&peps,&pT);
  
    pcs2    /= cs20;
    peps    /= eps0;
    pT      /= T0;
    pdummy0  = -1.0;
    pdummy1  = -1.0;
    // DEBUGGING
    double leps = peps * eps0;
    double ltcs2,ltp,ltT,ldummy;

    eos_a(ltau,leps,&ltcs2,&ltp,&ltT,&ldummy);

    ltcs2 /= cs20;
    ltp   /= p0;
    ltT   /= T0;

    if (   (fabs(pcs2 - ltcs2)/(0.5*(pcs2+ltcs2)) > 5.0e-7)
           || (fabs(pp   - ltp  )/(0.5*(pp  +ltp  )) > 5.0e-7)
           || (fabs(pT   - ltT)  /(0.5*(pT  +ltT  )) > 5.0e-7))
      {
        cerr << "Sun_init_primeos::operator(): " << endl;
        cerr << "   (tau,p  ) ==> (cs2,eps,T): (" << ptau << "," << pp   << ") ==> ("
             << pcs2 << "," << peps << "," << pT << ")" << endl;
        cerr << "   (tau,eps) ==> (cs2,p  ,T): (" << ptau << "," << peps << ") ==> ("
             << ltcs2 << "," << ltp << "," << ltT << ")" << endl;
      }
  
    peps *= eps0;
    pT *= T0_sun;
  }

} // end namespace Mhd 


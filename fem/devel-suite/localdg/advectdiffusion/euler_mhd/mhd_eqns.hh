#ifndef MHD_H_INCLUDED
#define MHD_H_INCLUDED

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "mhd_fluxes.hh"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

using namespace std;

namespace Mhd {
  
typedef Flux<Mhd> flux_t;
typedef Flux_dw<Mhd>     flux_dw_t;
typedef Flux_lf<Mhd>     flux_lf_t;
typedef Flux_dwc<Mhd>    flux_dwc_t;
typedef Flux_hlle<Mhd>   flux_hlle_t;
typedef Flux_hllem<Mhd>  flux_hllem_t;
typedef Flux_vfroe<Mhd>  flux_vfroe_t;
typedef Flux_hllemc<Mhd> flux_hllemc_t;


class MhdSolver {
public:
  typedef double Vec9[9];
  typedef double VEC[8];
  typedef double VEC1D[8];
  typedef double MAT[8][8];
private:
  static void (*ipflux_1d)(const VEC1D, const VEC1D, VEC1D *, double);
 
  double max(const double p1, const double p2) {return (p1 >= p2) ? p1 : p2;}
  double min(const double p1, const double p2) {return (p1 <= p2) ? p1 : p2;}

  double p_1d(const VEC1D pu,double gamma) {
    double lp,lrho,lrhoux,lrhouy,lrhouz,lBx,lBy,lBz,lrhoE;
    lrho   = pu[0];
    lrhoux = pu[1];
    lrhouy = pu[2];
    lrhouz = pu[3];
    lBx    = pu[4];
    lBy    = pu[5];
    lBz    = pu[6];
    lrhoE  = pu[7];
    assert(lrho > 0.0);
    lp  = (gamma - 1.0) 
      * (lrhoE - 0.5 * (  lrhoux * lrhoux
			  + lrhouy * lrhouy
			  + lrhouz * lrhouz ) / lrho
	 - ( lBy * lBy + lBz * lBz ) / ( 8.0 * M_PI ));
    assert(lp > 0.0);
    return lp;
  }

  double dt_local(const VEC1D prdate,const VEC1D prdatn,double gamma) {
    double lval;
    
    assert(gamma>1.0); 
    
    double lrhoe = prdate[0], lrhon = prdatn[0];
    assert(lrhoe > 0.0 && lrhon > 0.0);         
    
    double luxe = prdate[1] / lrhoe, luxn = prdatn[1] / lrhon;
    double lBxe = prdate[4], lBxn = prdatn[4];
    double lBye = prdate[5], lByn = prdatn[5];
    double lBze = prdate[6], lBzn = prdatn[6];
    double lpe  = p_1d(prdate,gamma), lpn  = p_1d(prdatn,gamma);
    
    lval =   max(fabs(luxe),fabs(luxn))
      + sqrt((  (  max(lBxe * lBxe, lBxn * lBxn)
		   + max(lBye * lBye, lByn * lByn)
		   + max(lBze * lBze, lBzn * lBzn) ) / (4.0 * M_PI)
		+ gamma * max(lpe,lpn) ) / min(lrhoe,lrhon));
    
    return lval;
  }
  double relax (const Vec9 &pdatae, const Vec9 &pdatan,const double (&x)[3], 
		Vec9 &pret) {
    double lgamma1,lGamma1;
    double lrepsl,lrepsr,lreps1l,lreps1r,lreps2l,lreps2r;
    cons_t ldatae(EOS_CONS,EOS_PRIM,pdatae);
    cons_t ldatan(EOS_CONS,EOS_PRIM,pdatan);
    
    VEC1D  ldat1de,ldat1dn,lret;
    double lBx;
    
    /* Relaxation (Korrektur fuer allg. Zustandsgleichungen) */
    lgamma1 = max(ldatae.gamma(),ldatan.gamma());
    lGamma1 = max(1.+EOS_CONS->dpdeps(ldatae)*ldatae.tau(),
		  1.+EOS_CONS->dpdeps(ldatan)*ldatan.tau());
    lgamma1 = max(lgamma1,lGamma1);
    lrepsl  =   ldatae[7]
      - 0.5/ldatae[0] * (  ldatae[1]*ldatae[1]
			   + ldatae[2]*ldatae[2]
			   + ldatae[3]*ldatae[3])
      - 1./(8.*M_PI) * (  ldatae[4]*ldatae[4]
			  + ldatae[5]*ldatae[5]
			  + ldatae[6]*ldatae[6]);
    lrepsr  =   ldatan[7]
      - 0.5/ldatan[0] * (  ldatan[1]*ldatan[1]
			   + ldatan[2]*ldatan[2]
			   + ldatan[3]*ldatan[3])
      - 1./(8.*M_PI) * (  ldatan[4]*ldatan[4]
			  + ldatan[5]*ldatan[5]
			  + ldatan[6]*ldatan[6]);
    lreps1l = ldatae.p()/(lgamma1-1.0);
    lreps1r = ldatan.p()/(lgamma1-1.0);
    lreps2l = lrepsl - lreps1l;
    lreps2r = lrepsr - lreps1r;
    ldatae[7] -= lreps2l;
    ldatan[7] -= lreps2r;
    /* END: Relaxation */
    /* Magnetfeld Korrektur */
    lBx = 0.5 * (ldatae[4] + ldatan[4]);
    ldat1de[0] = ldatae[0];
    ldat1de[1] = ldatae[1];
    ldat1de[2] = ldatae[2];
    ldat1de[3] = ldatae[3];
    ldat1de[4] = ldatae[4];
    ldat1de[5] = ldatae[5];
    ldat1de[6] = ldatae[6];
    ldat1de[7] = ldatae[7] - (ldatae[4] * ldatae[4]) / (8.0 * M_PI);;
    ldat1dn[0] = ldatan[0];
    ldat1dn[1] = ldatan[1];
    ldat1dn[2] = ldatan[2];
    ldat1dn[3] = ldatan[3];
    ldat1dn[4] = ldatan[4];
    ldat1dn[5] = ldatan[5];
    ldat1dn[6] = ldatan[6];
    ldat1dn[7] = ldatan[7] - (ldatan[4] * ldatan[4]) / (8.0 * M_PI);;
    /* END: Magnetfeld Korrektur */
    (*ipflux_1d)(ldat1de,ldat1dn,&lret,lgamma1);
    /* Magnetfeld Korrektur */
    lret[1]    -= (lBx * lBx) / (8.0 * M_PI);
    /* END: Magnetfeld Korrektur */
    /* Relaxation (Korrektur fuer allg. Zustandsgleichungen) */
    if (lret[0]>0.0)
      lret[7]+=lret[0]/ldat1de[0]*lreps2l;
    else
      lret[7]+=lret[0]/ldat1dn[0]*lreps2r;
    /* END: Relaxation */
    pret[0] = lret[0];
    pret[1] = lret[1];
    pret[2] = lret[2];
    pret[3] = lret[3];
    pret[4] = lret[4];
    pret[5] = lret[5];
    pret[6] = lret[6];
    pret[7] = lret[7];
    pret[8] = 0.0;
    return dt_local(ldat1de,ldat1dn,lgamma1);
  }
  double ipflux(const Vec9 &pdatae, const Vec9 &pdatan, const double (&x)[3], 
			   Vec9 &pret) {
    VEC1D  ldat1de,ldat1dn,lret;
    double lBxe,lBxn,lBx;
    ldat1de[0] = pdatae[0];
    ldat1de[1] = pdatae[1];
    ldat1de[2] = pdatae[2];
    ldat1de[3] = pdatae[3];
    ldat1de[4] = pdatae[4];
    ldat1de[5] = pdatae[5];
    ldat1de[6] = pdatae[6];
    ldat1de[7] = pdatae[7];
    ldat1dn[0] = pdatan[0];
    ldat1dn[1] = pdatan[1];
    ldat1dn[2] = pdatan[2];
    ldat1dn[3] = pdatan[3];
    ldat1dn[4] = pdatan[4];
    ldat1dn[5] = pdatan[5];
    ldat1dn[6] = pdatan[6];
    ldat1dn[7] = pdatan[7];
    lBxe = ldat1de[4];
    lBxn = ldat1dn[4];
    lBx  = 0.5 * (lBxe + lBxn);
    ldat1de[7] -= (lBxe * lBxe) / (8.0 * M_PI);
    ldat1dn[7] -= (lBxn * lBxn) / (8.0 * M_PI);
    ldat1de[4]=lBx;
    ldat1dn[4]=lBx;
    (*ipflux_1d)(ldat1de,ldat1dn,&lret,gamma_id());
    lret[1]    -= (lBx * lBx) / (8.0 * M_PI);
    pret[0] = lret[0];
    pret[1] = lret[1];
    pret[2] = lret[2];
    pret[3] = lret[3];
    pret[4] = lret[4];
    pret[5] = lret[5];
    pret[6] = lret[6];
    pret[7] = lret[7];
    pret[8] = 0.0;
    return dt_local(ldat1de,ldat1dn,gamma_id());
  }
public:
    class Eosmode
    {
      public:
        typedef enum {me_default=0,me_ideal,me_waals,me_osborne,
                      me_tmv,me_nconv,me_file,me_func,me_adtab,me_sun} meos_t;
        typedef enum {mt_none=0,mt_tau,mt_rho} mtab_t;

        meos_t eos,adtab_eos,func_eos;
        mtab_t tabmode;
        char *constabfile,*primtabfile;
        int firstdim,epsdim;
        int adtab_constau,adtab_primtau,adtab_dim,adtab_addim;
        int adtab_minlevel,adtab_maxlevel;
        double firstmin,firstmax,epsmin,epsmax,gamma_id,R_id;
        double adtab_mintau,adtab_maxtau,adtab_minp,adtab_maxp;
        double adtab_mineps,adtab_maxeps,adtab_tol;
        double adtab_idealgamma,adtab_idealR,func_idealgamma,func_idealR;

        Eosmode()
          : eos(me_default), adtab_eos(me_default), func_eos(me_default),
            tabmode(mt_none), constabfile(0), primtabfile(0),
            firstdim(0), epsdim(0), adtab_constau(0), adtab_primtau(0),
            adtab_dim(2), adtab_addim(2), adtab_minlevel(0),
            firstmin(1.0), firstmax(2.0),
            epsmin(1.0), epsmax(2.0), gamma_id(-1.0), R_id(1.0),
            adtab_mintau(0.1), adtab_maxtau(1.0), adtab_minp(1.0),
            adtab_maxp(10.0), adtab_mineps(2.0), adtab_maxeps(20.0),
            adtab_idealgamma(-1.0), adtab_idealR(1.0),
            func_idealgamma(-1.0), func_idealR(1.0)
          {}         
    };
  MhdSolver(MhdSolver::Eosmode::meos_t peos,double pgamma=1.4,double pR=1.0) {
    MhdSolver::FLUXCHOICE = MhdSolver::mf_dw;
    MhdSolver::Eosmode leosmode;
    if (peos==MhdSolver::Eosmode::me_ideal) { // Perfektes Gas 
      assert(pgamma > 1.0);
      leosmode.eos = MhdSolver::Eosmode::me_ideal;
      leosmode.gamma_id = pgamma;
      leosmode.R_id     = pR;
    } else 
      leosmode.eos = peos;
    MhdSolver::init_eos(leosmode);
    MhdSolver::init(MhdSolver::FLUXCHOICE);
  }
  double operator()(const Vec9 &pdatae, const Vec9 &pdatan, const double (&x)[3], 
		    Vec9 &pret) {
    if (MhdSolver::FLUXCHOICE < MhdSolver::first_rgflux) {
      if (eos_is_ideal())
	return ipflux(pdatae,pdatan,x,pret);
      else
	return relax(pdatae,pdatan,x,pret);
    }
    else abort();
  }
  double operator()(const double (&pdatae)[5], const double (&pdatan)[5], const double (&x)[3], double (&pret)[5]) {
    Vec9 ldatae=
      {pdatae[0],pdatae[1],pdatae[2],pdatae[3],0.,0.,0.,pdatae[4],0.};
    Vec9 ldatan=
      {pdatan[0],pdatan[1],pdatan[2],pdatan[3],0.,0.,0.,pdatan[4],0.};
    Vec9 lret;
    double dt=(*this)(ldatae,ldatan,x,lret);
    pret[0]=lret[0],pret[1]=lret[1],pret[2]=lret[2],pret[3]=lret[3],pret[4]=lret[7];
    return dt;
  }

  double p(const Vec9& pdat) const {
    cons_t cons(EOS_CONS, EOS_PRIM,pdat);
    return cons.p();
  }

  double p(const double (&pdat)[5]) const {
    Vec9 ldat={pdat[0],pdat[1],pdat[2],pdat[3],0.,0.,0.,pdat[4],0.};
    return p(ldat);
  }

  /************************************************************************************/
 public:
  static const int first_rgflux = 10;
  typedef enum {mf_bct=0,mf_dw,mf_hlle,mf_hllemg,mf_hlleml,mf_roe,
		mf_rgdw=first_rgflux,mf_rgdwc,mf_rghlle,mf_rghllem,
		mf_rghllemc,mf_rglf,mf_rgvfroe} mflux_t;
  
  static const double s2;
  static mflux_t FLUXCHOICE;
  static double GAMMA_ID,R_ID;
  static eos_cons_t *EOS_CONS;
  static eos_prim_t *EOS_PRIM;
  
  static double gamma_id() {
    assert(GAMMA_ID > 1.0);
    return GAMMA_ID;
  }
  
  static int eos_is_ideal() {
    return (GAMMA_ID > 1.0);
  }

  static void init_eos(const Eosmode&);
  public :
  typedef mflux_t flux1d_t ;
  static const double MhdEPS;
  static const int ewave;
  static double maxtime;
  static void (*initialize)(const double px, VEC1D *pu);
  static int docheck;
  static int tstep_auto;
  static int use_hll_fix;
  static int verbose_mode;
  static int use_point_fix;
  static int use_entropy_fix;
  static int use_velocity_fix;
  static int check_rmv_errors;
  static int use_reference_fix;
  static int check_2D_conservation;
  static int hll_use_roe_mean_values;
  static double visc;
  static char algorithm[20];
  static int undo;
  static int hllcount;
  static int point_fixes;
  static int entropy_fixes;
  static int velocity_fixes;
  static int reference_fixes;
  static int conservation_errors_2D;
  static double max_check_rmv_error;
  static double max_conservation_error_2D;
  static void (*init_ipflux)(void);
  static flux_t *rgflux_1d;
  static void init(flux1d_t) ;
};


/***************************************************************************
     functions and numerical fluxes
****************************************************************************/

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/
/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class T, const int vsize> class MhdVector
*******************************************************************************/

template <class T, const int vsize>
inline void MhdVector<T,vsize>::init(const T pvec[vsize])
{
  for (int i=0;i<vsize;i++)
    vec[i] = pvec[i];
}

template <class T, const int vsize>
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator=(const vec_t &pvec)
{
  if (this != &pvec)
    init(pvec);
  return *this;
}

template <class T, const int vsize> 
inline T& MhdVector<T,vsize>::operator[](int pidx)
{
  assert((0 <= pidx) && (pidx < vsize));
  return vec[pidx];
}

template <class T, const int vsize> 
inline T MhdVector<T,vsize>::operator[](int pidx) const
{
  assert((0 <= pidx) && (pidx < vsize));
  return vec[pidx];
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator+=(const vec_t &pvec)
{
  for (int i=0;i<vsize;i++)
   (*this)[i] += pvec[i];
  return *this;
}

template <class T, const int vsize>
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator-=(const vec_t &pvec)
{
  for (int i=0;i<vsize;i++)
   (*this)[i] -= pvec[i];
  return *this;
}

template <class T, const int vsize>
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator*=(const mat_t &pmat)
{
  vec_t hvec;

  for (int j=0;j<vsize;j++)
  {
    for (int i=0;i<vsize;i++)
      hvec[j] += (*this)[i] * pmat(i,j);
  }

  init(hvec);

  return *this;
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator*=(const T pt)
{
  for (int i=0;i<vsize;i++)
    (*this)[i] *= pt;
  return *this;
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t& MhdVector<T,vsize>::operator/=(const T pt)
{
  assert(pt != (T)0);
  for (int i=0;i<vsize;i++)
    (*this)[i] /= pt;
  return *this;
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t MhdVector<T,vsize>::operator+(const vec_t &pvec)
const
{
  vec_t ret(*this);
  ret += pvec;
  return ret;
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t MhdVector<T,vsize>::operator-(const vec_t &pvec)
const
{
  vec_t ret(*this);
  ret -= pvec;
  return ret;
}

template <class T, const int vsize>
inline typename MhdVector<T,vsize>::vec_t MhdVector<T,vsize>::operator*(const mat_t &pmat)
const
{
  vec_t ret(*this);
  ret *= pmat;
  return ret;
}

template <class T, const int vsize>
inline double MhdVector<T,vsize>::operator*(const vec_t &pvec)
const
{
  double d = 0.0;

  for (int i=0;i<vsize;i++)
    d += (*this)[i] * pvec[i];

  return d;
}

template <class T, const int vsize>
inline typename MhdVector<T,vsize>::vec_t MhdVector<T,vsize>::operator*(const T pt)
const
{
  vec_t ret(*this);
  ret *= pt;
  return ret;
}

template <class T, const int vsize> 
inline typename MhdVector<T,vsize>::vec_t MhdVector<T,vsize>::operator/(const T pt)
const
{
  vec_t ret(*this);
  assert(pt != (T)0);
  ret /= pt;
  return ret;
}

template <class T, const int vsize> 
inline void MhdVector<T,vsize>::write(ostream &pout) const
{
  pout << '(';
  for (int i=0;i<vsize;i++)
  {
    pout << (*this)[i];
    if (i<(vsize-1))
      pout << ',';
  }
  pout << ')';
}

/*******************************************************************************
   template <class T, const int msize> class Matrix
*******************************************************************************/

template <class T, const int msize>
inline void Matrix<T,msize>::init(const T pmat[msize][msize])
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) = pmat[i][j];
}

template <class T, const int msize>
inline void Matrix<T,msize>::init(const vec_t pmat[msize])
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) = pmat[i][j];
}

template <class T, const int msize>
inline void Matrix<T,msize>::clear()
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) = (T)0;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator=(const mat_t &pmat)
{
  if (this != &pmat)
    init(pmat);

  return *this;
}

template <class T, const int msize>
inline T& Matrix<T,msize>::operator()(int pi,int pj)
{
  assert((0 <= pi) && (pi < msize));
  assert((0 <= pj) && (pj < msize));
  return mat[pi][pj];
}

template <class T, const int msize>
inline T Matrix<T,msize>::operator()(int pi,int pj) const
{
  assert((0 <= pi) && (pi < msize));
  assert((0 <= pj) && (pj < msize));
  return mat[pi][pj];
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator+=(const mat_t &pmat)
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) += pmat(i,j);

  return *this;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator-=(const mat_t &pmat)
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) -= pmat(i,j);

  return *this;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator*=(const mat_t &pmat)
{
  mat_t hmat;

  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      for (int k=0;k<msize;k++)
        hmat(i,j) += (*this)(i,k) * pmat(k,j);

  init(hmat);

  return *this;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::vec_t Matrix<T,msize>::operator*=(const vec_t &pvec)
{
  vec_t ret;

  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      ret[i] += (*this)(i,j) * pvec[j];

  return ret;
}


template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator*=(const T pt)
{
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) *= pt;

  return *this;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t& Matrix<T,msize>::operator/=(const T pt)
{
  assert(pt != (T)0);
  for (int i=0;i<msize;i++)
    for (int j=0;j<msize;j++)
      (*this)(i,j) /= pt;

  return *this;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t Matrix<T,msize>::operator+(const mat_t &pmat)
const
{
  mat_t ret(*this);
  ret += pmat;
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t Matrix<T,msize>::operator-(const mat_t &pmat)
const
{
  mat_t ret(*this);
  ret -= pmat;
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t Matrix<T,msize>::operator*(const mat_t &pmat)
const
{
  mat_t ret(*this);
  ret *= pmat;
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::vec_t Matrix<T,msize>::operator*(const vec_t &pvec)
const
{
  mat_t hmat(*this);
  vec_t ret = (hmat *= pvec);
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t Matrix<T,msize>::operator*(const T pt)
const
{
  mat_t ret(*this);
  ret *= pt;
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::mat_t Matrix<T,msize>::operator/(const T pt)
const
{
  mat_t ret(*this);
  assert(pt != (T)0);
  ret /= pt;
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::vec_t Matrix<T,msize>::row(int pidx)
const
{
  vec_t ret;
  assert((0 <= pidx) && (pidx < msize));
  for (int j=0;j<msize;j++)
    ret[j] = (*this)(pidx,j);
  return ret;
}

template <class T, const int msize>
inline typename Matrix<T,msize>::vec_t Matrix<T,msize>::col(int pidx)
const
{
  vec_t ret;
  assert((0 <= pidx) && (pidx < msize));
  for (int i=0;i<msize;i++)
    ret[i] = (*this)(i,pidx);
  return ret;
}

template <class T, const int msize>
inline void Matrix<T,msize>::write(ostream &pout) const
{
  for (int i=0;i<msize;i++)
    pout << row(i) << endl;
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                            OUTPUT OPERATORS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

template <class T, const int vsize>
ostream& operator<<(ostream &pout, MhdVector<T,vsize> &pvar)
{
  pvar.write(pout);
  return pout;
}

template <class T, const int msize>
ostream& operator<<(ostream &pout, Matrix<T,msize> &pmat)
{
  pmat.write(pout);
  return pout;
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                            INPUT OPERATORS                             ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

template <class T, const int vsize>
istream& operator>>(istream &pin, MhdVector<T,vsize> &pvar)
{
  for (int j=0;j<pvar.size();j++)
    pin >> pvar[j];
  return pin;
}



/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons
*******************************************************************************/

template <const int size>
inline Eos_cons<size>::Eos_cons(int pdpdeps_avail)
: dpdeps_avail(pdpdeps_avail)
{
  assert((dpdeps_avail == 0) || (dpdeps_avail == 1));
}

template <const int size>
double Eos_cons<size>::gamma(const var_t &pvar) const
{
  return (cs2(pvar) / (pvar.tau() * p(pvar)));
}

template <const int size>
double Eos_cons<size>::dpdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();

  assert(ltau > 0.0);

  return (p(pvar) * dpdeps(pvar) - cs2(pvar) / (ltau * ltau));
}

template <const int size>
double Eos_cons<size>::dpdeps(const var_t &pvar) const
{
  cerr << "ERROR: Eos_cons::dpdeps() is not implemented" << endl
       << "       for the equation of state selected!" << endl << endl;
  abort();
}

/*******************************************************************************
   template <const int size> class Eos_prim
*******************************************************************************/

template <const int size>
inline Eos_prim<size>::Eos_prim(int pdeps_avail) : deps_avail(pdeps_avail)
{
  assert((pdeps_avail == 0) || (pdeps_avail == 2));
}

template <const int size>
double Eos_prim<size>::gamma(const var_t &pvar) const
{
  return (cs2(pvar) / (pvar.tau() * pvar.p()));
}

template <const int size>
double Eos_prim<size>::depsdp(const var_t &pvar) const
{
  cerr << "ERROR: Eos_prim::depsdp() is not implemented" << endl
       << "       for the equation of state selected!" << endl << endl;
  abort();
}

template <const int size>
double Eos_prim<size>::depsdtau(const var_t &pvar) const
{
  cerr << "ERROR: Eos_prim::depsdtau() is not implemented" << endl
       << "       for the equation of state selected!" << endl << endl;
  abort();
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Eos_thin_entry
*******************************************************************************/

inline Eos_thin_entry::~Eos_thin_entry()
{
  if (next_table) 
    delete next_table;
}

/*******************************************************************************
   template <const int size> class Eos_full_entry
*******************************************************************************/

template <const int size>
inline Eos_full_entry<size>::Eos_full_entry(const double (& pval)[size])
: Eos_thin_entry()
{
  for (int i=0;i<size;i++)
    value[i] = pval[i];
}

template <const int size>
inline double Eos_full_entry<size>::get_value(int pidx) const
{
  assert((0 <= pidx) && (pidx < size));

  return value[pidx];
}

/*******************************************************************************
   class Eos_table
*******************************************************************************/

inline double Eos_table::max(double pval1, double pval2) const
{
  return ((pval1 > pval2) ? pval1 : pval2);
}

inline int Eos_table::refine() const
{
  int ret = 0;
 
  if (get_level() < minlevel)
    ret = 1;
  else if (get_level() < maxlevel)
  {
    int i;
    double val_sw[5],val_se[5],val_nw[5],val_ne[5],val_c[5];
    double pnt_sw[2],pnt_se[2],pnt_nw[2],pnt_ne[2],pnt_c[2];
    double lam1,lam2,ival,err=0.0;

    pnt_sw[0] = (first_is_tau ? min1 : 1.0/max1);
    pnt_sw[1] = min2;

    pnt_se[0] = (first_is_tau ? max1 : 1.0/min1);
    pnt_se[1] = min2;

    pnt_nw[0] = (first_is_tau ? min1 : 1.0/max1);
    pnt_nw[1] = max2;

    pnt_ne[0] = (first_is_tau ? max1 : 1.0/min1);
    pnt_ne[1] = max2;

    pnt_c[0]  = (first_is_tau ? 0.5 * (min1 + max1)
                              : 2.0 / (min1 + max1));
    pnt_c[1]  = 0.5 * (min2 + max2);

    table_init(pnt_sw[0],pnt_sw[1],
               val_sw[0],val_sw[1],val_sw[2],val_sw[3],val_sw[4]);
    table_init(pnt_se[0],pnt_se[1],
               val_se[0],val_se[1],val_se[2],val_se[3],val_se[4]);
    table_init(pnt_nw[0],pnt_nw[1],
               val_nw[0],val_nw[1],val_nw[2],val_nw[3],val_nw[4]);
    table_init(pnt_ne[0],pnt_ne[1],
               val_ne[0],val_ne[1],val_ne[2],val_ne[3],val_ne[4]);
    table_init(pnt_c[0], pnt_c[1],
               val_c[0], val_c[1], val_c[2], val_c[3], val_c[4]);

    lam1 = 0.5;
    lam2 = 0.5;

    for (i=0;i<3+additional_values;i++)
    {
      ival =   (val_ne[i] + val_sw[i] - val_se[i] - val_nw[i]) * lam1 * lam2
             + (val_se[i] - val_sw[i]) * lam1
             + (val_nw[i] - val_sw[i]) * lam2
             + val_sw[i];
      assert(val_c[i] > 0.0);
      err = max(err,fabs(ival - val_c[i]) / val_c[i]);
    }

    if (err > tol)
      ret = 1;
  }

  return ret;
}

inline double Eos_table::table(int pcomp,
                               int pidx1,int pidx2) const
{
  double lret;

  assert((0 <= pcomp) && (pcomp < 3+additional_values));
  assert((0 <= pidx1) && (pidx1 < dim1));
  assert((0 <= pidx2) && (pidx2 < dim2));
  assert(entry_type == et_full);

  switch (additional_values)
  {
    case 0:
      lret = ((Eos_full_entry<3> *)eos_entry[pidx1*dim2+pidx2])
             ->get_value(pcomp);
      break;
    case 1:
      lret = ((Eos_full_entry<4> *)eos_entry[pidx1*dim2+pidx2])
             ->get_value(pcomp);
      break;
    case 2:
      lret = ((Eos_full_entry<5> *)eos_entry[pidx1*dim2+pidx2])
             ->get_value(pcomp);
      break;
    default:
      cerr << "ERROR (Eos_table::table): illegal value of"
           << " \"additional_values\"!" << endl;
      abort();
  }

  return lret;
}

/*******************************************************************************
   template <const int size> class Eos_cons_adtab
*******************************************************************************/

template <const int size>
inline int Eos_cons_adtab<size>::get_table_value(int pcomp,
                                                 double ptau,
                                                 double peps,
                                                 double &pret) const
{
  double first = (first_is_tau ? ptau : 1.0/ptau);

  return eos_table->get_table_value(pcomp,first,peps,pret);
}

template <const int size>
inline int Eos_cons_adtab<size>::approximate_dtau(int pcomp,
                                                  double ptau,
                                                  double peps,
                                                  double &pret) const
{
  return eos_table->get_table_derivative(Eos_table::dt_tau,
                                         pcomp,ptau,peps,pret);
}

template <const int size>
inline int Eos_cons_adtab<size>::approximate_deps(int pcomp,
                                                  double ptau,
                                                  double peps,
                                                  double &pret) const
{
  return eos_table->get_table_derivative(Eos_table::dt_second,
                                         pcomp,ptau,peps,pret);
}

template <const int size>
Eos_cons_adtab<size>::Eos_cons_adtab(int pfirst_is_tau,
                                     int pminlevel, int pmaxlevel,
                                     int pdim1, int pdim2,
                                     int paddim1, int paddim2,
                                     double pmin_first, double pmax_first,
                                     double pmin_eps, double pmax_eps,
                                     double ptol, Init_eos &ptable_init)
: Eos_cons<size>(ptable_init.additional_values()),
  first_is_tau(pfirst_is_tau), maxlevel(pmaxlevel), range_warning(1)
{
  assert(pmaxlevel >= pminlevel);
  assert(pminlevel >= 0);
  assert(pdim1 > 1);
  assert(pdim2 > 1);
  assert(paddim1 > 1);
  assert(paddim2 > 1);
  assert(pmin_first < pmax_first);
  assert(pmin_eps   < pmax_eps  );
  assert(ptol >= 0.0);

  eos_table = new Eos_table(first_is_tau,0,pminlevel,pmaxlevel,
                            pdim1,pdim2,paddim1,paddim2,
                            pmin_first,pmax_first,pmin_eps,pmax_eps,
                            ptol,ptable_init);
}

template <const int size>
Eos_cons_adtab<size>::~Eos_cons_adtab()
{
  int i;
  char gnuname[100];

  for (i=0;i<=maxlevel;i++)
  {
    //ostrstream ost(gnuname,100);
    //ost.fill('0');
    //ost << "adtab_cons" << setw(2) << i << ".gnu" << '\0';

    eos_table->tableinfo(i,gnuname);
  }

  delete eos_table;
}

template <const int size>
double Eos_cons_adtab<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2;
  int lflag = get_table_value(0,ltau,leps,lcs2);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::cs2(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_adtab<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lp;
  int lflag = get_table_value(1,ltau,leps,lp);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::p(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(lp > 0.0);

  return lp;
}

template <const int size>
double Eos_cons_adtab<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double ldpdeps;
  int lflag;

  if (this->dpdeps_available())
    lflag = get_table_value(3,ltau,leps,ldpdeps);
  else
    lflag = approximate_deps(1,ltau,leps,ldpdeps);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::dpdeps(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(ldpdeps > 0.0);

  return ldpdeps;
}

template <const int size>
double Eos_cons_adtab<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lT;
  int lflag = get_table_value(2,ltau,leps,lT);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::T(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(lT > 0.0);

  return lT;
}

/*******************************************************************************
   template <const int size> class Eos_prim_adtab
*******************************************************************************/

template <const int size>
inline int Eos_prim_adtab<size>::get_table_value(int pcomp,
                                                 double ptau,
                                                 double pp,
                                                 double &pret) const
{
  double first = (first_is_tau ? ptau : 1.0/ptau);

  return eos_table->get_table_value(pcomp,first,pp,pret);
}

template <const int size>
inline int Eos_prim_adtab<size>::approximate_dp(int pcomp,
                                                double ptau,
                                                double pp,
                                                double &pret) const
{
  return eos_table->get_table_derivative(Eos_table::dt_second,
                                         pcomp,ptau,pp,pret);
}

template <const int size>
inline int Eos_prim_adtab<size>::approximate_dtau(int pcomp,
                                                  double ptau,
                                                  double pp,
                                                  double &pret) const
{
  return eos_table->get_table_derivative(Eos_table::dt_tau,
                                         pcomp,ptau,pp,pret);
}

template <const int size>
Eos_prim_adtab<size>::Eos_prim_adtab(int pfirst_is_tau,
                                     int pminlevel, int pmaxlevel,
                                     int pdim1, int pdim2,
                                     int paddim1, int paddim2,
                                     double pmin_first, double pmax_first,
                                     double pmin_p, double pmax_p,
                                     double ptol, Init_eos &ptable_init)
: Eos_prim<size>(ptable_init.additional_values()),
  first_is_tau(pfirst_is_tau), maxlevel(pmaxlevel), range_warning(1) 
{
  assert(pmaxlevel >= pminlevel);
  assert(pminlevel >= 0);
  assert(pdim1 > 1);
  assert(pdim2 > 1);
  assert(paddim1 > 1);
  assert(paddim2 > 1);
  assert(pmin_first < pmax_first);
  assert(pmin_p     < pmax_p    );
  assert(ptol >= 0.0);

  eos_table = new Eos_table(first_is_tau,0,pminlevel,pmaxlevel,
                            pdim1,pdim2,paddim1,paddim2,
                            pmin_first,pmax_first,pmin_p,pmax_p,
                            ptol,ptable_init);
}

template <const int size>
Eos_prim_adtab<size>::~Eos_prim_adtab()
{
  int i;
  char gnuname[100];

  for (i=0;i<=maxlevel;i++)
  {
    //ostrstream ost(gnuname,100);
    //ost.fill('0');
    //ost << "adtab_prim" << setw(2) << i << ".gnu" << '\0';

    eos_table->tableinfo(i,gnuname);
  }

  delete eos_table;
}

template <const int size>
double Eos_prim_adtab<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2;
  int lflag = get_table_value(0,ltau,lp,lcs2);

  if ((lflag) && (range_warning))
    cerr << "Eos_prim_adtab::cs2(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_adtab<size>::eps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double leps;
  int lflag = get_table_value(1,ltau,lp,leps);

  if ((lflag) && (range_warning))
    cerr << "Eos_prim_adtab::eps(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(leps > 0.0);

  return leps;
}

template <const int size>
double Eos_prim_adtab<size>::depsdp(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double ldepsdp;
  int lflag;

  if (this->deps_available())
    lflag = get_table_value(3,ltau,lp,ldepsdp);
  else
    lflag = approximate_dp(1,ltau,lp,ldepsdp);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::depsdp(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(ldepsdp > 0.0);

  return ldepsdp;
}

template <const int size>
double Eos_prim_adtab<size>::depsdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double ldepsdtau;
  int lflag;

  if (this->deps_available())
    lflag = get_table_value(4,ltau,lp,ldepsdtau);
  else
    lflag = approximate_dtau(1,ltau,lp,ldepsdtau);

  if ((lflag) && (range_warning))
    cerr << "Eos_cons_adtab::depsdtau(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(ldepsdtau > 0.0);
  return ldepsdtau;
}

template <const int size>
double Eos_prim_adtab<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lT;
  int lflag = get_table_value(2,ltau,lp,lT);

  if ((lflag) && (range_warning))
    cerr << "Eos_prim_adtab::T(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range on level "
         << lflag-1 << ")!" << endl;

  assert(lT > 0.0);

  return lT;
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Wrap_eos_init_conseos
*******************************************************************************/

template <const int size>
inline Wrap_eos_init_conseos<size>
::Wrap_eos_init_conseos(eos_t *peos, conv_eos_t *pconv_eos)
: Init_eos(peos->dpdeps_available()), eos(peos), conv_eos(pconv_eos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int size>
void Wrap_eos_init_conseos<size>::operator()(double ptau, double peps,
                                             double &pcs2, double &pp,
                                             double &pT, double &pdpdeps,
                                             double &pdummy) const
{
  cons_t cons(eos,conv_eos);

  for (int i=0;i<size;i++)
    cons[i] = 0.0;

  assert(ptau > 0.0);

  cons[0]      = 1.0 /ptau;
  cons[size-1] = peps/ptau;

  pcs2 = cons.cs2();
  pp   = cons.p();
  pT   = cons.T();

  if (eos->dpdeps_available())
    pdpdeps = eos->dpdeps(cons);
  else
    pdpdeps = -1.0;

  pdummy = -1.0;
}

/*******************************************************************************
   class Wrap_eos_init_primeos
*******************************************************************************/

template <const int size>
inline Wrap_eos_init_primeos<size>
::Wrap_eos_init_primeos(eos_t *peos,conv_eos_t *pconv_eos)
: Init_eos(peos->deps_available()), eos(peos), conv_eos(pconv_eos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int size>
void Wrap_eos_init_primeos<size>::operator()(double ptau, double pp,
                                             double &pcs2, double &peps,
                                             double &pT, double &pdepsdp,
                                             double &pdepsdtau) const
{
  prim_tp_t prim(eos,conv_eos);

  for (int i=0;i<size;i++)
    prim[i] = 0.0;

  prim[0]      = ptau;
  prim[size-1] = pp;

  pcs2 = prim.cs2();
  peps = prim.eps();
  pT   = prim.T();

  if (eos->deps_available())
  {
    pdepsdp   = eos->depsdp(prim);
    pdepsdtau = eos->depsdtau(prim);
  }
  else
  {
    pdepsdp   = -1.0;
    pdepsdtau = -1.0;
  }
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_file
*******************************************************************************/

template <const int size>
Eos_cons_file<size>::Eos_cons_file(const char* pfname) : range_warning(1)
{
  int i,j,k;
  ifstream in(pfname);

  if (!in)
  {
    cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
         << " (could not open file \""
         << pfname << "\")!" << endl;
    abort();
  }
  else
  {
    char s[100];

    cout << "Eos_cons_file::Eos_cons_file(): reading file \""
         << pfname << "\"." << endl;

    in >> s;
    if ((strcmp(s,"rho") != 0) && (strcmp(s,"tau") != 0))
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (only tables in rho or tau and eps possible)!"
           << endl;
      in.close();
      abort();
    }
    else
      if (strcmp(s,"rho") == 0)
        first_is_tau = 0;
      else
        first_is_tau = 1;

    in >> min_first;
    in >> max_first;
    in >> dim_first;

    assert(min_first > 0.0);
    assert(max_first > min_first);
    assert(dim_first > 1);

    step_first = (max_first - min_first) / ((double)dim_first - 1.0);

    in >> s;
    if (strcmp(s,"eps") != 0)
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (only tables in rho or tau and eps possible)!"
           << endl;
      in.close();
      abort();
    }

    in >> min_eps;
    in >> max_eps;
    in >> dim_eps;

    assert(min_eps > 0.0);
    assert(max_eps > min_eps);
    assert(dim_eps > 1);

    step_eps = (max_eps - min_eps) / ((double)dim_eps - 1.0);

    in >> s;
    if (strcmp(s,"cs2") != 0)
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (first table value has to be cs2)!"
           << endl;
      in.close();
      abort();
    }

    in >> s;
    if (strcmp(s,"p") != 0)
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (second table value has to be p)!"
           << endl;
      in.close();
      abort();
    }

    in >> s;
    if (strcmp(s,"T") != 0)
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (third table value has to be T)!"
           << endl;
      in.close();
      abort();
    }

    for (i=0;i<3;i++)
      eostable[i] = new double [dim_first*dim_eps];

    for (i=0;i<dim_first;i++)
      for (j=0;j<dim_eps;j++)
        for (k=0;k<3;k++)
          in >> eostable[k][i*dim_eps+j];

    if (!in.good())
    {
      cerr << "Eos_cons_file::Eos_cons_file(): ERROR"
           << " (unexpected end of file \""
           << pfname << "\")!" << endl;
      in.close();
      abort();
    }
  }
}

template <const int size>
inline double Eos_cons_file<size>::table(int pcomp,
                                         int pidx_first,int pidx_eps) const
{
  assert((0 <= pcomp)      && (pcomp      < 3        ));
  assert((0 <= pidx_first) && (pidx_first < dim_first));
  assert((0 <= pidx_eps  ) && (pidx_eps   < dim_eps  ));

  return eostable[pcomp][pidx_first*dim_eps+pidx_eps];
}

template <const int size>
inline int Eos_cons_file<size>::get_table_value(int pcomp,
                                                double ptau,
                                                double peps,
                                                double &pret) const
{
  double lam_first,lam_eps,val_sw,val_se,val_nw,val_ne;
  double first  = ((first_is_tau) ? ptau : (1.0/ptau));
  int idx_first = (int)((first - min_first) / step_first);
  int idx_eps   = (int)((peps  - min_eps  ) / step_eps  );
  int ret = 1;

  // possibly hit regions (table covers region 0):
  //
  //        6  |:   1   :|  5
  //           |:.......:|
  //        ---|---------|---
  //        '''|:''''''''|:''
  //        4  |:   0   :|: 2 
  //        ...|:........|:..
  //        ---|---------|---
  //           |:''''''':|
  //        7  |:   3   :|  8
  //

  if (   (idx_first < 0) || (idx_first > dim_first-1)
      || (idx_eps   < 0) || (idx_eps   > dim_eps-1  ))
  {
    // value not in table range
    ret = 0;
  }

  if (idx_first < 0)                 // ranges 7,3,8
  {
    idx_first = 0;
    if (idx_eps < 0)                    // range 7
      idx_eps = 0;
    else if (idx_eps >= dim_eps-1)      // range 8
      idx_eps = dim_eps-2;
  }
  else if (idx_first >= dim_first-1) // ranges 6,1,5
  {
    idx_first = dim_first-2;
    if (idx_eps < 0)                    // range 6
      idx_eps = 0;
    else if (idx_eps >= dim_eps-1)      // range 5
      idx_eps = dim_eps-2;
  }
  else if (idx_eps < 0)              // range 4
    idx_eps = 0;
  else if (idx_eps >= dim_eps-1)     // range 2
    idx_eps = dim_eps-2;

  lam_first = (first - (double)idx_first * step_first - min_first) / step_first;
  lam_eps   = (peps  - (double)idx_eps   * step_eps   - min_eps  ) / step_eps;

  val_sw = table(pcomp,idx_first,idx_eps);
  val_se = table(pcomp,idx_first+1,idx_eps);
  val_nw = table(pcomp,idx_first,idx_eps+1);
  val_ne = table(pcomp,idx_first+1,idx_eps+1);

  pret =   (val_ne + val_sw - val_se - val_nw) * lam_first * lam_eps
         + (val_se - val_sw) * lam_first
         + (val_nw - val_sw) * lam_eps
         + val_sw;

  return ret;  
}

template <const int size>
double Eos_cons_file<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2;

  if ((!get_table_value(0,ltau,leps,lcs2)) && (range_warning))
    cerr << "Eos_cons_file::cs2(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range)!"
         << endl;

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_file<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lp;

  if ((!get_table_value(1,ltau,leps,lp)) && (range_warning))
    cerr << "Eos_cons_file::p(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range)!"
         << endl;

  assert(lp > 0.0);

  return lp;
}

template <const int size>
double Eos_cons_file<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lT;

  if ((!get_table_value(2,ltau,leps,lT)) && (range_warning))
    cerr << "Eos_cons_file::T(): WARNING (value [tau,eps] = [" 
         << ltau << "," << leps << "] not in table range)!"
         << endl;

  assert(lT > 0.0);

  return lT;
}

/*******************************************************************************
   template <const int size> class Eos_prim_file
*******************************************************************************/

template <const int size>
Eos_prim_file<size>::Eos_prim_file(const char* pfname) : range_warning(1)
{
  int i,j,k;
  ifstream in(pfname);

  if (!in)
  {
    cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
         << " (could not open file \""
         << pfname << "\")!" << endl;
    abort();
  }
  else
  {
    char s[100];

    cout << "Eos_prim_file::Eos_prim_file(): reading file \""
         << pfname << "\"." << endl;

    in >> s;
    if ((strcmp(s,"rho") != 0) && (strcmp(s,"tau") != 0))
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (only tables in rho or tau and p possible)!"
           << endl;
      in.close();
      abort();
    }
    else
      if (strcmp(s,"rho") == 0)
        first_is_tau = 0;
      else
        first_is_tau = 1;

    in >> min_first;
    in >> max_first;
    in >> dim_first;

    assert(min_first > 0.0);
    assert(max_first > min_first);
    assert(dim_first > 1);

    step_first = (max_first - min_first) / ((double)dim_first - 1.0);

    in >> s;
    if (strcmp(s,"p") != 0)
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (only tables in rho or tau and p possible)!"
           << endl;
      in.close();
      abort();
    }

    in >> min_p;
    in >> max_p;
    in >> dim_p;

    assert(min_p > 0.0);
    assert(max_p > min_p);
    assert(dim_p > 1);

    step_p = (max_p - min_p) / ((double)dim_p - 1.0);

    in >> s;
    if (strcmp(s,"cs2") != 0)
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (first table value has to be cs2)!"
           << endl;
      in.close();
      abort();
    }

    in >> s;
    if (strcmp(s,"eps") != 0)
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (second table value has to be eps)!"
           << endl;
      in.close();
      abort();
    }

    in >> s;
    if (strcmp(s,"T") != 0)
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (third table value has to be T)!"
           << endl;
      in.close();
      abort();
    }

    for (i=0;i<3;i++)
      eostable[i] = new double [dim_first*dim_p];

    for (i=0;i<dim_first;i++)
      for (j=0;j<dim_p;j++)
        for (k=0;k<3;k++)
          in >> eostable[k][i*dim_p+j];

    if (!in.good())
    {
      cerr << "Eos_prim_file::Eos_prim_file(): ERROR"
           << " (unexpected end of file \""
           << pfname << "\")!" << endl;
      in.close();
      abort();
    }
  }
}

template <const int size>
inline double Eos_prim_file<size>::table(int pcomp,
                                         int pidx_first,int pidx_p) const
{
  assert((0 <= pcomp)      && (pcomp      < 3        ));
  assert((0 <= pidx_first) && (pidx_first < dim_first));
  assert((0 <= pidx_p    ) && (pidx_p     < dim_p    ));

  return eostable[pcomp][pidx_first*dim_p+pidx_p];
}

template <const int size>
inline int Eos_prim_file<size>::get_table_value(int pcomp,
                                                double ptau,
                                                double pp,
                                                double &pret) const
{
  double lam_first,lam_p,val_sw,val_se,val_nw,val_ne;
  double first  = ((first_is_tau) ? ptau : (1.0/ptau));
  int idx_first = (int)((first - min_first) / step_first);
  int idx_p     = (int)((pp    - min_p    ) / step_p  );
  int ret = 1;

  // possibly hit regions (table covers region 0):
  //
  //        6  |:   1   :|  5
  //           |:.......:|
  //        ---|---------|---
  //        '''|:''''''''|:''
  //        4  |:   0   :|: 2 
  //        ...|:........|:..
  //        ---|---------|---
  //           |:''''''':|
  //        7  |:   3   :|  8
  //

  if (   (idx_first < 0) || (idx_first > dim_first-1)
      || (idx_p     < 0) || (idx_p     > dim_p-1    ))
  {
    // value not in table range
    ret = 0;
  }

  if (idx_first < 0)                 // ranges 7,3,8
  {
    idx_first = 0;
    if (idx_p < 0)                      // range 7
      idx_p = 0;
    else if (idx_p >= dim_p-1)          // range 8
      idx_p = dim_p-2;
  }
  else if (idx_first >= dim_first-1) // ranges 6,1,5
  {
    idx_first = dim_first-2;
    if (idx_p < 0)                      // range 6
      idx_p = 0;
    else if (idx_p >= dim_p-1)          // range 5
      idx_p = dim_p-2;
  }
  else if (idx_p < 0)                   // range 4
    idx_p = 0;
  else if (idx_p >= dim_p-1)            // range 2
    idx_p = dim_p-2;

  lam_first = (first - (double)idx_first * step_first - min_first) / step_first;
  lam_p     = (pp    - (double)idx_p     * step_p     - min_p    ) / step_p;

  val_sw = table(pcomp,idx_first,idx_p);
  val_se = table(pcomp,idx_first+1,idx_p);
  val_nw = table(pcomp,idx_first,idx_p+1);
  val_ne = table(pcomp,idx_first+1,idx_p+1);

  pret =   (val_ne + val_sw - val_se - val_nw) * lam_first * lam_p
         + (val_se - val_sw) * lam_first
         + (val_nw - val_sw) * lam_p
         + val_sw;

  return ret;  
}

template <const int size>
double Eos_prim_file<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2;

  if ((!get_table_value(0,ltau,lp,lcs2)) && (range_warning))
    cerr << "Eos_prim_file::cs2(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range)!"
         << endl;

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_file<size>::eps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double leps;

  if ((!get_table_value(1,ltau,lp,leps)) && (range_warning))
    cerr << "Eos_prim_file::eps(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range)!"
         << endl;

  assert(leps > 0.0);

  return leps;
}

template <const int size>
double Eos_prim_file<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lT;

  if ((!get_table_value(2,ltau,lp,lT)) && (range_warning))
    cerr << "Eos_prim_file::T(): WARNING (value [tau,p] = [" 
         << ltau << "," << lp << "] not in table range)!"
         << endl;

  assert(lT > 0.0);

  return lT;
}
/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_func
*******************************************************************************/

template <const int size>
inline Eos_cons_func<size>::Eos_cons_func(Init_eos *pinit_eos)
: Eos_cons<size>(pinit_eos->additional_values()),  init_eos(pinit_eos) {}

template <const int size>
double Eos_cons_func<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2,lp,lT,ldpdeps,ldummy;

  (*init_eos)(ltau,leps,lcs2,lp,lT,ldpdeps,ldummy);

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_func<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2,lp,lT,ldpdeps,ldummy;

  (*init_eos)(ltau,leps,lcs2,lp,lT,ldpdeps,ldummy);

  assert(lp > 0.0);

  return lp;
}

template <const int size>
double Eos_cons_func<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2,lp,lT,ldpdeps,ldummy;

  if (this->dpdeps_available())
  {
    (*init_eos)(ltau,leps,lcs2,lp,lT,ldpdeps,ldummy);

    assert(ldpdeps > 0.0);
  }
  else
  {
    cerr << "ERROR (Eos_cons_func::dpdeps): dpdeps is not available!" << endl;
    abort();
  }

  return ldpdeps;
}

template <const int size>
double Eos_cons_func<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lcs2,lp,lT,ldpdeps,ldummy;

  (*init_eos)(ltau,leps,lcs2,lp,lT,ldpdeps,ldummy);

  assert(lT > 0.0);

  return lT;
}

/*******************************************************************************
   template <const int size> class Eos_prim_func
*******************************************************************************/

template <const int size>
inline Eos_prim_func<size>::Eos_prim_func(Init_eos *pinit_eos)
: Eos_prim<size>(pinit_eos->additional_values()), init_eos(pinit_eos) {}

template <const int size>
double Eos_prim_func<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2,leps,lT,ldepsdp,ldepsdtau;

  (*init_eos)(ltau,lp,lcs2,leps,lT,ldepsdp,ldepsdtau);

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_func<size>::eps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2,leps,lT,ldepsdp,ldepsdtau;

  (*init_eos)(ltau,lp,lcs2,leps,lT,ldepsdp,ldepsdtau);

  assert(leps > 0.0);

  return leps;
}

template <const int size>
double Eos_prim_func<size>::depsdp(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2,leps,lT,ldepsdp,ldepsdtau;

  if (this->deps_available())
  {
    (*init_eos)(ltau,lp,lcs2,leps,lT,ldepsdp,ldepsdtau);

    assert(ldepsdp > 0.0);
  }
  else
  {
    cerr << "ERROR (Eos_prim_func::depsdp): depsdp is not available!" << endl;
    abort();
  }

  return ldepsdp;
}

template <const int size>
double Eos_prim_func<size>::depsdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2,leps,lT,ldepsdp,ldepsdtau;

  if (this->deps_available())
  {
    (*init_eos)(ltau,lp,lcs2,leps,lT,ldepsdp,ldepsdtau);

    assert(ldepsdtau > 0.0);
  }
  else
  {
    cerr << "ERROR (Eos_prim_func::depsdtau): depsdtau is not available!"
         << endl;
    abort();
  }

  return ldepsdtau;
}

template <const int size>
double Eos_prim_func<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lp   = pvar.p();
  double lcs2,leps,lT,ldepsdp,ldepsdtau;

  (*init_eos)(ltau,lp,lcs2,leps,lT,ldepsdp,ldepsdtau);

  assert(lT > 0.0);

  return lT;
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_ideal
*******************************************************************************/

template <const int size>
inline Eos_cons_ideal<size>::Eos_cons_ideal(double pgamma, double pR)
: Eos_cons<size>(1), gamma(pgamma), R(pR)
{
  assert(gamma > 1.0);
  assert(R     > 0.0);
}

template <const int size>
double Eos_cons_ideal<size>::cs2(const var_t &pvar) const
{
  double lp   = p(pvar);
  double lcs2 = gamma * lp * pvar.tau();

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_ideal<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();

  assert(ltau > 0.0);
  assert(leps > 0.0);

  return ( (gamma - 1.0) * leps / ltau );
}

template <const int size>
double Eos_cons_ideal<size>::dpdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();

  assert(ltau > 0.0);

  return ( - (gamma - 1.0) * leps / (ltau * ltau) );
}

template <const int size>
double Eos_cons_ideal<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();

  assert(ltau > 0.0);

  return ( (gamma - 1.0) / ltau );
}

template <const int size>
double Eos_cons_ideal<size>::T(const var_t &pvar) const
{
  return ( pvar.tau() * pvar.p() / R );
}

/*******************************************************************************
   template <const int size> class Eos_prim_ideal
*******************************************************************************/

template <const int size>
inline Eos_prim_ideal<size>::Eos_prim_ideal(double pgamma, double pR)
: Eos_prim<size>(2), gamma(pgamma), R(pR)
{
  assert(gamma > 1.0);
  assert(R     > 0.0);
}

template <const int size>
double Eos_prim_ideal<size>::cs2(const var_t &pvar) const
{
  double lcs2 = gamma * pvar.p() * pvar.tau();

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_ideal<size>::eps(const var_t &pvar) const
{
  return ( pvar.p() * pvar.tau() / (gamma - 1.0) );  
}

template <const int size>
double Eos_prim_ideal<size>::depsdp(const var_t &pvar) const
{
  return ( pvar.tau() / (gamma - 1.0) );
}

template <const int size>
double Eos_prim_ideal<size>::depsdtau(const var_t &pvar) const
{
  return ( pvar.p() / (gamma - 1.0) );
}

template <const int size>
double Eos_prim_ideal<size>::dcs2dp(const var_t &pvar) const
{
  return ( gamma * pvar.tau() );
}

template <const int size>
double Eos_prim_ideal<size>::dcs2dtau(const var_t &pvar) const
{
  return ( gamma * pvar.p() );
}

template <const int size>
double Eos_prim_ideal<size>::T(const var_t &pvar) const
{
  return ( pvar.tau() * pvar.p() / R );
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_nconv
*******************************************************************************/

template <const int size>
inline double Eos_cons_nconv<size>::dpdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lc   = (1.0 - gamma)/(ltau * ltau);
  double lret;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = (((5.0 * c5 * ltau + 4.0 * c4)
                             * ltau + 3.0 * c3)
                             * ltau + 2.0 * c2)
                             * ltau + c1;
    lret = (alpha * lnc + (1.0 - alpha) * lc) * pvar.eps();
  }
  else
  {
    lret = lc * pvar.eps();
  }
  
  return lret;
}

template <const int size>
inline double Eos_cons_nconv<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lc   = (gamma - 1.0)/ltau;
  double lret;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = ((((c5 * ltau + c4)
                        * ltau + c3)
                        * ltau + c2)
                        * ltau + c1)
                        * ltau + c0;
    lret = (alpha * lnc + (1.0 - alpha) * lc);
  }
  else
  {
    lret = lc;
  }
  
  return lret;
}

template <const int size>
inline Eos_cons_nconv<size>::Eos_cons_nconv(double palpha)
: alpha(palpha), gamma(1.4), mintaunc(1.0), maxtaunc(5.5),
  c0(1.418181816), c1(-1.894214877), c2(1.175957924),
  c3(-0.3444027039), c4(0.04688204343), c5(-0.002404207354)
{
  assert((0.0 <= alpha) && (alpha <= 1.0));
}

template <const int size>
double Eos_cons_nconv<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();

  return (- ltau * ltau * (dpdtau(pvar) - p(pvar) * dpdeps(pvar)));
}

template <const int size>
double Eos_cons_nconv<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lc   = (gamma - 1.0) / ltau;
  double lret;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = ((((c5 * ltau + c4)
                        * ltau + c3)
                        * ltau + c2)
                        * ltau + c1)
                        * ltau + c0;
    lret = (alpha * lnc + (1.0 - alpha) * lc) * leps;
  }
  else
  {
    lret = lc * leps;
  }
  
  return lret;
}

/*******************************************************************************
   template <const int size> class Eos_prim_nconv
*******************************************************************************/

template <const int size>
inline double Eos_prim_nconv<size>::dpdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lc   = (1.0 - gamma)/(ltau * ltau);
  double lret;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = (((5.0 * c5 * ltau + 4.0 * c4)
                             * ltau + 3.0 * c3)
                             * ltau + 2.0 * c2)
                             * ltau + c1;
    lret = (alpha * lnc + (1.0 - alpha) * lc) * eps(pvar);
  }
  else
  {
    lret = lc * eps(pvar);
  }
  
  return lret;
}

template <const int size>
inline double Eos_prim_nconv<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lc   = (gamma - 1.0)/ltau;
  double lret;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = ((((c5 * ltau + c4)
                        * ltau + c3)
                        * ltau + c2)
                        * ltau + c1)
                        * ltau + c0;
    lret = (alpha * lnc + (1.0 - alpha) * lc);
  }
  else
  {
    lret = lc;
  }
  
  return lret;
}

template <const int size>
inline Eos_prim_nconv<size>::Eos_prim_nconv(double palpha)
: alpha(palpha), gamma(1.4), mintaunc(1.0), maxtaunc(5.5),
  c0(1.418181816), c1(-1.894214877), c2(1.175957924),
  c3(-0.3444027039), c4(0.04688204343), c5(-0.002404207354)
{
  assert((0.0 <= alpha) && (alpha <= 1.0));
}

template <const int size>
double Eos_prim_nconv<size>::cs2(const var_t &pvar) const
{
  double ltau = pvar.tau();

  return (- ltau * ltau * (dpdtau(pvar) - pvar.p() * dpdeps(pvar)));
}

template <const int size>
double Eos_prim_nconv<size>::eps(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double lc   = (gamma - 1.0) / ltau;
  double ldiv;

  if ((mintaunc <= ltau) && (ltau <= maxtaunc))
  {
    double lnc = ((((c5 * ltau + c4)
                        * ltau + c3)
                        * ltau + c2)
                        * ltau + c1)
                        * ltau + c0;
    ldiv = alpha * lnc + (1.0 - alpha) * lc;
  }
  else
  {
    ldiv = lc;
  }

  assert(fabs(ldiv) > 1e-12);

  return (pvar.p() / ldiv);
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_osborne
*******************************************************************************/

template <const int size>
inline Eos_cons_osborne<size>::Eos_cons_osborne()
: tau0(100.0), a1(3.84e-4), a2(1.756e-3), b0(1.312e-2), b1(6.265e-2),
  b2(0.2133), c0(0.5132), c1(0.6761), phi0(2.0e-2)
{}

template <const int size>
double Eos_cons_osborne<size>::cs2(const var_t &pvar) const
{
  double lp   = p(pvar);
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lzet = tau0 / ltau - 1.0;
  double let0 = leps / tau0;
  double lpdt = - tau0 / (ltau * ltau * (let0 + phi0))
                  * (  (a1 + 2.0 * a2 * lzet)
                      + let0 * (b1 + 2.0 * b2 * lzet + let0 * c1));
  double lpde =   (- lp + b0 + lzet * (b1 + b2 * lzet)
                   + 2.0 * let0 * (c0 + c1 * lzet))
                / (leps + tau0 * phi0);
  double lcs2 = - ltau * ltau * (lpdt - lp * lpde);

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_osborne<size>::p(const var_t &pvar) const
{
  double lzet = tau0 / pvar.tau() - 1.0;
  double let0 = pvar.eps() / tau0;
  double lp   =   (  lzet * (a1 + a2 * lzet)
                   + let0 * (b0 + lzet * (b1 + b2 * lzet)
                              + let0 * (c0 + c1 * lzet)))
                / (let0 + phi0);

  return lp;
}

/*******************************************************************************
   template <const int size> class Eos_prim_osborne
*******************************************************************************/

template <const int size>
inline Eos_prim_osborne<size>::Eos_prim_osborne()
: tau0(100.0), a1(3.84e-4), a2(1.756e-3), b0(1.312e-2), b1(6.265e-2),
  b2(0.2133), c0(0.5132), c1(0.6761), phi0(2.0e-2)
{}

template <const int size>
double Eos_prim_osborne<size>::cs2(const var_t &pvar) const
{
  double lp   = pvar.p();
  double ltau = pvar.tau();
  double leps = eps(pvar);
  double lzet = tau0 / ltau - 1.0;
  double let0 = leps / tau0;
  double lpdt = - tau0 / (ltau * ltau * (let0 + phi0))
                  * (  (a1 + 2.0 * a2 * lzet)
                      + let0 * (b1 + 2.0 * b2 * lzet + let0 * c1));
  double lpde =   (- lp + b0 + lzet * (b1 + b2 * lzet)
                   + 2.0 * let0 * (c0 + c1 * lzet))
                / (leps + tau0 * phi0);
  double lcs2 = - ltau * ltau * (lpdt - lp * lpde);

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_osborne<size>::eps(const var_t &pvar) const
{
  double lp   = pvar.p();
  double lzet = tau0 / pvar.tau() - 1.0;
  double lxi  = b0 + lzet * (b1 + b2 * lzet);
  double lsqr = sqrt(  (lxi - lp) * (lxi - lp)
                     - 4.0 * (c0 + c1 * lzet)
                           * (lzet * (a1 + a2 * lzet) - phi0 * lp));
  double leps = tau0 * (lp - lxi + lsqr) / (2.0 * (c0 + c1 * lzet));

  assert(leps > 0.0);
  assert(tau0 * (lp - lxi - lsqr) / (2.0 * (c0 + c1 * lzet)) <= 0.0);

  return leps;  
}
/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_tmv
*******************************************************************************/

template <const int size>
inline Eos_cons_tmv<size>::Eos_cons_tmv()
: r(287.086), cvtr(r / 0.4), thetavib(1000.0), alpha(r), newtoneps(1e-8)
{}

template <const int size>
inline double Eos_cons_tmv<size>::f(double pT, double peps) const
{
  return cvtr * pT + (alpha * thetavib / (exp(thetavib/pT) - 1.0)) - peps;
}

template <const int size>
inline double Eos_cons_tmv<size>::df(double pT) const
{
  double le = exp(thetavib/pT);

  return cvtr + (  alpha * thetavib * thetavib * le 
                 / ((le - 1.0) * (le - 1.0) * pT * pT));
}

template <const int size>
double Eos_cons_tmv<size>::cs2(const var_t &pvar) const
{
  double lT = T(pvar);

  return r * lT * (1.0 + r  / df(lT));
}

template <const int size>
double Eos_cons_tmv<size>::p(const var_t &pvar) const
{
  double lT = T(pvar);

  return r * lT / pvar.tau();;
}

template <const int size>
double Eos_cons_tmv<size>::T(const var_t &pvar) const
{
  double leps = pvar.eps();
  double lT = leps / cvtr;

  while (fabs(f(lT,leps)) > newtoneps)
    lT -= f(lT,leps) / df(lT);

  return lT;
}

/*******************************************************************************
   template <const int size> class Eos_prim_tmv
*******************************************************************************/

template <const int size>
inline double Eos_prim_tmv<size>::deps(double pT) const
{
  double le = exp(thetavib/pT);

  return cvtr + (  alpha * thetavib * thetavib * le 
                 / ((le - 1.0) * (le - 1.0) * pT * pT));
}

template <const int size>
inline Eos_prim_tmv<size>::Eos_prim_tmv()
: r(287.086), cvtr(r / 0.4), thetavib(1000.0), alpha(r)
{}

template <const int size>
double Eos_prim_tmv<size>::cs2(const var_t &pvar) const
{
  double lT = pvar.p() * pvar.tau() / r;

  return r * lT * (1.0 + r  / deps(lT));
}

template <const int size>
double Eos_prim_tmv<size>::eps(const var_t &pvar) const
{
  double lT = pvar.p() * pvar.tau() / r;

  return cvtr * lT + alpha * thetavib / (exp(thetavib/lT) - 1.0);   
}

template <const int size>
double Eos_prim_tmv<size>::T(const var_t &pvar) const
{
  return ( pvar.p() * pvar.tau() / r );
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_waals
*******************************************************************************/

template <const int size>
double Eos_cons_waals<size>::cs2(const var_t &pvar) const
{
  double lp   = p(pvar);
  double ltau = pvar.tau();
  double lcs2 = - 2.0 * a / ltau 
                + (lp * ltau * ltau + a) / (ltau - b) * (1.0 + R / cv); 

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_cons_waals<size>::p(const var_t &pvar) const
{
  double ltau = pvar.tau();
  double leps = pvar.eps();
  double lT   = (leps - eps0 + a / ltau) / cv;
  double lp   = R * lT / (ltau - b) - a / (ltau * ltau);

  assert(ltau > 0.0);
  assert(leps > 0.0);

  return lp;
}

template <const int size>
double Eos_cons_waals<size>::dpdtau(const var_t &pvar) const
{
  double ltau = pvar.tau();

  assert(ltau > 0.0);

  return ( - R * a / (cv * (ltau - b) * ltau * ltau)
           - ( pvar.p() + a / (ltau * ltau)) / (ltau - b)
           + 2.0 * a / (ltau * ltau * ltau) );
}

template <const int size>
double Eos_cons_waals<size>::dpdeps(const var_t &pvar) const
{
  double ltau = pvar.tau();

  assert(ltau > 0.0);

  return ( R / (cv * (ltau - b)) );
}

template <const int size>
double Eos_cons_waals<size>::T(const var_t &pvar) const
{
  return ( (pvar.eps() - eps0 + a / pvar.tau()) / cv );
}

/*******************************************************************************
   template <const int size> class Eos_prim_waals
*******************************************************************************/

template <const int size>
double Eos_prim_waals<size>::cs2(const var_t &pvar) const
{
  double lp   = pvar.p();
  double ltau = pvar.tau();
  double lcs2 = - 2.0 * a / ltau 
                + (lp * ltau * ltau + a) / (ltau - b) * (1.0 + R / cv); 

  assert(lcs2 > 0.0);

  return lcs2;
}

template <const int size>
double Eos_prim_waals<size>::eps(const var_t &pvar) const
{
  double lp   = pvar.p();
  double ltau = pvar.tau();
  double lT   = (lp + a / (ltau * ltau)) * (ltau - b) / R;
  double leps = eps0 + cv * lT - a / ltau;

  assert(leps > 0.0);

  return leps;  
}

template <const int size>
double Eos_prim_waals<size>::depsdp(const var_t &pvar) const
{
  return ( cv * (pvar.tau() - b) / R );
}

template <const int size>
double Eos_prim_waals<size>::depsdtau(const var_t &pvar) const
{
  double lp   = pvar.p();
  double ltau = pvar.tau();

  return (   (cv / R) * (   lp + a / (ltau * ltau)
                          - 2.0 * a * (ltau - b) / (ltau * ltau * ltau) )
           + a / (ltau * ltau) );
}

template <const int size>
double Eos_prim_waals<size>::dcs2dp(const var_t &pvar) const
{
  double ltau = pvar.tau();

  return ( (1.0 + R / cv) * ltau * ltau / (ltau - b) );
}

template <const int size>
double Eos_prim_waals<size>::dcs2dtau(const var_t &pvar) const
{
  double lp   = pvar.p();
  double ltau = pvar.tau();

  return (     (1.0 + R / cv)
             * (2.0 * ltau * lp * (ltau - b) - (lp * ltau * ltau + a) )
             / ((ltau - b) * (ltau - b))
           + 2.0 * a / (ltau * ltau) );
}

template <const int size>
double Eos_prim_waals<size>::T(const var_t &pvar) const
{
  double ltau = pvar.tau();

  return ( ((pvar.p() + a / (ltau * ltau)) * (ltau - b)) / R );
}


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

inline double Sun_eos::SIGN(double pa, double pb) const
{
  return ((pb >= 0.0) ? fabs(pa) : -fabs(pa));
}

inline Sun_eos::Sun_eos()
: boltz(1.3807e-16), amass(1.66057e-24), chi(2.178e-11), cii(2.415e+15),
  hmu(1.0078), mu1(1.297), rho0(2.31200e-06), T0(14400.), Rgask(8.314e7),
  mu0(1.057), p0(rho0*T0*Rgask/mu0), eps0(p0/rho0), cs20(eps0),
  TOL(3.0e-8), ITMAX(100)
{}



/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int vsize> class Variable
*******************************************************************************/

template <const int vsize>
inline double& Variable<vsize>::operator[](int pidx)
{
  return val[pidx];
}

template <const int vsize>
inline double Variable<vsize>::operator[](int pidx) const
{
  return val[pidx];
}

template <const int vsize>
inline void Variable<vsize>::write(ostream &pout) const
{
  val.write(pout);  
}

/*******************************************************************************
   template <const int tsize> class Cons
*******************************************************************************/

template <const int tsize>
inline Cons<tsize>::Cons(eos_t *peos, conv_eos_t *pceos)
: eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Cons<tsize>::Cons(conv_eos_t *pceos, eos_t *peos)
: eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Cons<tsize>::Cons(eos_t *peos, conv_eos_t *pceos, const vec_t &pval)
: Variable<tsize>(pval), eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Cons<tsize>::Cons(eos_t *peos, conv_eos_t *pceos, const dvec_t &pval)
: Variable<tsize>(pval), eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline double Cons<tsize>::p() const
{
  double lp = eos->p(*this);

  assert(lp > 0.0);

  return lp;
}

template <const int tsize>
inline void Cons<tsize>::dWdU_t(mat_t &pret) const
{
  mat_t dmat;

  dWdU(dmat);

  for (int i=0;i<pret.size();i++)
    for (int j=0;j<pret.size();j++)
      pret(i,j) = dmat(j,i);
}

/*******************************************************************************
   template <const int tsize> class Prim
*******************************************************************************/

template <const int tsize>
inline Prim<tsize>::Prim(eos_t *peos, conv_eos_t *pceos)
: eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Prim<tsize>::Prim(conv_eos_t *pceos, eos_t *peos)
: eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Prim<tsize>::Prim(eos_t *peos, conv_eos_t *pceos, const vec_t &pval)
: Variable<tsize>(pval), eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

template <const int tsize>
inline Prim<tsize>::Prim(eos_t *peos, conv_eos_t *pceos, const dvec_t &pval)
: Variable<tsize>(pval), eos(peos), conv_eos(pceos)
{
  assert(eos);
  assert(conv_eos);
}

/*******************************************************************************
   template <const int tsize> class Prim_rp
*******************************************************************************/

template <const int tsize>
inline Prim_rp<tsize>::Prim_rp(eos_t *peos, conv_eos_t *pceos, const vec_t &pval)
: Prim<tsize>(peos,pceos,pval)
{
}

template <const int tsize>
inline Prim_rp<tsize>::Prim_rp(eos_t *peos, conv_eos_t *pceos, const dvec_t &pval)
: Prim<tsize>(peos,pceos,pval)
{
}

/*******************************************************************************
   template <const int tsize> class Prim_tp
*******************************************************************************/

template <const int tsize>
inline Prim_tp<tsize>::Prim_tp(eos_t *peos, conv_eos_t *pceos, const vec_t &pval)
: Prim<tsize>(peos,pceos,pval)
{
}

template <const int tsize>
inline Prim_tp<tsize>::Prim_tp(eos_t *peos, conv_eos_t *pceos, const dvec_t &pval)
: Prim<tsize>(peos,pceos,pval)
{
}

/*******************************************************************************
   template <class VAR> class VariableTop
*******************************************************************************/

template <class VAR>
inline VariableTop<VAR>::VariableTop(const var_t &pvar)
: VAR(pvar.get_eos(),pvar.get_conv_eos())
{
  init(pvar);
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(const cvar1_t &pcvar1)
: VAR(pcvar1.get_conv_eos(),pcvar1.get_eos())
{
  init(pcvar1);
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(const cvar2_t &pcvar2)
: VAR(pcvar2.get_conv_eos(),pcvar2.get_eos())
{
  init(pcvar2);
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(eos_t *peos,
                                     conv_eos_t *pceos,
                                     const vec_t &pvec)
: VAR(peos,pceos,pvec)
{
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(eos_t *peos,
                                     conv_eos_t *pceos,
                                     const dvec_t &pdvec)
: VAR(peos,pceos,pdvec)
{
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(eos_t *peos,
                                     conv_eos_t *pceos,
                                     const var_t &pvar)
: VAR(peos,pceos)
{
  init(pvar);
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(eos_t *peos,
                                     conv_eos_t *pceos,
                                     const cvar1_t &pcvar1)
: VAR(peos,pceos)
{
  init(pcvar1);
}

template <class VAR>
inline VariableTop<VAR>::VariableTop(eos_t *peos,
                                     conv_eos_t *pceos,
                                     const cvar2_t &pcvar2)
: VAR(peos,pceos)
{
  init(pcvar2);
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator=(const var_t &pvar)
{
  if (this != &pvar)
    init(pvar);
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator=(const cvar1_t &pcvar1)
{
  init(pcvar1);
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator=(const cvar2_t &pcvar2)
{
  init(pcvar2);
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator+=(const var_t &pvar)
{
  this->val += pvar.val;
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator-=(const var_t &pvar)
{
  this->val -= pvar.val;
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator*=(double pd)
{
  this->val *= pd;
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t& VariableTop<VAR>::operator/=(double pd)
{
  assert(pd != 0.0);
  this->val /= pd;
  return *this;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t VariableTop<VAR>::operator+(const var_t &pvar)
const
{
  var_t ret(*this);
  ret += pvar;
  return ret;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t VariableTop<VAR>::operator-(const var_t &pvar)
const
{
  var_t ret(*this);
  ret -= pvar;
  return ret;
}

template <class VAR>
inline double VariableTop<VAR>::operator*(const var_t &pvar)
const
{
  return (this->val * pvar.val);
}

template <class VAR>
inline typename VariableTop<VAR>::var_t VariableTop<VAR>::operator*(double pd)
const
{
  var_t ret(*this);
  ret *= pd;
  return ret;
}

template <class VAR>
inline typename VariableTop<VAR>::var_t VariableTop<VAR>::operator/(double pd)
const
{
  var_t ret(*this);
  assert(pd != 0.0);
  ret /= pd;
  return ret;
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                            OUTPUT OPERATORS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

template <const int vsize>
inline ostream& operator<<(ostream &pout, const typename Variable<vsize>::var_t &pvar)
{
  pvar.write(pout);
  return pout;
}

template <const int tsize>
inline ostream& operator<<(ostream &pout, const typename Cons<tsize>::var_t &pvar)
{
  pout << "Cons<" << tsize << "> = ";
  pvar.write(pout);
  return pout;
}

template <const int tsize>
inline ostream& operator<<(ostream &pout, const typename Prim_rp<tsize>::var_t &pvar)
{
  pout << "Prim_rp<" << tsize << "> = ";
  pvar.write(pout);
  return pout;
}

template <const int tsize>
inline ostream& operator<<(ostream &pout, const typename Prim_tp<tsize>::var_t &pvar)
{
  pout << "Prim_tp<" << tsize << "> = ";
  pvar.write(pout);
  return pout;
}

template <class VAR>
inline ostream& operator<<(ostream &pout, const typename VariableTop<VAR>::var_t &pvar)
{
  pvar.write(pout);
  return pout;
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class CONS> class Equations
*******************************************************************************/

template <class CONS> 
inline double Equations<CONS>::min(double pd1, double pd2) const
{
  return (pd1 <= pd2) ? pd1 : pd2;
}

template <class CONS>
inline double Equations<CONS>::max(double pd1, double pd2) const
{
  return (pd1 >= pd2) ? pd1 : pd2;
}

template <class CONS>
inline void Equations<CONS>::r_mat(const cons_t &pvar, mat_t &pret) const
{
  mat_t rmat,dmat;
  prim_t pprim(pvar);

  pprim.dUdW(dmat);

  r_mat(pprim,rmat);

  pret = dmat * rmat;
}

template <class CONS>
inline void Equations<CONS>::l_mat(const cons_t &pvar, mat_t &pret) const
{
  mat_t lmat,dmat;
  prim_t pprim(pvar);

  pvar.dWdU(dmat);

  l_mat(pprim,lmat);

  pret = lmat * dmat;
}

template <class CONS>
inline void Equations<CONS>::grad_lambda(int pidx,
                                         const cons_t &pvar,
                                         vec_t &pret) const
{
  mat_t gmat;

  grad_lambda_vec(pvar,gmat);

  pret = gmat.row(pidx);
}

template <class CONS>
inline void Equations<CONS>::r(int pidx, const cons_t &pvar, vec_t &pret) const
{
  mat_t dmat;
  prim_t pprim(pvar);
  vec_t rvec;

  pprim.dUdW(dmat);

  r(pidx,pprim,rvec);

  pret = dmat * rvec;
}

template <class CONS>
inline double Equations<CONS>::lambda(int pidx, const cons_t &pvar) const
{
  prim_t pprim(pvar);

  return lambda(pidx,pprim);
}


/*******************************************************************************
   class Cons_mhd
*******************************************************************************/

inline Cons_mhd::Cons_mhd(eos_t *peos,
                                  conv_eos_t *pceos,
                                  const vec_t &pvec)
: Cons<8>(peos,pceos,pvec)
{
}

inline Cons_mhd::Cons_mhd(eos_t *peos,
                                  conv_eos_t *pceos,
                                  const dvec_t &pdvec)
: Cons<8>(peos,pceos,pdvec)
{
}

inline double Cons_mhd::min(double pd1, double pd2) const
{
  return (pd1 <= pd2) ? pd1 : pd2;
}

inline double Cons_mhd::max(double pd1, double pd2) const
{
  return (pd1 >= pd2) ? pd1 : pd2;
}

inline void Cons_mhd::init(const var_t &pvar)
{
  if (this != &pvar)
    for (int i=0;i<size();i++)
      (*this)[i] = pvar[i];
}

inline void Cons_mhd::init(const cvar1_t &pcvar1)
{
  (*this)[0] = pcvar1.rho();
  (*this)[1] = pcvar1.rhoux();
  (*this)[2] = pcvar1.rhouy();
  (*this)[3] = pcvar1.rhouz();
  (*this)[4] = pcvar1.Bx();
  (*this)[5] = pcvar1.By();
  (*this)[6] = pcvar1.Bz();
  (*this)[7] = pcvar1.rhoe();
}

inline void Cons_mhd::init(const cvar2_t &pcvar2)
{
  (*this)[0] = pcvar2.rho();
  (*this)[1] = pcvar2.rhoux();
  (*this)[2] = pcvar2.rhouy();
  (*this)[3] = pcvar2.rhouz();
  (*this)[4] = pcvar2.Bx();
  (*this)[5] = pcvar2.By();
  (*this)[6] = pcvar2.Bz();
  (*this)[7] = pcvar2.rhoe();
}

inline double Cons_mhd::rho() const
{
  assert(val[0] > 0.0);
  return val[0];
}

inline double Cons_mhd::rhoux() const
{
  return val[1];
}

inline double Cons_mhd::rhouy() const
{
  return val[2];
}

inline double Cons_mhd::rhouz() const
{
  return val[3];
}

inline double Cons_mhd::Bx() const
{
  return val[4];
}

inline double Cons_mhd::By() const
{
  return val[5];
}

inline double Cons_mhd::Bz() const
{
  return val[6];
}

inline double Cons_mhd::rhoe() const
{
  assert(val[7] > 0.0);
  return val[7];
}

inline double Cons_mhd::ux() const
{
  return (rhoux() / rho());
}

inline double Cons_mhd::uy() const
{
  return (rhouy() / rho());
}

inline double Cons_mhd::uz() const
{
  return (rhouz() / rho());
}

inline double Cons_mhd::e() const
{
  return (rhoe() / rho());
}

inline double Cons_mhd::vsqrt() const
{
  return sqrt(max((va2() + cs2()) * (va2() + cs2()) - 4.0 * vax2() * cs2(),
                  0.0));
}

inline double Cons_mhd::vf2() const
{
  return 0.5 * ((va2() + cs2()) + vsqrt());
}

inline double Cons_mhd::vs2() const
{
  return max(0.5 * ((va2() + cs2()) - vsqrt()),0.0);
}

inline double Cons_mhd::vf() const
{
  return sqrt(vf2());
}

inline double Cons_mhd::vs() const
{
  return sqrt(vs2());
}

/*******************************************************************************
   class Prim_rp_mhd
*******************************************************************************/

inline Prim_rp_mhd::Prim_rp_mhd(eos_t *peos,
                                        conv_eos_t *pceos,
                                        const vec_t &pvec)
: Prim_rp<8>(peos,pceos,pvec)
{
}

inline Prim_rp_mhd::Prim_rp_mhd(eos_t *peos,
                                        conv_eos_t *pceos,
                                        const dvec_t &pdvec)
: Prim_rp<8>(peos,pceos,pdvec)
{
}

inline double Prim_rp_mhd::min(double pd1, double pd2) const
{
  return (pd1 <= pd2) ? pd1 : pd2;
}

inline double Prim_rp_mhd::max(double pd1, double pd2) const
{
  return (pd1 >= pd2) ? pd1 : pd2;
}

inline void Prim_rp_mhd::init(const var_t &pvar)
{
  if (this != &pvar)
    for (int i=0;i<size();i++)
      (*this)[i] = pvar[i];
}

inline void Prim_rp_mhd::init(const cvar1_t &pcvar1)
{
  (*this)[0] = pcvar1.rho();
  (*this)[1] = pcvar1.ux();
  (*this)[2] = pcvar1.uy();
  (*this)[3] = pcvar1.uz();
  (*this)[4] = pcvar1.Bx();
  (*this)[5] = pcvar1.By();
  (*this)[6] = pcvar1.Bz();
  (*this)[7] = pcvar1.p();
}

inline void Prim_rp_mhd::init(const cvar2_t &pcvar2)
{
  (*this)[0] = pcvar2.rho();
  (*this)[1] = pcvar2.ux();
  (*this)[2] = pcvar2.uy();
  (*this)[3] = pcvar2.uz();
  (*this)[4] = pcvar2.Bx();
  (*this)[5] = pcvar2.By();
  (*this)[6] = pcvar2.Bz();
  (*this)[7] = pcvar2.p();
}

inline double Prim_rp_mhd::rho() const
{
  assert(val[0] > 0.0);
  return val[0];
}

inline double Prim_rp_mhd::rhoux() const
{
  return (rho() * ux());
}

inline double Prim_rp_mhd::rhouy() const
{
  return (rho() * uy());
}

inline double Prim_rp_mhd::rhouz() const
{
  return (rho() * uz());
}

inline double Prim_rp_mhd::Bx() const
{
  return val[4];
}

inline double Prim_rp_mhd::By() const
{
  return val[5];
}

inline double Prim_rp_mhd::Bz() const
{
  return val[6];
}

inline double Prim_rp_mhd::rhoe() const
{
  return (rho() * (eps() + 0.5 * u2()) + B2() / (8.0 * M_PI));
}

inline double Prim_rp_mhd::ux() const
{
  return val[1];
}

inline double Prim_rp_mhd::uy() const
{
  return val[2];
}

inline double Prim_rp_mhd::uz() const
{
  return val[3];
}

inline double Prim_rp_mhd::e() const
{
  return (rhoe() / rho());
}

inline double Prim_rp_mhd::vsqrt() const
{
  return sqrt(max((va2() + cs2()) * (va2() + cs2()) - 4.0 * vax2() * cs2(),
                  0.0));
}

inline double Prim_rp_mhd::vf2() const
{
  return 0.5 * ((va2() + cs2()) + vsqrt());
}

inline double Prim_rp_mhd::vs2() const
{
  return max(0.5 * ((va2() + cs2()) - vsqrt()),0.0);
}

inline double Prim_rp_mhd::vf() const
{
  return sqrt(vf2());
}

inline double Prim_rp_mhd::vs() const
{
  return sqrt(vs2());
}

/*******************************************************************************
   class Prim_tp_mhd
*******************************************************************************/

inline Prim_tp_mhd::Prim_tp_mhd(eos_t *peos,
                                        conv_eos_t *pceos,
                                        const vec_t &pvec)
: Prim_tp<8>(peos,pceos,pvec)
{
}

inline Prim_tp_mhd::Prim_tp_mhd(eos_t *peos,
                                        conv_eos_t *pceos,
                                        const dvec_t &pdvec)
: Prim_tp<8>(peos,pceos,pdvec)
{
}

inline double Prim_tp_mhd::min(double pd1, double pd2) const
{
  return (pd1 <= pd2) ? pd1 : pd2;
}

inline double Prim_tp_mhd::max(double pd1, double pd2) const
{
  return (pd1 >= pd2) ? pd1 : pd2;
}

inline void Prim_tp_mhd::init(const var_t &pvar)
{
  if (this != &pvar)
    for (int i=0;i<size();i++)
      (*this)[i] = pvar[i];
}

inline void Prim_tp_mhd::init(const cvar1_t &pcvar1)
{
  (*this)[0] = pcvar1.tau();
  (*this)[1] = pcvar1.ux();
  (*this)[2] = pcvar1.uy();
  (*this)[3] = pcvar1.uz();
  (*this)[4] = pcvar1.Bx();
  (*this)[5] = pcvar1.By();
  (*this)[6] = pcvar1.Bz();
  (*this)[7] = pcvar1.p();
}

inline void Prim_tp_mhd::init(const cvar2_t &pcvar2)
{
  (*this)[0] = pcvar2.tau();
  (*this)[1] = pcvar2.ux();
  (*this)[2] = pcvar2.uy();
  (*this)[3] = pcvar2.uz();
  (*this)[4] = pcvar2.Bx();
  (*this)[5] = pcvar2.By();
  (*this)[6] = pcvar2.Bz();
  (*this)[7] = pcvar2.p();
}

inline double Prim_tp_mhd::rho() const
{
  return (1.0 / tau());
}

inline double Prim_tp_mhd::rhoux() const
{
  return (rho() * ux());
}

inline double Prim_tp_mhd::rhouy() const
{
  return (rho() * uy());
}

inline double Prim_tp_mhd::rhouz() const
{
  return (rho() * uz());
}

inline double Prim_tp_mhd::Bx() const
{
  return val[4];
}

inline double Prim_tp_mhd::By() const
{
  return val[5];
}

inline double Prim_tp_mhd::Bz() const
{
  return val[6];
}

inline double Prim_tp_mhd::rhoe() const
{
  return (  rho() * (eps() + 0.5 * u2()) + B2() / (8.0 * M_PI));
}

inline double Prim_tp_mhd::ux() const
{
  return val[1];
}

inline double Prim_tp_mhd::uy() const
{
  return val[2];
}

inline double Prim_tp_mhd::uz() const
{
  return val[3];
}

inline double Prim_tp_mhd::e() const
{
  return (rhoe() / rho());
}

inline double Prim_tp_mhd::vsqrt() const
{
  return sqrt(max((va2() + cs2()) * (va2() + cs2()) - 4.0 * vax2() * cs2(),
                  0.0));
}

inline double Prim_tp_mhd::vf2() const
{
  return 0.5 * ((va2() + cs2()) + vsqrt());
}

inline double Prim_tp_mhd::vs2() const
{
  return max(0.5 * ((va2() + cs2()) - vsqrt()),0.0);
}

inline double Prim_tp_mhd::vf() const
{
  return sqrt(vf2());
}

inline double Prim_tp_mhd::vs() const
{
  return sqrt(vs2());
}

/*******************************************************************************
   class Mhd
*******************************************************************************/

inline void Mhd::r_mat(const cons_t &pvar, mat_t &pret) const
{
  Equations< VariableTop<Cons_mhd> >::r_mat(pvar,pret);
}

inline void Mhd::l_mat(const cons_t &pvar, mat_t &pret) const
{
  Equations< VariableTop<Cons_mhd> >::l_mat(pvar,pret);
}

inline void Mhd::r(int pidx, const cons_t &pvar, vec_t &pret) const
{
  Equations< VariableTop<Cons_mhd> >::r(pidx,pvar,pret);
}

inline double Mhd::lambda(int pidx, const cons_t &pvar) const
{
  return Equations< VariableTop<Cons_mhd> >::lambda(pidx,pvar);
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                            OUTPUT OPERATORS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

inline ostream& operator<<(ostream &pout, const Cons_mhd::var_t &pcons)
{
  pout << "(rho,rhoux,rhouy,rhouz,Bx,By,Bz,rhoe) = ";
  pcons.write(pout);
  return pout;
}

inline ostream& operator<<(ostream &pout, const Prim_rp_mhd::var_t &pprim)
{
  pout << "(rho,ux,uy,uz,Bx,By,Bz,p) = ";
  pprim.write(pout);
  return pout;
}

inline ostream& operator<<(ostream &pout, const Prim_tp_mhd::var_t &pprim)
{
  pout << "(tau,ux,uy,uz,Bx,By,Bz,p) = ";
  pprim.write(pout);
  return pout;
}

} // end namespace Mhd 
#include "mhd_eqns.cc"
#include "mhd_fluxes.cc"

#endif  // MHD_H_INCLUDED


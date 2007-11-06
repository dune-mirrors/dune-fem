#ifndef __HEADER_MHD_FLUX
#define __HEADER_MHD_FLUX
using namespace std;

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Mhd {

template <const int vsize> class Variable;
template <class T, const int msize> class Matrix;

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 *****************************************************************************/

/******************************************************************************
   template <class T, const int vsize> class Vector
******************************************************************************/

template <class T, const int vsize> class MhdVector
{
protected:
  typedef MhdVector <T,vsize> vec_t;
  typedef Matrix <T,vsize> mat_t;

  T vec[vsize];

  inline void init(const T [vsize]);
  inline void init(const vec_t &pvec) {init(pvec.vec);}
public:
  void clear() {for (int i=0;i<vsize;i++) vec[i]=(T)0;}

  MhdVector() {clear();}
  MhdVector(const T pvec[vsize]) {init(pvec);}
  MhdVector(const vec_t &pvec) {init(pvec.vec);}
  int size() const {return vsize;}
  inline vec_t& operator=(const vec_t&);
  inline T& operator[](int);
  inline T operator[](int) const;
  inline vec_t& operator+=(const vec_t&);
  inline vec_t& operator-=(const vec_t&);
  inline vec_t& operator*=(const mat_t&);
  inline vec_t& operator*=(const T);
  inline vec_t& operator/=(const T);
  inline vec_t operator+(const vec_t&) const;
  inline vec_t operator-(const vec_t&) const;
  inline vec_t operator*(const mat_t&) const;
  inline double operator*(const vec_t&) const;
  inline vec_t operator*(const T) const;
  inline vec_t operator/(const T) const;
  inline void write(ostream &) const;

  friend class Variable <vsize>;
};

/*******************************************************************************
   template <class T, const int msize> class Matrix
*******************************************************************************/

template <class T, const int msize> class Matrix
{
protected:
  typedef MhdVector <T,msize> vec_t;
  typedef Matrix <T,msize> mat_t;

  T mat[msize][msize];

  inline void init(const T [msize][msize]);
  inline void init(const vec_t [msize]);
  inline void init(const mat_t &pmat) {init(pmat.mat);}
public:
  inline void clear();

  Matrix() {clear();}
  Matrix(const T pmat[msize][msize]) {init(pmat);}
  Matrix(const vec_t pmat[msize]) {init(pmat);}
  Matrix(const mat_t &pmat) {init(pmat.mat);}
  int size() const {return msize;}
  inline mat_t& operator=(const mat_t&);
  inline T& operator()(int,int);
  inline T operator()(int,int) const;
  inline mat_t& operator+=(const mat_t&);
  inline mat_t& operator-=(const mat_t&);
  inline mat_t& operator*=(const mat_t&);
  inline vec_t  operator*=(const vec_t&);
  inline mat_t& operator*=(const T);
  inline mat_t& operator/=(const T);
  inline mat_t operator+(const mat_t&) const;
  inline mat_t operator-(const mat_t&) const;
  inline mat_t operator*(const mat_t&) const;
  inline vec_t operator*(const vec_t&) const;
  inline mat_t operator*(const T) const;
  inline mat_t operator/(const T) const;
  inline vec_t row(int) const;
  inline vec_t col(int) const;
  inline void write(ostream &) const;
};

template <const int size> class Cons;
template <const int size> class Prim;

#ifdef _IPGAS
extern double global_gamma_id;
#endif

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class VAR> class Eos
*******************************************************************************/

template <class VAR> class Eos
{
protected:
  typedef VAR var_t;

public:
  virtual ~Eos() {}
  virtual double cs2(const var_t &) const = 0;
  double cs(const var_t &pvar) const {return sqrt(cs2(pvar));}
};

/*******************************************************************************
   template <const int size> class Eos_cons
*******************************************************************************/

template <const int size> class Eos_cons : public Eos < Cons<size> >
{
protected:
  typedef typename Eos < Cons<size> > :: var_t var_t;

  const int dpdeps_avail;
public:
  inline Eos_cons(int = 0);
  virtual double p(const var_t &) const = 0;
  double gamma(const var_t &) const;
  int dpdeps_available() const { return dpdeps_avail; }

  /* extensions for cws and conservative solvers */

  virtual double dpdtau(const var_t &) const;
  virtual double dpdeps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const {assert(0 > 1); abort();}
};

/*******************************************************************************
   template <const int size> class Eos_prim
*******************************************************************************/

template <const int size> class Eos_prim : public Eos < Prim<size> >
{
protected:
  typedef typename Eos < Prim<size> > :: var_t var_t;

  const int deps_avail;
public:
  inline Eos_prim(int = 0);
  virtual double eps(const var_t &) const = 0;
  double gamma(const var_t &) const;
  int deps_available() const { return deps_avail; }

  /* extensions for cws and conservative solvers */

  virtual double depsdp(const var_t &)   const;
  virtual double depsdtau(const var_t &) const;

  /* extensions for cws */

  virtual double dcs2dp(const var_t &)   const {assert(0 > 1); abort();}
  virtual double dcs2dtau(const var_t &) const {assert(0 > 1); abort();}

  /* extensions for 2d applications */

  virtual double T(const var_t &) const {assert(0 > 1); abort();}
};
/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int vsize> class Variable
*******************************************************************************/

template <const int vsize> class Variable
{
protected:
  typedef double dvec_t[vsize];
  typedef Variable<vsize> var_t;
  typedef MhdVector<double,vsize> vec_t;
  typedef Matrix<double,vsize> mat_t;

  vec_t val;

  Variable() {}
  Variable(const vec_t& pvar) {val.init(pvar.vec);}
  Variable(const dvec_t &pval) {val.init(pval);}
public:
  virtual ~Variable() {}
  void clear() {val.clear();}
  void init(const vec_t &pvar) {val.init(pvar.vec);}
  void init(const dvec_t &pval) {val.init(pval);}
  vec_t vec() const {return val;}
  virtual double tau() const = 0;
  int size() const {return vsize;}
  inline double& operator[](int);
  inline double operator[](int) const;
  inline void write(ostream&) const;
};

/*******************************************************************************
   template <const int tsize> class Cons
*******************************************************************************/

template <const int tsize> class Cons : public Variable <tsize>
{
protected:
  typedef Cons <tsize> var_t;
  typedef Eos_cons <tsize> eos_t;
  typedef Eos_prim <tsize> conv_eos_t;
  typedef typename Variable <tsize> :: vec_t vec_t;
  typedef typename Variable <tsize> :: dvec_t dvec_t;
  typedef typename Variable <tsize> :: mat_t mat_t;

  eos_t *eos;
  conv_eos_t *conv_eos;

  Cons() : eos((eos_t*)0), conv_eos((conv_eos_t*)0) {}
  inline Cons(eos_t*, conv_eos_t*);
  inline Cons(conv_eos_t*, eos_t*);
  inline Cons(eos_t*, conv_eos_t*, const vec_t&);
  inline Cons(eos_t*, conv_eos_t*, const dvec_t&);
public:
  eos_t *get_eos() const {return eos;}
  conv_eos_t *get_conv_eos() const {return conv_eos;}
  virtual double eps() const = 0;
  inline double p() const;
  double cs() const {return eos->cs(*this);}
  double cs2() const {return eos->cs2(*this);}
  double gamma() const {return eos->gamma(*this);}

  /* extensions for cws */

  virtual void gradp(vec_t&) const {cout<<"vars.h1\n"<<flush;abort();}
  virtual void dWdU(mat_t&) const {cout<<"vars.h2\n"<<flush;abort();}
  inline void dWdU_t(mat_t&) const;

  /* extensions for 2d applications */

  double T() const {return eos->T(*this);}
};

/*******************************************************************************
   template <const int tsize> class Prim
*******************************************************************************/

template <const int tsize> class Prim : public Variable <tsize>
{
protected:
  typedef Prim <tsize> var_t;
  typedef Eos_prim <tsize> eos_t;
  typedef Eos_cons <tsize> conv_eos_t;
  typedef typename Variable <tsize> :: vec_t vec_t;
  typedef typename Variable <tsize> :: dvec_t dvec_t;
  typedef typename Variable <tsize> :: mat_t mat_t;

  eos_t *eos;
  conv_eos_t *conv_eos;

  Prim() : eos((eos_t*)0), conv_eos((conv_eos_t*)0) {}
  inline Prim(eos_t*, conv_eos_t*);
  inline Prim(conv_eos_t*, eos_t*);
  inline Prim(eos_t*, conv_eos_t*, const vec_t&);
  inline Prim(eos_t*, conv_eos_t*, const dvec_t&);
public:
  eos_t *get_eos() const {return eos;}
  conv_eos_t *get_conv_eos() const {return conv_eos;}
  virtual double p() const = 0;
  double eps() const {return eos->eps(*this);}
  double cs() const {return eos->cs(*this);}
  double cs2() const {return eos->cs2(*this);}
  double gamma() const {return eos->gamma(*this);}

  /* extensions for 2d applications */

  double T() const {return eos->T(*this);}
};

/*******************************************************************************
   template <const int size> class Prim_rp
*******************************************************************************/

template <const int tsize> class Prim_rp : public Prim <tsize>
{
protected :
  typedef typename Prim<tsize> :: vec_t vec_t;
  typedef typename Prim<tsize> :: dvec_t dvec_t;
  typedef typename Prim<tsize> :: mat_t mat_t;
  typedef typename Prim<tsize> :: eos_t eos_t;
  typedef typename Prim<tsize> :: conv_eos_t conv_eos_t;

  Prim_rp() : Prim<tsize>(){}
  Prim_rp(eos_t *peos, conv_eos_t *pceos) : Prim<tsize>(peos,pceos) {}
  Prim_rp(conv_eos_t *pceos, eos_t *peos) : Prim<tsize>(pceos,peos) {}
  inline Prim_rp(eos_t*, conv_eos_t*, const vec_t&);
  inline Prim_rp(eos_t*, conv_eos_t*, const dvec_t&);
};

/*******************************************************************************
   template <const int size> class Prim_tp
*******************************************************************************/

template <const int tsize> class Prim_tp : public Prim <tsize>
{
protected :
  typedef typename Prim<tsize> :: vec_t vec_t;
  typedef typename Prim<tsize> :: dvec_t dvec_t;
  typedef typename Prim<tsize> :: mat_t mat_t;
  typedef typename Prim<tsize> :: eos_t eos_t;
  typedef typename Prim<tsize> :: conv_eos_t conv_eos_t;

  Prim_tp() : Prim<tsize>(){}
  Prim_tp(eos_t *peos, conv_eos_t *pceos) : Prim<tsize>(peos,pceos) {}
  Prim_tp(conv_eos_t *pceos, eos_t *peos) : Prim<tsize>(pceos,peos) {}
  inline Prim_tp(eos_t*, conv_eos_t*, const vec_t&);
  inline Prim_tp(eos_t*, conv_eos_t*, const dvec_t&);

  /* extensions for cws */

public:
  virtual void gradcs2(vec_t&) const {cout<<"vars.h3\n"<<flush; abort();}
  virtual void gradeps(vec_t&) const {cout<<"vars.h4\n"<<flush; abort();}
  virtual void dUdW(mat_t&)    const {cout<<"vars.h5\n"<<flush; abort();}
};

/*******************************************************************************
   template <class VAR> class VariableTop
*******************************************************************************/

template <class VAR> class VariableTop : public VAR
{
protected:
  VariableTop() : VAR() {}
public:
  typedef VariableTop <VAR> var_t;
  typedef VariableTop <typename VAR::cvar1_t> cvar1_t;
  typedef VariableTop <typename VAR::cvar2_t> cvar2_t;
  typedef VariableTop <typename VAR::prim_t>  prim_t;
  typedef typename VAR::vec_t  vec_t;
  typedef typename VAR::mat_t  mat_t;
  typedef typename VAR::dvec_t dvec_t;
  typedef typename VAR::eos_t      eos_t;
  typedef typename VAR::conv_eos_t conv_eos_t;
 
  VariableTop(eos_t *peos, conv_eos_t *pceos) : VAR(peos,pceos) {}
  inline VariableTop(const var_t&);
  inline explicit VariableTop(const cvar1_t&);
  inline explicit VariableTop(const cvar2_t&);
  inline VariableTop(eos_t*, conv_eos_t*, const vec_t&);
  inline VariableTop(eos_t*, conv_eos_t*, const dvec_t&);
  inline VariableTop(eos_t*, conv_eos_t*, const var_t&);
  inline VariableTop(eos_t*, conv_eos_t*, const cvar1_t&);
  inline VariableTop(eos_t*, conv_eos_t*, const cvar2_t&);
  inline var_t& operator=(const var_t&);
  inline var_t& operator=(const cvar1_t&);
  inline var_t& operator=(const cvar2_t&);
  inline var_t& operator+=(const var_t&);
  inline var_t& operator-=(const var_t&);
  inline var_t& operator*=(double);
  inline var_t& operator/=(double);
  inline var_t operator+(const var_t&) const;
  inline var_t operator-(const var_t&) const;
  inline double operator*(const var_t&) const;
  inline var_t operator*(double) const;
  inline var_t operator/(double) const;
};



/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Init_eos
*******************************************************************************/

class Init_eos
{
protected:
  const int add_values;
public:
  Init_eos(int padd_values = 0) : add_values(padd_values) {}
  virtual ~Init_eos() {}
  int additional_values() const { return add_values; }
  virtual void operator()(double,double,
                          double&,double&,double&,double&,double&) const = 0;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class CONS> class Equations
*******************************************************************************/

template <class CONS> class Equations
{
public:
  typedef CONS cons_t;
  typedef typename CONS::cvar1_t cvar1_t;
  typedef typename CONS::cvar2_t cvar2_t;
  typedef typename CONS::prim_t  prim_t;
  typedef typename CONS::vec_t   vec_t;
  typedef typename CONS::mat_t   mat_t;
protected:
  inline double min(double, double) const;
  inline double max(double, double) const;
public:
  virtual ~Equations() {}
  virtual int ewave() const = 0;
  virtual int nonlinear(int) const = 0;
  virtual void f(const cons_t&, cons_t&) const = 0;
  virtual void f(const prim_t&, cons_t&) const = 0;
  virtual double lambda(int, const prim_t&) const = 0;
  virtual void lambda_vec(const prim_t&, vec_t&) const = 0;
  virtual void r(int, const prim_t&, vec_t&) const = 0;
  virtual void r_mat(const prim_t&, mat_t&) const = 0;
  virtual void l(int, const prim_t&, vec_t&) const = 0;
  virtual void l_mat(const prim_t&, mat_t&) const = 0;
  virtual double dt_local(double, double,
                          const cons_t&, const cons_t&) const = 0;
  virtual double indicator(bool, int, const cons_t&, const cons_t&) const = 0;

  /* extensions for flux_dwc */

  inline void r_mat(const cons_t&, mat_t&) const;
  inline void l_mat(const cons_t&, mat_t&) const;

  /* extensions for cws (cws also requires l_mat(const cons_t&, mat_t&)) */

  virtual void Df(const cons_t&, mat_t&) const {assert(0>1); abort();}
  virtual void grad_lambda_vec(const cons_t&, mat_t&) const {assert(0>1); abort();}
  inline void grad_lambda(int, const cons_t&, vec_t&) const;
  inline void r(int, const cons_t&, vec_t&) const;
  inline double lambda(int, const cons_t&) const;
};

class Cons_mhd;
class Prim_rp_mhd;
class Prim_tp_mhd;

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Cons_mhd
*******************************************************************************/

class Cons_mhd : public Cons <8>
{
public:
  typedef Cons_mhd    var_t;
  typedef Prim_rp_mhd cvar1_t;
  typedef Prim_tp_mhd cvar2_t;
  typedef Prim_tp_mhd prim_t;
protected:
  inline double min(double, double) const;
  inline double max(double, double) const;
  inline void init(const var_t&);
  inline void init(const cvar1_t&);
  inline void init(const cvar2_t&);

  Cons_mhd() : Cons<8>(){}
  Cons_mhd(eos_t *peos, conv_eos_t *pceos) : Cons<8>(peos,pceos){}
  Cons_mhd(conv_eos_t *pceos, eos_t *peos) : Cons<8>(pceos,peos){}
  inline Cons_mhd(eos_t*, conv_eos_t*, const vec_t&);
  inline Cons_mhd(eos_t*, conv_eos_t*, const dvec_t&);
public:
  virtual double tau() const;
  virtual double eps() const;
  inline double rho() const;
  inline double rhoux() const;
  inline double rhouy() const;
  inline double rhouz() const;
  inline double Bx() const;
  inline double By() const;
  inline double Bz() const;
  inline double rhoe() const;
  inline double ux() const;
  inline double uy() const;
  inline double uz() const;
  inline double e() const;
  inline double& rho()  {return val[0];}
  inline double& rhoux(){return val[1];}
  inline double& rhouy(){return val[2];}
  inline double& rhouz(){return val[3];}
  inline double& Bx()   {return val[4];}
  inline double& By()   {return val[5];}
  inline double& Bz()   {return val[6];}
  inline double& rhoe() {return val[7];}
  inline double u2() const    {return ux() * ux() + uy() * uy() + uz() * uz();}
  inline double B2() const    {return Bx() * Bx() + By() * By() + Bz() * Bz();}
  inline double Bu() const    {return Bx() * ux() + By() * uy() + Bz() * uz();}
  inline double vax2() const  {return Bx() * Bx() / (4.0 * M_PI * rho());}
  inline double va2() const   {return B2() / (4.0 * M_PI * rho());}
  inline double vsqrt() const;
  inline double vf2() const;
  inline double vs2() const;
  inline double vax() const   {return sqrt(vax2());}
  inline double va() const    {return sqrt(va2());}
  inline double vf() const;
  inline double vs() const;

  /* extensions for cws */

  virtual void gradp(vec_t&) const;
  virtual void dWdU(mat_t&) const;
};

/*******************************************************************************
   class Prim_rp_mhd
*******************************************************************************/

class Prim_rp_mhd : public Prim_rp <8>
{
public:
  typedef Prim_rp_mhd var_t;
  typedef Cons_mhd    cvar1_t;
  typedef Prim_tp_mhd cvar2_t;
  typedef Prim_tp_mhd prim_t;
protected:
  inline double min(double, double) const;
  inline double max(double, double) const;
  inline void init(const var_t&);
  inline void init(const cvar1_t&);
  inline void init(const cvar2_t&);

  Prim_rp_mhd() : Prim_rp<8>(){}
  Prim_rp_mhd(eos_t *peos, conv_eos_t *pceos) : Prim_rp<8>(peos,pceos){}
  Prim_rp_mhd(conv_eos_t *pceos, eos_t *peos) : Prim_rp<8>(pceos,peos){}
  inline Prim_rp_mhd(eos_t*, conv_eos_t*, const vec_t&);
  inline Prim_rp_mhd(eos_t*, conv_eos_t*, const dvec_t&);
public:
  virtual double tau() const;
  virtual double p() const;
  inline double rho() const;
  inline double rhoux() const;
  inline double rhouy() const;
  inline double rhouz() const;
  inline double Bx() const;
  inline double By() const;
  inline double Bz() const;
  inline double rhoe() const;
  inline double ux() const;
  inline double uy() const;
  inline double uz() const;
  inline double e() const;
  inline double& rho() {return val[0];}
  inline double& ux()  {return val[1];}
  inline double& uy()  {return val[2];}
  inline double& uz()  {return val[3];}
  inline double& Bx()  {return val[4];}
  inline double& By()  {return val[5];}
  inline double& Bz()  {return val[6];}
  inline double& p()   {return val[7];}
  inline double u2() const    {return ux() * ux() + uy() * uy() + uz() * uz();}
  inline double B2() const    {return Bx() * Bx() + By() * By() + Bz() * Bz();}
  inline double Bu() const    {return Bx() * ux() + By() * uy() + Bz() * uz();}
  inline double vax2() const  {return Bx() * Bx() / (4.0 * M_PI * rho());}
  inline double va2() const   {return B2() / (4.0 * M_PI * rho());}
  inline double vsqrt() const;
  inline double vf2() const;
  inline double vs2() const;
  inline double vax() const   {return sqrt(vax2());}
  inline double va() const    {return sqrt(va2());}
  inline double vf() const;
  inline double vs() const;
};

/*******************************************************************************
   class Prim_tp_mhd
*******************************************************************************/

class Prim_tp_mhd : public Prim_tp <8>
{
public:
  typedef Prim_tp_mhd var_t;
  typedef Cons_mhd    cvar1_t;
  typedef Prim_rp_mhd cvar2_t;
  typedef Prim_tp_mhd prim_t;
protected:
  inline double min(double, double) const;
  inline double max(double, double) const;
  inline void init(const var_t&);
  inline void init(const cvar1_t&);
  inline void init(const cvar2_t&);

  Prim_tp_mhd() : Prim_tp<8>(){}
  Prim_tp_mhd(eos_t *peos, conv_eos_t *pceos) : Prim_tp<8>(peos,pceos){}
  Prim_tp_mhd(conv_eos_t *pceos, eos_t *peos) : Prim_tp<8>(pceos,peos){}
  inline Prim_tp_mhd(eos_t*, conv_eos_t*, const vec_t&);
  inline Prim_tp_mhd(eos_t*, conv_eos_t*, const dvec_t&);
public:
  virtual double tau() const;
  virtual double p() const;
  inline double rho() const;
  inline double rhoux() const;
  inline double rhouy() const;
  inline double rhouz() const;
  inline double Bx() const;
  inline double By() const;
  inline double Bz() const;
  inline double rhoe() const;
  inline double ux() const;
  inline double uy() const;
  inline double uz() const;
  inline double e() const;
  inline double& tau() {return val[0];}
  inline double& ux()  {return val[1];}
  inline double& uy()  {return val[2];}
  inline double& uz()  {return val[3];}
  inline double& Bx()  {return val[4];}
  inline double& By()  {return val[5];}
  inline double& Bz()  {return val[6];}
  inline double& p()   {return val[7];}
  inline double u2() const    {return ux() * ux() + uy() * uy() + uz() * uz();}
  inline double B2() const    {return Bx() * Bx() + By() * By() + Bz() * Bz();}
  inline double Bu() const    {return Bx() * ux() + By() * uy() + Bz() * uz();}
  inline double vax2() const  {return Bx() * Bx() / (4.0 * M_PI * rho());}
  inline double va2() const   {return B2() / (4.0 * M_PI * rho());}
  inline double vsqrt() const;
  inline double vf2() const;
  inline double vs2() const;
  inline double vax() const   {return sqrt(vax2());}
  inline double va() const    {return sqrt(va2());}
  inline double vf() const;
  inline double vs() const;

  /* extensions for cws */

  virtual void gradcs2(vec_t&) const;
  virtual void gradeps(vec_t&) const;
  virtual void dUdW(mat_t&)    const;
};

/*******************************************************************************
   class Mhd
*******************************************************************************/

class Mhd : public Equations< VariableTop<Cons_mhd> >
{
public:
  typedef VariableTop<Prim_rp_mhd> print_t;
  typedef enum {mi_none=0,mi_rho=1,mi_rhoux=2,mi_rhouy=4,mi_rhouz=8,
                mi_bx=16,mi_by=32,mi_bz=64,mi_rhoe=128,mi_ux=256,
                mi_uy=512,mi_uz=1024,mi_e=2048,mi_p=4096} mindicator_t;
protected:
  const double MhdEPS;
public:
  Mhd() : MhdEPS(1e-12) {}
  virtual ~Mhd() {}
  virtual int ewave() const {return 3;}
  virtual int nonlinear(int) const;
  virtual void f(const cons_t&, cons_t&) const;
  virtual void f(const prim_t&, cons_t&) const;
  virtual double lambda(int, const prim_t&) const;
  virtual void lambda_vec(const prim_t&, vec_t&) const;
  virtual void r(int, const prim_t&, vec_t&) const;
  virtual void r_mat(const prim_t&, mat_t&) const;
  virtual void l(int, const prim_t&, vec_t&) const;
  virtual void l_mat(const prim_t&, mat_t&) const;
  virtual double dt_local(double, double,
                          const cons_t&, const cons_t&) const;
  virtual double indicator(bool, int, const cons_t&, const cons_t&) const;

  /* extensions for flux_dwc */

  inline void r_mat(const cons_t&, mat_t&) const;
  inline void l_mat(const cons_t&, mat_t&) const;

  /* extensions for cws (cws also requires l_mat(const cons_t&, mat_t&)) */

  virtual void Df(const cons_t&, mat_t&) const;
  virtual void grad_lambda_vec(const cons_t&, mat_t&) const;
  inline void r(int, const cons_t&, vec_t&) const;
  inline double lambda(int, const cons_t&) const;
};

typedef VariableTop<Cons_mhd>    cons_t;
typedef VariableTop<Prim_rp_mhd> prim_rp_t;
typedef VariableTop<Prim_tp_mhd> prim_tp_t;

typedef Mhd eqns_t;
typedef Eos_cons<8> eos_cons_t;
typedef Eos_prim<8> eos_prim_t;

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Eos_thin_entry
*******************************************************************************/

class Eos_table;

class Eos_thin_entry
{
protected:
  Eos_table *next_table;

public:
  Eos_thin_entry() : next_table(0) {}
  inline ~Eos_thin_entry();

  Eos_table *get_table() { return next_table; }
  void attach_table(Eos_table *ptable) { next_table = ptable; }
};

/*******************************************************************************
   template <const int size> class Eos_full_entry
*******************************************************************************/

template <const int size> class Eos_full_entry : public Eos_thin_entry
{
protected:
  double value[size];

public:
  inline Eos_full_entry(const double (&)[size]);

  inline double get_value(int) const;
};

/*******************************************************************************
   class Eos_table
*******************************************************************************/

class Eos_table
{
protected:
  typedef enum {et_thin,et_full} entry_type_t;

  entry_type_t entry_type;
  const int additional_values;
  int first_is_tau,level,minlevel,maxlevel,dim1,dim2,addim1,addim2;
  double min1,max1,step1,min2,max2,step2,tol;
  Init_eos &table_init;
  Eos_thin_entry **eos_entry;
  inline double max(double,double) const;
  inline int refine() const;
  inline double table(int,int,int) const;
  void tableinfo(ostream&) const;
  void tableinfo(int,ostream&) const;
public:
  typedef enum {dt_tau,dt_second} derivative_type_t;

  Eos_table(int,int,int,int,int,int,int,int,
            double,double,double,double,double,Init_eos&);
  ~Eos_table();

  int get_level() const { return level; }
  int get_table_value(int,double,double,double&) const;
  int get_table_derivative(derivative_type_t,int,double,double,double&) const;
  void tableinfo(int,const char*) const;
};

/*******************************************************************************
   template <const int size> class Eos_cons_adtab
*******************************************************************************/

template <const int size> class Eos_cons_adtab : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;
 
  int first_is_tau,maxlevel,range_warning;
  Eos_table *eos_table;

  inline int get_table_value(int,double,double,double&) const;
  inline int approximate_dtau(int,double,double,double&) const;
  inline int approximate_deps(int,double,double,double&) const;
public:
  Eos_cons_adtab(int,int,int,int,int,int,int,
                 double,double,double,double,double,Init_eos&);
  virtual ~Eos_cons_adtab();
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double dpdeps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_adtab
*******************************************************************************/

template <const int size> class Eos_prim_adtab : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;
 
  int first_is_tau,maxlevel,range_warning;
  Eos_table *eos_table;

  inline int get_table_value(int,double,double,double&) const;
  inline int approximate_dp(int,double,double,double&) const;
  inline int approximate_dtau(int,double,double,double&) const;
public:
  Eos_prim_adtab(int,int,int,int,int,int,int,
                 double,double,double,double,double,Init_eos&);
  virtual ~Eos_prim_adtab();
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double depsdp(const var_t &) const;
  virtual double depsdtau(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_file
*******************************************************************************/

template <const int size> class Eos_cons_file : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;

  int range_warning,first_is_tau,dim_first,dim_eps;
  double min_first,max_first,step_first,min_eps,max_eps,step_eps;
  double *eostable[3];

  inline double table(int,int,int) const;
  inline int get_table_value(int,double,double,double&) const;
public:
  Eos_cons_file(const char*);
  virtual ~Eos_cons_file() {for (int i=2;i>=0;i--) delete [] eostable[i];}
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_file
*******************************************************************************/

template <const int size> class Eos_prim_file : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  int range_warning,first_is_tau,dim_first,dim_p;
  double min_first,max_first,step_first,min_p,max_p,step_p;
  double *eostable[3];

  inline double table(int,int,int) const;
  inline int get_table_value(int,double,double,double&) const;
public:
  Eos_prim_file(const char*);
  virtual ~Eos_prim_file() {for (int i=2;i>=0;i--) delete [] eostable[i];}
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_func
*******************************************************************************/

template <const int size> class Eos_cons_func : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;

  Init_eos *init_eos;

public:
  inline Eos_cons_func(Init_eos *);
  virtual ~Eos_cons_func() {}
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double dpdeps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_func
*******************************************************************************/

template <const int size> class Eos_prim_func : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  Init_eos *init_eos;

public:
  inline Eos_prim_func(Init_eos *);
  virtual ~Eos_prim_func() {}
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double depsdp(const var_t &)   const;
  virtual double depsdtau(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_ideal
*******************************************************************************/

template <const int size> class Eos_cons_ideal : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;

  double gamma,R;
public:
  inline Eos_cons_ideal(double, double lR = 1.0);
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double dpdtau(const var_t &) const;
  virtual double dpdeps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_ideal
*******************************************************************************/

template <const int size> class Eos_prim_ideal : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  double gamma,R;
public:
  inline Eos_prim_ideal(double, double lR = 1.0);
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double depsdp(const var_t &)   const;
  virtual double depsdtau(const var_t &) const;

  /* extensions for cws */

  virtual double dcs2dp(const var_t &)   const;
  virtual double dcs2dtau(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_nconv
*******************************************************************************/

template <const int size> class Eos_cons_nconv : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;
 
  double alpha, gamma, mintaunc, maxtaunc, c0, c1, c2, c3, c4, c5;

  inline double dpdtau(const var_t &) const;
  inline double dpdeps(const var_t &) const;
public:
  inline Eos_cons_nconv(double palpha=1.0);
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_nconv
*******************************************************************************/

template <const int size> class Eos_prim_nconv : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  double alpha, gamma, mintaunc, maxtaunc, c0, c1, c2, c3, c4, c5;

  inline double dpdtau(const var_t &) const;
  inline double dpdeps(const var_t &) const;
public:
  inline Eos_prim_nconv(double palpha=1.0);
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_osborne
*******************************************************************************/

template <const int size> class Eos_cons_osborne : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;

  double tau0,a1,a2,b0,b1,b2,c0,c1,phi0;
public:
  inline Eos_cons_osborne();
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_osborne
*******************************************************************************/

template <const int size> class Eos_prim_osborne : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  double tau0,a1,a2,b0,b1,b2,c0,c1,phi0;
public:
  inline Eos_prim_osborne();
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Eos_tabular
*******************************************************************************/

class Eos_tabular
{
public:
  virtual ~Eos_tabular() {}
  virtual double p(double tau,double eps) = 0;
  virtual double eps(double tau,double p) = 0;
  virtual double dp_deps(double tau,double eps) = 0;
  virtual double dp_dtau(double tau,double eps) = 0;
  virtual double cs2(double tau,double eps) = 0;
};

/*******************************************************************************
   class Eostab_cons
*******************************************************************************/

class Eostab_cons : public eos_cons_t
{
  Eos_tabular &ptable;
public:
  virtual double cs2(const var_t &u) const 
  {return ptable.cs2(u.tau(),u.eps());}
  virtual double p(const var_t &u) const
  {return ptable.p(u.tau(),u.eps());}
  Eostab_cons(Eos_tabular &pptable) : ptable(pptable) {}
};

/*******************************************************************************
   class Eostab_prim
*******************************************************************************/

class Eostab_prim : public eos_prim_t
{
  Eos_tabular &ptable;
public:
  virtual double cs2(const var_t &v) const
  {return ptable.cs2(v.tau(),ptable.eps(v.tau(),v.p()));}
  virtual double eps(const var_t &v) const
  {return ptable.eps(v.tau(),v.p());}
  Eostab_prim(Eos_tabular &pptable) : ptable(pptable) {}
};

/*******************************************************************************
   class Eos_tabular_tau
*******************************************************************************/

class Eos_tabular_tau : public Eos_tabular 
{
  double **_ptable;
  double taumin,taumax,epsmin,epsmax,htau,heps;
  int dimtau,dimeps;
  inline double ptable(int t,int e)
  {
    assert(0<=t && t<dimtau);
    assert(0<=e && e<dimeps);
    return _ptable[t][e];
  };
  inline void findintable(double tau,double eps,
                          int &t,int &e,
                          double &psw,double &pse,double &pnw,double &pne,
                          double &lambdat,double &lambdae);
public:
  virtual double p(double tau,double eps);
  virtual double eps(double tau,double p);
  virtual double dp_deps(double tau,double eps);
  virtual double dp_dtau(double tau,double eps);
  virtual double cs2(double tau,double eps);
  Eos_tabular_tau(eos_cons_t *peos_c,eos_prim_t *peos_p,
                  double ptaumin,double ptaumax,
                  double pepsmin,double pepsmax,
                  int pdimtau,int pdimeps);
  ~Eos_tabular_tau()
  {
    for (int t=0;t<dimtau;t++) delete [] _ptable[t];
    delete [] _ptable;
  }
};

/*******************************************************************************
   class Eos_tabular_rho
*******************************************************************************/

class Eos_tabular_rho : public Eos_tabular
{
  double **_ptable;
  double rhomin,rhomax,epsmin,epsmax,hrho,heps;
  int dimrho,dimeps;
  inline double ptable(int r,int e)
  {
    assert(0<=r && r<dimrho);
    assert(0<=e && e<dimeps);
    return _ptable[r][e];
  };
  inline void findintable(double rho,double eps,
                          int &r,int &e,
                          double &psw,double &pse,double &pnw,double &pne,
                          double &lambdar,double &lambdae);
public:
  virtual double p(double tau,double eps);
  virtual double eps(double tau,double p);
  virtual double dp_deps(double tau,double eps);
  virtual double dp_dtau(double tau,double eps);
  virtual double cs2(double tau,double eps);
  Eos_tabular_rho(eos_cons_t *peos_c,eos_prim_t *peos_p,
                  double ptaumin,double ptaumax,
                  double pepsmin,double pepsmax,
                  int pdimtau,int pdimeps);
  ~Eos_tabular_rho()
  {
    for (int r=0;r<dimrho;r++) delete [] _ptable[r];
    delete [] _ptable;
  }
};



/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_tmv
*******************************************************************************/

template <const int size> class Eos_cons_tmv : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;

  double r,cvtr,thetavib,alpha,newtoneps;

  inline double f(double, double) const;
  inline double df(double) const;
public:
  inline Eos_cons_tmv();
  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_tmv
*******************************************************************************/

template <const int size> class Eos_prim_tmv : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;

  double r,cvtr,thetavib,alpha;

  inline double deps(double) const;
public:
  inline Eos_prim_tmv();
  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <const int size> class Eos_cons_waals
*******************************************************************************/

template <const int size> class Eos_cons_waals : public Eos_cons <size>
{
protected:
  typedef typename Eos_cons <size> :: var_t var_t;
public:
  const double a,b,R,cv,eps0;

  Eos_cons_waals() : Eos_cons<size>(1), a(1684.54), b(0.001692),
                     R(461.5), cv(1401.88), eps0(0.0) {}

  virtual double cs2(const var_t &) const;
  virtual double p(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double dpdtau(const var_t &) const;
  virtual double dpdeps(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};

/*******************************************************************************
   template <const int size> class Eos_prim_waals
*******************************************************************************/

template <const int size> class Eos_prim_waals : public Eos_prim <size>
{
protected:
  typedef typename Eos_prim <size> :: var_t var_t;
public:
  const double a,b,R,cv,eps0;

  Eos_prim_waals() : Eos_prim<size>(2), a(1684.54), b(0.001692),
                     R(461.5), cv(1401.88), eps0(0.0) {}

  virtual double cs2(const var_t &) const;
  virtual double eps(const var_t &) const;

  /* extensions for cws and conservative solvers */

  virtual double depsdp(const var_t &)   const;
  virtual double depsdtau(const var_t &) const;

  /* extensions for cws */

  virtual double dcs2dp(const var_t &)   const;
  virtual double dcs2dtau(const var_t &) const;

  /* extensions for 2d applications */

  virtual double T(const var_t &) const;
};


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Sun_eos
*******************************************************************************/

class Sun_eos
{
protected:
  inline double SIGN(double,double) const;

public:
  const long double boltz,amass,chi,cii,hmu,mu1;
  const double rho0,T0,Rgask,mu0,p0,eps0,cs20;
  const double TOL;
  const int ITMAX;

  inline Sun_eos();

  float zbrent(int,float,float,float,float,float) const;
  void eos_a(float,float,double*,double*,double*,double*) const;
  void eos_b(float,float,double*,double*,double*) const;
  float ener_a(float,float,float) const;
  float ener_b(float,float,float) const;
};

/*******************************************************************************
   class Sun_init_conseos
*******************************************************************************/

class Sun_init_conseos : public Sun_eos, public Init_eos
{
public:
  Sun_init_conseos() : Init_eos(1) {};
  virtual void operator()(double,double,
                          double&,double&,double&,double&,double&) const;
};

/*******************************************************************************
   class Sun_init_primeos
*******************************************************************************/

class Sun_init_primeos : public Sun_eos, public Init_eos
{
public:
  Sun_init_primeos() {};
  virtual void operator()(double,double,
                          double&,double&,double&,double&,double&) const;
};

/********************************************************/
/********************************************************/


template <const int size> class Wrap_eos_init_conseos;
template <const int size> class Wrap_eos_init_primeos;
typedef Wrap_eos_init_conseos<8> wrap_eos_init_conseos_t;
typedef Wrap_eos_init_primeos<8> wrap_eos_init_primeos_t;

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   class Wrap_eos_init_conseos
*******************************************************************************/

template <const int size>
class Wrap_eos_init_conseos : public Init_eos
{
protected:
  typedef Eos_cons<size> eos_t;
  typedef Eos_prim<size> conv_eos_t;
  eos_t *eos;
  conv_eos_t *conv_eos;

public:
  inline Wrap_eos_init_conseos(eos_t*, conv_eos_t*);
  virtual void operator()(double,double,
                          double&,double&,double&,double&,double&) const;
};

/*******************************************************************************
   class Wrap_eos_init_primeos
*******************************************************************************/

template <const int size>
class Wrap_eos_init_primeos : public Init_eos
{
protected:
  typedef Eos_prim<size> eos_t;
  typedef Eos_cons<size> conv_eos_t;
  eos_t *eos;
  conv_eos_t *conv_eos;

public:
  inline Wrap_eos_init_primeos(eos_t*, conv_eos_t*);
  virtual void operator()(double,double,
                          double&,double&,double&,double&,double&) const;
};


/********************************************************/
/********************************************************/

typedef Eos_cons_tmv<8> eos_cons_tmv_t;
typedef Eos_prim_tmv<8> eos_prim_tmv_t;
typedef Eos_cons_file<8> eos_cons_file_t;
typedef Eos_prim_file<8> eos_prim_file_t;
typedef Eos_cons_func<8> eos_cons_func_t;
typedef Eos_prim_func<8> eos_prim_func_t;
typedef Eos_cons_adtab<8> eos_cons_adtab_t;
typedef Eos_prim_adtab<8> eos_prim_adtab_t;
typedef Eos_cons_ideal<8> eos_cons_ideal_t;
typedef Eos_prim_ideal<8> eos_prim_ideal_t;
typedef Eos_cons_waals<8> eos_cons_waals_t;
typedef Eos_prim_waals<8> eos_prim_waals_t;
typedef Eos_cons_nconv<8> eos_cons_nconv_t;
typedef Eos_prim_nconv<8> eos_prim_nconv_t;
typedef Eos_cons_osborne<8> eos_cons_osborne_t;
typedef Eos_prim_osborne<8> eos_prim_osborne_t;


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux
*******************************************************************************/

template <class EQNS> class Flux
{
protected:
  typedef typename EQNS::cons_t cons_t;
  typedef typename EQNS::prim_t prim_t;
  typedef typename EQNS::vec_t  vec_t;
  typedef typename EQNS::mat_t  mat_t;

  const double eps;
  const double del;

  EQNS *eqns;

  Flux() : eps(1e-12), del(5e-12) {eqns = new EQNS();}
  double min(double pl, double pr) const {return (pl <= pr) ? pl : pr;}
  double max(double pl, double pr) const {return (pl >= pr) ? pl : pr;}
  inline double phi(double) const;
  inline double psi(double) const;
  double ce(const prim_t&, const prim_t&) const;
  double ch(const prim_t&, const prim_t&, const prim_t&) const;
public:
  virtual ~Flux() {delete eqns;}
  inline const EQNS &get_eqns() const {return *eqns;}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const = 0;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux
*******************************************************************************/


template <class EQNS>
void Flux<EQNS>::operator()(const cons_t &pconsl,
                            const cons_t &pconsr,
                            cons_t &pret,
                            double pdt,double phl,double phr) const
{
  prim_t priml(pconsl);
  prim_t primr(pconsr);

  (*this)(priml,primr,pret,pdt,phl,phr);
}

template <class EQNS>
inline double Flux<EQNS>::phi(double px) const
{
  double lret;

  if (px < eps)
    lret = 0.0;
  else if (px < del)
    lret = (px - eps) / (del - eps);
  else
    lret = 1.0;

  return lret;
}

template <class EQNS>
inline double Flux<EQNS>::psi(double px) const
{
  double lret;

  if (px < -del)
    lret = 1.0;
  else if (px < -eps)
    lret = 0.5 * (px + del) / (eps - del) + 1.0; 
  else if (px < eps)
    lret = 0.5;
  else if (px < del)
    lret = 0.5 * (px - del) / (eps - del);
  else
    lret = 0.0;

  return lret;
}

template <class EQNS>
double Flux<EQNS>::ce(const prim_t &ppriml, const prim_t &pprimr) const
{
  double lret = 0.0;
  vec_t lambdal,lambdar;

  eqns->lambda_vec(ppriml,lambdal);
  eqns->lambda_vec(pprimr,lambdar);

  for (int j=0;j<lambdal.size();j++)
    if ((eqns->nonlinear(j)) && (lambdal[j] < 0.0) && (0.0 < lambdar[j]))
      lret += lambdar[j] - lambdal[j]; 

  return lret;
}

template <class EQNS>
double Flux<EQNS>::ch(const prim_t &ppriml,
                      const prim_t &pprimr,
                      const prim_t &pprimmv) const
{
  int m = pprimmv.size()-1;
  int e = eqns->ewave();

  double lvmax_mv = eqns->lambda(m,pprimmv) - eqns->lambda(e,pprimmv);

  return max(0.0,(eqns->lambda(e,pprimr) - eqns->lambda(e,ppriml)) - lvmax_mv);
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hlle
*******************************************************************************/

template <class EQNS> class Flux_hlle : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;
public:
  Flux_hlle() : Flux<EQNS>() {}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hlle
*******************************************************************************/

template <class EQNS>
void Flux_hlle<EQNS>::operator()(const prim_t &ppriml,
                                 const prim_t &pprimr,
                                 cons_t &pret,
                                 double pdt,double phl,double phr) const
{
  cons_t consl(ppriml);
  cons_t consr(pprimr);

  (*this)(consl,consr,pret,pdt,phl,phr);
}

template <class EQNS>
void Flux_hlle<EQNS>::operator()(const cons_t &pconsl,
                                 const cons_t &pconsr,
                                 cons_t &pret,
                                 double pdt,double phl,double phr) const
{
  /* check if equations of state are valid */

  assert(pconsl.get_eos() == pconsr.get_eos());
  assert(pconsl.get_conv_eos() == pconsr.get_conv_eos());
  assert(pret.get_eos() == pconsr.get_eos());
  assert(pret.get_conv_eos() == pconsr.get_conv_eos());

  /* calculate jump */

  cons_t diff(pconsr - pconsl);

  /* calculate mean value */

  prim_t mval((pconsl + pconsr) * 0.5);

  /* calculate blm,brp */

  prim_t priml(pconsl), primr(pconsr);
  int m = pconsl.size()-1;

  double blm = min(0.0,min(this->eqns->lambda(0,priml),
                           this->eqns->lambda(0,mval)));
  double brp = max(0.0,max(this->eqns->lambda(m,primr),
                           this->eqns->lambda(m,mval)));
  
  /* calculate flux */

  cons_t fl(pconsl),fr(pconsr);

  this->eqns->f(pconsl,fl);
  this->eqns->f(pconsr,fr);

  assert(brp - blm > 1e-12);

  pret = ((fl * brp - fr * blm) + diff * brp * blm) / (brp - blm);
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_dwc
*******************************************************************************/

template <class EQNS> class Flux_dwc : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;

  const double eps;
  const double visc;
  Flux_hlle<EQNS> *flux_hlle;
public:
  Flux_dwc() : Flux<EQNS>(), eps(1e-12), visc(0.05)
  {
    flux_hlle = new Flux_hlle<EQNS>;
  }
  virtual ~Flux_dwc() {delete flux_hlle;}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double,double,double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double,double,double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_dwc
*******************************************************************************/

template <class EQNS>
void Flux_dwc<EQNS>::operator()(const prim_t &ppriml,
                                const prim_t &pprimr,
                                cons_t &pret,
                                double pdt,double phl,double phr) const
{
  abort();
}

template <class EQNS>
void Flux_dwc<EQNS>::operator()(const cons_t &pconsl,
                                const cons_t &pconsr,
                                cons_t &pret,
                                double pdt,double phl,double phr) const
{
  int k;

  /* check if equations of state are valid */

  assert(pconsl.get_eos() == pconsr.get_eos());
  assert(pconsl.get_conv_eos() == pconsr.get_conv_eos());
  assert(pret.get_eos() == pconsr.get_eos());
  assert(pret.get_conv_eos() == pconsr.get_conv_eos());

  /* calculate mean value */

  prim_t priml(pconsl),primr(pconsr);
  cons_t mval((pconsl + pconsr) * 0.5);
  prim_t mpval(mval);

  /* perform hll-fix */

  double lphi = phi(ch(priml,primr,mpval));

  if (lphi > eps)
    (*flux_hlle)(pconsl,pconsr,pret,pdt,phl,phr);
  else
    pret.clear();

  if (lphi < 1.0 - eps)
    {
      pret *= lphi;

      /*** calculate native DW-flux ***/

      /* calculate jump */

      cons_t diff(pconsl - pconsr);

      /* calculate upwind states */

      vec_t lambda;
      vec_t *v0 = new vec_t[lambda.size()];

      this->eqns->lambda_vec(mpval,lambda);

      for (k=0;k<lambda.size();k++)
        {
          double lpsi = psi(lambda[k]);
          v0[k] = (pconsl * (1.0 - lpsi) + pconsr * lpsi).vec();
        }

      /* calculate right hand side of linear system */

      mat_t lmat;
      vec_t rhs;

      this->eqns->l_mat(mval,lmat);

      for (k=0;k<lambda.size();k++)
        rhs[k] = lmat.row(k) * v0[k];

      /* solve linear system */

      mat_t rmat;

      this->eqns->r_mat(mval,rmat);

      cons_t ref(priml.get_conv_eos(),priml.get_eos(),rmat*rhs);

      /* calculate native flux */

      cons_t flux(pconsl);

      this->eqns->f(ref,flux);

      cons_t natflux(flux + diff * min(visc*sqrt(diff*diff),1.0));

      /*** calculate DW-flux ***/

      pret += natflux * (1.0 - lphi);

      /* cleaning up */

      delete [] v0;
    }
}

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_dw
*******************************************************************************/

template <class EQNS> class Flux_dw : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;

  const double eps;
  const double visc;
  Flux_hlle<EQNS> *flux_hlle;
public:
  Flux_dw() : Flux<EQNS>(), eps(1e-12), visc(0.05)
  {
    flux_hlle = new Flux_hlle<EQNS>;
  }
  virtual ~Flux_dw() {delete flux_hlle;}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double,double,double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double,double,double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_dw
*******************************************************************************/

template <class EQNS>
void Flux_dw<EQNS>::operator()(const prim_t &ppriml,
                               const prim_t &pprimr,
                               cons_t &pret,
                               double pdt,double phl,double phr) const
{
  abort();
}

template <class EQNS>
void Flux_dw<EQNS>::operator()(const cons_t &pconsl,
                               const cons_t &pconsr,
                               cons_t &pret,
                               double pdt,double phl,double phr) const
{
  int k;

  /* check if equations of state are valid */

  assert(pconsl.get_eos() == pconsr.get_eos());
  assert(pconsl.get_conv_eos() == pconsr.get_conv_eos());
  assert(pret.get_eos() == pconsr.get_eos());
  assert(pret.get_conv_eos() == pconsr.get_conv_eos());

  /* calculate mean value */

  prim_t priml(pconsl),primr(pconsr);
  prim_t mval((priml + primr) * 0.5);

  /* perform hll-fix */

  double lphi = phi(ch(priml,primr,mval));

  if (lphi > eps)
    (*flux_hlle)(pconsl,pconsr,pret,pdt,phl,phr);
  else
    pret.clear();

  if (lphi < 1.0 - eps)
    {
      pret *= lphi;

      /*** calculate native DW-flux ***/

      /* calculate jump */

      cons_t diff(pconsl - pconsr);

      /* calculate upwind states */

      vec_t lambda;
#ifdef _MHD2D
      vec_t v0[8];
#else
      vec_t *v0 = new vec_t[lambda.size()];
#endif

      this->eqns->lambda_vec(mval,lambda);

      for (k=0;k<lambda.size();k++)
        {
          double lpsi = psi(lambda[k]);
          v0[k] = priml.vec() * (1.0 - lpsi) + primr.vec() * lpsi;
        }

      /* calculate right hand side of linear system */

      mat_t lmat;
      vec_t rhs;

      this->eqns->l_mat(mval,lmat);

      for (k=0;k<lambda.size();k++)
        rhs[k] = lmat.row(k) * v0[k];

      /* solve linear system */

      mat_t rmat;

      this->eqns->r_mat(mval,rmat);

      prim_t ref(priml.get_eos(),priml.get_conv_eos(),rmat*rhs);

      /* calculate native flux */

      cons_t flux(pconsl);

      this->eqns->f(ref,flux);

      cons_t natflux(flux + diff * min(visc*sqrt(diff*diff),1.0));

      /*** calculate DW-flux ***/

      pret += natflux * (1.0 - lphi);

      /* cleaning up */

#ifndef _MHD2D
      delete [] v0;
#endif
    }
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hllemc
*******************************************************************************/

template <class EQNS> class Flux_hllemc : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;
public:
  Flux_hllemc() : Flux<EQNS>() {}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hllemc
*******************************************************************************/

template <class EQNS>
void Flux_hllemc<EQNS>::operator()(const prim_t &ppriml,
                                   const prim_t &pprimr,
                                   cons_t &pret,
                                   double pdt,double phl,double phr) const
{
  cons_t consl(ppriml);
  cons_t consr(pprimr);

  (*this)(consl,consr,pret,pdt,phl,phr);
}

template <class EQNS>
void Flux_hllemc<EQNS>::operator()(const cons_t &pconsl,
                                   const cons_t &pconsr,
                                   cons_t &pret,
                                   double pdt,double phl,double phr) const
{
  int i;

  /* check if equations of state are valid */

  assert(pconsl.get_eos() == pconsr.get_eos());
  assert(pconsl.get_conv_eos() == pconsr.get_conv_eos());
  assert(pret.get_eos() == pconsr.get_eos());
  assert(pret.get_conv_eos() == pconsr.get_conv_eos());

  /* calculate jump */

  cons_t diff(pconsr - pconsl);

  /* calculate mean value */

  cons_t cmval((pconsl + pconsr) * 0.5);
  prim_t pmval(cmval);

  /* calculate eigenvectors for mean value */

  mat_t rmat,lmat;
  vec_t lambdamv;

  this->eqns->r_mat(cmval,rmat);
  this->eqns->l_mat(cmval,lmat);
  this->eqns->lambda_vec(pmval,lambdamv);

  /* decompose jump into wave contributions */

  vec_t alpha = lmat * diff.vec();

  /* calculate speeds */

  vec_t vel;
  int ewave = this->eqns->ewave();

  for (i=0;i<ewave;i++)
    vel[i] = lambdamv[ewave] - lambdamv[i];
  vel[ewave] = fabs(lambdamv[ewave]);

  /* calculate anti-diffusion coefficients */

  vec_t delta;
  double phi1 = phi(vel[0] - vel[1]);
  double phi2 = 0.0;

  phi2 = phi(cmval.B2());

  double div = vel[0] + (1.0 - phi2) * vel[ewave];
  
  for (i=1;i<ewave;i++)
    delta[i] = phi1 * vel[0] / (div + vel[i]);
  delta[ewave] = phi1 * vel[0] / (vel[0] + vel[ewave]); 

  /* calculate anti-diffusion term */

  vec_t adt;
  int m = pconsl.size()-1;

  adt = rmat.col(ewave) * alpha[ewave] * delta[ewave];
  for (i=1;i<ewave;i++)
    adt += (rmat.col(i) * alpha[i] + rmat.col(m-i) * alpha[m-i]) * delta[i];

  cons_t cadt(diff.get_eos(),diff.get_conv_eos(),adt);

  /* calculate blm,brp */

  prim_t priml(pconsl), primr(pconsr);

  double blm = min(0.0,min(this->eqns->lambda(0,priml),lambdamv[0]));
  double brp = max(0.0,max(this->eqns->lambda(m,primr),lambdamv[m]));
  
  /* calculate flux */

  cons_t fl(pconsl),fr(pconsr);

  this->eqns->f(pconsl,fl);
  this->eqns->f(pconsr,fr);

  assert(brp - blm > 1e-12);

  pret = ((fl * brp - fr * blm) + (diff - cadt) * brp * blm) / (brp - blm);
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hllem
*******************************************************************************/

template <class EQNS> class Flux_hllem : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;
public:
  Flux_hllem() : Flux<EQNS>() {}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_hllem
*******************************************************************************/

template <class EQNS>
void Flux_hllem<EQNS>::operator()(const prim_t &ppriml,
                                  const prim_t &pprimr,
                                  cons_t &pret,
                                  double pdt,double phl,double phr) const
{
  int i;

  /* check if equations of state are valid */

  assert(ppriml.get_eos() == pprimr.get_eos());
  assert(ppriml.get_conv_eos() == pprimr.get_conv_eos());
  assert(pret.get_eos() == pprimr.get_conv_eos());
  assert(pret.get_conv_eos() == pprimr.get_eos());

  /* calculate jump */

  prim_t diff(pprimr - ppriml);

  /* calculate mean value */

  prim_t mval((ppriml + pprimr) * 0.5);

  /* calculate eigenvectors for mean value */

  mat_t rmat,lmat;
  vec_t lambdamv;

  this->eqns->r_mat(mval,rmat);
  this->eqns->l_mat(mval,lmat);
  this->eqns->lambda_vec(mval,lambdamv);

  /* decompose jump into wave contributions */

  vec_t alpha = lmat * diff.vec();

  /* calculate speeds */

  vec_t vel;
  int ewave = this->eqns->ewave();

  for (i=0;i<ewave;i++)
    vel[i] = lambdamv[ewave] - lambdamv[i];
  vel[ewave] = fabs(lambdamv[ewave]);

  /* calculate anti-diffusion coefficients */

  vec_t delta;
  double phi1 = phi(vel[0] - vel[1]);
  double phi2 = 0.0;

  phi2 = phi(mval.B2());

  double div = vel[0] + (1.0 - phi2) * vel[ewave];
  
  for (i=1;i<ewave;i++)
    delta[i] = phi1 * vel[0] / (div + vel[i]);
  delta[ewave] = phi1 * vel[0] / (vel[0] + vel[ewave]); 

  /* calculate anti-diffusion term */

  vec_t adt;
  int m = ppriml.size()-1;

  adt = rmat.col(ewave) * alpha[ewave] * delta[ewave];
  for (i=1;i<ewave;i++)
    adt += (rmat.col(i) * alpha[i] + rmat.col(m-i) * alpha[m-i]) * delta[i];

  /* calculate blm,brp */

  double blm = min(0.0,min(this->eqns->lambda(0,ppriml),lambdamv[0]));
  double brp = max(0.0,max(this->eqns->lambda(m,pprimr),lambdamv[m]));
  
  /* calculate flux */

  cons_t fl(ppriml.get_conv_eos(),ppriml.get_eos());
  cons_t fr(pprimr.get_conv_eos(),pprimr.get_eos());

  this->eqns->f(ppriml,fl);
  this->eqns->f(pprimr,fr);

  mat_t dUdW;
  mval.dUdW(dUdW);
  cons_t corr(ppriml.get_conv_eos(),ppriml.get_eos(),dUdW*(diff.vec()-adt));

  assert(brp - blm > 1e-12);

  pret = ((fl * brp - fr * blm) + corr * brp * blm) / (brp - blm);
}

template <class EQNS>
void Flux_hllem<EQNS>::operator()(const cons_t &pconsl,
                                  const cons_t &pconsr,
                                  cons_t &pret,
                                  double pdt,double phl,double phr) const
{
  prim_t priml(pconsl);
  prim_t primr(pconsr);

  (*this)(priml,primr,pret,pdt,phl,phr);
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_lf
*******************************************************************************/

template <class EQNS> class Flux_lf : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;
public:
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_lf
*******************************************************************************/

template <class EQNS>
void Flux_lf<EQNS>::operator()(const prim_t &ppriml,
                               const prim_t &pprimr,
                               cons_t &pret,
                               double pdt,double phl,double phr) const
{
  cons_t consl(ppriml);
  cons_t consr(pprimr);

  (*this)(consl,consr,pret,pdt,phl,phr);
}

template <class EQNS>
void Flux_lf<EQNS>::operator()(const cons_t &pconsl,
                               const cons_t &pconsr,
                               cons_t &pret,
                               double pdt,double phl,double phr) const
{
  /* check if equations of state are valid */

  assert(pconsl.get_eos() == pconsr.get_eos());
  assert(pconsl.get_conv_eos() == pconsr.get_conv_eos());
  assert(pret.get_eos() == pconsr.get_eos());
  assert(pret.get_conv_eos() == pconsr.get_conv_eos());

  /* calculate flux */
  double lambda=pdt/(0.5*(phl+phr));
  cons_t lfl(pconsl),lfr(pconsr);
  this->eqns->f(pconsl,lfl);
  this->eqns->f(pconsr,lfr);
  pret=(lfl+lfr)*0.5 + (pconsl-pconsr)*0.5/lambda;
}


/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                           CLASS DEFINITIONS                            ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_vfroe
*******************************************************************************/

template <class EQNS> class Flux_vfroe : public Flux <EQNS>
{
protected:
  typedef typename Flux <EQNS> :: cons_t cons_t;
  typedef typename Flux <EQNS> :: prim_t prim_t;
  typedef typename Flux <EQNS> :: vec_t vec_t;
  typedef typename Flux <EQNS> :: mat_t mat_t;

  const double eps;
  Flux_hlle<EQNS> *flux_hlle;
public:
  Flux_vfroe() : Flux<EQNS>(), eps(1e-12) {flux_hlle = new Flux_hlle<EQNS>;}
  virtual ~Flux_vfroe() {delete flux_hlle;}
  virtual void operator()(const prim_t&, const prim_t&, cons_t&,
                          double, double, double) const;
  virtual void operator()(const cons_t&, const cons_t&, cons_t&,
                          double, double, double) const;
};

/******************************************************************************
 ******************************************************************************
 ***                                                                        ***
 ***                         CLASS IMPLEMENTATIONS                          ***
 ***                                                                        ***
 ******************************************************************************
 ******************************************************************************/

/*******************************************************************************
   template <class EQNS> class Flux_vfroe
*******************************************************************************/

template <class EQNS>
void Flux_vfroe<EQNS>::operator()(const prim_t &ppriml,
                                  const prim_t &pprimr,
                                  cons_t &pret,
                                  double pdt,double phl,double phr) const
{
  int i,k,done;
  vec_t intvec,alpha,lambda;
  mat_t lmat,rmat;

  /* check if equations of state are valid */

  assert(ppriml.get_eos() == pprimr.get_eos());
  assert(ppriml.get_conv_eos() == pprimr.get_conv_eos());
  assert(pret.get_eos() == pprimr.get_conv_eos());
  assert(pret.get_conv_eos() == pprimr.get_eos());

  /* calculate mean value */

  prim_t mval((ppriml + pprimr) * 0.5);

  /* perform hll- and entropy-fix */

  double lphi = phi(ce(ppriml,pprimr)+ch(ppriml,pprimr,mval));

  if (lphi > eps)
    (*flux_hlle)(ppriml,pprimr,pret,pdt,phl,phr);
  else
    pret.clear();

  if (lphi < 1.0 - eps)
    {
      pret *= lphi;

      /*** calculate native VFROE-flux ***/

      /* calculate jump */

      prim_t diff(pprimr - ppriml);

      /* calculate eigenvalues, eigenvectors and wave strenghts */

      this->eqns->lambda_vec(mval,lambda);
      this->eqns->l_mat(mval,lmat);
      this->eqns->r_mat(mval,rmat);
      alpha = lmat * diff.vec();

      /* determine state on cell interface */

      assert(lambda.size() >= 3);

      k = -1;
      done = 0;
      while((!done) && (k < lambda.size()-1))
        if (lambda[k+1] <= eps)
          k++;
        else
          done = 1;

      if (k > lambda.size()/2 - 1)
        {
          intvec = pprimr.vec();
          for (i=lambda.size()-1;i>k;i--)
            intvec -= rmat.col(i) * alpha[i];
        }
      else
        {
          intvec = ppriml.vec();
          for (i=0;i<=k;i++)
            intvec += rmat.col(i) * alpha[i]; 
        }

      prim_t intstate(ppriml.get_eos(),ppriml.get_conv_eos(),intvec);
  
      /* calculate native flux */

      cons_t natflux(ppriml.get_conv_eos(),ppriml.get_eos());

      this->eqns->f(intstate,natflux);

      /*** calculate VFROE-flux ***/

      pret += natflux * (1.0 - lphi);
    }
}

template <class EQNS>
void Flux_vfroe<EQNS>::operator()(const cons_t &pconsl,
                                  const cons_t &pconsr,
                                  cons_t &pret,
                                  double pdt,double phl,double phr) const
{
  Flux<EQNS>::operator()(pconsl,pconsr,pret,pdt,phl,phr);
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

// ** NEWLY INSERTED FROM MHD_EQNS_HH **


} // end namespace Mhd 
#endif

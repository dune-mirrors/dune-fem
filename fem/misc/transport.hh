#ifndef __TRANSPORT_HH__
#define __TRANSPORT_HH__

// Dune header includes
#include <dune/common/function.hh>
#include <dune/common/fmatrix.hh>

using namespace Dune;

#include <math.h>

//**************************************************************************
//! Describes the Problem physically  
//!
//!  We want to solve the problem 
//
//!  d_t B(c) + div ( f(u c) - grad ( D (u) c ) = q(c) 
//
//!  where c is the unknown and u a given velocity, B ( Retardation )
//!  a given function which can be nonlinear and D a given dispersion tensor. 
//!  q is a right hand side reaction term and f(u c) is the so called flux
//!  function.
//
//**************************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
class LinearTransport 
{
public:
  // usefull typedefs 
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType FieldType;

  enum { dimRange = RangeType::dimension };
  
protected:  
  typedef FieldVector < FieldType , dim > FieldVectorType;

  // flux function, describing the problem

  // the corresponding numerical flux function
  struct FluxFunction  
  {
    //! evaluate Function 
    FieldType evaluate ( FieldVectorType & velo , FieldVectorType &normal,
            FieldType u , const FieldType & df_u ) const;  
    
    //! evaluate first derivative 
    void jacobian ( FieldVectorType & velo , FieldVectorType & normal,
            RangeType & u , RangeType & ret ) const;  
  };

private:
  // the initial data and exact solution 
  class ExactSolution : public Function < FunctionSpaceType , ExactSolution > 
  {
    public:
      // Constructor 
      ExactSolution (FunctionSpaceType &f);
     
      // evaluate the function
      void evaluate (const DomainType & x, DomainFieldType time, RangeType & ret) const ;
     
      // evaluate the function a time zero
      void evaluate (const DomainType & x, RangeType & ret) const; 
  };

public:
  //! Constructor 
  LinearTransport  (FunctionSpaceType &f);

  //! Destructor 
  ~LinearTransport ();
  
  //! initial data of the problem  
  const ExactSolution& initialData () const; 

  //! return true if exact Solution for problem exists 
  bool hasExactSolution() const;

  //! return exact Solution 
  const ExactSolution& exactSolution () const; 
    
  //! velocity of the transport problem
  void v ( const FieldVectorType &x, FieldVectorType & velo) const;

  //! flux function , returns local delta t ,
  double flux ( FieldVectorType  & velo , FieldVectorType & normal, 
          RangeType & uEl, RangeType & uNeigh, RangeType & ret ) const;
  
  //! return true if dispersion tensor is constant
  bool constDispersion () const;

  //! dispersion  coefficient
  template <class LocalFunctionType> 
  void D ( LocalFunctionType &lf, const FieldVectorType &x, FieldMatrix<FieldType,dim,dim> &disp) const;

  //! dispersion  coefficient
  template <class LocalFunctionType> 
  void q ( LocalFunctionType &lf, const FieldVectorType &x, 
           RangeType &ret) const;
  
  //! do the boundary 
  void boundary ( const FieldVectorType &x, const RangeType & u, RangeType & ret ) const;

  //! if B(u) is not the identity of u then return true
  //! and we have to evaluate the first derivative of B
  bool hasRetardation () const;

  //! returns one over jacobian of Reatardation 
  template <class LocalFunctionIteratorType> 
  void RetardationJacobian ( LocalFunctionIteratorType &lf, RangeType & ret ) const;

  void setTime ( const FieldType t) const;
  const FieldType getTime () const;
  
protected:
  mutable FieldType time_;
  
  //! the flux function corresponding to this problem
  const FluxFunction fluxFunc_;
  
  // diffusion coefficient 
  mutable double D_;

  // temporay data for storage of the flux derivative 
  RangeType * upwind_;

private:  
  // exact Solution 
  const ExactSolution exSol_;
}; // end class LinearTransport

// Implementation
template <int dim, int dimworld, class FunctionSpaceType> 
inline LinearTransport<dim,dimworld,FunctionSpaceType>::LinearTransport
(FunctionSpaceType &f ) : 
  time_ (0.0) , fluxFunc_ () , D_ (0.1) , upwind_ ( 0 ) , exSol_ ( f ) 
{
  upwind_ = new RangeType ();
}

//! Destructor 
template <int dim, int dimworld, class FunctionSpaceType> 
inline LinearTransport<dim,dimworld,FunctionSpaceType>::
~LinearTransport () 
{
  if( upwind_ ) delete upwind_;
}
  
//! initial data of the problem  
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename LinearTransport<dim,dimworld,FunctionSpaceType>::ExactSolution& 
LinearTransport<dim,dimworld,FunctionSpaceType>::initialData () const 
{ 
  return exSol_;   
}

  //! return true if exact Solution for problem exists 
template <int dim, int dimworld, class FunctionSpaceType> 
inline bool LinearTransport<dim,dimworld,FunctionSpaceType>::hasExactSolution() const 
{ 
  return true; 
}

//! return exact Solution 
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename LinearTransport<dim,dimworld,FunctionSpaceType>::ExactSolution& 
LinearTransport<dim,dimworld,FunctionSpaceType>::exactSolution () const 
{
  return exSol_; 
}
    
//! velocity of the transport problem
template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>:: 
v ( const FieldVectorType &x, FieldVectorType & velo) const
{
  velo = 0.0;
  velo[0] = 1.0;
}

//! flux function , returns local delta t ,
template <int dim, int dimworld, class FunctionSpaceType> 
inline double LinearTransport<dim,dimworld,FunctionSpaceType>::flux 
( FieldVectorType & velo , FieldVectorType & normal, 
  RangeType & uEl, RangeType & uNeigh, RangeType & retFlux ) const 
{
  enum { dimrange = RangeType::dimension };

  retFlux = 0.0;

  RangeType &upwind = *upwind_; 

  fluxFunc_.jacobian(velo,normal,uEl, upwind );
  double dt = upwind_->infinity_norm();

  for(int l=0; l<dimrange; l++)
  {
    if(upwind[l] > 0.0) 
    {
      retFlux[l] += fluxFunc_.evaluate(velo,normal,uEl[l],upwind[l]);
      dt = upwind.infinity_norm();
    }
    // * else statement? only upwind or what?
  } 

  fluxFunc_.jacobian(velo,normal,uNeigh, upwind );
  for(int l=0; l<dimrange; l++)
  {
    if(upwind[l] < 0.0)
    {
      retFlux[l] += fluxFunc_.evaluate(velo,normal,uNeigh[l],upwind[l]);
      dt = std::min(upwind.infinity_norm(),dt);
    }
  }
  return dt;
}

//! return true if dispersion  coefficient is constant
template <int dim, int dimworld, class FunctionSpaceType> 
inline bool LinearTransport<dim,dimworld,FunctionSpaceType>::
constDispersion () const 
{
  return true;
}


//! dispersion  coefficient
template <int dim, int dimworld, class FunctionSpaceType> 
template <class LocalFunctionType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::
D ( LocalFunctionType &lf, const FieldVectorType &x, FieldMatrix<FieldType,dim,dim> &disp) const 
{
  // rotating pulse 
  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++)
      if(j==i) disp(j,i) = D_;
      else disp(j,i) = 0.0;
}

//! dispersion  coefficient
template <int dim, int dimworld, class FunctionSpaceType> 
template <class LocalFunctionType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::
q ( LocalFunctionType &lf, const FieldVectorType &x, 
    RangeType & ret ) const 
{
  ret = lf[0];
  return;
}

//! do the boundary 
template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::
boundary ( const FieldVectorType &x, const RangeType & u , RangeType & ret ) const 
{
  // return u, for no flow boundary 
  ret = u;
  return;
}

//! if B(u) is not the identity of u then return true
//! and we have to evaluate the first derivative of B
template <int dim, int dimworld, class FunctionSpaceType> 
inline bool LinearTransport<dim,dimworld,FunctionSpaceType>::
hasRetardation () const 
{
  return false;
}

//! returns one over jacobian of B prime 
template <int dim, int dimworld, class FunctionSpaceType> 
template <class LocalFunctionIteratorType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>:: 
RetardationJacobian ( LocalFunctionIteratorType &lf , RangeType & ret  ) const 
{
  ret = 1.0;
}

template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>:: 
setTime ( const FieldType time ) const 
{
  time_ = time;
}

template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename LinearTransport<dim,dimworld,FunctionSpaceType>::FieldType  
LinearTransport<dim,dimworld,FunctionSpaceType>::getTime () const 
{
  return time_;
}

//*******************************************************************
//
// class ExactSolution 
//
//*******************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
inline LinearTransport<dim,dimworld,FunctionSpaceType>::ExactSolution::
ExactSolution ( FunctionSpaceType &f) : 
 Function < FunctionSpaceType , ExactSolution > ( f )  
{
}  

template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::ExactSolution::
evaluate (const DomainType & x , DomainFieldType time, RangeType & ret) const 
{
  ret = 0.0;
  return;
}

template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::ExactSolution::
evaluate (const DomainType & x1 , RangeType & ret) const 
{
  evaluate (x1,0.0,ret);
}

// evaluate flux function 
template <int dim, int dimworld, class FunctionSpaceType> 
inline typename LinearTransport<dim,dimworld,FunctionSpaceType>::FieldType 
LinearTransport<dim,dimworld,FunctionSpaceType>::FluxFunction::
evaluate ( FieldVectorType & velo , FieldVectorType &normal,
           FieldType u , const FieldType & df_u ) const  
{
  // we have a linear function
  return u * df_u;
}

//! evaluate first derivative 
template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearTransport<dim,dimworld,FunctionSpaceType>::FluxFunction::
jacobian ( FieldVectorType & velo , FieldVectorType & normal,
           RangeType & u , RangeType & ret ) const  
{
  // each component of ret is set to the value of the scalar product 
  ret = velo * normal;
  return;
}


#endif


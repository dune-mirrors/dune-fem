#ifndef __RIEMANNPROBLEM_CC__
#define __RIEMANNPROBLEM_CC__

using namespace Dune;

#include "../../misc/transport.hh"

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
class RiemannProblem : public LinearTransport < dim ,dimworld, FunctionSpaceType >
{
private:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType FieldType;
  
  typedef typename LinearTransport
    <dim,dimworld,FunctionSpaceType>::FieldVectorType FieldVectorType;
  // the initial data and exact solution 
  class Schock : public Function < FunctionSpaceType , Schock > 
  {
    //! Schock velocity 
    DomainType schock_;
    public:
      // Constructor 
      Schock (FunctionSpaceType &f);
     
      // evaluate the function
      void evaluate (const DomainType & x, DomainFieldType time, RangeType & ret) const ;
     
      // evaluate the function a time zero
      void evaluate (const DomainType & x, RangeType & ret) const; 
  };

public:
  //! Constructor 
  RiemannProblem  (FunctionSpaceType &f);

  //! initial data of the problem  
  const Schock& initialData () const; 

  //! return exact Solution 
  const Schock& exactSolution () const; 
    
  //! velocity of the transport problem
  bool kreis ( const FieldVectorType &x) const;
  
  //! velocity of the transport problem
  void v ( const FieldVectorType &x, FieldVectorType & velo) const;
  
  //! do the boundary 
  void boundary ( const FieldVectorType &x, const RangeType & u, RangeType & ret ) const;
  
private:  
  // exact Solution 
  const Schock schock_;
}; // end class RiemannProblem


//**********************************************************************
//
//  Implementation 
//
//**********************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
inline RiemannProblem<dim,dimworld,FunctionSpaceType>::RiemannProblem
(FunctionSpaceType &f ) : LinearTransport <dim,dimworld,FunctionSpaceType>(f) 
  , schock_ ( f )
{
}
  
//! initial data of the problem  
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename RiemannProblem<dim,dimworld,FunctionSpaceType>::Schock& 
RiemannProblem<dim,dimworld,FunctionSpaceType>::initialData () const 
{ 
  return schock_;   
}

//! return exact Solution 
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename RiemannProblem<dim,dimworld,FunctionSpaceType>::Schock& 
RiemannProblem<dim,dimworld,FunctionSpaceType>::exactSolution () const 
{
  return schock_; 
}
    
//! velocity of the transport problem
template <int dim, int dimworld, class FunctionSpaceType> 
inline void RiemannProblem<dim,dimworld,FunctionSpaceType>:: 
v ( const FieldVectorType &x, FieldVectorType & velo) const
{
  velo = 0.0;
  velo[0] = 1.0;
}

//! velocity of the transport problem
template <int dim, int dimworld, class FunctionSpaceType> 
inline bool RiemannProblem<dim,dimworld,FunctionSpaceType>:: 
kreis ( const FieldVectorType &x) const
{
  FieldVectorType z (15.0);
  z[0] = 150.0;
  z -= x;
  if( std::abs(x.two_norm()) < 5.)
    return true;

  return false;
}

//! do the boundary 
template <int dim, int dimworld, class FunctionSpaceType> 
inline void RiemannProblem<dim,dimworld,FunctionSpaceType>::
boundary ( const FieldVectorType &x, const RangeType & u , RangeType & ret ) const 
{
  // return u, for no flow boundary 
  ret = u;
  return;
}

//*******************************************************************
//
// class Schock 
//
//*******************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
inline RiemannProblem<dim,dimworld,FunctionSpaceType>::Schock::
Schock ( FunctionSpaceType &f) : 
 Function < FunctionSpaceType , Schock > ( f )  
 , schock_ (0.0)  
{
  schock_[0] = 1.0;
}  

template <int dim, int dimworld, class FunctionSpaceType> 
inline void RiemannProblem<dim,dimworld,FunctionSpaceType>::Schock::
evaluate (const DomainType & x , DomainFieldType time, RangeType & ret) const 
{
  ret = -1.0;
  if(x[0] - time < - 0.5 ) ret = 1.0;
  //ret = 0.5;
  return;
}

template <int dim, int dimworld, class FunctionSpaceType> 
inline void RiemannProblem<dim,dimworld,FunctionSpaceType>::Schock::
evaluate (const DomainType & x1 , RangeType & ret) const 
{
  evaluate (x1,0.0,ret);
}

#endif

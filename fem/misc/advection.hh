#ifndef __ADVECTION_HH__
#define __ADVECTION_HH__

// Include general headers
#include <iostream>

// Include local headers
#include "transport.hh"

using namespace Dune;

namespace Siemens {
  
  template <int dim, int dimworld, class FunctionSpaceType>
  class LinearAdvection : 
    public LinearTransport<dim, dimworld, FunctionSpaceType>
  {
    // useful typedefs 
    typedef typename FunctionSpaceType::Range RangeType;
    typedef typename FunctionSpaceType::Domain DomainType;
    typedef typename FunctionSpaceType::DomainField DomainFieldType;

    enum { dimRange = RangeType::dimension };
  private:  
    // flux function, describing the problem
    typedef typename FunctionSpaceType::RangeField FieldType;
    typedef typename LinearTransport<dim,dimworld,FunctionSpaceType>::FieldVectorType FieldVectorType;
    
    // the initial data and exact solution 
    class LinAdSolution : 
      public Function < FunctionSpaceType , LinAdSolution > 
    {
    public:
      // Constructor 
      LinAdSolution (FunctionSpaceType &f);
     
      // ??
      //      bool kreis(const FieldVectorType& v) const;

      // evaluate the function
      void evaluate (const DomainType & x,
                     DomainFieldType time, 
                     RangeType & ret) const ;
      
      // evaluate the function a time zero
      void evaluate (const DomainType & x, 
                     RangeType & ret) const; 
      
    };

  public:
    //! Constructor 
    LinearAdvection(FunctionSpaceType &f);
    
    //! Destructor 
    ~LinearAdvection();
    
    //! initial data of the problem  
    const LinAdSolution& initialData () const; 
    
    //! initial data of the problem  
    const LinAdSolution& exactSolution () const; 
    
    //! velocity of the transport problem
    void v ( const FieldVectorType &x, FieldVectorType & velo) const;
    
  //! dispersion  coefficient
    template <class LocalFunctionType> 
    void q ( LocalFunctionType &lf, const FieldVectorType &x, 
             RangeType &ret) const;
    
    //! do the boundary 
    void boundary ( const FieldVectorType &x, const RangeType & u, RangeType & ret ) const;
    
  private:  
    //! the flux function corresponding to this problem
    const LinAdSolution linAd_;
  }; // end class LinearAdvection


//**********************************************************************
//
//  Implementation 
//
//**********************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
inline LinearAdvection<dim,dimworld,FunctionSpaceType>::LinearAdvection(FunctionSpaceType &f) : 
  LinearTransport<dim,dimworld,FunctionSpaceType> (f),
  linAd_( f ) 
{
  this->D_ = 0.0;
}

//! initial data of the problem  
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename LinearAdvection<dim,dimworld,FunctionSpaceType>::LinAdSolution& 
LinearAdvection<dim,dimworld,FunctionSpaceType>::initialData () const 
{ 
  return linAd_;   
}

//! return exact Solution 
template <int dim, int dimworld, class FunctionSpaceType> 
inline const typename LinearAdvection<dim,dimworld,FunctionSpaceType>::LinAdSolution& 
LinearAdvection<dim,dimworld,FunctionSpaceType>::exactSolution () const 
{
  return linAd_; 
}
    
//! velocity of the transport problem
template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearAdvection<dim,dimworld,FunctionSpaceType>:: 
v ( const FieldVectorType &x, FieldVectorType & velo) const
{
  assert(dim == 2);
  velo[0] = 1.0;
  velo[1] = 0.0;
  return;
}

// what the heck!?
//template <int dim, int dimworld, class FunctionSpaceType> 
//inline void LinearAdvection<dim,dimworld,FunctionSpaceType>:: 
//kreis(const FieldVectorType& x) const
//{
//return true;
//}

//! dispersion  coefficient
template <int dim, int dimworld, class FunctionSpaceType> 
template <class LocalFunctionType> 
inline void LinearAdvection<dim,dimworld,FunctionSpaceType>::
q ( LocalFunctionType &lf, const FieldVectorType &x, RangeType & ret ) const 
{
  ret = 0.0;  
}

//! do the boundary 
template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearAdvection<dim,dimworld,FunctionSpaceType>::
boundary ( const FieldVectorType &x, const RangeType & u , RangeType & ret ) const 
{
  linAd_.evaluate(x, this->time_, ret);
}

//*******************************************************************
//
// class LinAdSolution 
//
//*******************************************************************
template <int dim, int dimworld, class FunctionSpaceType> 
inline LinearAdvection<dim,dimworld,FunctionSpaceType>::LinAdSolution::
LinAdSolution ( FunctionSpaceType &f ) : 
 Function < FunctionSpaceType , LinAdSolution > ( f ) 
{
}  

template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearAdvection<dim,dimworld,FunctionSpaceType>::LinAdSolution::
evaluate (const DomainType & x , DomainFieldType time, RangeType & ret) const 
{
  ret = (x[0] < time) ? 1.0 : 0.0;
}

template <int dim, int dimworld, class FunctionSpaceType> 
inline void LinearAdvection<dim,dimworld,FunctionSpaceType>::LinAdSolution::
evaluate (const DomainType & x1 , RangeType & ret) const 
{
  ret = 0.0;
}


} // end namespace Siemens

#endif

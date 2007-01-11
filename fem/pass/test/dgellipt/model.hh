
/***********************************************************************************************
 Model is the interface class used for the discretization of the problem
***********************************************************************************************/
#ifndef __LDGEXAMPLE_MODEL_HH__
#define __LDGEXAMPLE_MODEL_HH__

using namespace Dune;

#include <math.h>

#include "problem.cc"

// include function space
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/discretefunction/common/function.hh>

namespace LDGExample { 

template <typename Field, class GridImp, int dimR>
struct ModelParam
{
  enum { dim      = GridImp::dimension };
  enum { dimworld = GridImp::dimensionworld };
  
  enum { dimRange  = 1 };
  enum { dimDomain = dimworld };

  typedef Field FieldType;
  
  typedef GridImp GridType;
};

template <class Grid>
  class ModelTraits {
  public:
    enum { dimRange2 = 1 };
    enum { dimRange1= dimRange2*Grid::dimensionworld };
    typedef Grid GridType;
    enum { dimDomain = GridType::dimensionworld };
    enum { dimRange = dimRange2 }; 
    enum { dimGradRange = dimRange1 };
    typedef FieldVector<double, dimDomain> DomainType;
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
    typedef FieldVector<double,dimRange> RangeType;
    typedef FieldVector<double,dimGradRange> GradientType;
    typedef FieldMatrix<double,dimDomain+1,dimDomain> FluxRangeType;
    typedef FieldMatrix<double,dimGradRange,dimDomain> DiffusionRangeType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
   };

//**********************************************************************
//
// capillary pressure model 
//
//**********************************************************************
template <class ModelParamType> 
class Model
{
public:
  enum {dim      = ModelParamType :: dim };
  enum {dimworld = ModelParamType :: dimworld };

  enum {dimRange  = ModelParamType :: dimRange };
  enum {dimDomain = ModelParamType :: dimDomain };
  enum {dimGradRange = dimRange * dimDomain }; 

  //! boundary types that may be used by the model
  enum BndType { Dirichlet, Neumann, NoFlow, OutFlow };

  typedef typename ModelParamType :: GridType          GridType; 
  typedef typename ModelParamType :: FieldType         FieldType; 

  typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
  typedef FunctionSpace < FieldType , FieldType, dim , dimGradRange > GradFuncSpaceType;

  typedef typename FuncSpaceType :: RangeFieldType RangeFieldType;
  
  typedef typename FuncSpaceType::RangeType                RangeType;
  typedef typename FuncSpaceType::DomainType               DomainType;
  typedef typename FuncSpaceType::JacobianRangeType        JacobianRangeType;
  
  typedef typename GradFuncSpaceType::RangeType            GradRangeType;
  typedef typename GradFuncSpaceType::DomainType           GradDomainType;
  typedef typename GradFuncSpaceType::JacobianRangeType    GradJacobianRangeType;
  
  typedef typename GridType::template Codim<0>::Entity Entity;

public:
  struct Traits
  {
    typedef typename ModelParamType :: GridType GridType;
    enum {dimRange = ModelParamType :: dimRange };
    enum {dimDomain = dim };
    enum {dimGradRange = dimDomain * dimRange };
    typedef FieldVector< FieldType, dimDomain-1 > FaceDomainType;
    typedef typename GridType::template Codim<0>::Entity Entity;

    typedef FunctionSpace < FieldType , FieldType, dim , dimRange > FuncSpaceType;
    typedef typename FuncSpaceType::RangeType          RangeType;
    typedef typename FuncSpaceType::DomainType         DomainType;
    typedef typename FuncSpaceType::JacobianRangeType  JacobianRangeType;
    typedef FieldMatrix<double,dimDomain+1,dimDomain> FluxRangeType;
  };

public:
  //! constructor 
  Model() 
    : funcSpace_() // fake id = 1
    , initial_ ( funcSpace_ )
    , rhsData_ ( funcSpace_ )
  {
  }

  template <class IntersectionIterator, class FaceDomainType> 
  bool hasBoundaryValue(const IntersectionIterator& it,
             const double time,
             const FaceDomainType& x) const 
  {
    return true;
  }

  template <class IntersectionIterator, class FaceDomainType> 
  void boundaryValue(const IntersectionIterator& it,
           const double time,
           const FaceDomainType& x,
           const RangeType& uLeft,
           RangeType& uRight) const 
  {
    uRight= 0.0;
  }
  
  template <class EntityType , class FaceDomainType> 
  void diffusion(const EntityType & en, const double time ,
                 const FaceDomainType &x ,
                 GradRangeType & ret) const
  {
    ret = exactFactor();
    return ;
  }

  template <class EntityType , class FaceDomainType> 
  void diffusion(const EntityType & en, const double time ,
                 const FaceDomainType &x, 
                 GradJacobianRangeType & ret) const
  {
    ret = 0.0; 
    for(int i=0; i<dimGradRange; ++i) 
    {
      ret[i][i] = exactFactor();
    }
    return ;
  }

  template <class EntityType , class FaceDomainType> 
  void diffusion(const EntityType & en, const double time ,
                 const FaceDomainType &x, 
                 JacobianRangeType & ret) const
  {
    ret = exactFactor(); 
    return ;
  }


  template <class EntityType, class ArgumentTuple,
           class RanType>
  void source(EntityType& en, double time, const DomainType& x,
              const ArgumentTuple& u,
              RanType& s) const
  {
    assert(false);
  }


  template <class EntityType, class DomType, class RanType, 
            class DiffusionRangeType>
  void gradient(const EntityType& en,
             const double time,
             const DomType& x,
             const RanType& s,
             DiffusionRangeType& a) const
  {
    enum { rows = DiffusionRangeType :: rows };
    enum { cols = DiffusionRangeType :: cols };

    assert( (int) rows == (int) dimDomain );
    assert( (int) cols == (int) dimDomain );

    a=0.0;
    for (int i=0; i<dimDomain; ++i)
    {
      a[i][i] = -s[0];
    }
  }
  
  //! interface methods for the initial data
  class InitialData : public Function < FuncSpaceType , InitialData > {
  public:  
    InitialData ( FuncSpaceType &f)
      : Function < FuncSpaceType , InitialData > ( f ) { } ; 

    void evaluate (const DomainType & x , RangeType & ret) const 
    { 
      ret = 1.0;
      for(int i=0; i<DomainType::dimension; i++)
        ret *= ( x[i] - SQR(x[i]) );
      return;
    }

    void evaluate (const DomainType & x , double time, RangeType & ret) const 
    { evaluate (x,ret); return;};  

  };
  typedef InitialData InitialDataType;
  
  const InitialData & initialData () const {
    return initial_;
  } 
    
  //! interface methods for the initial data
  class RHSData : public Function < FuncSpaceType , RHSData > {
  public:  
    RHSData ( FuncSpaceType &f)
      : Function < FuncSpaceType , RHSData > ( f ) { } ; 

    void evaluate (const DomainType & x , RangeType & ret) const 
    { 
      enum { dimR = RangeType :: dimension };
      ret = 0.0;
      ret[dimR - 1] = rhsFunction( &x[0] );
      return;
    }

    void evaluate (const DomainType & x , double time, RangeType & ret) const 
    { evaluate (x,ret); return;};  

  };
  typedef RHSData RHSDataType;
  
  const RHSData & rhsData () const {
    return rhsData_;
  } 
    
private:
  FuncSpaceType funcSpace_; 
  InitialDataType initial_;
  RHSDataType rhsData_;
};

} // end namespace LDGExample 
#endif

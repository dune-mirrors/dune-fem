
/***********************************************************************************************
 Model is the interface class used for the discretization of the problem
***********************************************************************************************/
#ifndef __LDGEXAMPLE_MODEL_HH__
#define __LDGEXAMPLE_MODEL_HH__

using namespace Dune;

#include <math.h>
#include <dune/fem/io/file/asciiparser.hh>

#include "problem.cc"
#include "benchmark.cc"

// include function space
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/function.hh>

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

  typedef FunctionSpace < typename GridType :: ctype , FieldType, dim , dimRange > FuncSpaceType;
  typedef FunctionSpace < typename GridType :: ctype , FieldType, dim , dimGradRange > GradFuncSpaceType;

  typedef typename FuncSpaceType :: RangeFieldType RangeFieldType;
  typedef typename FuncSpaceType :: DomainFieldType DomainFieldType;
  
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

  typedef DataFunctionIF<dim,DomainFieldType,RangeFieldType> DataFunctionType; 
  DataFunctionType * createData(const std::string & paramFile) const 
  {
    int problem = 1;
    readParameter(paramFile,"Problem",problem);

    double shift = 2.0;
    readParameter(paramFile,"GlobalShift",shift);

    double factor = 1.0;
    readParameter(paramFile,"Factor",factor);

    if( problem == 0 ) 
    {
      return new BenchMark_1<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 1 ) 
    {
      return new BenchMark_1_2<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 2 ) 
    {
      return new BenchMark_2<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 3 ) 
    {
      return new BenchMark_3<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 4 ) 
    {
      return new BenchMark_4<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 5 ) 
    {
      return new BenchMark_5<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 6 ) 
    {
      return new BenchMark_6<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 7 ) 
    {
      return new BenchMark_7<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 8 ) 
    {
      return new SinSin<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 9 ) 
    {
      return new CosCos<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 10 ) 
    {
      return new CastilloProblem<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 11 ) 
    {
      return new InSpringingCorner<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }
    if( problem == 12 ) 
    {
      return new RiviereProblem<dim,DomainFieldType,RangeFieldType> (shift,factor);
    }

    return 0;
  }
public:
  //! constructor 
  Model(std::string paramFile) 
    : funcSpace_() // fake id = 1
    , gradFuncSpace_() // fake id = 1
    , problem_( createData ( paramFile ) )
    , solution_( funcSpace_, problem() )
    , gradient_( gradFuncSpace_, problem() )
    , rhsData_ ( funcSpace_ , problem() )
  {
  }
  
  ~Model() 
  {
    delete problem_;
  }

  const DataFunctionType& problem() const 
  { 
    assert( problem_ );
    return *problem_;
  }
  
  bool boundaryValue(const DomainType & p, 
                     RangeType& val) const 
  {
    return problem().boundaryDataFunction(&p[0],val[0]);
  }
  
  void neumann(const DomainType & p, 
               JacobianRangeType & grad) const 
  {
    return problem().neumann(&p[0],&grad[0][0]);
  }
  
  template <class EntityType , class FaceDomainType, class ReturnType> 
  void diffusion(const EntityType & en, const double time ,
                 const FaceDomainType &x, 
                 ReturnType & ret) const
  {
    DomainType point = en.geometry().global(x);
    double k[dim][dim];
    problem().K(&point[0],k);
    for(int i=0; i<dim; ++i)
      for(int j=0; j<dim; ++j)
        ret[i][j] = k[i][j];
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
  class RHSData : public Function < FuncSpaceType , RHSData > 
  {
    enum { dimDomain = FuncSpaceType :: dimDomain };
    typedef DataFunctionIF<dimDomain,DomainFieldType,RangeFieldType> DataFunctionType;

    const DataFunctionType& func_;
  public:  
    RHSData ( FuncSpaceType &f, const DataFunctionType & data)
      : Function < FuncSpaceType , RHSData > ( f ) 
      , func_(data)
    {}  

    void evaluate (const DomainType & x , RangeType & ret) const 
    { 
      enum { dimR = RangeType :: dimension };
      ret = 0.0;
      ret[dimR - 1] = func_.rhs( &x[0] );
      return;
    }

    void evaluate (const DomainType & x , double time, RangeType & ret) const 
    { 
      evaluate (x,ret); return;
    } 
  };
  typedef RHSData RHSDataType;
  const RHSData & rhsData () const {
    return rhsData_;
  } 

  //! the exact solution to the problem for EOC calculation 
  class ExactSolution : public Function < FuncSpaceType , ExactSolution >
  {
    typedef typename FuncSpaceType::RangeType RangeType;
    typedef typename FuncSpaceType::RangeFieldType RangeFieldType;
    typedef typename FuncSpaceType::DomainType DomainType;

    enum { dimDomain = FuncSpaceType :: dimDomain };
    typedef DataFunctionIF<dimDomain,DomainFieldType,RangeFieldType> DataFunctionType;

    const DataFunctionType& data_;
  public:
    ExactSolution (FuncSpaceType &f, const DataFunctionType& d) 
    : Function < FuncSpaceType , ExactSolution > ( f ) 
    , data_(d)   
    {}

    //! see problem.cc 
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      ret = data_.exact( &x[0] );
    }
    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };

  typedef ExactSolution ExactSolutionType;
  const ExactSolutionType& solution() const { return solution_; }

  //! the exact solution to the problem for EOC calculation 
  template <class FuncSPCType>
  class ExactGradient : public Function < FuncSPCType , ExactGradient<
                        FuncSPCType > >
  {
    typedef typename FuncSPCType::RangeType RangeType;
    typedef typename FuncSPCType::RangeFieldType RangeFieldType;
    typedef typename FuncSPCType::DomainType DomainType;
    enum { dimDomain = FuncSpaceType :: dimDomain };
    typedef DataFunctionIF<dimDomain,DomainFieldType,RangeFieldType> DataFunctionType;

    const DataFunctionType& data_;
  public:
    ExactGradient (const FuncSPCType &f,
        const DataFunctionType& d)
      : Function < FuncSPCType , ExactGradient< FuncSPCType > > ( f ) 
      , data_(d)
      {}

    void evaluate (const DomainType & x , RangeType & ret) const
    {
      data_.gradExact( &x[0], &ret[0] );
    }

    void evaluate (const DomainType & x , double time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };
  typedef ExactGradient<GradFuncSpaceType> ExactGradientType;
  const ExactGradientType& gradient() const { return  gradient_; }

private:
  FuncSpaceType funcSpace_; 
  GradFuncSpaceType gradFuncSpace_; 
  DataFunctionType * problem_;
  ExactSolutionType solution_;
  ExactGradientType gradient_;
  RHSDataType rhsData_;
};

} // end namespace LDGExample 
#endif

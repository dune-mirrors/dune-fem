#ifndef DUNE_LEGENDREDGBASEFUNCTIONS_HH
#define DUNE_LEGENDREDGBASEFUNCTIONS_HH

#include <dune/common/misc.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>

#include <dune/fem/space/dgspace/legendrepoly.hh>

namespace Dune
{
  
  namespace Fem {  

  //Template Meta Programm for evaluating tensorproduct polynomial in arbitrary dimensions
  template<int dim,int i,int PolOrd>
  class Eval
  {
  public:
    template <class FieldVectorType>
    static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int idx)
    {
      assert( x.dimension == dim );
      const int num = idx % (PolOrd+1);
      // LegendrePoly lp=LegendrePoly(num);
      return lp.evaluate(num,x[i-1]) * Eval<dim, i-1, PolOrd>::apply(lp,x,(idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class Eval<dim,0,PolOrd>
  {
  public:
    template <class FieldVectorType>
    static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int idx)
    {
      return 1.0;
    }
  };

  //Template Meta Programm for evaluating the partial derivative of tensorproduct polynomials in arbitrary dimensions
  template<int dim,int i,int PolOrd>
  class EvalD{
  public:
    template <class FieldVectorType>
    static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int j, const int idx)
    {
      const int num=idx%(PolOrd+1);
      if( (i-1) != j )
        return lp.evaluate(num,x[i-1]) * EvalD<dim,i-1,PolOrd>::apply(lp,x, j, (idx-num)/(PolOrd+1));
      else
        return lp.jacobian(num,x[i-1]) * EvalD<dim,i-1,PolOrd>::apply(lp,x, j, (idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class EvalD<dim,0,PolOrd>{
  public:
    template <class FieldVectorType>
    static double apply(const LegendrePoly& lp, const FieldVectorType& x, const int j, const int idx)
    {
      return 1.0;
    }
  };


  //! number of legendre base functions for given polord and dim 
  template <int p, int dim>
  struct NumLegendreBaseFunctions
  {
    enum { numBaseFct = Power_m_p<p+1,dim>::power };
  };

  template < class FunctionSpaceType, int polOrd>
  class LegendreDGBaseFunction :
    public BaseFunctionInterface<FunctionSpaceType>
  {
  protected:
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    enum{ dim=DomainType::dimension };

  public:
    LegendreDGBaseFunction( const int baseNum ) 
      : lp_( LegendrePoly :: instance() ),
        baseNum_( baseNum )
    {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }
    
    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = Eval<dim,dim,polOrd>::apply( lp_, x, baseNum_ );
    }

    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = EvalD<dim,dim,polOrd>::apply( lp_, x, diffVariable[0], baseNum_);
    
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      assert(false); // Not implemented
      abort();
    }

    static int numBaseFunctions() 
    {
      return NumLegendreBaseFunctions<polOrd,dim>::numBaseFct;
    }

  protected:
    // reference to legendre polynomials 
    const LegendrePoly& lp_;

    // my number of basis function
    const int baseNum_;
  };

  template <class ScalarFunctionSpaceImp, int polOrd>
  class LegendreDGBaseFunctionFactory : 
    public BaseFunctionFactory<ScalarFunctionSpaceImp> 
  {
  public:
    // true if hierarchical ordering of basis functions 
    // should be used 
    static const bool hierarchical = true ;

    // Add compile time checker: only scalar functions allowed
    typedef ScalarFunctionSpaceImp FunctionSpaceType;
    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;
  public:
    //! constructor creating basis function factory 
    LegendreDGBaseFunctionFactory(const GeometryType& geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>( geo ) {}

    virtual BaseFunctionType* baseFunction(int i) const 
    {
      // only for cubes we have LegendreBaseFunctions 
      if( ! this->geometry().isCube() )
      {
        DUNE_THROW(NotImplemented,"LegendreBaseFunctions only implemented for cubes!");
      }
      return new LegendreDGBaseFunction<FunctionSpaceType ,polOrd> ( i );
    }
    
    virtual int numBaseFunctions() const 
    {
      return NumLegendreBaseFunctions<polOrd,FunctionSpaceType::dimDomain>::numBaseFct;
    }
  };

  } // end namespace Fem 

} // end namespace Dune
#endif

#ifndef DUNE_LEGENDREDGBASEFUNCTIONS_HH
#define DUNE_LEGENDREDGBASEFUNCTIONS_HH

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/common/misc.hh>
// Local includes
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include "legendrepoly.hh"
namespace Dune {
  
  //Template Meta Programm for evaluating tensorproduct polynomial in arbitrary dimensions
  template<int dim,int i,int PolOrd>
  class Eval{
  public:
    static double apply(const LegendrePoly& lp,FieldVector<double,dim> x,int idx){
      int num=idx%(PolOrd+1);
      // LegendrePoly lp=LegendrePoly(num);
      return lp.eval(num,x[i-1])*Eval<dim,i-1,PolOrd>::apply(lp,x,(idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class Eval<dim,0,PolOrd>{
  public:
   static double apply(const LegendrePoly& lp,FieldVector<double,dim> x,int idx){
     return 1.0;
   }
  };
//Template Meta Programm for evaluating the partial derivative of tensorproduct polynomials in arbitrary dimensions
template<int dim,int i,int PolOrd>
  class EvalD{
  public:
    static double apply(const LegendrePoly& lp,FieldVector<double,dim> x,int j,int idx){
      int num=idx%(PolOrd+1);
      if((i-1)!=j)
	return lp.eval(num,x[i-1])*EvalD<dim,i-1,PolOrd>::apply(lp,x,j,(idx-num)/(PolOrd+1));
      else
	return lp.eval1(num,x[i-1])*EvalD<dim,i-1,PolOrd>::apply(lp,x,j,(idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class EvalD<dim,0,PolOrd>{
  public:
    static double apply(const LegendrePoly& lp,FieldVector<double,dim> x,int j,int idx){
      return 1.0;
    }
  };

  //! number of legendre base functions for given polord and dim 
  template <int p, int dim>
  struct NumLegendreBaseFunctions
  {
    enum { numBaseFct = Power_m_p<p+1,dim>::power };
  };

  namespace {
    const LegendrePoly legPoly;
  }

  template <class FunctionSpaceType,int polOrd>
  class LegendreDGBaseFunction :
    public BaseFunctionInterface<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;
   
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    enum{dim=DomainType::dimension};

  public:
    LegendreDGBaseFunction(int baseNum) :
      baseNum_(baseNum),
      lp(legPoly)
    {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }
    
    ~LegendreDGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = Eval<dim,dim,polOrd>::apply(lp,x,baseNum_);
   
    }

    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = EvalD<dim,dim,polOrd>::apply(lp,x,diffVariable[0],baseNum_);
    
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
      abort();
    }

    static int numBaseFunctions() 
    {
      return NumLegendreBaseFunctions<polOrd,dim>::numBaseFct;
    }
    const LegendrePoly& lp;
  };

  template <class ScalarFunctionSpaceImp, int polOrd>
  class LegendreDGBaseFunctionFactory : 
    public BaseFunctionFactory<ScalarFunctionSpaceImp> 
  {
  public:
    // Add compile time checker: only scalar functions allowed

    typedef ScalarFunctionSpaceImp FunctionSpaceType;
    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;
  public:
    LegendreDGBaseFunctionFactory(GeometryType geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>(geo) {}

    virtual BaseFunctionType* baseFunction(int i) const 
    {
      // only for cubes we have LegendreBaseFunctions 
      if( ! this->geometry().isCube() )
      {
        DUNE_THROW(NotImplemented,"LegendreBaseFunctions only implemented for cubes!");
      }
      return new LegendreDGBaseFunction<FunctionSpaceType ,polOrd> (i);
    }
    
    virtual int numBaseFunctions() const 
    {
      return NumLegendreBaseFunctions<polOrd,FunctionSpaceType::dimDomain>::numBaseFct;
    }
  };

} // end namespace Dune
#endif

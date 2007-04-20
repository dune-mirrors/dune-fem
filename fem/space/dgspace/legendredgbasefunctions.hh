#ifndef DUNE_LEGENDREDGBASEFUNCTIONS_HH
#define DUNE_LEDENDREDGBASEFUNCTIONS_HH

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/common/misc.hh>
// Local includes
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include "legendrepoly.hh"
namespace Dune {
  
  typedef int deriType;

  //Template Meta Programm for evaluating tensorproduct polynomial in arbitrary dimensions
  template<int dim,int i,int PolOrd>
  class Eval{
  public:
    static double apply(FieldVector<double,dim> x,int idx){
      int num=idx%(PolOrd+1);
      LegendrePoly  lp=LegendrePoly(num);
      return lp.eval(x[i-1])*Eval<dim,i-1,PolOrd>::apply(x,(idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class Eval<dim,0,PolOrd>{
  public:
   static double apply(FieldVector<double,dim> x,int idx){
     return 1.0;
   }
  };
//Template Meta Programm for evaluating the partial derivative of tensorproduct polynomials in arbitrary dimensions
template<int dim,int i,int PolOrd>
  class EvalD{
  public:
    static double apply(FieldVector<double,dim> x,int j,int idx){
      int num=idx%(PolOrd+1);
      LegendrePoly  lp=LegendrePoly(num);
      if((i-1)!=j)
	return lp.eval(x[i-1])*EvalD<dim,i-1,PolOrd>::apply(x,j,(idx-num)/(PolOrd+1));
      else
	return lp.eval1(x[i-1])*EvalD<dim,i-1,PolOrd>::apply(x,j,(idx-num)/(PolOrd+1));
    }
  };

  template<int dim,int PolOrd>
  class EvalD<dim,0,PolOrd>{
  public:
    static double apply(FieldVector<double,dim> x,int j,int idx){
      return 1.0;
    }
  };
 

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
      baseNum_(baseNum) {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }
    
    ~LegendreDGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = Eval<dim,dim,polOrd>::apply(x,baseNum_);
   
    }

    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = EvalD<dim,dim,polOrd>::apply(x,diffVariable[0],baseNum_);
    
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
      abort();
    }
    static int numBaseFunctions() 
    {
      return Power_m_p<polOrd+1,dim>::power;
    }

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
      GeometryType type = this->geometry();
      assert(type.isCube());
	
      return new LegendreDGBaseFunction<FunctionSpaceType ,polOrd>(i);
      
      DUNE_THROW(NotImplemented, 
                 "The chosen geometry type is not implemented");
      return 0;
    }
    
    virtual int numBaseFunctions() const 
    {
      return Power_m_p<polOrd+1,FunctionSpaceType::DimDomain>::power;
    }
  };

} // end namespace Dune

//#include "legendre_imp.cc"
#endif

#ifndef DUNE_BASEFUNCTIONSTUB_HH
#define DUNE_BASEFUNCTIONSTUB_HH

#include <cmath>

// Dune includes
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>

namespace Dune {

  typedef int deriType;

  template <class FunctionSpaceImp>
  class BaseFunctionStub : 
    public BaseFunctionInterface<FunctionSpaceImp>
  {
  public:
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;

  public:
    BaseFunctionStub(int baseNum) :
      baseNum_(static_cast<double>(baseNum+1))
    {}

    ~BaseFunctionStub() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = 0;
      std::cout << "eval(" << baseNum_ << "): x = " << x << std::endl;
      for (int i = 0; i < DomainType::size; ++i) {
        phi += std::pow(x[i], i);
      }
      phi *= baseNum_;
      std::cout << "eval = " << phi << std::endl;
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const 
    {
      switch (diffVar[0]) {
      case 0:
        phi = 0.;
        break;
      case 1:
        phi = 1.;
        break;
      case 2:
        phi = 2.*x[2];
        break;
      default:
        assert(false);
      }
      phi *= baseNum_;
    }
        
    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }
  
  private:
    double baseNum_;
  };

  template <class FunctionSpaceImp>
  class BaseFunctionStubFactory :
    public BaseFunctionFactory<FunctionSpaceImp>
  {
  public:
    typedef BaseFunctionInterface<FunctionSpaceImp> BaseFunctionType;

  public:
    BaseFunctionStubFactory(GeometryType geo) :
      BaseFunctionFactory<FunctionSpaceImp>(geo)
    {}

    virtual BaseFunctionType* baseFunction(int i) const 
    {
      return new BaseFunctionStub<FunctionSpaceImp>(i);
    }

    virtual int numBaseFunctions() const {
      return 2;
    }

  private:
  };
  
}

#endif

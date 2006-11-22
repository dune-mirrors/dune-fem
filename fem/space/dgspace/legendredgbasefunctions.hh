#ifndef DUNE_LEGENDREDGBASEFUNCTIONS_HH
#define DUNE_LEDENDREDGBASEFUNCTIONS_HH

// Dune includes
#include <dune/grid/common/grid.hh>

// Local includes
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>

namespace Dune {
  
  typedef int deriType;

  //! Wrapper interface for DG base functions
  template <class FunctionSpaceType>
  class LegendreDGBaseFunctionWrapper {
  protected:
    LegendreDGBaseFunctionWrapper() {}
    virtual ~LegendreDGBaseFunctionWrapper() {}
    
    enum { dimDomain = FunctionSpaceType::DimDomain };

    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    static int numBaseFunctions(int polOrder) 
    {
      switch (dimDomain) {
      case 2:
        return (polOrder+1)*(polOrder+1);
      case 3:
        return (polOrder+1)*(polOrder+1)*(polOrder+1);
      default:

        DUNE_THROW(NotImplemented, "DGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }
   
    double eval_quadrilateral_2d_l (int j,int i, const DomainType & xi ) const;
    double eval_hexahedron_3d_l (int j,int i, const DomainType & xi ) const;
    void grad_quadrilateral_2d_l (int j,int i, const DomainType & xi,
                                 JacobianRangeType & grad ) const;
    void grad_hexahedron_3d_l (int j,int i, const DomainType & xi,
             JacobianRangeType & grad ) const;

  }; // end class DGBaseFunctionWrapper

  //! Base class for DG base functions
  template <class FunctionSpaceType, GeometryIdentifier::IdentifierType ElType, int polOrd>
  class LegendreDGBaseFunction;

  //! Specialisation for quadrilaterals 
  template <class FunctionSpaceType, int polOrd>
  class LegendreDGBaseFunction<FunctionSpaceType,
        GeometryIdentifier::Quadrilateral , polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private LegendreDGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;
    
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    LegendreDGBaseFunction(int baseNum) :
      baseNum_(baseNum) {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~LegendreDGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_quadrilateral_2d_l(polOrd,baseNum_,x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_quadrilateral_2d_l(polOrd,baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }
    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
      abort();
    }

    static int numBaseFunctions() {
      return LegendreDGBaseFunctionWrapper<FunctionSpaceType>::numBaseFunctions(static_cast<int>(polOrd));
    }
  }; // end class DGBaseFunction<FunctionSpaceType, quadrilateral, polOrd>


  //! Specialisation for hexahedrons 
  template <class FunctionSpaceType, int polOrd>
  class LegendreDGBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private LegendreDGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;
    //static const bool leg_=false;
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    LegendreDGBaseFunction(int baseNum) :
      //DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~LegendreDGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      
      phi = this->eval_hexahedron_3d_l(polOrd,baseNum_,x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_hexahedron_3d_l(polOrd,baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }
    
    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
      abort();
    }

    static int numBaseFunctions() 
    {
      return LegendreDGBaseFunctionWrapper<FunctionSpaceType>:: numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class LegendreDGBaseFunction<FunctionSpaceType, hexahedron, polOrd>

  /*
  template <class FunctionSpaceImp, GeometryIdentifier::IdentifierType elType, int polOrd>
  class LegendreDGFastBaseFunctionSet :
    public FastBaseFunctionSet<FunctionSpaceImp>
  {
    enum { dimRange = FunctionSpaceImp::DimRange };
   
    typedef LegendreDGBaseFunction<
      FunctionSpaceImp, elType, polOrd> MyBaseFunctionType;

  public:
    LegendreDGFastBaseFunctionSet(FunctionSpaceImp& spc) :
      FastBaseFunctionSet<FunctionSpaceImp>(spc, MyBaseFunctionType::numBaseFunctions()) 
    {
      assert(dimRange == 1); // works only for scalar spaces
      int numBaseFct = MyBaseFunctionType::numBaseFunctions();
      this->setNumOfDiffFct(numBaseFct);

      for (int i = 0; i < numBaseFct; ++i) {
        this->setBaseFunctionPointer(i, new MyBaseFunctionType(i));
      }
    }

    virtual ~LegendreDGFastBaseFunctionSet() {}
  };
  */

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
      if(type.isCube())
      {
        if(type.isQuadrilateral())
          return new LegendreDGBaseFunction<FunctionSpaceType, GeometryIdentifier::Quadrilateral,polOrd>(i);
        if(type.isHexahedron())
          return new LegendreDGBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd>(i);
      } 
      DUNE_THROW(NotImplemented, 
                 "The chosen geometry type is not implemented");
      return 0;
    }

    virtual int numBaseFunctions() const {
      switch (FunctionSpaceType::DimDomain) {
      case 2:
        return  (polOrd+1)*(polOrd+1);
      case 3:
        return ((polOrd+1)*(polOrd+1)*(polOrd+1));
      default:
        DUNE_THROW(NotImplemented, 
                   "LegendreDGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      abort();
      return -1;
    }
  };

} // end namespace Dune

#include "legendre_imp.cc"
#endif

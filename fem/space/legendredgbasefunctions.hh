#ifndef DUNE_DGBASEFUNCTIONS_HH
#define DUNE_DGBASEFUNCTIONS_HH

// Dune includes
#include <dune/fem/common/basefunctions.hh>
#include <dune/fem/common/basefunctionfactory.hh>
#include <dune/fem/common/fastbase.hh>
#include <dune/grid/common/grid.hh>

// Local includes

using namespace Dune;

namespace Dune {
  
  typedef int deriType;

  //! Wrapper interface for DG base functions
  template <class FunctionSpaceType>
  class DGBaseFunctionWrapper {
  protected:
    DGBaseFunctionWrapper() {}
    virtual ~DGBaseFunctionWrapper() {}
    
    enum { dimDomain = FunctionSpaceType::DimDomain };

    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    static int numBaseFunctions(int polOrder) {
      switch (dimDomain) {
      case 2:
        return (polOrder + 2) * (polOrder + 1) / 2;
      case 3:
        return ((polOrder+1)*(polOrder+2)*(2*polOrder+3)/6 +
                         (polOrder+1)*(polOrder+2)/2)/2;
      default:
        DUNE_THROW(NotImplemented, "DGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }

    static int numBaseFunctions_l(int polOrder) {
      switch (dimDomain) {
      case 2:
        return (polOrder + 1) * (polOrder + 1);
      case 3:
        return (polOrder+1)*(polOrder+1)*(polOrder+1);
      default:

        DUNE_THROW(NotImplemented, "DGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }
    
    double eval_triangle_2d (int i, const DomainType & xi ) const;
    double eval_quadrilateral_2d (int polOrd,int i, const DomainType & xi ) const;//for legendrePolys the polOrd is needed to construct the right polynomial
    double eval_quadrilateral_2d (int i, const DomainType & xi ) const;
    double eval_tetrahedron_3d (int i, const DomainType & xi ) const;
    double eval_pyramid_3d (int i, const DomainType & xi ) const;
    double eval_prism_3d (int i, const DomainType & xi ) const;
    double eval_hexahedron_3d (int polOrd,int i, const DomainType & xi ) const;
    double eval_hexahedron_3d (int i, const DomainType & xi ) const;
    
    void grad_triangle_2d (int i, const DomainType & xi,
                                    JacobianRangeType & grad ) const;
    void grad_quadrilateral_2d (int polOrd,int i, const DomainType & xi,
                                 JacobianRangeType & grad ) const;
    void grad_quadrilateral_2d (int i, const DomainType & xi,
                                 JacobianRangeType & grad ) const;
    void grad_tetrahedron_3d (int i, const DomainType & xi,
                                  JacobianRangeType & grad ) const;
    void grad_pyramid_3d (int i, const DomainType & xi,
                             JacobianRangeType & grad ) const;
    void grad_prism_3d (int i, const DomainType & xi,
			JacobianRangeType & grad ) const;
    void grad_hexahedron_3d (int i, const DomainType & xi,
                              JacobianRangeType & grad ) const;
    void grad_hexahedron_3d (,int i, const DomainType & xi,
                              JacobianRangeType & grad ) const;
  


//legendre evaluate
    double eval_quadrilateral_2d_l (int j,int i, const DomainType & xi ) const;
    double eval_hexahedron_3d_l (int j,int i, const DomainType & xi ) const;
    void grad_quadrilateral_2d_l (int j,int i, const DomainType & xi,
                                 JacobianRangeType & grad ) const;
    void grad_hexahedron_3d_l (int j,int i, const DomainType & xi,
			       JacobianRangeType & grad ) const;






  }; // end class DGBaseFunctionWrapper

  //! Base class for DG base functions
  template <class FunctionSpaceType, GeometryType ElType, int polOrd>
  class LegendreDGBaseFunction;



  //! Specialisation for quadrilaterals 
template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, quadrilateral, polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;
    
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
     
	phi = this->eval_quadrilateral_2d_l(polOrd,baseNum_,x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;

      i
	  this->grad_quadrilateral_2d_l(polOrd,baseNum_, x, tmp);
	  phi = tmp[0][diffVar[0]];

    }
    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      
	return DGBaseFunctionWrapper<FunctionSpaceType>::numBaseFunctions_l(static_cast<int>(polOrd));
    }
  }; // end class DGBaseFunction<FunctionSpaceType, quadrilateral, polOrd>



  



//! Specialisation for hexahedrons 
  template <class FunctionSpaceType, int polOrd>
  class legendreDGBaseFunction<FunctionSpaceType, hexahedron, polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;
    static const bool leg_=false;
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      
	pfi = this->eval_hexahedron_3d_l(polOrd,baseNum_,x);
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
    }

    static int numBaseFunctions() {
      if(!leg_)
	return DGBaseFunctionWrapper<FunctionSpaceType>:: numBaseFunctions(static_cast<int>(polOrd));
      else
	return DGBaseFunctionWrapper<FunctionSpaceType>:: numBaseFunctions_l(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, hexahedron, polOrd>


 





  template <class FunctionSpaceImp, GeometryType elType, int polOrd>
  class DGFastBaseFunctionSet :
    public FastBaseFunctionSet<FunctionSpaceImp>
  {
    enum { dimRange = FunctionSpaceImp::DimRange };
   
    typedef LegendreDGBaseFunction<
      FunctionSpaceImp, elType, polOrd> DGBaseFunctionType;

  public:
    DGFastBaseFunctionSet(FunctionSpaceImp& spc) :
      FastBaseFunctionSet<FunctionSpaceImp>(spc, DGBaseFunctionType::numBaseFunctions()) 
    {
      assert(dimRange == 1); // works only for scalar spaces
      int numBaseFct = DGBaseFunctionType::numBaseFunctions();
      this->setNumOfDiffFct(numBaseFct);

      for (int i = 0; i < numBaseFct; ++i) {
        this->setBaseFunctionPointer(i, new DGBaseFunctionType(i));
      }
    }

    virtual ~DGFastBaseFunctionSet() {}
      
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
    LegendreBaseFunctionFactory(GeometryType geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>(geo)
    {}

    virtual BaseFunctionType* baseFunction(int i) const {
      switch (this->geometry()) {
      
      case cube:
        if (FunctionSpaceType::DimDomain == 2) {
          return new LegendreBaseFunction<FunctionSpaceType,quadrilateral,polOrd>(i);
        }
        else {
          return new LegendreBaseFunction<FunctionSpaceType, hexahedron, polOrd>(i);
        }
      
      case quadrilateral:
        return new LegendreBaseFunction<FunctionSpaceType, quadrilateral, polOrd>(i);
      
      case hexahedron:
        return new LegendreBaseFunction<FunctionSpaceType, hexahedron, polOrd>(i);
      
      default:
        DUNE_THROW(NotImplemented, 
                   "The chosen geometry type is not implemented");
      }
      return 0;
    }

    virtual int numBaseFunctions() const {
      switch (FunctionSpaceType::DimDomain) {
      case 2:
        return (polOrd + 1) * (polOrd + 1);
      case 3:
        return ((polOrd+1)*(polOrd+1)*(polOrd+1));
      default:
        DUNE_THROW(NotImplemented, 
                   "DGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }
    
  };




#include "orthonormalbase_mod.cc"
#include "quadraevaluate.cc"
} // end namespace Dune

#endif

#ifndef DUNE_DGBASEFUNCTIONS_HH
#define DUNE_DGBASEFUNCTIONS_HH

//- Dune includes 
#include <dune/grid/common/grid.hh>

//- Local includes
#include "../../basefunctions/common/basefunctions.hh"
#include "../../basefunctions/common/basefunctionfactory.hh"
#include "../../basefunctions/common/fastbase.hh"

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
      case 1:
	return (polOrder + 1);
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
    double eval_line(int i,const DomainType & xi ) const;
    double eval_triangle_2d (int i, const DomainType & xi ) const;
    double eval_quadrilateral_2d (int i, const DomainType & xi ) const;
    double eval_tetrahedron_3d (int i, const DomainType & xi ) const;
    double eval_pyramid_3d (int i, const DomainType & xi ) const;
    double eval_prism_3d (int i, const DomainType & xi ) const;
    double eval_hexahedron_3d (int i, const DomainType & xi ) const;
   
    void grad_line(int i, const DomainType & xi,
                             JacobianRangeType & grad ) const;
    void grad_triangle_2d (int i, const DomainType & xi,
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
  }; // end class DGBaseFunctionWrapper

  //! Base class for DG base functions
  template <class FunctionSpaceType, GeometryIdentifier::IdentifierType ElType, int polOrd>
  class DGBaseFunction;

  //! Specialisation for triangles
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Line, polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;

  private:
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
      phi = this->eval_line(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_line(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, triangle, polOrd>

  //! Specialisation for triangles
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Triangle, polOrd> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType>
  {
  private:
    //- Local data
    int baseNum_;

  private:
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
      phi = this->eval_triangle_2d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_triangle_2d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, triangle, polOrd>

  //! Specialisation for quadrilaterals
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Quadrilateral, polOrd> :
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
      phi = this->eval_quadrilateral_2d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_quadrilateral_2d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }
  }; // end class DGBaseFunction<FunctionSpaceType, quadrilateral, polOrd>

  //! Specialisation for tetrahedrons
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Tetrahedron, polOrd> :
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
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_tetrahedron_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_tetrahedron_3d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, tetrahedron, polOrd>

  //! Specialisation for pyramids
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Pyramid, polOrd> :
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
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_pyramid_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_pyramid_3d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, pyramid, polOrd>

  //! Specialisation for prisms
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Prism, polOrd> :
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
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_prism_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_prism_3d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]];
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, prism, polOrd>

  //! Specialisation for hexahedrons
  template <class FunctionSpaceType, int polOrd>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd> :
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
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_hexahedron_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVar,
                          const DomainType& x, RangeType& phi) const {
      JacobianRangeType tmp;
      this->grad_hexahedron_3d(baseNum_, x, tmp);
      phi = tmp[0][diffVar[0]]; 
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVar,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, hexahedron, polOrd>
  
  template <class FunctionSpaceImp, GeometryIdentifier::IdentifierType elType, int polOrd>
  class DGFastBaseFunctionSet :
    public FastBaseFunctionSet<FunctionSpaceImp>
  {
    enum { dimRange = FunctionSpaceImp::DimRange };
   
    typedef DGBaseFunction<
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
  class DiscontinuousGalerkinBaseFunctionFactory : 
    public BaseFunctionFactory<ScalarFunctionSpaceImp> 
  {
  public:
    // Add compile time checker: only scalar functions allowed

    typedef ScalarFunctionSpaceImp FunctionSpaceType;
    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;
  public:
    DiscontinuousGalerkinBaseFunctionFactory(GeometryType geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>(geo)
    {}

    virtual BaseFunctionType* baseFunction(int i) const {
      switch (GeometryIdentifier::fromGeo(this->geometry())) 
      {
	case GeometryIdentifier::Line:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Line, polOrd>(i);
        case GeometryIdentifier::Triangle:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Triangle, polOrd>(i);
        case GeometryIdentifier::Quadrilateral:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Quadrilateral, polOrd>(i);
        case GeometryIdentifier::Tetrahedron:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Tetrahedron, polOrd>(i);
        case GeometryIdentifier::Pyramid:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Pyramid, polOrd>(i);
        case GeometryIdentifier::Prism:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Prism, polOrd>(i);
        case GeometryIdentifier::Hexahedron:
          return new DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd>(i);
        default:
          DUNE_THROW(NotImplemented, 
                     "The chosen geometry type is not implemented");
      }
      return 0;
    }

    virtual int numBaseFunctions() const {
      switch (FunctionSpaceType::DimDomain) {
      case 1:
	return (polOrd + 1);
      case 2:
        return (polOrd + 2) * (polOrd + 1) / 2;
      case 3:
        return ((polOrd+1)*(polOrd+2)*(2*polOrd+3)/6 +
                         (polOrd+1)*(polOrd+2)/2)/2;
      default:
        DUNE_THROW(NotImplemented, 
                   "DGBaseFunctionWrapper only supports 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }
    
  };

} // end namespace Dune

#include "orthonormalbase_mod.cc"

#endif

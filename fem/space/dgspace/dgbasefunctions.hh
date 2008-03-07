#ifndef DUNE_DGBASEFUNCTIONS_HH
#define DUNE_DGBASEFUNCTIONS_HH

//- Dune includes 
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

//- Local includes
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include "orthonormalbase_mod.hh"

namespace Dune {
  
  template <int polOrder, int dimDomain> 
  struct DGNumberOfBaseFunctions
  {
    enum { numBaseFunctions = 0 };
  };
  
  //! number of base functions for dimDomain = 1 
  template <int polOrder> 
  struct DGNumberOfBaseFunctions<polOrder,1> 
  {
    enum { numBaseFunctions = polOrder + 1 };
  };
  
  //! number of base functions for dimDomain = 2 
  template <int polOrder> 
  struct DGNumberOfBaseFunctions<polOrder,2> 
  {
    enum { numBaseFunctions = (polOrder + 2) * (polOrder + 1) / 2};
  };

  //! number of base functions for dimDomain = 3 
  template <int polOrder> 
  struct DGNumberOfBaseFunctions<polOrder,3> 
  {
    enum { numBaseFunctions = ((polOrder+1)
            *(polOrder+2)*(2*polOrder+3)/6 
              + (polOrder+1)*(polOrder+2)/2)/2 };
  };

  //! Wrapper interface for DG base functions
  template <class FunctionSpaceType>
  class DGBaseFunctionWrapper 
  {
  protected:
    DGBaseFunctionWrapper() {}
    virtual ~DGBaseFunctionWrapper() {}
    
    enum { dimDomain = FunctionSpaceType::DimDomain };

    // temporary variable for grad evaluation 
    mutable double grad_[dimDomain];

    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  public:  
    /** \brief return number of base functions depending on polynomial order
               and dimension of domain */
    static int numBaseFunctions(const int polOrder, const int dimension) 
    {
      switch (dimension) 
      {
        case 1:
          return (polOrder + 1);
        case 2:
          return (polOrder + 2) * (polOrder + 1) / 2;
        case 3:
          return ((polOrder+1)*(polOrder+2)*(2*polOrder+3)/6 +
                         (polOrder+1)*(polOrder+2)/2)/2;
        default:
          DUNE_THROW(NotImplemented, "DGBaseFunctionWrapper only supports 1D, 2D and 3D Domain");
      }
      assert(false); // can't get here!
      return -1;
    }
  
  protected:
    // numBaseFunctions for derived classes 
    static int numBaseFunctions(int polOrder) 
    {
      return numBaseFunctions(polOrder, dimDomain);
    }
    ////////////////////////////
    /// 1d functions 
    ////////////////////////////
    // eval function
    RangeFieldType eval_line(const int i,const DomainType & xi ) const
    {
      return OrthonormalBase_1D::eval_line(i,&xi[0]); 
    }

    // eval gradient 
    RangeFieldType grad_line(const int i, const int comp, const DomainType & xi) const
    {
      OrthonormalBase_1D::grad_line(i,&xi[0],
                                     grad_);
      return grad_[comp];
    }

    ///////////////////////////////////
    //  2d functions 
    //////////////////////////////////
    RangeFieldType eval_triangle_2d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_2D::eval_triangle_2d(i,&xi[0]);
    }
    
    RangeFieldType eval_quadrilateral_2d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_2D::eval_quadrilateral_2d(i,
                                                       &xi[0]);  
    }
    
    RangeFieldType grad_triangle_2d (const int i, const int comp, const DomainType & xi) const
    {
      OrthonormalBase_2D::grad_triangle_2d(i,&xi[0],grad_);
      return grad_[comp];
    }
    
    RangeFieldType grad_quadrilateral_2d (const int i, const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_2D::grad_quadrilateral_2d(i, &xi[0], grad_); 
      return grad_[comp];
    }

    //////////////////////////////////////
    //  3d functions 
    //////////////////////////////////////
    RangeFieldType eval_tetrahedron_3d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_tetrahedron_3d(i,&xi[0]); 
    }
    
    RangeFieldType eval_pyramid_3d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_pyramid_3d(i,&xi[0]); 
    }
    
    RangeFieldType eval_prism_3d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_prism_3d(i,&xi[0]); 
    }
    
    RangeFieldType eval_hexahedron_3d (const int i, const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_hexahedron_3d(i,&xi[0]); 
    }
   
    RangeFieldType grad_tetrahedron_3d (const int i, const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_3D::grad_tetrahedron_3d(i,&xi[0],grad_);
      return grad_[comp];
    }
    
    RangeFieldType grad_pyramid_3d (const int i, const int comp, const DomainType & xi) const
    {
      OrthonormalBase_3D::grad_pyramid_3d(i,&xi[0], grad_); 
      return grad_[comp];
    }

    RangeFieldType grad_prism_3d (const int i, const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_3D::grad_prism_3d(i,&xi[0], grad_ ); 
      return grad_[comp];
    }

    RangeFieldType grad_hexahedron_3d (const int i, const int comp, const DomainType & xi) const
    {
      OrthonormalBase_3D::grad_hexahedron_3d(i,&xi[0], grad_ ); 
      return grad_[comp];
    }

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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

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

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_line(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_line(baseNum_, diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

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

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_triangle_2d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_triangle_2d(baseNum_, diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    int baseNum_;

  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) 
    {
      // Check if base number is valid
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      phi = this->eval_quadrilateral_2d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      phi = this->grad_quadrilateral_2d(baseNum_, diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_tetrahedron_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_tetrahedron_3d(baseNum_, diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_pyramid_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_pyramid_3d(baseNum_, diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_prism_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_prism_3d(baseNum_, diffVariable[0], x );
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
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
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  private:
    //- Local data
    const int baseNum_;

  public:
    DGBaseFunction(int baseNum) :
      DGBaseFunctionWrapper<FunctionSpaceType>(),
      baseNum_(baseNum) {
      assert(baseNum_ >= 0 && baseNum_ < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::DimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_hexahedron_3d(baseNum_, x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_hexahedron_3d(baseNum_,diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, hexahedron, polOrd>
  
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

    virtual BaseFunctionType* baseFunction(int i) const 
    {
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

    virtual int numBaseFunctions() const 
    {
      // call numBaseFunctions from Wrapper class depending on 
      // geometry's dimension 
      return DGBaseFunctionWrapper<FunctionSpaceType> 
                :: numBaseFunctions( polOrd , this->geometry().dim() );
    }
    
  };

} // end namespace Dune
#endif

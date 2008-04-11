#ifndef DUNE_DGBASEFUNCTIONS_HH
#define DUNE_DGBASEFUNCTIONS_HH

//- Dune includes 
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

//- Local includes
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/common/functionspace.hh>
#include "orthonormalbase_mod.hh"

namespace Dune
{
  
  // Number of DG Base Functions
  // ---------------------------

  template< int polOrder, int dimDomain >
  struct DGNumberOfBaseFunctions;
  
  template< int polOrder >
  struct DGNumberOfBaseFunctions< polOrder, 1 >
  {
    enum { numBaseFunctions = polOrder + 1 };
  };
  
  template< int polOrder >
  struct DGNumberOfBaseFunctions< polOrder, 2 >
  {
    enum { numBaseFunctions = (polOrder + 2) * (polOrder + 1) / 2};
  };

  template< int polOrder >
  struct DGNumberOfBaseFunctions< polOrder, 3 >
  {
    enum { numBaseFunctions = ((polOrder+1)
            *(polOrder+2)*(2*polOrder+3)/6 
              + (polOrder+1)*(polOrder+2)/2)/2 };
  };



  //! Wrapper interface for DG base functions
  template <class FunctionSpaceType, int baseNum >
  class DGBaseFunctionWrapper 
  {
  protected:
    DGBaseFunctionWrapper() {}
    virtual ~DGBaseFunctionWrapper() {}
    
    enum { dimDomain = FunctionSpaceType::dimDomain };

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
    RangeFieldType eval_line(const DomainType & xi ) const
    {
      return OrthonormalBase_1D::eval_line(baseNum,&xi[0]); 
    }

    // eval gradient 
    RangeFieldType grad_line(const int comp, const DomainType & xi) const
    {
      OrthonormalBase_1D::grad_line(baseNum,&xi[0],grad_);
      return grad_[comp];
    }

    ///////////////////////////////////
    //  2d functions 
    //////////////////////////////////
    RangeFieldType eval_triangle_2d (const DomainType & xi ) const
    {
      return OrthonormalBase_2D::eval_triangle_2d(baseNum,&xi[0]);
    }
    
    RangeFieldType eval_quadrilateral_2d (const DomainType & xi ) const
    {
      return OrthonormalBase_2D::eval_quadrilateral_2d(baseNum,&xi[0]);  
    }
    
    RangeFieldType grad_triangle_2d (const int comp, const DomainType & xi) const
    {
      OrthonormalBase_2D::grad_triangle_2d(baseNum,&xi[0],grad_);
      return grad_[comp];
    }
    
    RangeFieldType grad_quadrilateral_2d (const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_2D::grad_quadrilateral_2d(baseNum, &xi[0], grad_); 
      return grad_[comp];
    }

    //////////////////////////////////////
    //  3d functions 
    //////////////////////////////////////
    RangeFieldType eval_tetrahedron_3d (const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_tetrahedron_3d(baseNum,&xi[0]); 
    }
    
    RangeFieldType eval_pyramid_3d (const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_pyramid_3d(baseNum,&xi[0]); 
    }
    
    RangeFieldType eval_prism_3d (const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_prism_3d(baseNum,&xi[0]); 
    }
    
    RangeFieldType eval_hexahedron_3d (const DomainType & xi ) const
    {
      return OrthonormalBase_3D::eval_hexahedron_3d(baseNum,&xi[0]); 
    }
   
    RangeFieldType grad_tetrahedron_3d (const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_3D::grad_tetrahedron_3d(baseNum,&xi[0],grad_);
      return grad_[comp];
    }
    
    RangeFieldType grad_pyramid_3d (const int comp, const DomainType & xi) const
    {
      OrthonormalBase_3D::grad_pyramid_3d(baseNum,&xi[0], grad_); 
      return grad_[comp];
    }

    RangeFieldType grad_prism_3d (const int comp, const DomainType & xi) const 
    {
      OrthonormalBase_3D::grad_prism_3d(baseNum,&xi[0], grad_ ); 
      return grad_[comp];
    }

    RangeFieldType grad_hexahedron_3d (const int comp, const DomainType & xi) const
    {
      OrthonormalBase_3D::grad_hexahedron_3d(baseNum,&xi[0], grad_ ); 
      return grad_[comp];
    }

  }; // end class DGBaseFunctionWrapper

  //! Base class for DG base functions
  template <class FunctionSpaceType, GeometryIdentifier::IdentifierType ElType, int polOrd, int baseNum>
  class DGBaseFunction;

  //! Specialisation for triangles
  template <class FunctionSpaceType, int polOrd, int baseNum>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Line, polOrd, baseNum> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType, baseNum>
  {
  private:
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      // Check if base number is valid
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_line(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_line(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType, baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, triangle, polOrd>

  //! Specialisation for triangles
  template <class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Triangle, polOrd, baseNum > :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType, baseNum>
  {
  private:
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      // Check if base number is valid
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_triangle_2d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_triangle_2d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, triangle, polOrd>

  //! Specialisation for quadrilaterals
  template <class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Quadrilateral, polOrd, baseNum > :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType, baseNum>
  {
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      // Check if base number is valid
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_quadrilateral_2d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_quadrilateral_2d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }
  }; // end class DGBaseFunction<FunctionSpaceType, quadrilateral, polOrd>

  //! Specialisation for tetrahedrons
  template <class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Tetrahedron, polOrd, baseNum> :
    public BaseFunctionInterface<FunctionSpaceType>,
    private DGBaseFunctionWrapper<FunctionSpaceType, baseNum>
  {
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_tetrahedron_3d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_tetrahedron_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, tetrahedron, polOrd>

  //! Specialisation for pyramids
  template <class FunctionSpaceType, int polOrd, int baseNum>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Pyramid, polOrd, baseNum> :
    public BaseFunctionInterface<FunctionSpaceType>,
     private DGBaseFunctionWrapper<FunctionSpaceType,baseNum>
  {
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_pyramid_3d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_pyramid_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, pyramid, polOrd>

  //! Specialisation for prisms
  template <class FunctionSpaceType, int polOrd, int baseNum>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Prism, polOrd, baseNum> :
    public BaseFunctionInterface<FunctionSpaceType>,
     private DGBaseFunctionWrapper<FunctionSpaceType,baseNum>
  {
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_prism_3d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_prism_3d(diffVariable[0], x );
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, prism, polOrd>

  //! Specialisation for hexahedrons
  template <class FunctionSpaceType, int polOrd, int baseNum>
  class DGBaseFunction<FunctionSpaceType, GeometryIdentifier::Hexahedron, polOrd, baseNum > :
    public BaseFunctionInterface<FunctionSpaceType>,
     private DGBaseFunctionWrapper<FunctionSpaceType,baseNum>
  {
    //- Local typedefs
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  
  public:
    DGBaseFunction() :
      DGBaseFunctionWrapper<FunctionSpaceType,baseNum>()
    {
      assert(baseNum >= 0 && baseNum < numBaseFunctions());
      // Only for scalar function spaces
      assert(FunctionSpaceType::dimRange == 1);
    }

    ~DGBaseFunction() {}

    virtual void evaluate(const FieldVector<deriType, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_hexahedron_3d(x);
    }
    
    virtual void evaluate(const FieldVector<deriType, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_hexahedron_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<deriType, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
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
    enum {dimDomain = FunctionSpaceType::dimDomain};
  public:
    DiscontinuousGalerkinBaseFunctionFactory(GeometryType geo) :
      BaseFunctionFactory<ScalarFunctionSpaceImp>(geo)
    {}

    virtual BaseFunctionType* baseFunction(int i) const 
    {
      // see orthonormalbase_mod.cc for max number of base functions 
      // enum { maxBaseFunctions = 164 };
      enum { maxBaseFunctions = DGNumberOfBaseFunctions<polOrd,dimDomain>::numBaseFunctions };
      switch (GeometryIdentifier::fromGeo(this->geometry())) 
      {
        case GeometryIdentifier::Line:
          return BaseFunctionCreator<GeometryIdentifier::Line, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Triangle:
          return BaseFunctionCreator<GeometryIdentifier::Triangle, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Quadrilateral:
          return BaseFunctionCreator<GeometryIdentifier::Quadrilateral, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Tetrahedron:
          return BaseFunctionCreator<GeometryIdentifier::Tetrahedron, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Pyramid:
          return BaseFunctionCreator<GeometryIdentifier::Pyramid, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Prism:
          return BaseFunctionCreator<GeometryIdentifier::Prism, maxBaseFunctions, 0>::create(i);
        case GeometryIdentifier::Hexahedron:
          return BaseFunctionCreator<GeometryIdentifier::Hexahedron, maxBaseFunctions, 0>::create(i);
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
      return DGBaseFunctionWrapper<FunctionSpaceType,0> 
                :: numBaseFunctions( polOrd , this->geometry().dim() );
    }

  protected:
    template <GeometryIdentifier::IdentifierType GeomId, int maxBaseNum, int baseNum> 
    struct BaseFunctionCreator
    {
      // create base function if i equals baseNum 
      static BaseFunctionType* create(const int i) 
      {
        if ( i == baseNum ) 
        {
          return new DGBaseFunction<FunctionSpaceType, GeomId, polOrd, baseNum> ();
        }
        else 
        {
          return BaseFunctionCreator<GeomId,maxBaseNum,baseNum+1>::create( i );
        }
      }
    };
    
    template <GeometryIdentifier::IdentifierType GeomId, int maxBaseNum> 
    struct BaseFunctionCreator<GeomId,maxBaseNum,maxBaseNum>
    {
      enum { baseNum = maxBaseNum }; 
      // create base function if i equals baseNum 
      static inline BaseFunctionType* create(const int i) 
      {
        if ( i == baseNum ) 
        {
          return new DGBaseFunction<FunctionSpaceType, GeomId, polOrd, baseNum> ();
        }
        else 
        {
          DUNE_THROW(NotImplemented,"BaseFunction " << baseNum << " not found");
          return 0;
        }
      }
    };  
  };



#ifndef COMPILE_LIB
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 0 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > > 
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > > 
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > > 
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 1, 1 >, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 1, 1 > > 
  {
    typedef FunctionSpace< double, double, 1, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 0 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > > 
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > > 
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > > 
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 2, 1 >, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 2, 1 > > 
  {
    typedef FunctionSpace< double, double, 2, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 0 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > > 
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 1 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > > 
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };

  
  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 2 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > > 
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };


  template<>
  class DiscontinuousGalerkinBaseFunctionFactory
    < FunctionSpace< double, double, 3, 1 >, 3 >
  : public BaseFunctionFactory< FunctionSpace< double, double, 3, 1 > > 
  {
    typedef FunctionSpace< double, double, 3, 1 > ScalarFunctionSpace;
    
  public:
    DiscontinuousGalerkinBaseFunctionFactory( GeometryType geometry )
    : BaseFunctionFactory< ScalarFunctionSpace >( geometry )
    {}

    virtual BaseFunctionInterface< ScalarFunctionSpace > *
    baseFunction ( int i ) const;

    virtual int numBaseFunctions () const;
  };
#endif

} // end namespace Dune
#endif

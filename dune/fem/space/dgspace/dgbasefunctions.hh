#ifndef DUNE_DGBASEFUNCTIONS_HH
#define DUNE_DGBASEFUNCTIONS_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>
#include <dune/fem/space/dgspace/orthonormalbase_mod.hh>

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

  //- do not use this elsewhere
  struct DGBaseId 
  {
     enum { simplexId, cubeId, prismId, pyramidId }; 
  };

  //! Wrapper interface for DG base functions
  template <class FunctionSpaceType, int baseNum >
  class DGBaseFunctionWrapper 
  {
  protected:
    DGBaseFunctionWrapper() {}
    virtual ~DGBaseFunctionWrapper() {}
    
    enum { dimDomain = FunctionSpaceType::dimDomain };

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
      double grad[dimDomain];
      OrthonormalBase_1D::grad_line(baseNum,&xi[0],grad);
      return grad[comp];
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
      double grad[dimDomain];
      OrthonormalBase_2D::grad_triangle_2d(baseNum,&xi[0],grad);
      return grad[comp];
    }

    RangeFieldType hess_triangle_2d (const int comp, const DomainType & xi) const
    {
      double h[3];
      OrthonormalBase_2D::hess_triangle_2d(baseNum,&xi[0],h);
      return h[comp];
    }
    
    RangeFieldType grad_quadrilateral_2d (const int comp, const DomainType & xi) const 
    {
      double grad[dimDomain];
      OrthonormalBase_2D::grad_quadrilateral_2d(baseNum, &xi[0], grad); 
      return grad[comp];
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
      double grad[dimDomain];
      OrthonormalBase_3D::grad_tetrahedron_3d(baseNum,&xi[0],grad);
      return grad[comp];
    }
    
    RangeFieldType grad_pyramid_3d (const int comp, const DomainType & xi) const
    {
      double grad[dimDomain];
      OrthonormalBase_3D::grad_pyramid_3d(baseNum,&xi[0], grad); 
      return grad[comp];
    }

    RangeFieldType grad_prism_3d (const int comp, const DomainType & xi) const 
    {
      double grad[dimDomain];
      OrthonormalBase_3D::grad_prism_3d(baseNum,&xi[0], grad ); 
      return grad[comp];
    }

    RangeFieldType grad_hexahedron_3d (const int comp, const DomainType & xi) const
    {
      double grad[dimDomain];
      OrthonormalBase_3D::grad_hexahedron_3d(baseNum,&xi[0], grad ); 
      return grad[comp];
    }

  }; // end class DGBaseFunctionWrapper



  //! Base class for DG base functions
  template< class FunctionSpaceType, int dimDomain, int basicGeo, int polOrd, int baseNum >
  class DGBaseFunction;

  //! Specialisation for lines 
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction< FunctionSpaceType, 1, DGBaseId::simplexId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_line(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_line(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType, baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, line, polOrd>

  //! Specialisation for triangles
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction< FunctionSpaceType, 2, DGBaseId::simplexId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_triangle_2d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_triangle_2d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      int c=diffVariable[0]+diffVariable[1];
      phi = this->hess_triangle_2d(c, x);
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, triangle, polOrd>

  //! Specialisation for quadrilaterals
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, 2, DGBaseId::cubeId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_quadrilateral_2d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_quadrilateral_2d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
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
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction< FunctionSpaceType, 3, DGBaseId::simplexId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_tetrahedron_3d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_tetrahedron_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
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
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction< FunctionSpaceType, 3, DGBaseId::pyramidId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType,baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_pyramid_3d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_pyramid_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, pyramid, polOrd>

  //! Specialisation for prisms
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, 3, DGBaseId::prismId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const {
      phi = this->eval_prism_3d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_prism_3d(diffVariable[0], x );
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      DUNE_THROW(NotImplemented,"hessian not implemented!");
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, prism, polOrd>

  //! Specialisation for hexahedrons
  template< class FunctionSpaceType, int polOrd, int baseNum >
  class DGBaseFunction<FunctionSpaceType, 3, DGBaseId::cubeId, polOrd, baseNum >
  : public BaseFunctionInterface< FunctionSpaceType >,
    private DGBaseFunctionWrapper< FunctionSpaceType, baseNum >
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

    virtual void evaluate(const FieldVector<int, 0>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->eval_hexahedron_3d(x);
    }
    
    virtual void evaluate(const FieldVector<int, 1>& diffVariable,
                          const DomainType& x, RangeType& phi) const 
    {
      phi = this->grad_hexahedron_3d(diffVariable[0], x);
    }

    virtual void evaluate(const FieldVector<int, 2>&diffVariable,
                          const DomainType& x, RangeType& phi) const {
      assert(false); // Not implemented
    }

    static int numBaseFunctions() {
      return DGBaseFunctionWrapper<FunctionSpaceType,baseNum>::
        numBaseFunctions(static_cast<int>(polOrd));
    }

  }; // end class DGBaseFunction<FunctionSpaceType, hexahedron, polOrd>
  
  template< class ScalarFunctionSpaceImp, int polOrd >
  class DiscontinuousGalerkinBaseFunctionFactory
  : public BaseFunctionFactory< ScalarFunctionSpaceImp >
  {
    typedef BaseFunctionFactory< ScalarFunctionSpaceImp > BaseType;
  public:
    // Add compile time checker: only scalar functions allowed

    typedef ScalarFunctionSpaceImp FunctionSpaceType;
    typedef BaseFunctionInterface<FunctionSpaceType> BaseFunctionType;

    static const int dimDomain = FunctionSpaceType::dimDomain;

  public:
    DiscontinuousGalerkinBaseFunctionFactory ( const GeometryType &type )
    : BaseType( type )
    {}

    using BaseType::geometry;

    virtual BaseFunctionType* baseFunction ( const int i ) const
    {
      const int maxBaseFunctions = DGNumberOfBaseFunctions< polOrd, dimDomain >::numBaseFunctions;

      const GeometryType geo = geometry();
      assert( (int) geo.dim() == dimDomain );

      if( dimDomain == 1 )
        return BaseFunctionCreator< 1, DGBaseId::simplexId, maxBaseFunctions, 0 >::create( i );
      else if( dimDomain == 2 )
      {
        if( geo.isSimplex() ) 
        {
          return BaseFunctionCreator< 2, DGBaseId::simplexId, maxBaseFunctions, 0 >::create( i );
        } 
        else 
        {
          assert( geo.isCube() );
          return BaseFunctionCreator< 2, DGBaseId::cubeId, maxBaseFunctions, 0 >::create( i );
        }
      }
      else if( dimDomain == 3 )
      {
        if( geo.isSimplex() )
        {
          return BaseFunctionCreator< 3, DGBaseId::simplexId, maxBaseFunctions, 0 >::create( i );
        }
        else if( geo.isCube() )
        {
          return BaseFunctionCreator< 3, DGBaseId::cubeId, maxBaseFunctions, 0 >::create( i );
        }
        else if( geo.isPrism() ) 
        {
          return BaseFunctionCreator< 3, DGBaseId::prismId, maxBaseFunctions, 0 >::create( i );
        }
        else if( geo.isPyramid() )
        {
          return BaseFunctionCreator< 3, DGBaseId::pyramidId, maxBaseFunctions, 0 >::create( i );
        }
        else 
          DUNE_THROW( InvalidStateException, "Invalid geometry type." );
      }
      else
        DUNE_THROW( NotImplemented, "DGBaseFunction only implemented for dim=1,2,3." );
    }

    virtual int numBaseFunctions() const 
    {
      // call numBaseFunctions from Wrapper class depending on 
      // geometry's dimension 
      return DGBaseFunctionWrapper<FunctionSpaceType,0> 
                :: numBaseFunctions( polOrd , this->geometry().dim() );
    }

  protected:
    template< int dimD, int basicGeo, int maxBaseNum, int baseNum >
    struct BaseFunctionCreator
    {
      // create base function if i equals baseNum 
      static BaseFunctionType* create ( const int i )
      {
        if( i == baseNum ) 
          return new DGBaseFunction< FunctionSpaceType, dimD, basicGeo, polOrd, baseNum >();
        else 
          return BaseFunctionCreator< dimD, basicGeo, maxBaseNum, baseNum+1 >::create( i );
      }
    };
    
    template< int dimD, int basicGeo, int maxBaseNum >
    struct BaseFunctionCreator< dimD, basicGeo, maxBaseNum, maxBaseNum >
    {
      static const int baseNum = maxBaseNum;

      // create base function if i equals baseNum 
      static inline BaseFunctionType* create( const int i )
      {
        if( i == baseNum )
          return new DGBaseFunction< FunctionSpaceType, dimD, basicGeo, polOrd, baseNum >();
        else 
          DUNE_THROW( NotImplemented, "BaseFunction " << baseNum << " not found." );
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
#endif // #ifndef COMPILE_LIB

} // end namespace Dune

#endif // #ifndef DUNE_DGBASEFUNCTIONS_HH

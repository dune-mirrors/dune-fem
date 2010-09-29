#ifndef DUNE_FEM_VECTORIALBASEFUNCTIONSET_HH
#define DUNE_FEM_VECTORIALBASEFUNCTIONSET_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>


#ifdef USE_BASEFUNCTIONSET_CODEGEN 
#define USE_BASEFUNCTIONSET_OPTIMIZED
#endif

#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
#include <dune/fem/space/basefunctions/codegen.hh>
#endif

#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
#include <dune/fem/space/basefunctions/evaluatecaller.hh>
#endif

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------
  template< class FunctionSpace, template< class > class Storage >
  class VectorialBaseFunctionSet;



  // VectorialBaseFunctionSetTraits
  // ------------------------------

  template< class FunctionSpace, template< class > class Storage >
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpace FunctionSpaceType;
    typedef VectorialBaseFunctionSet< FunctionSpace, Storage > BaseFunctionSetType;
  };



  // VectorialBaseFunctionSet
  // ------------------------

  /** \class VectorialBaseFunctionSet
   *  \ingroup BaseFunction
   *  \brief special base function set taking advantage of the vectorial
   *         structure of base functions
   *
   *  This base function can be used in conjunction with scalar basefunctions
   *  \f$ \phi_i \f$ which are extended to vectorial base functions like 
   *  \f$ \Phi_j = \phi_i e_k \f$, where \f$ e_k = [ \delta_{ik} ]_i \f$.
   *
   *  \tparam  FunctionSpace  analytical function space
   *  \tparam  Storage        base function storage (either CachingStorage or SimpleStorage)
   */
  template< class FunctionSpace, template< class > class Storage >
  class VectorialBaseFunctionSet
  : public BaseFunctionSetDefault< VectorialBaseFunctionSetTraits< FunctionSpace, Storage > >
  {
    typedef VectorialBaseFunctionSet< FunctionSpace, Storage > ThisType;
    typedef BaseFunctionSetDefault< VectorialBaseFunctionSetTraits< FunctionSpace, Storage > > BaseType;

  public:
    typedef VectorialBaseFunctionSetTraits< FunctionSpace, Storage > Traits;

    typedef FunctionSpace FunctionSpaceType;

    static const int dimDomain = FunctionSpaceType::dimDomain;
    static const int dimRange  = FunctionSpaceType::dimRange;

    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
    
    typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType::JacobianRangeType ScalarJacobianRangeType;

    typedef Storage< ScalarFunctionSpaceType > StorageType;

    typedef RangeFieldType DofType;
    
  private:  
    typedef BaseFunctionFactory< ScalarFunctionSpaceType > FactoryType;
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;

  public:
    //! Constructor
    explicit VectorialBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory ),
      util_( dimRange )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      // add my dimrange 
      Fem::CodegenInfo::instance().addDimRange( dimRange );
#endif
    }

    const StorageType &storage () const { return storage_; }

    using BaseType::evaluate;
    using BaseType::evaluateSingle;
    using BaseType::evaluateGradientSingle;
    using BaseType::jacobian;

  private:
    VectorialBaseFunctionSet ( const ThisType & );
    const ThisType &operator= ( const ThisType & );
    
  public:
    GeometryType geometryType () const { return storage_.geometryType(); }
    
    int numBaseFunctions () const { return dimRange*numDifferentBaseFunctions(); }
    int numDifferentBaseFunctions () const { return storage_.numBaseFunctions(); }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int,const FieldVector<deriType,diffOrd>&,const PointType&,RangeType&) const */
    template< int diffOrd, class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< deriType, diffOrd > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const
    {
      ScalarRangeType tmp;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      evaluateScalar( scalarBaseFunction, diffVariable, x, tmp );

      phi = 0;
      phi[ util_.component( baseFunction ) ] = tmp[ 0 ];
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */
    template< class PointType >
    void jacobian ( const int baseFunction,
                    const PointType &x, 
                    JacobianRangeType &phi ) const
    {
      ScalarJacobianRangeType tmp;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      jacobianScalar( scalarBaseFunction, x, tmp );

      phi = 0;
      phi[ util_.component( baseFunction )] = tmp[ 0 ];
    }

    template< class PointType >
    void jacobianScalar ( const int baseFunction,
                          const PointType &x,
                          ScalarJacobianRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.jacobian( baseFunction, x, phi );
    }
    
    template< int diffOrder, class PointType, class LocalDofVectorType >
    void evaluateAll ( const FieldVector< int, diffOrder > &diffVariable,
                       const PointType &x,
                       const LocalDofVectorType& dofs, 
                       RangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0, iR = 0; i < numScalarBase ; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, diffVariable, x, phi ); 
        for( int r = 0; r < dimRange; ++r , ++iR )
        { 
          ret[ r ] += phi[ 0 ] * dofs[ iR ];
        }
      }
    }

    template< class PointType, class LocalDofVectorType >
    void evaluateAll ( const PointType &x,
                       const LocalDofVectorType& dofs, 
                       RangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions();
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi ); 
        for( int r = 0; r < dimRange; ++r, ++iR )
        { 
          ret[ r ] += phi[ 0 ] * dofs[ iR ];
        }
      }
    }

    template< class PointType, class RangeVectorType >
    void evaluateAll ( const PointType &x, RangeVectorType &ret ) const
    {
      const int numScalarBase = numDifferentBaseFunctions();
      ret.resize( dimRange*numScalarBase );
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi ); 
        for( int r = 0; r < dimRange; ++r, ++iR )
        { 
          ret[ iR ] = 0;
          ret[ iR ][ r ] = phi[ 0 ];
        }
      }
    }

    template< class PointType, 
              class GeometryJacobianInverseType,
              class LocalDofVectorType, 
              class GlobalJacobianRangeType>
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       const LocalDofVectorType& dofs, 
                       GlobalJacobianRangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions();
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarJacobianRangeType gradPhiRef;
        // get type of scalar global jacobian range 
        typename GlobalJacobianRangeType::row_type gradPhi;

        jacobianScalar( i, x, gradPhiRef );

        gjit.mv( gradPhiRef[ 0 ], gradPhi );

        for( int r = 0; r < dimRange; ++r, ++iR )
          ret[ r ].axpy( dofs[ iR ], gradPhi );
      }
    }
   
    template< class PointType, class GeometryJacobianInverseType,
              class GlobalJacobianRangeVectorType >
    void jacobianAll ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       GlobalJacobianRangeVectorType &ret ) const
    {
      const int numScalarBase = numDifferentBaseFunctions();
      ret.resize( dimRange*numScalarBase );
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarJacobianRangeType gradPhiRef;
        jacobianScalar( i, x, gradPhiRef );

        const int iR0 = iR++;
        ret[ iR0 ] = 0;
        gjit.mv( gradPhiRef[ 0 ], ret[ iR0 ][ 0 ] );

        for( int r = 1; r < dimRange; ++r, ++iR )
        {
          ret[ iR ] = 0;
          ret[ iR ][ r ] = ret[ iR0 ][ 0 ];
        }
      }
    }

    template< class PointType, class LocalDofVectorType >
    void axpy ( const PointType &x,
                const RangeType &rangeFactor,
                LocalDofVectorType& dofs ) const
    {
      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0, iR = 0; i < numScalarBase ; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi ); 
        for( int r = 0; r < dimRange; ++r , ++iR )
        { 
          dofs[ iR ] += phi[ 0 ] * rangeFactor[ r ];
        }
      }
    }
    
    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit, 
                const GlobalJacobianRangeType &jacFactor,
                LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacFactor[ r ], jacFactorInv[ r ] );

      const int numScalarBase = numDifferentBaseFunctions();
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarJacobianRangeType grad;
        jacobianScalar( i, x, grad );
        for( int r = 0; r < dimRange; ++r, ++iR )
          dofs[ iR ] += grad[ 0 ] * jacFactorInv[ r ];
      }
    }
 

    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    void axpy ( const PointType &x,
                const GeometryJacobianInverseType& gjit, 
                const RangeType& rangeFactor,
                const GlobalJacobianRangeType &jacFactor,
                LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      for( int r = 0; r < dimRange; ++r )
        gjit.mtv( jacFactor[ r ], jacFactorInv[ r ] );

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0, iR = 0; i < numScalarBase; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi );
        ScalarJacobianRangeType grad;
        jacobianScalar( i, x, grad );
        for( int r = 0; r < dimRange; ++r, ++iR  )
          dofs[ iR ] += (phi[ 0 ] * rangeFactor[ r ]) + (grad[ 0 ] * jacFactorInv[ r ]);
      }
    }
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    RangeFieldType evaluateSingle ( const int baseFunction,
                                    const PointType &x,
                                    const RangeType &psi ) const
    {
      ScalarRangeType phi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );

      evaluateScalar( scalarBaseFunction, x, phi );
      return psi[ util_.component( baseFunction ) ] * phi[ 0 ];
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateGradientSingle(const int baseFunction,const EntityType &entity,const PointType &x,const JacobianRangeType &psi) const */
    template< class EntityType, class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    RangeFieldType evaluateGradientSingle( const int baseFunction,
                                           const EntityType &entity,
                                           const PointType &x,
                                           const JacobianRangeType &psi ) const
    {
      typedef typename EntityType :: Geometry GeometryType;
      typedef FieldMatrix< typename GeometryType :: ctype,
                           GeometryType :: mydimension,
                           GeometryType :: mydimension >
        GeometryJacobianType;

      const GeometryType &geometry = entity.geometry();
      const GeometryJacobianType &jacobianInverseTransposed
        = geometry.jacobianInverseTransposed( coordinate( x ) );
      
      ScalarJacobianRangeType gradPhi;
      const int scalarBaseFunction = util_.containedDof( baseFunction );
      jacobianScalar( scalarBaseFunction, x, gradPhi );
    
      DomainType gradScaled( 0 );
      jacobianInverseTransposed.umv( gradPhi[ 0 ], gradScaled );
      return gradScaled * psi[ util_.component( baseFunction ) ];
    }

    /////////////////////////////////////////////////////////////////////////
    //
    //  evaluate and store results in a vector 
    //
    /////////////////////////////////////////////////////////////////////////
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
    template <class BaseFunctionSet, class Geometry, int dimRange, int numRows, int numCols>
    struct EvaluateRanges 
    {
      template< class QuadratureType, 
                class RangeVectorType,
                class LocalDofVectorType,
                class RangeFactorType>
      static void eval( const QuadratureType& quad,
                        const RangeVectorType& rangeStorage,
                        const LocalDofVectorType& dofs,
                        RangeFactorType &rangeFactors)
      {
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateRanges< "
                  << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
      }
    };

    template <class BaseFunctionSet, int dimRange, int numRows, int numCols>
    struct EvaluateRanges<BaseFunctionSet, Fem :: EmptyGeometry, dimRange, numRows, numCols > 
    {
      template< class QuadratureType, 
                class RangeVectorType,
                class LocalDofVectorType,
                class RangeFactorType>
      static void eval( const QuadratureType& quad,
                        const RangeVectorType& rangeStorage,
                        const LocalDofVectorType& dofs,
                        RangeFactorType &rangeFactors)
      {
#ifndef USE_BASEFUNCTIONSET_CODEGEN
        BaseFunctionSet::evaluateRanges( quad, rangeStorage, dofs, rangeFactors, numRows, numCols );
#else 
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateRanges< "
                  << "EmptyGeo, " << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
#endif
      }
    };

#endif

    template< class QuadratureType, 
              class LocalDofVectorType,
              class RangeFactorType>
    void evaluateRanges ( const QuadratureType& quad,
                          const LocalDofVectorType& dofs,
                          RangeFactorType &rangeVector) const 
    {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      typedef Fem :: EvaluateCallerInterfaceTraits< 
          QuadratureType, RangeFactorType, LocalDofVectorType > Traits;
      typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

      // get base function evaluate caller (calls evaluateRanges) 
      const BaseEvaluationType& baseEval = 
        BaseEvaluationType::storage( *this, storage_.getRangeStorage( quad ), quad );

      baseEval.evaluateRanges( quad, dofs, rangeVector );
#else 
      evaluateRanges( quad, storage_.getRangeStorage( quad ), 
          dofs, rangeVector, quad.nop(), numDifferentBaseFunctions() );
#endif
    }

    template< class QuadratureType, 
              class LocalDofVectorType,
              class RangeVectorType,
              class RangeFactorType>
    static void 
    evaluateRanges ( const QuadratureType& quad,
                     const RangeVectorType& rangeStorage,
                     const LocalDofVectorType& dofs,
                     RangeFactorType &rangeVector,
                     const unsigned int numRows,
                     const unsigned int numCols )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      static std::map< const size_t, size_t > storedRows; 
      typename std::map< const size_t, size_t > :: iterator it = storedRows.find( numRows );
      if( it != storedRows.end() ) 
      {
        Fem::CodegenInfo::instance().notify( (*it).second );
      }
      else 
      {
        storedRows[ numRows ] = 
          Fem::CodegenInfo::instance().addEntry( "evalranges", 
              ( storedRows.size() == 1 ), Fem :: CodeGeneratorType :: evaluateCodegen, dimDomain, dimRange, numRows, numCols );
        std::cout << "Generate code evalranges for (" << dimRange << "," << numRows << "," << numCols << ")" << std::endl;
      }
#endif

      assert( (int) numCols * dimRange == dofs.numDofs() );

      assert( rangeStorage.size() >= (int) numRows );
      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = quad.cachingPoint( row ); 

        assert( rangeStorage.size() > (int) baseRow );
        assert( rangeStorage[ baseRow ].size() >= (int) numCols );

        RangeType& result = rangeVector[ row ]; 
        result = 0;

        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          const ScalarRangeType& phi = rangeStorage[ baseRow ][ col ];
          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            result[ r ] +=  dofs[ colR ] * phi[ 0 ];
          }
        }
      }
    }

    ////////////////////////////////////////////////////////////////////
    //
    //  --evaluateJacobians 
    //
    ////////////////////////////////////////////////////////////////////
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
    template <class BaseFunctionSet, class Geometry, 
              int dimRange, int numRows, int numCols>
    struct EvaluateJacobians 
    {
      template< class QuadratureType, 
                class JacobianRangeVectorType,
                class JacobianRangeFactorType,
                class LocalDofVectorType >
      static void eval( const QuadratureType& quad,
                        const Geometry& geometry,
                        const JacobianRangeVectorType& jacobianStorage,
                        const LocalDofVectorType& dofs,
                        JacobianRangeFactorType &jacFactors)
      {
#ifndef USE_BASEFUNCTIONSET_CODEGEN
        BaseFunctionSet :: 
          evaluateJacobians( quad, geometry, jacobianStorage, dofs, jacFactors, 
                             jacFactors[ 0 ], numRows, numCols );
#else 
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateJacobians< "
                  << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
#endif
      }
    };

    template <class BaseFunctionSet,
              int dimRange, int numRows, int numCols>
    struct EvaluateJacobians< BaseFunctionSet, Fem :: EmptyGeometry, dimRange, numRows, numCols > 
    {
      template< class QuadratureType, 
                class JacobianRangeVectorType,
                class JacobianRangeFactorType,
                class LocalDofVectorType >
      static void eval( const QuadratureType&,
                        const Fem :: EmptyGeometry&,
                        const JacobianRangeVectorType&,
                        const LocalDofVectorType&,
                        const JacobianRangeFactorType& )
      {
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::evaluateJacobians< "
                  << "EmptyGeo, " << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
      }
    };
#endif // endif USE_BASEFUNCTIONSET_OPTIMIZED 

    template< class QuadratureType, 
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    void evaluateJacobians ( const QuadratureType& quad,
                             const Geometry& geometry,
                             const LocalDofVectorType& dofs,
                             JacobianRangeVectorType &jacVector ) const 
    {
      assert( jacVector.size() > 0 );
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      typedef Fem :: EvaluateCallerInterfaceTraits< QuadratureType, 
              JacobianRangeVectorType, LocalDofVectorType, Geometry >  Traits;
      typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

      // get base function evaluate caller (calls axpyRanges) 
      const BaseEvaluationType& baseEval = 
        BaseEvaluationType::storage( *this, storage_.getJacobianStorage( quad ), quad );

      // call appropriate axpyJacobian method 
      baseEval.evaluateJacobians( quad, geometry, dofs, jacVector );
#else 
      evaluateJacobians( quad, geometry, 
                         storage_.getJacobianStorage( quad ),
                         dofs, jacVector, jacVector[ 0 ], 
                         quad.nop(), numDifferentBaseFunctions() );
#endif
    }

    template< class QuadratureType, 
              class Geometry,
              class JacobianRangeVectorType,
              class LocalDofVectorType,
              class JacobianVectorType,
              class GlobalJacobianRangeType >
    static void 
    evaluateJacobians ( const QuadratureType& quad,
                        const Geometry& geometry,
                        const JacobianRangeVectorType& jacobianStorage,
                        const LocalDofVectorType& dofs,
                        JacobianVectorType &jacVector,
                        const GlobalJacobianRangeType&,
                        const unsigned int numRows,
                        const unsigned int numCols )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      static std::map< const size_t , size_t > storedRows; 
      typename std::map< const size_t , size_t > :: iterator it = storedRows.find( numRows );
      if( it != storedRows.end() ) 
      {
        Fem::CodegenInfo::instance().notify( (*it).second );
      }
      else 
      {
        storedRows[ numRows ] = 
          Fem::CodegenInfo::instance().addEntry( "evaljacobians", 
              ( storedRows.size() == 1 ), Fem :: CodeGeneratorType :: evaluateJacobiansCodegen, dimDomain, dimRange, numRows, numCols );
        std::cout << "Generate code evaljacobians for (" << numRows << "," << numCols << ")" << std::endl;
      }
#endif

      assert( (int) numCols * dimRange == dofs.numDofs() );
      assert( jacobianStorage.size() >= (int)numRows );

      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = quad.cachingPoint( row ); 

        // if geometry has non-affine mapping we need to update jacobian inverse
        typedef typename Geometry::Jacobian GeometryJacobianType;
        const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );

        assert( jacobianStorage.size() > (int) baseRow );
        assert( jacobianStorage[ baseRow ].size() >= (int) numCols );

        GlobalJacobianRangeType& result = jacVector[ row ]; 
        result = 0;

        // get type of scalar global jacobian 
        // (which is one row of the GlobalJacobianRangeType)
        typedef typename GlobalJacobianRangeType :: row_type JacobianRangeType;
        JacobianRangeType gradPhi;

        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          gjit.mv( jacobianStorage[ baseRow ][ col ][ 0 ], gradPhi );

          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            result[ r ].axpy( dofs[ colR ], gradPhi );
          }
        }
      }
    }
   
    /////////////////////////////////////////////////////////////
    //
    //  --axpyRanges -- add a vector of ranges to the dof vector  
    //
    /////////////////////////////////////////////////////////////
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
    template <class BaseFunctionSet, class Geometry, 
              int dimRange, int numRows, int numCols>
    struct AxpyRanges 
    {
      template< class QuadratureType, 
                class RangeVectorType,
                class RangeFactorType,
                class LocalDofVectorType >
      static void axpy( const QuadratureType& quad,
                        const RangeVectorType& rangeStorage,
                        const RangeFactorType &rangeFactors,
                        LocalDofVectorType& dofs)
      {
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyRanges <"
                  << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
      }
    };

    template <class BaseFunctionSet,
              int dimRange, int numRows, int numCols>
    struct AxpyRanges<BaseFunctionSet, Fem :: EmptyGeometry, dimRange, numRows, numCols>  
    {
      template< class QuadratureType, 
                class RangeVectorType,
                class RangeFactorType,
                class LocalDofVectorType >
      static void axpy( const QuadratureType& quad,
                        const RangeVectorType& rangeStorage,
                        const RangeFactorType &rangeFactors,
                        LocalDofVectorType& dofs)
      {
#ifndef USE_BASEFUNCTIONSET_CODEGEN
        BaseFunctionSet::axpyRanges( quad, rangeStorage, rangeFactors, dofs, numRows, numCols );
#else 
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyRanges <"
                  << dimRange << " , " << numRows << " , " << numCols << " >!" << std::endl;
        abort();
#endif
      }
    };

#endif // endif USE_BASEFUNCTIONSET_OPTIMIZED 

    template< class QuadratureType, 
              class RangeFactorType,
              class LocalDofVectorType >
    void axpyRanges ( const QuadratureType& quad,
                      const RangeFactorType &rangeFactors,
                      LocalDofVectorType& dofs ) const
    {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      typedef Fem :: EvaluateCallerInterfaceTraits< 
          QuadratureType, RangeFactorType, LocalDofVectorType > Traits;
      typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

      // get base function evaluate caller (calls axpyRanges) 
      const BaseEvaluationType& baseEval = 
        BaseEvaluationType::storage( *this, storage_.getRangeStorage( quad ), quad );

      // call appropriate axpyRanges method 
      baseEval.axpyRanges( quad, rangeFactors, dofs );
#else 
      axpyRanges( quad, storage_.getRangeStorage( quad ), 
          rangeFactors, dofs, quad.nop(), numDifferentBaseFunctions() );
#endif
    }

    template< class QuadratureType, 
              class RangeVectorType,
              class RangeFactorType,
              class LocalDofVectorType >
    static void axpyRanges ( const QuadratureType& quad,
                             const RangeVectorType& rangeStorage,
                             const RangeFactorType &rangeFactors,
                             LocalDofVectorType& dofs,
                             const unsigned int numRows, 
                             const unsigned int numCols )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      static std::map< const size_t , size_t > storedRows; 
      typename std::map< const size_t , size_t > :: iterator it = storedRows.find( numRows );
      if( it != storedRows.end() ) 
      {
        Fem::CodegenInfo::instance().notify( (*it).second );
      }
      else 
      {
        storedRows[ numRows ] = 
          Fem::CodegenInfo::instance().addEntry( "axpyranges", 
              ( storedRows.size() == 1 ), Fem :: CodeGeneratorType :: axpyCodegen, dimDomain, dimRange, numRows, numCols );
        std::cout << "Generate code axpyranges for (" << numRows << "," << numCols << ")" << std::endl;
      }
#endif

      assert( numRows == quad.nop() );
      assert( (int) numCols * dimRange == dofs.numDofs() );
      assert( rangeStorage.size() >= (int) numRows );

      // this way the codes seems faster 
      for( size_t col = 0, colR = 0; col < numCols; ++col ) 
      {
        for( int r = 0; r < dimRange; ++r , ++colR ) 
        {
          for( size_t row = 0; row < numRows ; ++row )
          {
            dofs[ colR ] += rangeStorage[ quad.cachingPoint( row ) ][ col ][ 0 ] * rangeFactors[ row ][ r ];
          }
        }
      }
    }

    ///////////////////////////////////////////////////////////
    //  applyAxpy Jacobian 
    ///////////////////////////////////////////////////////////
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
    template <class BaseFunctionSet, class Geometry, 
              int dimRange, int numRows, int numCols>
    struct AxpyJacobians 
    {
      template< class QuadratureType, 
                class JacobianRangeVectorType,
                class JacobianRangeFactorType,
                class LocalDofVectorType >
      static void axpy( const QuadratureType& quad,
                        const Geometry& geometry,
                        const JacobianRangeVectorType& jacobianStorage,
                        const JacobianRangeFactorType &jacFactors,
                        LocalDofVectorType& dofs)
      {
        BaseFunctionSet :: 
          axpyJacobians( quad, geometry, jacobianStorage, jacFactors, dofs, 
                         jacFactors[ 0 ], numRows, numCols );
      }
    };

    template <class BaseFunctionSet,
              int dimRange, int numRows, int numCols>
    struct AxpyJacobians< BaseFunctionSet, Fem :: EmptyGeometry, dimRange, numRows, numCols > 
    {
      template< class QuadratureType, 
                class JacobianRangeVectorType,
                class JacobianRangeFactorType,
                class LocalDofVectorType >
      static void axpy( const QuadratureType&,
                        const Fem :: EmptyGeometry&,
                        const JacobianRangeVectorType&,
                        const JacobianRangeFactorType &,
                        LocalDofVectorType&)
      {
        std::cerr << "ERROR: wrong code generated for VectorialBaseFunctionSet::axpyJacobians" << std::endl;
        abort();
      }
    };

#ifdef USE_BASEFUNCTIONSET_CODEGEN
    #include <autogeneratedcode.hh>
#endif
#endif // endif USE_BASEFUNCTIONSET_OPTIMIZED 

    template< class QuadratureType, 
              class Geometry,
              class JacobianVectorType,
              class LocalDofVectorType >
    void axpyJacobians ( const QuadratureType& quad,
                         const Geometry& geometry,
                         const JacobianVectorType &jacVector,
                         LocalDofVectorType& dofs) const 
    {
#ifdef USE_BASEFUNCTIONSET_OPTIMIZED
      typedef Fem :: EvaluateCallerInterfaceTraits< QuadratureType, 
              JacobianVectorType, LocalDofVectorType, Geometry >  Traits;
      typedef Fem :: EvaluateCallerInterface< Traits > BaseEvaluationType;

      // get base function evaluate caller (calls axpyRanges) 
      const BaseEvaluationType& baseEval = 
        BaseEvaluationType::storage( *this, storage_.getJacobianStorage( quad ), quad );

      // call appropriate axpyJacobian method 
      baseEval.axpyJacobians( quad, geometry, jacVector, dofs );
#else 
      assert( jacVector.size() > 0 );
      axpyJacobians( quad, geometry,
                     storage_.getJacobianStorage( quad ),
                     jacVector, dofs, jacVector[ 0 ], 
                     quad.nop(), numDifferentBaseFunctions() );
#endif
    }

    template< class QuadratureType, 
              class Geometry,
              class JacobianRangeVectorType,
              class JacobianVectorType,
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    static void axpyJacobians ( const QuadratureType& quad,
                                const Geometry& geometry,
                                const JacobianRangeVectorType& jacobianStorage,
                                const JacobianVectorType &jacVector,
                                LocalDofVectorType& dofs,
                                const GlobalJacobianRangeType&,
                                const unsigned int numRows, 
                                const unsigned int numCols )
    {
#ifdef BASEFUNCTIONSET_CODEGEN_GENERATE
      static std::map< const size_t , size_t > storedRows; 
      typename std::map< const size_t , size_t > :: iterator it = storedRows.find( numRows );
      if( it != storedRows.end() ) 
      {
        Fem::CodegenInfo::instance().notify( (*it).second );
      }
      else 
      {
        storedRows[ numRows ] = 
          Fem::CodegenInfo::instance().addEntry( "axpyjacobians", 
              ( storedRows.size() == 1 ), Fem :: CodeGeneratorType :: axpyJacobianCodegen, dimDomain, dimRange, numRows, numCols );
        std::cout << "Generate code axpyjacobians for (" << numRows << "," << numCols << ")" << std::endl;
      }
#endif
      assert( (int ) numCols * dimRange == dofs.numDofs() );
      assert( jacobianStorage.size() >= (int) numRows );

      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = quad.cachingPoint( row ); 
        assert( jacobianStorage.size() > (int)baseRow );
        assert( jacobianStorage[ baseRow ].size() >= (int)numCols );

        typedef typename Geometry::Jacobian GeometryJacobianType;
        const GeometryJacobianType& gjit = geometry.jacobianInverseTransposed( quad.point( row ) );

        JacobianRangeType jacFactorInv;

        // multiply jacobian factor with geometry inverse 
        for( int r = 0; r < dimRange; ++r )
          gjit.mtv( jacVector[ row ][ r ], jacFactorInv[ r ] );

        for( size_t col = 0, colR = 0; col < numCols; ++col ) 
        {
          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            dofs[ colR ] += jacobianStorage[ baseRow ][ col ][ 0 ] * jacFactorInv[ r ];
          }
        }
      }
    }

  protected:
    template< int diffOrd, class PointType >
    void evaluateScalar ( const int baseFunction,
                          const FieldVector< int, diffOrd > &diffVariable,
                          const PointType &x,
                          ScalarRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }
    
    template< class PointType >
    void evaluateScalar ( const int baseFunction,
                          const PointType &x,
                          ScalarRangeType &phi ) const
    {
      FieldVector< int, 0 > diffVar;
      evaluateScalar( baseFunction, diffVar, x, phi );
    }

  protected:
    StorageType storage_;
    PointBasedDofConversionUtility< dimRange > util_;
  };

} // end namespace Dune

#endif // #ifndef DUNE_FEM_VECTORIALBASEFUNCTIONSET_HH

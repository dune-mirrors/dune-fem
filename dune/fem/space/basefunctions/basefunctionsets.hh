#ifndef DUNE_BASEFUNCTIONSETS_HH
#define DUNE_BASEFUNCTIONSETS_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/basefunctions/basefunctioninterface.hh>
#include <dune/fem/space/basefunctions/basefunctionfactory.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

namespace Dune
{
  
  /** \addtogroup BaseFunction
   *  \{
   */
  
  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet;



  //! Traits class for standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StandardBaseFunctionSetTraits
  {
    //! Export function space type
    typedef FunctionSpaceImp FunctionSpaceType;
    //! Type of the base function storage policy
    typedef StorageImp<FunctionSpaceType> StorageType;
    //! Exact type of the base function
    typedef StandardBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    //! Factory type for the corresponding base functions (polymorphic)
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;
  };



  /** \class StandardBaseFunctionSet
   *  \brief standard base function set
   */
  template< class FunctionSpaceImp, template< class > class StorageImp >
  class StandardBaseFunctionSet
  : public BaseFunctionSetDefault
    < StandardBaseFunctionSetTraits< FunctionSpaceImp, StorageImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

    typedef StandardBaseFunctionSetTraits< FunctionSpaceType, StorageImp > Traits;

  private:
    typedef StandardBaseFunctionSet< FunctionSpaceType, StorageImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;
    
  protected:  
    typedef typename Traits::StorageType StorageType; 
  public:
    enum { dimRange = FunctionSpaceType :: dimRange };  
      
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: RangeFieldType DofType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    typedef typename Traits :: FactoryType FactoryType;

  public:
    // use evaluate of default implementation 
    using BaseType :: evaluate;
    using BaseType :: jacobian;

  public:
    //! Constructor
    inline explicit StandardBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory )
    {
    }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions () const
    {
      return storage_.numBaseFunctions();
    }
    
    /** \copydoc Dune::BaseFunctionSetInterface::geometryType */
    inline GeometryType geometryType () const
    {
      return storage_.geometryType();
    }
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */ 
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const
    {
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }

    /** \copydoc Dune::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */ 
    template< class PointType >
    inline void jacobian ( const int baseFunction,
                           const PointType &x,
                           JacobianRangeType &phi ) const
    {
      storage_.jacobian( baseFunction, x, phi );
    }

  private:
    StandardBaseFunctionSet( const StandardBaseFunctionSet& );

  protected:
    StorageType storage_;
  };


  //- VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VectorialBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };



  //! \brief Special base function implementation that takes advantage
  //! of the vectorial structure of the base functions.
  //! This base function can be used in conjunction with scalar basefunctions
  //! \f$ \phi_i \f$ which are extended to vectorial base functions like 
  //! \f$ \Phi_j = \phi_i e_k \f$, where \f$ e_k = [ \delta_{ik} ]_i \f$.
  template< class FunctionSpaceImp, template< class > class StorageImp >
  class VectorialBaseFunctionSet
  : public BaseFunctionSetDefault
    < VectorialBaseFunctionSetTraits< FunctionSpaceImp, StorageImp > >
  {
  public:
    typedef FunctionSpaceImp FunctionSpaceType;

    typedef VectorialBaseFunctionSetTraits< FunctionSpaceType, StorageImp > Traits;
    
  private:
    typedef VectorialBaseFunctionSet< FunctionSpaceType, StorageImp > ThisType;
    typedef BaseFunctionSetDefault< Traits > BaseType;
    
  public:
    typedef typename FunctionSpaceType :: ScalarFunctionSpaceType
      ScalarFunctionSpaceType;
    
  private:  
    typedef BaseFunctionFactory< ScalarFunctionSpaceType > FactoryType;
    typedef typename FactoryType :: BaseFunctionType BaseFunctionType;
    typedef StorageImp< ScalarFunctionSpaceType > StorageType;

  public:
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

    enum
    {
      dimDomain = FunctionSpaceType :: dimDomain,
      dimRange = FunctionSpaceType :: dimRange
    };
    
    typedef RangeFieldType DofType;
    
    typedef typename ScalarFunctionSpaceType :: RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType :: JacobianRangeType
      ScalarJacobianRangeType;
    
  protected:
    StorageType storage_;
    PointBasedDofConversionUtility< dimRange > util_;

  public:
    using BaseType :: evaluate;
    using BaseType :: evaluateSingle;
    using BaseType :: evaluateGradientSingle;
    using BaseType :: jacobian;
       
  public:
    //! Constructor
    inline explicit VectorialBaseFunctionSet ( const FactoryType &factory )
    : storage_( factory ),
      util_( FunctionSpaceType::dimRange )
    {
    }

  private:
    VectorialBaseFunctionSet ( const ThisType & );
    
  public:
    inline GeometryType geometryType () const
    { 
      return storage_.geometryType(); 
    }
    
    inline int numBaseFunctions () const
    {
      return numDifferentBaseFunctions() * dimRange;
    }

    inline int numDifferentBaseFunctions () const
    {
      return storage_.numBaseFunctions();
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< int, diffOrd > &diffVariable,
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
    inline void jacobian ( const int baseFunction,
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
    inline void jacobianScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarJacobianRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.jacobian( baseFunction, x, phi );
    }
    
    template< int diffOrder, class PointType, class LocalDofVectorType >
    inline void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                           const PointType &x,
                           const LocalDofVectorType& dofs, 
                           RangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase ; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, diffVariable, x, phi ); 
        for( int r = 0; r < dimRange; ++r )
        { 
          ret[ r ] += phi[ 0 ] * dofs[ util_.combinedDof( i , r ) ];
        }
      }
    }

    template< class PointType, class LocalDofVectorType >
    inline void evaluate ( const PointType &x,
                           const LocalDofVectorType& dofs, 
                           RangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase ; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi ); 
        for( int r = 0; r < dimRange; ++r )
        { 
          ret[ r ] += phi[ 0 ] * dofs[ util_.combinedDof( i , r ) ];
        }
      }
    }

    template< class PointType, 
              class GeometryJacobianInverseType,
              class LocalDofVectorType, 
              class GlobalJacobianRangeType>
    inline void jacobian ( const PointType &x,
                           const GeometryJacobianInverseType& gjit, 
                           const LocalDofVectorType& dofs, 
                           GlobalJacobianRangeType &ret ) const
    {
      ret = 0;

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase ; ++i )
      {
        ScalarJacobianRangeType gradPhiRef, gradPhi;
        jacobianScalar( i, x, gradPhiRef );

        FieldMatrixHelper :: multiply( gjit, gradPhiRef[ 0 ], gradPhi[ 0 ] );

        for( int r = 0; r < dimRange; ++r )
          ret[ r ].axpy( dofs[ util_.combinedDof( i , r ) ], gradPhi[ 0 ] );
      }
    }
    
    template< class PointType, class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const RangeType &rangeFactor,
                       LocalDofVectorType& dofs ) const
    {
      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase ; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi ); 
        for( int r = 0; r < dimRange; ++r )
        { 
          dofs[ util_.combinedDof( i , r ) ] += phi[ 0 ] * rangeFactor[ r ];
        }
      }
    }
    
    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       const GlobalJacobianRangeType &jacFactor,
                       LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      FieldMatrixHelper :: multiply( jacFactor, gjit, jacFactorInv );

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase; ++i )
      {
        ScalarJacobianRangeType grad;
        jacobianScalar( i, x, grad );
        for( int r = 0; r < dimRange; ++r )
          dofs[ util_.combinedDof( i , r ) ] += grad[ 0 ] * jacFactorInv[ r ];
      }
    }
 

    template< class PointType, class GeometryJacobianInverseType, 
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    inline void axpy ( const PointType &x,
                       const GeometryJacobianInverseType& gjit, 
                       const RangeType& rangeFactor,
                       const GlobalJacobianRangeType &jacFactor,
                       LocalDofVectorType& dofs ) const
    {
      GlobalJacobianRangeType jacFactorInv;
      FieldMatrixHelper :: multiply( jacFactor, gjit, jacFactorInv );

      const int numScalarBase = numDifferentBaseFunctions ();
      for( int i = 0; i < numScalarBase; ++i )
      {
        ScalarRangeType phi;
        evaluateScalar( i, x, phi );
        ScalarJacobianRangeType grad;
        jacobianScalar( i, x, grad );
        for( int r = 0; r < dimRange; ++r )
          dofs[ util_.combinedDof( i , r ) ] += (phi[ 0 ] * rangeFactor[ r ]) + (grad[ 0 ] * jacFactorInv[ r ]);
      }
    }
 
    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    DUNE_VERSION_DEPRECATED(1,2,remove)
    inline RangeFieldType evaluateSingle ( const int baseFunction,
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
    inline RangeFieldType evaluateGradientSingle( const int baseFunction,
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

#if 1
    /////////////////////////////////////////////////////////////////////////
    //
    //  evaluate and store results in a vector 
    //
    /////////////////////////////////////////////////////////////////////////
    template< class QuadratureType, 
              class LocalDofVectorType,
              class RangeVectorType>
    inline void 
    evaluateRanges ( const QuadratureType& quad,
                     const LocalDofVectorType& dofs,
                     RangeVectorType &rangeVector) const 
    {
#ifdef DUNE_FEM_BASEFUNC_USE_SSE
      typedef typename StorageType :: RangeCacheMatrixType RangeCacheMatrixType;
      const RangeCacheMatrixType& baseFunctionMatrix = storage_.getRangeMatrix( quad );

      const bool faceCachingQuad = (QuadratureType :: codimension == 1) && 
        Conversion< QuadratureType, CachingInterface > :: exists ;

      if( faceCachingQuad ) 
      {
        // apply mapping of quadrature points 
        RangeCacheMatrixType mat = baseFunctionMatrix.getQuadMat( quad );
        mat.mv( dofs, rangeVector, dimRange );  
      }
      else 
      {
        baseFunctionMatrix.mv( dofs, rangeVector, dimRange );
      }
#else 
      typedef typename StorageType :: RangeVectorType RangeVectorType;
      const RangeVectorType& baseFctStorage = storage_.getRangeStorage( quad );

      const size_t numRows = quad.nop();
      const int numDiffBase = numDifferentBaseFunctions();
      assert( numDiffBase * dimRange == dofs.numDofs() );

      assert( baseFctStorage.size() >= numRows );
      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = storage_.applyCaching( quad , row ); 

        assert( baseFctStorage.size() > baseRow );
        assert( baseFctStorage[ baseRow ].size() >= numDiffBase );

        RangeType& result = rangeVector[ row ]; 
        result = 0;

        for( size_t col = 0; col < numDiffBase; ++col ) 
        {
          const ScalarRangeType& phi = baseFctStorage[ baseRow ][ col ];
          size_t colR = util_.combinedDof( col, 0 );
          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            result[ r ] +=  dofs[ colR ] * phi[ 0 ];
          }
        }
      }
#endif
    }
#endif

    template< class QuadratureType, 
              class Geometry,
              class LocalDofVectorType,
              class JacobianRangeVectorType>
    inline void 
    evaluateJacobians ( const QuadratureType& quad,
                        const Geometry& geometry,
                        const LocalDofVectorType& dofs,
                        JacobianRangeVectorType &jacVector ) const 
    {
      assert( jacVector.size() > 0 );
      evaluateJacobians( quad, geometry, dofs, jacVector, jacVector[ 0 ] );
    }

    template< class QuadratureType, 
              class Geometry,
              class LocalDofVectorType,
              class JacobianVectorType,
              class GlobalJacobianRangeType >
    inline void 
    evaluateJacobians ( const QuadratureType& quad,
                        const Geometry& geometry,
                        const LocalDofVectorType& dofs,
                        JacobianVectorType &jacVector,
                        const GlobalJacobianRangeType& ) const 
    {
      typedef typename StorageType :: JacobianRangeVectorType JacobianRangeVectorType;
      const JacobianRangeVectorType& baseFctStorage = storage_.getJacobianStorage( quad );

      const size_t numRows = quad.nop();
      const size_t numDiffBase = numDifferentBaseFunctions();
      assert( numDiffBase * dimRange == dofs.numDofs() );
      assert( baseFctStorage.size() >= numRows );

      const bool affineGeometry = geometry.affine();
      typedef typename Geometry :: JacobianTransposed GeometryJacobianType;
      GeometryJacobianType gjitTmp ;

      const GeometryJacobianType& gjit = 
        ( affineGeometry ) ? 
        geometry.jacobianInverseTransposed( quad.point( 0 )) :
        gjitTmp ;


      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = storage_.applyCaching( quad , row ); 

        // if geometry has non-affine mapping we need to update jacobian inverse
        if( ! affineGeometry ) 
          gjitTmp = geometry.jacobianInverseTransposed( quad.point( row ) );

        assert( baseFctStorage.size() > baseRow );
        assert( baseFctStorage[ baseRow ].size() >= numDiffBase );

        GlobalJacobianRangeType& result = jacVector[ row ]; 
        result = 0;

        // work to do here, need scalar global jacobian 
        ScalarJacobianRangeType gradPhi;

        for( size_t col = 0; col < numDiffBase; ++col ) 
        {
          FieldMatrixHelper :: multiply(gjit,
                                        baseFctStorage[ baseRow ][ col ][ 0 ], 
                                        gradPhi[ 0 ] );

          size_t colR = util_.combinedDof( col, 0 );
          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            result[ r ].axpy( dofs[ colR ], gradPhi[ 0 ] );
          }
        }
      }
    }
   
    ////////////////////////////////////////////////////
    //  axpyRanges 
    ////////////////////////////////////////////////////
    template< class QuadratureType, 
              class RangeVectorType,
              class LocalDofVectorType >
    inline void axpyRanges ( const QuadratureType& quad,
                             const RangeVectorType &rangeFactors,
                             LocalDofVectorType& dofs ) const
    {
#ifdef DUNE_FEM_BASEFUNC_USE_SSE
      typedef typename StorageType :: RangeCacheMatrixType RangeCacheMatrixType;
      const RangeCacheMatrixType& baseFunctionMatrix = storage_.getRangeMatrix( quad );


      const bool faceCachingQuad = (QuadratureType :: codimension == 1) && 
        Conversion< QuadratureType, CachingInterface > :: exists ;

      if( faceCachingQuad ) 
      {
        RangeCacheMatrixType mat = baseFunctionMatrix.getQuadMat( quad );
        // apply mapping of quadrature points 
        mat.umtv( rangeFactors, dofs , dimRange );
      }
      else 
        baseFunctionMatrix.umtv( rangeFactors, dofs , dimRange );

#else 
      typedef typename StorageType :: RangeVectorType RangeVectorType;
      const RangeVectorType& baseFctStorage = storage_.getRangeStorage( quad );

      const size_t numRows = quad.nop();
      const size_t numDiffBase = numDifferentBaseFunctions();
      assert( numDiffBase * dimRange == dofs.numDofs() );

      assert( baseFctStorage.size() >= numRows );
      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = storage_.applyCaching( quad , row ); 

        assert( baseFctStorage.size() > baseRow );
        assert( baseFctStorage[ baseRow ].size() >= numDiffBase );
        const RangeType& factor = rangeFactors[ row ]; 

        for( size_t col = 0; col < numDiffBase; ++col ) 
        {
          const ScalarRangeType& phi = baseFctStorage[ baseRow ][ col ];
          size_t colR = util_.combinedDof( col, 0 );
          for( int r = 0; r < dimRange; ++r , ++colR ) 
          {
            dofs[ colR ] += phi[ 0 ] * factor[ r ];
          }
        }
      }
#endif
    }

    ///////////////////////////////////////////////////////////
    //  applyAxpy Jacobian 
    ///////////////////////////////////////////////////////////
    template< class QuadratureType, 
              class Geometry,
              class JacobianVectorType,
              class LocalDofVectorType >
    inline void axpyJacobians ( const QuadratureType& quad,
                                const Geometry& geometry,
                                const JacobianVectorType &jacVector,
                                LocalDofVectorType& dofs) const 
    {
      assert( jacVector.size() > 0 );
      axpyJacobians( quad, geometry,
                     jacVector, dofs, jacVector[ 0 ] );
    }

    template< class QuadratureType, 
              class Geometry,
              class JacobianVectorType,
              class GlobalJacobianRangeType,
              class LocalDofVectorType >
    inline void axpyJacobians ( const QuadratureType& quad,
                                const Geometry& geometry,
                                const JacobianVectorType &jacVector,
                                LocalDofVectorType& dofs,
                                const GlobalJacobianRangeType& ) const
    {
      typedef typename StorageType :: JacobianRangeVectorType JacobianRangeVectorType;
      const JacobianRangeVectorType& baseFctStorage = storage_.getJacobianStorage( quad );


      const size_t numRows = quad.nop();
      const size_t numDiffBase = numDifferentBaseFunctions();
      assert( numDiffBase * dimRange == dofs.numDofs() );

      typedef typename Geometry :: JacobianTransposed GeometryJacobianType;
      GlobalJacobianRangeType jacFactorInv;

      const bool affineGeometry = geometry.affine();
      GeometryJacobianType gjitTmp ;

      const GeometryJacobianType& gjit = 
        ( affineGeometry ) ? 
        geometry.jacobianInverseTransposed( quad.point( 0 )) :
        gjitTmp ;


      assert( baseFctStorage.size() >= numRows );
      for( size_t row = 0; row < numRows ; ++row )
      {
        const size_t baseRow = storage_.applyCaching( quad , row ); 
        assert( baseFctStorage.size() > baseRow );
        assert( baseFctStorage[ baseRow ].size() >= numDiffBase );

        // if geometry has non-affine mapping we need to update jacobian inverse
        if( ! affineGeometry ) 
          gjitTmp = geometry.jacobianInverseTransposed( quad.point( row ));

        // multiply jacobian factor with geometry inverse 
        FieldMatrixHelper :: multiply( jacVector[ row ], 
                                       gjit,
                                       jacFactorInv );

        for( size_t col = 0; col < numDiffBase; ++col ) 
        {
          size_t colR = util_.combinedDof( col, 0 );
          for( int r = 0; r < dimRange; ++r, ++colR ) 
          {
            dofs[ colR ] += baseFctStorage[ baseRow ][ col ][ 0 ] * jacFactorInv[ r ];
          }
        }
      }
    }

  protected:
    template< int diffOrd, class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const FieldVector< int, diffOrd > &diffVariable,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numDifferentBaseFunctions()) );
      storage_.evaluate( baseFunction, diffVariable, x, phi );
    }
    
    template< class PointType >
    inline void evaluateScalar ( const int baseFunction,
                                 const PointType &x,
                                 ScalarRangeType &phi ) const
    {
      FieldVector< int, 0 > diffVar;
      evaluateScalar( baseFunction, diffVar, x, phi );
    }
  };

  /** \} */

} // end namespace Dune

#endif

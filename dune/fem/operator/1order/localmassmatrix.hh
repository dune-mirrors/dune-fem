#ifndef DUNE_FEM_LOCALMASSMATRIX_HH
#define DUNE_FEM_LOCALMASSMATRIX_HH

//- dune-common includes
#include <dune/common/dynmatrix.hh>
#include <dune/common/fmatrix.hh>

//- dune-geometry includes
#include <dune/geometry/typeindex.hh>

//- dune-fem includes
#include <dune/fem/common/memory.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/misc/checkgeomaffinity.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/allgeomtypes.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    /*! @addtogroup PassHyp
     *  ** @{
    */

    /** \brief Local Mass Matrix inversion implementation, select the correct method in your
        implementation */
    template< class DiscreteFunctionSpace, class VolumeQuadrature >
    class LocalMassMatrixImplementation
    {
      typedef LocalMassMatrixImplementation< DiscreteFunctionSpace, VolumeQuadrature > ThisType;

    public:
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType ctype;
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType RangeFieldType;
      typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

      enum { dimRange = DiscreteFunctionSpaceType :: dimRange };
      enum { localBlockSize = DiscreteFunctionSpaceType :: localBlockSize };
      enum { dgNumDofs = localBlockSize };

      typedef Dune::FieldMatrix< ctype, dgNumDofs, dgNumDofs > DGMatrixType;
      typedef Dune::FieldVector< ctype, dgNumDofs >            DGVectorType;

      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

      typedef typename DiscreteFunctionSpaceType :: IndexSetType IndexSetType;
      typedef typename IndexSetType :: IndexType   IndexType;

      typedef typename DiscreteFunctionSpaceType :: BasisFunctionSetType BasisFunctionSetType;

      typedef typename GridPartType :: GridType GridType;
      typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;
      typedef typename EntityType :: Geometry  Geometry;

      typedef VolumeQuadrature VolumeQuadratureType;

      typedef Fem::GeometryAffinityCheck<VolumeQuadratureType>  GeometryAffinityCheckType;

      //! is true if grid is structured grid
      enum { StructuredGrid = Dune::Capabilities::isCartesian< GridType >::v };

      typedef AllGeomTypes< typename GridPartType :: IndexSetType,GridType> GeometryInformationType;
      typedef typename GeometryInformationType :: DomainType DomainType;

      // use dynamic matrix from dune-common
      typedef Dune::DynamicMatrix< RangeFieldType > MatrixType;
      typedef Dune::DynamicVector< RangeFieldType > VectorType;

    protected:
      std::shared_ptr< const DiscreteFunctionSpaceType > spc_;
      const IndexSetType& indexSet_;

      GeometryInformationType geoInfo_;
      const std::function<int(const int)> volumeQuadratureOrder_;
      const bool affine_;

      mutable DGMatrixType dgMatrix_;
      mutable DGVectorType dgX_, dgRhs_;

      // use dynamic vector from dune-common

      mutable VectorType rhs_, row_;
      mutable MatrixType matrix_;

      mutable std::vector< RangeType > phi_;
      mutable std::vector< RangeType > phiMass_;

      typedef std::pair< std::unique_ptr< MatrixType >, std::unique_ptr< VectorType > > MatrixPairType;
      typedef std::map< const int, MatrixPairType > MassMatrixStorageType;
      typedef std::vector< MassMatrixStorageType > LocalInverseMassMatrixStorageType;

      mutable LocalInverseMassMatrixStorageType localInverseMassMatrix_;

      // index of entity from index set, don't setup mass matrix for the same entity twice
      mutable IndexType lastEntityIndex_;
      mutable unsigned int lastTopologyId_ ;
      // sequence number (obtained from DofManager via the space)
      mutable int sequence_;

      struct NoMassDummyCaller
      {
        typedef Dune::FieldMatrix< ctype, dimRange, dimRange > MassFactorType;

        // return false since we don;t have a mass term
        bool hasMass () const { return false; }

        void mass ( const EntityType &, const VolumeQuadratureType &, int, MassFactorType & ) const
        {}
      };


      bool checkDiagonalMatrix( const MatrixType& matrix ) const
      {
        const int rows = matrix.rows();
        for( int r=0; r<rows; ++r )
        {
          for( int c=0; c<r; ++c ) // the mass matrix is symmetric
          {
            // if we find one off diagonal non-zero return false
            if( std::abs(matrix[r][c]) > 1e-12 )
              return false;
          }
        }
        return true;
      }

      template< class BasisFunctionSet >
      MatrixPairType&
      getLocalInverseMassMatrix ( const EntityType &entity, const Geometry &geo,
                                  const BasisFunctionSet &basisSet, int numBasisFct ) const
      {
        const GeometryType geomType = geo.type();
        typedef typename MassMatrixStorageType::iterator iterator;
        MassMatrixStorageType &massMap = localInverseMassMatrix_[ GlobalGeometryTypeIndex::index( geomType ) ];

        auto it = massMap.find( numBasisFct );
        if( it == massMap.end() )
        {
          std::pair< iterator, bool > insertPair = massMap.insert( std::make_pair( numBasisFct, MatrixPairType(nullptr,nullptr) ) );
          it = insertPair.first;
          insertPair.first->second.first.reset( new MatrixType( numBasisFct, numBasisFct, 0.0 ));
          MatrixType& matrix = insertPair.first->second.first.operator *();
          VolumeQuadratureType volQuad( entity, volumeQuadratureOrder( entity ) );
          buildMatrixNoMassFactor( entity, geo, basisSet, volQuad, numBasisFct, matrix, false );
          try {
            matrix.invert();
          }
          catch ( Dune::FMatrixError &e )
          {
            std::cerr << "Matrix is singular:" << std::endl << matrix << std::endl;
            std::terminate();
          }
          const bool diagonal = checkDiagonalMatrix( matrix );
          // store diagonal if matrix is diagonal
          if( diagonal )
          {
            insertPair.first->second.second.reset( new VectorType( matrix.rows() ) );
            VectorType& diag = insertPair.first->second.second.operator *();
            const int rows = matrix.rows();
            for( int row=0; row<rows; ++row )
            {
              diag[ row ] = matrix[ row ][ row ];
            }
          }
        }

        return it->second;
      }

      template< class MassCaller, class BasisFunctionSet >
      MatrixType &getLocalInverseMassMatrixDefault ( MassCaller &caller, const EntityType &entity,
                                                     const Geometry &geo, const BasisFunctionSet &basisSet ) const
      {
        const int numDofs = basisSet.size();
        // if sequence changed or entity index changed: recompute mass matrix
        if( entityHasChanged( entity ) || (numDofs != int( matrix_.rows())) )
        {
          // resize temporary memory if necessary
          if( numDofs != int( matrix_.rows() ) )
            matrix_.resize( numDofs, numDofs );

          buildMatrix( caller, entity, geo, basisSet, numDofs, matrix_ );
          matrix_.invert();
        }

        // make sure that rhs_ has the correct size
        assert( numDofs == int( matrix_.rows() ) );
        return matrix_;
      }

      // return number of max non blocked dofs
      int maxNumDofs () const
      {
        return space().blockMapper().maxNumDofs() * localBlockSize;
      }

      //! return appropriate quadrature order, default is 2 * order()
      int maxVolumeQuadratureOrder () const
      {
        return volumeQuadratureOrder_( space().order() );
      }

    public:
      //! return appropriate quadrature order, default is 2 * order(entity)
      int volumeQuadratureOrder ( const EntityType &entity ) const
      {
        return volumeQuadratureOrder_( space().order( entity ) );
      }

      //! constructor taking space and volume quadrature order
      explicit LocalMassMatrixImplementation ( const DiscreteFunctionSpaceType &spc, int volQuadOrd )
        : LocalMassMatrixImplementation( spc, [volQuadOrd](const int order) { return volQuadOrd; } )
      {}

      //! constructor taking space and volume quadrature order
      explicit LocalMassMatrixImplementation ( const DiscreteFunctionSpaceType &spc,
              std::function<int(const int)> volQuadOrderFct = [](const int order) { return Capabilities::DefaultQuadrature< DiscreteFunctionSpaceType >::volumeOrder(order); } )
        : spc_( referenceToSharedPtr( spc ) )
        , indexSet_( space().indexSet() )
        , geoInfo_( indexSet_ )
        , volumeQuadratureOrder_ ( volQuadOrderFct )
        , affine_ ( setup() )
        , rhs_(), row_(), matrix_()
        , phi_( maxNumDofs() )
        , phiMass_( maxNumDofs() )
        , localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) )
        , lastEntityIndex_( -1 )
        , lastTopologyId_( ~0u )
        , sequence_( -1 )
      {}

      //! copy constructor
      LocalMassMatrixImplementation ( const ThisType &other )
      : spc_(other.spc_),
        indexSet_( space().indexSet() ),
        geoInfo_( indexSet_ ),
        volumeQuadratureOrder_( other.volumeQuadratureOrder_ ),
        affine_( other.affine_ ),
        rhs_( other.rhs_ ), row_( other.row_ ), matrix_( other.matrix_ ),
        phi_( other.phi_ ),
        phiMass_( other.phiMass_ ),
        localInverseMassMatrix_( GlobalGeometryTypeIndex :: size( GridType::dimension ) ),
        lastEntityIndex_( other.lastEntityIndex_ ),
        lastTopologyId_( other.lastTopologyId_ ),
        sequence_( other.sequence_ )
      {}

    public:
      //! returns true if geometry mapping is affine
      bool affine () const { return affine_; }

      //! return mass factor for diagonal mass matrix
      double getAffineMassFactor(const Geometry& geo) const
      {
        return geoInfo_.referenceVolume( geo.type() ) / geo.volume();
      }

      template< class BasisFunctionSet >
      bool checkInterpolationBFS(const BasisFunctionSet &bfs) const
      {
        static const int quadPointSetId = SelectQuadraturePointSetId< VolumeQuadratureType >::value;
        // if BasisFunctionSet does not have an static int member called pointSetId then this will be -1
        static const int basePointSetId = detail::SelectPointSetId< BasisFunctionSet >::value;
        // for Lagrange-type basis evaluated on interpolation points
        // this is the Kronecker delta, so the mass matrix is diagonal even
        // on non affine grids
        if constexpr ( quadPointSetId == basePointSetId )
        {
          const unsigned int numShapeFunctions = bfs.size() / dimRange;
          return VolumeQuadratureType( bfs.entity(), volumeQuadratureOrder( bfs.entity() ) )
                     .isInterpolationQuadrature(numShapeFunctions);
        }
        return false;
      }

      //! apply local dg mass matrix to local function lf
      //! using the massFactor method of the caller
      template< class MassCaller, class BasisFunctionSet, class LocalFunction >
      void applyInverse ( MassCaller &caller, const EntityType &entity, const BasisFunctionSet &basisFunctionSet, LocalFunction &lf ) const
      {
        Geometry geo = entity.geometry();
        if( ( affine() || geo.affine() || checkInterpolationBFS(basisFunctionSet) )
            && !caller.hasMass() )
          applyInverseLocally( entity, geo, basisFunctionSet, lf );
        else
          applyInverseDefault( caller, entity, geo, basisFunctionSet, lf );
      }

      //! apply local dg mass matrix to local function lf
      //! using the massFactor method of the caller
      template< class MassCaller, class LocalFunction >
      void applyInverse ( MassCaller &caller, const EntityType &entity, LocalFunction &lf ) const
      {
        applyInverse( caller, entity, lf.basisFunctionSet(), lf );
      }

      //! apply local dg mass matrix to local function lf without mass factor
      template< class LocalFunction >
      void applyInverse ( const EntityType &entity, LocalFunction &lf ) const
      {
        NoMassDummyCaller caller;
        applyInverse( caller, entity, lf.basisFunctionSet(), lf );
      }
      template< class BasisFunctionSet, class LocalFunction >
      void applyInverse ( const EntityType &entity, const BasisFunctionSet &basisFunctionSet, LocalFunction &lf ) const
      {
        NoMassDummyCaller caller;
        applyInverse( caller, entity, basisFunctionSet, lf );
      }


      //! apply local dg mass matrix to local function lf without mass factor
      template< class LocalFunction >
      void applyInverse ( LocalFunction &lf ) const
      {
        applyInverse( lf.entity(), lf );
      }

      //! compute localMatrix * M^-1
      template< class LocalMatrix >
      void rightMultiplyInverse ( LocalMatrix &localMatrix ) const
      {
        const EntityType &entity = localMatrix.rangeEntity();
        Geometry geo = entity.geometry();
        if( ( affine() || geo.affine() || checkInterpolationBFS(localMatrix.rangeBasisFunctionSet())) )
          rightMultiplyInverseLocally( entity, geo, localMatrix );
        else
          rightMultiplyInverseDefault( entity, geo, localMatrix );
      }

      //! compute M^-1 * localMatrix
      template< class LocalMatrix >
      void leftMultiplyInverse ( LocalMatrix &localMatrix ) const
      {
        const EntityType &entity = localMatrix.domainEntity();
        Geometry geo = entity.geometry();
        if( ( affine() || geo.affine() || checkInterpolationBFS(localMatrix.rangeBasisFunctionSet())) )
          leftMultiplyInverseLocally( entity, geo, localMatrix );
        else
          leftMultiplyInverseDefault( entity, geo, localMatrix );
      }

      const DiscreteFunctionSpaceType &space () const { return *spc_; }

      /////////////////////////////////////////////
      // end of public methods
      /////////////////////////////////////////////

    protected:
      ///////////////////////////////////////////////////////////
      //  applyInverse for DG spaces
      ///////////////////////////////////////////////////////////
      template< class MassCaller, class BasisFunctionSet, class LocalFunction >
      void applyInverseDgOrthoNormalBasis ( MassCaller &caller, const EntityType &entity,
                                            const BasisFunctionSet &basisFunctionSet, LocalFunction &lf ) const
      {
        Geometry geo = entity.geometry();
        assert( dgNumDofs == lf.size() );

        // affine_ can be a static information
        const bool isAffine = affine() || geo.affine();
        // make sure that for affine grids the geometry info is also correct
        assert( affine() ? geo.affine() : true );

        // in case of affine mappings we only have to multiply with a factor
        if( isAffine && !caller.hasMass() )
        {
          const double massVolInv = getAffineMassFactor( geo );

          // apply inverse mass matrix
          for( int l = 0; l < dgNumDofs; ++l )
            lf[ l ] *= massVolInv;
        }
        else
        {
          // copy local function to right hand side
          for( int l = 0; l < dgNumDofs; ++l )
            dgRhs_[ l ] = lf[ l ];

          buildMatrix( caller, entity, geo, basisFunctionSet, dgNumDofs, dgMatrix_ );
          dgMatrix_.solve( dgX_, dgRhs_ );

          // copy back to local function
          for( int l = 0; l < dgNumDofs; ++l )
            lf[ l ] = dgX_[ l ];
        }
      }

      //! compute localMatrix * M^-1
      template< class LocalMatrix >
      void rightMultiplyInverseDgOrthoNormalBasis ( LocalMatrix &localMatrix ) const
      {
        const EntityType &entity = localMatrix.rangeEntity();
        Geometry geo = entity.geometry();
        assert( dgNumDofs == localMatrix.columns() );

        // in case of affine mappings we only have to multiply with a factor
        if( affine() || geo.affine() )
        {
          localMatrix.scale( getAffineMassFactor( geo ) );
        }
        else
        {
          NoMassDummyCaller caller;
          buildMatrix( caller, entity, geo, localMatrix.rangeBasisFunctionSet(), dgNumDofs, dgMatrix_ );
          dgMatrix_.invert();

          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            for( int j = 0; j < dgNumDofs; ++j )
              dgRhs_[ j ] = localMatrix.get( i, j );
            dgMatrix_.mtv( dgRhs_, dgX_ );
            for( int j = 0; j < dgNumDofs; ++j )
              localMatrix.set( i, j, dgX_[ j ] );
          }
        }
      }

      //! compute M^-1 * localMatrix
      template< class LocalMatrix >
      void leftMultiplyInverseDgOrthoNormalBasis ( LocalMatrix &localMatrix ) const
      {
        const EntityType &entity = localMatrix.domainEntity();
        Geometry geo = entity.geometry();
        assert( dgNumDofs == localMatrix.columns() );

        // in case of affine mappings we only have to multiply with a factor
        if( affine() || geo.affine() )
        {
          localMatrix.scale( getAffineMassFactor( geo ) );
        }
        else
        {
          NoMassDummyCaller caller;
          buildMatrix( caller, entity, geo, localMatrix.domainBasisFunctionSet(), dgNumDofs, dgMatrix_ );
          dgMatrix_.invert();

          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            for( int j = 0; j < dgNumDofs; ++j )
              dgRhs_[ j ] = localMatrix.get( i, j );
            dgMatrix_.mv( dgRhs_, dgX_ );
            for( int j = 0; j < dgNumDofs; ++j )
              localMatrix.set( i, j, dgX_[ j ] );
          }
        }
      }

      //! returns true if the entity has been changed
      bool entityHasChanged( const EntityType& entity ) const
      {
        // don't compute matrix new for the same entity
        const int currentSequence   = space().sequence();
        const unsigned int topologyId = entity.type().id();
        const IndexType entityIndex = indexSet_.index( entity ) ;

        // check whether sequence has been updated
        if( sequence_ != currentSequence ||
            lastEntityIndex_ != entityIndex ||
            lastTopologyId_  != topologyId )
        {
          // update identifiers
          lastEntityIndex_ = entityIndex ;
          sequence_        = currentSequence;
          lastTopologyId_  = topologyId ;

          return true ;
        }
        else
          // the entity did not change
          return false ;
      }

      ///////////////////////////////////////////////////////////
      //  standard applyInverse method
      ///////////////////////////////////////////////////////////
      //! apply local mass matrix to local function lf
      //! using the massFactor method of the caller
      template< class MassCaller, class BasisFunctionSet, class LocalFunction >
      void applyInverseDefault ( MassCaller &caller, const EntityType &entity,
                                 const Geometry &geo, const BasisFunctionSet &basisFunctionSet, LocalFunction &lf ) const
      {
        // get local inverted mass matrix
        MatrixType &invMassMatrix
          = getLocalInverseMassMatrixDefault ( caller, entity, geo, basisFunctionSet );

        // copy local function to right hand side
        const int numDofs = lf.size();
        rhs_.resize( numDofs );
        for( int l = 0; l < numDofs; ++l )
          rhs_[ l ] = lf[ l ];

        // apply inverse to right hand side and store in lf
        multiply( numDofs, invMassMatrix, rhs_, lf );
      }

      //! compute localMatrix * M^-1
      template< class LocalMatrix >
      void rightMultiplyInverseDefault ( const EntityType &entity, const Geometry &geo, LocalMatrix &localMatrix ) const
      {
        NoMassDummyCaller caller;
        MatrixType &invMassMatrix
          = getLocalInverseMassMatrixDefault ( caller, entity, geo, localMatrix.rangeBasisFunctionSet() );

        const int cols = localMatrix.columns();
        rhs_.resize( cols );
        row_.resize( cols );

        const int rows = localMatrix.rows();
        for( int i = 0; i < rows; ++i )
        {
          // get i-th row from localMatrix
          for( int j = 0; j < cols; ++j )
            rhs_[ j ] = localMatrix.get( i, j );

          // multiply with all columns of inverse mass matrix
          invMassMatrix.mtv( rhs_, row_ );

          // store as i-th row in localMatrix
          for( int j = 0; j < cols; ++j )
            localMatrix.set( i, j, row_[ j ] );
        }
      }

      //! compute M^-1 * localMatrix
      template< class LocalMatrix >
      void leftMultiplyInverseDefault ( const EntityType &entity, const Geometry &geo, LocalMatrix &localMatrix ) const
      {
        NoMassDummyCaller caller;
        MatrixType &invMassMatrix
          = getLocalInverseMassMatrixDefault ( caller, entity, geo, localMatrix.rangeBasisFunctionSet() );

        const int cols = localMatrix.columns();
        rhs_.resize( cols );
        row_.resize( cols );

        const int rows = localMatrix.rows();
        for( int i = 0; i < rows; ++i )
        {
          // get i-th column from localMatrix
          for( int j = 0; j < cols; ++j )
            rhs_[ j ] = localMatrix.get( j, i );

          // multiply with all rows in inverse mass matrix
          invMassMatrix.mv( rhs_, row_ );

          // store as i-th column in localMatrix
          for( int j = 0; j < cols; ++j )
            localMatrix.set( j, i, row_[ j ] );
        }
      }


      ///////////////////////////////////////////////////////////
      //  local applyInverse method for affine geometries
      ///////////////////////////////////////////////////////////
      //! apply local mass matrix to local function lf
      template< class BasisFunctionSet, class LocalFunction >
      void applyInverseLocally ( const EntityType &entity,
                                 const Geometry &geo, const BasisFunctionSet &basisFunctionSet, LocalFunction &lf ) const
      {
        const int numDofs = lf.size();

        // get local inverted mass matrix
        MatrixPairType& matrixPair =
          getLocalInverseMassMatrix( entity, geo, basisFunctionSet, numDofs );

        // if diagonal exists then matrix is in diagonal form
        if( matrixPair.second )
        {
          const VectorType& diagonal = *matrixPair.second;
          assert( int(diagonal.size()) == numDofs );

          VolumeQuadratureType volQuad( entity, volumeQuadratureOrder( entity ) );

          const int nop = volQuad.nop();
          assert(nop*dimRange == numDofs);
          for( int l=0, qt = 0; qt < nop; ++qt )
          {
            const auto intel = geo.integrationElement( volQuad.point(qt) );
            for (int r = 0; r < dimRange; ++r, ++l )
            {
              lf[ l ] *= diagonal[ l ] / intel;
            }
          }
        }
        else
        {
          const double massVolInv = getAffineMassFactor( geo );
          // copy local function to right hand side
          // and apply inverse mass volume fraction
          rhs_.resize( numDofs );
          for( int l = 0; l < numDofs; ++l )
            rhs_[ l ] = lf[ l ] * massVolInv;

          const MatrixType& invMassMatrix = *matrixPair.first;
          // apply inverse local mass matrix and store in lf
          multiply( numDofs, invMassMatrix, rhs_, lf );
        }
      }

      template <class LocalMatrix>
      const VectorType&
      setupInverseDiagonal( const EntityType &entity, const Geometry &geo,
                            const VectorType& refElemDiagonal,
                            LocalMatrix &localMatrix ) const
      {
        const int cols = localMatrix.columns();

        assert( int(refElemDiagonal.size()) == cols );

        VolumeQuadratureType volQuad( entity, volumeQuadratureOrder( entity ) );

        VectorType& elementDiagonal = rhs_;
        elementDiagonal.resize( cols );

        const int nop = volQuad.nop();
        assert(nop*dimRange == cols);
        for( int l = 0, qt = 0; qt < nop; ++qt )
        {
          const auto intel = geo.integrationElement( volQuad.point(qt) );
          for (int r = 0; r < dimRange; ++r,++l )
          {
            elementDiagonal[ l ] = refElemDiagonal[ l ] / intel;
          }
        }
        return elementDiagonal;
      }

      template< class LocalMatrix >
      void rightMultiplyInverseLocally ( const EntityType &entity, const Geometry &geo, LocalMatrix &localMatrix ) const
      {
        const int cols = localMatrix.columns();
        MatrixPairType& matrixPair =
          getLocalInverseMassMatrix( entity, geo, localMatrix.rangeBasisFunctionSet(), cols );

        // if diagonal exists then matrix is in diagonal form
        // stored as inverse on the reference element
        if( matrixPair.second )
        {
          const VectorType& elementDiagonal =
            setupInverseDiagonal( entity, geo, *matrixPair.second, localMatrix );

          row_.resize( cols );
          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            // get i-th row from localMatrix
            // and multiply with diagonal of inverse mass matrix
            for( int j = 0; j < cols; ++j )
              row_[ j ] = elementDiagonal[ j ] * localMatrix.get( i, j );

            // store as i-th row in localMatrix
            for( int j = 0; j < cols; ++j )
              localMatrix.set( i, j, row_[ j ] );
          }
        }
        else
        {
          const MatrixType &invMassMatrix = *matrixPair.first;

          const double massVolInv = getAffineMassFactor( geo );

          rhs_.resize( cols );
          row_.resize( cols );

          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            // get i-th row from localMatrix
            // and multiply with diagonal of inverse mass matrix
            for( int j = 0; j < cols; ++j )
              rhs_[ j ] = localMatrix.get( i, j ) * massVolInv;

            // multiply with all columns of inverse mass matrix
            invMassMatrix.mtv( rhs_, row_ );

            // store as i-th row of localMatrix
            for( int j = 0; j < cols; ++j )
              localMatrix.set( i, j, row_[ j ] );
          }
        }
      }

      //! compute M^-1 * localMatrix
      template< class LocalMatrix >
      void leftMultiplyInverseLocally ( const EntityType &entity, const Geometry &geo, LocalMatrix &localMatrix ) const
      {
        const int cols = localMatrix.columns();
        MatrixPairType& matrixPair =
          getLocalInverseMassMatrix( entity, geo, localMatrix.rangeBasisFunctionSet(), cols );

        // if diagonal exists then matrix is in diagonal form
        // stored as inverse on the reference element
        if( matrixPair.second )
        {
          const VectorType& elementDiagonal =
            setupInverseDiagonal( entity, geo, *matrixPair.second, localMatrix );

          row_.resize( cols );
          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            // get i-th column from localMatrix
            // and multiply with diagonal of inverse mass matrix
            for( int j = 0; j < cols; ++j )
              row_[ j ] = elementDiagonal[ j ] * localMatrix.get( j, i );

            // store as i-th column of localMatrix
            for( int j = 0; j < cols; ++j )
              localMatrix.set( j, i, row_[ j ] );
          }
        }
        else
        {
          const MatrixType &invMassMatrix = *matrixPair.first;

          const double massVolInv = getAffineMassFactor( geo );

          rhs_.resize( cols );
          row_.resize( cols );

          const int rows = localMatrix.rows();
          for( int i = 0; i < rows; ++i )
          {
            // get i-th column from localMatrix
            for( int j = 0; j < cols; ++j )
              rhs_[ j ] = localMatrix.get( j, i ) * massVolInv;

            // apply to all rows of inverse mass matrix
            invMassMatrix.mv( rhs_, row_ );

            // store as i-th column of localMatrix
            for( int j = 0; j < cols; ++j )
              localMatrix.set( j, i, row_[ j ] );
          }
        }
      }


      //! setup and return affinity
      bool setup () const
      {
        // for structured grids this is always true
        if( StructuredGrid )
          return true;

        // get types for codim 0
        const std::vector<Dune::GeometryType>& geomTypes = geoInfo_.geomTypes(0);

        // for simplices we also have affine mappings
        if( (geomTypes.size() == 1) && geomTypes[0].isSimplex() )
        {
          return true;
        }

        // otherwise use geometry affinity
        return false ;
      }

      //! build local mass matrix
      template< class MassCaller, class Matrix >
      void buildMatrix ( MassCaller &caller, const EntityType &entity,
                         const Geometry &geo, const BasisFunctionSetType &set,
                         std::size_t numDofs, Matrix &matrix ) const
      {
        assert( numDofs == set.size() );

        // clear matrix
        matrix = 0;

        // create quadrature
        VolumeQuadratureType volQuad( entity, volumeQuadratureOrder( entity ) );

        if( caller.hasMass() )
          buildMatrixWithMassFactor( caller, entity, geo, set, volQuad, numDofs, matrix );
        else
          buildMatrixNoMassFactor( entity, geo, set, volQuad, numDofs, matrix );
      }

      //! build local mass matrix with mass factor
      template <class Matrix>
      void buildMatrixNoMassFactor(
                       const EntityType& en,
                       const Geometry& geo,
                       const BasisFunctionSetType& set,
                       const VolumeQuadratureType& volQuad,
                       const int numDofs,
                       Matrix& matrix,
                       const bool applyIntegrationElement = true ) const
      {
        const int volNop = volQuad.nop();
        for(int qp=0; qp<volNop; ++qp)
        {
          // calculate integration weight
          const double intel = ( applyIntegrationElement ) ?
              ( volQuad.weight(qp) * geo.integrationElement(volQuad.point(qp)) ) : volQuad.weight(qp) ;

          // eval basis functions
          set.evaluateAll(volQuad[qp], phi_);

          for(int m=0; m<numDofs; ++m)
          {
            const RangeType& phi_m = phi_[m];
            const ctype val = intel * (phi_m * phi_m);
            matrix[m][m] += val;

            for(int k=m+1; k<numDofs; ++k)
            {
              const ctype val = intel * (phi_m * phi_[k]);
              matrix[m][k] += val;
              matrix[k][m] += val;
            }
          }
        }
      }

      //! build local mass matrix with mass factor
      template <class MassCallerType, class Matrix>
      void buildMatrixWithMassFactor(
                       MassCallerType& caller,
                       const EntityType& en,
                       const Geometry& geo,
                       const BasisFunctionSetType& set,
                       const VolumeQuadratureType& volQuad,
                       const int numDofs,
                       Matrix& matrix) const
      {
        typedef typename MassCallerType :: MassFactorType MassFactorType;
        MassFactorType mass;

        const int volNop = volQuad.nop();
        for(int qp=0; qp<volNop; ++qp)
        {
          // calculate integration weight
          const double intel = volQuad.weight(qp)
             * geo.integrationElement(volQuad.point(qp));

          // eval basis functions
          set.evaluateAll( volQuad[qp], phi_);

          // call mass factor
          caller.mass( en, volQuad, qp, mass);

          // apply mass matrix to all basis functions
          for(int m=0; m<numDofs; ++m)
          {
            mass.mv( phi_[m], phiMass_[m] );
          }

          // add values to matrix
          for(int m=0; m<numDofs; ++m)
          {
            for(int k=0; k<numDofs; ++k)
            {
              matrix[m][k] += intel * (phiMass_[m] * phi_[k]);
            }
          }
        }
      }

      // implement matvec with matrix (mv of densematrix is too stupid)
      template <class Matrix, class Rhs, class X>
      void multiply( const int size,
                     const Matrix& matrix,
                     const Rhs& rhs,
                     X& x ) const
      {
        assert( (int) matrix.rows() == size );
        assert( (int) matrix.cols() == size );
        assert( (int) rhs.size() == size );

        for( int row = 0; row < size; ++ row )
        {
          RangeFieldType sum = 0;
          // get matrix row
          typedef typename Matrix :: const_row_reference  MatRow;
          MatRow matRow = matrix[ row ];

          // multiply row with right hand side
          for( int col = 0; col < size; ++ col )
          {
            sum += matRow[ col ] * rhs[ col ];
          }

          // set to result to result vector
          x[ row ] = sum;
        }
      }
    };



    // LocalMassMatrix
    // ---------------

    /** \brief Local Mass Matrix for arbitrary spaces */
    template< class DiscreteFunctionSpace, class VolumeQuadrature >
    class LocalMassMatrix
    : public LocalMassMatrixImplementation< DiscreteFunctionSpace, VolumeQuadrature >
    {
      typedef LocalMassMatrixImplementation< DiscreteFunctionSpace, VolumeQuadrature > BaseType;

    public:
      // copy base class constructors
      using BaseType :: LocalMassMatrixImplementation;
    };



    ///////////////////////////////////////////////////////////////////
    //
    //  DG LocalMassMatrix Implementation
    //
    ///////////////////////////////////////////////////////////////////

    /** \brief DG Local Mass Matrix for arbitrary spaces */
    template< class DiscreteFunctionSpace, class VolumeQuadrature >
    class LocalMassMatrixImplementationDgOrthoNormal
    : public LocalMassMatrixImplementation< DiscreteFunctionSpace, VolumeQuadrature >
    {
      typedef LocalMassMatrixImplementation< DiscreteFunctionSpace, VolumeQuadrature > BaseType;

    public:
      typedef typename BaseType :: EntityType  EntityType;

      // copy base class constructors
      using BaseType :: LocalMassMatrixImplementation;

      //! apply local dg mass matrix to local function lf
      //! using the massFactor method of the caller
      template <class MassCallerType, class BasisFunction, class LocalFunctionType>
      void applyInverse(MassCallerType& caller,
                        const EntityType& en,
                        const BasisFunction &basisFunction,
                        LocalFunctionType& lf) const
      {
        BaseType :: applyInverseDgOrthoNormalBasis( caller, en, basisFunction, lf );
      }
      template <class MassCallerType, class LocalFunctionType>
      void applyInverse(MassCallerType& caller,
                        const EntityType& en,
                        LocalFunctionType& lf) const
      {
        BaseType :: applyInverseDgOrthoNormalBasis( caller, en, lf.basisFunctionSet(), lf );
      }

      //! apply local dg mass matrix to local function lf without mass factor
      template <class LocalFunctionType>
      void applyInverse(const EntityType& en,
                        LocalFunctionType& lf) const
      {
        typename BaseType :: NoMassDummyCaller caller;
        applyInverse(caller, en, lf );
      }

      template <class BasisFunction, class LocalFunctionType>
      void applyInverse(const EntityType& en,
                        const BasisFunction &basisFunction,
                        LocalFunctionType& lf) const
      {
        typename BaseType :: NoMassDummyCaller caller;
        applyInverse(caller, en, basisFunction, lf );
      }

      //! apply local dg mass matrix to local function lf without mass factor
      template< class LocalFunction >
      void applyInverse ( LocalFunction &lf ) const
      {
        applyInverse( lf.entity(), lf.basisFunctionSet(), lf );
      }

      //! compute localMatrix * M^-1
      template< class LocalMatrix >
      void rightMultiplyInverse ( LocalMatrix &localMatrix ) const
      {
        BaseType::rightMultiplyInverseDgOrthoNormalBasis( localMatrix );
      }

      //! compute M^-1 * localMatrix
      template< class LocalMatrix >
      void leftMultiplyInverse ( LocalMatrix &localMatrix ) const
      {
        BaseType::leftMultiplyInverseDgOrthoNormalBasis( localMatrix );
      }
    };

//! @}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALMASSMATRIX_HH

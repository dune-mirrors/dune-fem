#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

//- Dune includes
#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  //! \brief The Laplace operator
  template< class DiscreteFunctionImp, class TensorImp >
  class LaplaceFEOp
  : public Operator< typename DiscreteFunctionImp :: RangeFieldType,
                     typename DiscreteFunctionImp :: RangeFieldType,
                     DiscreteFunctionImp,
                     DiscreteFunctionImp >
  {
  public:
    //! type of discrete functions
    typedef DiscreteFunctionImp DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType :: FunctionSpaceType
      DiscreteFunctionSpaceType;
        
    //! type of jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;
    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;
    //! type of range vectors
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };
        
    //! type of tensor
    typedef TensorImp TensorType;

    //! type of system matrix
    typedef SparseRowMatrix< RangeFieldType > MatrixType;

    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
            
  private:
    typedef LaplaceFEOp< DiscreteFunctionType, TensorType > ThisType;

  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

    //! pointer to the system matrix 
    mutable MatrixType *matrix_;
 
    //! flag indicating whether the system matrix has been assembled
    mutable bool matrix_assembled_;
      
    TensorType *const stiffTensor_;

  private:
    mutable JacobianRangeType grad;

    enum { maxBaseFunctions = 100 };
    mutable JacobianRangeType mygrad[ maxBaseFunctions ];
 
  public:
    //! constructor
    LaplaceFEOp( const DiscreteFunctionSpaceType &discreteFunctionSpace )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      matrix_( 0 ),
      matrix_assembled_( false ),
      stiffTensor_( 0 )
    {
    }
        
    //! constructor
    LaplaceFEOp( TensorType &stiffTensor,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      matrix_( 0 ),
      matrix_assembled_( false ),
      stiffTensor_( &stiffTensor )
    { 
    }

    virtual ~LaplaceFEOp ()
    {
      if( matrix_ != 0 )
        delete matrix_;
    }

    //! \brief apply the operator
    virtual void operator() ( const DiscreteFunctionType &u, 
                              DiscreteFunctionType &w ) const 
    {
      systemMatrix().apply( u, w );
    }
  
    /*! \brief obtain a reference to the system matrix 
     *
     *  The assembled matrix is returned. If the system matrix has not been
     *  assembled, yet, the assembly is performed.
     *
     *  \returns a reference to the system matrix
     */
    MatrixType &systemMatrix () const
    {
      if( !matrix_assembled_ )
        assemble();
      return *matrix_;
    }

    //! print the system matrix into a stream
    void print ( std :: ostream out = std :: cout ) const 
    {
      systemMatrix().print( out );
    }

    const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    /*! 
     *   assemble: perform grid-walkthrough and assemble global matrix
     * 
     *   If the matrix storage is 
     *   not allocated, new storage is allocated by newEmptyMatrix.
     *   the begin and end iterators are determined and the assembling
     *   of the global matrix initiated by call of assembleOnGrid and 
     *   bndCorrectOnGrid. The assemled flag is set. 
     */
    void assemble () const 
    {
      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      // if the matrix has not been allocated, allocate it.
      if( this->matrix_ == 0 )
      {
        const int size = dfSpace.size();
        const int numNonZero = 8 * (1 << dimension) * polynomialOrder;
        matrix_ = new MatrixType( size, size, numNonZero );
        assert( matrix_ != 0 );
      }

      matrix_->clear();
      assembleOnGrid();
      boundaryCorrectOnGrid();

      matrix_assembled_ = true;
    }

  protected:
    /*! perform grid walkthrough and assemble matrix
     *
     *  For each element, the local element matrix is determined into the
     *  given local matrix storage and distributed into the global matrix.
     *  Distribution is performed by an add(row,col,val) method on the 
     *  global matrix class.
     */
    void assembleOnGrid () const
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
  
      FieldMatrix< RangeFieldType, maxBaseFunctions, maxBaseFunctions > matrix;

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
        assembleOnEntity( *it, matrix );
    }

    /*! perform matrix assemble for one entity
     *
     *  \param[in] entity entity for current local update
     *  \param{in] localMatrix local matrix storate
     */
    template< class EntityType, class LocalMatrixImp >
    void assembleOnEntity ( const EntityType &entity, 
                            LocalMatrixImp &localMatrix ) const
    {
      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      //const BaseFunctionSetType baseSet = dfSpace.baseFunctionSet( entity );
      //const int numBaseFunctions = baseSet.numBaseFunctions();
      
      // obtain local matrix
      const unsigned int size
        = assembleLocalMatrix( entity, localMatrix );

      for( unsigned int i = 0; i < size; ++i )
      { 
        const unsigned int row = dfSpace.mapToGlobal( entity , i );
        for( unsigned int j = 0; j < size; ++j )
        {
          const unsigned int col = dfSpace.mapToGlobal( entity, j );
          matrix_->add( row, col, localMatrix[ i ][ j ] );
        }
      }
    }

    /*! assemble the local matrix
     *
     *  In the Finite Element Method, most base functions are zero on a given
     *  entity. To build the matrix on one entity, we only need to store a very
     *  limited amount of matrix entries, the socalled local matrix.
     */
    template< class  EntityType, class LocalMatrixType >
    unsigned int assembleLocalMatrix ( const EntityType &entity,
                                       LocalMatrixType &matrix ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType baseSet
        = dfSpace.baseFunctionSet( entity );
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
            
      assert( numBaseFunctions <= maxBaseFunctions );
      for( unsigned int i = 0; i < numBaseFunctions; ++i )
        for( unsigned int j = 0; j <= i; ++j ) 
          matrix[ j ][ i ] = 0;

      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const unsigned int numQuadraturePoints = quad.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt ) {
        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimension, dimension > &inv
                    = geometry.jacobianInverseTransposed( quad.point( pt ) );
        const double volume = geometry.integrationElement( quad.point( pt ) );

        for( unsigned int i = 0; i < numBaseFunctions; ++i ) {
          baseSet.jacobian( i, quad, pt, mygrad[ i ] ); 
      
          // multiply with transpose of jacobian inverse 
          mygrad[ i ][ 0 ] = FMatrixHelp :: mult( inv, mygrad[ i ][ 0 ] );
        }
        
        RangeFieldType weight = quad.weight( pt ) * volume;
        if( stiffTensor_ ) {
          RangeType phi;
          stiffTensor_->evaluate( geometry.global( quad.point( pt ) ), phi );
          weight *= phi[ 0 ];
        }
        for( unsigned int i = 0; i < numBaseFunctions; ++i ) 
          for ( unsigned int j = 0; j <= i; ++j ) 
            matrix[ j ][ i ] += (mygrad[ i ][ 0 ] * mygrad[ j ][ 0 ]) * weight;
      }
      
      // symmetrize matrix
      for( unsigned int i = 0; i < numBaseFunctions; ++i ) 
        for( unsigned int j = numBaseFunctions; j > i; --j ) 
          matrix[ j ][ i ] = matrix[ i ][ j ];

      return numBaseFunctions;
    }

    void boundaryCorrectOnGrid () const
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        if( it->hasBoundaryIntersections() )
          boundaryCorrectOnEntity( *it );
      }
    }

    /*! treatment of Dirichlet-DoFs for one entity
     *
     *   delete rows for dirichlet-DoFs, setting diagonal element to 1.
     *
     *   \note A LagrangeDiscreteFunctionSpace is implicitly assumed.
     *
     *   \param[in]  entity  entity to perform Dirichlet treatment on
     */
    template< class EntityType >
    void boundaryCorrectOnEntity ( const EntityType &entity ) const
    {
      typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
        LagrangePointSetType;

      enum { faceCodim = 1 };
      typedef typename GridPartType :: IntersectionIteratorType
        IntersectionIteratorType;
      typedef typename LagrangePointSetType :: template Codim< faceCodim >
                                            :: SubEntityIteratorType
        FaceDofIteratorType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      const GridPartType &gridPart = dfSpace.gridPart();

      const LagrangePointSetType &lagrangePointSet
        = dfSpace.lagrangePointSet( entity );
 
      IntersectionIteratorType it = gridPart.ibegin( entity );
      const IntersectionIteratorType endit = gridPart.iend( entity );
      for( ; it != endit ; ++it ) {
        if( !it.boundary() )
          continue;
       
        const int face = it.numberInSelf();
        FaceDofIteratorType faceIt
          = lagrangePointSet.template beginSubEntity< faceCodim >( face );
        const FaceDofIteratorType faceEndIt
          = lagrangePointSet.template endSubEntity< faceCodim >( face );
        for( ; faceIt != faceEndIt; ++faceIt ) {
          const unsigned int dof
            = dfSpace.mapToGlobal( entity, *faceIt );
          matrix_->kroneckerKill( dof, dof );
        }
      }
    }
  };



  // Assembler for right rand side
  // note: this is not an L2-Projection
  template< class DiscreteFunctionImp >
  class RightHandSideAssembler
  {
  public:
    typedef DiscreteFunctionImp DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType :: LocalFunctionType
      LocalFunctionType;
  
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;

    enum { dimension = GridType :: dimension };
  
  public:
    // discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    static void assemble( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;
      //typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;
      
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();
  
      discreteFunction.clear(); //discreteFunction auf Null setzen

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        const GeometryType &geometry = entity.geometry(); //Referenz auf Geometrie
      
        LocalFunctionType localFunction = discreteFunction.localFunction( entity ); 
        const BaseFunctionSetType baseFunctionSet //BaseFunctions leben immer auf Refernzelement!!!
          = discreteFunctionSpace.baseFunctionSet( entity ); 

        CachingQuadrature< GridPartType, 0 > quadrature( entity, polOrd ); //0 --> codim 0
        const int numQuadraturePoints = quadrature.nop();
        for( int qP = 0; qP < numQuadraturePoints; ++qP )
        {
          const double det
            = geometry.integrationElement( quadrature.point( qP ) );
          const double factor = det * quadrature.weight( qP );

          RangeType phi;
          function.evaluate( geometry.global( quadrature.point( qP ) ), phi );
            
          const int numDofs = localFunction.numDofs(); //Dofs = Freiheitsgrade (also die Unbekannten)
          for( int i = 0; i < numDofs; ++i )
          {
            RangeType psi; //R"uckgabe-Funktionswerte
        
            baseFunctionSet.evaluate( i, quadrature, qP, psi ); //i = i'te Basisfunktion; qP Quadraturpunkt
            localFunction[ i ] += factor * (phi * psi);
          }
        }
      }
    }
  };

} // end namespace 

#endif

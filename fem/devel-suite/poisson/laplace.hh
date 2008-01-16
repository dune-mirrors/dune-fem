#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

#ifdef ENABLE_TIMING
#include <time.h>
#endif

//- Dune includes
#include <dune/common/fmatrix.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  //! \brief The Laplace operator
  template< class DiscreteFunction, class MatrixObject, class Tensor >
  class LaplaceFEOp
  : public Operator< typename DiscreteFunction :: RangeFieldType,
                     typename DiscreteFunction :: RangeFieldType,
                     DiscreteFunction,
                     DiscreteFunction >
  {
  public:
    //! type of discrete functions
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of the system matrix object
    typedef MatrixObject MatrixObjectType;
    
    //! type of tensor
    typedef Tensor TensorType;

    //! type of this LaplaceFEOp
    typedef LaplaceFEOp< DiscreteFunctionType, MatrixObjectType, TensorType >
      LaplaceFEOpType;

  private:
    typedef LaplaceFEOpType ThisType;

  public:
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
        
    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 > QuadratureType;

    typedef typename MatrixObjectType :: LocalMatrixType LocalMatrixType;
    
  private:
    enum { maxBaseFunctions = 100 };
    
  protected:
    class Assembler;
   
  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

    //! pointer to the system matrix 
    mutable MatrixObjectType matrixObject_;
 
    //! flag indicating whether the system matrix has been assembled
    mutable bool matrix_assembled_;
      
    TensorType *const stiffTensor_;

  private:
    mutable JacobianRangeType grad;

    mutable JacobianRangeType mygrad[ maxBaseFunctions ];
 
  public:
    //! constructor
    inline explicit LaplaceFEOp( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
      matrix_assembled_( false ),
      stiffTensor_( 0 )
    {
    }
        
    //! constructor
    inline LaplaceFEOp( TensorType &stiffTensor,
                        const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
      matrix_assembled_( false ),
      stiffTensor_( &stiffTensor )
    { 
    }

  private:
    // prohibit copying
    inline LaplaceFEOp ( const ThisType & );

  public:
    //! \brief apply the operator
    virtual void operator() ( const DiscreteFunctionType &u, 
                              DiscreteFunctionType &w ) const 
    {
      //systemMatrix().apply( u, w );
      systemMatrix().multOEM( u.leakPointer(), w.leakPointer() );
    }
  
    /*! \brief obtain a reference to the system matrix 
     *
     *  The assembled matrix is returned. If the system matrix has not been
     *  assembled, yet, the assembly is performed.
     *
     *  \returns a reference to the system matrix
     */
    MatrixObjectType &systemMatrix () const
    {
      if( !matrix_assembled_ )
        assemble();
      return matrixObject_;
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

      LagrangeMatrixSetup stencil;
      matrixObject_.reserve( stencil );

#if 0
      // if the matrix has not been allocated, allocate it.
      if( this->matrix_ == 0 )
      {
        const int size = dfSpace.size();
        const int numNonZero = 8 * (1 << dimension) * polynomialOrder;
        matrix_ = new MatrixType( size, size, numNonZero );
        assert( matrix_ != 0 );
      }
#endif

#ifdef ENABLE_TIMING
      time_t starttime = time( NULL );
#endif
      
      matrixObject_.clear();
      Assembler a( *this );
      dfSpace.forEach( a );
      //assembleOnGrid();
      boundaryCorrectOnGrid();

#ifdef ENABLE_TIMING
      time_t endtime = time( NULL );
      std :: cout << "Time to assemble matrix: " << (endtime - starttime) << std :: endl;
#endif

      matrix_assembled_ = true;
    }

  protected:
    /*! assemble the local matrix
     *
     *  In the Finite Element Method, most base functions are zero on a given
     *  entity. To build the matrix on one entity, we only need to store a very
     *  limited amount of matrix entries, the socalled local matrix.
     */
    template< class EntityType, class LocalMatrixType >
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
          baseSet.jacobian( i, quad[ pt ], mygrad[ i ] ); 
      
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

        LocalMatrixType localMatrix = matrixObject_.localMatrix( entity, entity );
       
        const int face = it.numberInSelf();
        FaceDofIteratorType faceIt
          = lagrangePointSet.template beginSubEntity< faceCodim >( face );
        const FaceDofIteratorType faceEndIt
          = lagrangePointSet.template endSubEntity< faceCodim >( face );
        for( ; faceIt != faceEndIt; ++faceIt )
          localMatrix.unitRow( *faceIt );
#if 0
        {
          const unsigned int dof = dfSpace.mapToGlobal( entity, *faceIt );
          //matrix_->kroneckerKill( dof, dof );
          matrix_->unitRow( dof );
        }
#endif
      }
    }
  };



  template< class DiscreteFunction, class MatrixObject, class Tensor >
  class LaplaceFEOp< DiscreteFunction, MatrixObject, Tensor > :: Assembler
  {
  protected:
    const LaplaceFEOpType &feop_;

    mutable FieldMatrix< RangeFieldType, maxBaseFunctions, maxBaseFunctions >
      elementMatrix_;
#if 0
    mutable StandardArray< unsigned int, maxBaseFunctions >
      localMap_;
#endif

  public:
    inline Assembler( const LaplaceFEOpType &feop )
    : feop_( feop )
    {
    }
    
    template< class EntityType >
    inline void operator() ( const EntityType &entity ) const
    {
      //const DiscreteFunctionSpaceType &dfSpace = feop_.discreteFunctionSpace();
    
      // obtain local matrix
      const unsigned int size
        = feop_.assembleLocalMatrix( entity, elementMatrix_ );

#if 0
      for( unsigned int i = 0; i < size; ++i )
        localMap_[ i ] = dfSpace.mapToGlobal( entity, i );
#endif
      LocalMatrixType localMatrix = feop_.matrixObject_.localMatrix( entity, entity );
      
      for( unsigned int i = 0; i < size; ++i )
      { 
        for( unsigned int j = 0; j < size; ++j )
          localMatrix.add( i, j, elementMatrix_[ i ][ j ] );
//          feop_.matrix_->add( localMap_[ i ], localMap_[ j ], localMatrix_[ i ][ j ] );
      }
    }
  };




  // Assembler for right rand side
  // note: this is not an L2-Projection
  template< class DiscreteFunction >
  class RightHandSideAssembler
  {
  public:
    typedef DiscreteFunction DiscreteFunctionType;
    
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

    typedef CommunicationManager< DiscreteFunctionSpaceType >
      CommunicationManagerType;
  
  public:
    // discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    static void assemble( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;

      // We use a caching quadrature for codimension 0 entities
      typedef CachingQuadrature< GridPartType, 0 > QuadratureType;
      
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();
  
      discreteFunction.clear(); //discreteFunction auf Null setzen

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        // *it gives a reference to the current entity
        const EntityType &entity = *it;

        // obtain a reference to the entity's geometry
        const GeometryType &geometry = entity.geometry();
      
        // obtain BaseFunctionSet for the entity
        // note that base functions are always defined on the reference geometry
        LocalFunctionType localFunction = discreteFunction.localFunction( entity ); 
        const BaseFunctionSetType baseFunctionSet
          = discreteFunctionSpace.baseFunctionSet( entity ); 

        // obtain number of DoFs (degrees of freedom, the unknowns)
        const int numDofs = localFunction.numDofs();

        QuadratureType quadrature( entity, polOrd );
        const int numQuadraturePoints = quadrature.nop();
        for( int qP = 0; qP < numQuadraturePoints; ++qP )
        {
          const double det
            = geometry.integrationElement( quadrature.point( qP ) );
          const double factor = det * quadrature.weight( qP );

          RangeType phi;
          function.evaluate( geometry.global( quadrature.point( qP ) ), phi );
            
          for( int i = 0; i < numDofs; ++i )
          {
            RangeType psi;
        
            // evaluate the i-th base function in the quadrature point qp
            // the result is stored in psi
            baseFunctionSet.evaluate( i, quadrature[ qP ], psi );
            localFunction[ i ] += factor * (phi * psi);
          }
        }
      }

      CommunicationManagerType communicate( discreteFunctionSpace );
      communicate.exchange( discreteFunction );
    }
  };

} // end namespace 

#endif

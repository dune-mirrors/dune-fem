#ifndef DUNE_FEM_P12DSPACE_TEST_LAPLACE_HH
#define DUNE_FEM_P12DSPACE_TEST_LAPLACE_HH

//- Dune includes
#include <dune/common/fmatrix.hh>
#include <dune/common/timer.hh>

#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

namespace Dune
{

  //! \brief The Laplace operator
  template< class DiscreteFunction, class MatrixTraits, class Tensor >
  class LaplaceFEOp
  : public Operator< typename DiscreteFunction :: RangeFieldType,
                     typename DiscreteFunction :: RangeFieldType,
                     DiscreteFunction,
                     DiscreteFunction >,
    public OEMSolver::PreconditionInterface
  {
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
  public:
    //! type of discrete functions
    typedef DiscreteFunction                                         DiscreteFunctionType;

    //! type of tensor
    typedef Tensor                                                   TensorType;

    //! type of this LaplaceFEOp
    typedef LaplaceFEOp< DiscreteFunctionType, MatrixTraits,
                         TensorType >                                LaplaceFEOpType;

  private:
    typedef LaplaceFEOpType                                          ThisType;

  public:
    //! type of discrete function space
    typedef typename DiscreteFunctionType
              :: DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;
    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType     RangeFieldType;
       
  protected:
    //! type of jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType  JacobianRangeType;
    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType
              :: BaseFunctionSetType                                 BaseFunctionSetType;

  public:
    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    //! type of grid
    typedef typename DiscreteFunctionSpaceType :: GridType           GridType;

    //! polynomial order of base functions
    enum { polynomialOrder = DiscreteFunctionSpaceType :: polynomialOrder };

    //! The grid's dimension
    enum { dimension = GridType :: dimension };
        
    //! type of quadrature to be used
    typedef CachingQuadrature< GridPartType, 0 >                     QuadratureType;

    typedef typename MatrixTraits
              :: template MatrixObject
                    < LagrangeMatrixTraits< MatrixTraits > >
                :: MatrixObjectType                                  MatrixObjectType;

    typedef typename MatrixObjectType :: LocalMatrixType             LocalMatrixType;
    typedef typename MatrixObjectType :: PreconditionMatrixType      PreconditionMatrixType;
    typedef typename MatrixObjectType :: MatrixType                  MatrixType;

  public:
    // types for boundary treatment
    // ----------------------------
    typedef typename DiscreteFunctionSpaceType :: MapperType         MapperType;
/*    typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType >       SlaveDofsType;
 *    typedef typename SlaveDofsType :: SingletonKey                   SlaveDofsKeyType;
 *    typedef SingletonList< SlaveDofsKeyType, SlaveDofsType >         SlaveDofsProviderType;*/
    
  protected:
    class Assembler;

    typedef DofManager< GridType >                                   DofManagerType;
    typedef DofManagerFactory< DofManagerType >                      DofManagerFactoryType;
   
  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;
/*    const DofManagerType &dofManager_;*/

    //! pointer to the system matrix
    mutable MatrixObjectType matrixObject_;
 
    //! flag indicating whether the system matrix has been assembled
/*    mutable int sequence_;*/
      
    TensorType *const stiffTensor_;

/*    SlaveDofsType *const slaveDofs_;*/

  private:
    mutable JacobianRangeType grad;

  public:
    //! constructor
    explicit LaplaceFEOp( const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
/*      sequence_( -1 ),*/
      stiffTensor_( 0 )//,
/*      slaveDofs_( getSlaveDofs( discreteFunctionSpace_ ) )*/
    { }
        
    //! constructor
    LaplaceFEOp( TensorType &stiffTensor, const DiscreteFunctionSpaceType &dfSpace )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
/*      sequence_( -1 ),*/
      stiffTensor_( &stiffTensor )//,
/*      slaveDofs_( getSlaveDofs( discreteFunctionSpace_ ) )*/
    { }

  private:
    // prohibit copying
    LaplaceFEOp ( const ThisType & );

  public:
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
    MatrixObjectType &systemMatrix () const
    {
      // if stored sequence number it not equal to the one of the
      // dofManager (or space) then the grid has changed 
      // and matrix has to be assembled new 
/*      if( sequence_ != dofManager_.sequence() )*/
/*        assemble();*/

      return matrixObject_;
    }

    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType &preconditionMatrix () const
    {
      return systemMatrix().preconditionMatrix();
    }

    //! return true if preconditioning is enabled
    bool hasPreconditionMatrix () const
    {
      return matrixObject_.hasPreconditionMatrix();
    }

    //! print the system matrix into a stream
    void print ( std :: ostream & out = std :: cout ) const 
    {
      systemMatrix().matrix().print( out );
    }

    const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    /** \brief perform a grid walkthrough and assemble the global matrix */
    void assemble () const 
    {
      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      // reserve memory for matrix 
      matrixObject_.reserve();

      // create timer (also stops time)
      Timer timer;

      // clear matrix 
      matrixObject_.clear();

      // create local matrix assembler
      Assembler a( *this );

      // apply local matrix assembler on each element 
      dfSpace.forEach( a );

      // correct matrix due to boundary conditions 
      boundaryCorrectOnGrid();

      // get elapsed time 
      const double assemblyTime = timer.elapsed();
      // in verbose mode print times 
      std :: cout << "Time to assemble matrix: " << assemblyTime << "s" << std :: endl;

      // get grid sequence number from space (for adaptive runs)
/*      sequence_ = dofManager_.sequence();*/
    }

  protected:
    // evaluate of stiffness tensor 
    template< class Point >
    RangeFieldType stiffTensor ( const Point &x ) const
    {
      if( stiffTensor_ != 0 )
      {
        typename DiscreteFunctionSpaceType :: RangeType phi;
        stiffTensor_->evaluate( x, phi );
        return phi[ 0 ];
      }
      else
        return RangeFieldType( 1 );
    }

    //! adjust matrix due to boundary conditions
    void boundaryCorrectOnGrid () const
    {
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      typedef typename IteratorType::Entity EntityType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      const IteratorType end = dfSpace.end();
      for( IteratorType it = dfSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        // if entity has boundary intersections 
        if( entity.hasBoundaryIntersections() )
          boundaryCorrectOnEntity( entity );
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
      typedef typename DiscreteFunctionSpaceType::LocalFiniteElementType LocalFiniteElementType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();
      const LocalFiniteElementType &fem        = dfSpace.localFiniteElement();

      // get local matrix from matrix object 
      LocalMatrixType localMatrix = matrixObject_.localMatrix( entity, entity );

      // This is a hack: it returns a vector of ones and zeros where a one
      // marks a boundary dof
      std::vector< typename DiscreteFunctionSpaceType::RangeType > dofOut;
      BoundaryCheck< EntityType > bc( entity );
      fem.localInterpolation().interpolate( bc, dofOut );

      for( size_t i = 0; i < dofOut.size(); ++i )
      {
        if( dofOut[ i ][ 0 ] == 1 )
        {
          // old version of bnd treatment
          // localMatrix.unitRow( i );
 
          localMatrix.clearRow(i);
          localMatrix.set(i, i, 1.);
        }
      }
    }
  };

  //! Laplace Local Assembler 
  template< class DiscreteFunction, class MatrixObject, class Tensor >
  class LaplaceFEOp< DiscreteFunction, MatrixObject, Tensor > :: Assembler
  {
  protected:
    const LaplaceFEOpType &feop_;

    mutable DynamicArray< JacobianRangeType > gradCache_;
    mutable RangeFieldType weight_;
 
  public:
    explicit Assembler ( const LaplaceFEOpType &feop )
    : feop_( feop ),
      gradCache_( feop_.discreteFunctionSpace().mapper().maxNumDofs() )
    {}

    //! add scalar product of cached gradients to local matrix
    void updateLocalMatrix ( LocalMatrixType &localMatrix ) const
    {
      const unsigned int columns = localMatrix.columns();
      for( unsigned int i = 0; i < columns; ++i )
      {
        localMatrix.add( i, i, weight_ *
                         (gradCache_[ i ][ 0 ] * gradCache_[ i ][ 0 ]  ) );
        for ( unsigned int j = 0; j < i; ++j )
        {
          const RangeFieldType value
            = weight_ * (gradCache_[ i ][ 0 ] * gradCache_[ j ][ 0 ]);
          localMatrix.add( j, i, value );
          localMatrix.add( i, j, value );
        }
      }
    }

    //! assemble local matrix for given entity 
    template< class EntityType >
    void operator() ( const EntityType &entity ) const
    {
      typedef typename EntityType :: Geometry                        GeometryType;

      assert( entity.partitionType() != GhostEntity );

      const GeometryType &geometry = entity.geometry();

      // get local matrix from matrix object
      LocalMatrixType localMatrix
        = feop_.matrixObject_.localMatrix( entity, entity );
      
      // get base function set 
      const BaseFunctionSetType &baseSet = localMatrix.domainBaseFunctionSet();

      // get number of local base functions 
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
            
      // create quadrature 
      QuadratureType quadrature( entity, 2 * (polynomialOrder - 1) );

      // loop over all quadrature points 
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        // get local coordinate of quadrature point 
        const typename QuadratureType :: CoordinateType &x
          = quadrature.point( pt );

        // get jacobian inverse transposed 
        const FieldMatrix< double, dimension, dimension > &inv
          = geometry.jacobianInverseTransposed( x );

        // for all base functions evaluate the graddient 
        // on quadrature point pt and apply jacobian inverse 
        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          baseSet.jacobian( i, quadrature[ pt ], gradCache_[ i ] );
          gradCache_[ i ][ 0 ] = FMatrixHelp :: mult( inv, gradCache_[ i ][ 0 ] );
        }
        
        // evaluate integration weight 
        weight_ = quadrature.weight( pt ) * geometry.integrationElement( x )
                  * feop_.stiffTensor( geometry.global( x ) );
        
        // add scalar product of gradients to local matrix 
        updateLocalMatrix( localMatrix );
      }
    }
  };

  // Assembler for right rand side
  // note: this is not an L2-Projection
  template< class DiscreteFunction >
  class RightHandSideAssembler
  {
  public:
    typedef DiscreteFunction                                         DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType
              :: DiscreteFunctionSpaceType                           DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType :: LocalFunctionType       LocalFunctionType;
  
    typedef typename DiscreteFunctionSpaceType
              :: BaseFunctionSetType                                 BaseFunctionSetType;
    typedef typename DiscreteFunctionSpaceType :: RangeType          RangeType;
    typedef typename DiscreteFunctionSpaceType :: GridPartType       GridPartType;
    typedef typename GridPartType :: GridType                        GridType;

    enum { dimension = GridType :: dimension };

  public:
    //! discreteFunction is an output parameter (kind of return value)
    template< int polOrd, class FunctionType >
    static void assemble( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction )
    {
      typedef typename DiscreteFunctionSpaceType :: IteratorType     IteratorType;
      typedef typename IteratorType :: Entity                        EntityType;
      typedef typename EntityType :: Geometry                        GeometryType;

      // We use a caching quadrature for codimension 0 entities
      typedef CachingQuadrature< GridPartType, 0 >                   QuadratureType;
      
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();
  
      // set discreteFunction to zero
      discreteFunction.clear(); 

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
          // evaluate right hand side function
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

      // communicate data (for parallel runs)
      discreteFunctionSpace.communicate( discreteFunction );
    }
  };

} // end namespace 
#endif //DUNE_FEM_P12DSPACE_TEST_LAPLACE_HH


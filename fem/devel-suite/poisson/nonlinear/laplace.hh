#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

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
  template< class DiscreteFunction, class MatrixTraits, class Tensor, class MassFunction >
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
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of tensor
    typedef Tensor TensorType;

    //! type of mass function
    typedef MassFunction MassFunctionType;

    //! type of this LaplaceFEOp
    typedef LaplaceFEOp< DiscreteFunctionType, MatrixTraits, TensorType, MassFunctionType >
      LaplaceFEOpType;

  private:
    typedef LaplaceFEOpType ThisType;

  public:
    //! type of discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! field type of range
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType
      RangeFieldType;
       
  protected:
    //! type of jacobian
    typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
      JacobianRangeType;
    //! type of the base function set
    typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType
      BaseFunctionSetType;

  public:
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

    typedef typename MatrixTraits
      :: template  MatrixObject< LagrangeMatrixTraits< MatrixTraits > >
      :: MatrixObjectType MatrixObjectType;

    typedef typename MatrixObjectType :: LocalMatrixType LocalMatrixType;
    typedef typename MatrixObjectType :: PreconditionMatrixType PreconditionMatrixType;
    typedef typename MatrixObjectType :: MatrixType MatrixType;

    // types for boundary treatment
    // ----------------------------
    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType > SlaveDofsType;
    typedef typename SlaveDofsType :: SingletonKey SlaveDofsKeyType;
    typedef SingletonList< SlaveDofsKeyType, SlaveDofsType >
      SlaveDofsProviderType;
    
  protected:
    class Assembler;
   
  protected:
    const DiscreteFunctionSpaceType &discreteFunctionSpace_;

    //! pointer to the system matrix 
    mutable MatrixObjectType matrixObject_;
 
    //! flag indicating whether the system matrix has been assembled
    mutable bool matrix_assembled_;
      
    TensorType *const stiffTensor_;

    SlaveDofsType *const slaveDofs_;

    mutable MassFunctionType  massFunction_;
   
#ifdef ENABLE_TIMING
    Timer timer_;
#endif

  private:
    mutable JacobianRangeType grad;

  public:
    //! constructor
    inline explicit LaplaceFEOp( const DiscreteFunctionSpaceType &dfSpace, const MassFunctionType &massFunction )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
      matrix_assembled_( false ),
      stiffTensor_( 0 ),
      slaveDofs_( getSlaveDofs( discreteFunctionSpace_ ) ),
      massFunction_( massFunction )
    {
    }
        
    //! constructor
    inline LaplaceFEOp( TensorType &stiffTensor,
                        const DiscreteFunctionSpaceType &dfSpace, const MassFunctionType &massFunction  )
    : discreteFunctionSpace_( dfSpace ),
      matrixObject_( discreteFunctionSpace_, discreteFunctionSpace_ ),
      matrix_assembled_( false ),
      stiffTensor_( &stiffTensor ),
      slaveDofs_( getSlaveDofs( discreteFunctionSpace_ ) ),
      massFunction_( massFunction )
    { 
    }

  private:
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    // prohibit copying
    inline LaplaceFEOp ( const ThisType & );

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
      if( !matrix_assembled_ )
        assemble();
      return matrixObject_;
    }


    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType &preconditionMatrix () const
    {
      return systemMatrix().preconditionMatrix();
    }

    //! return true if preconditioning is enabled
    bool hasPreconditionMatrix () const { return matrixObject_.hasPreconditionMatrix(); }


    //! print the system matrix into a stream
    void print ( std :: ostream out = std :: cout ) const 
    {
      systemMatrix().print( out );
    }

    const DiscreteFunctionSpaceType &discreteFunctionSpace () const
    {
      return discreteFunctionSpace_;
    }

    /** \brief perform a grid walkthrough and assemble the global matrix */
    void assemble () const 
    {
      const DiscreteFunctionSpaceType &dfSpace = discreteFunctionSpace();

      matrixObject_.reserve();

#ifdef ENABLE_TIMING
      timer_.reset();
#endif
      
      matrixObject_.clear();
      Assembler a( *this, massFunction_ );
      dfSpace.forEach( a );

      boundaryCorrectOnGrid();

#ifdef ENABLE_TIMING
      std :: cout << "Time to assemble matrix: " << timer_.elapsed() << std :: endl;
#endif

      matrix_assembled_ = true;
    }

  protected:
    template< class Point >
    inline RangeFieldType stiffTensor ( const Point &x ) const
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

      SlaveDofsType &slaveDofs = this->slaveDofs();
      int numSlaveDofs = slaveDofs.size();
 
      IntersectionIteratorType it = gridPart.ibegin( entity );
      const IntersectionIteratorType endit = gridPart.iend( entity );
      for( ; it != endit ; ++it )
      {
        if( !it->boundary() )
          continue;

        LocalMatrixType localMatrix = matrixObject_.localMatrix( entity, entity );
       
        const int face = it->numberInSelf();
        FaceDofIteratorType faceIt
          = lagrangePointSet.template beginSubEntity< faceCodim >( face );
        const FaceDofIteratorType faceEndIt
          = lagrangePointSet.template endSubEntity< faceCodim >( face );
        for( ; faceIt != faceEndIt; ++faceIt )
        {
          const int localDof = *faceIt;
          localMatrix.unitRow( localDof );

          const int globalDof = dfSpace.mapToGlobal( entity, localDof );
          for( int i = 0; i < numSlaveDofs; ++i )
          {
            if( globalDof == slaveDofs[ i ] )
              localMatrix.set( localDof, localDof, 0 );
          }
        }
      }
    }

  protected:
    inline static SlaveDofsType *getSlaveDofs ( const DiscreteFunctionSpaceType &space )
    {
      SlaveDofsKeyType key( space, space.mapper() );
      return &(SlaveDofsProviderType :: getObject( key ));
    }

    inline SlaveDofsType &slaveDofs () const
    {
      slaveDofs_->rebuild();
      return *slaveDofs_;
    }
  };



  template< class DiscreteFunction, class MatrixObject, class Tensor, class MassFunction >
  class LaplaceFEOp< DiscreteFunction, MatrixObject, Tensor, MassFunction > :: Assembler
  {
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  protected:
    const LaplaceFEOpType &feop_;
    const MassFunctionType &massFunction_;
    mutable DynamicArray< JacobianRangeType > gradCache_;
    mutable RangeFieldType weight_;
 
  public:
    inline explicit Assembler ( const LaplaceFEOpType &feop, const MassFunctionType &massFunction )
    : feop_( feop ),
      massFunction_( massFunction),
      gradCache_( feop_.discreteFunctionSpace().mapper().maxNumDofs() )
    {}

    inline void updateLocalMatrix ( LocalMatrixType &localMatrix, RangeFieldType localweight_) const
    {
      const unsigned int columns = localMatrix.columns();
      for( unsigned int i = 0; i < columns; ++i )
      {
        localMatrix.add( i, i, localweight_ * (gradCache_[ i ][ 0 ] * gradCache_[ i ][ 0 ]  ) );
        for ( unsigned int j = 0; j < i; ++j )
        {
          const RangeFieldType value
            = localweight_ * (gradCache_[ i ][ 0 ] * gradCache_[ j ][ 0 ]);
          localMatrix.add( j, i, value );
          localMatrix.add( i, j, value );
        }
      }
    } 

template< class EntityType >
inline void assemlbeMassMatrix(const EntityType &entity, LocalMatrixType &localMatrix) const
 {
   const typename EntityType :: Geometry &geometry = entity.geometry();
   const BaseFunctionSetType &baseSet = localMatrix.domainBaseFunctionSet();
   QuadratureType massQuadrature( entity, 3 * polynomialOrder  );
   const unsigned int massNumQuadraturePoints = massQuadrature.nop(); 
   RangeType phi, psi,psi2;
   

      for( unsigned int qP = 0; qP < massNumQuadraturePoints; ++qP )
      {
        weight_ = massQuadrature.weight( qP ) * geometry.integrationElement( massQuadrature.point( qP ) );

	massFunction_.evaluate(geometry.global( massQuadrature.point( qP )), phi); 

        const unsigned int columns = localMatrix.columns();
        for( unsigned int i = 0; i < columns; ++i )
        {
          baseSet.evaluate( i, massQuadrature[ qP ], psi );
          localMatrix.add( i, i, weight_ * (psi * psi  )*phi);

          for ( unsigned int j = 0; j < i; ++j )
          {
            baseSet.evaluate( j, massQuadrature[ qP ], psi2 );
            const RangeFieldType value = weight_ * (psi * psi2)* phi;
            localMatrix.add( j, i, value );
            localMatrix.add( i, j, value );
           }
         }
       }
}

template< class EntityType >
inline void assemlbeStiffMatrix(const EntityType &entity, LocalMatrixType &localMatrix) const
{
      const typename EntityType :: Geometry &geometry = entity.geometry();
      const BaseFunctionSetType &baseSet = localMatrix.domainBaseFunctionSet();
      const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
            
      QuadratureType quadrature( entity, 2 * (polynomialOrder - 1) );
      const unsigned int numQuadraturePoints = quadrature.nop();
      for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
        const typename QuadratureType :: CoordinateType &x
          = quadrature.point( pt );

        const FieldMatrix< double, dimension, dimension > &inv
          = geometry.jacobianInverseTransposed( x );
        for( unsigned int i = 0; i < numBaseFunctions; ++i )
        {
          baseSet.jacobian( i, quadrature[ pt ], gradCache_[ i ] );
          gradCache_[ i ][ 0 ] = FMatrixHelp :: mult( inv, gradCache_[ i ][ 0 ] );
        }
//              weight_ = quadrature.weight( pt ) * geometry.integrationElement( x )
//                   * feop_.stiffTensor( geometry.global( x ) );
//         updateLocalMatrix( localMatrix );
        weight_ = quadrature.weight( pt ) * geometry.integrationElement( x );
        RangeFieldType localweight_ = weight_* feop_.stiffTensor( geometry.global( x ) );
        updateLocalMatrix( localMatrix, localweight_ );
     }

}

    template< class EntityType >
    inline void operator() ( const EntityType &entity ) const
    {
      assert( entity.partitionType() == InteriorEntity );

//       const typename EntityType :: Geometry &geometry = entity.geometry();

       LocalMatrixType localMatrix
         = feop_.matrixObject_.localMatrix( entity, entity );

      
//       const BaseFunctionSetType &baseSet = localMatrix.domainBaseFunctionSet();
//       const unsigned int numBaseFunctions = baseSet.numBaseFunctions();
//             
//       QuadratureType quadrature( entity, 2 * (polynomialOrder - 1) );
//       const unsigned int numQuadraturePoints = quadrature.nop();
//       for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
//       {
//         const typename QuadratureType :: CoordinateType &x
//           = quadrature.point( pt );
// 
//         const FieldMatrix< double, dimension, dimension > &inv
//           = geometry.jacobianInverseTransposed( x );
//         for( unsigned int i = 0; i < numBaseFunctions; ++i )
//         {
//           baseSet.jacobian( i, quadrature[ pt ], gradCache_[ i ] );
//           gradCache_[ i ][ 0 ] = FMatrixHelp :: mult( inv, gradCache_[ i ][ 0 ] );
//         }
// //              weight_ = quadrature.weight( pt ) * geometry.integrationElement( x )
// //                   * feop_.stiffTensor( geometry.global( x ) );
// //         updateLocalMatrix( localMatrix );
//         weight_ = quadrature.weight( pt ) * geometry.integrationElement( x );
//         RangeFieldType localweight_ = weight_* feop_.stiffTensor( geometry.global( x ) );
//         updateLocalMatrix( localMatrix, localweight_ );
//      }
assemlbeStiffMatrix(entity, localMatrix);
assemlbeMassMatrix(entity, localMatrix);
//     QuadratureType massQuadrature( entity, 3 * polynomialOrder  );
//     const unsigned int massNumQuadraturePoints = massQuadrature.nop(); 
//     RangeType phi, psi,psi2;
//  phi= massFunction_;
// 
//       for( unsigned int qP = 0; qP < massNumQuadraturePoints; ++qP )
//       {
// //        const typename QuadratureType :: CoordinateType &x
// //           = quadrature.point( qP );   
// 
//         const unsigned int columns = localMatrix.columns();
//       for( unsigned int i = 0; i < columns; ++i )
//       {
//         baseSet.evaluate( i, massQuadrature[ qP ], psi );
//         localMatrix.add( i, i, weight_ * (psi * psi  )*phi);
//         for ( unsigned int j = 0; j < i; ++j )
//         {
//           baseSet.evaluate( j, massQuadrature[ qP ], psi2 );
//           const RangeFieldType value
//             = weight_ * (psi * psi2)* phi;
//           localMatrix.add( j, i, value );
//           localMatrix.add( i, j, value );
//         }
//       }
//       }


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

      // communicate data 
      discreteFunctionSpace.communicate( discreteFunction );
    }
  }; 

} // end namespace 

#endif

#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

//- Dune includes
#include <dune/fem/quadrature/quadrature.hh>

//- local includes
#include "feop.hh"

namespace Dune 
{

  /** \brief The Laplace operator
   */
  template< class DiscreteFunctionImp, class TensorImp >
  class LaplaceFEOp
  : public FEOp< DiscreteFunctionImp,
                 SparseRowMatrix< typename DiscreteFunctionImp :: RangeFieldType >,
                 LaplaceFEOp< DiscreteFunctionImp, TensorImp > >
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
    typedef FEOp< DiscreteFunctionType, MatrixType, ThisType > BaseType;

  public:
    //! Operation mode for Finite Element Operator
    typedef typename BaseType :: OpMode OpMode;
 
  private:
    mutable JacobianRangeType grad;
    mutable JacobianRangeType othGrad;

    enum { maxnumOfBaseFct = 100 };
    mutable JacobianRangeType mygrad[ maxnumOfBaseFct ];
        
  private:
    //DiscreteFunctionType *stiffFunction_;
    TensorType *stiffTensor_;
       
  public:
    //! constructor
    LaplaceFEOp( const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode )
    : BaseType( discreteFunctionSpace, opMode ),
      //stiffFunction_( NULL ),
      stiffTensor_( NULL )
    {
    }
        
        #if 0
        //! constructor
        LaplaceFEOp( const DiscreteFunctionType &stiff,
                     const DiscreteFunctionSpaceType &discreteFunctionSpace,
                     OpMode opMode )
          : BaseType( discreteFunctionSpace, opMode ),
            stiffFunction_( &stiff ),
            stiffTensor_( NULL )
        {
        }
        #endif
        
    //! constructor
    LaplaceFEOp( TensorType &stiff,
                 const DiscreteFunctionSpaceType &discreteFunctionSpace,
                 OpMode opMode )
    : BaseType( discreteFunctionSpace, opMode ),
      //stiffFunction_( NULL ),
      stiffTensor_( &stiff )
    { 
    }
        
    //! Returns the actual matrix if it is assembled
    const MatrixType* getMatrix () const
    {
      assert( this->matrix_ );
      return this->matrix_;
    }
        
    //! Creates a new empty matrix
    MatrixType* newEmptyMatrix () const 
    {
      return new MatrixType( this->functionSpace_.size(),
                             this->functionSpace_.size(), 
                             15 * dimension );
    }
        
    //! Prepares the local operator before calling apply()
    void prepareGlobal ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest )
    {
      this->arg_  = &arg;
      this->dest_ = &dest;
      // assert( this->arg_ != 0 );
      // assert(this->dest_ != 0);
      this->dest_.clear();
    }
        
    //! return the matrix entr that belong to basis function i and j 
    template< class EntityType >
    double getLocalMatrixEntry( const EntityType &entity,
                                const int i,
                                const int j ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseSet
        = discreteFunctionSpace.baseFunctionSet( entity );

      double val = 0;
      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) { 
        baseSet.jacobian( i, quad, pt, grad );

        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimension, dimension > &inv
          = geometry.jacobianInverseTransposed( quad.point( pt ) );
        const double vol = geometry.integrationElement( quad.point( pt ) );
            
        // multiply with transpose of jacobian inverse 
        grad[ 0 ] = FMatrixHelp :: mult( inv, grad[ 0 ] );
                    
        if( i != j ) {
          baseSet.jacobian( j, quad, pt, othGrad );
                            
          // multiply with transpose of jacobian inverse 
          othGrad[ 0 ] = FMatrixHelp :: multTransposed( inv, othGrad[ 0 ] );
          val += (grad[ 0 ] * othGrad[ 0 ] ) * quad.weight( pt ) * vol;
        } else
          val += ( grad[ 0 ] * grad[ 0 ] ) * quad.weight( pt ) * vol;
      }
      return val;
    }
        
    //! ???
    template< class  EntityType, class LocalMatrixType >
    void getLocalMatrix( const EntityType &entity,
                         const int matrixSize,
                         LocalMatrixType &matrix ) const
    {
      typedef typename EntityType :: Geometry GeometryType;

      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = this->functionSpace_;
      const GeometryType &geometry = entity.geometry();

      const BaseFunctionSetType &baseSet
        = discreteFunctionSpace.baseFunctionSet( entity );
            
      assert( matrixSize <= maxnumOfBaseFct );
      for( int i = 0; i < matrixSize; ++i )
        for( int j = 0; j <= i; ++j ) 
          matrix[ j ][ i ] = 0;

      QuadratureType quad( entity, 2 * (polynomialOrder - 1) );
      const int numQuadraturePoints = quad.nop();
      for( int pt = 0; pt < numQuadraturePoints; ++pt ) {
        // calc Jacobian inverse before volume is evaluated 
        const FieldMatrix< double, dimension, dimension > &inv
                    = geometry.jacobianInverseTransposed( quad.point( pt ) );
        const double volume = geometry.integrationElement( quad.point( pt ) );

        for( int i = 0; i < matrixSize; ++i ) {
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
        for( int i = 0; i < matrixSize; ++i ) 
          for ( int j = 0; j <= i; ++j ) 
            matrix[ j ][ i ] += (mygrad[ i ][ 0 ] * mygrad[ j ][ 0 ]) * weight;
      }
      
      // symmetrize matrix
      for( int i = 0; i < matrixSize; ++i ) 
        for( int j = matrixSize; j > i; --j ) 
          matrix[ j ][ i ] = matrix[ i ][ j ];
    }
  }; // end class



  // Assebmler for right rand side
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
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      typedef typename EntityType :: Geometry GeometryType;
      
      const DiscreteFunctionSpaceType &discreteFunctionSpace
        = discreteFunction.space();
  
      discreteFunction.clear(); //discreteFunction auf Null setzen

      const IteratorType endit = discreteFunctionSpace.end();
      for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
      {
        //it* Pointer auf ein Element der Entity
        const GeometryType &geometry = (*it).geometry(); //Referenz auf Geometrie
      
        LocalFunctionType localFunction = discreteFunction.localFunction( *it ); 
        const BaseFunctionSetType &baseFunctionSet //BaseFunctions leben immer auf Refernzelement!!!
          = discreteFunctionSpace.baseFunctionSet( *it ); 

        CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrd ); //0 --> codim 0
        const int numDofs = localFunction.numDofs(); //Dofs = Freiheitsgrade (also die Unbekannten)
        for( int i = 0; i < numDofs; ++i )
        {
          RangeType phi, psi; //R"uckgabe-Funktionswerte
        
          const int numQuadraturePoints = quadrature.nop();
          for( int qP = 0; qP < numQuadraturePoints; ++qP )
          {
            const double det
              = geometry.integrationElement( quadrature.point( qP ) );
            function.evaluate( geometry.global( quadrature.point( qP ) ), phi );
            baseFunctionSet.evaluate( i, quadrature, qP, psi ); //i = i'te Basisfunktion; qP Quadraturpunkt
            localFunction[ i ] += det * quadrature.weight( qP ) * (phi * psi);
          }
        }
      }
    }
  };

} // end namespace 

#endif


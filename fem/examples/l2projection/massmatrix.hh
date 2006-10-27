#ifndef DUNE_MASSMATRIX_HH
#define DUNE_MASSMATRIX_HH

#include <dune/fem/examples/poisson/feop.hh>
#include <dune/fem/operator/feop/spmatrix.hh>
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune {
  
  /** \brief The mass matrix
   *
   * \tparam polOrd The quadrature order
   *
   * \todo Giving the quadrature order as a template parameter is
   * a hack.  It would be better to determine the optimal order automatically.
   */
  template <class DiscFunctionType, class TensorImp, int polOrd>
  class MassMatrixFEOp :
    public FEOp<DiscFunctionType,
                        SparseRowMatrix<double>,
                        MassMatrixFEOp<DiscFunctionType,TensorImp, polOrd> >
    {
      
    typedef TensorImp TensorType;
    typedef typename DiscFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FEOp<DiscFunctionType,
                      SparseRowMatrix<double>,
     MassMatrixFEOp<DiscFunctionType, TensorType, polOrd> >::OpMode OpMode;

    
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    
    //! Quadrature
    Quadrature < typename FunctionSpaceType::RangeFieldType,
               GridType::dimension> quad; 
    
  public:

    //! Returns the actual matrix if it is assembled
    /** \todo Should this be in a base class? */
    const SparseRowMatrix<double>* getMatrix() const {
      assert(this->matrix_);
      return this->matrix_;
    }

    //! Constructor
    MassMatrixFEOp( const typename DiscFunctionType::FunctionSpaceType &f, OpMode opMode ):
      FEOp<DiscFunctionType,SparseRowMatrix<double>,
             MassMatrixFEOp<DiscFunctionType,TensorImp, polOrd> >( f, opMode ) ,
          quad( f.begin()->geometry().type() , polOrd )
    {
    }
        
    //! ???
    SparseRowMatrix<double>* newEmptyMatrix( ) const {
      return new SparseRowMatrix<double>( this->functionSpace_.size () , 
                                          this->functionSpace_.size () , 
                                          30);
    }
        
    //! ???
    template <class EntityType>
    double getLocalMatrixEntry( EntityType &entity, const int i, const int j ) const 
    {
      typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      const int dim = GridType::dimension;
      double val = 0;
      
      FieldVector<double, dim> tmp(1.0);
      
      const BaseFunctionSetType & baseSet 
        = this->functionSpace_.getBaseFunctionSet( entity );
      
      RangeType v1 (0.0);
      RangeType v2 (0.0);

      for ( int pt=0; pt < quad.nop(); pt++ ) 
      {      
          const double det = entity.geometry().integrationElement( quad.point(pt) ); 
          baseSet.eval(i,quad,pt,v1); 
          baseSet.eval(j,quad,pt,v2); 
          val += ( v1 * v2 ) * det * quad.weight(pt);
        }
      
      return val;
    }
    
    //! ???
    template < class  EntityType, class MatrixType>
    void getLocalMatrix( EntityType &entity, const int matSize, MatrixType& mat) const 
    {
      int i,j;
      typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      const int dim = GridType::dimension;
      const BaseFunctionSetType & baseSet 
        = this->functionSpace_.getBaseFunctionSet( entity );
      
      /** \todo What's the correct type here? */
      static FieldVector<double, dim> tmp(1.0);
      
      static RangeType v[30];
      // Check magic constant. Otherwise program will fail in loop below
      assert(matSize <= 30); 

      for(i=0; i<matSize; i++) 
        for (j=0; j<=i; j++ ) 
          mat[j][i]=0.0;

      for ( int pt=0; pt < quad.nop(); pt++ )
        {
          const double det = entity.geometry().integrationElement( quad.point(pt) ); 
          
          for(i=0; i<matSize; i++) 
            baseSet.eval(i,quad,pt,v[i]); 
          
          for(i=0; i<matSize; i++) 
            for (j=0; j<=i; j++ ) 
              mat[j][i] += ( v[i] * v[j] ) * det * quad.weight(pt);
        }

      for(i=0; i<matSize; i++) 
        for (j=matSize; j>i; j--) 
          mat[j][i] = mat[i][j];
    }

  protected:
    DiscFunctionType *_tmp;
  };

} // namespace Dune

#endif

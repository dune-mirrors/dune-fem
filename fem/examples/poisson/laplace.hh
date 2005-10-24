#ifndef DUNE_LAPLACE_HH
#define DUNE_LAPLACE_HH

#include <dune/fem/feop/spmatrix.hh>
#include <dune/fem/feop.hh>
#include <dune/quadrature/fixedorder.hh>

namespace Dune 
{

    /** \brief The Laplace operator
     */
    template <class DiscFunctionType, class TensorType, int polOrd>
    class LaplaceFEOp : 
        public FEOp<DiscFunctionType,
                    SparseRowMatrix<double>,
                    LaplaceFEOp<DiscFunctionType,TensorType, polOrd> > {
        
        //! The corresponding function space type
        typedef typename DiscFunctionType::FunctionSpaceType FunctionSpaceType;

        //! The grid
        typedef typename FunctionSpaceType::GridType GridType;

        //! The grid's dimension
        enum { dim = GridType::dimension };

        //! ???
        typedef typename FunctionSpaceType::JacobianRangeType JacobianRange;

        //! ???
        typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
        
        //! ???
        typedef typename FEOp<DiscFunctionType,
                              SparseRowMatrix<double>,
                              LaplaceFEOp<DiscFunctionType, TensorType, polOrd> >::OpMode OpMode;
        
        
        mutable JacobianRange grad;
        mutable JacobianRange othGrad;


        enum { maxnumOfBaseFct = 100 };
        mutable JacobianRange mygrad[maxnumOfBaseFct];
        
    public:

        //! ???
        FixedOrderQuad < typename FunctionSpaceType::RangeFieldType, typename
                   FunctionSpaceType::DomainType , polOrd > quad;
        
        //! ???
        DiscFunctionType *stiffFunktion_;

        //! ???
        TensorType *stiffTensor_;
        
        //! ???
        LaplaceFEOp( const typename DiscFunctionType::FunctionSpaceType &f, OpMode opMode ): 
            FEOp<DiscFunctionType,
                 SparseRowMatrix<double>,
                 LaplaceFEOp<DiscFunctionType,TensorType, polOrd> >( f, opMode ) , 
            quad ( *(f.begin() )), stiffFunktion_(NULL), stiffTensor_(NULL)
        {
        }
        
        //! ???
        LaplaceFEOp( const DiscFunctionType &stiff, const typename DiscFunctionType::FunctionSpaceType &f, OpMode opMode ): 
            FEOp<DiscFunctionType,
            SparseRowMatrix<double>,
            LaplaceFEOp<DiscFunctionType,TensorType, polOrd> >( f, opMode ) ,
            quad ( *(f.begin())), stiffFunktion_(&stiff), stiffTensor_(NULL)
        { 
        }
        
        //! ???
        LaplaceFEOp( TensorType &stiff, const typename DiscFunctionType::FunctionSpaceType &f, OpMode opMode ): 
            FEOp<DiscFunctionType,SparseRowMatrix<double>,LaplaceFEOp<DiscFunctionType,TensorType, polOrd> >( f, opMode ) ,
            quad ( *(f.begin())), stiffFunktion_(NULL), stiffTensor_(&stiff)
        { 
        }
        
        //! Returns the actual matrix if it is assembled
        const SparseRowMatrix<double>* getMatrix() const {
            assert(this->matrix_);
            return this->matrix_;
        }
        
        //! Creates a new empty matrix
        SparseRowMatrix<double>* newEmptyMatrix( ) const 
        {
            return new SparseRowMatrix<double>( this->functionSpace_.size () , 
                                                this->functionSpace_.size () , 
                                                15 * dim);
        }
        
        //! Prepares the local operator before calling apply()
        void prepareGlobal ( const DiscFunctionType &arg, DiscFunctionType &dest ) 
        {
            this->arg_  = &arg;
            this->dest_ = &dest;
            assert(this->arg_ != 0); assert(this->dest_ != 0);
            dest.clear();
        }
        
        //! return the matrix entr that belong to basis function i and j 
        template <class EntityType>
        double getLocalMatrixEntry( EntityType &entity, const int i, const int j ) const 
        {
            enum { dim = GridType::dimension };
            typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
            
            const BaseFunctionSetType & baseSet = this->functionSpace_.getBaseFunctionSet( entity );
            
            double val = 0.;
            for ( int pt=0; pt < quad.nop(); pt++ ) 
                {      
                    baseSet.jacobian(i,quad,pt,grad); 
                    
                    // calc Jacobian inverse before volume is evaluated 
                    const FieldMatrix<double,dim,dim>& inv = entity.geometry().jacobianInverseTransposed(quad.point(pt));
                    const double vol = entity.geometry().integrationElement(quad.point(pt)); 
            
                    // multiply with transpose of jacobian inverse 
                    grad[0] = FMatrixHelp :: mult ( inv, grad[0] );
                    
                    if( i != j ) 
                        {
                            baseSet.jacobian(j,quad,pt,othGrad); 
                            
                            // multiply with transpose of jacobian inverse 
                            othGrad[0] = FMatrixHelp :: multTransposed ( inv, othGrad[0] );
                            val += ( grad[0] * othGrad[0] ) * quad.weight( pt ) * vol;
                        }
                    else 
                        {
                            val += ( grad[0] * grad[0] ) * quad.weight( pt ) * vol;
                        }
                }

            return val;
        }
        
        //! ???
        template < class  EntityType, class MatrixType>
        void getLocalMatrix( EntityType &entity, const int matSize, MatrixType& mat) const {
            enum { dim = GridType::dimension };

            typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
            
            const BaseFunctionSetType & baseSet = this->functionSpace_.getBaseFunctionSet( entity );
            
            assert( matSize <= maxnumOfBaseFct );
            
            for(int i=0; i<matSize; i++) 
                for (int j=0; j<=i; j++ ) 
                    mat[j][i]=0.0;

            for ( int pt=0; pt < quad.nop(); pt++ ) {

                // calc Jacobian inverse before volume is evaluated 
                const FieldMatrix<double,dim,dim>& inv = entity.geometry().jacobianInverseTransposed(quad.point(pt));
                const double vol = entity.geometry().integrationElement(quad.point(pt));

                for(int i=0; i<matSize; i++) {
      
                    baseSet.jacobian(i,quad,pt,mygrad[i]); 
      
                    // multiply with transpose of jacobian inverse 
                    mygrad[i][0] = FMatrixHelp :: mult ( inv,mygrad[i][0] );
                }
                    
                typename FunctionSpaceType::RangeType ret;
                    
                if(stiffTensor_){
                    stiffTensor_->evaluate(entity.geometry().global(quad.point(pt)),ret);
                    ret[0] *= quad.weight( pt );
                    for(int i=0; i<matSize; i++) 
                        for (int j=0; j<=i; j++ ) 
                            mat[j][i] += ( mygrad[i][0] * mygrad[j][0] ) * ret[0] * vol;
                }
                else{
                    for(int i=0; i<matSize; i++) 
                        for (int j=0; j<=i; j++ ) 
                            mat[j][i] += ( mygrad[i][0] * mygrad[j][0] ) * quad.weight( pt ) * vol;
                }
            }
            
            for(int i=0; i<matSize; i++) 
                for (int j=matSize; j>i; j--) 
                    mat[j][i] = mat[i][j];
            
        }
        
    }; // end class
    
} // end namespace 

#endif


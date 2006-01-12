#ifndef __DUNE_FEM_CC__
#define __DUNE_FEM_CC__

#include "gausspoints.hh"

#include <dune/fem/feoperator.hh>
#include <dune/fem/feop/spmatrix.hh>
#include <dune/fem/inverseoperators.hh>

#include <dune/quadrature/fixedorder.hh>

using namespace Adi;

namespace Dune {

// projection of the rhs 
template <class DiscreteFunctionType> 
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  typedef typename FunctionSpaceType::GridType GridType;
  typedef typename FunctionSpaceType::RangeType Range;
  typedef typename GridType::template Codim<0>::LevelIterator LevelIterator;

public:
  //! a small helper for project lagrange
template <class FunctionType>
  void assignLagrange(FunctionType& f, 
                      LocalFunctionType& lf, 
                      LevelIterator& it,
                      int polOrd, Int2Type<2>) {
    GaussPoints gausspts;
    std::vector<double> gp;

    for (int i = 0; i < polOrd; ++i) {
      gp.push_back(gausspts.point(polOrd, i));
    }

    for (int l = 0; l < polOrd*polOrd; ++l) {  
      const int i = l%polOrd;
      const int j = l/polOrd;
      
      FieldVector<double, 2> local;
      local[0] = gp[i];
      local[1] = gp[j];
      
      lf[l] = f(it->geometry().global(local));
    }
  }
  
  template <class FunctionType>
  void assignLagrange(FunctionType& f,
                      LocalFunctionType& lf, 
                      LevelIterator& it,
                      int polOrd, Int2Type<3>) {
    GaussPoints gausspts;
    std::vector<double> gp;
    
    for (int i = 0; i < polOrd; ++i) {
      gp.push_back(gausspts.point(polOrd, i));
    }

    for (int l = 0; l < polOrd*polOrd*polOrd; ++l) {
      const int i = l%polOrd;
      const int j = (l/polOrd)%polOrd;
      const int k = l/(polOrd*polOrd);
      
      FieldVector<double, 3> local;
      local[0] = gp[i];
      local[1] = gp[j];
      local[2] = gp[k];
      
      lf[l] = f(it->geometry().global(local));
    }
  }
  

  //! A simple projection for the lagrange test case
  template <class FunctionType>
  void projectLagrange(FunctionType& f,
                       DiscreteFunctionType& df,
                       int level) {
    enum { dimRange = FunctionSpaceType::DimRange };
    enum { dimDomain = FunctionSpaceType::DimDomain };
    assert(dimRange == 1);
  
    const int polOrd = FunctionSpaceType::PolOrd;

    df.clear();
    const GridType& grid = df.getFunctionSpace().getGrid();

    // Temporaries
    LocalFunctionType lf = df.newLocalFunction();
    Range res;
    Range phi;

    LevelIterator endit = grid.template lend<0>(grid.maxlevel());
    for (LevelIterator it = grid.template lbegin<0>(grid.maxlevel()); it != endit; ++it){
      df.localFunction(*it, lf);

      assignLagrange(f, lf, it, polOrd, Int2Type<dimDomain>() );
    }
  }

  //! An optimized L2 projection for orthonormal DG base functions
  template <int polOrd, class FunctionType>
  void projectDG(const FunctionType& f,
                 DiscreteFunctionType& df,
                 int level) {
    //- Local typedefs
    typedef FixedOrderQuad<
      typename FunctionSpaceType::RangeFieldType,
      typename FunctionSpaceType::DomainType, 
      polOrd> QuadType;

    enum { dimRange = FunctionSpaceType::DimRange };

    //- Actual code
    df.clear();
    const GridType& grid = df.getFunctionSpace().getGrid();

    // Temporaries
    LocalFunctionType lf = df.newLocalFunction();
    Range res;
    Range phi;

    LevelIterator endit = grid.template lend<0>(level);
    for (LevelIterator it = grid.template lbegin<0>(level); it != endit; ++it){
      QuadType quad(*it);
      df.localFunction(*it, lf);
      BaseFunctionSetType& baseSet =
        df.getFunctionSpace().getBaseFunctionSet(*it);

      for (int i = 0; i < baseSet.numBaseFunctions(); ++i) {
        double sum(0.0);
        double norm = 0.0;
        for (int l = 0; l < quad.nop(); ++l) {
          double res = f(it->geometry().global(quad.point(l)));      
          baseSet.eval(i, quad, l, phi);
          double ds = it->geometry().integrationElement(quad.point(l));
          
          res *= ds*quad.weight(l)*phi[0];
          sum += res;

          norm += ds*quad.weight(l)*phi[0]*phi[0];
        } // end quadrature loop

        lf[i] = sum/norm;
      /*
      for (int i = 0; i < baseSet.getNumberOfBaseFunctions(); ++i) {
        Range sum(0.0);
        double norm = 0.0;
        for (int l = 0; l < quad.nop(); ++l) {
          f.evaluate(it->geometry().global(quad.point(l)), res);
          baseSet.eval(i, quad, l, phi);
          double ds = it->geometry().integrationElement(quad.point(l));
          
          res *= ds*quad.weight(l)*phi[0];
          sum += res;

          norm += ds*quad.weight(l)*phi[0]*phi[0];
        } // end quadrature loop

        for (int j = 0; j < dimRange; ++j) {
          lf[i*dimRange + j] = sum[j]/norm;
        }
      */
      } // end base function loop
    }
  }

  //! A local L2 projection
  //! \warning: Only works for dimRange == 1
  template <int polOrd, class FunctionType>
  void projectLocal(FunctionType& f,
                    DiscreteFunctionType& df, 
                    int level) {
    //- Local typedefs
    typedef FixedOrderQuad<
      typename FunctionSpaceType::RangeFieldType,
      typename FunctionSpaceType::DomainType, 
      polOrd> QuadType;
    // * Hack!
    const int nBaseFct = 27;

    typedef FieldVector<double, nBaseFct> VectorType;
    typedef FieldMatrix<double, nBaseFct, nBaseFct> MatrixType;

    // * do we need that?
    df.clear();
    const GridType& grid = df.getFunctionSpace().getGrid();

    // Temporaries
    Range res;
    Range phi_i;
    Range phi_j;
    VectorType tmp;

    LevelIterator endit = grid.template lend<0>(level);
    for (LevelIterator it = grid.template lbegin<0>(level); it != endit; ++it){
    assert(df.getFunctionSpace().getBaseFunctionSet(*it).numBaseFunctions() 
           == nBaseFct);

      QuadType quad(*it);
      LocalFunctionType lf = df.localFunction(*it);
      const BaseFunctionSetType& baseSet =
        df.getFunctionSpace().getBaseFunctionSet(*it);

      MatrixType mat(0.0);
      VectorType rhs(0.0);
      for (int l = 0; l < quad.nop(); ++l) {
          // Get values at integration points
        f.evaluate(it->geometry().global(quad.point(l)), res);
        double ds = it->geometry().integrationElement(quad.point(l));
 
        // Distribute the stuff on the different degrees of freedom
        for (int i = 0; i < nBaseFct; ++i) {
          baseSet.eval(i, quad, l, phi_i);
          // rhs
          rhs[i] += ds*quad.weight(l)*res[0]*phi_i[0];

          // mat
          for (int j = 0; j < nBaseFct; ++j) {
            baseSet.eval(j, quad, l, phi_j);
            mat[i][j] += ds*quad.weight(l)*phi_i[0]*phi_j[0];
          }
        }
      } // end quadrature loop

      //mat.invert();
      //mat.umv(rhs, tmp);
      mat.solve(tmp, rhs);
      for (int i = 0; i < nBaseFct; ++i) {
        lf[i] = tmp[i];
      }
    } // end grid loop 
  }

  template <int polOrd, class FunctionType> 
  void project (FunctionType &f, DiscreteFunctionType &discFunc, int level)
  {
    const typename DiscreteFunctionType::FunctionSpace 
        & functionSpace_= discFunc.getFunctionSpace();  
  
    discFunc.clear();
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename GridType::template Codim<0>::LevelIterator LevelIterator;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      

    GridType & grid = functionSpace_.getGrid();

    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    LevelIterator it = grid.template lbegin<0> ( level );
    LevelIterator endit = grid.template lend<0> ( level );
    FixedOrderQuad <typename FunctionSpaceType::RangeFieldType,
              typename FunctionSpaceType::DomainType , polOrd > quad ( *it );
              
    for( ; it != endit ; ++it)
    {
      LocalFuncType lf = discFunc.localFunction( *it ); 
      double det = (*it).geometry().integrationElement(quad.point(0));
      
      const typename FunctionSpaceType::BaseFunctionSetType & set = 
            functionSpace_.getBaseFunctionSet(*it);

      for(int i=0; i<lf.numDofs(); i++)
      {
        for(int qP = 0; qP < quad.nop(); qP++)
        {
          f.evaluate((*it).geometry().global( quad.point(qP) ), ret);
          set.eval(i,quad,qP,phi);
          lf[i] += det * quad.weight(qP) * (ret * phi);
        }
      }
    }
  }

  template <int polOrd, class FunctionType> 
  void lumpi (int level, FunctionType &f, DiscreteFunctionType &discFunc, double time = 0.0)
  {
    const typename DiscreteFunctionType::DiscreteFunctionSpaceType
        & functionSpace_= discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename GridType::template Codim<0>::LeafIterator LeafIterator;
    typedef typename GridType::template Codim<0>::LevelIterator LevelIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename DiscreteFunctionType::
        LocalFunctionType LocalFuncType;
      
    const GridType & grid = functionSpace_.grid();

    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    discFunc.clear();

    LevelIterator endit = grid.template lend<0> ( level );
    LevelIterator it = grid.template lbegin<0> ( level ); 
    
    FixedOrderQuad < typename FunctionSpaceType::RangeFieldType,
    typename FunctionSpaceType::DomainType , polOrd > quad ( *it );
    
    double sum;
    double intWeight;
    
    for( ; it != endit ; ++it)
    {
      EntityType & en = (*it);
      LocalFuncType lf = discFunc.localFunction ( en );

      int numDof = lf.numDofs ();  
      for(int i=0; i<numDof; i++)
      {
        sum = 0.0;
        lf[i] = 0.0;
        for(int qP = 0; qP < quad.nop(); qP++)
        {
          intWeight = quad.weight(qP) * en.geometry().integrationElement(quad.point(qP));
          f.evaluate(en.geometry().global(quad.point(qP)),time,ret);
          lf[i] += intWeight * ret[i];
          sum += intWeight;
        }
        lf[i] /= sum;
      }
    }
  }
};

// used for calculation of the initial values 
template <class DiscreteFunctionType> 
class LagrangeInterpolation
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <int polOrd, class FunctionType> 
  void interpol (int level, FunctionType &f, DiscreteFunctionType &discFunc)
  {
    const typename DiscreteFunctionType::FunctionSpace 
        & functionSpace_= discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      

    typename FunctionSpaceType::RangeType ret (0.0);

    typedef typename FunctionSpaceType::IteratorType IteratorType;
    IteratorType endit = functionSpace_.end(); 
    for(IteratorType it = functionSpace_.begin(); it != endit ; ++it)
    {
      LocalFuncType lf = discFunc.localFunction(*it); 

      int numDof = lf.numDofs ();  
      for(int i=0; i<numDof; i++)
      {
        f.evaluate(it->geometry()[i],ret);
        lf[i] = ret[0];
      }
    }
  }
};


//************************************************************************
//
// Local Operators 
//
//************************************************************************
template <class DiscFunctionType>
class ConstLocalFEOp :
  public FiniteElementOperator<DiscFunctionType,
                               SparseRowMatrix<double>,
                               ConstLocalFEOp<DiscFunctionType> > {

  typedef typename FiniteElementOperator<DiscFunctionType,
               SparseRowMatrix<double>,
               ConstLocalFEOp<DiscFunctionType> >::OpMode OpMode;
               
public:
  ConstLocalFEOp( const typename DiscFunctionType::FunctionSpace &f, double localmatrix[4][4], OpMode opMode )://= ON_THE_FLY ) :
    FiniteElementOperator<DiscFunctionType,SparseRowMatrix<double>,ConstLocalFEOp<DiscFunctionType> >( f, opMode ) {
    
    for ( int i=0; i<4; i++ ) {
      for ( int j=0; j<4; j++ ) {
        localmatrix_[i][j] = localmatrix[i][j];
      }
    }
  }

  SparseRowMatrix<double>* newEmptyMatrix( ) const {
    return new SparseRowMatrix<double>( this->functionSpace_.size ( this->functionSpace_.getGrid().maxlevel() ) , 
          this->functionSpace_.size ( this->functionSpace_.getGrid().maxlevel() ) , 
          10, 0.0 );
  }

  template <class EntityType>
  double getLocalMatrixEntry( EntityType &entity, const int i, const int j ) const {
    return localmatrix_[i][j];
  }

  template < class  EntityType, class MatrixType>
  void getLocalMatrix( EntityType &entity, const int matSize, MatrixType& mat) const {
    for(int i=0; i<matSize; i++) 
      for (int j=0; j<matSize; j++ ) 
        mat(i,j) =  localmatrix_[i][j]; 
  
    return;
  }

protected:
  double localmatrix_[4][4];
};




// calculates || u-u_h ||_L2
template <class DiscreteFunctionType> 
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <int polOrd, class FunctionType> 
  double norm (int level, FunctionType &f, DiscreteFunctionType &discFunc,
      double time)
  {
    const typename DiscreteFunctionType::FunctionSpace 
        & functionSpace_= discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename GridType::template Traits<0>::LevelIterator LevelIterator;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    
    typename FunctionSpaceType::RangeType ret (0.0);
    typename FunctionSpaceType::RangeType phi (0.0);

    double sum = 0.0;
    LevelIterator endit = grid.template lend<0> ( level );
    LevelIterator it = grid.template lbegin<0> ( level );
    //FaceCenterQuad < typename FunctionSpaceType::RangeField,
    //typename FunctionSpaceType::Domain > // , polOrd > 
    //  quad ( *it );
    
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    IteratorType endit = functionSpace_.end(); 
    IteratorType it = functionSpace_.begin();

    assert(it != endit );
    FixedOrderQuad < typename FunctionSpaceType::RangeFieldType,
               typename FunctionSpaceType::DomainType , polOrd > quad ( *it );
    
    for(; it != endit ; ++it)
    {
      double det = (*it).geometry().integrationElement(quad.point(0));
      LocalFuncType lf = discFunc.localFunction(*it); 
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate((*it),quad,qP,phi);
        sum += det * quad.weight(qP) * SQR(ret[0] - phi[0]);
      }
    }
    return sqrt(sum);
  }
};

} // end namespace 

#endif


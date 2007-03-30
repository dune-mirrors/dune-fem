#ifndef HDIV_PROJECTION_CC 
#define HDIV_PROJECTION_CC 

//- Dune includes
#include <dune/grid/common/referenceelements.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

//- Dune-fem includes 
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>

namespace Dune {

template <class DiscreteFunctionType>
class HdivTest
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeType RangeType; 
  typedef typename FunctionSpaceType::DomainType DomainType; 
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType; 
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType; 
  typedef typename FunctionSpaceType::GridPartType GridPartType; 
  
    enum { dimR = 1 };
    enum { dimD = FunctionSpaceType::DimDomain - 1 };
    enum { polOrdN = FunctionSpaceType::polOrd };

  typedef FunctionSpace<DomainFieldType,RangeFieldType,dimD,dimR> FaceSpaceType;
  //typedef FunctionSpace<DomainFieldType,RangeFieldType,FunctionSpaceType::DimDomain,1> FaceSpaceType;
  typedef FunctionSpace<DomainFieldType,RangeFieldType,FunctionSpaceType::DimDomain,dimR> ElSpaceType;

  typedef DiscontinuousGalerkinSpace<FaceSpaceType,GridPartType, polOrdN > FaceDiscreteSpaceType; 
  //typedef LagrangeDiscreteFunctionSpace<FaceSpaceType,GridPartType, polOrdN > FaceDiscreteSpaceType; 
  enum { gradPolOrd = ((polOrdN - 1) < 0) ? 0 : (polOrdN - 1) };
  typedef DiscontinuousGalerkinSpace<ElSpaceType,GridPartType, gradPolOrd, CachingStorage > ElementGradientSpaceType; 
  enum { bubblePolOrd = ((polOrdN - 2) < 0) ? 0 : (polOrdN - 2) };
  typedef DiscontinuousGalerkinSpace<ElSpaceType,GridPartType, bubblePolOrd , CachingStorage > ElementDiscreteSpaceType; 

 public:
  static double localMassConserve(const DiscreteFunctionType &discFunc, int polyOrder = -1 ) 
  {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 
   
    enum { dim = GridType::dimension };
    
    const FunctionSpaceType& space =  discFunc.space();
    const GridPartType & gridPart = space.gridPart();
    int polOrd = (polyOrder <0) ? (2 * space.order() + 2) : polyOrder;
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    RangeType ret (0.0);
    RangeType neighRet (0.0);
    
    // Get quadraturae rule
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
     
    typedef CachingQuadrature <GridPartType , 1> FaceQuadratureType; 

    double sum = 0.0;

    //const LocalIdSetType & idSet = space.grid().localIdSet();
    
    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      EntityType & en = *it;
      const LocalFuncType lf = discFunc.localFunction(en);
     
      double localValue = 0.0;
      
      IntersectionIteratorType endnit = gridPart.iend(en);
      for(IntersectionIteratorType nit = gridPart.ibegin(en);
          nit != endnit; ++nit)
      {
        // only interior faces are considered 
        if(nit.neighbor())
        {
          EntityPointerType neighEp = nit.outside();
          EntityType&            nb = *neighEp;
          
          {
            FaceQuadratureType faceQuadInner(gridPart, nit, polOrd, FaceQuadratureType::INSIDE);
            FaceQuadratureType faceQuadOuter(gridPart, nit, polOrd, FaceQuadratureType::OUTSIDE);

            const LocalFuncType neighLf = discFunc.localFunction(nb);
            
            const int quadNop = faceQuadInner.nop();
            for (int l = 0; l < quadNop ; ++l)
            {
              DomainType normal = 
                nit.integrationOuterNormal(faceQuadInner.localPoint(l));

              lf.evaluate(faceQuadInner,l,ret);
              neighLf.evaluate(faceQuadOuter,l,neighRet);
              
              ret += neighRet;
              ret *= 0.5;
              
              double val = ret * normal; 
              val *= faceQuadInner.weight(l);

              localValue += val;
            }
          }
        }

        // only interior faces are considered 
        if(nit.boundary())
        {
          FaceQuadratureType faceQuadInner(gridPart, nit, polOrd, FaceQuadratureType::INSIDE);

          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l)
          {
            DomainType normal = 
              nit.integrationOuterNormal(faceQuadInner.localPoint(l));

            lf.evaluate(faceQuadInner,l,ret);
            
            double val = ret * normal; 
            val *= faceQuadInner.weight(l);

            localValue += val;
          }
        }
      } // end of intersection iterator 
      sum += std::abs(localValue);
    }

    return sum;
  }
  
  void bubbleQuad(const DomainType & point, RangeType & result) const 
  {
    double ret = 1.0;
    for(int i=0; i<DomainType::dimension; ++i)
    {
      ret *= 4.0 * point[i] * (1.0-point[i]);
    }
    result = ret;
  }
  

  // bubble function for simplex 
  double bubbleSimplex(const DomainType & point) const 
  {
    double ret = 1.0;
    double lambdaNull = 1.0;
    for(int i=0; i<DomainType::dimension; ++i) 
    {
      ret *= point[i];
      lambdaNull -= point[i];
    }
    ret *= lambdaNull;
    return ret; 
  }
  
  // bubble function for simplex 
  double gradBubbleSimplex(int deri, const DomainType & point) const 
  {
    assert( deri >= 0 );
    assert( deri < DomainType::dimension );
      
    double ret = 1.0;
    double lambdaNull = 1.0;

    for(int i=0; i<DomainType::dimension; ++i) 
    {
      if( i != deri )  ret *= point[i];
      lambdaNull -= point[i];
    }
    ret *= lambdaNull;
    
    double first = -1.0;
    for(int i=0; i<DomainType::dimension; ++i)
    {
      first *= point[i]; 
    }
    return (ret+first); 
  }
  
  void gradientBubbleSimplex(const DomainType & point, DomainType & result) const 
  {
    for(int i=0; i<DomainType::dimension; ++i) 
    {
      result[i] = gradBubbleSimplex(i,point); 
    }
  }

  template <class BaseFunctionSetType, class QuadratureType, class JacoRangeType>
  void evalBubbleBase(const BaseFunctionSetType & bSet, int func, QuadratureType & quad,
      int quadPoint, JacoRangeType & result) const 
  {
    const DomainType & point = quad.point(quadPoint);
    
    bSet.jacobian(func,quad,quadPoint,result);
    double bub = bubbleSimplex(point);
    result *= bub; 
    
    typedef typename ElementDiscreteSpaceType :: RangeType ElRangeType; 
    ElRangeType bTmp; 
    bSet.eval(func,quad,quadPoint,bTmp);
    RangeType tmp; 
    gradientBubbleSimplex(point,tmp);

    tmp *= bTmp[0];
    result[0] += tmp;
  }

  // only works for 2d right now 
  void curl(const DomainType & arg, DomainType & dest) const 
  {
    dest[0] =  arg[1];
    dest[1] = -arg[0];
  }

  void evalBubbleFunc(int comp, const DomainType & point, 
                      RangeType & result) const 
  {
    assert( comp == 0 );
    gradientBubbleSimplex(point,result);
  }

  template <class DiscreteFunctionSpaceImp, class EntityType, 
            class LocalFunctionType, class MatrixType, class VectorType> 
  void bubblePart(const DiscreteFunctionSpaceImp & space, 
                  EntityType & en, 
                  const LocalFunctionType & uLF, const int startRow , 
                  MatrixType & matrix, VectorType& rhs ) const 
  {
    typedef DiscreteFunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType; 

    typedef typename ElementDiscreteSpaceType::JacobianRangeType JacobianRangeType; 
    typedef typename FunctionSpaceType :: DomainType DomainType;
    
    enum { dim = GridType::dimension };
    
    const BaseFunctionSetType & bSet = space.baseFunctionSet(en); 
    int polOrd = 2 * space.order() + 1;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    const int localRows = bSet.numBaseFunctions(); 
    const int cols = uLF.numDofs();

    QuadratureType quad (en,polOrd);
    DomainType result;
    JacobianRangeType val;
    DomainType bVal; 
    DomainType aVal; 

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      const double intel = quad.weight(l) * 
        en.geometry().integrationElement(quad.point(l));
      
      uLF.evaluate(quad,l,result);
      for(int i=0; i<localRows; ++i)     
      {
        // we might have other row 
        int row = startRow + i;
        
        evalBubbleBase(bSet,i,quad,l,val);
        curl(val[0],aVal);

        double r = aVal * result; 
        r *= intel; 
        rhs[row] += r; 

        // for cols make matrix 
        for(int j=0; j<cols; ++j)     
        {
          uLF.baseFunctionSet().eval(j,quad,l,bVal);
          double t = aVal * bVal;
          t *= intel; 
          matrix[row][j] += t;
        }
      }
    }
  }

  template <class EntityType, class LocalFunctionType,
            class MatrixType, class VectorType> 
  void gradientPart(const ElementGradientSpaceType & space, 
                    EntityType & en, 
                    const LocalFunctionType & uLF, const int startRow , 
                    MatrixType& matrix, VectorType& rhs ) const 
  {
    typedef typename ElementGradientSpaceType::BaseFunctionSetType BaseFunctionSetType;
    const BaseFunctionSetType & bSet = space.baseFunctionSet(en); 
    int polOrd = 2 * space.order() + 1;
    
    typedef CachingQuadrature <GridPartType , 0> QuadratureType; 

    const int localRows = bSet.numBaseFunctions(); 
    const int cols = uLF.numDofs();

    QuadratureType quad (en,polOrd);

    RangeType result;
    RangeType uPhi;

    typedef typename ElementGradientSpaceType::JacobianRangeType GradJacobianRangeType; 
    GradJacobianRangeType gradPhi;

    const int quadNop = quad.nop();
    for (int l = 0; l < quadNop ; ++l)
    {
      const double intel = quad.weight(l) * 
        en.geometry().integrationElement(quad.point(l));
      
      uLF.evaluate(quad,l,result);
      for(int i=0; i<localRows; ++i)     
      {
        // we might have other row 
        int row = startRow + i;
        
        bSet.jacobian(i,quad,l,gradPhi);
        double uDGVal = result * gradPhi[0];
        rhs[row] += uDGVal * intel;  

        // for cols make matrix 
        for(int j=0; j<cols; ++j)     
        {
          uLF.baseFunctionSet().eval(j,quad,l,uPhi);
          double val = uPhi * gradPhi[0]; 
          matrix[row][j] += val * intel;
        }
      }
    }
  }

  //template <class AdaptationType>
  void project(const DiscreteFunctionType &uDG,
               DiscreteFunctionType & velo) const 
  {
    typedef typename DiscreteFunctionType::Traits::DiscreteFunctionSpaceType FunctionSpaceType;
    enum { localBlockSize = FunctionSpaceType::localBlockSize };

    typedef typename FunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
    typedef typename FunctionSpaceType::Traits::GridType GridType;
    typedef typename FunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename FunctionSpaceType::Traits::IteratorType Iterator;
    typedef typename GridType :: template Codim<0> :: Entity EntityType;
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType; 
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType; 
    typedef typename FunctionSpaceType :: DomainType DomainType;
    
    enum { dim = GridType::dimension };
    typedef typename GridType :: ctype coordType;
    
    const FunctionSpaceType& space = uDG.space();
    GridPartType & gridPart = const_cast<GridPartType &> (space.gridPart());
    int polOrd = 2 * space.order() + 2;

    // only working for polOrd = 1 at the moment 
    //assert( space.order() == 1 );
    
    FaceDiscreteSpaceType faceSpace(gridPart);
    ElementDiscreteSpaceType elSpace(gridPart);
    ElementGradientSpaceType gradSpace(gridPart);

    typedef typename FaceDiscreteSpaceType :: BaseFunctionSetType FaceBSetType  ; 
    typedef typename FaceDiscreteSpaceType :: RangeType FaceRangeType; 
    FaceRangeType faceVal;
    
    typedef typename ElementDiscreteSpaceType :: BaseFunctionSetType ElementBaseSetType  ; 
    typedef typename ElementGradientSpaceType :: BaseFunctionSetType GradientBaseSetType  ; 
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
    
    Iterator start = space.begin();
    // if empty grid, do nothing 
    if( start == space.end() ) return; 

    // only working for spaces with one element type 
    assert( space.multipleGeometryTypes() == false );

    // colums are dofs of searched function 
    LocalFuncType lf = uDG.localFunction(*start); 
    const int numDofs = lf.numDofs();

    //const ElementBaseSetType & elSet = elSpace.baseFunctionSet(*start);
    //const int numBubbleDofs = elSet.numBaseFunctions();
  
    const GradientBaseSetType & gradSet = gradSpace.baseFunctionSet(*start);
    const int numGradDofs = gradSet.numBaseFunctions();
  
#if SPASTGRID 
    GeometryType geoType (GeometryType::cube,dim-1);
    const FaceBSetType & faceSet =
      faceSpace.subBaseFunctionSet( geoType, true ); 
#else 
    const FaceBSetType & faceSet =
      faceSpace.baseFunctionSet( *((*start).template entity<1> (0)) );
#endif
    
    const int numFaceDofs = faceSet.numBaseFunctions();
    
    const ReferenceElement< coordType, dim > & refElem =
        ReferenceElements< coordType, dim >::general(start->geometry().type());

    // get number of faces 
    const int overallFaceDofs = numFaceDofs * refElem.size(1);

    //const int rows = (space.order() <= 1) ? (overallFaceDofs) : (overallFaceDofs + numGradDofs + numBubbleDofs);
    const int rows = (space.order() <= 1) ? (overallFaceDofs) : (overallFaceDofs + numGradDofs);
    
    // number of columns 
    //const int cols = numDofs; 

    //std::cout << numFaceDofs << "numFaceDofs \n";
    //std::cout << space.order() << " Order \n";
    //std::cout << overallFaceDofs << " faceDofs \n";
    //std::cout << numGradDofs << " gDofs \n";
    //std::cout << numDofs << " columns \n";

    //std::cout << numBubbleDofs << " bDofs \n";
    // should be equal 
    //std::cout << rows << " r | c " << cols << "\n";
    //if( rows != cols )
    //assert( (rows == cols) ? 1 : (std::cout << rows << " r | c " << cols <<"\n",0) );
    /*
    if( rows != cols ) 
    {
      std::cout << rows << " r | c " << cols <<" in: " << __FILE__ << "\n";
      abort();
    }
    */

    //std::cout << numDofs << " cols " << fD << " numFace \n";

    RangeType ret (0.0);
    RangeType neighRet (0.0);
    RangeType uPhi (0.0);

    // matrix type 
    //typedef DenseMatrix<double> MatrixType; 
    typedef FieldMatrix<RangeFieldType,localBlockSize,localBlockSize> MatrixType;
    //MatrixType matrix(rows,cols);
    MatrixType matrix;
      
    // matrix type 
    //MatrixType fakeMatrix(cols,cols);
      
    /*
    std::vector<double> rhs(rows,0.0);
    std::vector<double> fakeRhs(numDofs,0.0);
    std::vector<double> x(numDofs,0.0);
    */
    typedef FieldVector<RangeFieldType,localBlockSize> VectorType; 
    VectorType rhs(0.0);
    VectorType x(0.0);

    // Get quadraturae rule
    typedef CachingQuadrature <GridPartType , 1> FaceQuadratureType; 

    /*
    DuneODE::SolverInterfaceImpl<MatrixType> op(matrix,numDofs);
    DuneODE::SolverInterfaceImpl<MatrixType> fakeOp(fakeMatrix,numDofs);
    DuneODE::BICGSTAB solver(DuneODE::Communicator::instance());
    
    double epsilon = 1e-16;
    solver.set_max_number_of_iterations(2*numDofs);
    solver.set_tolerance(epsilon);
    */

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // reset values 
      matrix = 0.0;
      rhs = 0.0;
      
      EntityType & en = *it;
      LocalFuncType veloLF = velo.localFunction(en);

      const LocalFuncType uLF = uDG.localFunction(en);
      const BaseFunctionSetType & bSet = uLF.baseFunctionSet(); 
      
      IntersectionIteratorType endnit = gridPart.iend(en);
      for(IntersectionIteratorType nit = gridPart.ibegin(en);
          nit != endnit; ++nit)
      {
        FaceQuadratureType faceQuadInner(gridPart, nit, polOrd, FaceQuadratureType::INSIDE);
        const FaceBSetType & faceSet =
          faceSpace.subBaseFunctionSet(nit.intersectionGlobal());
       
        int firstRow = nit.numberInSelf() * numFaceDofs;
        
        // only interior faces are considered 
        if(nit.neighbor())
        {
          EntityPointerType neighEp = nit.outside();
          EntityType&            nb = *neighEp;
    
          {
            FaceQuadratureType faceQuadOuter(gridPart, nit, polOrd, FaceQuadratureType::OUTSIDE);
            const LocalFuncType uNeighLf = uDG.localFunction(nb);
            
            const int quadNop = faceQuadInner.nop();
            for (int l = 0; l < quadNop ; ++l)
            {
              DomainType unitNormal = 
                nit.integrationOuterNormal(faceQuadInner.localPoint(l));

              const double faceVol = unitNormal.two_norm();
              unitNormal *= 1.0/faceVol;
                           
              const double intel = faceVol * faceQuadInner.weight(l);

              uLF.evaluate(faceQuadInner,l,ret);
              uNeighLf.evaluate(faceQuadOuter,l,neighRet);

              // evaluate uDG 
              ret += neighRet;
              ret *= 0.5;

              double val = ret * unitNormal; 
              val *= intel;
             
              for(int j=0; j<numFaceDofs; ++j)
              {
                int row = firstRow + j; 
                faceSet.eval(j,faceQuadInner.localPoint(l), faceVal);
                rhs[row] += val*faceVal[0];

                for(int i=0; i<numDofs; ++i) 
                {
                  bSet.eval(i,faceQuadInner,l,uPhi); 
                  double prod = uPhi * unitNormal;
                  prod *= intel * faceVal[0];
                  matrix[row][i] += prod;      
                }
              }
            }
          }
        }
       
        // only interior faces are considered 
        if(nit.boundary())
        {
          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l)
          {
            DomainType unitNormal = 
              nit.integrationOuterNormal(faceQuadInner.localPoint(l));

            const double faceVol = unitNormal.two_norm();
            unitNormal *= 1.0/faceVol;
                         
            const double intel = faceVol * faceQuadInner.weight(l);

            uLF.evaluate(faceQuadInner,l,ret);

            double val = ret * unitNormal; 
            val *= intel;
           
            for(int j=0; j<numFaceDofs; ++j)
            {
              int row = firstRow + j; 
              faceSet.eval(j,faceQuadInner.localPoint(l), faceVal);
              rhs[row] += val*faceVal[0];

              for(int i=0; i<numDofs; ++i) 
              {
                bSet.eval(i,faceQuadInner,l,uPhi); 
                double prod = uPhi * unitNormal;
                prod *= intel * faceVal[0];
                matrix[row][i] += prod;      
              }
            }
          }
        }
      }

      // if we have element parts 
      if( overallFaceDofs < rows ) 
      {
        // add part of gradient 
        int gradStartRow = overallFaceDofs;
        gradientPart(gradSpace, en, uLF, gradStartRow, matrix,rhs);
        
        // add part of bubble space 
        int bubbleStartRow = numFaceDofs + numGradDofs; 
        bubblePart(elSpace,en,uLF,bubbleStartRow, matrix,rhs);
      }
      
      /*
      if( rows != cols )
      {
        matrix.multTransposed(rhs,fakeRhs);
        fakeMatrix.multiply_AT_A(matrix);
        
        solver.solveOEM(fakeOp, &x[0], &fakeRhs[0] ); 
      }
      else 
      {
        solver.solveOEM(op, &x[0], &rhs[0] ); 
      }
      */
      matrix.invert();

      // apply M^-1 * rhs 
      x = 0.0;
      matrix.umv(rhs,x);
      
      for(int i=0; i<numDofs; ++i)
      {
        veloLF[i] = x[i];
      }
    }
  }
};

} // end namespace 
#endif


/**************************************************************************
**       Title: examples of elliptic models 
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: implementation of default elliptic models for 
**              use in solvers such as FEOp. A model is a class providing
**              all data functions of a problem. The model is assumed to
**              be instantiated once and passed to solvers. It is depending
**              on a traits class providing typedefs.
**
**************************************************************************/
#ifndef DUNE_ELLIPTICMODEL_HH
#define DUNE_ELLIPTICMODEL_HH

#include <config.h>
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/fem/operator/elementintegratortraits.hh>

namespace Dune
{
  
/*======================================================================*/
/*!
 *  \class PoissonModel
 *  \brief The PoissonModel class provides a default model
 *         for an elliptic problem to be handled by FEOp
 *
 *  The model problem simply is  - div grad u = 2(x(1-x)+ y(1-y))
 *  Homogeneous Dirichlet boundary values. If the problem is solved on
 *  the unit square, the solution should be u = x(1-x)y(1-y)
 *
 *  All types are extracted from the DefaultElementMatrixTraits class
 *  additionally, the model contains member variables and methods.
 */
/*======================================================================*/

class PoissonModel
{
public:  
  typedef DefaultElementIntegratorTraits TraitsType;
  typedef TraitsType::EntityType EntityType;
  typedef TraitsType::ElementQuadratureType ElementQuadratureType;
  typedef TraitsType::IntersectionQuadratureType 
                   IntersectionQuadratureType;
  typedef TraitsType::RangeType RangeType;
  typedef TraitsType::DomainType DomainType;
  typedef TraitsType::JacobianRangeType JacobianRangeType;
  typedef TraitsType::DiscreteFunctionSpaceType 
                   DiscreteFunctionSpaceType;
  
//! constructor with functionspace argument such that the space and the 
//! grid is available
  PoissonModel(DiscreteFunctionSpaceType& fuspace)
            : fuspace_(fuspace) 
        {
          // currently implementation is fitted to 2 dimensions
          assert(dimworld==2);
        }    
  
//! give access to the functionspace
  DiscreteFunctionSpaceType& discreteFunctionSpace() 
        {
          return fuspace_;
        }   

//! return boundary type of a boundary point p used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline TraitsType::BoundaryType boundaryType(EntityType& en, 
                                               QuadratureType& quad, int p)
        {
#ifdef FORCENEUMANN
#warning "POISSONMODEL HAS ONLY NEUMANN-BOUNDARY!"
          return TraitsType::Neumann;
#else          
          return TraitsType::Dirichlet;
#endif // FORCENEUMANN         
        }

//! determine dirichlet value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void dirichletValues(EntityType& en, QuadratureType& quad, int p, 
                              RangeType& ret)
        {
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
          //ret[0] = glob[0]+glob[1]; // 1.0;
          //ret[0] = 0.0;
          ret[0] = 1.0;
          for(int i=0; i<DomainType::dimension; i++)
              ret[0] *= ( glob[i] - SQR(glob[i]) );
    //    ret[0] += x[0]+x[1]; 
    // ret += 1.0;   // add dirichlet-values!!

        }

//! determine neumann value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void neumannValues(EntityType& en, QuadratureType& quad, int p, 
                            RangeType& ret)
        {
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
          //ret[0] = glob[0]+glob[1]; // 1.0;
          //ret[0] = 0.0;
          ret[0] = 0.0;
          for(int i=0; i<DomainType::dimension; i++)
              ret[0] += - ( glob[i] - SQR(glob[i]) );
        }

//! determine robin value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void robinValues(EntityType& en, QuadratureType& quad, int p, 
                          RangeType& ret)
        {
          ret[0] = 0.0;
        }

//! determine mass (i.e. value of the function c) in a domain point used in a 
//! quadrature
  template <class EntityType, class QuadratureType>  
  inline void mass(EntityType& en, QuadratureType& quad, int p, RangeType& ret)
        {
          ret[0] = 0.0;
        }

//! determine source (i.e. value of the function f) in a domain point used in 
//! a quadrature. Here simply f(x,y) = 2 (x(1-x) + y(y-x))
  template <class EntityType, class QuadratureType>  
  inline void source(EntityType& en, QuadratureType& quad, int p, 
                     RangeType& ret)
        {
          const DomainType& glob = en.geometry().global(quad.point(p));          
          ret[0] = 2*(glob[0]*(1-glob[0]) + glob[1]*(1-glob[1]));
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! diffflux = stiffness * grad(phi) 
  template <class EntityType, class QuadratureType>  
  inline void diffusiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   JacobianRangeType& gradphi, 
                   JacobianRangeType& ret)
        {
          // - laplace phi = -div (1 * grad phi - 0 * phi)
          ret = gradphi;          
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! convectiveFlux =  - velocity * phi 
  template <class EntityType, class QuadratureType>  
  inline void convectiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   RangeType& phi, 
                   DomainType& ret)
        {
          // - laplace phi = -div (1 * grad phi - 0 * phi)
          ret = 0.0; // should fill whole vector!!          
        }

  //! the coefficient for robin boundary condition
  inline double alpha()
        { 
          return 1.0;
        }

private:
  DiscreteFunctionSpaceType& fuspace_;  
};  // end of PoissonModel class


/*======================================================================*/
/*!
 *  \class PoissonExactSolution
 *  \brief The class provides the exact solution for the model given by 
 *         the PoissonModel class
 *
 *  The function represents u = x(1-x)y(1-y), which is the solution of
 *  the model problem  - div grad u = 2(x(1-x)+ y(1-y)) on
 *  the unit square with homogeneous Dirichlet boundary values. 
 *  The function can be used for EOC calculation
 */
/*======================================================================*/

  
class PoissonExactSolution : 
    public Function < DefaultElementIntegratorTraits::FunctionSpaceType , PoissonExactSolution > 
{
  typedef DefaultElementIntegratorTraits::FunctionSpaceType FunctionSpaceType;
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;
public:
  PoissonExactSolution (FunctionSpaceType &f) : 
          Function < FunctionSpaceType , PoissonExactSolution > ( f ) {}
 
  //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
  void evaluate (const DomainType & x , RangeType & ret) const
  {
    ret = 1.0;
    for(int i=0; i<DomainType::dimension; i++)
        ret *= ( x[i] - SQR(x[i]) );
    //    ret[0] += x[0]+x[1]; 
    // ret += 1.0;   // add dirichlet-values!!
  }

  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};


/*======================================================================*/
/*!
 *  \class Elliptic2dModel
 *  \brief The Elliptic2dModel class provides a complete model
 *         for an elliptic problem to be handled by FEOp
 *
 *  The model problem is defined on the 2D unit square with
 *        - div ( a grad u - b u) + cu = f
 *            u = g_D on Dirichlet-boundary (upper and lower edge)
 *      (a grad u - bu ) n = g_N on Neuman boundary (left edge)
 *      (a grad u - bu ) n + alpha u= g_R on Robin boundary (right edge)
 *
 *  The data functions are parametrized by nonnegative scalar q,r,s for
 *  switching on/off certain contributions and regulating the stability
 *     alpha = 1
 *     a = [(1+q), -q; -q,  (1+q)]
 *     b = [1; 1] * s * y
 *     c = xy * r
 *     f = 2q + r* x^2y^2 + r* x^2y + s * (x+y+y^2 + 2xy) 
 *     g_D = xy+x
 *     g_N = -(1+q)(y+1)
 *     g_R = 2 - s*y^2 + (2+q-s)* y
 *
 *  The solution is simply u = xy+x
 *
 *  All types are extracted from the DefaultElementMatrixTraits class
 *  additionally, the model contains member variables and methods.
 */
/*======================================================================*/

class Elliptic2dModel
{
public:  
  typedef DefaultElementIntegratorTraits TraitsType;
  typedef TraitsType::EntityType EntityType;
  typedef TraitsType::ElementQuadratureType ElementQuadratureType;
  typedef TraitsType::IntersectionQuadratureType 
                   IntersectionQuadratureType;
  typedef TraitsType::RangeType RangeType;
  typedef TraitsType::DomainType DomainType;
  typedef TraitsType::JacobianRangeType JacobianRangeType;
  typedef TraitsType::DiscreteFunctionSpaceType 
                   DiscreteFunctionSpaceType;

  static const double eps = 1e-10;
  static const double q = 1.0; // q>0 => non-unit diffusivity
  static const double r = 1.0; // r>0 => mass term activated 
  static const double s = 1.0; // 0.001; // s>0 => convective term activated
  
//  typedef typename TraitsType::EntityType EntityType;
//  typedef typename TraitsType::ElementQuadratureType ElementQuadratureType;
//  typedef typename TraitsType::IntersectionQuadratureType 
//                   IntersectionQuadratureType;
//  typedef typename TraitsType::RangeType RangeType;
//  typedef typename TraitsType::DomainType DomainType;
//  typedef typename TraitsType::JacobianRangeType JacobianRangeType;
//  typedef typename TraitsType::DiscreteFunctionSpaceType 
//                   DiscreteFunctionSpaceType;
  
//! constructor with functionspace argument such that the space and the 
//! grid is available
  Elliptic2dModel(DiscreteFunctionSpaceType& fuspace)
            : fuspace_(fuspace) 
        {
          // currently implementation is fitted to 2 dimensions
          assert(dimworld==2);
        }    
  
//! give access to the functionspace
  DiscreteFunctionSpaceType& discreteFunctionSpace() 
        {
          return fuspace_;
        }   

//! return boundary type of a boundary point p used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline TraitsType::BoundaryType boundaryType(EntityType& en, 
					       QuadratureType& quad, int p)
        {
#ifdef FORCENEUMANN
#warning "ELLIPTIC2DMODEL HAS ONLY NEUMANN-BOUNDARY!"
          return TraitsType::Neumann;	    
#else
#if 1 // Full problem: Dirichlet, Neuman and Dirichlet-boundary
//#if 0 // temporary check: only Dirichlet Bnd
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
	  if (glob[0]<eps)
	    return TraitsType::Neumann;
	  else if (glob[0]>1-eps)
	    return TraitsType::Robin;
	  else 
#endif 
	    return TraitsType::Dirichlet;	    
#endif // FORCENEUMANN
        }

//! determine dirichlet value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void dirichletValues(EntityType& en, QuadratureType& quad, int p, 
                              RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
    ret[0] = glob[0] * ( 1.0 + glob[1]);
  }
  
//! determine neumann value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void neumannValues(EntityType& en, QuadratureType& quad, int p, 
                            RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
//    assert(glob[0]<eps);
    ret[0] = - (1+q) * (glob[1] +1);
  }

//! determine robin value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void robinValues(EntityType& en, QuadratureType& quad, int p, 
                          RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
    assert(glob[0]>1.0-eps);
    ret[0] = 2.0 - s * SQR(glob[1]) + (2+q-s) * glob[1];
  }
  
//! determine mass (i.e. value of the function c) in a domain point used in a 
//! quadrature
  template <class EntityType, class QuadratureType>  
  inline void mass(EntityType& en, QuadratureType& quad, int p, RangeType& ret)
        {
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
	  ret[0] = glob[0]*glob[1] * r;
        }

//! determine source (i.e. value of the function f) in a domain point used in 
//! a quadrature. 
  template <class EntityType, class QuadratureType>  
  inline void source(EntityType& en, QuadratureType& quad, int p, 
                     RangeType& ret)
        {
          const DomainType& glob = en.geometry().global(quad.point(p));     
          ret[0] = 2.0 * q
              + r * SQR(glob[0])* SQR(glob[1])  
              + r * SQR(glob[0]) * glob[1]
              + s * ( glob[1] + glob[0] 
                      + SQR(glob[1]) + 2 * glob[0]*glob[1]);         
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! diffflux = stiffness * grad(phi) 
  template <class EntityType, class QuadratureType>  
  inline void diffusiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   JacobianRangeType& gradphi, 
                   JacobianRangeType& ret)
        {
//          ret = gradphi;
          ret[0][0] = (1+q) * gradphi[0][0] - q    * gradphi[0][1];
          // the following cumbersome multiplication with 1.0 is required, 
          // otherwise compile error: missing reference to 
          //Dune::Elliptic2dModel::q
          ret[0][1] = - 1.0 * q * gradphi[0][0]  + (1+q) * gradphi[0][1];
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! convectiveFlux =  - velocity * phi 
  template <class EntityType, class QuadratureType>  
  inline void convectiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   RangeType& phi, 
                   DomainType& ret)
        {
          const DomainType& glob = en.geometry().global(quad.point(p));     
	  ret[0] = - glob[1] * s * phi[0];
	  ret[1] = - glob[1] * s * phi[0];
        }

  //! the coefficient for robin boundary condition
  inline double alpha()
        { 
          return 1.0;
        }

private:
  DiscreteFunctionSpaceType& fuspace_;  
};  // end of Elliptic2dModel class

/*======================================================================*/
/*!
 *  \class Elliptic2dExactSolution
 *  \brief The class provides the exact solution for the model given by 
 *         the Elliptic2dModel class
 *
 *  The function represents u = x y + x, which is the solution of
 *  the model problem on the unit square with inhomogeneous Dirichlet 
 *  Neumann and Dirichlet boundary values. 
 *  Function can be used for EOC calculation
 */
/*======================================================================*/

class Elliptic2dExactSolution : 
    public Function < DefaultElementIntegratorTraits::FunctionSpaceType , Elliptic2dExactSolution > 
{
  typedef DefaultElementIntegratorTraits::FunctionSpaceType FunctionSpaceType;
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;
public:
  Elliptic2dExactSolution (FunctionSpaceType &f) : 
          Function < FunctionSpaceType , Elliptic2dExactSolution > ( f ) {}
 
  void evaluate (const DomainType & x , RangeType & ret) const
  {
    ret = x[0]*x[1] + x[0];
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
}; // end class Elliptic2dExactSolution

/*======================================================================*/
/*!
 *  \class Elliptic3dModel
 *  \brief The Elliptic3dModel class provides a complete model
 *         for an elliptic problem to be handled by FEOp
 *
 *  The model problem is defined on the 3D unit cube with
 *        - div ( a grad u - b u) + cu = f
 *      (a grad u - bu ) n = g_N on Neuman boundary (x=0 face)
 *      (a grad u - bu ) n + alpha u= g_R on Robin boundary (x=1 face)
 *            u = g_D on Dirichlet-boundary (remaining 4 faces)
 *
 *  The data functions are
 *     alpha = 1
 *     a = [3 -1 -1; -1 3 -1; -1 -1 3]
 *     b = [1; 1; 1] y
 *     c = xy
 *     f = 2z+3y+3x+y^2z+2xyz+xy^2+x^2y^2z+x^2y
 *     g_D = xyz+x
 *     g_N = -3yz-3
 *     g_R = 4yz-2y-z+4-y^2z
 *
 *  The solution is simply u = xyz+x
 *
 *  All types are extracted from the DefaultElementMatrixTraits class
 *  additionally, the model contains member variables and methods.
 */
/*======================================================================*/

class Elliptic3dModel
{
public:  
  typedef DefaultElementIntegratorTraits TraitsType;
  typedef TraitsType::EntityType EntityType;
  typedef TraitsType::ElementQuadratureType ElementQuadratureType;
  typedef TraitsType::IntersectionQuadratureType 
                   IntersectionQuadratureType;
  typedef TraitsType::RangeType RangeType;
  typedef TraitsType::DomainType DomainType;
  typedef TraitsType::JacobianRangeType JacobianRangeType;
  typedef TraitsType::DiscreteFunctionSpaceType 
                   DiscreteFunctionSpaceType;

  static const  double eps = 1e-10;

//  typedef typename TraitsType::EntityType EntityType;
//  typedef typename TraitsType::ElementQuadratureType ElementQuadratureType;
//  typedef typename TraitsType::IntersectionQuadratureType 
//                   IntersectionQuadratureType;
//  typedef typename TraitsType::RangeType RangeType;
//  typedef typename TraitsType::DomainType DomainType;
//  typedef typename TraitsType::JacobianRangeType JacobianRangeType;
//  typedef typename TraitsType::DiscreteFunctionSpaceType 
//                   DiscreteFunctionSpaceType;
  
//! constructor with functionspace argument such that the space and the 
//! grid is available
  Elliptic3dModel(DiscreteFunctionSpaceType& fuspace)
            : fuspace_(fuspace) 
        {
          // currently implementation is fitted to 3 dimensions
          assert(dimworld==3);
        }    
  
//! give access to the functionspace
  DiscreteFunctionSpaceType& discreteFunctionSpace() 
        {
          return fuspace_;
        }   

//! return boundary type of a boundary point p used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline TraitsType::BoundaryType boundaryType(EntityType& en, 
					       QuadratureType& quad, int p)
        {
#ifdef FORCENEUMANN
#warning "ELLIPTIC3DMODEL HAS ONLY NEUMANN-BOUNDARY!"
          return TraitsType::Neumann;
#else
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
	  if (glob[0]<eps)
	    return TraitsType::Neumann;
	  else if (glob[0]>1-eps)
	    return TraitsType::Robin;
	  else 
	    return TraitsType::Dirichlet;	    
#endif // FORCENEUMANN
        }

//! determine dirichlet value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void dirichletValues(EntityType& en, QuadratureType& quad, int p, 
                              RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
    ret[0] = glob[0] * ( 1.0 + glob[1]*glob[2]);
  }
  
//! determine neumann value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void neumannValues(EntityType& en, QuadratureType& quad, int p, 
                            RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
    assert(glob[0]<eps);
    ret[0] = -3.0 * glob[1]*glob[2] - 3.0;
  }

//! determine robin value in a boundary point used in a quadrature
  template <class EntityType, class QuadratureType>  
  inline void robinValues(EntityType& en, QuadratureType& quad, int p, 
                          RangeType& ret)
  {
    const DomainType& glob = en.geometry().global(quad.point(p)); 
    assert(glob[0]>1.0-eps);
    ret[0] = 4 * glob[1] * glob[2] - 2* glob[1] - glob[2] + 
        4.0 - SQR(glob[1])*glob[2];
  }
  
//! determine mass (i.e. value of the function c) in a domain point used in a 
//! quadrature
  template <class EntityType, class QuadratureType>  
  inline void mass(EntityType& en, QuadratureType& quad, int p, RangeType& ret)
        {
	  const DomainType& glob = en.geometry().global(quad.point(p)); 
	  ret[0] = glob[0]*glob[1];
        }

//! determine source (i.e. value of the function f) in a domain point used in 
//! a quadrature. 
  template <class EntityType, class QuadratureType>  
  inline void source(EntityType& en, QuadratureType& quad, int p, 
                     RangeType& ret)
        {
          const DomainType& glob = en.geometry().global(quad.point(p));     
          ret[0] = 2 * glob[2] + 3* glob[1] + 3 * glob[0] + 
              SQR(glob[1])*glob[2] + 2* glob[0] * glob[1]* glob[2] + 
              glob[0] * SQR(glob[1]) + SQR(glob[0]*glob[1]) * glob[2] + 
              SQR(glob[0]) * glob[1];
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! diffflux = stiffness * grad(phi) 
  template <class EntityType, class QuadratureType>  
  inline void diffusiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   JacobianRangeType& gradphi, 
                   JacobianRangeType& ret)
        {
          ret[0][0] =   3* gradphi[0][0] - gradphi[0][1]   - gradphi[0][2];
          ret[0][1] =   - gradphi[0][0]  +3* gradphi[0][1] - gradphi[0][2];
          ret[0][2] =   - gradphi[0][0] - gradphi[0][1] + 3* gradphi[0][2];
        }

//! no direct access to stiffness and velocity, but whole flux, i.e.
//! convectiveFlux =  - velocity * phi 
  template <class EntityType, class QuadratureType>  
  inline void convectiveFlux(EntityType& en, QuadratureType& quad, int p, 
                   RangeType& phi, 
                   DomainType& ret)
        {
          const DomainType& glob = en.geometry().global(quad.point(p));     
	  ret[0] = - glob[1]*phi[0];
	  ret[1] = - glob[1]*phi[0];
	  ret[2] = - glob[1]*phi[0];
        }

  //! the coefficient for robin boundary condition
  inline double alpha()
        { 
          return 1.0;
        }

private:
  DiscreteFunctionSpaceType& fuspace_;  
};  // end of Elliptic3dModel class

/*======================================================================*/
/*!
 *  \class Elliptic3dExactSolution
 *  \brief The class provides the exact solution for the model given by 
 *         the Elliptic2dModel class
 *
 *  The function represents u = x y z + x, which is the solution of
 *  the model problem on the unit cube with inhomogeneous Dirichlet 
 *  Neumann and Dirichlet boundary values. 
 *  Function can be used for EOC calculation
 */
/*======================================================================*/

class Elliptic3dExactSolution : 
    public Function < DefaultElementIntegratorTraits::FunctionSpaceType , Elliptic3dExactSolution > 
{
  typedef DefaultElementIntegratorTraits::FunctionSpaceType FunctionSpaceType;
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;
public:
  Elliptic3dExactSolution (FunctionSpaceType &f) : 
          Function < FunctionSpaceType , Elliptic3dExactSolution > ( f ) {}
 
  void evaluate (const DomainType & x , RangeType & ret) const
  {
    ret = x[0]*x[1]*x[2] + x[0];
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
}; // end class Elliptic3dExactSolution

}; // end namespace Dune

#endif

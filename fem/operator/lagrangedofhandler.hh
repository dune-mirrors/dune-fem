/**************************************************************************
**       Title: Class collecting and giving access to information about 
**              LagrangePoints in a Lagrange basis
**              A corresponding Interpolation is provided
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: Class, which is only a beginning, but can be developed further
**
**************************************************************************/

#ifndef DUNE_LAGRANGEDOFHANDLER_HH
#define DUNE_LAGRANGEDOFHANDLER_HH

#include <config.h>
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/fem/misc/entityfunction.hh>

namespace Dune 
{

/*======================================================================*/
/*!
 *  \class LagrangeDofHandler
 *  \brief The LagrangeDofHandler class returns local coordinates of a 
 *         Lagrange-point by a point() method.
 *
 *  This is an adapter, initialized with the corresponding space and entity, 
 *  which returns the local coordinates of a Lagrange node by a point() 
 *  method. It provides access to DOF-numbers located on faces of the entity.
 *
 *  The point() method is required, as the data functions in problem models 
 *  only allow access in a "quadrature-style" call, e.g. 
 *
 *      model.dirichletValues(entity, quadrature, point, ret)
 *
 *  So in particular evaluations in Dirichlet-Dofs require a "quadrature"-like 
 *  access wrapper, so quadrature in this call can be replaces by in 
 *  instance of the lagrangeDofHandler. Currently this class is used in 
 *  elementintegrators.hh and feop.hh
 *
 *  Currently it is only implemented trivially for order 1 elements assuming
 *  a correspondence between DOFs and element vertices, in particular 
 *  identical local enumeration. 
 */
/*======================================================================*/

template <class DiscreteFunctionSpaceType>
class LagrangeDofHandler
{ 
public:
  typedef typename 
      DiscreteFunctionSpaceType::GridType::template Codim<0>::Entity 
      EntityType;
  typedef typename EntityType::Geometry GeometryImp;
  typedef typename EntityType::ctype coordType;
  
  enum { dim = EntityType :: dimension };
  
/*======================================================================*/
/*! 
 *   constructor: initialization with fspace and entity
 *
 *   \param fspace the lagrange-function space for which DOFs are handled
 *
 *   \param entity the entity for which the handler operates
 */
/*======================================================================*/

  LagrangeDofHandler(const DiscreteFunctionSpaceType& fspace, 
                     const EntityType& entity):
          fspace_(fspace),
          entity_(entity), 
          geo_(entity_.geometry()), 
          geotype_(geo_.type()),
          refElem_(refElemCont_(geotype_))
        {
          // currently only implementation for order-1 base functions is
          // realized with correspondence of grid-points to lagrange-points
          assert(fspace.polynomOrder() == 1);
//          geo_= entity_.geometry();
//          geotype_ = geo_.type();
//          refElem_ = refElemCont_(geotype_);
        };
  
/*======================================================================*/
/*! 
 *   point: point access in a quadrature-style fashion without 
 *          specifying the entity!
 *
 *   \param p integer specifying the Lagrange-Point number within the entity 
 *            the number is assumed to be identical with local-function 
 *            dof index.
 *
 *   \return the local coordinates of the Lagrange-Point
 */
/*======================================================================*/

//  const FieldVector<typename EntityType::ctype, EntityType::dimension > 
  const FieldVector<typename EntityType::ctype, 
                    EntityType::Geometry::mydimension > 
  point(int p) const
        {
          //const  FieldVector<typename EntityType::ctype, 
          //    EntityType::Geometry::mydimension > 
          //    ret =  geo_.local(geo_[p]);
          //std::cout << "original point: " <<  geo_.local(geo_[p]);
          //std::cout << "new point     : " <<  ret;
          return geo_.local(geo_[p]);
        };

/*======================================================================*/
/*! 
 *   numDofsOnFace: determine number of Lagrange-nodes on a face.
 *
 *   the method can be used for determining basis functions, which have a 
 *   support on the face.
 *
 *   \param faceNum local number of face in entity
 *
 *   \param faceCodim codimension of face in entity
 *
 *   \return number of Lagrange-Nodes on the specified face
 */
/*======================================================================*/

  int numDofsOnFace(int faceNum, int faceCodim) const
        {
          return refElem_.size( faceNum, faceCodim , dim );
        };
  
/*======================================================================*/
/*! 
 *   entityDofNum: determine DOF number in entity for given DOF number on face 
 *
 *   determine DOF number in entity for given DOF number on face 
 *   By iterating over 0<=faceDofNum<numDofsOnFace(...) all number of 
 *   base functions with support on the face are accessed.
 *
 *   \param faceNum local number of face in entity
 *   \param faceCodim codimension of face in entity
 *   \param faceDofNum number of DOF in face (0...numDofsOnFace()) 
 *
 *   \return the DOF number in the entity, e.g. for indexing localFunctions.
 */
/*======================================================================*/

  int entityDofNum(int faceNum, int faceCodim , int faceDofNum ) const
        {
          return refElem_.subEntity(faceNum, faceCodim , 
                                    faceDofNum , dim);
        };
  
private:

  //! reference to the (Lagrange-)functionspace used in initialization
  const DiscreteFunctionSpaceType& fspace_;  
  //! reference to the entity used in initialization
  const EntityType& entity_;  
  //! the referencecontainer for the element
  ReferenceElementContainer<coordType,dim> refElemCont_;
  //! The geometry of the element
  const GeometryImp& geo_;
  //! the geometrytype of the element
  const GeometryType& geotype_;
  //! the referenceelement for the entity
  const ReferenceElement<coordType,dim>& refElem_;  

}; // end class LagrangeDofHandler

/*======================================================================*/
/*!
 *  \class LagrangeInterpolator 
 *  \brief The LagrangeInterpolator class performs lagrange
 *         interpolation of an analytical function (Dune::Function) 
 *
 *  could also be implemented as a singleton, as only functionality is 
 *  provided by a template method
 */
/*======================================================================*/
  
  class LagrangeInterpolator
  {
    public:
    //! constructor
    LagrangeInterpolator()
        {
        }
    
    //! destructor
    ~LagrangeInterpolator()
        {
        }
    
    //! interpolation routine
/*======================================================================*/
/*! 
 *   interpolFunction: interpolate analytical function into Lagrange Discrete
 *                     function
 *
 *   This method can be used, if a function is given, which can be evaluated
 *   globally by an evaluate(glob, ret) method. In constrast: Functions, which
 *   are based on discrete functions cannot be used by this method, but must
 *   use interpolateEntityFunction.
 *
 *   nothing else than calling interpolateEntityFunction with suitable
 *   Wrapper around FunctionType is performed
 *
 *   \param f analytical function (derived from Dune::Function)
 *
 *   \param discFunc reference to discrete lagrange-function
 */
/*======================================================================*/

    template <class FunctionType, class DiscreteFunctionType> 
        void interpolFunction (FunctionType &f, 
                               DiscreteFunctionType &discFunc)
          {
            // these types are not available!!
            //                typename DiscreteFunctionType::GridType::EntityType, 
            //                typename DiscreteFunctionType::GridPartType::EntityType, 
            //            typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType::EntityType, 
            typedef typename DiscreteFunctionType::GridType::template Codim<0>::Entity   EntityType;
            EntityFunctionAdapter<
                EntityType,
                FunctionType> enf(f);
            interpolEntityFunction 
                (enf, discFunc);            
          };
    
/*======================================================================*/
/*! 
 *   interpolEntityFunction: interpolate entity-function into Lagrange 
 *                 Discrete function
 *
 *   interpolation of EntityFunctions, which cannot be evaluated globally, but 
 *   locally. By this, functions, which depend on discretefunction can also be 
 *   interpolated, etc.
 *
 *   \param f an the entityfunction to be interpolated 
 *
 *   \param discFunc reference to the generated discrete lagrange-function
 */
/*======================================================================*/
    
    template < class EntityFunctionType, 
              class DiscreteFunctionType> 
    void interpolEntityFunction (EntityFunctionType &f, 
                                 DiscreteFunctionType &discFunc)
          {
            typedef typename DiscreteFunctionType::FunctionSpaceType 
                DiscreteFunctionSpaceType;
            
            typedef typename DiscreteFunctionSpaceType::IteratorType Iterator;
            typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
            typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
            typedef typename DiscreteFunctionType::LocalFunctionType 
                LocalFunctionType;
            
          const DiscreteFunctionSpaceType
              & functionSpace = discFunc.getFunctionSpace();  
          
          discFunc.clear();      
          
          RangeType ret (0.0);
          
          Iterator it = functionSpace.begin();
          Iterator endit = functionSpace.end();
          
          for( ; it != endit ; ++it)
          {
            LagrangeDofHandler<DiscreteFunctionSpaceType> 
                dofHandler(functionSpace, *it);
            
            // initialize functions
            LocalFunctionType lf = discFunc.localFunction( *it ); 
            f.init(*it);
            
            //DomainType local(0.0);
            for(int i=0; i<lf.numDofs(); i++)
            {
              // f.evaluate(*it, dofHandler, i, ret);
              f.evaluate(dofHandler.point(i), ret);
              //local = dofHandler.point(i);            
              //DomainType glob = (*it).geometry().global(local);            
              //f.evaluate(glob, ret);
              lf[i] = ret[0];
            }
          }
        }

  }; // end class LagrangeInterpolator
  
}; // end namespace dune

#endif

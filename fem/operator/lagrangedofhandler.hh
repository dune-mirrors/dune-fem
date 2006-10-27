/**************************************************************************
**       Title: Class collecting and giving access to information about 
**              LagrangeDofs in a Lagrange basis
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
 *
 *   \return the local coordinates of the Lagrange-Point
 */
/*======================================================================*/

  const FieldVector<typename EntityType::ctype, EntityType::dimension >& 
  point(int p)
        {
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
 *   \param faceCoDim codimension of face in entity
 *
 *   \return number of Lagrange-Nodes on the specified face
 */
/*======================================================================*/

  int numDofsOnFace(int faceNum, int faceCodim)
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
 *
 *   \param faceCoDim codimension of face in entity
 *
 *   \param faceDofNum number of DOF in face (0...numDofsOnFace()) 
 *
 *   \return the DOF number in the entity, e.g. for indexing localFunctions.
 */
/*======================================================================*/

  int entityDofNum(int faceNum, int faceCodim , int faceDofNum )
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

}; // end namespace dune

#endif

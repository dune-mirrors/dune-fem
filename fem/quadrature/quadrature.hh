#ifndef DUNE_FEM_QUADRATURE_HH
#define DUNE_FEM_QUADRATURE_HH

//#include <vector>
#include <cassert>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>

#include <dune/fem/storage/array.hh>

// quadrature storage classes 
#include <dune/fem/quadrature/quadprovider.hh>
#include <dune/fem/quadrature/defaultquadratures.hh>

namespace Dune
{

  /** \addtogroup Quadrature 
   *
   *  In dune-fem, quadratures are a set of quadrature points and corresponding
   *  weights.
   *
   *  \remark To get an overview of a quadrature, see Quadrature.
   *
   *  \remark The quadratures usually used are ElementQuadrature and
   *  CachingQuadrature.
   */


  
  /** \class   QuadraturePointWrapper
   *  \ingroup Quadrature
   *  \brief   wrapper for a (Quadrature,int) pair
   */
  template< class Quadrature >
  class QuadraturePointWrapper
  {
  public:
    typedef Quadrature QuadratureType;

  private:
    typedef QuadraturePointWrapper< QuadratureType > ThisType;

  protected:
    const QuadratureType &quadrature_;
    const unsigned int quadPoint_;

  public:
    inline QuadraturePointWrapper ( const QuadratureType &quadrature,
                                    const unsigned int quadPoint )
    : quadrature_( quadrature ),
      quadPoint_( quadPoint )
    {
    }

  private:
    // forbid assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline const QuadratureType &quadrature () const
    {
      return quadrature_;
    }

    inline unsigned int point () const
    {
      return quadPoint_;
    }
  };



  template< class Point >
  struct QuadraturePointWrapperHelper
  {
    typedef Point PointType;
    typedef Point DomainType;

    inline static const DomainType &point ( const PointType &x )
    {
      return x;
    }
  };



  template< class Quadrature >
  struct QuadraturePointWrapperHelper< QuadraturePointWrapper< Quadrature > >
  {
    typedef Quadrature QuadratureType;
    typedef QuadraturePointWrapper< QuadratureType > PointType;
    typedef typename QuadratureType :: CoordinateType DomainType;

    inline static const DomainType &point ( const PointType &x )
    {
      return x.quadrature().point( x.point() );
    }
  };



  /** \brief   extract the real coordinate from a point
   *  \relates Dune::QuadraturePointWrapper
   *
   *  This function unwraps a possibly \ref Dune::QuadraturePointWrapper
   *  "wrapped quadrature point". If the point is not a wrapped quadrature
   *  point, the point itself is returned. This makes it possible to have
   *  one evaluation method for the two different kinds of evaluation points.
   *
   *  \param[in]  x  possibly wrapped point
   *  
   *  \returns a reference to the actual point
   */
  template< class Point >
  inline const typename QuadraturePointWrapperHelper< Point > :: DomainType &
    coordinate ( const Point &x )
  {
    return QuadraturePointWrapperHelper< Point > :: point( x );
  }



  /** \class IntegrationPointList
   *  \ingroup Quadrature
   *  \brief actual interface class for integration point lists
   *
   *  IntegrationPointList is a proxy for the actual implementations of the
   *  integration point lists. During construction, the IntegrationPointList
   *  object is configured with an appropriate implementation object from the
   *  QuadratureProvider (monostate pattern).
   *
   *  The design goal is minimization of construction time. The actual
   *  implementation can be created once and reused whenever it is needed.
   *  Moreover, this layout insulates the user from all initialization and
   *  storage stuff.
   *
   *  \note The difference between integration point lists and quadratures is
   *        that quadratures have weights.
   */
  template< typename FieldImp, int dim,
            template< class, int > class IntegrationTraits >
  class IntegrationPointList 
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

  private:
    typedef IntegrationPointList< FieldType, dimension, IntegrationTraits > ThisType;

    typedef IntegrationTraits< FieldType, dimension > Traits;

    typedef QuadratureProvider< FieldType, dimension, IntegrationTraits >
      QuadratureProviderType;

  public:
    //! type of integration point list implementation 
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;
    
    //! type of coordinate
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;
    
    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  protected:
    const IntegrationPointListType &ipList_;

  public:
    /** \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  order     order of the requested quadrature
     */
    inline IntegrationPointList ( const GeometryType &geometry,
                                  int order )
    : ipList_( QuadratureProviderType :: getQuadrature( geometry, order ) )
    {
    }

    /** \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  elementGeometry  geometry type of element that resulting 
     *              quadrature is used for (in case of face quadratures)
     *  \param[in]  order     order of the requested quadrature
     */
    inline IntegrationPointList ( const GeometryType &geometry,
                                  const GeometryType &elementGeometry,
                                  int order )
    : ipList_( QuadratureProviderType :: getQuadrature( geometry, elementGeometry, order ) )
    {
    }

    /** \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for testing purposes.
     *
     *  \param[in]  ipList  implementation of the integration point list
     */
    inline IntegrationPointList ( const IntegrationPointListType &ipList )
    : ipList_( ipList )
    {
    }

    /** \brief copy constructor
     *  
     *  \param[in]  org  integration point list to be copied
     */ 
    inline IntegrationPointList ( const IntegrationPointList &org )
    : ipList_( org.ipList_ )
    {
    }

    const QuadraturePointWrapperType operator[] ( unsigned int i ) const
    {
      return QuadraturePointWrapperType( *this, i );
    }

    /** \brief obtain a reference the actual implementation
     * 
     *  \returns a reference to the implementation of this integration point
     *           list
     */
    const IntegrationPointListType &ipList () const
    {
      return ipList_;
    }

    /** \brief obtain the number of integration points
     *
     *  \returns number of integration points within this list
     */
    int nop () const
    {
      return ipList_.nop();
    }
    
    /** \brief obtain coordinates of i-th integration point
     *
     *  This method returns a reference to the coordinates of the i-th
     *  integration point for 0 <= i < nop(). The integration point is given
     *  in local coordinates, i.e., coordinates with respect to the reference
     *  element.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns reference to i-th integration point
     */
    const CoordinateType &point ( size_t i ) const
    {
      return ipList_.point( i );
    }

    /** \brief obtain the identifier of the integration point list
     * 
     *  The identifier of an integration point list must be globally unique.
     *  Even integration point lists for different dimensions must have
     *  different identifiers.
     *
     *  \note Quadratures are considered distinct if they differ in one of the
     *        following points: geometry type, order, dimension or implementation.
     *
     *  \returns globally unique identifier of the integration point list
     */
    size_t id () const
    {
      return ipList_.id();
    }
    
    /** \brief obtain order of the integration point list
     *
     *  The order of a quadrature is the maximal polynomial degree that is
     *  guaranteed to be integrated exactly by the quadrature.
     *
     *  In case of an integration point list, the definition of this value is
     *  left to the implementor.
     *  
     *  \note Calling this method yields a virtual function call, so do not
     *        call this method unnecessarily.
     *
     *  \returns the order of the integration point list
     */
    int order () const
    {
      return ipList_.order();
    }
    
    /** \brief obtain GeometryType for this integration point list
     *
     *  Integration point lists are specified in local coordinates, i.e.,
     *  coordinates with respect to the reference element. Hence, each 
     *  integration point list is only valid for one type of geometry, i.e.,
     *  for one reference element. The type can be retrieved via this method.
     *
     *  \note Calling this method yields a virtual function call, so do not
     *        call this method unnecessarily.
     *
     *  \returns GeometryType for this integration point list
     */
    GeometryType geometry () const
    {
      return ipList_.geometry();
    }
  };



  /** \class Quadrature
   *  \ingroup Quadrature
   *  \brief actual interface class for quadratures
   *
   *  IntegrationPointList is a proxy for the actual implementations of the
   *  integration point lists. During construction, the IntegrationPointList
   *  object is configured with an appropriate implementation object from the
   *  QuadratureProvider (monostate pattern).
   *
   *  The design goal is minimization of construction time. The actual
   *  implementation can be created once and reused whenever it is needed.
   *  Moreover, this layout insulates the user from all initialization and
   *  storage stuff.
   *
   *  \note The difference between integration point lists and quadratures is
   *        that quadratures have weights.
   */
  template< class FieldImp, int dim,
            template< class, int > class QuadratureTraits = DefaultQuadratureTraits >
  class Quadrature
  : public IntegrationPointList< FieldImp, dim, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

  private:
    typedef Quadrature< FieldType, dimension, QuadratureTraits > ThisType;
    typedef IntegrationPointList< FieldType, dimension, QuadratureTraits > BaseType;

    typedef QuadratureTraits< FieldType, dimension > Traits;
    
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits >
      QuadratureProviderType;

  public:
    using BaseType :: ipList;

    //! type of the implementation (this must actually be a quadrature implementation)
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;

    //! type of local coordinate vectors
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! to be revised, look at caching quad 
    enum { codimension = 0 };

  public:
    /** \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  order     order of the requested quadrature
     */
    inline Quadrature( const GeometryType &geometry, int order )
    : BaseType( geometry, order )
    {
    }

    /** \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     * 
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometry  geometry type of the requested quadrature
     *  \param[in]  elementGeometry  geometry type of element that resulting 
     *              quadrature is used for (in case of face quadratures)
     *  \param[in]  order     order of the requested quadrature
     *
     *  \note This is a specialized constructor for constructing 
     *  face quadratures for UGGrid.
     */
    inline Quadrature ( const GeometryType &geometry,
                        const GeometryType &elementGeometry,
                        int order )
    : BaseType( geometry, elementGeometry, order )
    {
    }

    /** \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for testing purposes.
     *
     *  \param[in]  ipList  implementation of the integration point list
     */
    inline explicit Quadrature( const IntegrationPointListType& ipList )
    : BaseType( ipList )
    {
    }

    /** \brief copy constructor
     *  
     *  \param[in]  org  quadrature to be copied
     */ 
   //! Copy constructor
    inline Quadrature( const Quadrature &org )
    : BaseType( org )
    {
    }

    /** \brief obtain weight of i-th integration point
     *
     *  This method returns the weight of the i-th integration point for
     *  0 <= i < nop() within the quadrature.
     *
     *  \note The integration point can be obtained via the point() method.
     *
     *  \note The quadrature weights sum up to the volume of the reference
     *        element.
     *
     *  \param[in]  i  number of the integration point, 0 <= i < nop()
     *
     *  \returns weight of the i-th integration point
     */
    const FieldType &weight( size_t i ) const
    {
      return ipList().weight( i );
    }
  };

}

#endif

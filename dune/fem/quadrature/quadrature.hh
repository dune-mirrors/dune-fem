#ifndef DUNE_FEM_QUADRATURE_HH
#define DUNE_FEM_QUADRATURE_HH

#include <cassert>
#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/common/coordinate.hh>
#include <dune/fem/common/utility.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/version.hh>
#include <dune/fem/quadrature/quadprovider.hh>
#include <dune/fem/quadrature/defaultquadratures.hh>

namespace Dune
{
  namespace Fem
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

    typedef typename QuadratureType::CoordinateType CoordinateType;
    typedef typename QuadratureType::RealType RealType;
    typedef typename QuadratureType::LocalCoordinateType LocalCoordinateType;

  private:
    typedef QuadraturePointWrapper< QuadratureType > ThisType;

  protected:
    const QuadratureType &quadrature_;
    const unsigned int index_;

  public:
    QuadraturePointWrapper ( const QuadratureType &quadrature, unsigned int index )
      : quadrature_( quadrature ), index_( index )
    {}

    ThisType &operator= ( const ThisType & ) = delete;

    inline const QuadratureType &quadrature () const
    {
      return quadrature_;
    }

    unsigned int index () const { return index_; }
    const CoordinateType &position () const { return quadrature().point( index() ); }
    const RealType &weight () const { return quadrature().weight( index() ); }
    const LocalCoordinateType &localPosition () const { return quadrature().localPoint( index() ); }
  };



  /** \brief   extract the real coordinate from a point
   *  \relates Dune::Fem::QuadraturePointWrapper
   *
   *  This function unwraps a possibly \ref Dune::Fem::QuadraturePointWrapper
   *  "wrapped quadrature point". If the point is not a wrapped quadrature
   *  point, the point itself is returned. This makes it possible to have
   *  one evaluation method for the two different kinds of evaluation points.
   *
   *  \param[in]  x  possibly wrapped point
   *
   *  \returns a reference to the actual point
   */
  template< class Quadrature >
  static inline typename Quadrature::CoordinateType
  coordinate ( const QuadraturePointWrapper< Quadrature > &x )
  {
    return x.position();
  }

  /**
   * \class   QuadraturePointIterator
   * \ingroup Quadrature
   * \brief   iterator over quadrature points
   */
  template< class Quadrature >
  class QuadraturePointIterator
  {
    typedef QuadraturePointIterator< Quadrature > ThisType;

  public:
    typedef std::forward_iterator_tag iterator_category;
    typedef QuadraturePointWrapper< Quadrature > value_type;
    typedef std::ptrdiff_t difference_type;
    typedef Envelope< QuadraturePointWrapper< Quadrature > > pointer;
    typedef QuadraturePointWrapper< Quadrature > reference;

    QuadraturePointIterator () noexcept = default;

    QuadraturePointIterator ( const Quadrature &quadrature, std::size_t point ) noexcept
      : quadrature_( &quadrature ), point_( point  )
    {}

    reference operator* () const noexcept { return value_type( quadrature(), point_ ); }
    pointer operator-> () const noexcept { return pointer( value_type( quadrature(), point_ ) ); }

    bool operator== ( const ThisType &other ) const noexcept { return (point_ == other.point_); }
    bool operator!= ( const ThisType &other ) const noexcept { return (point_ != other.point_); }

    ThisType &operator++ () noexcept { ++point_; return *this; }
    ThisType operator++ ( int ) noexcept { ThisType copy( *this ); ++(*this); return copy; }

    const Quadrature &quadrature () const noexcept { assert( quadrature_ ); return *quadrature_; }

  protected:
    const Quadrature *quadrature_ = nullptr;
    std::size_t point_ = 0;
  };



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

    typedef IntegrationTraits< FieldType, dimension > Traits;

  private:
    typedef IntegrationPointList< FieldType, dimension, IntegrationTraits > ThisType;

    typedef QuadratureProvider< FieldType, dimension, IntegrationTraits >
      QuadratureProviderType;

  public:
    //! type of integration point list implementation
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;

    //! type of coordinate
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! type of key to identify quadrature on user side (default the order of the quadrature)
    typedef typename Traits :: QuadratureKeyType  QuadratureKeyType;

    typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;

    typedef std::shared_ptr< const IntegrationPointListType >  IntegrationPointListStorageType;

    //! to be revised, look at caching quad
    enum { codimension = 0 };

  protected:
    // when the integration point list is obtained from the
    // QuadratureProvider then it should not be deleted
    struct NoDelete
    {
      void operator()( const IntegrationPointListType* ) {}
    };

    IntegrationPointListStorageType ipListPtr_;

  public:
    /** \brief create a quadrature for a given geometry type and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     *
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometryType  geometry type of the requested quadrature
     *  \param[in]  order         order of the requested quadrature
     */
    inline IntegrationPointList ( const GeometryType &geometryType,
                                  const QuadratureKeyType& quadKey )
    : ipListPtr_( &QuadratureProviderType :: getQuadrature( geometryType, quadKey ), NoDelete() )
    {
    }

    /** \brief create a quadrature for a given geometry type and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     *
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometryType     geometry type of the requested quadrature
     *  \param[in]  elementGeometry  geometry type of element that resulting
     *              quadrature is used for (in case of face quadratures)
     *  \param[in]  order            order of the requested quadrature
     */
    inline IntegrationPointList ( const GeometryType &geometryType,
                                  const GeometryType &elementGeometry,
                                  const QuadratureKeyType& quadKey )
    : ipListPtr_( &QuadratureProviderType :: getQuadrature( geometryType, elementGeometry, quadKey ), NoDelete() )
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
    inline IntegrationPointList ( const IntegrationPointListType& ipList )
    : ipListPtr_( &ipList, NoDelete() )
    {
    }

    /** \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for agglomeration quadratures
     *
     *  \param[in]  ipListPtr  implementation of the integration point list
     */
    inline IntegrationPointList ( const IntegrationPointListStorageType& ipListPtr )
    : ipListPtr_( ipListPtr )
    {
    }

    /** \brief copy constructor
     *
     *  \param[in]  org  integration point list to be copied
     */
    inline IntegrationPointList ( const IntegrationPointList &org )
    : ipListPtr_( org.ipListPtr_ )
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
      assert( ipListPtr_ );
      return *ipListPtr_;
    }

    /** \brief obtain the number of integration points
     *
     *  \returns number of integration points within this list
     */
    int nop () const
    {
      return ipList().nop();
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
      return ipList().point( i );
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
      return ipList().id();
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
      return ipList().order();
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
    GeometryType geometryType () const
    {
      return ipList().geometryType();
    }

    /** \brief returns list of element interpolation points for a given face quadrature
      */
    auto interpolationPoints( const int reqDim ) const
    {
      return ipList().interpolationPoints( reqDim );
    }

    /** \brief return true if quadrature is also a set of interpolation points
     * for the given shape functions */
    bool isFaceInterpolationQuadrature( const size_t numShapeFunctions ) const
    {
      return ipList().isFaceInterpolationQuadrature( numShapeFunctions );
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

    static const unsigned int dimension = dim ;

    typedef QuadratureTraits< FieldType, dimension > Traits;

  private:
    typedef Quadrature< FieldType, dimension, QuadratureTraits > ThisType;
    typedef IntegrationPointList< FieldType, dimension, QuadratureTraits > BaseType;

    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits >
      QuadratureProviderType;

  public:
    using BaseType :: ipList;

    typedef typename BaseType::IntegrationPointListStorageType  IntegrationPointListStorageType;

    //! type of the implementation (this must actually be a quadrature implementation)
    typedef typename Traits :: IntegrationPointListType IntegrationPointListType;

    //! type of local coordinate vectors
    typedef typename IntegrationPointListType :: CoordinateType CoordinateType;

    //! type of key to identify the quadrature
    typedef typename Traits :: QuadratureKeyType  QuadratureKeyType;

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
     *  \param[in]  geometryType  geometry type of the requested quadrature
     *  \param[in]  key           key to identify the quadrature (default = order)
     */
    inline Quadrature( const GeometryType &geometryType, const QuadratureKeyType &key )
    : BaseType( geometryType, key )
    {
    }

    /** \brief create a quadrature for a given geometry and order
     *
     *  This constructor creates a quadrature for the specified geometry which
     *  is capable of integrating polynoms up the given order exactly.
     *
     *  \note The order of the quadrature may be higher than the requested one.
     *
     *  \param[in]  geometryType     geometry type of the requested quadrature
     *  \param[in]  elementGeometry  geometry type of element that resulting
     *              quadrature is used for (in case of face quadratures)
     *  \param[in]  key              key to identify the quadrature (default = order)
     *
     *  \note This is a specialized constructor for constructing
     *  face quadratures for UGGrid.
     */
    inline Quadrature ( const GeometryType &geometryType,
                        const GeometryType &elementGeometry,
                        const QuadratureKeyType &key )
    : BaseType( geometryType, elementGeometry, key )
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

    /** \brief create an integration point list from an implementation
     *
     *  This constructor creates an integration point list from a given
     *  implementation.
     *
     *  \note This constructor is provided mainly for agglomeration quadratures.
     *
     *  \param[in]  ipListPtr  shared_ptr of implementation of the integration point list
     */
    inline explicit Quadrature( const IntegrationPointListStorageType& ipListPtr )
    : BaseType( ipListPtr )
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


  /** \class SelectQuadraturePointSetId
   *
   * Select point set id from quadrature if available, otherwise set default
   * value (-Dune::QuadratureType::size)
   *
   */
  template <class Quadrature>
  using SelectQuadraturePointSetId = detail::SelectPointSetId< Quadrature, -Dune::QuadratureType::size >;

  } //end namespace Fem

} //end namespace Dune

#endif

#ifndef DUNE_QUADPROVIDER_HH
#define DUNE_QUADPROVIDER_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/quadrature/quadratureimp.hh>
#include <dune/fem/quadrature/idprovider.hh>

namespace Dune
{

  /*! \class QuadCreator
   *  \ingroup Quadrature
   *  \brief the actual quadrature storage
   *
   *  QuadCreator is a utility class providing the actual quadrature storage.
   *
   *  The template argument is used to distinguish classes for different geometry
   *  types (maybe GeometryType :: BasicType would be a better choice).
   */
  template< unsigned int dummy >
  class QuadCreator
  {
  private:
    //! class holding vector with pointer to quadrature objects 
    template< class QuadImp >
    class QuadratureStorage  
    {
    public:
      typedef QuadImp QuadType;

      typedef QuadType *QuadPtr;

    protected:
      DynamicArray< QuadPtr > storage_;
      //std :: vector< QuadPtr > storage_;
      
    public:   
      inline QuadratureStorage ()
      : storage_( QuadType :: maxOrder() + 1, (QuadPtr)0 )
      {
      }
      
      inline ~QuadratureStorage () 
      {
        for( unsigned int i = 0; i < storage_.size(); ++i )
          delete storage_[ i ];
      }

      inline QuadImp &getQuadrature( const GeometryType &geometry,
                                     unsigned int order )
      {
        assert( order < storage_.size() );
        
        QuadPtr& quadPtr = storage_[ order ];
        if( quadPtr == 0 )
          quadPtr = new QuadImp( geometry, order, IdProvider :: instance().newId() );
        return *quadPtr;
      }      
    }; // end class QuadratureStorage 
  
  public:
    /*! \brief provide quadrature
     *
     *  \param[in]  geometry  type of geometry, the quadrature is requested for
     *  \param[in]  order     minimal order of the requested quadrature
     */
    template< class QuadImp >
    static const QuadImp &provideQuad( const GeometryType &geometry,
                                       unsigned int order )
    {
      static QuadratureStorage< QuadImp > storage;
      return storage.getQuadrature( geometry, order );
    }
  };



  /*! \class QuadratureProvider
   *  \ingroup Quadrature
   *  \brief provide a single instance pool of quadratures
   *
   *  QuadratureProvider follows the monostate pattern. It provides a single
   *  point of access (and storage) for the actual implementation of 
   *  quadratures. Hence, the expensive creations of quadratures should be
   *  reduced to a minimum.
   *
   *  There are the following specializations:
   *  - QuadratureProvider<FieldImp,0,QuadratureTraits>
   *  - QuadratureProvider<FieldImp,1,QuadratureTraits>
   *  - QuadratureProvider<FieldImp,2,QuadratureTraits>
   *  - QuadratureProvider<FieldImp,3,QuadratureTraits>
   */
  template< typename FieldImp, int dim, template< class, int > class QuadratureTraits >
  class QuadratureProvider
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = dim };

  private:
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits > ThisType;
    
    typedef QuadratureTraits< FieldType, dimension > QuadratureTraitsType;

  public:
    //! type for cube quadrature
    typedef typename QuadratureTraitsType :: CubeQuadratureType CubeQuadratureType;

    //! type of integration point list implementation
    typedef typename QuadratureTraitsType :: IntegrationPointListType
      IntegrationPointListType;

    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          int order )
    {
      assert( geometry.isCube() );
      return QuadCreator< 0 > :: template provideQuad< CubeQuadratureType > ( geometry, order );
    }
    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          const GeometryType &elementGeometry,
                                                          int order )
    {
      return getQuadrature( geometry, order );
    }
    
  private:
    // forbid creation
    QuadratureProvider();
    
    // forbid copying
    QuadratureProvider( const ThisType& );
   
    // forbid assignment
    QuadratureProvider &operator=( const ThisType& );
  };



  /** \copydoc Dune::QuadratureProvider */
  template< typename FieldImp, template< class, int > class QuadratureTraits >
  class QuadratureProvider< FieldImp, 0, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = 0 };

  private:
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits > ThisType;

    typedef QuadratureTraits< FieldType, dimension > QuadratureTraitsType;

  public:
    //! type of point quadrature
    typedef typename QuadratureTraitsType :: PointQuadratureType PointQuadratureType;

    //! type of integration point list implementation
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;

  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex() );
      assert( order >= 0 );
      //return QuadCreator< 0 > :: template provideQuad< PointQuadratureType > ( geometry, order );
      static PointQuadratureType quad( geometry, 
                                       order, 
                                       IdProvider ::instance().newId()); 
      return quad; 
    }

    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          const GeometryType &elementGeometry,
                                                          int order )
    {
      return getQuadrature(geometry, order);
    }

  private:
    // forbid creation
    QuadratureProvider();
    
    // forbid copying
    QuadratureProvider( const ThisType& );
   
    // forbid assignment
    QuadratureProvider &operator=( const ThisType& );
  }; 
  

  
  /** \copydoc Dune::QuadratureProvider */
  template< class FieldImp, template< class, int > class QuadratureTraits >
  class QuadratureProvider< FieldImp, 1, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = 1 };

  private:
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits > ThisType;

    typedef QuadratureTraits< FieldType, dimension > QuadratureTraitsType;

  public:
    //! type of line quadrature
    typedef typename QuadratureTraitsType :: LineQuadratureType LineQuadratureType;

    //! type of integration point list implementation
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;

  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex() );
      assert( order >= 0 );
      return QuadCreator< 0 > :: template provideQuad< LineQuadratureType > ( geometry, order );
    }

    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          const GeometryType &elementGeometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex() );
      assert( order >= 0 );
      // we need here to distinguish between the basic types 
      // otherwise the this won't work for UGGrid 
      return ( elementGeometry.basicType() == GeometryType :: simplex ) ? 
        QuadCreator< 0 > :: template provideQuad< LineQuadratureType > ( geometry, order ) :
        QuadCreator< 1 > :: template provideQuad< LineQuadratureType > ( geometry, order ) ;
    }

  private:
    // forbid creation
    QuadratureProvider();
    
    // forbid copying
    QuadratureProvider( const ThisType& );
   
    // forbid assignment
    QuadratureProvider &operator=( const ThisType& );
  }; 



  /** \copydoc Dune::QuadratureProvider */
  template< class FieldImp, template< class, int > class QuadratureTraits >
  class QuadratureProvider< FieldImp, 2, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = 2 };

  private:
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits > ThisType;

    typedef QuadratureTraits< FieldType, dimension > QuadratureTraitsType;

  public:
    //! type of simplex quadrature
    typedef typename QuadratureTraitsType :: SimplexQuadratureType SimplexQuadratureType;
    //! type of cube quadrature
    typedef typename QuadratureTraitsType :: CubeQuadratureType CubeQuadratureType;

    //! type of integration point list implementation
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;

  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex() );
      assert( order >= 0 );

      if( geometry.isSimplex() ) 
      {
        return QuadCreator< 0 > :: 
          template provideQuad< SimplexQuadratureType > ( geometry, order ); 
      }
      else 
      {
        return QuadCreator< 1 > :: 
          template provideQuad< CubeQuadratureType >    ( geometry, order ) ;
      }
    }

    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          const GeometryType &elementGeometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex() );
      assert( order >= 0 );

      // if geometry is simplex return simplex quadrature 
      if ( geometry.isSimplex() ) 
      {
        // check element geometry to provide quadratures with different ids 
        if( elementGeometry.isSimplex() )
          return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isCube() )
          return QuadCreator< 1 > :: template provideQuad< SimplexQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isPrism() )
          return QuadCreator< 2 > :: template provideQuad< SimplexQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isPyramid() )
          return QuadCreator< 3 > :: template provideQuad< SimplexQuadratureType > ( geometry, order ) ;
        else 
          DUNE_THROW( RangeError, "Element type not available for dimension 3" );
      }
      else 
      {
        // return cube quadrature 
        // check element geometry to provide quadratures with different ids 
        if( elementGeometry.isSimplex() )
          return QuadCreator< 4 > :: template provideQuad< CubeQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isCube() )
          return QuadCreator< 5 > :: template provideQuad< CubeQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isPrism() )
          return QuadCreator< 6 > :: template provideQuad< CubeQuadratureType > ( geometry, order ) ;
        else if( elementGeometry.isPyramid() )
          return QuadCreator< 7 > :: template provideQuad< CubeQuadratureType > ( geometry, order ) ;
        else 
          DUNE_THROW( RangeError, "Element type not available for dimension 3" );
      }

      DUNE_THROW( RangeError, "Element type not available for dimension 2" );
      // dummy return
      return QuadCreator< 0 > :: 
        template provideQuad< SimplexQuadratureType >( geometry, 0 );
    }

  private:
    // forbid creation
    QuadratureProvider();
    
    // forbid copying
    QuadratureProvider( const ThisType& );
   
    // forbid assignment
    QuadratureProvider &operator=( const ThisType& );
  };


  
  /** \copydoc Dune::QuadratureProvider */
  template< class FieldImp, template< class, int > class QuadratureTraits >
  class QuadratureProvider< FieldImp, 3, QuadratureTraits >
  {
  public:
    typedef FieldImp FieldType;

    enum { dimension = 3 };

  private:
    typedef QuadratureProvider< FieldType, dimension, QuadratureTraits > ThisType;

    typedef QuadratureTraits< FieldType, dimension > QuadratureTraitsType;

  public:
    //! type of simplex quadrature
    typedef typename QuadratureTraitsType :: SimplexQuadratureType SimplexQuadratureType;
    //! type of cube quadrature
    typedef typename QuadratureTraitsType :: CubeQuadratureType CubeQuadratureType;
    //! type of prims quadrature
    typedef typename QuadratureTraitsType :: PrismQuadratureType PrismQuadratureType;
    //! type of pyramid quadrature
    typedef typename QuadratureTraitsType :: PyramidQuadratureType PyramidQuadratureType;

    //! type of integration point list implementation
    typedef typename QuadratureTraitsType :: IntegrationPointListType IntegrationPointListType;

  public:
    //! Access to the quadrature implementations.
    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          int order )
    {
      assert( geometry.isCube() || geometry.isSimplex()
              || geometry.isPrism() || geometry.isPyramid() );
      assert( order >= 0 );
      
      if( geometry.isSimplex() )
        return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType >
          ( geometry, order );
      if( geometry.isCube() )
        return QuadCreator< 1 > :: template provideQuad< CubeQuadratureType >
          ( geometry, order );
      
      if( geometry.isPrism() )
        return QuadCreator< 2 > :: template provideQuad< PrismQuadratureType >
          ( geometry, order );
      if( geometry.isPyramid() )
        return QuadCreator< 3 > :: template provideQuad< PyramidQuadratureType >
          ( geometry, order );
      
      DUNE_THROW( RangeError, "Element type not available for dimension 3" );
      // dummy return
      return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType >
        ( geometry, 0 );
    }

    static const IntegrationPointListType &getQuadrature( const GeometryType &geometry,
                                                          const GeometryType &elementGeometry,
                                                          int order )
    {
      DUNE_THROW( RangeError, "QuadProvider::getQuadrature not implemented for 3d face quadratures!" );
      // dummy return
      return QuadCreator< 0 > :: template provideQuad< SimplexQuadratureType >
        ( geometry, 0 );
    }
  private:
    // forbid creation
    QuadratureProvider();
    
    // forbid copying
    QuadratureProvider( const ThisType& );
   
    // forbid assignment
    QuadratureProvider &operator=( const ThisType& );
  };

} // end namespace Dune 
#endif

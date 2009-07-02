#ifndef DUNE_REFERENCEELEMENTPROVIDER_HH
#define DUNE_REFERENCEELEMENTPROVIDER_HH

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{

  // ReferenceElementProviderBase
  // ----------------------------

  template< class Grid, int codim >
  struct ReferenceElementProviderBase
  {
    typedef Grid GridType;

    typedef typename GridType::ctype FieldType;
    static const int dimension = GridType::dimension;
    static const int codimension = codim;

    typedef GenericReferenceElement< FieldType, dimension-codimension >
      ReferenceElementType;

  protected:
    typedef GenericReferenceElementContainer< FieldType, dimension-codimension >
      ContainerType;

    ReferenceElementProviderBase ()
    : container_( ContainerType::instance() )
    {}

    const ContainerType &container () const
    {
      return container_;
    }

  private:
    const ContainerType &container_;
  };



  // ReferenceElementProvider
  // ------------------------

  template< class Grid, int codim = 0 >
  class ReferenceElementProvider
  : public ReferenceElementProviderBase< Grid, codim >
  {
    typedef ReferenceElementProviderBase< Grid, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container() ( object.type() );
    }
  };


  template< class Grid, int codim >
  class ReferenceElementProvider< const Grid, codim >
  : public ReferenceElementProvider< Grid, codim >
  {};



  // ReferenceElementProvider for AlbertaGrid
  // ----------------------------------------

#if HAVE_ALBERTA
  template< int dim, int dimworld >
  class AlbertaGrid;

  template< int dim, int dimworld, int codim >
  class ReferenceElementProvider< AlbertaGrid< dim, dimworld >, codim >
  : public ReferenceElementProviderBase< AlbertaGrid< dim, dimworld >, codim >
  {
    typedef ReferenceElementProviderBase< AlbertaGrid< dim, dimworld >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().simplex();
    }
  };
#endif // #if HAVE_ALBERTA



  // ReferenceElementProvider for AluConformGrid
  // -------------------------------------------

#if HAVE_ALUGRID
  template< int dim, int dimworld >
  class AluConformGrid;

  template< int dim, int dimworld, int codim >
  class ReferenceElementProvider< AluConformGrid< dim, dimworld >, codim >
  : public ReferenceElementProviderBase< AluConformGrid< dim, dimworld >, codim >
  {
    typedef ReferenceElementProviderBase< AluConformGrid< dim, dimworld >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().simplex();
    }
  };
#endif // #if HAVE_ALUGRID



  // ReferenceElementProvider for AluSimplexGrid
  // -------------------------------------------

#if HAVE_ALUGRID
  template< int dim, int dimworld >
  class AluSimplexGrid;

  template< int dim, int dimworld, int codim >
  class ReferenceElementProvider< AluSimplexGrid< dim, dimworld >, codim >
  : public ReferenceElementProviderBase< AluSimplexGrid< dim, dimworld >, codim >
  {
    typedef ReferenceElementProviderBase< AluSimplexGrid< dim, dimworld >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().simplex();
    }
  };
#endif // #if HAVE_ALUGRID



  // ReferenceElementProvider for AluCubeGrid
  // ----------------------------------------

#if HAVE_ALUGRID
  template< int dim, int dimworld >
  class AluCubeGrid;

  template< int dim, int dimworld, int codim >
  class ReferenceElementProvider< AluCubeGrid< dim, dimworld >, codim >
  : public ReferenceElementProviderBase< AluCubeGrid< dim, dimworld >, codim >
  {
    typedef ReferenceElementProviderBase< AluCubeGrid< dim, dimworld >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().cube();
    }
  };
#endif // #if HAVE_ALUGRID



  // ReferenceElementProvider for OneDGrid
  // -------------------------------------

  class OneDGrid;

  template< int codim >
  class ReferenceElementProvider< OneDGrid, codim >
  : public ReferenceElementProviderBase< OneDGrid, codim >
  {
    typedef ReferenceElementProviderBase< OneDGrid, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().cube();
    }
  };



  // ReferenceElementProvider for SGrid
  // ----------------------------------

  template< int dim, int dimworld >
  class SGrid;

  template< int dim, int dimworld, int codim >
  class ReferenceElementProvider< SGrid< dim, dimworld >, codim >
  : public ReferenceElementProviderBase< SGrid< dim, dimworld >, codim >
  {
    typedef ReferenceElementProviderBase< SGrid< dim, dimworld >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().cube();
    }
  };


  // ReferenceElementProvider for YaspGrid
  // -------------------------------------

  template< int dim >
  class YaspGrid;

  template< int dim, int codim >
  class ReferenceElementProvider< YaspGrid< dim >, codim >
  : public ReferenceElementProviderBase< YaspGrid< dim >, codim >
  {
    typedef ReferenceElementProviderBase< YaspGrid< dim >, codim > BaseType;

  public:
    typedef typename BaseType::ReferenceElementType ReferenceElementType;

    template< class Object >
    const ReferenceElementType &operator() ( const Object &object ) const
    {
      return BaseType::container().cube();
    }
  };

}

#endif // #ifndef DUNE_REFERENCEELEMENTPROVIDER_HH

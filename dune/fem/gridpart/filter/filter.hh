#ifndef DUNE_FEM_GRIDPART_FILTER_FILTER_HH
#define DUNE_FEM_GRIDPART_FILTER_FILTER_HH

//- system includes
#include <algorithm>
#include <cassert>
#include <vector>

//- dune-common includes
#include <dune/common/bartonnackmanifcheck.hh>


namespace Dune
{

  namespace Fem
  {

    // forward declarations
    // --------------------

    template< class > struct DefaultFilterTraits;
    template< class > class FilterInterface;
    template< class > class FilterDefaultImplementation;


    // DefaultFilterTraits
    // -------------------

    //! \brief type definitions
    template < class FilterImp >
    struct DefaultFilterTraits
    {
      //! \brief filter type
      typedef FilterImp FilterType;

      //! \brief entity types
      template < int cd >
      struct Codim
      {
        //! \brief entity type for given codimension
        typedef typename FilterType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief entity type for codimension 0
      typedef typename Codim< 0 >::EntityType EntityType;

    }; // end DefaultFilterTraits



    // FilterInterface
    // ---------------

    /** \ingroup FilterGridPart
     *  \brief   Interface class for filter to use with a Dune::FilteredGridPart
     */
    template< class FilterTraits >
    class FilterInterface
    {
      typedef FilterInterface< FilterTraits > ThisType;

      friend class FilterDefaultImplementation< FilterTraits >;

      typedef FilterTraits Traits;

    public:
      //! \brief type of the filter implementation
      typedef typename Traits :: FilterType FilterType;

      //! \brief entity types
      template< int cd >
      struct Codim
      {
        typedef typename Traits::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of entity with codim=0
      typedef typename Codim< 0 >::EntityType EntityType;

    private:
      FilterInterface () = default;

      FilterInterface ( const ThisType & ) = default;
      FilterInterface ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;
      ThisType &operator= ( ThisType && ) = default;

    public:
      //! \brief returns true if the given entity of the pointer in the domain
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().contains( entity ) );
        return asImp().contains< cd >( entity );
      }

      //! \brief returns true if the given entity of the pointer in the domain
      template< class Entity >
      bool contains ( const Entity & entity ) const
      {
        enum { cc = Entity::codimension };
        CHECK_INTERFACE_IMPLEMENTATION( asImp().contains< cc >( entity ) );
        return asImp().contains< cc >( entity );
      }

      //! \brief returns true if an intersection is interior
      //         (allows boundarys within a given domain)
      template< class Intersection >
      bool interiorIntersection ( const Intersection &intersection ) const
      {
        return asImp().interiorIntersection( intersection );
      }

      //! \brief returns true if an intersection is a boundary intersection
      template< class Intersection >
      bool intersectionBoundary( const Intersection &intersection ) const
      {
        return asImp().intersectionBoundary( intersection );
      }

      //! \brief returns the boundary id for an intersection
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection &intersection ) const
      {
        return asImp().intersectionBoundaryId( intersection );
      }

      //! \brief returns true if for an intersection a neighbor exsits
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection &intersection ) const
      {
        return asImp().intersectionNeighbor( intersection );
      }

    protected:
      FilterType &asImp ()
      {
        return static_cast< FilterType & >( *this );
      }

      const FilterType &asImp () const
      {
        return static_cast< const FilterType & >( *this );
      }
    };


    // FilterDefaultImplementation
    // ---------------------------

    template< class FilterTraits >
    class FilterDefaultImplementation
    : public FilterInterface< FilterTraits >
    {
      typedef FilterDefaultImplementation< FilterTraits > ThisType;
      typedef FilterInterface< FilterTraits > BaseType;

    public:
      //! \brief type of the filter implementation
      typedef typename BaseType::FilterType FilterType;

      //! \brief entity types
      template< int cd >
      struct Codim
      {
        //! \brief type of codim cd
        typedef typename BaseType::template Codim< cd >::EntityType EntityType;
      };

      //! \brief type of codim 0 entity
      typedef typename BaseType::EntityType EntityType;

    protected:
      using BaseType::asImp;

      FilterDefaultImplementation () = default;

      FilterDefaultImplementation ( const ThisType & ) = default;
      FilterDefaultImplementation ( ThisType && ) = default;

      ThisType &operator= ( const ThisType & ) = default;
      ThisType &operator= ( ThisType && ) = default;

    public:
      using BaseType::contains;

      //! \brief default implementation returns contains from neighbor
      template< class Intersection >
      bool interiorIntersection( const Intersection &intersection ) const
      {
        typedef typename Intersection::Entity EntityType;
        const EntityType outside(intersection.outside());
        return asImp().contains( outside );
      }

      //! \brief returns true if the given entity of the pointer in the domain
      template< int cd >
      bool contains ( const typename Codim< cd >::EntityType & ) const;

      //! \brief returns true if an intersection is a boundary intersection
      template< class Intersection >
      bool intersectionBoundary( const Intersection & ) const;

      //! \brief returns the boundary id for an intersection
      template< class Intersection >
      int intersectionBoundaryId ( const Intersection & ) const;

      //! \brief returns true if for an intersection a neighbor exsits
      template< class Intersection >
      bool intersectionNeighbor ( const Intersection & ) const;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_FILTER_HH

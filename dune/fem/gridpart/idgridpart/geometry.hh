#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_GEOMETRY_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_GEOMETRY_HH

#include <type_traits>

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< int, int, class > class IdGeometry;
    template< int, int, class > class IdLocalGeometry;


    // IdBasicGeometry
    // ---------------

    template< class Traits >
    struct IdBasicGeometry
    {
      typedef typename Traits::HostGeometryType HostGeometryType;

      static const int mydimension = HostGeometryType::mydimension;
      static const int coorddimension = HostGeometryType::coorddimension;

      typedef typename HostGeometryType::ctype ctype;
      typedef FieldVector< ctype, mydimension > LocalVector;
      typedef FieldVector< ctype, coorddimension > GlobalVector;

      typedef typename HostGeometryType::JacobianTransposed JacobianTransposed;
      typedef typename HostGeometryType::JacobianInverseTransposed JacobianInverseTransposed;

      IdBasicGeometry ( const HostGeometryType &hostGeometry )
      : hostGeometry_( hostGeometry )
      {}

      operator bool () const { return bool( hostGeometry_ ); }

      GeometryType type () const { return hostGeometry_.type(); }
      bool affine () const { return hostGeometry_.affine(); }

      int corners () const { return hostGeometry_.corners(); }
      GlobalVector corner ( const int i ) const { return hostGeometry_.corner( i ); }
      GlobalVector center () const { return hostGeometry_.center(); }

      GlobalVector global ( const LocalVector &local ) const { return hostGeometry_.global( local ); }
      LocalVector local ( const GlobalVector &global ) const { return hostGeometry_.local( global ); }

      ctype integrationElement ( const LocalVector &local ) const { return hostGeometry_.integrationElement( local ); }
      ctype volume () const { return hostGeometry_.volume(); }

      JacobianTransposed jacobianTransposed ( const LocalVector &local ) const
      {
        return hostGeometry_.jacobianTransposed( local );
      }

      JacobianInverseTransposed jacobianInverseTransposed ( const LocalVector &local ) const
      {
        return hostGeometry_.jacobianInverseTransposed( local );
      }

    private:
      HostGeometryType hostGeometry_;
    };



    // IdGeometryTraits
    // ----------------

    template< int mydim, class GridFamily >
    struct IdGeometryTraits
    {
      typedef typename std::remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      static const int dimension = HostGridPartType::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;

      typedef typename HostGridPartType::template Codim< codimension >::GeometryType HostGeometryType;
    };



    // IdGeometry
    // ----------

    template< int mydim, int cdim, class GridFamily >
    class IdGeometry
    : public IdBasicGeometry< IdGeometryTraits< mydim, GridFamily > >
    {
      typedef IdBasicGeometry< IdGeometryTraits< mydim, GridFamily > > Base;

    public:
      typedef typename Base::HostGeometryType HostGeometryType;

      IdGeometry ()
      {}

      IdGeometry ( const HostGeometryType &hostGeometry )
      : Base( hostGeometry )
      {}
    };



    // IdLocalGeometryTraits
    // ---------------------

    template< int mydim, class GridFamily >
    struct IdLocalGeometryTraits
    {
      typedef typename std::remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      static const int dimension = HostGridPartType::dimension;
      static const int mydimension = mydim;
      static const int codimension = dimension - mydimension;

      typedef typename HostGridPartType::template Codim< codimension >::LocalGeometryType HostGeometryType;
    };



    // IdLocalGeometry
    // -------------------

    template< int mydim, int cdim, class GridFamily >
    class IdLocalGeometry
    : public IdBasicGeometry< IdLocalGeometryTraits< mydim, GridFamily > >
    {
      typedef IdBasicGeometry< IdLocalGeometryTraits< mydim, GridFamily > > Base;

    public:
      typedef typename Base::HostGeometryType HostGeometryType;

      IdLocalGeometry ()
      {}

      IdLocalGeometry ( const HostGeometryType &hostGeometry )
      : Base( hostGeometry )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_GEOMETRY_HH

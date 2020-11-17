#ifndef DUNE_FEM_TWISTPROVIDER_HH
#define DUNE_FEM_TWISTPROVIDER_HH

//- System includes
#include <cassert>

#include <map>
#include <memory>
#include <vector>

//- Dune includes
#include <dune/geometry/referenceelements.hh>

#include <dune/fem/quadrature/quadrature.hh>

#include <dune/fem/storage/singleton.hh>

//- Local includes
#include "pointmapper.hh"
#include "topology.hh"

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    template <class ct, int dim>
    class TwistProvider;
    template <class ct, int dim>
    class TwistMapperCreator;
    template <class ct, int dim>
    class TwistStorage;


    //! \brief Identifies quadrature points on faces with twists
    //! For a given quadrature type and a face with a given twist the TwistMapper
    //! provides a mapping from the quadrature point number on the twisted face
    //! to the quadrature point number on the (untwisted) reference face. (It
    //! removes the twist from the quadrature, so to say.) This is needed in
    //! unstructured grids when a quadrature point on a reference element's
    //! face needs to be transformed to a quadrature point in the reference
    //! element itself.
    //!
    //! The PointMapper objects are filled by the TwistMapperCreator and its
    //! subclasses.


    //! \brief Helper class which stores information about twist mapping for a
    //! given quadrature id.
    template <class ct, int dim>
    class TwistStorage
    {
      typedef CachingTraits<ct, dim> Traits;
    public:
      typedef typename Traits::PointType PointType;
      typedef typename Traits::PointVectorType PointVectorType;
      typedef typename Traits::MapperType MapperType;

    public:
      //! Constructor
      //! \param minTwist Minimal possible twist for the face
      //! \param maxTwist Maximal possible twist for the face + 1
      explicit TwistStorage(int minTwist, int maxTwist);

      //! Add a new mapper for a given twist
      void addMapper(const MapperType& mapper, int twist);

      //! Add a point (in the case of asymmetric quadratures)
      size_t addPoint(const PointType& points);

      //! Access to a mapper
      const MapperType& getMapper(int twist) const;

      //! Access to the points for all twists (in the case of symmetric
      //! quadratures, the points are identical with the quadrature points).
      const PointVectorType& getPoints() const;

      //! Minimal twist
      int minTwist() const { return minTwist_; }

      //! Maximal twist + 1
      int maxTwist() const { return maxTwist_; }

    private:
      typedef typename Traits::MapperVectorType MapperVectorType;

    private:
      MapperVectorType mappers_;
      PointVectorType points_;

      int minTwist_;
      int maxTwist_;
    };



    // TiwstProvider
    // -------------

    //! \brief Access point for PointMapper objects with twist information
    //! PointMapper objects get created once and are reused as often as needed.
    //! The TwistProvider serves in this context as the single point of access
    //! which is responsible for the creation and management of these objects.
    //! TwistProvider follows the monostate pattern.
    template <class ct, int dim>
    class TwistProvider
    {
      typedef CachingTraits<ct, dim> Traits;
    public:
      //! Generic quadrature type
      typedef typename Traits::QuadratureType QuadratureType;
      //! Storage for the mappings of a specific quadrature id, alongside the
      //! resulting caching points (which may differ from the quadrature points
      //! in the case of asymmetric quadratures)
      typedef typename TwistMapperCreator<ct, dim>::TwistStorageType TwistStorageType;

    public:
      //! Delivers the PointMapper object for quadrature quad and twist twist.
      static const TwistStorageType& getTwistStorage(const QuadratureType& quad);

    private:
      typedef std::vector< const TwistStorageType* > MapperContainerType;
      typedef typename MapperContainerType::iterator IteratorType;

    private:
      // singleton class holding map with storages
      class MapperContainer
      {
        // instance of map
        MapperContainerType mappers_;

      public:
        //! constructor
        MapperContainer() : mappers_(100, (TwistStorageType*) 0)
        {}

        //! destructor
        ~MapperContainer()
        {
          IteratorType endit = mappers_.end();
          for(IteratorType it = mappers_.begin(); it != endit; ++it)
          {
            delete (*it);
          }
        }

        friend class Dune::Fem::Singleton< MapperContainerType >;
        //! return reference to mappers
        static MapperContainerType& instance()
        {
          return Singleton< MapperContainer > :: instance().mappers_;
        }
      };
    };



    // TwistMapperStrategy
    // -------------------

    //! This class factors out all geometry dependent stuff in a strategy class
    template <class ct, int dim>
    struct TwistMapperStrategy
    {
      typedef FieldMatrix<ct, dim+1, dim> MatrixType;

    public:
      TwistMapperStrategy(int minTwist, int maxTwist) :
        minTwist_(minTwist),
        maxTwist_(maxTwist)
      {}

      //! virtual desctructor because of virtual functions
      virtual ~TwistMapperStrategy () = default;

      virtual MatrixType buildTransformationMatrix ( int twist ) const = 0;

      int minTwist() const { return minTwist_; }
      int maxTwist() const { return maxTwist_; }

    private:
      int minTwist_;
      int maxTwist_;
    };



    // TwistMapperCreator
    // ------------------

    //! Helper class for TwistProvider which takes care of the creation process
    template <class ct, int dim>
    class TwistMapperCreator
    {
      typedef CachingTraits<ct, dim> Traits;
    public:
      typedef typename Traits::QuadratureType QuadratureType;
      typedef typename Traits::PointType PointType;
      typedef typename Traits::MapperType MapperType;
      typedef FieldVector<ct, dim+1> CoordinateType;
      typedef TwistStorage<ct, dim> TwistStorageType;

    public:
      //! Constructor
      TwistMapperCreator(const QuadratureType& quad);

      //! Destructor
      ~TwistMapperCreator();

      //! Create the actual mapper for a given twist
      const TwistStorageType* createStorage() const;

      //! Lowest possible twist for the quadrature's geometry
      int minTwist() const {
        return helper_->minTwist();
      }

      //! Largest possible twist + 1 for the quadrature's geometry
      int maxTwist() const {
        return helper_->maxTwist();
      }

    private:
      typedef typename TwistMapperStrategy<ct, dim>::MatrixType MatrixType;

    private:
      TwistMapperCreator(const TwistMapperCreator&);
      TwistMapperCreator& operator=(const TwistMapperCreator&);

    private:
      const QuadratureType& quad_;
      TwistMapperStrategy<ct, dim>* helper_;

      static const ct eps_;
    };



    // PointTwistMapperStrategy
    // ------------------------

    //! Implements the creator's functionality that depends on the underlying
    //! geometry. This is the special implementation for points.
    template <class ct, int dim>
    class PointTwistMapperStrategy final
      : public TwistMapperStrategy<ct, dim>
    {
      typedef TwistMapperStrategy<ct, dim> BaseType;
      typedef typename BaseType::MatrixType MatrixType;

    public:
      PointTwistMapperStrategy(GeometryType geo)
        : BaseType( 0, 1 )
      {
        assert(dim == 0);
      }

      virtual MatrixType buildTransformationMatrix ( int twist ) const
      {
        const auto &refElement = Dune::ReferenceElements< ct, dim >::cube();

        assert( twist == 0 );
        MatrixType mat;
        mat[ 0 ] = refElement.position( 0, dim );
        return mat;
      }
    };



    // LineTwistMapperStrategy
    // -----------------------

    //! Implements the creator's functionality that depends on the underlying
    //! geometry. This is the special implementation for line.
    template <class ct, int dim>
    class LineTwistMapperStrategy final
      : public TwistMapperStrategy<ct, dim>
    {
      typedef TwistMapperStrategy<ct, dim> BaseType;
      typedef typename BaseType::MatrixType MatrixType;

    public:
      LineTwistMapperStrategy(GeometryType geo)
        : BaseType( 0, 2 )
      {
        assert(dim == 1);
      }

      MatrixType buildTransformationMatrix ( int twist ) const override
      {
        const auto &refElement = Dune::ReferenceElements< ct, dim >::cube();

        assert( (twist == 0) || (twist == 1) );
        MatrixType mat;
        mat[ twist   ] = refElement.position( 0, dim );
        mat[ 1-twist ] = refElement.position( 1, dim );
        return mat;
      }
    };



    // TriangleTwistMapperStrategy
    // ---------------------------

    //! Implements the creator's functionality that depends on the underlying
    //! geometry. This is the special implementation for triangle.
    template <class ct, int dim>
    class TriangleTwistMapperStrategy final
      : public TwistMapperStrategy<ct, dim>
    {
      typedef TwistMapperStrategy<ct, dim> BaseType;
      typedef typename BaseType::MatrixType MatrixType;

    public:
      TriangleTwistMapperStrategy(GeometryType geo)
        : BaseType( -3, 3 )
      {
        assert(dim == 2);
      }

      virtual MatrixType buildTransformationMatrix ( int twist ) const override
      {
        typedef Dune::Fem::FaceTopologyMapping<tetra> FaceTopo;

        const auto &refElement = Dune::ReferenceElements< ct, dim >::simplex();

        MatrixType mat( ct( 0 ) );
        for (int idx = 0; idx < dim+1; ++idx)
          mat[idx] = refElement.position( FaceTopo::twistedDuneIndex(idx, twist), dim); // dim == codim here
        return mat;
      }
    };



    // QuadrilateralTwistMapperStrategy
    // --------------------------------

    //! Implements the creator's functionality that depends on the underlying
    //! geometry. This is the special implementation for quadrilaterals
    template <class ct, int dim>
    class QuadrilateralTwistMapperStrategy final
      : public TwistMapperStrategy<ct, dim>
    {
      typedef TwistMapperStrategy<ct, dim> BaseType;
      typedef typename BaseType::MatrixType MatrixType;

    public:
      QuadrilateralTwistMapperStrategy(GeometryType geo)
        : BaseType( -4, 4 )
      {
        assert(dim == 2);
      }

      MatrixType buildTransformationMatrix ( int twist ) const override
      {
        typedef Dune::Fem::FaceTopologyMapping<hexa> FaceTopo;

        const auto &refElement = Dune::ReferenceElements< ct, dim >::cube();

        MatrixType mat( ct( 0 ) );
        for (int idx = 0; idx < dim+1; ++idx)
          mat[idx] = refElement.position( FaceTopo::twistedDuneIndex(idx, twist), dim); // dim == codim here
        return mat;
      }
    };

  } // namespace Fem

} // namespace Dune

#include "twistprovider.cc"
#endif // #ifndef DUNE_FEM_TWISTPROVIDER_HH

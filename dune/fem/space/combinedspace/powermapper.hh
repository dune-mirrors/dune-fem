#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH
#define DUNE_FEM_COMBINEDSPACE_MAPPER_HH

#include <dune/common/math.hh>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  namespace Fem
  {

    // PowerMapper
    // -----------

    template< class Mapper, int N >
    class PowerMapper;


    // PowerMapperTraits
    // -----------------

    template< class Mapper, int N >
    struct PowerMapperTraits
    {
      typedef Mapper MapperType;

      // type of codim 0 elements
      typedef typename MapperType::ElementType ElementType;

      // type of Size's
      typedef typename MapperType::SizeType SizeType;

      // type of dofmapper
      typedef PowerMapper< MapperType, N >  DofMapperType;
    };



    /** \class PowerMapper
     *  \ingroup PowerSpace
     *  \brief DofMapper for the PowerSpace
     */
    template< class Mapper, int N >
    class PowerMapper
      : public AdaptiveDofMapper< PowerMapperTraits< Mapper, N > >
    {
      typedef PowerMapper< Mapper, N > ThisType;
      typedef AdaptiveDofMapper< PowerMapperTraits< Mapper, N > > BaseType;

    public:
      typedef typename BaseType::Traits Traits;

      //! Type of Element (e.g. Codim 0 Entiy)
      typedef typename BaseType::ElementType ElementType;

      //! Type of indices
      typedef typename Traits::SizeType SizeType;

    protected:
      // funtor wrapper
      template< class Functor >
      struct FunctorWrapper
      {
        FunctorWrapper ( SizeType offset, Functor functor )
          : offset_( offset ),
            functor_( functor )
        {}

        template< class GlobalKey >
        void operator() ( int localBlock, const GlobalKey &globalKey )
        {
          int localDof = localBlock*numComponents;
          for( int component = 0; component < numComponents; ++component, ++localDof )
            functor_( localDof, globalKey + component * offset_ );
        }

        template< class GlobalKey >
        void operator() ( const GlobalKey &globalKey )
        {
          for( int component = 0; component < numComponents; ++component )
            functor_( globalKey + component * offset_ );
        }

      private:
        SizeType offset_;
        Functor functor_;
      };

      typedef typename Traits::MapperType MapperType;

      static const int numComponents = N;

    public:

      //! Constructor
      PowerMapper ( MapperType &mapper )
        : mapper_( mapper ),
          offset_( mapper.size() ),
          oldOffset_( -1 )
      {}

      //! Total number of degrees of freedom
      int size () const {  return mapper().size() * numComponents; }

      /** \copydoc Dune::DofMapper::contains( const int codim ) const */
      bool contains ( const int codim ) const
      {
        return mapper().contains( codim );
      }

      /** \copydoc Dune::DofMapper::fixedDataSize( const int codim ) const */
      bool fixedDataSize ( const int codim ) const
      {
        return mapper().fixedDataSize( codim );
      }

      /** \copydoc DofMapper::mapEach */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        FunctorWrapper< Functor > wrapper( offset_, f );
        mapper().mapEach( element, wrapper );
      }

      /** \copydoc Dune::DofMapper::mapEntityDofToGlobal(const Entity &entity,const int localDof) const */
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        FunctorWrapper< Functor > wrapper( offset_, f );
        mapper().mapEachEntityDof( entity, wrapper );
      }

      /** \copydoc Dune::DofMapper::maxNumDofs() const */
      int maxNumDofs () const { return mapper().maxNumDofs() * numComponents; }

      /** \copydoc Dune::DofMapper::numDofs(const ElementType &element) const */
      int numDofs ( const ElementType &element ) const { return mapper().numDofs( element ) * numComponents; }

      /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
      template< class Entity >
      int numEntityDofs ( const Entity &entity ) const { return mapper().numEntityDofs( entity ) * numComponents; }

      //! return new index in dof array
      int newIndex ( const int hole, const int block ) const
      {
        const int numContainedBlocks = mapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffset = mapper().newIndex( hole, containedBlock );
        return containedOffset + component * offset_;
      }

      //! return old index in dof array of given index ( for dof compress )
      int oldIndex ( const int hole, const int block ) const
      {
        const int numContainedBlocks = mapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffset = mapper().oldIndex( hole, containedBlock );
        return containedOffset + component * oldOffset_;
      }

      //! return number of holes in the data
      int numberOfHoles ( const int block ) const { return mapper().numberOfHoles( block % mapper().numBlocks() ); }

      //! returnn number of mem blocks
      int numBlocks () const { return mapper().numBlocks() * numComponents; }

      //! return current old offset of block
      int oldOffSet ( const int block ) const
      {
        const int numContainedBlocks = mapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffset = mapper().oldOffSet( containedBlock );
        return containedOffset + component * oldOffset_;
      }

      //! return current offset of block
      int offSet ( const int block ) const
      {
        const int numContainedBlocks = mapper().numBlocks();
        const int containedBlock     = block % numContainedBlocks;
        const int component          = block / numContainedBlocks;

        const int containedOffset = mapper().offSet( containedBlock );
        return containedOffset + component * offset_;
      }

      //! return true if compress will affect data (always true for this mapper)
      bool consecutive () const { return true; }

      template< class Entity >
      void insertEntity ( const Entity &entity ) { update(); }

      template< class Entity >
      void removeEntity ( const Entity &entity ) { }

      void resize () { update(); }

      bool compress () { update(); return true; }

      template< class StreamTraits >
      void write ( OutStreamInterface< StreamTraits > &out ) const {}

      template< class StreamTraits >
      void read ( InStreamInterface< StreamTraits > &in )
      {
        update();
      }

      void backup () const {}
      void restore () {}

    protected:
      // store sizes for old-new index mapping
      void update ()
      {
        oldOffset_ = offset_;
        offset_ = mapper().size();
      }

      MapperType &mapper () { return mapper_; }
      const MapperType &mapper () const { return mapper_; }

    private:
      MapperType &mapper_;
      SizeType offset_, oldOffset_;
    };

    /** @} **/

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_COMBINEDSPACE_MAPPER_HH

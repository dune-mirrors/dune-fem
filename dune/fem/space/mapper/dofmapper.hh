#ifndef DUNE_FEM_DOFMAPPER_HH
#define DUNE_FEM_DOFMAPPER_HH

//- Dune includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/space/mapper/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    /** @addtogroup DofMapper

        Dof Mapper are special implementations to provide the mapToGlobal
        feature of the spaces. Furthermore, during adaptation the mapper
        provide information about holes in the data vectors.

        \remarks
        The DofMapper interface is described by the class
        DofMapper.

      @{
     */

    //-------------------------------------------------------------------------
    //-
    //-  --MapperInterface
    //-
    //-------------------------------------------------------------------------
    /** \brief
       Interface for calculating the size of a function space for a grid on a
       specified level.
       Furthermore the local to global mapping of dof number is done.
       Also during grid adaptation this mapper knows about old and new indices
       of entities.
    */
    template< class DofMapperTraits >
    class DofMapper
    : public Fem::BartonNackmanInterface< DofMapper< DofMapperTraits >,
                                     typename DofMapperTraits :: DofMapperType >
    {
      typedef DofMapper< DofMapperTraits > ThisType;
      typedef Fem::BartonNackmanInterface< ThisType, typename DofMapperTraits::DofMapperType > BaseType;

    public:
      typedef DofMapperTraits Traits;

      //! type of the DofMapper implementation
      typedef typename Traits::DofMapperType DofMapperType;

      //! type of codimension 0 entities
      typedef typename Traits::ElementType ElementType;

      //! type of size integer
      typedef typename Traits::SizeType SizeType;

      typedef ElementType EntityType;

    protected:
      using BaseType::asImp;

    public:
      /** \brief  return number of dofs for special function space and grid on
                  specified level  */
      SizeType size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
        return asImp().size();
      }

      /** \brief returns true if DoFs for given codimension exist

          \param[in]  codim   codimension to check

          \returns true if DoFs for codimension exist
       */
      bool contains ( const int codim ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().contains( codim ) );
        return asImp().contains( codim );
      }

      /** \brief Check, whether the data in a codimension has fixed size */
      bool fixedDataSize ( const int codim ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().fixedDataSize( codim ) );
        return asImp().fixedDataSize( codim );
      }

      /** \brief map each local DoF number to a global key

          \param[in]  element  element, the DoFs belong to
          \param[in]  f        functor to call for each DoF

          The functor has to be a copyable object satisfying the following
          interface:
          \code
          struct Functor
          {

            // application operator
            template< class GlobalKey >
            void operator() ( const int localDoF, const GlobalKey &globalDoF );
          };
          \endcode

          For each DoF to be mapped, this method will call the application operator
          once.

          \note There is no guarantee on the order, in which the functor is applied.
          \note The global key has to be compatible with the Dof storage.
       */
      template< class Functor >
      void mapEach ( const ElementType &element, Functor f ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().mapEach( element, f ) );
      }

      /** \brief map each local DoF number to a global key

          \param[in]  entity   entity, the DoFs belong to
          \param[in]  f        functor to call for each DoF

          The functor has to be a copyable object satisfying the following
          interface:
          \code
          struct Functor
          {

            // application operator
            template< class GlobalKey >
            void operator() ( const int localDoF, const GlobalKey &globalKey );
          };
          \endcode

          For each DoF to be mapped, this method will call the application operator
          once.

          \note There is no guarantee on the order, in which the functor is applied.
          \note The global key has to be compatible with the Dof storage.
       */
      template< class Entity, class Functor >
      void mapEachEntityDof ( const Entity &entity, Functor f ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().mapEachEntityDof( entity, f ) );
      }

      /** \brief obtain maximal number of DoFs on one entity
       */
      int maxNumDofs () const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().maxNumDofs());
        return asImp().maxNumDofs();
      }

      /** \brief obtain number of DoFs on an entity
       *
          \param[in]  element  entity of codimension 0

          \returns number of DoFs on the entity
       */
      SizeType numDofs ( const ElementType &element ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().numDofs( element ) );
        return asImp().numDofs( element );
      }

      /** \brief obtain number of DoFs actually belonging to an entity

          In contrast to numDofs, this method returns the number of DoFs actually
          associated with an entity (usually a subentity). We have the following
          relation for an entity \f$E\f$ of codimension 0:
          \f[
          \mathrm{numDofs}( E ) = \sum_{e \subset E} \mathrm{numEntityDofs}( e ),
          \f]
          where \f$\subset\f$ denotes the subentity relation.

          \param[in]  entity  entity of codimension

          \returns number of DoFs on the entity
       */
      template< class Entity >
      SizeType numEntityDofs ( const Entity &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().numEntityDofs( entity ) );
        return asImp().numEntityDofs( entity );
      }

      /**
       * \brief update DoF mapping after grid modification
       *
       * \note This method may not have any semantic side effects and may
       *       may be called arbitrarily often.
       *
       * \note If the update is expensive, implementors might choose to check
       *       the sequence number of the underlying grid part.
       **/
      void update ()
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().update() );
      }
    };



    //-------------------------------------------------------------------------
    //-
    //-  --AdaptiveMapperInterface
    //-
    //-------------------------------------------------------------------------
    /** \brief
       Extended interface for adaptive DoF mappers
    */
    template< class DofMapperTraits >
    class AdaptiveDofMapper
    : public DofMapper< DofMapperTraits >
    {
      typedef DofMapper< DofMapperTraits > BaseType;

    protected:
      using BaseType::asImp;

    public:
      //! type of size integer
      typedef typename BaseType::SizeType SizeType;

      //! at the moment this should be similar to SizeType
      typedef SizeType GlobalKeyType;

      /** \brief return number of holes for data block */
      SizeType numberOfHoles ( const int block ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().numberOfHoles(block));
        return asImp().numberOfHoles(block);
      }

      /** \brief return old index of hole for data block (with resprect to new offset) */
      GlobalKeyType oldIndex ( const int hole, const int block ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().oldIndex(hole,block));
        return asImp().oldIndex(hole,block);
      }

      /** \brief return new index of hole for data block (with resprect to new offset) */
      GlobalKeyType newIndex ( const int hole, const int block ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().newIndex(hole,block));
        return asImp().newIndex(hole,block);
      }

      /** \brief return true if compress will affect data */
      bool consecutive () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().consecutive() );
        return asImp().consecutive();
      }

      /** \brief return old offsets for given block */
      SizeType oldOffSet ( const int block ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().oldOffSet(block));
        return asImp().oldOffSet(block);
      }

      /** \brief return current offsets for given block */
      SizeType offSet ( const int block ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().offSet(block));
        return asImp().offSet(block);
      }

      /** \brief return number of supported blocks */
      SizeType numBlocks () const
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().numBlocks());
        return asImp().numBlocks();
      }

      /**
       * \brief update DoF mapping after grid modification
       *
       * Adaptive DoF mappers are considered to be always up to date and this
       * method does nothing.
       **/
      void update () {}
    };

    /// @}

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_HH

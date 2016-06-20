#ifndef DUNE_FEM_COMBINEDDOFSTORAGE_HH
#define DUNE_FEM_COMBINEDDOFSTORAGE_HH

//- local includes
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/storage/subvector.hh>

namespace Dune
{

  namespace Fem
  {

    //! Utility class that helps in the transformation between dofs in the
    //! combined space and its enclosed spaces
    template< class ContainedMapper, int N, DofStoragePolicy policy >
    class CombinedDofConversionUtility;


    //! does the same as DofConversionUtility<PointBased>, just other
    //! construtor
    template< class ContainedMapper , int N >
    class CombinedDofConversionUtility< ContainedMapper, N, PointBased >
    : public PointBasedDofConversionUtility< N >
    {
      typedef PointBasedDofConversionUtility< N >  BaseType;

    public:
      typedef ContainedMapper ContainedMapperType;

      CombinedDofConversionUtility ( const ContainedMapperType & mapper, const int numComponents )
      : BaseType( numComponents )
      {}
    };

    //! Specialisation for VariableBased approach
    template< class ContainedMapper, int N >
    class CombinedDofConversionUtility< ContainedMapper, N, VariableBased >
    {
    public:
      typedef ContainedMapper ContainedMapperType;

      /** \brief constructor
       *
       *  \param[in]  mapper  mapper of the contained space
       *  \param[in]  size    number of global DoFs per component
       */
      CombinedDofConversionUtility ( const ContainedMapperType &mapper, int size )
      : mapper_( mapper )
      {}

      //! Find out what type of policy this is.
      static DofStoragePolicy policy ()
      {
        return VariableBased;
      }

      //! Set new size after adaptation.
      void newSize ( int size )
      {}

      //! Component which the actual base function index gives a contribution
      //! \return is in range {0, dimRange-1}
      int component ( int combinedIndex ) const
      {
        return combinedIndex / containedSize();
      }

      //! Number of the (scalar) base function belonging to base function index
      int containedDof ( int combinedIndex ) const
      {
        return combinedIndex % containedSize();
      }

      //! Reverse operation of containedDof, component
      //! i == combinedDof(containedDof(i), component(i))
      int combinedDof ( int containedIndex, int component ) const
      {
        return containedIndex + (component * containedSize());
      }

    protected:
      const ContainedMapperType &mapper_;

      int containedSize () const
      {
        return mapper_.size();
      }
    };



    template< class MapperImp, int N, DofStoragePolicy policy  >
    class CombinedSubMapper
    : public Fem :: IndexMapperInterface< CombinedSubMapper< MapperImp, N, policy > >
    {
      typedef CombinedSubMapper< MapperImp, N , policy > ThisType;

    public:
      // original mapper
      typedef MapperImp ContainedMapperType;
      typedef CombinedDofConversionUtility< ContainedMapperType, N, policy > DofConversionType;

      CombinedSubMapper ( const ContainedMapperType& mapper, unsigned int component )
      : mapper_( mapper ),
        component_( component ),
        utilGlobal_(mapper_, policy  == PointBased ? N : mapper.size() )
      {
        assert(component_ < N);
      }

      CombinedSubMapper(const ThisType& ) = default;
      CombinedSubMapper(ThisType&& ) = default;
      ThisType& operator=(const ThisType& ) = delete;
      ThisType& operator=(ThisType&& ) = delete;

      //! Total number of degrees of freedom
      unsigned int size () const
      {
        return mapper_.size();
      }

      unsigned int range () const
      {
        return size() * N;
      }

      unsigned int operator[] ( unsigned int index ) const
      {
        utilGlobal_.newSize( mapper_.size() );
        return utilGlobal_.combinedDof(index, component_);
      }

    private:
      const ContainedMapperType& mapper_;
      const unsigned int component_;
      mutable DofConversionType utilGlobal_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMBINEDDOFSTORAGE_HH

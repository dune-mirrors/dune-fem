#ifndef DUNE_COMBINEDDOFSTORAGE_HH
#define DUNE_COMBINEDDOFSTORAGE_HH

//- local includes 
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/storage/subarray.hh>

namespace Dune
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
  public:
    typedef ContainedMapper ContainedMapperType;

  private:
    typedef PointBasedDofConversionUtility< N >  BaseType;

  public:
    inline CombinedDofConversionUtility ( const ContainedMapperType & mapper,
                                          const int numComponents )
    : BaseType( numComponents )
    {}
  };

  //! Specialisation for VariableBased approach
  template< class ContainedMapper, int N >
  class CombinedDofConversionUtility< ContainedMapper, N, VariableBased >
  {
  public:
    typedef ContainedMapper ContainedMapperType;

  protected:
    const ContainedMapperType &mapper_;

  public:
    /** \brief constructor
     *  
     *  \param[in]  mapper  mapper of the contained space
     *  \param[in]  size    number of global DoFs per component
     */
    inline CombinedDofConversionUtility ( const ContainedMapperType &mapper,
                                          int size )
    : mapper_( mapper )
    {}

    //! Find out what type of policy this is.
    inline static DofStoragePolicy policy ()
    {
      return VariableBased;
    }

    //! Set new size after adaptation.
    inline void newSize ( int size )
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
    inline int containedSize () const
    {
      return mapper_.size();
    }
  };


  
  template< class MapperImp, int N, DofStoragePolicy policy  >
  class CombinedSubMapper
  : public IndexMapperInterface< CombinedSubMapper< MapperImp, N, policy > >
  {
  public:
    //- original mapper 
    typedef MapperImp ContainedMapperType;

  private:
    typedef CombinedSubMapper< ContainedMapperType, N , policy > ThisType;

  public:
    typedef CombinedDofConversionUtility< ContainedMapperType, N, policy >   DofConversionType;

  public:
    //- Public methods
    CombinedSubMapper ( const CombinedMapperType &mapper,
                        const unsigned int component )
    : mapper_( mapper ),
      component_( component ),
      utilGlobal_(mapper_,
                  policy  == PointBased ? 
                  N :  mapper.size()/N )
    {
      assert(component_ < N);
    }

    CombinedSubMapper(const ThisType& other) :
      mapper_(other.mapper_),
      component_(other.component_),
      utilGlobal_(other.utilGlobal_) 
    { 
      assert(component_ < N);
    }

    //! Total number of degrees of freedom
    inline unsigned int size () const {
      return mapper_.size();
    }
    inline unsigned int range () const {
      return size() * N;
    }
    inline const unsigned int operator[] ( unsigned int index ) const 
    {
      utilGlobal_.newSize( mapper_.size() );
      return utilGlobal_.combinedDof(index, component_);
    }

  private:
    ThisType& operator=(const ThisType& other);
    //- Data members
    const ContainedMapperType& mapper_;
    const unsigned int component_;
    mutable DofConversionType utilGlobal_;
  };

} // end namespace Dune

#endif

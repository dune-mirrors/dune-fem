#ifndef DUNE_FEM_REFERENCEBLOCKVECTOR_HH
#define DUNE_FEM_REFERENCEBLOCKVECTOR_HH

#include <vector>

#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

namespace Dune {
namespace Fem {

  /** \class ReferenceBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< typename F, unsigned int BlockSize >
  class ReferenceBlockVector : public MutableBlockVector< std::vector< F >, BlockSize >
  {
    typedef MutableBlockVector< std::vector< F >, BlockSize > BaseType;

  public:
    typedef typename BaseType::SizeType  SizeType;

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  The dofs are not initialized.
     *
     *  \param[in]  size         Number of blocks
     */
    explicit ReferenceBlockVector ( SizeType size )
    : BaseType( size )
    {}

  };

} // namespace Fem
} // namespace Dune

#endif // DUNE_FEM_REFERENCEBLOCKVECTOR_HH

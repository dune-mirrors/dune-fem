// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_REFERENCEBLOCKVECTOR_HH
#define DUNE_FEM_REFERENCEBLOCKVECTOR_HH

#include <algorithm>
#include <cassert>
#include <vector>

#include <dune/fem/function/blockvectors/simpleblockvector.hh>

namespace Dune {
namespace Fem {

  // Forward declaration
  template< typename F, unsigned int BlockSize >
  class ReferenceBlockVectorBlock;

  /** \class ReferenceBlockVector
  *   \brief This is the reference implementation of a block vector as it is expected
  *      as the second template parameter to Dune::Fem::BlockVectorDiscreteFunction
  *
  *   \tparam  F           The ground fields. All dofs are elements of this field.
  *   \tparam  BlockSize   Size of the blocks
  */
  template< typename F, unsigned int BlockSize >
  class ReferenceBlockVector : public SimpleBlockVector< std::vector< F >, BlockSize >
  {
    typedef SimpleBlockVector< std::vector< F >, BlockSize > BaseType;

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

    /** \brief Constructor; use this to create a block vector with 'size' blocks.
     *
     *  All the dofs are set to 'initialValue'.
     *
     *  \param[in]  size          Number of blocks
     *  \param[in]  initialValue  This is the value to which each dof is set
     */
    ReferenceBlockVector ( SizeType size, const F& initialValue )
    : BaseType( size, initialValue )
    {}
  };

} // namespace Fem
} // namespace Dune

#endif // DUNE_FEM_REFERENCEBLOCKVECTOR_HH

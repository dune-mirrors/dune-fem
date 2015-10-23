#ifndef DUNE_FEM_RANGEGENERATORS_HH
#define DUNE_FEM_RANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>

namespace Dune
{
  namespace Fem
  {
    /**
     * \brief Iterator ranges for entities and DOFs to support iteration with range-based for loops.
     *
     * <h2>Range-based for loop</h2>
     *
     * Range-based for loops are availabe in GCC 4.6+, Clang 3.2+ and Intel ICC 13+.
     *
     * For further details, see e.g. http://en.cppreference.com/w/cpp/language/range-for.
     *
     * <h2>Entities</h2>
     *
     * Assuming you have a DiscreteFunction `df`, you can iterate over the grid entities like this:
     *
     * \code
     * // iterate over entities
     * for (auto&& entity : entities(df))
     * {
     *   // do stuff
     * }
     * \endcode
     *
     * <h2>DOFs</h2>
     *
     * If you want to iterate over the DOFs of a DiscreteFunction `df`, you can use the function dofs()
     * to obtain a range that is suitable for a range-based for loop. As before you can do like this:
     *
     * \code
     * for (auto&& dof : dofs(df))
     * {
     *   // do stuff
     * }
     * \endcode
     *
     */

    //! Iterates over all entities.
    /**
     * This functions returns an object representing the range of entities. The main purpose of this function
     * is to enable iteration over those entities by means of a range-based for loop.
     *
     * \since      GCC 4.6
     * \param df   a DiscreteFunction.
     * \returns    an unspecified object that is guaranteed to fulfil the interface
     *             of IteratorRange and that can be iterated over using a range-based
     *             for loop.
     */
    template<typename DF>
    inline IteratorRange<typename DF::DiscreteFunctionSpaceType::IteratorType> entities(const DF& df)
    {
      typedef IteratorRange<typename DF::DiscreteFunctionSpaceType::IteratorType> ReturnType;
      return ReturnType(df.space().begin(),df.space().end());
    }
    //! \}

    //! Iterates over all DOFs.
    /**
     * This functions returns an object representing the range of DOFs. The main purpose of this function
     * is to enable iteration over those DOFs by means of a range-based for loop.
     *
     * \since      GCC 4.6
     * \param df   a DiscreteFunction.
     * \returns    an unspecified object that is guaranteed to fulfil the interface
     *             of IteratorRange and that can be iterated over using a range-based
     *             for loop.
     */
    template<typename DF>
    inline IteratorRange<typename DF::DofIteratorType> dofs(DF& df)
    {
      typedef IteratorRange<typename DF::DofIteratorType> ReturnType;
      return ReturnType(df.dbegin(),df.dend());
    }
    //! \}

    //! Iterates over all DOFs.
    /**
     * This functions returns an object representing the range of DOFs. The main purpose of this function
     * is to enable iteration over those DOFs by means of a range-based for loop.
     *
     * \since      GCC 4.6
     * \param df   a DiscreteFunction.
     * \returns    an unspecified object that is guaranteed to fulfil the interface
     *             of IteratorRange and that can be iterated over using a range-based
     *             for loop.
     */
    template<typename DF>
    inline IteratorRange<typename DF::ConstDofIteratorType> dofs(const DF& df)
    {
      typedef IteratorRange<typename DF::ConstDofIteratorType> ReturnType;
      return ReturnType(df.dbegin(),df.dend());
    }
    //! \}

  } // end namespace Fem

} // end namespace Dune

#endif // DUNE_FEM_RANGEGENERATORS_HH

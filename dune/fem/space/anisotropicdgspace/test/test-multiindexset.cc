#include <config.h>

// C++ includes
#include <iostream>

// dune-fem includes
#include <dune/fem/space/anisotropicdgspace/multiindexset.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Provide a test for class MultiIndexSet
*/


int main ( int argc, char **argv )
{
  // get dimension and order
  const int dimension = Dune::GridSelector::GridType::dimension;
  const int order = POLORDER;

  // create multi index set
  typedef AnisotropicDG::MultiIndexSet< dimension, order > MultiIndexSet;
  MultiIndexSet multiIndexSet;

  // get size of index set
  std::size_t size = multiIndexSet.size();

  // iterate over all indices and print on konsole
  std::cout << "Printing all " << size << " indices..." << std::endl;
  typedef MultiIndexSet::IteratorType IteratorType;
  std::size_t count = 0;
  const IteratorType end = multiIndexSet.end();
  for( IteratorType it = multiIndexSet.begin(); it != end; ++it )
  {
    // get multi index
    MultiIndexSet::MultiIndexType multiIndex = *it;

    // check whether index is valid
    if( !multiIndexSet.contains( multiIndex ) )
    {
      std::cerr << "Error, multiIndex not contained in index set" << std::endl;
      return 1;
    }

    // print index
    std::cout << "  " << multiIndex << std::endl;
    ++count;
  }

  // check number of indices returned by iterator
  if( count != size )
  {
    std::cerr << "Error, size() and number of indices returned do not coincide" << std::endl;
    return 1;
  }

  std::cout << "...done" << std::endl;

  return 0;
}

// vim: set expandtab ts=2 sw=2 sts=2:
#include <config.h>

#include <iostream>
#include <string>

#include <dune/fem/petsc/common/petsccommon.hh>
#include <dune/fem/petsc/discretefunction/petscdiscretefunction.hh>

#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>

#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/space/common/adaptmanager.hh>

/*
 * global constants/variables
 */
#if defined POLORDER
const int polOrder = POLORDER;
#else
const int polOrder = 1;
#endif

/*
 * global types
 */
typedef Dune::GridSelector::GridType GridType;
typedef Dune::LeafGridPart< GridType > GridPartType;
typedef Dune::FunctionSpace< double, double, GridType::dimensionworld, 1 > FunctionSpaceType;
typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpaceType;
typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > PetscDiscreteFunctionType;

//typedef Dune::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType1;
typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType1;

class VariableFilenameParameter
: public Dune::DataOutputParameters 
{

public:
  VariableFilenameParameter ( const std::string &prefix ) 
  : prefix_( prefix )
  {}

  virtual std::string prefix () const
  {
    return prefix_;
  }
  
private:
  std::string prefix_;
};

template< typename DF >
void setToRank ( DF &dFunction, const std::string &msg ) 
{
  typedef Dune::tuple< DF* > IOTupleType;
  typedef Dune::DataOutput< GridType, IOTupleType > DataOutputType;
  IOTupleType dataTup ( &dFunction );
  DataOutputType dataOutput( dFunction.gridPart().grid(), dataTup, VariableFilenameParameter( msg ) );

  std::cout << std::endl;
  std::cout << dFunction.gridPart().comm().rank() << ", size() of " << msg << ": " << dFunction.size() << std::endl;
  /*std::cout << dFunction.gridPart().comm().rank() << ", size/localBlockSize for " << 
            msg << ": " << dFunction.size()/DF::localBlockSize << std::endl;*/
  size_t counter = 0;
  for( typename DF::DofIteratorType it = dFunction.dbegin(); it != dFunction.dend(); ++it ) 
  {
    *it = dFunction.gridPart().comm().rank() + 1;
    ++counter;
  }
  std::cout << dFunction.gridPart().comm().rank() << ", " << msg << "s dof iterator walked over " 
            << counter << " dofs\n";
  dataOutput.writeData( 0, "before" );
  dFunction.communicate();
  dataOutput.writeData( 1, "after" );

  /*std::cout << msg << ", reading on rk " << dFunction.gridPart().comm().rank() << " after communication:\n";
  counter = 0;
  for( typename DF::DofIteratorType it = dFunction.dbegin(); it != dFunction.dend(); ++it ) 
  {
    std::cout << counter++ << " >> " << *it << std::endl;
  }*/
}

int main ( int argc, char ** argv ) 
{
  try
  {
    Dune::MPIManager::initialize( argc, argv );
    Dune::Petsc::PetscInitialize( &argc, &argv, static_cast< char* >( 0 ), static_cast< char* >( 0 ) );

    Dune::Parameter::append( argc, argv );
    Dune::Parameter::append( "parameter_commtest" );

    const std::string gridFileName = Dune::Parameter::getValue< std::string >( "gridfilename" );
    Dune::GridPtr< GridType > gridPtr( "domain_commtest.dgf" );

    GridPartType gridPart( *gridPtr );
    DiscreteFunctionSpaceType dfSpace( gridPart, Dune::All_All_Interface );
    //GridType &grid = gridPart.grid();


    // create the other discrete function, do stuff with it..
    DiscreteFunctionType1 dFunction1( "dFunction1", dfSpace );
    dFunction1.clear();
    dFunction1.communicate();
    setToRank ( dFunction1, "other" );


    // create the PETSc discrete function
    PetscDiscreteFunctionType dFunctionPetsc( "PETScDiscreteFunction", dfSpace ); 
    dFunctionPetsc.clear();
    dFunctionPetsc.communicate();
    setToRank( dFunctionPetsc, "PETSc" );

    Dune::Parameter::write( "parameter_commtest.log" );

    Dune::Petsc::PetscFinalize();

    return 0;
  }
  catch( const Dune::Exception &exception )
  {
    std::cerr << exception << std::endl;
    return 1;
  }

  // control should not reach here
  abort();
  return 0;
}

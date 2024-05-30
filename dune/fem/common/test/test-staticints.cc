#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/common/staticlistofint.hh>
#include <dune/fem/io/parameter.hh>

LIST_OF_INT(Solvers,
            cg = 1,
            bicgstab = 2,
            gmres = 3 );

int main (int argc, char **argv)
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  Dune::Fem::Parameter::append("fem.solver.method","bicgstab");

  std::string prefix ("fem.solver");
  int mth = Solvers::to_id( Dune::Fem::Parameter::getEnum( prefix + ".method", Solvers::names() ) );
  std::string name = Solvers::to_string( mth );
  std::cout << "Solver id   " << mth << std::endl;
  std::cout << "Solver name " << name << std::endl;

  // we should get bicgstab
  assert( mth == 2 );
  assert( name == std::string("bicgstab"));
  return 0;
}

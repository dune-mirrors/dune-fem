#include <config.h>
#include <dune/common/fvector.hh>
#include <iostream>
#include <dune/grid/alugrid.hh>
typedef Dune::ALUSimplexGrid<3,3> GridType;
typedef GridType :: Codim<0> :: LeafIterator LeafIterator;
/****************************************************************
 Type IDiscFuncType 
 Oeffentliches Interface:
   typedef  LocalDiscFuncInterfaceType
   typedef  EntityType
   enum     dimDomain
   typedef  DomainType
   method   double& operator[](int)
   method   double operator[](int) const
   method   int size() const 
   method   LocalDiscFuncInterfaceType* operator[](const EntityType&)
   method   double evaluate(const DomainType&) const 
   method   int map(const EntityType&,int) const 

 Local Function Interface (erfuellt von IDiscFuncType::LocalDiscFuncInterfaceType
   typedef EntityType
   enum    dimDomain 
   typedef DomainType
   method  double evaluate(const DomainType&)
   method  int size()
   method  double& operator[](int i)
   method  double operator[](int i) const 

****************************************************************/
#include "discfuncinterface.hh"
/****************************************************************
 ****************************************************************/
void testinterface(IDiscFuncType& df) {
  const int size = df.size();       
  for (int k=0;k<100000;++k) {
    for (int i=0;i<size;++i) {   // warum nicht for (int i=0;i<df.size();++i)?
      df[i] *= -1.0;
    }
  }
}
void testlocalinterface(GridType& grid,IDiscFuncType& df) {
  Dune::FieldVector<double,3> x[10];
  x[0][0] = 0 , x[0][1] = 0 , x[0][2] = 0;
  x[1][0] = 1 , x[1][1] = 0 , x[1][2] = 0;
  x[2][0] = 0 , x[2][1] = 1 , x[2][2] = 0;
  x[3][0] = 0 , x[3][1] = 0 , x[3][2] = 1;
  x[4][0] = 0.5 , x[4][1] = 0   , x[4][2] = 0;
  x[5][0] = 0   , x[5][1] = 0.5 , x[5][2] = 0;
  x[6][0] = 0.5 , x[6][1] = 0.5 , x[6][2] = 0;
  x[7][0] = 0   , x[7][1] = 0   , x[7][2] = 0.5;
  x[8][0] = 0   , x[8][1] = 0.5 , x[8][2] = 0.5;
  x[9][0] = 0.5 , x[9][1] = 0.5 , x[9][2] = 0.5;
  LeafIterator endit = grid.leafend();
  for (int k=0;k<100;++k) {
    for(LeafIterator it = grid.leafbegin();it != endit ; ++it ) {
      IDiscFuncType::LocalDiscFuncInterfaceType* ldf = df[*it];
      for (int i=0;i<10;++i) {
	(*ldf).evaluate(x[i]);
      }
      delete ldf;
    }
  }
}
// ***************************************************
double f(const Dune::FieldVector<double,3>& x) {
  return x[0]*x[1]*x[2];
}
int main(int argc, char ** argv, char ** envp) {
  std::cout << "HALLO" << std::endl;
  // Gitter anlegen, verfeinern und diskrete Funktion initialisieren
  GridType grid("alu-testgrid.tetra");
  grid.globalRefine(3);
  ContDiscFunc df(grid);         
  LeafIterator endit = grid.leafend();
  for(LeafIterator it = grid.leafbegin();it != endit ; ++it ) {
    IDiscFuncType::LocalDiscFuncInterfaceType* ldf = df[*it];
    for (int i=0;i<(*ldf).size();++i) {
      (*ldf)[i] = f(it->geometry()[i]);
    }
    delete ldf;
  }
  // Effizienz der Interfaces testen (kein Test der Interface effektivitaet!)
  testinterface(df);
  testlocalinterface(grid,df);
}

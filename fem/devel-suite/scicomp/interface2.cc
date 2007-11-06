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
   method   LocalDiscFuncInterfaceType operator[](const EntityType&)
   method   double evaluate(const DomainType&) const 
   method   int map(const EntityType&,int) const 

 Type Local FuncInterface
   typedef EntityType
   enum    dimDomain 
   typedef DomainType
   method  double evaluate(const DomainType&)
   method  int size()
   method  double& operator[](int i)
   method  double operator[](int i) const 

****************************************************************/
#include "discfuncinterface2.hh"
/****************************************************************
 ****************************************************************/
void testinterface(IDiscFuncType& df) {
  const int size = df.size();       // sicherste Variante!
  for (int k=0;k<100000;++k) {
    //for (int i=0;i<df.size();++i) {
    for (int i=0;i<size;++i) {
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
      IDiscFuncType::LocalDiscFuncInterfaceType  ldf = df[*it];
      for (int i=0;i<10;++i) {
	ldf.evaluate(x[i]);
      }
    }
  }
}
// ***************************************************
double f(const Dune::FieldVector<double,3>& x) {
  return x[0]*x[1]*x[2];
}
int main(int argc, char ** argv, char ** envp) {
  std::cout << "HALLO" << std::endl;
  GridType grid("alu-testgrid.tetra");
  grid.globalRefine(3);
  ContDiscFunc df(grid);         // Eine Moegliche Implementierung
  LeafIterator endit = grid.leafend();
  for(LeafIterator it = grid.leafbegin();it != endit ; ++it ) {
    IDiscFuncType::LocalDiscFuncInterfaceType ldf = df[*it];
    for (int i=0;i<ldf.size();++i) {
      ldf[i] = f(it->geometry()[i]);
    }
  }
  testinterface(df);
  testlocalinterface(grid,df);
  // df.disp();
}

/******************************************************************
 * Vorteil:
 * Compiler kann viel mehr optimieren - vor allem inlining. Interface
 * Methoden (sollten) komplett wegfallen, d.h. man erhaellt direkt
 * die gewuenschte Methode und das bei der Uebersetzung.
 * Nachteil:
 * 
 ******************************************************************/


// 0)   naive                     144.736u 0.651s 2:26.82  
// 1)   wrapper                   104.532u 0.387s 1:45.94       1/0  >-25%
// 1a)  wrapper mit size Aufruf   128.434u 0.460s 2:10.25      1a/1  >+20%
// 2)   BN                         70.590u 0.398s 1:11.84       2/1  >-30%
// 2a)  BN (mit Caching)           35.586u 0.306s 0:36.29      2a/1  >-65%
//                                                             2a/2  >-50%
// 3)   BN/Engine                  33.975u 0.329s 0:34.70   



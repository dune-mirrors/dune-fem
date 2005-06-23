#ifndef __MAKEGRID_CC__
#define __MAKEGRID_CC__

#include <cassert> 
//***********************************************************
//
//  --AlbertaGrid
//
//***********************************************************
#if AGRID 
// * change here
typedef AlbertaGrid<DIM,DIM_OF_WORLD> GrdType;

template <class GType>
class ParaCrit 
{
  typedef GType GridType;
  typedef typename GridType::template Traits<0>::Entity EntityType;
  enum { dim = EntityType::dimension };
  typedef BaryCenterQuad < double , FieldVector<double,dim> , 0 > QuadratureType;
  
  FieldVector<double,dim> bary_;
 
  GridType &grid_;
  QuadratureType quad_;
  
  // number of processors
  int procNum_;
  int myGridSize_;
  int count_;

  int proc_;
  
public:
  ParaCrit(GridType &grid , int procNum ) : grid_(grid), 
      quad_ ( *(grid.template lbegin<0> (0))), procNum_ (procNum) , 
      myGridSize_ ( grid.size( grid.maxlevel(), 0)) , count_ (0) , proc_(0)
      {
        if((myGridSize_ % procNum_) != 0)
          myGridSize_ = (int) myGridSize_ % (procNum_-1);
        
        myGridSize_ = myGridSize_/procNum_;
        std::cout << myGridSize_ << " Size of Partition \n\n";
      }
  
  void classify(EntityType &en)
  {
    if(procNum_ == 1) 
    {
      std::cout << "classify En " << en.global_index() << "\n";
      assert(proc_ < procNum_);
      en.partition( proc_ );
      count_++;
      if(count_ >= myGridSize_) 
      {
        proc_++;
        count_ = 0;
      }
    }
   
    if(procNum_ == 2)
    {
      bary_ = en.geometry().global(quad_.point(0));
      if(bary_[0] < 150)
      {
        //std::cout << "make partition 0 \n";
        en.partition( 0 );
        return;
      }
      if(bary_[0] < 300)
      {
        //std::cout << "make partition 1 \n";
        en.partition( 1 );
        return;
      }
    }

    if(procNum_ >2 )
    {
   
    bary_ = en.geometry().global(quad_.point(0));
    if(bary_[0] < 75)
    {
      en.partition( 0 );
      return;
    }
    if(bary_[0] < 150)
    {
      en.partition( 1 );
      return;
    }
    if(bary_[0] < 225)
    {
      en.partition( 2 );
      return;
    }
   
    if(bary_[0] < 300)
    {
      en.partition( procNum_-1 );
      return;
    }
    
    }
    /*
    bary_ = en.geometry().global(quad_.point(0));
    if(bary_(1) < 7.5)
    {
      en.partition( 0 );
      return;
    }
    if(bary_(1) < 15.0)
    {
      en.partition( 1 );
      return;
    }
    if(bary_(1) < 22.5)
    {
      en.partition( 2 );
      return;
    }
    
    if(bary_(1) < 30.0)
    {
      en.partition( 3 );
      return;
    }
    */

    /*
    bary_ = en.geometry().global(quad_.point(0));
    if(bary_(1) < 15.)
    {
      en.partition( 0 );
      return;
    }
    if(bary_(1) < 30.0)
    {
      en.partition( 1 );
      return;
    }
    */
  }
};

GrdType *
makeGrid (const char * gridName , int myRank, int mySize )
{
  GrdType * grd = 0;
#ifdef _PARALLEL_
  grd = new GrdType ( gridName , myRank );
  assert(grd);
  ParaCrit<GrdType> crit ( *grd , mySize );
  makeParallelGrid( *grd , crit );
  grd->createGhosts();
#else
  grd = new GrdType ( gridName );
  std::cout << "got here" << std::endl;
  assert(grd);
  std::cout << "too..." << std::endl;
#endif
  return grd;
}

#endif


//***********************************************************
//
//  --ALU3dGrid
//
//***********************************************************

#if BGRID
typedef ALU3dGrid<DIM,DIM_OF_WORLD,tetra> GrdType;

GrdType *
makeGrid (const char * gridName , int myRank, int mySize )
{
  GrdType * grd = 0;
#ifndef _ALU3DGRID_PARALLEL_
  grd = new GrdType ( gridName );
#else
  grd = new GrdType ( gridName , MPI_COMM_WORLD );
#endif
  assert(grd);
  grd->loadBalance();
  return grd; 
}
#endif


#if SGRID
typedef SGrid<DIM,DIM_OF_WORLD> GrdType;

GrdType *
makeGrid (const char * gridName , bool parallel )
{
  GrdType * grd = 0
#ifdef YGRID    
    Vec<DIM> lang = 30.0;
    lang(0) = 300.0;
    Vec<DIM,int> anz = 1;
    anz(0) = 10;
    Vec<DIM,bool> per = false;
    grd = new GrdType (MPI_COMM_WORLD, lang, anz, per , 1 );
#else 
    int n[DIM];
    double h[DIM];

#ifdef TWOPHASE 
    for(int i=0; i<DIM; i++) 
    {
      n[i] = 8; h[i] = 0.5;
    }
    n[1] = 16;
    h[1] = 1.0;
#else 
    n[0] = 10;
    h[0] = 300.0;
    for(int i=1; i<DIM; i++) 
    {
      n[i] = 1; h[i] = 30.0;
    }
#endif 
    std::cout << "Set adaptive to false, because using SGrid. \n";
    grd = new GrdType ((int *) &n, (double *) &h );
#endif
  assert(grd);  
  return grd;
}
    
#endif

#endif


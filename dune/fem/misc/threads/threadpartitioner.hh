#ifndef DUNE_FEM_THREADPARTITIONER_HH
#define DUNE_FEM_THREADPARTITIONER_HH

//- system includes 
#include <string>
#include <list>
#include <map>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#if HAVE_ALUGRID 
#include <dune/grid/alugrid/3d/alugrid.hh>
#elif HAVE_DUNE_ALUGRID
#include <dune/alugrid/3d/alugrid.hh>
#endif

#if ALU3DGRID_PARALLEL

#warning "Using the ThreadPartitioner"

#if HAVE_DUNE_ALUGRID 
namespace ALUGrid {
#else
namespace ALUGridSpace {
#endif

  // class fulfilling the ALUGrid communicator interface
  // but without any communication, this is needed to avoid 
  // communication during the call of DataBase :: repartition 
  class MpAccessSerial : public MpAccessGlobal
  {
  public:  
    MpAccessSerial() {} 

#if HAVE_DUNE_ALUGRID || defined ALUGRID_CONSTRUCTION_WITH_STREAMS 
    // this features is only available in the newer version 
    typedef MpAccessGlobal :: minmaxsum_t  minmaxsum_t;
#else 
    typedef double minmaxsum_t ;
#endif

    virtual int psize() const { return 1; }
    int myrank () const { return 0; }

    virtual int barrier () const { return psize(); }
    virtual bool gmax (bool value) const { return value; }
    virtual int gmax (int value) const { return value; }
    virtual int gmin ( int value ) const { return value; }
    virtual int gsum ( int value ) const { return value; }
    virtual long gmax ( long value ) const { return value; }
    virtual long gmin ( long value ) const { return value; }
    virtual long gsum ( long value ) const { return value; }
    virtual double gmax ( double value ) const { return value; }
    virtual double gmin ( double value ) const { return value; }
    virtual double gsum ( double value ) const { return value; }
    virtual void gmax (double* in, int length, double* out) const { std::copy(in, in+length, out ); }
    virtual void gmin (double* in, int length, double* out) const { std::copy(in, in+length, out ); }
    virtual void gsum (double* in, int length, double* out) const { std::copy(in, in+length, out ); }
    virtual void gmax (int*,int,int*) const {  }
    virtual void gmin (int*,int,int*) const {  }
    virtual void gsum (int*,int,int*) const {  }
    virtual minmaxsum_t minmaxsum( double value ) const { return minmaxsum_t( value ); }
    virtual std::pair<double,double> gmax (std::pair<double,double> value ) const { return value; }
    virtual std::pair<double,double> gmin (std::pair<double,double> value ) const { return value; }
    virtual std::pair<double,double> gsum (std::pair<double,double> value ) const { return value; }
    virtual void bcast(int*,int, int) const { }
    virtual void bcast(char*,int, int) const { }
    virtual void bcast(double*,int, int) const { }
    virtual int exscan( int value ) const { return 0; }
    virtual int scan( int value ) const { return value; }
    virtual std::vector < int > gcollect ( int value ) const { return std::vector<int> (psize(), value); }
    virtual std::vector < double > gcollect ( double value ) const { return std::vector<double> (psize(), value); }
    virtual std::vector < std::vector < int > > gcollect (const std::vector < int > & value) const 
    { 
      return std::vector < std::vector < int > > (psize(), value); 
    }
    virtual std::vector < ObjectStream > gcollect (const ObjectStream &, const std::vector<int>& ) const
    {
      return std::vector < ObjectStream > (psize()); 
    }
    virtual std::vector < std::vector < double > > gcollect (const std::vector < double > & value) const 
    { 
      return std::vector < std::vector < double > > (psize(), value); 
    }
    virtual std::vector < ObjectStream > gcollect (const ObjectStream &os) const { return std::vector < ObjectStream >(psize(),os); }

#if HAVE_ALUGRID
    //! return address of communicator (not optimal but avoid explicit MPI types here)
    virtual const CommIF* communicator() const { return 0; }
#endif
  };
}

namespace Dune {

template < class GridPartImp >
class ThreadPartitioner 
{

protected:  
  typedef GridPartImp GridPartType;
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;
  typedef typename LocalIdSetType :: IdType IdType;

  typedef typename GridPartType :: IndexSetType  IndexSetType;
  typedef typename IndexSetType :: IndexType IndexType;

protected:  
  typedef ALU3DSPACE LoadBalancer LoadBalancerType;
  typedef typename LoadBalancerType :: DataBase DataBaseType;

  // type of communicator interface 
  typedef ALU3DSPACE MpAccessLocal MPAccessInterfaceType;

  // type of communicator implementation 
  typedef ALU3DSPACE MpAccessSerial  MPAccessImplType;

  mutable MPAccessImplType  mpAccess_;

  DataBaseType db_;

  const GridPartType& gridPart_;
  const IndexSetType& indexSet_;

  typedef typename GridPartType :: template Codim<0> :: EntityType         EntityType;
  typedef typename GridPartType :: template Codim<0> :: EntityPointerType  EntityPointerType;

  // load balancer bounds 
  const double ldbOver_ ;
  const double ldbUnder_;
  const double cutOffFactor_;

  const int pSize_ ;
  int graphSize_; 

  std::vector< int > index_;
  int indexCounter_ ;

  std::vector< int > partition_;

public:
  enum Method { recursive = 0, // METIS_PartGraphRecursive,
                kway      = 1, // METIS_PartGraphKway
                sfc       = 2  // ALUGRID_SpaceFillingCurve 
  };

  /** \brief constructor 
      \param gridPart  grid part with set of entities that should be partitioned 
      \param pSize     number of partitions 
  */    
  ThreadPartitioner( const GridPartType& gridPart, 
                     const int pSize,
                     const double cutOffFactor = 1.0 )
    : mpAccess_(),
      db_ (),
      gridPart_( gridPart )
    , indexSet_( gridPart_.indexSet() )
    , ldbOver_(1.2)
    , ldbUnder_(0.0)
    , cutOffFactor_( std::max( double(1.0 - cutOffFactor), 0.0 ) )
    , pSize_( pSize )
    , graphSize_( pSize_ )
    , index_( indexSet_.size( 0 ), -1 )
    , indexCounter_( 0 )
  {
    calculateGraph( gridPart_ );
  }

protected:  
  //! create consecutive entity numbering on-the-fly 
  //! this is neccessary, because we might only be interating over a
  //! sub set of the given entities, and thus indices might be non-consecutive
  int getIndex( const size_t idx ) 
  {
    assert( idx < index_.size() );
    if( index_[ idx ] < 0 ) index_[ idx ] = indexCounter_ ++ ;
    return index_[ idx ] ;
  }

  //! access consecutive entity numbering (read-only)
  int getIndex( const size_t idx ) const 
  {
    assert( idx < index_.size() );
    return index_[ idx ] ;
  }

  //! get consecutive entity index, if not existing, it's created 
  int getIndex( const EntityType& entity ) 
  {
    return getIndex( indexSet_.index( entity ) );
  }

  //! get consecutive entity index (read-only)
  int getIndex( const EntityType& entity ) const 
  {
    return getIndex( indexSet_.index( entity ) );
  }

  void calculateGraph( const GridPartType& gridPart )
  {
    graphSize_ = 0;
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType Iterator;
    const Iterator end = gridPart.template end<0> ();
    const int cutOff = cutOffFactor_ * (indexSet_.size( 0 ) / pSize_) ;
    // create graph 
    for(Iterator it = gridPart.template begin<0> (); it != end; ++it )
    {
      const EntityType& entity = *it;
      assert( entity.partitionType() == InteriorEntity );
      ldbUpdateVertex ( entity, cutOff,
                        gridPart.ibegin( entity ),
                        gridPart.iend( entity ),
                        db_ );
    }
  }

  template <class IntersectionIteratorType>
  void ldbUpdateVertex ( const EntityType & entity, 
                         const int cutOff,
                         const IntersectionIteratorType& ibegin,
                         const IntersectionIteratorType& iend,
                         DataBaseType & db )
  {
    const int index = getIndex( entity );
    int weight = (index >= cutOff) ? 1 : 8; // a least weight 1 for macro element 

    {
      if( Fem::GridPartCapabilities::hasGrid< GridPartType >::v )
      {
        // calculate weight, which is number of children
        const int mxl = gridPart_.grid().maxLevel();
        if( mxl > entity.level() && ! entity.isLeaf() ) 
        {
          typedef typename EntityType :: HierarchicIterator HierIt; 
          const HierIt endit = entity.hend( mxl ); 
          for(HierIt it = entity.hbegin( mxl ); it != endit; ++it)
            ++weight;
        }
      }

#if HAVE_ALUGRID // old ALUGrid version 1.52
      double center[ 3 ] = { 0, 0, 0 };
      typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinate;
      const GlobalCoordinate barycenter = entity.geometry().center();
      for( int i = 0; i < GlobalCoordinate::dimension; ++i ) 
        center[ i ] = barycenter[ i ];
      db.vertexUpdate( typename LoadBalancerType::GraphVertex( index, weight, center ) );
#else // new dune-alugrid version
      db.vertexUpdate( typename LoadBalancerType::GraphVertex( index, weight ) );
#endif
      ++graphSize_; 
    }
    
    // set weight for faces (to be revised)
    updateFaces( entity, ibegin, iend, weight, db );   
  }

  template <class IntersectionIteratorType>
  void updateFaces(const EntityType& en,
                   IntersectionIteratorType nit,
                   const IntersectionIteratorType endit,
                   const int weight,
                   DataBaseType & db) 
  {
    for( ; nit != endit; ++nit )
    {
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      const IntersectionType& intersection = *nit;
      if(intersection.neighbor())
      {
        EntityPointerType ep = intersection.outside();
        const EntityType& nb = *ep;

        if( nb.partitionType() == InteriorEntity )
        {
          const int eid = getIndex( en );
          const int nid = getIndex( nb );
#if HAVE_DUNE_ALUGRID 
          // the newest ALU version only needs the edges to be inserted only once
          if( eid < nid ) 
#endif
          // the older version works with double insertion 
          // insert edges twice, with both orientations 
          // the ALUGrid partitioner expects it this way 
          {
            db.edgeUpdate ( typename LoadBalancerType :: GraphEdge ( eid, nid, weight
#if HAVE_DUNE_ALUGRID 
                  ,-1, -1
#endif
                  ) );
          }
        }
      }
    }
  }

public:
  /** \brief 
      \param method    partitioning method, available are:
                       - kway      = METIS_PartGraphKway 
                       - recursive = METIS_PartGraphRecursive (default)
                       - sfc       = space filling curve (only in dune-alugrid)
  */                     
  bool serialPartition( const Method method = recursive ) 
  {
    if( pSize_ > 1 ) 
    {
      // if the graph size is smaller then the number of partitions 
      // the distribution is easy to compute 
      if( graphSize_ <= pSize_ ) 
      {
        partition_.resize( graphSize_ );
        for( int i=0; i<graphSize_; ++ i ) 
          partition_[ i ] = i;
      }
      else 
      {
        // || HAVE_METIS 
        if( method == recursive ) 
          partition_ = db_.repartition( mpAccess_, DataBaseType :: METIS_PartGraphRecursive, pSize_ );
        else if( method == kway ) 
          partition_ = db_.repartition( mpAccess_, DataBaseType :: METIS_PartGraphKway, pSize_ );
        else if( method == sfc ) 
        {
#if HAVE_DUNE_ALUGRID
          partition_ = db_.repartition( mpAccess_, DataBaseType :: ALUGRID_SpaceFillingCurve, pSize_ );
#else   
          DUNE_THROW(InvalidStateException,"ThreadPartitioner::serialPartition: dune-alugrid not found, therefore not sfc partitioning");
#endif
        }
        else 
          DUNE_THROW(InvalidStateException,"ThreadPartitioner::serialPartition: wrong method");
        assert( int(partition_.size()) >= graphSize_ );
      }

      /*
      assert( partition_.size() > 0 );
      std::vector< int > counter( pSize_ , 0 );
      for( size_t i =0; i<partition_.size(); ++i) 
      {
        std::cout << "part[" << i << "] = " << partition_[ i ]  << endl;
        ++counter[  partition_[ i ]  ];
      }
      */
      return partition_.size() > 0;
    }
    else 
    {
      partition_.resize( indexSet_.size( 0 ) );
      for( size_t i =0; i<partition_.size(); ++i )
        partition_[ i ] = 0;
      return false ;
    }
  }

  int getDestination( const int idx ) const 
  {
    return db_.getDestination ( idx ) ;
  }

  std::set < int, std::less < int > > scan() const 
  {
    return db_.scan();
  }

  int getRank( const EntityType& entity ) const 
  {
    //std::cout << "partSize = " << partition_.size()  << " idx = " << getIndex( entity )  << std::endl;
    assert( (int) partition_.size() > getIndex( entity ) );
    return partition_[ getIndex( entity ) ];
  }

  bool validEntity( const EntityType& entity, const int rank ) const 
  {
    return getRank( entity ) == rank;
  }

};

} // end namespace Dune 
#else 
#warning "Parallel ALUGrid not available"
#endif // ALU3DGRID_PARALLEL
#endif // ifndef DUNE_FEM_THREADPARTITIONER_HH

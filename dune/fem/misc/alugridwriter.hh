#ifndef DUNE_ALUGRIDWRITER_HH
#define DUNE_ALUGRIDWRITER_HH

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/alugrid/3d/topology.hh>

namespace Dune {

template <class GridPartType>
class GlobalConsecutiveIndexSet
{
  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType  :: Traits :: LocalIdSet IdSetType;
  typedef typename IdSetType :: IdType IdType;

  enum { dimension = GridType :: dimension };

  const GridPartType& gridPart_;
  const GridType& grid_;
  const IdSetType& idSet_;

  typedef typename GridType :: template Codim< 0 >  :: Geometry :: GlobalCoordinate  CoordinateType;

  typedef std::pair< int , CoordinateType > VertexType ;
  typedef std::map< IdType, VertexType > IndexMapType;
  typedef typename IndexMapType :: iterator IndexMapIteratorType;
  mutable IndexMapType indices_;

  class DataHandle : 
    public CommDataHandleIF< 
     DataHandle, int > 
  {
    const IdSetType& idSet_;
    const int myRank_;
    IndexMapType& indices_;
  public:
    explicit DataHandle( const IdSetType& idSet, 
                         const int rank,
                         IndexMapType& indices ) 
      : idSet_( idSet ),
        myRank_( rank ),
        indices_( indices )
    {}


    //! returns true if combination is contained 
    bool contains ( int dim, int codim ) const
    {
      return dim == codim;
    }

    //! return whether we have a fixed size 
    bool fixedsize ( int dim, int codim ) const
    {
      return true;
    }

    //! return local dof size to be communicated 
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      return 2;
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void gather ( MessageBuffer &buffer,
                  const Entity &entity ) const
    {
      assert( (int) Entity :: codimension == (int) dimension );
      buffer.write( myRank_ );
      IndexMapIteratorType it = indices_.find( idSet_.id( entity ) );
      int index = ( it == indices_.end() ) ? -1 : (*it).second.first ;
      buffer.write( index );
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void scatter ( MessageBuffer &buffer,
                   const Entity &entity,
                   const size_t dataSize )
    {
      assert( (int) Entity :: codimension == (int) dimension );
      int rank;
      buffer.read( rank );
      int index ;
      buffer.read( index );
      if( myRank_ > rank && index > -1 ) 
      {
        const IdType id = idSet_.id( entity );
        IndexMapIteratorType it = indices_.find( id );
        if( it == indices_.end() ) 
        {
          VertexType vx( index, entity.geometry().corner(0));
          indices_[ id ] = vx;
        }
      }
    }
  };
public:
  explicit GlobalConsecutiveIndexSet( const GridPartType& gridPart ) 
   : gridPart_( gridPart ),
     grid_( gridPart.grid() ),
     idSet_( grid_.localIdSet() )
  {
    const int pSize = grid_.comm().size();
    const int myRank = grid_.comm().rank();
    int index = 0;
    for( int p = 0; p<pSize; ++p )
    {
      if( p == myRank ) 
      {
        typedef typename GridPartType :: template Codim< dimension > :: IteratorType
          IteratorType;
        const IteratorType endit = gridPart_.template end< dimension > ();
        for( IteratorType it = gridPart_.template begin< dimension > ();
             it != endit; ++it )
        {
          const IdType id = idSet_.id( *it );
          if( indices_.find( id ) == indices_.end() )
          {
            VertexType vx( index, it->geometry().corner(0) );
            indices_[ id ] = vx ;
            ++index;
          }
        }
      }

      std::cout << "P["<<myRank<< "] inddex = " << index << std::endl;
      // send current index number 
      grid_.comm().broadcast(&index, 1, p );
      std::cout << "P["<<myRank<< "] inddex = " << index << std::endl;

      {
        DataHandle dataHandle( idSet_, myRank, indices_ );
        gridPart_.communicate( dataHandle, 
                               InteriorBorder_InteriorBorder_Interface, 
                               ForwardCommunication );
      }
    }
  }

  void writeCoordinates ( std::ostream& out ) const 
  {
    out << indices_.size() << std::endl; 
    out.precision( 16 );
    out << std::scientific;
    IndexMapIteratorType end = indices_.end();
    for( IndexMapIteratorType it = indices_.begin(); it != end; ++it)
    {
      out << (*it).second.second << std::endl;
    }
  }

  void writeIndices ( std::ostream& out ) const 
  {
    IndexMapIteratorType end = indices_.end();
    for( IndexMapIteratorType it = indices_.begin(); it != end; ++it)
    {
      out << (*it).second.first << "  -1"  << std::endl;
    }
  }

  template <class EntityType>
  int index ( const EntityType& entity ) const 
  {
    assert( (int) EntityType :: codimension == (int) dimension );
    IndexMapIteratorType it = indices_.find( idSet_.id( entity ) );
    assert( it != indices_.end() );
    return (*it).second.first; 
  }

  int size() const 
  {
    return indices_.size();
  }
};

template <class GridPartType> 
class ALUGridWriter
{
  const GridPartType& gridPart_;

  typedef GlobalConsecutiveIndexSet < GridPartType > IndexSetType;

  typedef typename GridPartType :: GridType  GridType;

  IndexSetType indexSet_;

  enum { dimension = GridType :: dimension };

  typedef typename GridType :: template Codim< 0 > :: Entity  Entity;
protected:  
  ALUGridWriter( const GridPartType& gridPart ) 
    : gridPart_( gridPart ),
      indexSet_( gridPart_ )
  {
  }

  void write(const std::string& filename) const 
  {
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;
    const IteratorType endit = gridPart_.template end< 0 > ();
    IteratorType it = gridPart_.template begin< 0 > ();
    if( it == endit ) return;

    bool hexahedra = it->type().isHexahedron();
    if( ! hexahedra && ! it->type().isTetrahedron() ) 
    {
      DUNE_THROW(InvalidStateException,"Wrong geometry type");
    }

    int noElements = 0;
    for(; it != endit; ++it )
    {
      ++noElements;
    }

    std::stringstream filestr; 
    filestr << filename << "." << gridPart_.grid().comm().rank();

    std::ofstream file ( filestr.str().c_str() );
    if( ! file ) 
    {
      std::cerr << "ERROR: couldn't open file " << filestr.str() << std::endl;
      assert( false );
      abort();
    }

    const char* header = ( hexahedra ) ? "!Hexahedra" : "!Tetrahedra";
    // write header 
    file << header;
    file << "  ( noVertices = " << indexSet_.size();
    file << " | noElements = " << noElements << " )" << std :: endl;

    // write vertex coordinates  
    indexSet_.writeCoordinates( file );

    file << noElements << std::endl;
    if( hexahedra ) 
    {
      typedef ElementTopologyMapping< hexa > ElementTopo;
      writeElements< ElementTopo >( file );
      writeBoundaries< ElementTopo >( file );
    }
    else 
    {
      typedef ElementTopologyMapping< tetra > ElementTopo;
      writeElements< ElementTopo >( file );
      writeBoundaries< ElementTopo >( file );
    }

    // write global numbers of indices 
    indexSet_.writeIndices( file );
  }

  int getIndex( const Entity& entity, const int i ) const 
  {
    typedef typename GridType :: template Codim< dimension > :: EntityPointer  EntityPointer;
    EntityPointer vx = entity.template subEntity< dimension > ( i );
    return indexSet_.index( *vx );
  }
  
  template <class ElementTopo> 
  void writeElements( std::ostream& out ) const 
  {
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
    const IteratorType endit = gridPart_.template end< 0 > ();
    for(IteratorType it = gridPart_.template begin< 0 > ();
        it != endit; ++it )
    {
      const Entity& entity = *it ;
      const int nVx = entity.template count< dimension > ( );

      out << getIndex( entity, ElementTopo::dune2aluVertex( 0 ) );
      for(int i=1; i<nVx; ++i) 
      {
        out << "  " << getIndex( entity, ElementTopo::dune2aluVertex( i ));
      }
      out << std::endl;
    }
  }

  template <class ElementTopo> 
  void writeBoundaries( std::ostream& out ) const 
  {
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType IteratorType;
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType :: template Codim< 0 > :: Entity  Entity;
    const IteratorType endit = gridPart_.template end< 0 > ();
    int bndFaces = 0;
    for(IteratorType it = gridPart_.template begin< 0 > ();
        it != endit; ++it )
    {
      const Entity& entity = *it ;
      const IntersectionIteratorType endnit = gridPart_.iend( entity );
      for( IntersectionIteratorType nit = gridPart_.ibegin( entity );
           nit != endnit ; ++nit ) 
      {
        typedef typename IntersectionIteratorType :: Intersection Intersection;
        const Intersection& inter = *nit;
        if( inter.boundary() )
          ++bndFaces;
        else if( inter.neighbor() && 
                 inter.outside()->partitionType() != InteriorEntity ) 
         ++bndFaces;
      }
    }

    out << bndFaces << std::endl;

    for(IteratorType it = gridPart_.template begin< 0 > ();
        it != endit; ++it )
    {
      const Entity& entity = *it ;
      typedef typename GridType :: ctype coordType;
      const GenericReferenceElement< coordType, dimension > &refElem
         = GenericReferenceElements< coordType, dimension >::general( entity.type() );

      const IntersectionIteratorType endnit = gridPart_.iend( entity );
      for( IntersectionIteratorType nit = gridPart_.ibegin( entity );
           nit != endnit ; ++nit ) 
      {
        typedef typename IntersectionIteratorType :: Intersection Intersection;
        const Intersection& inter = *nit;
        int bndId = 0;
        if( inter.boundary() ) 
        {
          bndId = inter.boundaryId();
        }
        else if( inter.neighbor() && 
                 inter.outside()->partitionType() != InteriorEntity )
        {
          bndId = 111; 
        }

        if( bndId != 0 )
        {
          out << -bndId << "  ";
          const int duneFace  = inter.indexInInside();
          const int vxNr    = refElem.size( duneFace, 1, dimension );
          out << vxNr;

          std::vector< int > vertices( vxNr );
          const int aluFace = ElementTopo :: generic2aluFace( duneFace );
          for( int i=0; i<vxNr; ++i) 
          {
            const int j = ElementTopo :: faceVertex( aluFace, i );
            const int k = ElementTopo :: alu2genericVertex( j );
            typedef typename GridType :: template Codim< dimension > :: EntityPointer  EntityPointer;
            EntityPointer vx = entity.template subEntity< dimension > ( k );
            out << "  " << indexSet_.index( *vx );
          }
          out << std::endl;
        }
      }
    }
  }

public:  
  static void dumpMacroGrid(const GridPartType& gridPart, const std::string& filename)
  {
    ALUGridWriter< GridPartType > writer ( gridPart );
    writer.write( filename );
  }
};

}

#endif

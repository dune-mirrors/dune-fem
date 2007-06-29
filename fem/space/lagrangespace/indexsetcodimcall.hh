#ifndef DUNE_LAGRANGESPACE_INDEXSETCODIMCALL_HH
#define DUNE_LAGRANGESPACE_INDEXSETCODIMCALL_HH

#include "indexsetcodimcall.hh"

namespace Dune
{

  template< class GridPartImp >
  class IndexSetCodimCallInterface
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;

  private:
    typedef IndexSetCodimCallInterface< GridPartType > ThisType;
      
  public:
    virtual ~IndexSetCodimCallInterface ()
    {
    }
    
    virtual int subIndex ( const IndexSetType &indexSet,
                           const Entity0Type &entity,
                           int i ) const
    {
      DUNE_THROW( NotImplemented, "IndexSetCodimCallInterface :: "
                                  "subIndex is an abstract function." );
      return -1;
    }
  };
   

  
  template< class GridPartImp, unsigned int codim >
  class IndexSetCodimCall : public IndexSetCodimCallInterface< GridPartImp >
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;

  private:
    typedef IndexSetCodimCall< GridPartType, codim > ThisType;
    typedef IndexSetCodimCallInterface< GridPartType > BaseType;

  public:
    virtual int subIndex ( const IndexSetType& indexSet,
                           const Entity0Type &entity,
                           int i ) const
    {
      return indexSet.template subIndex< codim >( entity, i );
    }
  };



  template< class GridPartImp,
            unsigned int codim = GridPartImp :: GridType :: dimension >
  class IndexSetCodimCallFactory
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: GridType GridType;

    enum { dimension = GridType :: dimension };
    enum { codimension = codim };
    
    typedef IndexSetCodimCallInterface< GridPartType > IndexSetCodimCallType;

  public:
    static IndexSetCodimCallType* getCall( unsigned int i )
    {
      if( codimension == i )
        return new IndexSetCodimCall< GridPartType, codimension >();
      else
        return IndexSetCodimCallFactory< GridPartType, codimension - 1 >
               :: getCall( i );
    }

    static void getAllCalls( IndexSetCodimCallType *calls[ dimension+1 ] )
    {
      calls[ codimension ]
        = new IndexSetCodimCall< GridPartType, codimension >();
      IndexSetCodimCallFactory< GridPartType, codimension - 1 >
        :: getAllCalls( calls );
    }
  };


  
  template< class GridPartImp >
  class IndexSetCodimCallFactory< GridPartImp, 0 >
  {
  public:
    typedef GridPartImp GridPartType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: GridType GridType;
    
    enum { dimension = GridType :: dimension };
    enum { codimension = 0 };

    typedef IndexSetCodimCallInterface< GridPartType > IndexSetCodimCallType;

  public:
    static IndexSetCodimCallType* getCall( unsigned int i )
    {
      assert( i == 0 );
      return new IndexSetCodimCall< GridPartType, 0 >();
    }
    
    static void getAllCalls( IndexSetCodimCallType *calls[ dimension+1 ] )
    {
      calls[ codimension ]
        = new IndexSetCodimCall< GridPartType, codimension >();
    }
  };
 
}

#endif

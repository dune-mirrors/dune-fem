#ifndef DUNE_LAGRANGESPACE_INDEXSETCODIMCALL_HH
#define DUNE_LAGRANGESPACE_INDEXSETCODIMCALL_HH

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
                           int i ) const = 0;
  };
   

  
  template< class GridPartImp, unsigned int codim >
  class IndexSetCodimCall : public IndexSetCodimCallInterface< GridPartImp >
  {
  public:
    typedef GridPartImp GridPartType;

  private:
    typedef IndexSetCodimCall< GridPartType, codim > ThisType;

  public:
    typedef IndexSetCodimCallInterface< GridPartType > BaseType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;

  public:
    virtual int subIndex ( const IndexSetType& indexSet,
                           const Entity0Type &entity,
                           int i ) const
    {
      return indexSet.template subIndex< codim >( entity, i );
    }
  };
 
}

#endif

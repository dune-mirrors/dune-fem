#ifndef DUNE_FEM_CODIMMAP_HH
#define DUNE_FEM_CODIMMAP_HH

namespace Dune
{

  template< unsigned int nCodims, template< unsigned int > class CodimObjectImp >
  class CodimMap
  {
  public:
    enum { numCodims = nCodims };

  private:
    typedef CodimMap< numCodims, CodimObjectImp > ThisType;

    CompileTimeChecker< (numCodims > 0) > __numCodims_Must_Be_Positive__;

  public:
    typedef typename CodimObjectImp< 0 > :: BaseType CodimObjectBaseType;

  protected:
    CodimObjectBaseType *codimObjects_[ numCodims ];

  public:
    inline CodimMap ()
    {
      Factory< ThisType, numCodims - 1 > :: getObjects( *this );
    }
    
    inline ~CodimMap ()
    {
      for( unsigned int codim = 0; codim < numCodims; ++codim )
        delete codimObjects_[ codim ];
    }
    
    inline const CodimObjectBaseType &operator[] ( unsigned int codim ) const
    {
      assert( codim < numCodims );
      return *(codimObjects_[ codim ]);
    }

    inline CodimObjectBaseType &operator[] ( unsigned int codim )
    {
      assert( codim < numCodims );
      return *(codimObjects_[ codim ]);
    }

  protected:
    template< class CodimMapType,
              unsigned int codim >
    struct Factory
    {
      static void getObjects ( CodimMapType &map )
      {
        map.codimObjects_[ codim ] = new CodimObjectImp< codim >();
        assert( map.codimObjects_[ codim ] != 0 );
        Factory< CodimMapType, codim - 1 > :: getObjects( map );
      }
    };

    template< class CodimMapType >
    struct Factory< CodimMapType, 0 >
    {
      enum { codim = 0 };

      static void getObjects ( CodimMapType &map )
      {
        map.codimObjects_[ codim ] = new CodimObjectImp< codim >();
        assert( map.codimObjects_[ codim ] != 0 );
      }
    };

    template< class, unsigned int >
    friend struct Factory;
  };
  
}

#endif

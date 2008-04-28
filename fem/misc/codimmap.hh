#ifndef DUNE_FEM_CODIMMAP_HH
#define DUNE_FEM_CODIMMAP_HH

#include <dune/common/misc.hh>
#include <dune/fem/misc/metaprogramming.hh>

namespace Dune
{

  /** \class CodimMap
   *  \brief Simple array-like map for templates parameterized over codimensions
   *
   *  In some cases we need to call functions that depend on the codimension
   *  as template parameter. Now, if the codimension is actually a variable,
   *  the fastest way of accessing the functions is through an array mapping
   *  the codimension of the function. Building this array has to be done by
   *  template metaprogramming and hence is tedious and error prone. So, this
   *  is exactly what CodimMap does implicitly. Even the destruction of these
   *  objects is done implicitly be the destructor; you only need to destroy
   *  the CodimMap.
   *
   *  To use CodimMap you provide a class template (CodimObjectImp) depending
   *  only on the codimension as template parameter. All instances of this
   *  template have to have a common, virtual base class and a constructor
   *  requiring no arguments. CodimMap will then create an array of numCodims
   *  instances of CodimObjectImp. To the outside world it presents these
   *  instances as an array of objects of the common base class.
   *
   *  For example, a CodimMap could be created as follows:
   *  \code
   *    CodimMap< dimension+1, MyCodimObjectImp > map;
   *  \endcode
   *  Then, accessing the objects is similar to accessing an array:
   *  \code
   *    MyCodimBaseType &object = map[ codim ];
   *  \endcode
   *
   *  \param  nCodims         number of codimensions to create (usually
   *                          dimension+1)
   *  \param  CodimObjectImp  class template (parameterized by the codimension)
   *                          of which nCodims instances shall be created
   */
  template< unsigned int nCodims, template< unsigned int > class CodimObjectImp >
  class CodimMap
  {
  public:
    //! number of codimensions (= size of the array)
    enum { numCodims = nCodims };

  private:
    typedef CodimMap< numCodims, CodimObjectImp > ThisType;

    CompileTimeChecker< (numCodims > 0) > __numCodims_Must_Be_Positive__;

  public:
    //! common base type of all the instances of CodimObjectImp
    typedef typename CodimObjectImp< 0 > :: BaseType CodimObjectBaseType;

  protected:
    CodimObjectBaseType *codimObjects_[ numCodims ];

  private:
    template< unsigned int codim >
    struct MapConstructor
    {
      typedef ThisType &ArgumentType;

      inline static void apply ( ArgumentType map )
      {
        map.codimObjects_[ codim ] = new CodimObjectImp< codim >();
        assert( map.codimObjects_[ codim ] != 0 );
        // std :: cout << "Creating object " << codim << std :: endl;
      }
    };

  public:
    //! constructor building the CodimMap
    inline CodimMap ()
    {
      Loop< MetaSequence, MapConstructor, numCodims - 1 > :: apply( *this );
      //Factory< ThisType, numCodims - 1 > :: getObjects( *this );
    }
    
    //! destructor freeing all the instances of CodimObjectImp
    inline ~CodimMap ()
    {
      for( unsigned int codim = 0; codim < numCodims; ++codim )
        delete codimObjects_[ codim ];
    }
    
    /** \brief access the instance for one codimension
     *
     *  \param[in]  codim  codimension, the instance is to be obtained for
     *
     *  \returns instance of CodimObjectImp with the given codimension
     */
    inline const CodimObjectBaseType &operator[] ( unsigned int codim ) const
    {
      assert( codim < numCodims );
      return *(codimObjects_[ codim ]);
    }

    /** \brief access the instance for one codimension
     *
     *  \param[in]  codim  codimension, the instance is to be obtained for
     *
     *  \returns instance of CodimObjectImp with the given codimension
     */
    inline CodimObjectBaseType &operator[] ( unsigned int codim )
    {
      assert( codim < numCodims );
      return *(codimObjects_[ codim ]);
    }
 
#if 0
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
#endif
  };
  
}

#endif

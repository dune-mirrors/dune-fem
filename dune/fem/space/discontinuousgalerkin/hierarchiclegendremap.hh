#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HIERARCHICLEGENDREMAP_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HIERARCHICLEGENDREMAP_HH

namespace Dune
{
  namespace Fem
  {
    template<int polOrder, int dimension>
    class HierarchicLegendreMap
    {
      enum{numBaseFct=StaticPower<polOrder+1,dimension>::power};
      
      protected:
      
      typedef  array<int, numBaseFct> BaseNumberMapType;
    
      
      template<int pOrd, int d>
      struct DimLoop 
      {
        template < class A >
        static void loop( const int term, 
                          A& polOrders,
                          BaseNumberMapType& baseMap, 
                          int& oldBaseNum,
                          int& newBaseNum )
        {
          // go through all p because of old base counter
          for( int p=0; p<=pOrd; ++p) 
          {
            polOrders[ d ] = p ;
            DimLoop< pOrd, d-1> :: loop( term, polOrders, baseMap, oldBaseNum, newBaseNum );
          }
        }
      };

      template<int pOrd>
      struct DimLoop<pOrd, 0>
      {
        template< class A >
        static bool checkEntries( const A& array, const int p ) 
        {
          // at least one array entry has to be of order p
          // in addition no entry must be larger than order p 
          bool found = false ;
          for(typename A::size_type i=0; i<array.size(); ++i) 
          {
            if( array[ i ] == p ) found = true ;
            if( array[ i ] > p ) return false ;
          }
          return found ;
        }

        template< class A > 
        static void loop( const int term, 
                          A& polOrders,
                          BaseNumberMapType& baseMap, 
                          int& oldBaseNum,
                          int& newBaseNum )
        {
          // go through all p because of old base counter
          for( int p=0; p<=pOrd; ++p, ++oldBaseNum) 
          {
            polOrders[ 0 ] = p ;

            // at least one of the entries has to be of order p
            // non of the entries must be larger than p 
            if( checkEntries( polOrders, term ) )
            {
              assert( newBaseNum < int(numBaseFct) );
              baseMap[ newBaseNum ] = oldBaseNum;
              ++newBaseNum;
            }
          }
        }
      };
      
      // map for renumbering of basis functions to have a hierarchical basis 
      BaseNumberMapType baseFunctionMap_;
  

      public:
      //default constructor
      HierarchicLegendreMap()
      {
        baseFunctionMap_.fill(-1);
      
        // small array to store the curently used pol ords 
        array<int, dimension> polStorage ;
        polStorage.fill( -1 );

        int newBaseNum = 0;
        // check for all terms the number of base functions 
        for(int term=0; term <= polOrder; ++term )
         {
           // initialize oldBaseNum counter new for each term 
           int oldBaseNum = 0;
           // construct mapping for this term 
           DimLoop<polOrder, dimension-1> :: 
             loop( term, polStorage, baseFunctionMap_, oldBaseNum, newBaseNum );
         }
       
#ifndef NDEBUG 
        //std::cout << numBaseFct << std::endl;
        for( int i=0; i<numBaseFct; ++i ) 
        {
          std::cout << "base " << i << " is " << baseFunctionMap_[ i ] << std::endl;
          assert( baseFunctionMap_[ i ] >= 0 );
        }
#endif
 
      }
 

      const int operator[](int i) const {return baseFunctionMap_[i];}

    };
    

  
  }// namespace Fem

}//namespace Dune



#endif //DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HIERARCHICLEGENDREMAP_HH



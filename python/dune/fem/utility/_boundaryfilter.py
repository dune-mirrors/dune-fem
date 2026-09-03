# returns a boolean array contains whether or not an element belongs to a
# certain boundary
def boundaryFilter( gridView ):
    """
    Parameter:
        gridView     a grid view to create filter for

    Returns:
        returns a function that creates a boundary filter given a list of
        boundary ids.
    """
    import io, sys
    from numpy import array
    from dune.generator import algorithm

    code = """
    #include <dune/fem/misc/boundaryidprovider.hh>
    #include <dune/grid/common/rangegenerators.hh>

    template <class GridView>
    std::vector<int> boundaryId( const GridView& gv, const std::vector<int> & bndIds, const int domainId )
    {
      const auto& indexSet = gv.indexSet();
      std::vector< int > filter( indexSet.size(0), 0 );
      const bool allBoundaryIds = bndIds.empty();
      for( const auto& e : Dune::elements( gv, Dune::Partitions::interior ))
      {
        for( const auto& isec : Dune::intersections( gv, e ) )
        {
          if( isec.boundary() )
          {
            bool foundMatch = allBoundaryIds ;

            if( ! foundMatch )
            {
              const int myId = Dune::Fem::boundaryId( gv, isec );
              for( const auto& id : bndIds )
              {
                if( id == myId )
                {
                  foundMatch = true;
                  break;
                }
              }
            }

            if( foundMatch )
            {
              filter[ indexSet.index( e ) ] = domainId;
              break; // next element
            }
          }
        }
      }
      return filter;
    }
    """
    bndList = [int(0)]
    bndIdCheck = algorithm.load("boundaryId", io.StringIO(code), gridView, bndList, int(0) )

    def filterFct( boundaryIds = [], domainId = 1 ):
        """
        Parameter:
            boundaryIds  list containing one or more boundary ids. An empty list
                         refers to all boundary ids.
            domainId     label for the resulting domain (default is 1).

        Returns:
            returns a boolean array that contains which contains true of all
            elements with a corresponding boundary id from the given boundary ids.
        """
        assert isinstance(boundaryIds, (list,tuple))
        return bndIdCheck( gridView, boundaryIds, domainId )

    return filterFct

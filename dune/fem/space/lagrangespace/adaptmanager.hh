#ifndef DUNE_LAGRANGESPACE_ADAPTMANAGER_HH
#define DUNE_LAGRANGESPACE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/common/restrictprolonginterface.hh>

#include <dune/fem/space/lagrangespace/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/restrictprolong.hh>

namespace Dune
{

  /** @ingroup RestrictProlongImpl
   *
   *  \brief Restriction / prolongation operator for Lagrange discrete
   *         function spaces
   */
  template< class DF, class FS, class GP, int ord, template< class > class S >
  class RestrictProlongDefaultImplementation< DF, LagrangeDiscreteFunctionSpace< FS, GP, ord, S > >
  : public RestrictProlongInterfaceDefault< RestrictProlongTraits
      < RestrictProlongDefaultImplementation< DF, LagrangeDiscreteFunctionSpace< FS, GP, ord, S > >,
        typename GP::GridType::ctype > >
  {
    typedef RestrictProlongDefaultImplementation
      < DF, LagrangeDiscreteFunctionSpace< FS, GP, ord, S > >
      ThisType;
    typedef RestrictProlongInterfaceDefault< RestrictProlongTraits< ThisType, typename GP::GridType::ctype > >
      BaseType;

    typedef typename GP::GridType::ctype DomainFieldType;

  public:
    //! type of the discrete function
    typedef DF DiscreteFunctionType;

    //! type of the discrete function space
    typedef LagrangeDiscreteFunctionSpace< FS, GP, ord, S > DiscreteFunctionSpaceType;

  protected:
    using BaseType::entitiesAreCopies;

  public:
    //! type of the local functions
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

    //! type of the grid
    typedef typename DiscreteFunctionSpaceType::GridType Grid;

  public:
    //! constructor
    explicit
    RestrictProlongDefaultImplementation ( DiscreteFunctionType &discreteFunction )
    : discreteFunction_( discreteFunction )
    {}

    //! restrict data to the father
    template< class Entity >
    void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( !entitiesAreCopies( discreteFunction_.space().indexSet(), father, son ) )
      {
        LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
        LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

        localRestrictProlong_.restrictLocal( fatherFunction, sonFunction, initialize );
      }
    }

    //! prolong data to children
    template< class Entity >
    void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
    {
      if( !entitiesAreCopies( discreteFunction_.space().indexSet(), father, son ) )
      {
        LocalFunctionType fatherFunction = discreteFunction_.localFunction( father );
        LocalFunctionType sonFunction = discreteFunction_.localFunction( son );

        localRestrictProlong_.prolongLocal( fatherFunction, sonFunction );
      }
    }

    //! add discrete function to communicator 
    template< class Communicator >
    void addToList ( Communicator &comm )
    {
      // for Lagrange spaces this communication is not needed (since
      // data on ghosts is neglected 

      //  comm.addToList( discreteFunction_ );
    }

  private:
    DiscreteFunctionType &discreteFunction_;
    LagrangeLocalRestrictProlong< Grid, ord > localRestrictProlong_;
  };

}

#endif // #ifndef DUNE_LAGRANGESPACE_ADAPTMANAGER_HH

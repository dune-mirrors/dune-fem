#ifndef DUNE_FEMPY_PYYHON
#define DUNE_FEMPY_PYYHON

#include <dune/fempy/pybind11/pybind11.h>

#include <iostream>

#include <dune/common/std/utility.hh>

#include "pytuple.hh"

namespace PyDune
{
  // Always have void addToPython (...)
  // User defines as many functions of the type as he wants
  // template< class... Args >
  // void addToPython (pybind11::class_< T, Args... > &)
  // For this setting, return false if no standard constructor is to be added by default
  //
  // The user can also define a method
  // void addToPython( const char *) and in the python export section we
  // always call addToPython( modulename ) so that the user can add module
  // level stuff - return type is ignored for the module level call.
  bool addToPython(...) { return true; }

  // class taking a return type R of a method and provides information how
  // this type is exported to python - second argument is the possibly
  // required Container for keeping references alive. So the assumption is
  // that the type exposed to python is either
  // - identical to R
  // - Wrapper<R>
  // - Wrapper<R,Container>
  // When writting a wrapper class PyClass for a given dune class this clas
  // can be used for the return value R of a method on the dune class which
  // is to be exported on the wrapper class PyClass
  template < class R, class Container=void>
  struct ReturnValueConverter
  {
    // type used for exporting R to python
    typedef R type;
    // method used to do any required registration for the return type
    static void registerType(const char* name) {}
    // method used to convert the dune return to the C++ type exported to python
    // we assume that the C++ wrapper for R takes an instance of a
    // reference of r and some additional arguments (i.e. the container)
    template <class... ConstrArgs>
    static type convert(const R &r, ConstrArgs... args)
    { return r; }
  };

  // make sure that specialization work for references and const types
  template <class R, class Container>
  struct ReturnValueConverter<R&, Container> :
    public ReturnValueConverter<R, Container> {};
  template <class R, class Container>
  struct ReturnValueConverter<const R&, Container> :
    public ReturnValueConverter<R, Container> {};
  template <class R, class Container>
  struct ReturnValueConverter<const R, Container> :
    public ReturnValueConverter<R, Container> {};
  template <class R>
  struct ReturnValueConverter<R&, void> :
    public ReturnValueConverter<R> {};
  template <class R>
  struct ReturnValueConverter<const R&, void> :
    public ReturnValueConverter<R> {};
  template <class R>
  struct ReturnValueConverter<const R, void> :
    public ReturnValueConverter<R> {};

  // specialization of the ReturnValueConverter for tuple types
  template <class... TArgs, class Container>
  struct ReturnValueConverter<std::tuple<TArgs...>,Container>
  {
    template <class R>
    using Converter = ReturnValueConverter<Container, R>;
    typedef std::tuple<TArgs...> R;
    typedef std::tuple< typename Converter< TArgs >::type... > type;

    template <class... ConstrArgs>
    static type convert(const std::shared_ptr<Container> &c, const R &r, ConstrArgs... args)
    {
      return convert_impl( c, r, Dune::Std::index_sequence_for< TArgs... >(), args... );
    }
    static void registerType(const char* name)
    {
      TupleConverter<type>::create();
      registerType_(name, Dune::Std::index_sequence_for< TArgs... >() );
    }
    template< class... Args, class Method >
    static void def(pybind11::class_< Args... > &cls,
                   const char* name, Method method)
    {
      cls.def(name,method);
      registerType( (std::string(name)+"_sol").c_str() );
    }

    private:
    template< class... ConstrArgs,  std::size_t... i >
    static type convert_impl ( const std::shared_ptr<Container> &c, const R &r,
                               Dune::Std::index_sequence< i... >, ConstrArgs... args )
    {
      return std::make_tuple( Converter< typename std::tuple_element< i, R >::type >::
             convert( c, std::get< i >( r ), args... )... );
    }
    template< std::size_t... i >
    static void registerType_(const char *name, Dune::Std::index_sequence< i... > )
    {
      std::ignore = std::make_tuple( (Converter< typename std::tuple_element< i, R >::type >::registerType( (name+std::to_string( i )).c_str() ), i)... );
    }
  };
}
#endif // DUNE_FEMPY_PYYHON

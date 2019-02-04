# Master (will become release 2.7)

-  Added a new way to pass parameters to classes accepting an instance of
   the `ParameterReader`. The free standing function `parameterDict` in
   `dune/fem/io/parameter.hh` takes a prefix string and a variable number
   of parameters representing `key-value` pairs:
   ```
   int a;
   Dune::Fem::ParameterReader paraeters = Dune::Fem::paramDict("prefix",
           "key1",1.0, "key2","stringValue",
           "key3",[&a]() { return a; });
   ```
   The keys have to be of type `std::string` and the values can be either
   of type `std::string` or of a type which can be converted to a `std::string`
   using `std::to_string`. In addition `invocables` can be passed in
   returning one of the above. The result is a `ParameterReader` which
   returns the values for the given keys with the provided prefix, i.e.,
   `prefix.key1". If a querry for a key is made which is not one of the
   arguments then the querry is passed on to the global parameter
   container.

-  Overhaul of the way parameters could be passed into inverse operators:
   linear and nonlinear inverse operators should all now provide a type
   `SolverParameterType` which contains all methods to access user
   defineable parameters. This class should be constructable (implictely)
   from an instance of `ParameterReader`. The inverse operators
   should also now only have one constructor
   which takes an instance of `SolverParameterType`:
   ```
    explicit NewtonInverseOperator ( const ParameterType &parameter = Parameter::container() )
   ```
   All linear inverse operators should derive from the
   `InverseOperatorInterface` class and their `SolverParameterType` should
   derive from `SolverParametrs`.

   __Note__: using the `LocalParameter` approach to derive from the
   `SolverParameterType` to overload those methods will not always work as expected
   and should therefore not be used. To overload some given parameter the
   `parameterDict` provides similar functionality.

-  `reductionTol` and `absoluteTol` tolerance as solver parameters are
   deprecated. Use `tolerance` instead together with the `errormeassure`
   parameter.

-  `discreteFunction.communicate` and `linearOperator.communicate` are
   (nearly) deprecated, i.e., they shouldn't be called any more except in
   special circumtances.  Communication is taken care of by the local contributions.

-  The linear operators have a new method finalize which should be called
   after all changes to the matrix have been made. Although changing the
   matrix after this step could still be possible is some cases,
   no guarantee is given. To avoid issues with a premature finalization,
   two new methods `disableFinalize()` and `enableFinalize()` can be used -
   after `disableFinalize` is envoked any call to `finalize` on the linear
   operator will be ignored until `enableFinalize` is called.

    See dune-fem/dune-fem!290

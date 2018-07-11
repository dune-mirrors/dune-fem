# Release 2.6

## New discrete function spaces
- `TupleDiscreteFunctionSpace`
- spaces from `dune-localfunction` e.g. `RaviartThomas`
- p-adaptive discontinuous galerkin space

## Discrete functions and linear operator
- introduced a new concept of `bindable` local objects which are to
  replace the old `localFunction` concept

## Solvers
- unified some of the interfaces and replaced all the old `dune-fem`
  native solvers by a new `KrylovInverseOperator` class with a dynamic
  choice betweeen `cg/gmres/bicgstab` solvers

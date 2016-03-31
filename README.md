DUNE-FEM
========

[DUNE-FEM][0] is a [Distributed and Unified Numerics Environment][1]
module which defines interfaces for implementing discretization methods like Finite Element Methods (FEM)
and Finite Volume Methods (FV) and Discontinuous Galerkin Methods (DG).

If you need help, please ask on our [mailinglist][5]. Bugs can also be submitted
to the DUNE-FEM [bugtracker][6] instead.

Dependencies
------------

DUNE-FEM requires GCC (4.9+) or clang (3.5+) and depends on the following DUNE modules:

* [dune-common][10]

* [dune-geometry][11]

* [dune-grid][12]

The following DUNE modules are suggested:

* [dune-istl][13]

* [dune-localfunctions][14]

* [dune-alugrid][8]

* [dune-spgrid][9]

The following software is optional:

* [PETSc][3]

* [Eigen][4]

License
-------

The DUNE-FEM library, headers and test programs are free open-source software,
licensed under version 2 or later of the GNU General Public License.

See the file [LICENSE][7] for full copying permissions.

Installation
------------

For a full explanation of the DUNE installation process please read
the [installation notes][2].

 [0]: http://www.dune-project.org/fem/
 [1]: http://www.dune-project.org
 [2]: http://www.dune-project.org/doc/installation-notes.html
 [3]: http://www.mcs.anl.gov/petsc/
 [4]: http://eigen.tuxfamily.org
 [5]: http://lists.dune-project.org/mailman/listinfo/dune-fem
 [6]: http://gitlab.dune-project.org/dune-fem/dune-fem/issues
 [7]: LICENSE
 [8]: http://gitlab.dune-project.org/extensions/dune-alugrid
 [9]: http://gitlab.dune-project.org/extensions/dune-spgrid
 [10]: http://gitlab.dune-project.org/core/dune-common
 [11]: http://gitlab.dune-project.org/core/dune-geometry
 [12]: http://gitlab.dune-project.org/core/dune-grid
 [13]: http://gitlab.dune-project.org/core/dune-istl
 [14]: http://gitlab.dune-project.org/core/dune-localfunctions
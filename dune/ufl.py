from __future__ import absolute_import
import ufl

class Space(ufl.VectorElement):
    def __init__(self, dimDomain, dimRange):
        if isinstance(dimDomain, tuple):
            if len(dimDomain) != 2:
                raise Exception('dimDomain tuple must contain exactly two elements.')
            dimWorld = int(dimDomain[1])
            dimDomain = dimDomain[0]
        else:
            dimWorld = int(dimDomain)
        if dimDomain == 1:
            cell = ufl.Cell("interval", dimWorld)
        elif dimDomain == 2:
            cell = ufl.Cell("triangle", dimWorld)
        elif dimDomain == 3:
            cell = ufl.Cell("tetrahedron", dimWorld)
        else:
            raise NotImplementedError('dune.ufl.Space not implemented for dimension' + str(dimDomain) + '.')
        ufl.VectorElement.__init__(self, "Lagrange", cell, 1, dimRange)

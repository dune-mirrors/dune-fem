#!/usr/bin/python3

import pickle, numpy
import dune.generator
import ufl
from dune.grid import structuredGrid, gridFunction
from dune.fem.space import lagrange

class CheckPointer:
    def __init__(self,fileName,hgrid=None):
        if hgrid is None:
            self.file = open(fileName,'rb')
            self.hgrid = pickle.load(self.file)
        else:
            self.file = open(fileName,'wb')
            self.hgrid = hgrid
            pickle.dump(self.hgrid,self.file)
        self.items = []

    def hierarchicalGrid(self):
        return self.hgrid

    def add(self,item):
        self.items += [item]

    def backup(self):
        for i in self.items:
            if hasattr(i,"read"):
                pickle.dump(i.write(),self.file)
            else:
                pickle.dump(i,self.file)
    def restore(self):
        for i in self.items:
            if hasattr(i,"read"):
                i.read( pickle.load(self.file) )
            else:
                i = pickle.load( self.file ) # this doesn't work

def run(restore=False):
    if not restore:
        grid = structuredGrid([0,0],[1,1],[2,2])
        grid.hierarchicalGrid.globalRefine(2)
        checkPointer = CheckPointer("dumpA", grid.hierarchicalGrid)
    else:
        checkPointer = CheckPointer("dumpA")
        grid = checkPointer.hierarchicalGrid().leafView

    space = lagrange(grid)
    df = space.interpolate([0],name="test")
    checkPointer.add(df)

    spc = lagrange(grid, storage='istl')
    df_istl = spc.interpolate([0],name="test_istl")

    # test discrete function assignment
    df.assign( df_istl )
    df_istl.assign( df )

    if not restore:
        print("interpolating grid function")
        @gridFunction(grid)
        def gf(x): return x[0]*x[1]*(1-x[0])*(1-x[1])
        df.interpolate( gf )
    else:
        print("restoring discrete function")
        checkPointer.restore()

    df.plot()
    if not restore:
        checkPointer.backup()


if __name__ == "__main__":
    run(False)
    run(True)

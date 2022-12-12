import pickle
from dune.grid import structuredGrid
gridView = structuredGrid([0,0],[1,1],[2,2])
pickle.dump(gridView.hierarchicalGrid, open("structuredGrid.dbf",'wb'))

from dune.alugrid import aluConformGrid as leafGridView
gridView = leafGridView("sphere.dgf", dimgrid=2, dimworld=3)
pickle.dump(gridView.hierarchicalGrid, open("sphere.dbf",'wb'))

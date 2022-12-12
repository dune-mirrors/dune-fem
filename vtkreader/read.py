import dune.common.pickle
from dune.grid import Marker

with open("dump.dbf","rb") as f:
    dump = dune.common.pickle.load(f)
df = dump[2]
df.plot(level=3)
df.gridView.writeVTK("dump", pointdata=[df])
df.gridView.writeVTK("dump2", pointdata=[df], subsampling=2)

"""
# df.gridView.plot()
for i in range(4):
    df.gridView.hierarchicalGrid.mark(lambda e:
         Marker.refine if df.localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 4
         else Marker.keep)
    dune.fem.adapt([df])
    print(i,df.gridView.size(0))
df.gridView.plot()

import transform
l,_ = transform.error(df.gridView,[df])
l.plot()
"""

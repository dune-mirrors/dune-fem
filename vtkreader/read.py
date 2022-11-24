import dune.common.pickle
with open("dump.dbf","rb") as f:
    dump = dune.common.pickle.load(f)
df = dump[2]
df.gridView.writeVTK("dump", pointdata=[df])
df.gridView.writeVTK("dump2", pointdata=[df], subsampling=2)
import transform
l,_ = transform.error(df.gridView,[df])
l.plot()

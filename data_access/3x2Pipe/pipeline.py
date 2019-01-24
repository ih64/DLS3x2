import pandas as pd
import jpPipe

t = pd.read_hdf('../catalogs/MICE.hdf', '/data')
pipe = jpPipe.Pipe(t)
pipe.run()

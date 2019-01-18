import pandas as pd
import jpPipe

t = pd.read_hdf('../catalogs/buzzard_cutout.hdf', '/data')
pipe = jpPipe.Pipe(t)
pipe.run()

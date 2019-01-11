import pandas as pd
import jpPipe

#t = pd.read_hdf('/global/cscratch1/sd/ihasan/buzzard.hdf', '/data')
#t = pd.read_csv('../catalogs/F_all.tab', sep='\t')
t = pd.read_hdf('../catalogs/buzzard.hdf', '/data')
pipe = jpPipe.Pipe(t)
#pipe = jpPipe.Pipe()
pipe.run()

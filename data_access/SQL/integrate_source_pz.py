import pandas as pd
import numpy as np
import dask.dataframe as dd


scratch_path = '/global/cscratch1/sd/ihasan/source_pz.tab'
#scratch_path = '../catalogs/S1_testcase_pz.tab'
ddata = pd.read_csv(scratch_path, sep='\t')


def integrate(row, low, high):
    c_low = 'c'+str(low)
    c_high = 'c'+str(high - 1)
    total = np.trapz(row[c_low:c_high].values, x=np.arange(low,high,1))
    return total


w1 = ddata.apply(integrate, args=(40,60), axis=1)
w2 = ddata.apply(integrate, args=(60,80), axis=1)
w3 = ddata.apply(integrate, args=(80,100), axis=1)

ddata['w1'] = w1
ddata['w2'] = w2
ddata['w3'] = w3

#ddata.to_hdf('/global/cscratch1/sd/ihasan/output.hdf','/data')
ddata.to_hdf('source_pz.hdf','/data')


#w1 = ddata.apply(integrate, args=(40,60), axis=1, meta=('w1', 'f8'))
#w2 = ddata.apply(integrate, args=(60,80), axis=1, meta=('w2', 'f8'))
#w3 = ddata.apply(integrate, args=(80,100), axis=1, meta=('w3', 'f8'))
#w1 = w1.compute()
##w2 = w2.compute()
#w3 = w3.compute()

#ddata = ddata.reset_index().set_index('index')
#ddata = ddata.assign(w1=w1)
#ddata = ddata.assign(w2=w2)
#ddata = ddata.assign(w3=w3)
#ddata = ddata.compute()

#ddata.to_hdf('/global/cscratch1/sd/ihasan/source_pz.hdf', '/data') 
#ddata.to_parquet('/global/cscratch1/sd/ihasan/S1_testcase_pz.parquet')

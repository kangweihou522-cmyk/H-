"""Read saved result IDs as strings without changing numeric measurement columns."""
import pandas as pd


def read_frame(path, **kwargs):
    types=kwargs.get('dtype', {name:str for name in ('gene','id','symbol','entrez','term','source','target','sample','subject','group','module')})
    if kwargs.get('index_col') == 0:
        if isinstance(types,dict): types={0:str,**types}
    kwargs['dtype']=types
    return pd.read_csv(path, **kwargs)

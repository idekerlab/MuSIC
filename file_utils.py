import pandas as pd
import numpy as np
import pickle
import dill

def save_obj(obj, fname, method='pickle'):
    ''' Saving objects to designated filename in pickle format
    
    Args:
        obj: object that want to be saved
        fname: path to saved file
        method: {pickle, dill} specify package used for compressing
    '''
    with open(fname, 'wb') as f:
        if method == 'pickle':
            pickle.dump(obj, f)
        elif method == 'dill':
            dill.dump(obj, f)
        else:
            raise ValueError('Please select method from {pickle, dill}!')
    return

def load_obj(fname, method='pickle'):
    ''' Loading object that was saved in pickle format
    
    Args:
        fname: path to file
        method: {pickle, dill} specify package used for compressing
    '''
    with open(fname, 'rb') as f:
        if method == 'pickle':
            return pickle.load(f)
        elif method == 'dill':
            return dill.load(f)
        else:
            raise ValueError('Please select method from {pickle, dill}!')
        return
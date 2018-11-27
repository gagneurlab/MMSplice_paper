from collections import OrderedDict
import re
from pylab import find
import numpy as np
import pandas as pd
import warnings

def find_pos(longseq, shortseq, start=True):
    try:
        match = re.search(shortseq, longseq).span()
    except:
        # import pdb
        # pdb.set_trace()
        print(shortseq)
        if start:
            shortseq = shortseq[:5]
        else:
            shortseq = shortseq[-15:]
        match = re.search(shortseq, longseq).span()
        warnings.warn("Cannot match full length short sequence %s, used first 5 bases" % shortseq)
    if start:
        return match[0]
    else:
        return match[1]
        

def find_seq_diff_pos(wt_seq, mut_seq):
    """ Function to find the actual position of mutations based
    on the WT and MUT seq"""
    muts = np.zeros(len(wt_seq))
    for i in range(len(wt_seq)):
        try:
            muts[i] = (wt_seq[i]!=mut_seq[i])
        except IndexError:
            pass
    return find(muts)[0]

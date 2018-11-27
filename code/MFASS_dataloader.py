from kipoi.data import Dataset
from kipoi.metadata import GenomicRanges
import pandas as pd
import warnings

class MFASS_Exon(object):
    """ Class of MFASS exon
    Args: exon_cut_l: when extract exon feature, how many base pair to cut out at the begining of an exon
    exon_cut_r: when extract exon feature, how many base pair to cut out at the end of an exon 
    (cut out the part that is considered as acceptor site or donor site)
    acceptor_intron_cut: how many bp to cut out at the end of acceptor intron that consider as acceptor site
    donor_intron_cut: how many bp to cut out at the end of donor intron that consider as donor site
    intronl_len: intron length at the acceptor side
    intronr_len: intron length at the donor side
    acceptor_intron_len: what length in acceptor intron to consider for acceptor site model
    acceptor_exon_len: what length in acceptor exon to consider for acceptor site model
    donor_intron_len: what length in donor intron to consider for donor site model
    donor_exon_len: what length in donor exon to consider for donor site model
    """ 
        
    def __init__(self,
                 strand, 
                 intron1_len, 
                 exon_len, 
                 intron2_len, 
                 seq,
                 exon_cut_l=3,
                exon_cut_r=3,
                acceptor_intron_cut=20,
                donor_intron_cut=6,
                acceptor_intron_len=20,
                acceptor_exon_len=3,
                donor_exon_len=3,
                donor_intron_len=6):
        self.strand = strand
        self.exon_len = exon_len
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_intron_len = donor_intron_len
        self.donor_exon_len = donor_exon_len
        if strand in ['-', -1, '-1']:
            #seq = rev_seq(seq)
            self.intronl_len = intron2_len
            self.intronr_len = intron1_len
        else:
            self.intronl_len = intron1_len
            self.intronr_len = intron2_len            
        self.seq = seq
  
#     def get_acceptor_intron(self):
#         acceptor_intron = self.seq[:self.intronl_len-self.acceptor_intron_cut]
#         return acceptor_intron
    
#     def get_exon(self):
#         exon = self.seq[(self.intronl_len+self.exon_cut_l) : (self.intronl_len+self.exon_len-self.exon_cut_r)]
#         return exon
    
#     def get_donor_intron(self):
#         donor_intron = self.seq[-self.intronr_len+self.donor_intron_cut:]
#         return donor_intron
    
#     def get_acceptor_site(self):
#         acceptor = self.seq[self.intronl_len-self.acceptor_intron_len : self.intronl_len+self.acceptor_exon_len]
#         if acceptor[-self.acceptor_exon_len-2:-self.acceptor_exon_len] != "AG":
#             warnings.warn("None AG acceptor", UserWarning)
#         return acceptor
    
#     def get_donor_site(self):
#         donor = self.seq[-self.intronr_len-self.donor_exon_len : -self.intronr_len+self.donor_intron_len]
#         if donor[self.donor_exon_len:self.donor_exon_len+2] != "GT":
#             warnings.warn("None GT donor", UserWarning)
#         return donor
    
    def get_features(self):
        # get all features
        # acceptor_intron = self.get_acceptor_intron()
        # exon = self.get_exon()
        # donor_intron = self.get_donor_intron()
        # acceptor_site = self.get_acceptor_site()
        # donor_site = self.get_donor_site()
        # fts = {"acceptor_intron": acceptor_intron,
        #        "acceptor_site": acceptor_site,
        #        "exon": exon,
        #        "donor_site": donor_site,
        #        "donor_intron": donor_intron}
        # return fts
        x = {}
        x['seq'] = self.seq
        x['intronl_len'] = self.intronl_len
        x['intronr_len'] = self.intronr_len
        return x


class MFASS_Exon_Dataset(Dataset):
    def __init__(self, 
                 infile,
                 **kwargs):
        """ infile: process MFASS file 
        **kwargs: arguments feed to MFASS_Exon class
        """
        mfass = pd.read_table(infile)
        # alternative sequence
        self.exons_alt = mfass.apply(lambda x: MFASS_Exon(x['strand'], x['intron1_len'], 
                                 x['exon_len'], x['intron2_len'], 
                                 x['original_seq'], **kwargs), axis=1)
        # reference sequence
        self.exons_ref = mfass.apply(lambda x: MFASS_Exon(x['strand'], x['intron1_len'], 
                                 x['exon_len'], x['intron2_len'], 
                                 x['nat_seq'], **kwargs), axis=1)
        
    def __len__(self):
        return len(self.exons_alt)
    
    def __getitem__(self, idx):
        ref = self.exons_ref[idx]
        alt = self.exons_alt[idx]
        out = {}
        out['inputs'] = ref.get_features()
        out['inputs_mut'] = alt.get_features()
        return out
    
    
    
def rev_seq(seq):
    """Reverse the sequence and use the complement bases.
    This is for annotation on "-" strand.
    
    example: 
        rev_seq("atgc") 
        >>> 'gcat'
        rev_seq("akgc")
        >>> "k" is not valid base in:
        >>> dict_keys(['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'])
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return_seq = ""
    for base in seq:
        try:
            return_seq = complement[base] + return_seq
        except KeyError:
            print('"%s" is not valid base in:' %(base))
            print(complement.keys())
            raise
    return return_seq
    # from phylogeny_utilities.command_line_utils import *
import re, platform, time, os, sys, json
import numpy as np
import multiprocessing
# import matplotlib
# if platform.system()=='Linux':
#     matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import dendropy
# import phylogeny_utilities.common_vars as cv
# try:
#     import common_vars as cv
# except:
#     import phylogeny_utilities.common_vars as cv
from Bio import SeqIO

cog_lookup = {'COG0012':'Ribosome-binding ATPase YchF, GTP1/OBG family',
              'COG0016':'Phenylalanyl-tRNA synthetase alpha subunit','COG0018':'Arginyl-tRNA synthetase',
              'COG0048':'Ribosomal protein S12','COG0049':'Ribosomal protein S7','COG0052':'Ribosomal protein S2',
              'COG0080':'Ribosomal protein L11','COG0081':'Ribosomal protein L1',
              'COG0085':'DNA-directed RNA polymerase, beta subunit/140 kD subunit',
              'COG0087':'Ribosomal protein L3','COG0088':'Ribosomal protein L4',
              'COG0090':'Ribosomal protein L2','COG0091':'Ribosomal protein L22',
              'COG0092':'Ribosomal protein S3','COG0093':'Ribosomal protein L14','COG0094':'Ribosomal protein L5','COG0096':'Ribosomal protein S8','COG0097':'Ribosomal protein L6P/L9E','COG0098':'Ribosomal protein S5','COG0099':'Ribosomal protein S13','COG0100':'Ribosomal protein S11','COG0102':'Ribosomal protein L13','COG0103':'Ribosomal protein S9','COG0124':'Histidyl-tRNA synthetase','COG0172':'Seryl-tRNA synthetase','COG0184':'Ribosomal protein S15P/S13E','COG0185':'Ribosomal protein S19','COG0186':'Ribosomal protein S17','COG0197':'Ribosomal protein L16/L10AE','COG0200':'Ribosomal protein L15','COG0201':'Preprotein translocase subunit SecY','COG0202':'DNA-directed RNA polymerase, alpha subunit/40 kD subunit','COG0215':'Cysteinyl-tRNA synthetase','COG0256':'Ribosomal protein L18','COG0495':'Leucyl-tRNA synthetase','COG0522':'Ribosomal protein S4 or related protein','COG0525':'Valyl-tRNA synthetase','COG0533':'tRNA A37 threonylcarbamoyltransferase TsaD','COG0541':'Signal recognition particle GTPase','COG0552':'Signal recognition particle GTPase'}


codon_lookup = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                    'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
                    'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                    'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
                    'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
                    'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                    'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
                    'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'(stop)','TGA':'(stop)',
                    'TAG':'(stop)'}
codons_nostop = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                    'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
                    'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                    'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
                    'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
                    'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                    'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
                    'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'','TGA':'',
                    'TAG':''}



def read_from_fasta(file_path):
    """
    Reads from a fasta file and returns a dictionary where keys are taxon names and values are strings
    representing the sequence.
    :param file_path (string): The full system-readable path to the fasta file
    :return: fasta (dict)
    """
    output={}
    fasta=open(file_path,'r')
    first=True
    seq=''
    for l in fasta:
        if l[0]=='>':
            if first!=True:
                output[name]=seq
            else:
                first=False
            name=l[1:].strip()
            seq=''
        else:
            seq=seq + l.strip()
    output[name]=seq
    fasta.close()
    return output

def read_from_fastq(file_path):
    '''
    Basic FASTQ parser. Returns two dictionary objects. The first is a fasta-dict
    in the form {'<sequence_name>': <sequence String>}. The second has the quality
    scores in the form of a numpy array with data type np.uint8. The quality scores
    have had the bias subtracted (in this case 33).
    :param file_path:
    :return:
    '''
    fi = open(file_path,'r')
    ct=0
    fasta = {}
    quals = {}
    nm=''
    minlen=999999
    maxlen=0
    for ln in fi:
        if ct % 4 == 0:
            nm = ln.strip()[1:]
        if ct % 4 ==1:
            fasta[nm]=ln.strip()
        if ct % 4 ==2:
            assert ln[0]=='+'
        if ct % 4 ==3:
            quals[nm] = np.frombuffer(bytes(ln.strip(),'utf-8'),dtype=np.uint8) - 33
            dim=quals[nm].shape[0]
            minlen = min(minlen,dim)
            maxlen = max(maxlen,dim)
            nm=''
        ct += 1
    fi.close()
    print('max length: %s\tmin length: %s' % (maxlen,minlen))
    return fasta, quals

def write_to_fastq_subset(file_path, fastadict, qualsdict, subset_keys=None, quiet=False):
    seq_qual_keys = set(fastadict.keys()).intersection(set(qualsdict.keys()))
    nfa = len(fastadict)
    nqu = len(qualsdict)
    nfaqu = len(seq_qual_keys)
    if not quiet:
        print('# keys in... (fasta: %d), (quals: %d), (in common: %d)' % (nfa, nqu, nfaqu))

    if subset_keys is None:
        subset_keys_int = seq_qual_keys
        print('writing %d keys to file...' % nfaqu)
    else:
        nssk = len(subset_keys)
        subset_keys_int = seq_qual_keys.intersection(set(subset_keys))
        nsski = len(subset_keys_int)
        print('# subset keys in... (input: %d), (fasta/quals: %d), (in common: %d)' % (nssk, nfaqu, nsski))
        print('writing %d keys to file...' % nsski)

    fout = open(file_path, 'w')
    lct = 0

    for k in subset_keys_int:
        ct = fout.write('@%s\n%s\n+\n%s\n' % (k, fastadict[k], (qualsdict[k]+33).tobytes().decode('utf-8')) )
        lct += 1
        if lct % 100000 == 0:
            print('\rseqs done: %d' % lct, end = '')
    print('')
    fout.close()






def delete_taxa_and_remove_all_blank_columns_np(fasta_dict, subset_keys=None):
    '''

    :param fasta_dict:
    :return:
    '''
    tax, fanp = fasta_dict_to_nparray(fasta_dict)
    if subset_keys is not None:  # delete some taxa...
        remain_keys=list(set(fasta_dict.keys()).difference(set(subset_keys)))
        diff_keys=list(set(fasta_dict.keys()).difference(set(remain_keys)))
        diff_inds=list(map(lambda x: tax.index(x), diff_keys))
        diff_inds.sort()
        fanp = np.delete(fanp,diff_inds,0)
        tax = list(map(lambda x: tax[x], np.delete(np.arange(len(tax)),diff_inds,0)))

    # remove columns that are fully blank
    colsums=np.sum((fanp==45),0)
    nrows=fanp.shape[0]
    allblank=list(np.where(colsums==nrows)[0])
    fanp=np.delete(fanp,allblank,1)
    return tax,fanp



def remove_all_blank_columns_utils(fasta_dict,same_length_check=True):
    """
    Takes a dictionary representing a fasta file and removes any columns that are blank for all taxa. Data are
    assumed to be aligned starting with the first column in each string.

    NOTE: the operations in this function are in-place on the object provided by pythons nature, so while it
    returns a dictionary, catching the return value is not strictly necessary and the input will be
    modified after the fact.

    :param fasta_dict (dict): fasta dictionary (keys=taxon names, values=alignment strings)
    :param same_length_check (boolean) : OPTIONAL (default=True) If True, will throw an error if all sequences
        in fasta_dict are not the same length.
    :return: fasta_dict: dictionary with columns removed.
    """
    seqs_list = list(fasta_dict.values())
    num_seqs=len(seqs_list)
    seq_len = len(seqs_list[0])

    # check that all the sequences are the same length
    if same_length_check==True:
        for i in fasta_dict.values():
            if len(i) != seq_len:
                print ('The sequences were not all the same length.')
                return fasta_dict, -1

    # identify columns that are blank for every taxon
    all_blanks_list = []
    for i in range(seq_len):
        allblank=True
        for j in fasta_dict.values():
            if j[i]!='-':
                allblank=False
                break
        if allblank==True:
            all_blanks_list.append(i)

    newfasta={}
    for i in fasta_dict.keys():
        newfasta[i]=''

    non_blanks=list(set(range(seq_len)).difference(set(all_blanks_list)))
    # remove those columns (in place, so do it in reverse order)
    if len(all_blanks_list)>0:
        # all_blanks_list.sort(reverse=True)
        non_blanks.sort()

        for i in fasta_dict.keys():
            old=fasta_dict[i]
            new=''
            for j in non_blanks:
                new = new + old[j]
            newfasta[i]=new
    else:
        for i in fasta_dict.keys():
            newfasta[i]=fasta_dict[i]

    # return newfasta, all_blanks_list
    return newfasta

def read_fastq_to_qs_histogram(filepath, num_cols, num_rows):
    fastq = open(filepath,'r')
    qs = np.zeros((num_rows,num_cols),dtype=np.uint8)

    first=True
    seq=''
    ct = 0
    fastq.readline()
    for l in fastq:
        qv=np.frombuffer(bytes(l.strip(),'utf-8'),dtype=np.uint8) - 33
        qs[ct, :qv.shape[0]]=qv
        fastq.readline()
        ct += 1

    fastq.close()
    cts = np.ones((95,num_cols),dtype=np.int32)
    for c in range(num_cols):
        uq = np.unique(qs[:,c],return_counts=True)
        for j in range(uq[0].shape[0]):
            cts[uq[0][j]]=uq[1][j]
    return qs

def get_boxplot_stats(uqs):
    ncols = uqs.shape[1]
    nrows = uqs.shape[0]

    out=np.zeros((7,ncols),dtype=np.float64)
    # probs=uqs.cumsum(0)/np.sum(uqs,0)
    new = np.vstack((np.zeros(104, dtype=np.float64), uqs.cumsum(0)/np.sum(uqs,0)))
    probs = [.99,.95,.75,.5,.25,.05,.01]
    for c in range(ncols):
        out[:,c] = np.interp(probs, new[:,c], np.arange(new.shape[0]))

    return out



def list_from_stdin():
    l = []
    while True:
        a=input()
        if len(a)==0:
            break
        else:
            l.append(a)
    return l

def write_to_fasta(out_file_path,fasta_dict,subset_keys=None,raw=False,quiet=False):
    """
    Takes a fasta dictionary (keys=taxon names, values=alignment strings) and writes it to a file. Contains
    optional variables to specify only a subset of keys to write, or to write strings without blanks (assumed
    to be '-').
    :param out_file_path (string): system-readable path to write to
    :param fasta_dict (dict): fasta dictionary
    :param subset_keys: OPTIONAL iterator with values representing taxon names to include.
    :param raw: OPTIONAL (default: False) if TRUE, writes the "raw" (unaligned) sequence, which is the
    full sequence with blank ('-') characters removed.
    :return:
    """

    if subset_keys==None:
        mykeys=fasta_dict.keys()
    else:
        mykeys=list(set(subset_keys).intersection(fasta_dict.keys()))
        leftover_keys=list(set(subset_keys).difference(fasta_dict.keys()))
        if not quiet:
            print ('There were ' + str(len(leftover_keys)) + ' keys in the subset that were not in the original data.\n')


    fasta=open(out_file_path,'w')
    for i in mykeys:
        fasta.write('>'+i+'\n')
        if raw==False:
            fasta.write(fasta_dict[i]+'\n')
        else:
            fasta.write(fasta_dict[i].replace('-','')+'\n')
    fasta.close()
    if not quiet:
        print ('wrote file: ' + out_file_path + ' with the specified keys.')

# def write_to_fasta(out_file_path,fasta_dict,subset_keys=None,raw=False,quiet=False):
#     """
#     Takes a fasta dictionary (keys=taxon names, values=alignment strings) and writes it to a file. Contains
#     optional variables to specify only a subset of keys to write, or to write strings without blanks (assumed
#     to be '-').
#     :param out_file_path (string): system-readable path to write to
#     :param fasta_dict (dict): fasta dictionary
#     :param subset_keys: OPTIONAL iterator with values representing taxon names to include.
#     :param raw: OPTIONAL (default: False) if TRUE, writes the "raw" (unaligned) sequence, which is the
#     full sequence with blank ('-') characters removed.
#     :return:
#     """
#     if subset_keys==None:
#         mykeys=fasta_dict.keys()
#     else:
#         mykeys=list(set(subset_keys).intersection(fasta_dict.keys()))
#         leftover_keys=list(set(subset_keys).difference(fasta_dict.keys()))
#         if not quiet:
#             print ('There were ' + str(len(leftover_keys)) + ' keys in the subset that were not in the original data.\n')
#
#
#     fasta=open(out_file_path,'w')
#     for i in mykeys:
#         fasta.write('>'+i+'\n')
#         if raw==False:
#             fasta.write(fasta_dict[i]+'\n')
#         else:
#             fasta.write(fasta_dict[i].replace('-','')+'\n')
#     fasta.close()
#     if not quiet:
#         print ('wrote file: ' + out_file_path + ' with the specified keys.')


def read_from_aligned_phylip(file_path,as_nparray=False):
    f=open(file_path,'r')
    sz_ln = f.readline().strip().split(' ')

    #get size:
    rowsdone=False
    nrows = 0; ncols = 0;
    for i in range(len(sz_ln)):
        if sz_ln[i]!='' and rowsdone==False:
            nrows=int(sz_ln[i])
            rowsdone=True
        elif sz_ln[i]!='' and rowsdone==True:
            ncols = int(sz_ln[i])
            break

    assert ncols>0 and nrows>0
    arr = np.zeros((nrows,ncols),dtype=np.uint8)
    names = []

    #figure out sequence start positions
    ln = f.readline().strip()
    lctr = 0
    while ln[lctr]!=' ':
        lctr += 1
    while ln[lctr]==' ':
        lctr+=1
    names.append(ln[:lctr].strip())
    arr[0,:]=np.frombuffer(ln[lctr:(lctr+ncols)],np.uint8)

    for i in range(1,nrows):
        ln = f.readline().strip()
        names.append(ln[:lctr].strip())
        arr[i,:]=np.frombuffer(ln[lctr:(lctr+ncols)],np.uint8)

    if as_nparray==False:
        fasta={}
        for i in range(nrows):
            fasta[names[i]]=str(np.getbuffer(arr[i,:]))
        return fasta
    else:
        return names, arr

def get_fastadict_reverse_complement(fd):
    '''
    returns a fasta dictionary where every sequence is the reverse complement of the original.
    :param fd:
    :return:
    '''
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    fdrc = {}
    for i in fd.keys():
        fdrc[i] = str(Seq(fd[i],generic_dna).reverse_complement())
    return fdrc

def reverse_complement(seq):
    lkp = {'G':'C', 'g':'c', 'C':'G', 'c':'g', 'A':'T', 'a':'t', 'T':'A', 't':'a'}
    keys = set(list('GCATgcat'))
    a=''.join(map(lambda x: lkp[x] if x in keys else x, seq))
    return a[::-1]

def mask_fastadict(fasta, min_pct_nongap = 0.1):
    '''
    Takes a fasta dictionary and returns an equivalent version with the sequences masked
    based on 'min_pct_nongap'
    :param fasta:
    :param min_pct_nongap:
    :return:
    '''
    thresh = min_pct_nongap/100
    ntax = len(fasta.keys())
    ncols = len(fasta.values()[0])
    nparr = np.zeros((ntax,ncols),dtype=np.uint8)
    for i in range(ntax):
        seq = fasta[fasta.keys()[i]]
        nparr[i,:]=np.frombuffer(seq,np.uint8)

    # 45 is the uint8 code for the dash character '-':
    maskcols = np.where(np.sum((nparr!=45)*1,0).astype(np.float32)/nparr.shape[0]>thresh)
    newfasta = {}
    for i in range(ntax):
        k = fasta.keys()[i]
        newfasta[k] = str(np.getbuffer(nparr[i,maskcols]))
    return newfasta

def get_consensus_sequence(fasta, single_arb=True):
    '''
    Returns the consensus sequence from a multiple sequence alignment. If some oclumns have ties, reutnrs an
    arbitrary one if single_arb is True, and a list of all possibles if it is false. fasta can either be
    a list of aligned sequences or a fasta dict.
    :param fasta:
    :param single_arb:
    :return:
    '''
    from collections import Counter
    assert isinstance(fasta,list) or isinstance(fasta,dict)
    get_most_common = lambda test: [i[0] for i in test.most_common()[:sum(map(lambda x: 1 if x==test.most_common(1)[0][1] else 0, test.values()))]]
    # get_most_common = lambda test: (i[0] for i in test.most_common()[:sum(map(lambda x: 1 if x == test.most_common(1)[0][1] else 0, test.values()))])
    if isinstance(fasta,dict):
        l=len(next(iter(fasta.values())))
        assert l>2
        mostcomms=list(map(lambda x: get_most_common(Counter(map(lambda y: y[x], fasta.values()))), range(l)))
    else:
        # ***must be a list otherwise***
        l=len(fasta[0])
        assert l>2
        mostcomms = list(map(lambda x: get_most_common(Counter(map(lambda y: y[x], fasta))), range(l)))
    if single_arb:
        return ''.join(map(lambda x: x[0],mostcomms))
    else:
        f=[[a,b] for a in mostcomms[0] for b in mostcomms[1]]
        ct = 2
        while ct < len(mostcomms):
            f = [a + [b,] for a in f for b in mostcomms[ct]]
            ct += 1
        return list(map(lambda x: ''.join(x),f))


def get_avg_pdistance_of_fasta(fasta_path, getmax=False):
    f = read_from_fasta(fasta_path)
    taxn, fnp = fasta_dict_to_nparray(f)
    return get_avg_pdistance_of_nparray(fnp, getmax)

def get_avg_pdistance_of_list(seq_list, getmax=False):
    '''
    NOTE: all seqs in list must have the same length. If not, an error will arise.
    :param seq_list:
    :param getmax:
    :return:
    '''
    all_lens=list(set(map(len, seq_list)))
    assert (len(all_lens)==1), "sequences in seq_list are not all the same. some lengths are: %s" % all_lens[:min(5,len(all_lens))]
    fnp = np.zeros((len(seq_list), len(seq_list[0])), dtype=np.uint8)
    for i in range(len(seq_list)):
        fnp[i,:] = str2nparr(seq_list[i])
    return get_avg_pdistance_of_nparray(fnp)


def get_avg_pdistance_of_nparray(fnp, getmax=False, weighted=False):
    '''
    Computes average pairwise p-distance from an NP-array of uint8 representing an MSA.
    :param fnp:
    :param getmax: if True, return the max p-distance isntead of the average (default: False)
    :param weighted: if True, return the weighted p-distance instead of the raw average
    :return: (pd, equal_site_ct, common_site_ct, pair_ct)
    '''

    # f = read_from_fasta(fasta_path)
    # taxn, fnp = fasta_dict_to_nparray(f)

    maxpd = 0.
    run_tot = 0.
    run_ct = 0
    comm_sum=0
    same_sum=0
    for i in range(fnp.shape[0]):
        for j in range(i):
            comm = np.sum((fnp[i,:]!=45) & (fnp[j,:]!=45),dtype=np.float64)
            if comm ==0:
                continue
            same = np.sum((fnp[i,:]!=45) & (fnp[j,:]!=45) & (fnp[i,:]==fnp[j,:]),dtype=np.float64)
            run_tot+= 1.- same / comm
            run_ct += 1
            same_sum += same
            comm_sum += comm
            if (1. - same / comm > maxpd):
                maxpd = 1. - same / comm
            # print('i: %s\tj: %s\tcomm: %s\tsame: %s' % (i,j,comm,same))

    # print ('Avg P-Distance: %s' % pd)
    if getmax==False:
        if weighted:
            pd = 1.0 - same_sum / comm_sum
            return pd, int(same_sum), int(comm_sum), run_ct
        else:
            pd = run_tot / float(run_ct)
            return pd, int(same_sum), int(comm_sum), run_ct
    else:
        return maxpd

def compute_p_dist(np_seq1, np_seq2):
    comm = np.sum((np_seq1 != 45) & (np_seq2 != 45),dtype=np.float64)
    same = np.sum((np_seq1 != 45) & (np_seq2 != 45) & (np_seq1==np_seq2), dtype=np.float64)
    return (1-same/comm, int(comm), int(same))

def fasta_dict_to_nparray(fasta, taxnames=None, order='C'):
    '''
    Converts a dictionary representing an *aligned* fasta file (keys = seq names, values = seqs),
    into a numpy array of type np.uint8. Returns (taxa_names, array) where taxa_names is a list
    of names, in order, represented by the rows of the array. Should match the order of fasta.keys,
    but just in case...
    :param fasta:
    :param taxnames:
    :param order:
    :return:
    '''
    ntax = len(fasta.keys())
    ncols = max(map(len,fasta.values()))
    nparr = np.zeros((ntax,ncols),dtype=np.uint8,order=order)
    if taxnames is None:
        taxnames=[]
        # for i in range(ntax):
        i=0
        for k in fasta.keys():
            taxnames.append(k)
            seq = fasta[k]
            ls = len(seq)
            nparr[i,:ls]=np.frombuffer(bytes(seq,'utf8'),dtype=np.uint8)
            if ls < ncols:
                nparr[i,ls:]=45
            i+=1
    else:
        for i in range(len(taxnames)):
            seq = fasta[taxnames[i]]
            ls = len(seq)
            nparr[i,:ls]=np.frombuffer(bytes(seq,'utf8'),dtype=np.uint8)
            if ls < ncols:
                nparr[i,ls:]=45
    if order=='F':
        nparrF=np.array(nparr, copy=True, order='F')
        return taxnames, nparrF
    else:
        return taxnames, nparr

def get_fasta_nparray_sp_score(fa_arr, return_col_sps = False):
    '''
    Computes the sum-of-pairs score for a particular alignment. I.e. the number of pair-wise matching
    characters in an alignment.
    :param fa_arr:
    :return:
    '''
    from collections import Counter
    nseqs = fa_arr.shape[0]
    slen = fa_arr.shape[1]
    nc2 = lambda x: int(x*(x-1)/2)
    nc2_sum = lambda x: sum(map(nc2, x))

    sp_score = 0
    col_sps = []
    for c in range(slen):
        col_ctr = Counter(nparr2str(fa_arr[:,c]).replace('-',''))
        col_sps.append(nc2_sum(col_ctr.values()))
    sp_score = sum(col_sps)
    if not return_col_sps:
        return sp_score
    else:
        return sp_score, col_sps





def nparr2str(nparr):
    '''
    Converts a 1-d numpy array to a python string using utf-8 encoding. Ideally used on an
    array of type uint8
    :param nparr:
    :return:
    '''
    return nparr[:].tobytes().decode('utf-8')

def str2nparr(seq):
    '''
    Converts a python string to a 1d numpy array of type numpy.uint8.
    :param seq:
    :return:
    '''
    return np.frombuffer(bytes(seq, 'utf8'), dtype=np.uint8)

def pdist(s1, s2, all=False):
    '''
    returns the p-distance. Blanks are '-' ONLY.
    :param s1:
    :param s2:
    :return:
    '''
    assert len(s1)==len(s2), "lengths do not match"
    l=len(s1)
    a=np.zeros((2,l),dtype=np.uint8)
    a[0,]=str2nparr(s1)
    a[1,]=str2nparr(s2)
    c = np.sum((a[0,] != 45) & (a[1,] != 45))
    s = np.sum((a[0,] != 45) & (a[1,] != 45) & (a[0,]==a[1,]))
    if all:
        return 1.0 - float(s)/float(c), s, c
    else:
        return 1.0 - float(s)/float(c)


def compute_distinct_seq_aln_patterns(fasta_dict, return_numpy_array=False):
    '''
    Computes the number of distinct sequences and alignment patterns
    in an MSA, as RAxML might do to start up.
    :param fasta_dict:
    :return: tuple: (# seqs, # patterns) unless return_numpy_array is True,
    then it just returns the whole array (Fortran order!)
    '''
    nm2cls, eqcls = get_fasta_duplicate_datastruct(fasta_dict)
    # start by making a fasta_dict that dedups the sequences:
    fasm = dict(map(lambda x: (x['members'][0], x['seq']), eqcls.values()))
    n_seqs = len(fasm)
    tax, fasmarr = fasta_dict_to_nparray(fasm, order='F')
    seqs=[]
    for i in range(fasmarr.shape[1]):
        if np.sum(fasmarr[:,i]!=45)>0:
            seqs.append(nparr2str(fasmarr[:,i]))
    seqs_s=list(set(seqs))
    n_pats = len(seqs_s)
    if return_numpy_array:
        fasmarrF=np.zeros((n_seqs,n_pats),dtype=np.uint8, order='F')
        for i in range(n_pats):
            fasmarrF[:,i]=str2nparr(seqs_s[i])
        return fasmarrF
    else:
        return (n_seqs, n_pats)

def write_nparray_to_fasta(out_file_path, taxnames, fasta_nparr):
    f = open(out_file_path,'w')
    for i in range(len(taxnames)):
        c = f.write('>%s\n' % taxnames[i])
        # f.write('%s\n' % str(np.getbuffer(fasta_nparr[i,:])))
        c = f.write(fasta_nparr[i,:].tobytes().decode('utf-8') + '\n')
    # print ('wrote %s taxa names and sequences to the fasta file: %s' % (len(taxnames),out_file_path))
    f.close()


def convert_raxml_reduced_to_fasta(raxml_reduced, fasta_output_path):
    red_f = open(raxml_reduced,'r')
    lw=red_f.readline()

    fa={}
    for ln in red_f:
        seq = ln.strip().split(" ")
        fa[seq[0]]=seq[1]

    write_to_fasta(fasta_output_path,fa)


def get_min_max_avg_sequence_length(fasta_file):
    a = read_from_fasta(fasta_file)
    mylens = map(len,a.values())
    m1 = max(mylens)
    m2 = min(mylens)
    mavg = float(sum(mylens))/float(len(mylens))
    print ('%s --- min: %s, max: %s, avg: %.2f' % (fasta_file,m2, m1, mavg))


def sepp_json_to_tsv(in_path, out_path):
    if out_path[-4:] not in ('.tsv','.txt','.tab'):
        out_path += '.txt'
    inf=open(in_path,'r')
    sepp_json = json.load(inf)
    inf.close()
    outf = open(out_path,'w')
    headers = ['placements_ID','nm[0]','nm[1]']
    flds=['edge', 'logL', 'prob', 'distal', 'pendant']
    for i in range(7):
        headers+=list(map(lambda x: x+'_'+str(i),flds))
    headers_str='\t'.join(headers) + '\n'
    outf.write(headers_str)
    id = 0
    sub_str = '\t'.join(['%s',]*38) + '\n'
    for pl in sepp_json['placements']:
        all_pls = pl['p'][0]
        if len(pl['p'])>1:
            for pl_i in range(1,len(pl['p'])):
                all_pls += pl['p'][pl_i]
        if len(pl['p'])<7:
            for pl_i in range(7-len(pl['p'])):
                all_pls += [-1, 0.0, 0.0, 0.0, 0.0]
        for nm_id in range(len(pl['nm'])):
            rw = [id,pl['nm'][nm_id][0],pl['nm'][nm_id][1]]
            rw += all_pls
            outf.write(sub_str % tuple(rw))
            # print(len(rw))
            # print(rw)
            # print(sub_str % tuple(rw))
        id += 1
    outf.close()

def make_histogram_of_branch_lengths(tree_file, out_path):
    '''
    Makes a histogram of branch lengths over a tree.
    :param tree_file: Path to newick file representing a Tree.
    :param out_path: path to write PDF to. If this does not end with '.pdf', that extension is added.
    :return:
    '''
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tr=dendropy.Tree.get(path=tree_file,schema='newick')
    internals = [i.length for i in tr.postorder_internal_edge_iter(filter_fn=lambda x: x.length is not None)]
    leaves = [i.length for i in tr.leaf_edge_iter(filter_fn=lambda x: x.length is not None)]

    fig, axs = plt.subplots(2, 1, sharex=True, tight_layout=True)
    axs[0].hist(internals, bins=20)
    axs[0].set_xlabel('ln [internal]')
    axs[1].hist(internals, bins=20)
    axs[1].set_xlabel('ln [leaves]')

    if out_path[-4:]!='.pdf':
        out_path += '.pdf'
    plt.savefig(out_path)
    plt.clf()


def make_histogram_of_sequence_lengths(fasta_file,out_figure_path,subtitle=None):
    '''
    Does pretty much what it says. Saves the file in format dependent on the path, so best
    to end the path with '.pdf'
    :param fasta_file:
    :param out_figure_path:
    :return:
    '''
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    a = read_from_fasta(fasta_file)
    lsnp = np.array(list(map(len,a.values())))
    q99=np.quantile(lsnp,0.99)
    q01=np.quantile(lsnp,0.01)
    m=np.mean(lsnp)


    # n, bn, pat=plt.hist(lsnp,'auto',range=(min(q01,m*0.75),max(q99,m*1.25)))
    n, bn, pat = plt.hist(lsnp, 65, range=(min(q01, m * 0.75), max(q99, m * 1.25)))
    # print(bn.shape)
    plt.xlabel('# Base Pairs')
    plt.ylabel('count')
    plt.suptitle('Seq-Len Histogram: %s' % fasta_file)
    if subtitle is not None:
        plt.title(subtitle)
    plt.savefig(out_figure_path)
    plt.clf()

def make_fasta_with_clean_names(fastadict,outfasta,outnames):
    '''
    Creates a fasta file that has names as simple integer strings, starting at 1.
    :param fastadict:
    :param outfasta:
    :param outnames:
    :return:
    '''
    outn=open(outnames,'w')
    newfasta={}
    counter=1
    for i in fastadict.keys():
        outn.write(str(counter) + '\t' + i + '\n')
        newfasta[str(counter)]=fastadict[i]
        counter+=1
    write_to_fasta(outfasta,newfasta)

def read_from_fasta_dedup(file_path):
    """
    Reads from a fasta file and returns a dictionary where keys are taxon names and values are strings
    representing the sequence.
    :param file_path (string): The full system-readable path to the fasta file
    :return: fasta (dict)
    """
    output={}
    namects={}
    fasta=open(file_path,'r')
    first=True
    for l in fasta:
        if l[0]=='>':
            if first!=True:
                if name in namects.keys():
                    namects[name]+=1
                else:
                    namects[name]=0
                output[name+'-' +str(namects[name])]=seq
            else:
                first=False
            name=l[1:].strip()
            seq=''
        else:
            seq=seq + l.strip()
    output[name]=seq
    fasta.close()
    return output

def seq_length_data_for_histogram(a):
    '''
    returns a list of integer sequence lengths that can be written to a file and imported in R
    :param a:
    :return:
    '''
    return list(map(len, a.values()))

def bp_genbank_to_fasta(gbk_in):
    fasta_out = gbk_in[:-4] + '.fna'
    input_handle = open(gbk_in, "rU")
    output_handle = open(fasta_out, "w")

    sequences = SeqIO.parse(input_handle, "genbank")
    count = SeqIO.write(sequences, output_handle, "fasta")

    output_handle.close()
    input_handle.close()
    return count

def bp_genbank_get_CDS_dict(filename):
    '''
    Does some kind of parsing of a genbank file...
    :param filename:
    :return:
    '''
    # print filename
    p, f = os.path.split(filename)
    used_temp = False
    if f[-3:]=='.gz':
        temp_fn = '/dev/shm/' + f[:-3]
        used_temp = True
        os.system('gunzip -c %s > %s' % (filename,temp_fn))
        target_filename = temp_fn
    else:
        target_filename = filename
    seq = enumerate(SeqIO.parse(target_filename,"genbank"))
    # ind, rec = next(seq)
    # p,f = os.path.split(filename)
    cds = {}
    dna = {}
    name_pos = {}
    # gene_names = {}
    # gene_syns = {}
    no_translation_loci=[]
    length_errors=[]
    while True:
        try:
            ind, rec = next(seq)
        except:
            break
        # Getting Sequence ID
        seqid = rec.__dict__.get('id',None)

        # Getting NCBI Taxon ID
        taxid = None
        if rec.features[0].type == 'source' and 'db_xref' in rec.features[0].qualifiers.keys():
            for dx in rec.features[0].qualifiers['db_xref']:
                if dx[0:5]=='taxon':
                    taxid = dx.replace('taxon:','')
                    continue

        project_id = None
        for i in rec.dbxrefs:
            if i[:11]=='BioProject:':
                lr = len(i)
                project_id=rec.dbxrefs[0][-(lr - 11):]

        for i in rec.features:
            if i.type=='CDS':
                st = i.location.start + 0
                en = i.location.end + 0
                dir = i.location.strand

                id = None
                protein_id = None
                locus_tag = None
                if 'protein_id' in i.qualifiers.keys():
                    id = 'pt:' + i.qualifiers['protein_id'][0]
                    protein_id = i.qualifiers['protein_id'][0]
                if 'locus_tag' in i.qualifiers.keys():
                    locus_tag = i.qualifiers['locus_tag'][0]
                    if id is None:
                        id = 'lt:' + i.qualifiers['locus_tag'][0]
                else:
                    id = 'na:' + f + '_' + str(st) + '-' + str(en)
                    # id = 'file-' + str()

                gene_id = None
                dbx = i.qualifiers.get('db_xref',None)
                if dbx is not None:
                    for entry in dbx:
                        if entry[0:6].lower()=='geneid':
                            gene_id = entry[7:]


                # print '%s, %s, %s, %s' %(str(st), str(en), str(dir), id)
                try:
                    prot = i.qualifiers['translation'][0]
                except:
                    no_translation_loci.append(id)
                    continue
                if (en-st-3)/3!=len(prot):
                    length_errors.append(id)
                if dir == -1:
                    dnastr = rec.seq[st:en].reverse_complement()
                else:
                    dnastr = rec.seq[st:en]

                # if 'gene' in i.qualifiers.keys():
                #     gene_names[id]=i.qualifiers['gene'][0]
                # if 'gene_synonym' in i.qualifiers.keys():
                #     gene_syns[id] = i.qualifiers['gene_synonym'][0]
                g_nm = i.qualifiers.get('gene',[None,])[0]
                g_nm_sy = i.qualifiers.get('gene_synonym', [None, ])[0]
                cds[id] = prot
                dna[id] = dnastr

                fetch_mg_name_format = '%s.%s protein_id=\"%s\" gene_id=\"%s\" project_id=\"%s\"'
                fetchMG_name = fetch_mg_name_format % (taxid, locus_tag, protein_id, gene_id, project_id)
                name_pos[id] = (seqid, g_nm, st, en, dir, g_nm_sy, fetchMG_name)

    if used_temp:
        os.system('rm %s' % temp_fn)
    return cds, dna, name_pos, no_translation_loci, length_errors, f[:-3]

def fasta_dict_to_string(fasta_dict):
    a=''
    for i in fasta_dict.keys():
        a=a+'>'+i + '\n'
        a=a + fasta_dict[i] + '\n'
    return a

def run_fastSP_on_two_fastas(ref_file,est_file,out=None):
    import subprocess
    if out==None:
        outf=open('temp.txt','w')
        out='temp.txt'
    else:
        outf=open(out,'w')
    subprocess.call(["java","-jar","/home/mikenute/Phylolab/share/code/misc/FastSP_1.6.0.jar","-r",
        ref_file,"-e",est_file], stdout=outf,stderr=outf)
    outf.close()

    p,f=os.path.split(os.path.abspath(est_file))
    args=read_alignment_results_file(out,f)
    # os.unlink('temp.txt')
    return args

def read_alignment_results_file(filepath=None,optional_filename=None):
    pa, fi = os.path.split(filepath)
    myargs = {}
    if optional_filename!=None:
        fi=optional_filename
    myargs['file_name'] = fi

    # open the results file
    textfile = open(filepath, 'r')
    text = textfile.read()
    textfile.close()

    # regexes to pull the results
    file_regex = 'SP-Score (?P<sp>\d+\.\d+[E\-\d]*)[.$\n]*Modeler (?P<modeler>\d+\.\d+[E\-\d]*)[.$\n]*SPFN (?P<spfn>\-*\d+\.\d+[E\-\d]*)[.$\n]*SPFP (?P<spfp>\-*\d+\.\d+[E\-\d]*)[.$\n]*Compression (?P<comp>\d+\.\d+[E\-\d]*)[.$\n]*TC (?P<tc>\d+\.\d+[E\-\d]*)'
    file_regex_2 = 'MaxLenNoGap= (?P<maxlen>\d+).*NumSeq= (?P<numseq>\d+).*LenRef= (?P<lenref>\d+).*LenEst= (?P<lenest>\d+).*Cells= (?P<cells>\d+)'
    fileregs = re.compile(file_regex)
    fileregs2 = re.compile(file_regex_2)

    vals1list = ['sp', 'modeler', 'spfn', 'spfp', 'comp', 'tc']
    vals2list = ['maxlen', 'numseq', 'lenref', 'lenest', 'cells']
    vals1 = fileregs.search(text)
    vals2 = fileregs2.search(text)
    for i in vals1list:
        try:
            myargs[i] = vals1.group(i)
        except:
            myargs[i] = ''
    for i in vals2list:
        try:
            myargs[i] = vals2.group(i)
        except:
            myargs[i] = ''
    return myargs

def get_tree_error(ref_tree,est_tree):
    import dendropy
    rt=dendropy.Tree()
    rt.read_from_path(ref_tree,schema='newick')
    et=dendropy.Tree()
    et.read_from_path(est_tree, schema='newick')
    args={}
    b=rt.false_positives_and_negatives(et)
    num_splits=0
    for i in rt.get_edge_set():
        if i.is_internal():
            num_splits+=1
    num_splits_et=0
    for i in et.get_edge_set():
        if i.is_internal():
            num_splits_et+=1
    args['tree_fp']=b[0]
    args['tree_fn']=b[1]
    args['num_splits_rt']=num_splits
    args['num_splits_et']=num_splits_et

    return args

def make_alignment_header_and_entryline():
    vals1list = ['sp', 'modeler', 'spfn', 'spfp', 'comp', 'tc']
    vals2list = ['maxlen', 'numseq', 'lenref', 'lenest', 'cells']

    str_header='\t'.join(vals1list)+'\t'
    str_header=str_header + '\t'.join(vals2list) + '\n'

    str_line='%('+str_header.strip().replace('\t',')s\t%(') + ')s\n'
    return str_header, str_line

def write_results_lines_to_file(args_list,out_file):
    outf=open(out_file,'w')
    vals1list = ['sp', 'modeler', 'spfn', 'spfp', 'comp', 'tc']
    vals2list = ['maxlen', 'numseq', 'lenref', 'lenest', 'cells']
    tree_list = ['tree_fp','tree_fn','num_splits_rt','num_splits_et']
    str_header='file_name\t'+'\t'.join(vals1list)+'\t'
    str_header=str_header + '\t'.join(vals2list) + '\t' + '\t'.join(tree_list) + '\n'
    outf.write(str_header)

    str_line='%('+str_header.strip().replace('\t',')s\t%(') + ')s\n'
    # print str_line
    # print args_list[0]
    for i in args_list:
        outf.write(str_line % i)

    outf.close()

def open_json(fn):
    '''
    Wraps json.load(). I.e., opens a json file, converts it to the equivalent python object using json.load(), then
    closes the file.
    :param fn:
    :return:
    '''
    t=open(fn,'r')
    fs=json.load(t)
    t.close()
    return fs

def dna_to_protein(dna,verbose=False):
    '''
    Alias for dna_to_protein_fast...
    :param dna:
    :param verbose:
    :return:
    '''
    # cv.codon lookup table transcribed from wikipedia
    # cv.codon_lookup = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    #                 'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
    #                 'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    #                 'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
    #                 'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
    #                 'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    #                 'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
    #                 'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'(stop)','TGA':'(stop)',
    #                 'TAG':'(stop)'}
    k = int(len(dna)/3)
    # extras = len(dna)-k*3
    # if extras>0 and verbose:
    #     print ("The sequence given is not divisible by 3, there are %i extra nucleotides, which will be ignored" % extras)
    str_protein = dna_to_protein_fast(dna)
    # for i in range(k):
    #     codon = dna[(i*3):(i*3+3)].upper()
    #     str_protein = str_protein + codons_nostop[codon]
    # str_protein = str_protein.replace('(stop)','')
    return str_protein

def dna_to_protein_fast(s):
    '''
    Converts a DNA string to a protein string using the codon table from Wikipedia (as of 2017).
    :param s:
    :return:
    '''
    return ''.join(map(lambda x: codons_nostop.get(x,'X'), map(''.join, zip(*[iter(s.upper())]*3))))

def get_max_fasta_seqlen(fasta_file):
    '''
    Gets the dimensions of a fasta file
    :param fasta_file:
    :return: returns a tuple with (# of columns, # of rows) in
        the fasta.
    '''
    fasta=open(fasta_file,'r')
    first=True
    seq=''
    ml = 0
    seq_ct = 0
    for l in fasta:
        if l[0]=='>':
            seq_ct+=1
            if first!=True:
                ml = max(ml,len(seq))
            else:
                first=False
            name=l[1:].strip()
            seq=''
        else:
            seq=seq + l.strip()
    fasta.close()
    return (ml, seq_ct)

def get_fasta_duplicate_datastruct(fasta_dict, quiet=False):
    '''
    Returns a pair of dictionaries. The first is {orig_seq_id: seq_dedup_id} (i.e. equivalence class).
    The second is {seq_dedup_id: {'seq': seq_string, 'members': [orig_seq_id 1, ...., ], 'copynum': num_dupes}}
    :param fasta_dict:
    :return:
    '''
    seqs = list(set(fasta_dict.values()))
    seq_to_index = dict(zip(seqs, list(range(len(seqs)))))
    eq_classes = dict(zip(list(range(len(seqs))), map(lambda x: {'seq': seqs[x], 'members': [], 'copynum': 0}, range(len(seqs)))))
    seq_nm_to_class = dict(zip(fasta_dict.keys(), map(lambda x: seq_to_index[fasta_dict[x]], fasta_dict.keys())))

    orig_ct = 0
    for (k,v) in seq_nm_to_class.items():
        eq_classes[v]['members'].append(k)
        eq_classes[v]['copynum']+=1
        orig_ct += 1

    singleton_ct = sum(map(lambda x: 1 if x['copynum']<=1 else 0, eq_classes.values()))

    if not quiet:
        print("# seqs in original: %d" % orig_ct)
        print('# seqs deduped:     %d' % len(seqs))
        print('# singletons:       %d' % singleton_ct)
    return seq_nm_to_class, eq_classes

def get_fasta_deduped(fasta_dict):
    '''
    Just returns the same type of data structure but with no duplicate sequences. Keys assigned
    to sequences originally containing duplicates is done arbitrarily.
    :param fasta_dict:
    :return:
    '''
    n2c, eq = get_fasta_duplicate_datastruct(fasta_dict)
    c2n_dedup = {}
    for k,v in n2c.items():
        c2n_dedup[v]=k
    return dict(map(lambda x: (c2n_dedup[x], eq[x]['seq']), c2n_dedup.keys()))

def get_list_from_file(filepath):
    '''
    Converts a text file to a python list of strings, with each line becoming one item
    in the list.
    :param filepath:
    :return:
    '''
    myf=open(filepath,'r')
    ol=[]
    for i in myf:
        if i.strip()!='':
            ol.append(i.strip())
    myf.close()
    return ol

def write_list_to_file(mylist=None,filepath=None):
    '''
    Writes a list object to a file, one item per line.
    :param mylist: python list object
    :param filepath: file path for output
    :return:
    '''
    if mylist is None or filepath is None:
        print("usage: write_list_to_file(list,filepath)")
    myf=open(filepath,'w')
    for i in mylist:
        myf.write(str(i) + '\n')
    myf.close()

def write_dict_to_file(mydict, filepath, delimiter='\t'):
    '''
    Takes a dictionary object and writes it to a file with each entry one on line, in string form,
    separated by a delimiter.
    :param mydict:
    :param filepath:
    :param delimiter:
    :return:
    '''
    # if delimiter is None:
    #     delimiter = '\t'
    myf = open(filepath,'w')

    for k in mydict.keys():
        myf.write(k)
        myf.write(delimiter)
        myf.write(str(mydict[k]))
        myf.write("\n")
    myf.close()

# formats a byte as 8-digit binary from a bytearray. Index assumed to be 0 unless specified.
bin8 = lambda x: ''.join(['0',]*(8-(len(format(x,'b'))%8) ))+format(x,'b')

def get_file_md5_digest(file_path):
    import hashlib
    with open(file_path,'rb') as fb:
        hd=hashlib.md5(fb.read()).hexdigest()
    return hd

def alignment_stats_dir_to_tabd(aln_dir,outpath):
    '''
    Takes a directory full of alignment_stats files and converts it to a tab-delimited text spreadsheet
    :param aln_dir:
    :param outpath:
    :return:
    '''
    myf=open(outpath,'w')
    myf.write('File\tALIGNMENT_PERCENT_BLANKS_MARKER\tALIGNMENT_BLANKS_COUNT\tALIGNMENT_GAPS_COUNT\tALIGNMENT_ROWS_COUNT\tALIGNMENT_COLUMNS_COUNT\tALIGNMENT_AVERAGE_GAPS_PER_SEQUENCE\tALIGNMENT_AVERAGE_GAP_LENGTH\tALIGNMENT_STDDEV_GAP_LENGTH\tALIGNMENT_MEDIAN_GAP_LENGTH\tALIGNMENT_MNHD\tALIGNMENT_ANHD\n')

    ln='%(file)s\t%(APBM)s\t%(ABC)s\t%(AGC)s\t%(ARC)s\t%(ACC)s\t%(AAGPS)s\t%(AAGL)s\t%(ASGL)s\t%(AMGL)s\t%(AM)s\t%(AA)s\n'
    for i in os.listdir(aln_dir):
        a=read_alignment_stats_file(aln_dir + '/' + i)
        myf.write(ln % a)
    myf.close()

def read_alignment_stats_file(filepath):
    myf=open(filepath,'r')
    text=myf.read()
    myf.close()

    file_regex = 'ALIGNMENT_PERCENT_BLANKS_MARKER \| (?P<APBM>\d+\.\d+[E\-*\d]*)\nALIGNMENT_BLANKS_COUNT \| (?P<ABC>\d+\.\d+[E\-*\d]*)\nALIGNMENT_GAPS_COUNT \| (?P<AGC>\d+\.\d+[E\-*\d]*)\nALIGNMENT_ROWS_COUNT \| (?P<ARC>\d+\.\d+[E\-*\d]*)\nALIGNMENT_COLUMNS_COUNT \| (?P<ACC>\d+\.\d+[E\-*\d]*)\nALIGNMENT_AVERAGE_GAPS_PER_SEQUENCE \| (?P<AAGPS>\d+\.\d+[E\-*\d]*)\nALIGNMENT_AVERAGE_GAP_LENGTH \| (?P<AAGL>\d+\.\d+[E\-*\d]*)\nALIGNMENT_STDDEV_GAP_LENGTH \| (?P<ASGL>\d+\.\d+[E\-*\d]*)\nALIGNMENT_MEDIAN_GAP_LENGTH \| (?P<AMGL>\d+\.\d+[E\-*\d]*)\nALIGNMENT_MNHD \| (?P<AM>\d+\.\d+[E\-*\d]*)\nALIGNMENT_ANHD \| (?P<AA>\d+\.\d+[E\-*\d]*)\n'
    # file_regex = 'ALIGNMENT_PERCENT_BLANKS_MARKER\.*(?P<apbm>\d+\.\d+)\.*ALIGNMENT_BLANKS_COUNT\.*(?P<abc>\d+\.\d+*)\.*ALIGNMENT_GAPS_COUNT\.*(?P<agc>\d+\.\d+)'
    fileregs = re.compile(file_regex)
    vals1 = fileregs.search(text)
    try:
        args= {
            'file':filepath,
            'APBM': vals1.group('APBM'),
            'ABC': vals1.group('ABC'),
            'AGC': vals1.group('AGC'),
            'ARC': vals1.group('ARC'),
            'ACC': vals1.group('ACC'),
            'AAGPS': vals1.group('AAGPS'),
            'AAGL': vals1.group('AAGL'),
            'ASGL': vals1.group('ASGL'),
            'AMGL': vals1.group('AMGL'),
            'AM': vals1.group('AM'),
            'AA': vals1.group('AA')
        }
        return args
    except:
        print (filepath)
        args= {
            'file':filepath,
            'APBM': '',
            'ABC': '',
            'AGC': '',
            'ARC': '',
            'ACC': '',
            'AAGPS': '',
            'AAGL': '',
            'ASGL': '',
            'AMGL': '',
            'AM': '',
            'AA': ''
        }
        return args

def shrink_fasta_to_complement_of_another(bigfile,subtractfile,outfile):
    f1=read_from_fasta(bigfile)
    f2=read_from_fasta(subtractfile)
    newkeys=list(set(f1.keys()).difference(set(f2.keys())))
    write_to_fasta(outfile,f1,newkeys)

def get_dict_from_general_delimited_file(filename, key_col=0, val_col=1, delimiter='\t', header_rows=0):
    myf = open(filename, 'r')
    for r in range(header_rows):
        foo = myf.readline()
    args={}
    ct = 0
    for ln in myf:
        if len(ln.strip()) > 0:
            a = ln.strip().split(delimiter)
            args[a[key_col]] = a[val_col]
            ct += 1
            if ct % 100000 == 0:
                print('\r',end=''); print('%9d lines done' % ct, end = '')
    print('\n')
    myf.close()
    return args

    # args = dict(map(lambda x: (x[key_col], x[val_col]),
    #                     [ln.strip().split(delimiter) for ln in myf if len(ln.strip()) > 0]))

def get_dict_from_file_fast(filename, delimiter='\t', keysfirst = True):
    myf = open(filename,'r')
    args = {}
    if keysfirst:
        args = dict(map(lambda x: x.strip().split(delimiter),[ln for ln in myf if len(ln.strip())>0]))
    else:
        args = dict(map(lambda x: (x[1],x[0]), [ln.strip().split(delimiter) for ln in myf if len(ln.strip()) > 0]))
    # for i in myf:
    #     if len(i.strip())>0:
    #         a=i.strip().split(delimiter)
    #         args[a[0]]=a[1]
    myf.close()
    return args



def get_dict_from_tab_delimited_file(filename, keysfirst = True):
    myf = open(filename,'r')
    args = {}
    ct = 0
    if keysfirst==True:
        for i in myf:
            if len(i.strip())>0:
                a=i.strip().split('\t')
                args[a[0]]=a[1]
    else:
        for i in myf:
            if len(i.strip())>0:
                a=i.strip().split('\t')
                args[a[1]]=a[0]
                ct +=1
                if ct % 10000000 == 0:
                    print('%s lines done' % ct)
    myf.close()
    return args

def int_vs_node_brlens(fp):
    '''
    Takes a file path containing a newick tree and prints the Counts and, Average,
    Standard Deviation and C.V. of branch lengths for the internal and leaf
    branches separately.
    :param fp:      (string) path to file containing newick tree
    :return:        None
    '''
    tr = dendropy.Tree.get(path=fp,schema='newick',preserve_underscores=True)
    totlen=0
    totlen2=0
    ct=0
    for nd in tr.preorder_internal_node_iter():
        if nd.edge_length is not None:
            totlen2+=nd.edge_length**2
            totlen+=nd.edge_length
            ct+=1
    totlen_lf2=0
    totlen_lf=0
    ct_lf=0
    for nd in tr.leaf_node_iter():
        totlen_lf2+=nd.edge_length**2
        totlen_lf+=nd.edge_length
        ct_lf+=1
    int_avg = totlen / ct
    int_sd = (totlen2 / ct - (int_avg)**2)**0.5
    int_cv = int_avg / int_sd
    lf_avg = totlen_lf / ct_lf
    lf_sd = (totlen_lf2 / ct_lf - (lf_avg)**2)**0.5
    lf_cv = lf_avg / lf_sd
    print("                Ct           \tAvg Len      \tSD Len       \tCoef Var      ")
    print("Internal Nodes  %s         \t%.5f  \t%.5f  \t%.5f" % (ct, int_avg, int_sd, int_cv))
    print("Leaf Nodes      %s         \t%.5f  \t%.5f  \t%.5f" % (ct_lf, lf_avg, lf_sd, lf_cv))



if __name__=='__main__':
    if '-h' in sys.argv:
        help()
    elif '-pbs' in sys.argv:
        pbsargs=parse_pbs_args()
        make_pbs_file(*pbsargs)
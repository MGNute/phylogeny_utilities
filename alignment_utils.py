#-------------------------------------------------------------------------------
# Name:        alignment_utils.py
# Purpose:
#
# Author:      Michael
#
# Created:     10/02/2015
# Copyright:   (c) Michael 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import platform, os, io
import re
import dendropy
import numpy as np
# import pasta_output_readers as pouts
# import utilities



def fastsp_run_on_two_fastas(ref_file,est_file,out=None):
    import subprocess

    if out==None:
        out=open('temp.txt','w')
    else:
        out=open(out,'w')
    subprocess.call(["java","-jar","/home/mikenute/Phylolab/share/code/misc/FastSP_1.6.0.jar","-r",
        ref_file,"-e",est_file], stdout=out,stderr=out)
    out.close()
    out = open('temp.txt','r')
    text=out.read()
    print (text)
    out.close()

    # text=out.getvalue()
    # p,f=os.path.split(ref_file)
    # p2,f2=os.path.split(est_file)
    args=fastsp_read_results_string(text)
    args['ref_file']=ref_file
    args['est_file']=est_file
    os.unlink('temp.txt')
    return args

def fastsp_read_results_file(filepath=None,optional_filename=None):
    pa, fi = os.path.split(filepath)
    myargs = {}
    if optional_filename!=None:
        fi=optional_filename
    # myargs['file_name'] = fi

    # open the results file
    textfile = open(filepath, 'r')
    text = textfile.read()
    textfile.close()

    args=fastsp_read_results_string(text)
    args['file_name']=fi
    return args


def fastsp_read_results_string(text, ref_file=None, est_file=None):
    import re
    myargs = {}

    # open the results file
    # textfile = open(filepath, 'r')
    # text = textfile.read()
    # textfile.close()

    # regexes to pull the results
    file_regex = 'SP-Score (?P<sp>\d+\.\d+[E\-\d]*)[.\n]*Modeler (?P<modeler>\d+\.\d+[E\-*\d]*)[.\n]*SPFN (?P<spfn>\d+\.\d+[E\-\d]*)[.\n]*SPFP (?P<spfp>\d+\.\d+[E\-\d]*)[.\n]*Compression (?P<comp>\d+\.\d+[E\-\d]*)[.\n]*TC (?P<tc>\d+\.\d+[E\-\d]*)'
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

def fastsp_results_folder_to_tab_delimited(folder=None, output_file_path=None):
    # quick checks on folder string
    folder = folder.replace('\\', '/') + '/'
    folder = folder.replace('//', '/')

    # model output line:
    line_str = '%(file_name)s\t%(sp)s\t%(modeler)s\t%(spfn)s\t%(spfp)s\t%(comp)s\t%(tc)s\t%(maxlen)s\t%(numseq)s\t%(lenref)s\t%(lenest)s\t%(cells)s\t\n'

    outfile = open(output_file_path, 'w')
    outfile.write(
        'file_name\tSP-Score\tModeler\tSPFN\tSPFP\tCompression Ratio\tTC\tMaxLenNoGap\tNumSeq\tLenRef\tLenEst\tCells\t\n')

    for i in os.listdir(folder):
        myargs = {}
        myargs = fastsp_read_results_file(folder + i)
        outfile.write(line_str % myargs)

    outfile.close()

def mask_fastadict(fasta, min_pct_nongap = 0.1):
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

def read_from_fasta(file_path):
    '''
    Reads a fasta file into a dictionary where the keys are sequence names and values are strings representing the
    sequence.

    :param file_path:
    :return:
    '''
    output={}
    fasta=open(file_path,'r')
    first=True
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

def write_to_fasta(out_file_path=None,fasta_dict=None,subset_keys=None,raw=False):
    '''
    takes a dictionary where keys are the sequence names and values are the sequences, and writes that to a FASTA file.
    Optionally can include only a subset of the full set of taxa (defaults to all), and optionally can remove all gaps
    before writing.

    :param out_file_path:
    :param fasta_dict:
    :param subset_keys:
    :param raw:
    :return: None
    '''
    if out_file_path==None:
        a = '''
Arguments (in order):
    out_file_path:
    fasta_dict:
    subset_keys:
    raw:
        '''
        return None

    if subset_keys==None:
        mykeys=fasta_dict.keys()
    else:
        mykeys=subset_keys

    fasta=open(out_file_path,'w')
    for i in mykeys:
        fasta.write('>'+i+'\n')
        if raw==False:
            fasta.write(fasta_dict[i]+'\n')
        else:
            fasta.write(fasta_dict[i].replace('-','')+'\n')

    fasta.close()

class alignment_from_fasta:
    '''
    This is a clumsy class that represents a fasta file. It could easily be written as a smaller set of
    functions (and has been in another file) but this works for now.
    '''
    def __init__(self):
        self.taxa={}
        self.filename = None
        self.folder = None
        self.path = None

    def assign_fasta(self,file_path):
        self.sequence=self.read_from_fasta(file_path)
        self.folder, self.filename = os.path.split(file_path)
        self.path = file_path

    def read_from_fasta(self,file_path):
        output={}
        fasta=open(file_path,'r')
        first=True
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

    def write_to_fasta(self,out_file_path=None,fasta_dict=None,subset_keys=None,raw=False):
        if fasta_dict==None:
            fasta_dict=self.sequence

        if subset_keys==None:
            mykeys=fasta_dict.keys()
        else:
            mykeys=subset_keys

        fasta=open(out_file_path,'w')
        for i in mykeys:
            fasta.write('>'+i+'\n')
            if raw==False:
                fasta.write(fasta_dict[i]+'\n')
            else:
                fasta.write(fasta_dict[i].replace('-','')+'\n')

        fasta.close()

    def count_all_blank_columns(self):
        num_seqs=len(self.sequence.values())
        seq_len=len(self.sequence.values()[0])
        print (num_seqs)
        print (seq_len)

        record=io.StringIO.StringIO()
        record.write('line\tblank\toccupied\n')

        for i in self.sequence.keys():
            if len(self.sequence[i])!=seq_len:
                print (i + ' has sequence length different from the first')
                print (len(self.sequence[i]))

        for i in range(seq_len):
            all_blank=True
            blank_count=0
            occ_count=0
            for j in self.sequence.keys():
                if self.sequence[j][i]!='-':
                    occ_count+=1
                else:
                    blank_count +=1
            record.write('%s\t%s\t%s\n' % (i, blank_count, occ_count))

        return record

    def remove_all_blank_columns(self,same_length_check=True):
        num_seqs=len(self.sequence.values())
        seq_len = len(self.sequence.values()[0])

        if same_length_check==True:
            for i in self.sequence.values():
                assert len(i) == seq_len

        all_blanks_list = []
        for i in range(seq_len):
            allblank=True
            for j in self.sequence.values():
                if j[i]!='-':
                    allblank=False
                    break
            if allblank==True:
                all_blanks_list.append(i)

        if len(all_blanks_list)>0:
            all_blanks_list.sort(reverse=True)
            # for j in self.sequence.values():
            #     j = list(j)

            for i in all_blanks_list:
                for j in self.sequence.values():
                    j = j[0:i-1] + j[i:]
            #
            # for j in self.sequence.values():
            #     j = ''.join(j)




    def get_sequence_count(self):
        return len(self.sequence.keys())

    # def unmask_fasta_subset(self,fasta_dict,subset_keys):
    #     new_fasta={}
    #     vals=[]
    #     for i in subset_keys:
    #         new_fasta[i]=''
    #         vals.append(fasta_dict[i])
    #
    #     k=max(map(len,vals))
    #     for j in range(k):


def get_prefix(os=None):
    if os == None:
        os = platform.system()
    prefix_general_windows = 'C:/Users/Michael/Grad School Stuff/Research/Phylogenetics'
    prefix_general_linux = '/home/mikenute/Phylolab/share'
    if os == 'Windows':
        prefix_general = prefix_general_windows
    else:
        prefix_general = prefix_general_linux
    return prefix_general

def reduce_file(ref_file, in_file, reduced_file):
    # pref=get_prefix()
    # file1=pref+'/results/homfam/pasta_results/aln_full/homfam_aat__200.marker001.homfam_aat__raw.aln'
    file1=in_file
    # file_ref=pref+'/Data/homfam/aat/model/true.reduced.fasta'
    file_ref=ref_file
    # file_out=pref+'/results/homfam/pasta_results/aln_reduced/homfam_aat__200.marker001.homfam_aat.aln'
    file_out=reduced_file
    aln=alignment_from_fasta()
    fasta1=aln.read_from_fasta(file1)
    fasta_ref=aln.read_from_fasta(file_ref)

    a=set(fasta1.keys()).intersection(fasta_ref.keys())
    aln.write_to_fasta(file_out,fasta1,a)
    # for i in a:
    #     print i + ' - ' + str(len(fasta_ref[i])) + ' - ' + str(len(fasta1[i]))
    print ('done: ' + in_file)

def fastsp_get_empty_args():
    a= {
        'file_name':None,
        'sp':None,
        'modeler':None,
        'spfn':None,
        'spfp':None,
        'comp':None,
        'tc':None,
        'maxlen':None,
        'numseq':None,
        'lenref':None,
        'lenest':None,
        'cells':None
    }
    return a

if __name__ == '__main__':
    # zz_test_nams_fasta_files_for_blanks()
    myf=get_prefix()+'/code/sepp/test/unittest/data/upp/initial.fas'
    aln=alignment_from_fasta()
    aln.assign_fasta(myf)
    print (aln.get_sequence_count())
    # print aln.sequence.keys()







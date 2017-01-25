from startup import *
from command_line_utils import *
import re, platform, time, os
import numpy as np
import matplotlib.pyplot as plt
import dendropy
from Bio import SeqIO

def mask_fastadict(fasta, min_pct_nongap = 0.1):
    thresh = min_pct_nongap/100
    ntax = len(fasta.keys())
    ncols = len(fasta.values()[0])
    nparr = np.zeros((ntax,ncols),dtype=np.uint8)
    for i in range(ntax):
        seq = fasta[fasta.keys()[i]]
        nparr[i,:]=np.frombuffer(seq,np.uint8)

    # 45 is the uint8 code for the dash character '-':
    maskcols = np.where(np.sum((nparr<>45)*1,0).astype(np.float32)/nparr.shape[0]>thresh)
    newfasta = {}
    for i in range(ntax):
        k = fasta.keys()[i]
        newfasta[k] = str(np.getbuffer(nparr[i,maskcols]))
    return newfasta

def make_histogram_of_sequence_lengths(fasta_file,out_figure_path):
    a = read_from_fasta(fasta_file)
    lsnp = np.array(map(len,a.values()))
    plt.hist(lsnp,'auto')
    plt.xlabel('# Base Pairs')
    plt.ylabel('count')
    plt.title('Seq-Len Histogram: %s' % fasta_file)
    plt.savefig(out_figure_path)
    plt.clf()


def make_fasta_with_clean_names(fastadict,outfasta,outnames):
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
            if first<>True:
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
    L=[]
    for i in a.values():
        L.append(len(i))

    return L

def hist_of_seqs_by_name(seqs,i,fasta_dict):
    allkeys=[k for k in fasta_dict.keys() if k[0:10]==seqs[i]]
    hist={}
    for i in allkeys:
        l=len(fasta_dict[i])
        try:
            hist[l]+=1
        except:
            hist[l]=1
    return hist

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
    print filename
    seq = enumerate(SeqIO.parse(filename,"genbank"))
    ind, rec = seq.next()
    p,f = os.path.split(filename)
    cds = {}
    dna = {}
    st_en = {}
    no_translation_loci=[]
    length_errors=[]
    for i in rec.features:
        if i.type=='CDS':
            st = i.location.start + 0
            en = i.location.end + 0
            dir = i.location.strand

            if 'protein_id' in i.qualifiers.keys():
                id = 'pt:' + i.qualifiers['protein_id'][0]
            elif 'locus_tag' in i.qualifiers.keys():
                id = 'lt:' + i.qualifiers['locus_tag'][0]
            else:
                id = 'na:' + f + '_' + str(st) + '-' + str(en)
                # id = 'file-' + str()

            # print '%s, %s, %s, %s' %(str(st), str(en), str(dir), id)
            try:
                prot = i.qualifiers['translation'][0]
            except:
                no_translation_loci.append(id)
                continue
            if (en-st-3)/3<>len(prot):
                length_errors.append(id)
            if dir == -1:
                dnastr = rec.seq[st:en].reverse_complement()
            else:
                dnastr = rec.seq[st:en]
            cds[id] = prot
            dna[id] = dnastr
            st_en[id] = (st,en,dir)

    return cds, dna, st_en, no_translation_loci, length_errors, filename

    pass

def all_gbk_to_fna():
    from multiprocessing import Pool
    p=Pool(15)
    gbkl = get_list_from_file(kraken_file_list)
    p.map(bp_genbank_to_fasta,gbkl)
    # n = len(gbkl)
    # k=1
    # for i in gbkl:
    #     fnai = i[:-4] + '.fna'
    #     ct = bp_genbank_to_fasta(i,fnai)
    #     print "%i of %i\t-- %i lines\t-- %s" % (k,n,ct,i)
    #     k+=1


def fasta_dict_to_string(fasta_dict):
    a=''
    for i in fasta_dict.keys():
        a=a+'>'+i + '\n'
        a=a + fasta_dict[i] + '\n'
    return a

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

# def get_true_alignment_path(set):
#     pref=get_prefix()
#     return pref + '/Data/RoseDNA/100S1.tar/' + set + '/rose.aln.true.fasta'
#
# def get_true_tree_path(set):
#     pref=get_prefix()
#     return pref + '/Data/RoseDNA/100S1.tar/' + set + '/rose.tt'

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
    if optional_filename<>None:
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

def dna_to_protein(dna):
    # codon lookup table transcribed from wikipedia
    codon_lookup = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                    'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
                    'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                    'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
                    'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
                    'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                    'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
                    'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'(stop)','TGA':'(stop)',
                    'TAG':'(stop)'}
    k = int(len(dna)/3)
    extras = len(dna)-k*3
    if extras>0:
        print "The sequence given is not divisible by 3, there are %i extra nucleotides, which will be ignored" % extras
    str_protein = ''
    for i in range(k):
        codon = dna[(i*3):(i*3+3)].upper()
        str_protein = str_protein + codon_lookup[codon]
    str_protein = str_protein.replace('')
    return str_protein

def aligned_protein_to_nucleotides(prot,raw_dna):
    codon_lookup = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C',
                'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
                'CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','TTA':'L','TTG':'L','CTT':'L',
                'CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','ATG':'M','TTT':'F','TTC':'F',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'AGT':'S','AGC':'S','ACT':'T','ACC':'T','ACA':'T','ACG':'T','TGG':'W','TAT':'Y',
                'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'(stop)','TGA':'(stop)',
                'TAG':'(stop)'}
    #first check that they are the same:
    raw_prot=prot.replace('-','')
    prot_from_dna=dna_to_protein(raw_dna)
    prot_from_dna_minus_stop=prot_from_dna.replace('(stop)','')
    assert raw_prot==prot_from_dna_minus_stop, "Protein and DNA sequence do not match"
    # assert raw_prot==prot_from_dna[0:len(raw_prot)]




def get_list_from_file(filepath):
    myf=open(filepath,'r')
    ol=[]
    for i in myf:
        if i.strip()<>'':
            ol.append(i.strip())

    myf.close()
    return ol

def make_pbs_file(jobname='nute-job-id',tmstr="04:00:00",queue='secondary',nodes=None,filename="pbs-job-id.pbs",commandfile="",commandlist=None):
    if nodes==None:
        if queue=='secondary':
            nodes=12
        elif queue=='stat':
            nodes=16
        elif queue=='tallis':
            nodes=20
    if not isinstance(tmstr,str) and not isinstance(tmstr,unicode):
        tmstr=hours_to_wallclock_string(tmstr)

    cmdstr='''
#PBS -l walltime=XXXtimeXXX
#PBS -l nodes=1:ppn=XXXnodesXXX
#PBS -N XXXjobnameXXX
#PBS -q XXXqueueXXX
##PBS -j oe
###PBS -o XXXjobnameXXX.out
###PBS -e XXXjobnameXXX.err
#PBS -m be
#PBS -M mike.nute@gmail.com
#
#####################################

# Load Modules
module load java
module load python/2.7.8

# Change to the directory from which the batch job was submitted
cd $PBS_O_WORKDIR
    \n\n'''
    cmdstr=cmdstr.replace("XXXtimeXXX",str(tmstr))
    cmdstr=cmdstr.replace("XXXnodesXXX",str(nodes))
    cmdstr=cmdstr.replace("XXXqueueXXX",queue)
    cmdstr=cmdstr.replace("XXXjobnameXXX",jobname)

    # print commandfile

    if commandfile=="" or commandfile==None:
        if commandlist==None:
            for line in sys.stdin:
                cmdstr=cmdstr + line
        else:
            cmdstr = cmdstr + '\n'.join(commandlist)
    else:
        cf = open(commandfile,'r')
        for i in cf:
            cmdstr=cmdstr + i

    # cmdstr = cmdstr + commandString
    pbs=open(filename,'w')
    pbs.write(cmdstr)
    pbs.close()



def write_list_to_file(mylist,filepath):
    myf=open(filepath,'w')

    for i in mylist:
        myf.write(i + '\n')

    myf.close()

def get_two_trees(a,b):
    tax = dendropy.TaxonNamespace()

    tr1 = dendropy.Tree.get(path=a,schema='newick',rooting="force-unrooted",preserve_underscores=True,taxon_namespace=tax)
    tr2 = dendropy.Tree.get(path=b,schema='newick',rooting="force-unrooted",preserve_underscores=True,taxon_namespace=tax)
    return tr1, tr2, tax


def testbed():
    print 'this is a test'

def help():
    print '''
    Nute Phylogeny Utilities: Command Line Options:

    -pbs    Makes a pbs file according to the specifications
            given. Requires the following additional arguments:
            -jobname    name of the job
            -walltime   wall clock time for the job (hours)
            -queue      queue to submit to (default: secondary)
            -nodes      (default: based on queue)
                        number of nodes to request
            -pbsfile    path to the output pbs file
            -commands   file containing commands (optional,
                        default will take from stdin)
    '''


def parse_pbs_args():
    global jobname, tm
    jobname = sys.argv[sys.argv.index('-jobname') + 1]
    tm_ind = sys.argv.index('-walltime') + 1
    tm = float(sys.argv[tm_ind])
    timestr=hours_to_wallclock_string(tm)
    if '-queue' in sys.argv:
        queue_ind = sys.argv.index('-queue') + 1
        queue = sys.argv[queue_ind]
    else:
        queue = 'secondary'
    if '-nodes' in sys.argv:
        nodes_ind = sys.argv.index('-nodes') + 1
        nodes = sys.argv[nodes_ind]
    else:
        nodes = None
    outfile_ind = sys.argv.index('-pbsfile') + 1
    outfile = sys.argv[outfile_ind]
    if '-commands' in sys.argv:
        cmdfile = sys.argv[sys.argv.index('-commands') + 1]
    else:
        cmdfile=None
    # print 'command file:\t' + cmdfile
    # print '-commands' in sys.argv
    # print sys.argv.index('-commands')
    # print sys.argv[sys.argv.index('-commands') + 1]
    # print cmdfile
    return(jobname,timestr,queue,nodes,outfile,cmdfile)


def hours_to_wallclock_string(tm):
    import datetime
    h = int(tm)
    mf = 60 * (tm - h)
    m = int(mf)
    s = int(60 * (mf - m))
    timestr = str(h) + ":" + datetime.time(h, m, s).strftime('%M:%S')
    return timestr

def alignment_stats_dir_to_tabd(aln_dir,outpath):
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
        print filepath
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

def get_dict_from_tab_delimited_file(filename, keysfirst = True):
    myf = open(filename,'r')
    args = {}
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
    myf.close()
    return args

if __name__=='__main__':
    if '-h' in sys.argv:
        help()
    elif '-pbs' in sys.argv:
        pbsargs=parse_pbs_args()
        make_pbs_file(*pbsargs)
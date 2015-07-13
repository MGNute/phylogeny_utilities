__author__ = 'Michael'

import subprocess, sys
import multiprocessing
# from alignment_utils import get_prefix

def main(myfrom=None, myto=None, suffix='', inputlistfilename=None, proc=16):
    # a=read_list_from_file(get_prefix() + '/results/2015_04_20_baliphy_in_UPP/input_fasta_list.txt')
    if inputlistfilename==None:
        a=read_list_from_file()
    else:
        a=read_list_from_file(inputlistfilename)
    # print myfrom
    # print len(a[myfrom:myto])
    # print multiprocessing.cpu_count()

    p=multiprocessing.Pool(proc)
    args=zip(a[myfrom:myto],[suffix]*len(a[myfrom:myto]))
    p.map(run_baliphy_intermediate,args)
    pass

def run_baliphy_intermediate(args):
    run_baliphy(*args)

def run_baliphy(set, suffix=''):
    a=set.split('/')[1:4]
    # output_path='/scratch/users/nute2/results/baliphy_upp/' + '_'.join(a)+ '_' + suffix
    # print output_path
    output_path = '/projects/sciteam/jtr/work/nute/baliphy_upp/output/' + '_'.join(a)

    arg1='bali-phy'
    # arg2='~/data/rosedna/100S1/' + set + '.fasta'
    arg2=set
    # arg3='--config'
    # arg4='100S1_' + set + '_raw.txt'
    arg5= '--builtins-path'
    # arg6='/home/nute2/Phylolab/bali-phy/bali-phy-2.3.5/lib/bali-phy/'
    arg6 = '/projects/sciteam/jtr/programs/bali-phy/bali-phy-2.3.6/lib/bali-phy/'
    arg7='--modules-path'
    # arg8='/home/nute2/Phylolab/bali-phy/bali-phy-2.3.5/lib/bali-phy/modules/'
    arg8='/projects/sciteam/jtr/programs/bali-phy/bali-phy-2.3.6/lib/bali-phy/modules/'
    arg9='--name'
    arg10 = output_path
    # subprocess.call([arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10])
    # subprocess.call([arg1, arg2, arg5, arg6, arg7, arg8, arg9, arg10])
    subprocess.call([arg1, arg2, arg7, arg8, arg9, arg10])

def read_list_from_file(filename='/home/nute2/baliphy_upp_work/input_fasta_list.txt'):
    mylist = open(filename,'r')
    outlist=[]
    for i in mylist:
        outlist.append(i.strip())

    return outlist

def print_help():
    print """
    baliphy_raw_dispatcher.py dispatches baliphy to run baliphy on a set of files from a list.
    USAGE: python baliphy_raw_dispatcher [from] [to] (additional arguments)
        [from]: the entry of the list to start AFTER (so if you want to run baliphy on entries 17 to 32, [from] should be 16).
        [to]:   the entry of the list to stop ON (so in the previous example, this should be 32).

    Additional Optional Arguments:
        --suffix [suffix]:      a suffix to append to the output folders
        --filename [filename]:  the path to a specific list of inputs to use (defaults to 'input_fasta_list.txt')
        --proc [proc]:          number of processes to spawn. Ideally this should match the number of processors available.
        --help:     displays the help section (this text)
    """

if __name__=='__main__':
    if '--help' not in sys.argv:
        myfrom = int(sys.argv[1].strip())
        myto = int(sys.argv[2].strip())
        if '--suffix' in sys.argv:
            suffix = sys.argv[sys.argv.index('--suffix')+1]
        else:
            suffix= ''

        if '--filename' in sys.argv:
            filename = sys.argv[sys.argv.index('--filename')+1]

        if '--proc' in sys.argv:
            proc=int(sys.argv[sys.argv.index('--proc')+1])
        else:
            proc=16

        # print myfrom
        # print myto
        main(myfrom, myto, suffix,filename,proc)
    else:
        print_help()
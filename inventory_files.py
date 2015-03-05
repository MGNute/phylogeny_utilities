#-------------------------------------------------------------------------------
# Name:        inventory_files.py
# Purpose:     inventory all the results files from a large batch of pasta runs
#
# Author:      Michael
#
# Created:     28/01/2015
# Copyright:   (c) Michael 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

"""
Note that the variable 'resultspath' below must be set to give the path to the folder
jsut above where all the pasta-output folders are, and the argument is the folder to
inventory for pasta results. For example, in the case below, either "rosedna" or "rnasim"
could be passed as the first argument from the command line, depending on which one we 
want to list the outputs of.

\scratch    <--resultspath should point to this folder.
    \rosedna
        \rosedna_pastajob1
            ...
        \rosedna_pastajob2
            ...
        ...
    \rnasim
        \rnasim_pastajob1
            ...
        \rnasim_pastajob2
            ...
        ...
"""
import os, sys, re

def main(folder=None):
    rtime_match=re.compile('PASTA INFO: Total time spent: (?P<rtime>\d+\.\d+)[a-z]')

    resultspath='...SET ME...' + folder
    outfile=resultspath+'/output_files.txt'
    outf=open(outfile,'w')

    re_aln=re.compile('.marker001.')
    re_tree0=re.compile('_temp_iteration_0_tree')
    re_tree1=re.compile('_temp_iteration_1_tree')
    re_tree2=re.compile('_temp_iteration_2_tree')
    re_tree_initial=re.compile('_temp_iteration_initialsearch_')
    for rt, pa, fi in os.walk(resultspath):
        aln_name=''
        tree0=''
        tree1=''
        tree2=''
        tree_final=''
        for i in fi:
            if i[-8:]=='.out.txt':
##                print rt
                outf.write(rt + '\t' + i + '\t')
                print i
                myfile=open(rt+'/'+i)
                mystr=myfile.read()
                b=rtime_match.search(mystr)
                myfile.close()
                if b <> None:
                    outf.write(b.group('rtime'))
                outf.write('\t')
            if re_aln.search(i)<>None:
                aln_name=i
            elif re_tree0.search(i)<>None:
                tree0=i
            elif re_tree1.search(i)<>None:
                tree1=i
            elif re_tree2.search(i)<>None:
                tree2=i
            elif i[-4:]=='.tre' and re_tree_initial.search(i)==None:
                tree_final=i
        if aln_name<>'':
            outf.write(aln_name + '\t' + tree0 + '\t' + tree1 +'\t' + tree2 + '\t' + tree_final + '\n')
    outf.close()

if __name__ == '__main__':
    folder = sys.argv[1]
    main(folder)

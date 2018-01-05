import sys, collections
import utilities

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
    for l in fasta:
        if l[0]=='>':
            if first<>True:
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

def rename_fasta_to_numbers(fasta_path, name_map_path, new_fasta_path, name_prefix = ''):
    a = read_from_fasta(fasta_path)
    if name_prefix <> '':
        name_prefix = name_prefix + '_'
    nmap= open(name_map_path,'w')
    newf={}
    ct = 1
    for i in a.keys():
        newname = name_prefix + str(ct)
        newf[newname]=a[i]
        nmap.write(i + '\t' + newname + '\n')
        ct+=1

    nmap.close()
    write_to_fasta(new_fasta_path,newf)

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
            print 'There were ' + str(len(leftover_keys)) + ' keys in the subset that were not in the original data.\n'


    fasta=open(out_file_path,'w')
    for i in mykeys:
        fasta.write('>'+i+'\n')
        if raw==False:
            fasta.write(fasta_dict[i]+'\n')
        else:
            fasta.write(fasta_dict[i].replace('-','')+'\n')
    fasta.close()
    if not quiet:
        print 'wrote file: ' + out_file_path + ' with the specified keys.'


def remove_all_blank_columns(fasta_dict,same_length_check=True):
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
    num_seqs=len(fasta_dict.values())
    seq_len = len(fasta_dict.values()[0])

    # check that all the sequences are the same length
    if same_length_check==True:
        for i in fasta_dict.values():
            if len(i) <> seq_len:
                print 'The sequences were not all the same length.'
                return fasta_dict, -1

    # identify columns that are blank for every taxon
    all_blanks_list = []
    for i in range(seq_len):
        allblank=True
        for j in fasta_dict.values():
            if j[i]<>'-':
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

    return newfasta, all_blanks_list

def rename_fasta_by_name_map(inalign,outalign,namemap,keysfirst = True):
    nm = utilities.get_dict_from_tab_delimited_file(namemap,keysfirst)
    initalgn = read_from_fasta(inalign)
    outfasta = {}
    for i in initalgn.keys():
        nk = nm[i]
        outfasta[nk]=initalgn[i]
    write_to_fasta(outalign,outfasta)

def rename_tree_by_name_map(intree, outtree, namemap,keysfirst=True):
    import dendropy
    tr = dendropy.Tree.get(path=intree, schema="newick",preserve_underscores=True)
    nm = utilities.get_dict_from_tab_delimited_file(namemap, keysfirst)
    for i in tr.leaf_nodes():
        a=i.taxon.label
        i.taxon.label=unicode(nm[a])
    tr.write(path=outtree, schema="newick")
    pass

def shrink_one_fasta_to_match_another(big_fasta_path, smaller_fasta_path, out_path, small_fasta_is_list=False):
    """
    Takes one fasta file, removes all the taxa except those that are in another fasta file, removes resulting
    all-blank columns, and writes the new fasta to a new file.
    :param big_fasta_path: the larger (source) fasta file path
    :param smaller_fasta_path: the smaller (target) fasta file path
    :param out_path: the path of the output file (!!!opens as writeable without warning if it's overwriting
        another file!!!)
    :return:
    """
    big = read_from_fasta(big_fasta_path)
    out_dict = {}
    if small_fasta_is_list==False:
        small = read_from_fasta(smaller_fasta_path)
        new_keys = list(set(big.keys()).intersection(set(small.keys())))
    else:
        new_keys = utilities.get_list_from_file(smaller_fasta_path)

    for i in new_keys:
        out_dict[i]=big[i]

    out_dict, all_blanks_list = remove_all_blank_columns(out_dict)
    write_to_fasta(out_path,out_dict)

def split_fasta_into_parts(in_path, num_parts, out_folder = None):
    import math, os
    fa = read_from_fasta(in_path)

    if out_folder==None:
        path_template = in_path + '.part%s'
    else:
        fo, fi = os.path.split(in_path)
        path_template = os.path.join(out_folder,fi) + '.part%s'

    cuts = [i*len(fa)/num_parts for i in range(num_parts+1)]
    for i in range(num_parts):
        c1 = cuts[i]
        c2 = cuts[i+1]
        fanew = dict(fa.items()[c1:c2])
        new_path = path_template % i
        write_to_fasta(new_path,fanew)
        print 'wrote file: %s' % new_path
        del fanew

def split_fastq_into_fasta_and_quality(in_path, out_fasta_p, out_quality_p):
    myf = open(in_path,'r')
    out_fasta = open(out_fasta_p,'w')
    out_quality = open(out_quality_p,'w')

    ct = 0
    ln4 = []
    reads = 0

    for line in myf:
        ln4.append(line.strip())
        ct +=1
        if ct % 4 == 0:
            # fix the name
            nm = ln4[0].replace('@','')
            nm = nm.replace(':','_')
            nm = nm.replace(' ','--')

            assert ln4[2]=='+'
            out_fasta.write('>' + nm + '\n')
            out_quality.write('>' + nm + '\n')
            out_fasta.write(ln4[1] + '\n')
            out_quality.write(ln4[3] + '\n')

            ln4 = []
            reads += 1

    myf.close()
    out_fasta.close()
    out_quality.close()
    print "Done. Successfully converted fastq file %s into fasta file %s for sequences and %s for quality scores" % \
          (in_fastq, out_fasta_p, out_quality_p)

def dedupe_and_rename(in_aln, name_map, out_aln, out_mult_file=None):
    ia = read_from_fasta(in_aln)
    if out_mult is not None:
        mults = dict(collections.Counter(ia.values()))
    print "de-duping %s" % in_aln
    nm = {}
    # back_lkp = {}
    # fwd_lkp = {}
    a_len = len(ia.values())
    ia_seq_uq = list(set(ia.values()))
    a_len_uq = len(ia_seq_uq)
    print '\toriginal file has %s sequences and %s uniques' % (a_len, a_len_uq)
    back_lkp_it = []
    fwd_lkp_it = []

    print '\tmaking iterables'
    ct = 0
    done = 0
    for v in ia_seq_uq:
        newname = 's' + str(ct)
        back_lkp_it.append((v,newname))
        fwd_lkp_it.append((newname,v))
        ct += 1

    print '\tmaking dictionaries'
    fwd_lkp = dict(fwd_lkp_it)
    bck_lkp = dict(back_lkp_it)


    if out_mult is not None:
        out_mult_dict = dict.fromkeys(fwd_lkp.keys())
        for k in out_mult_dict.keys():
            out_mult_dict[k]=mults[fwd_lkp[k]]
        write_dict_to_file(out_mult_dict,out_mult_file)


    # a = len(ia.keys())
    # b = len(fwd_lkp.keys())
    # c = a-b
    # print "\tthe original had %s sequences, new has %s sequences. %s dupes" % (a,b,c)
    print '\twriting to fasta'
    write_to_fasta(out_aln,fwd_lkp)

    print '\tmaking name_map'
    myf = open(name_map,'w')
    mapct = 0
    for k in ia.keys():
        myf.write(k + '\t')
        myf.write(bck_lkp[ia[k]] + '\n')
        mapct +=1
    myf.close()
    print '\twrote %s values to the name map file %s' % (mapct, name_map)

def make_rc_file(fasta,out_path):
    fas = read_from_fasta(fasta)
    fasrc = utilities.get_fastadict_reverse_complement(fas)
    write_to_fasta(out_path,fasrc)

def remove_blanks_from_file(fasta,out_path):
    fas=read_from_fasta(fasta)
    fas2,ab = remove_all_blank_columns(fas)
    write_to_fasta(out_path,fas2)

if __name__=='__main__':
    if '-h' in sys.argv or '--help' in sys.argv:
        print '''
usage: python command_line_utils.py [...options...]

Options:
    --shrink-to-fit -b [big_fasta] -s [small_fasta] -o [out_path]
        --->Shrinks one fasta file to include only the taxa in the other file
    --shrink-alignment-to-fit-list -b [big_fasta] -s [small_list] -o [out_path]
        --->same as previous but assuming the smaller file is a list of taxon names instead of a fasta
    --extract-complement -b [big_fasta] -s [small_fasta] -o [out_path]
        --->Shrinks one fasta file to include only the taxa NOT in the other
    --remove-blank-columns -f [fasta] OPTIONAL: -o [out_path]
        --->Removes all of the columns from a fasta file that are exclusively blank
    --FastSPfold -f [folder_path] -o [output_file_path]
        --->Converts a folder full of FastSP results to a tab delimited file
    --AlignmentStatsFold -f [folder_path] -o [output_file_path]
        --->Converts a folder full of Alignment Statistics results to a tab delimited file
    --SanitizeNames -f [fasta_path] -o [output_fasta_path] -m [name_map_path] -p [prefix (opt)]
        --->Converts a fasta to another fasta with all taxa renamed as 'prefix_#'. Conversions are
            output to a tab-delimited file name_map_path
    --RenameAlignment -ia [in_alignment] -nm [name_map] -oa [out_alignment] OPTIONAL: -ks
        --->Converts a fasta to another fasta with all taxa renamed according to a tab-delimited name
            map.
            if -ks (keys-second) is given, name map file is assumed to be "value<tab>key", otherwise it is assumed
            to be "key<value>value" where 'key' is the name from the first fasta.
    --RenameTree -it [in_tree] -nm [name_map] -ot [out_tree] OPTIONAL: -ks
        --->Converts a Newick tree to another Newick tree with all taxa renamed according to a tab-
            delimited name-map.
            if -ks (keys-second) is given, name map file is assumed to be "value<tab>key", otherwise it is assumed
            to be "key<value>value" where 'key' is the name from the first tree.
    --SplitFasta -i [input_file] -p [# parts] OPTIONAL: -o [dump_folder]
        --->Cuts a file into '-p' parts. Each part has the identical path and filename, but iwth ".part#" appended.
            if '-o' is given, results will be written to that folder (with same names) instead of original.
    --FastqToFasta -i [input_file] -of [out_fasta] -oq [out_quality]
        --->Takes a Fastq file and parses it, writing it to two different files for the sequences and their quality
            scores. Each is (effectively) a fasta-formatted file. Also does some string cleanup on the sequence names t0
            use only dashes and underscores.
    --DedupeAndRename -ia [in_alignment] -nm [name_map] -oa [out_alignment] OPTIONAL: -mult [multiplicity_file]
        --->Takes a fasta file and removes all duplicate sequences from the file. Additionally renames the sequences
            with a sequenectial "s###" scheme and prints a name map to the file [name-map]. Optionally, if [multi-
            plicity_file] is given, a reference is written as a tab-delimited file with the new squence name and a
            count of occurrences in the input alignment.
    --MakeLengthHistogram -in [in_fasta] -out [output_image_path]
        --->Makes a histogram of the sequence lengths for the sequences in the input fasta file. Saves the file
            in a format that is dependent on the path, so best to end it with '.pdf'
    --MakeReverseComplement -in [in_fasta] -out [output_fasta]
        --->Makes a version of the same fasta with all the same sequence names, but with every string as its reverse
            complement.
        '''
        sys.exit(0)
    if sys.argv[1]=='--shrink-to-fit':
        big_fasta = sys.argv[sys.argv.index('-b')+1]
        small_fasta = sys.argv[sys.argv.index('-s')+1]
        out_path = sys.argv[sys.argv.index('-o')+1]
        print big_fasta
        print small_fasta
        print out_path
        shrink_one_fasta_to_match_another(big_fasta,small_fasta,out_path)
    elif sys.argv[1]=='--shrink-alignment-to-fit-list':
        big_fasta = sys.argv[sys.argv.index('-b')+1]
        small_fasta = sys.argv[sys.argv.index('-s')+1]
        out_path = sys.argv[sys.argv.index('-o')+1]
        print big_fasta
        print small_fasta
        print out_path
        shrink_one_fasta_to_match_another(big_fasta,small_fasta,out_path,True)
    elif sys.argv[1]=='--extract-complement':
        big_fasta = sys.argv[sys.argv.index('-b')+1]
        small_fasta = sys.argv[sys.argv.index('-s')+1]
        out_path = sys.argv[sys.argv.index('-o')+1]
        import utilities
        utilities.shrink_fasta_to_complement_of_another(big_fasta,small_fasta,out_path)
    elif sys.argv[1]=='--remove-blank-columns':
        fasta = sys.argv[sys.argv.index('-f')+1]
        if '-o' in sys.argv:
            out_path = sys.argv[sys.argv.index('-o')+1]
            print "out_path: %s" % out_path
        else:
            out_path = fasta
        remove_blanks_from_file(fasta,out_path)
    elif sys.argv[1]=='--FastSPfold':
        fold=sys.argv[sys.argv.index('-f')+1]
        outfile=sys.argv[sys.argv.index('-o')+1]
        print "converting FastSP folder %s to tab-delimited file %s" % (fold,outfile)
        # from phylogeny_utilities.alignment_utils import *
        from alignment_utils import *
        fastsp_results_folder_to_tab_delimited(fold,outfile)
    elif sys.argv[1]=='--AlignmentStatsFold':
        fold=sys.argv[sys.argv.index('-f')+1]
        outfile=sys.argv[sys.argv.index('-o')+1]
        print "converting Alignment Stats folder %s to tab-delimited file %s" % (fold,outfile)
        # from phylogeny_utilities.utilities import *
        from utilities import *
        alignment_stats_dir_to_tabd(fold,outfile)
    elif sys.argv[1]=='--SanitizeNames':
        infile=sys.argv[sys.argv.index('-f')+1]
        outfile=sys.argv[sys.argv.index('-o')+1]
        name_map = sys.argv[sys.argv.index('-m')+1]
        if '-p' in sys.argv:
            prefix = sys.argv[sys.argv.index('-p')+1]
        else:
            prefix = ''
        print "converting fasta %s to %s with sanitized names" % (infile,outfile)
        rename_fasta_to_numbers(infile,name_map,outfile,prefix)
    elif sys.argv[1]=='--RenameAlignment':
        # intree=sys.argv[sys.argv.index('-it')+1]
        inaln = sys.argv[sys.argv.index('-ia') + 1]
        # outtree = sys.argv[sys.argv.index('-ot') + 1]
        outaln = sys.argv[sys.argv.index('-oa') + 1]
        name_map = sys.argv[sys.argv.index('-nm') + 1]
        if '-ks' in sys.argv:
            rename_fasta_by_name_map(inaln, outaln, name_map, False)
        else:
            rename_fasta_by_name_map(inaln, outaln, name_map, True)
    elif sys.argv[1]=='--RenameTree':
        intree=sys.argv[sys.argv.index('-it')+1]
        # inaln = sys.argv[sys.argv.index('-ia') + 1]
        outtree = sys.argv[sys.argv.index('-ot') + 1]
        # outaln = sys.argv[sys.argv.index('-oa') + 1]
        name_map = sys.argv[sys.argv.index('-nm') + 1]
        if '-ks' in sys.argv:
            rename_tree_by_name_map(intree,outtree, name_map,False)
        else:
            rename_tree_by_name_map(intree, outtree, name_map, True)
    elif sys.argv[1]=='--SplitFasta':
        infasta = sys.argv[sys.argv.index('-i')+1]
        numparts = int(sys.argv[sys.argv.index('-p')+1])
        if '-o' in sys.argv:
            out_folder = sys.argv[sys.argv.index('-o')+1]
            split_fasta_into_parts(infasta,numparts,out_folder)
        else:
            split_fasta_into_parts(infasta,numparts)
    elif sys.argv[1]=='--FastqToFasta':
        in_fastq = sys.argv[sys.argv.index('-i')+1]
        out_fasta = sys.argv[sys.argv.index('-of') + 1]
        out_quality = sys.argv[sys.argv.index('-oq') + 1]
        split_fastq_into_fasta_and_quality(in_fastq, out_fasta, out_quality)
    elif sys.argv[1] == '--DedupeAndRename':
        inaln = sys.argv[sys.argv.index('-ia') + 1]
        outaln = sys.argv[sys.argv.index('-oa') + 1]
        name_map = sys.argv[sys.argv.index('-nm') + 1]
        if '-mult' not in sys.argv:
            dedupe_and_rename(inaln, name_map, outaln)
        else:
            multfile = sys.argv[sys.argv.index('-mult') + 1]
            dedupe_and_rename(inaln, name_map, outaln,multfile)
    elif sys.argv[1] == '--MakeLengthHistogram':
        inaln = sys.argv[sys.argv.index('-in') + 1]
        outpath = sys.argv[sys.argv.index('-out') + 1]
        make_histogram_of_sequence_lengths(inaln, outpath)
    elif sys.argv[1] == '--MakeReverseComplement':
        inaln = sys.argv[sys.argv.index('-in') + 1]
        outpath = sys.argv[sys.argv.index('-out') + 1]
        make_rc_file(inaln,outpath)

    else:
        print "No major option recognized. Check your \'--\' option and try again."
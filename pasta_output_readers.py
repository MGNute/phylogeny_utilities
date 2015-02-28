#-------------------------------------------------------------------------------
# Name:        pasta_output_readers
# Purpose:
#
# Author:      Michael
#
# Created:     10/02/2015
# Copyright:   (c) Michael 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import os, re, ConfigParser, platform, datetime

class pasta_job:
    """
    member properties:
        job_name        - used for pasta & pbs job name, etc...
        alignment_file  - path to file containing true alignment
        tree_file       - path to file containing true truee
        local_fasta_folder  - folder containing the raw .fasta files
        data_set_name   - name of the data set, used in various places (e.g. rosedna_MS)
        results_folder  - general folder hodling the results.
        pbs_folder      - folder to warehouse the pbs files until they are transferred
        cc_data_folder  - folder containing data files on the campus cluster

    class to handle the following tasks for each pasta_run:
        - identifying the location of the true tree and alignment file
        - creating the "raw.fasta" file for PASTA input
        - creating the pbs file to submit the job
        - creating the sftp command to transfer the raw files
    """

    def __init__(self,alignment_file=None,tree_file=None,job_name=None, subproblem_size=None, config_file_path=None):
        self.job_name=job_name
        self.alignment_file=alignment_file
        self.tree_file=tree_file
        self.subproblem_size=subproblem_size
        self.path_to_input_file=None
        self.prefix=None
        self.get_prefix()
        self.cfg=ConfigParser.ConfigParser()
        if config_file_path==None:
            self.cfg.read('pasta_batch_config.cfg')
        else:
            self.cfg.read(config_file_path)
        self.local_fasta_folder=self.cfg.get('pasta_batch_config','local_fasta_folder')
        self.data_set_name=self.cfg.get('pasta_batch_config','data_set_name')
        self.data_type=self.cfg.get('pasta_batch_config','data_type')
        self.results_folder=self.cfg.get('pasta_batch_config','results_folder')
        self.pbs_folder=self.cfg.get('pasta_batch_config','pbs_folder')
        self.cc_data_folder=self.cfg.get('pasta_batch_config','cc_data_folder')
        self.cc_pbs_folder=self.cfg.get('pasta_batch_config','cc_pbs_folder')
        self.cc_results_folder=self.cfg.get('pasta_batch_config','cc_results_folder')

    def get_prefix(self):
        if platform.system()=='Windows':
            self.prefix='C:/Users/Michael/Grad School Stuff/Research/Phylogenetics'
        elif platform.system()=='Linux':
            self.prefix='/home/mikenute/Phylolab/share'

    def set_alignment_file(self,alignment_file=None):
        self.alignment_file=alignment_file

    def set_tree_file(self,tree_file=None):
        self.tree_file=tree_file

    def set_all_parameters(self,alignment_file,tree_file, job_name, subproblem_size):
        self.alignment_file=alignment_file
        self.tree_file=tree_file
        self.job_name=job_name
        self.subproblem_size=subproblem_size

    def make_pbs_file(self, queue='secondary',append_to=None):
        # path_to_input_file=self.cc_data_folder + self.job_name+'_raw.fasta'
        if queue == 'secondary':
            outf=open(self.prefix + self.pbs_folder + self.job_name + '.pbs','w')
            tmpf=open(self.prefix + '/code/phylogeny_utilities/pasta_pbs_template.pbs','r')
            for line in tmpf:
                a=line.replace('PATH_TO_INPUT_FILE',self.path_to_input_file)
                a=a.replace('JOB_NAME',self.job_name)
                a=a.replace('PATH_TO_RESULTS_FOLDER',self.cc_results_folder)
                a=a.replace('TIMES_FILE_NAME',self.job_name + '_times.log')
                a=a.replace('DATA_TYPE',self.data_type)
                a=a.replace('SFTP_COMMAND_FILE',self.cfg.get('pasta_batch_config','cc_sftp_folder')+self.job_name + '_sftp.txt')
                a=a.replace('PASTA_TEMP_FOLDER',self.cfg.get('pasta_batch_config','cc_pasta_temp_folder') + self.job_name)
                a=a.replace('MAX_SUBPROBLEM_SIZE',str(self.subproblem_size))
                outf.write(a)
            outf.close()
            tmpf.close()
            return '# qsub '+self.cc_pbs_folder + self.job_name + '.pbs\n'
        if queue == 'stat':
            if append_to==None:
                outf=open(self.prefix + self.pbs_folder + self.job_name + '.pbs','w')
                tmpf=open('pasta_pbs_template.pbs','r')
                for line in tmpf:
                    a=line.replace('PATH_TO_INPUT_FILE',self.path_to_input_file)
                    a=a.replace('JOB_NAME',self.job_name)
                    a=a.replace('PATH_TO_RESULTS_FOLDER',self.cc_results_folder)
                    a=a.replace('TIMES_FILE_NAME',self.job_name + '_times.log')
                    a=a.replace('DATA_TYPE',self.data_type)
                    a=a.replace('SFTP_COMMAND_FILE',self.cfg.get('pasta_batch_config','cc_sftp_folder')+self.job_name + '_sftp.txt')
                    a=a.replace('PASTA_TEMP_FOLDER',self.cfg.get('pasta_batch_config','cc_pasta_temp_folder') + self.job_name)
                    a=a.replace('MAX_SUBPROBLEM_SIZE',str(self.subproblem_size))
                    a=a.replace('#PBS -q secondary','#PBS -q stat')
                    a=a.replace('#PBS -l walltime=4:00:00','#PBS -l walltime=48:00:00')
                    outf.write(a)
                outf.close()
                tmpf.close()
                return '# qsub '+self.cc_pbs_folder + self.job_name + '.pbs\n'
            else:
                outf=open(self.prefix + self.pbs_folder + append_to + '.pbs','a')
                tmpf=open('pasta_pbs_main_commands.txt','r')
                for line in tmpf:
                    a=line.replace('PATH_TO_INPUT_FILE',self.path_to_input_file)
                    a=a.replace('JOB_NAME',self.job_name)
                    a=a.replace('PATH_TO_RESULTS_FOLDER',self.cc_results_folder)
                    a=a.replace('TIMES_FILE_NAME',self.job_name + '_times.log')
                    a=a.replace('DATA_TYPE',self.data_type)
                    a=a.replace('SFTP_COMMAND_FILE',self.cfg.get('pasta_batch_config','cc_sftp_folder')+self.job_name + '_sftp.txt')
                    a=a.replace('PASTA_TEMP_FOLDER',self.cfg.get('pasta_batch_config','cc_pasta_temp_folder') + self.job_name)
                    a=a.replace('MAX_SUBPROBLEM_SIZE',str(self.subproblem_size))
                    outf.write(a)
                outf.close()
                tmpf.close()
                return ''

    def make_raw_fasta(self, file_name=None,override=False):
        if file_name==None:
            file_name_full = self.prefix +self.local_fasta_folder + self.job_name + '_raw.fasta'
            file_name = self.job_name
        else:
            file_name_full = self.prefix +self.local_fasta_folder + file_name + '_raw.fasta'

        self.path_to_input_file=self.cc_data_folder + file_name + '_raw.fasta'
        if os.path.isfile(file_name_full)==False or override==True:
            trueF = open(self.alignment_file,'r')
            rawF = open(file_name_full,'w')
            i=0
            for l in trueF:
                rawF.write(l.replace('-',''))
            trueF.close()
            rawF.close()

    def make_sftp_file(self,garbage=False):
        if garbage==True:
            myfile=self.prefix + self.cfg.get('pasta_batch_config','sftp_folder_garbage') + self.job_name + '_sftp.txt'
        else:
            myfile=self.prefix + self.cfg.get('pasta_batch_config','sftp_folder') + self.job_name + '_sftp.txt'
        outf=open(myfile,'w')
        dumploc=self.cfg.get('pasta_batch_config','cincy_dump_folder')
        outf.write('cd ' + dumploc + '\n')
        if garbage==False:
            outf.write('mkdir '+self.job_name + '\n')

        outf.write('cd '+self.job_name + '\n')
        outf.write('mput ' + self.cc_results_folder + self.job_name + '/*\n')
        outf.close()
        pass

    def make_job_list_tuple(self):
        out=(self.data_set_name, self.alignment_file, self.tree_file, self.job_name, self.subproblem_size)
        return out

class pastaRunsPreparer:
    def __init__(self,prefix_data=None):
        self.aln_tree_files_true=[]
        self.prefix_general_windows='C:/Users/Michael/Grad School Stuff/Research/Phylogenetics'
        self.prefix_general_linux='/home/mikenute/Phylolab/share'
        self.get_prefix_general()
        self.prefix_data=prefix_data

    def set_prefix_data(self,prefix_data):
        self.prefix_data=prefix_data

    def get_prefix_general(self):
        if platform.system()=='Windows':
            self.prefix_general=self.prefix_general_windows
        else:
            self.prefix_general=self.prefix_general_linux

    def walk_data_folder_for_fastas(self,regex_aln=None, regex_tree=None, regex_exclude=None):
##        folder=folder.replace(self.prefix_general,'')
##        folder=folder.replace(self.prefix_data,'')
        folder=self.prefix_general + self.prefix_data

        if regex_tree<>None:
            reg1=re.compile(regex_tree)
        else:
            reg1=None
        if regex_exclude<>None:
            reg2=re.compile(regex_exclude)
        else:
            excl=None
            reg2=None
        if regex_aln<>None:
            reg3=re.compile(regex_aln)

        for (i,j,k) in os.walk(folder):
            aln_tree_pair=None
            aln=''
            tree=''
            if reg2<>None:
                excl=reg2.search(i)
            if excl==None:
                for myfi in k:
                    if reg1<>None:
                        incl=reg1.search(myfi)
                        if incl<>None:
                            path=i+'/'+myfi
                            path=path.replace('\\','/')
                            tree=path
                    if reg3<>None:
                        incl2=reg3.search(myfi)
                        if incl2<>None:
                            path=i+'/'+myfi
                            path=path.replace('\\','/')
                            aln=path
            if aln<>'':
                aln_tree_pair=(aln,tree)
                self.aln_tree_files_true.append(aln_tree_pair)

def main2(cfg):
    prefix=get_prefix()
    # inputs - general
    prefix_data=cfg.get('pasta_batch_config','local_data_search_folder')
    runsf=open(prefix + cfg.get('pasta_batch_config','results_folder') + 'pasta_job_list.txt','w')
    aln_str=cfg.get('data_strings','aln_str')
    tree_str=cfg.get('data_strings','tree_str')
    excl_str=cfg.get('data_strings','exclude_str')
    reg1=re.compile(cfg.get('data_strings','reg1_expr'))
    try:
        reg2=re.compile(cfg.get('data_strings','reg2_expr'))
    except:
        reg2=None
    batchsubf=open(prefix + cfg.get('pasta_batch_config','results_folder') +'batch_submit_'+cfg.get('pasta_batch_config','data_set_name')+'.sh','w')
    jobname_prefix=cfg.get('pasta_batch_config','data_set_name')+'_'
    config_file_loc=cfg.get('pasta_batch_config','results_folder')+'pasta_batch_config.cfg'

    # intputs for Rosedna
##    prefix_data='/Data/rosedna/'
##    runsf=open(prefix + '/results/rosedna_MS/pasta_job_list.txt','w')
##    aln_str='rose.aln.true.fasta'
##    tree_str='rose.tt'
##    excl_str='1000L'
##    reg1=re.compile('/1000(?P<model>[MS][1234]).tar/')
##    reg2=re.compile('/(?P<sample>R\d+)/')
##    batchsubf=open(prefix + '/results/rosedna_MS/batch_submit_rosedna_MS.sh','w')
##    jobname_prefix='rose_'
##    config_file_loc='/results/rosedna_MS/pasta_batch_config.cfg'

    # intputs for homfam
##    prefix_data='/Data/homfam/'
##    runsf=open(prefix + '/results/homfam/pasta_job_list.txt','w')
##    aln_str='initial.fas'
##    tree_str=None
##    excl_str=None
##    reg1=re.compile('/homfam/(?P<model>.*?)/model/')
####    reg2=re.compile('/(?P<sample>R\d+)/')
##    samp=''
##    batchsubf=open(prefix + '/results/homfam/batch_submit_homfam.sh','w')
##    jobname_prefix='homfam_'
##    config_file_loc='/results/homfam/pasta_batch_config.cfg'

    # inputs for Indelible
##    prefix_data='/Data/indelible/'
##    runsf=open(prefix + '/results/indelible_v2/pasta_job_list.txt','w')
##    aln_str='true.fasta$'
##    tree_str='true.tt'
##    excl_str='homfam'
##    reg1=re.compile('/10000(?P<model>[MS][1234])/')
##    reg2=re.compile('/(?P<sample>\d+)/model/')
##    batchsubf=open(prefix + '/results/indelible_v2/batch_submit_indelible_v2_extras.sh','w')
##    jobname_prefix='indel2_'
##    config_file_loc='/results/indelible_v2/pasta_batch_config.cfg'

    # create new pastas object
    newpastas=pastaRunsPreparer(prefix_data)
    newpastas.walk_data_folder_for_fastas(aln_str,tree_str,excl_str)
    print len(newpastas.aln_tree_files_true)


    # logging
    runsf.write('data_set_name\talignment_file\ttree_file\tjob_name\tsubproblem_size\n')
    str_runs_list='%s\t%s\t%s\t%s\t%s\n'

    # config file for this batch
    mypastajob=pasta_job(config_file_path=prefix + config_file_loc)

    sizes=['200','100','50','25','10']

    for i in newpastas.aln_tree_files_true:
        try:
            a=reg1.search(i[0])
            mod=a.group('model')
            if reg2<>None:
                a=reg2.search(i[0])
                samp=a.group('sample')
            else:
                samp='1'
            file_name= jobname_prefix + mod + '_' + samp


            for j in sizes:
                job_name= jobname_prefix + mod+'_'+samp+'_' + j
                mypastajob.set_all_parameters(i[0],i[1],job_name,j)
                mypastajob.make_raw_fasta(file_name)
                qsub=mypastajob.make_pbs_file()

                #for indelible specifically:
    ##            if mod=='M2' and samp=='2':
    ##                qsub=mypastajob.make_pbs_file()
    ##            else:
    ##                append_to='indel2_M2_0_' + j
    ##                qsub=mypastajob.make_pbs_file('stat',append_to)
                mypastajob.make_sftp_file()
                out=mypastajob.make_job_list_tuple()

                runsf.write(str_runs_list % out)
                print job_name
                batchsubf.write(qsub)
        except:
            print 'no model or sample found for file: ' + i[0]
##        break
    sftp=open(prefix + cfg.get('pasta_batch_config','results_folder') + 'send_files_to_CC.sftp','w')
    sftp.write('lcd ' + prefix + cfg.get('pasta_batch_config','results_folder') + '\n')
    sftp.write('mput '+prefix + cfg.get('pasta_batch_config','results_folder') +'batch_submit_'+cfg.get('pasta_batch_config','data_set_name')+'.sh\n')
    sftp.write('cd '+ cfg.get('pasta_batch_config','cc_data_folder') +'\n')
    sftp.write('mput ' + prefix + cfg.get('pasta_batch_config','local_fasta_folder')+'*\n')
    sftp.write('cd '+ cfg.get('pasta_batch_config','cc_pbs_folder') +'\n')
    sftp.write('mput ' + prefix + cfg.get('pasta_batch_config','pbs_folder')+'*\n')
    sftp.write('cd '+ cfg.get('pasta_batch_config','cc_sftp_folder') +'\n')
    sftp.write('mput ' + prefix + cfg.get('pasta_batch_config','sftp_folder')+'*\n')
    sftp.close()


    runsf.close()

def collect_garbage():
    prefix=get_prefix()
    folder = '/results/indelible_v2/'
    cc_pbs_folder = '/home/nute2/pbs_secondarys/indelible_v2/'
    config_file_path = prefix+ folder + 'pasta_batch_config.cfg'
    pjob=pasta_job(config_file_path=config_file_path)

    # get completed files
    done_list=[]
    done_list_file=open(prefix + folder + 'output_files.txt','r')
    for line in done_list_file:
        a=line.split('\t')
        b=a[1].replace('.out.txt','')
        done_list.append(b)
    done_list_file.close()

    # get original files
    job_list = []
    job_list_dict={}
    job_list_file=open(prefix + folder + 'pasta_job_list.txt','r')
    for line in job_list_file:
        a=line.strip().split('\t')
        if a[3]<>'job_name':
            job_list.append(a[3])
            job_list_dict[a[3]]=a
    job_list_file.close()


    # new batch file
    bl=open(prefix + folder + 'batch_submit_garbage.sh','w')
    for i in job_list:
        if i not in done_list:
            bl.write('qsub '+ cc_pbs_folder + i + '.pbs\n')
            jobline=job_list_dict[i]
            size=i.split('_')[-1]
            pjob.set_all_parameters(jobline[1],jobline[2],jobline[3],size)
            pjob.make_sftp_file(True)

    bl.close()


##    alignment_file='/home/mikenute/Phylolab/share/Data\\RoseDNA\\1000M4.tar\\1000,1000,.00003,.005,medium_gap_pdf,GTR+second,5,2.0,1.0\\R19\\rose.aln.true.fasta'
##    alignment_file=alignment_file.replace('\\','/')
##    tree_file='/home/mikenute/Phylolab/share/Data\\RoseDNA\\1000M4.tar\\1000,1000,.00003,.005,medium_gap_pdf,GTR+second,5,2.0,1.0\\R19\\rose.tt'
##    tree_file=tree_file.replace('\\','/')
##    job_name='test1'
##    myjob=pasta_job(alignment_file,tree_file,job_name)
##    myjob.make_pbs_file()
##    myjob.make_raw_fasta()
##    myjob.make_sftp_file()

def main3(cfg):
    # main2(cfg)
   mypasta=groupOfCompletedPastaRuns(cfg)
   # mypasta.make_results_sftp_files()
   mypasta.make_tree_batch_file('get_tree_errors.bat')
   mypasta.make_alignment_batch_file('get_aln_errors.sh')
   import subprocess
   subprocess.call(['bash','get_aln_errors.sh'])
   # mypasta.assemble_all_results('R_data.txt')
   # mypasta.assemble_times_results()
##    collect_garbage()

class groupOfCompletedPastaRuns:
    def __init__(self, config=None):
        self.pf=get_prefix()
        self.prefixes={} #old I think
        self.job_list={}
        self.output_list={}

        if config<>None:
            self.cfg=config

        self.read_job_list()
        self.read_output_list()




        # specific locations of folders
##        if start_prefixes==None:
##            self.prefixes['tree_error_prefix']='results/indelible/error_calcs_tree/'
##            self.prefixes['alignment_error_prefix']='results/indelible/error_calcs_aln/'
##            self.prefixes['general_repo_prefix']='results/indelible/'
##        else:
##            for i in start_prefixes.keys():
##                self.prefixes[i]=start_prefixes[i]


    def get_prefix_commands(self,which_os=None):
        """
        returns a list containing [command to set environment variable,
        environment variable reference, folder string]
        """
##        print which_os
        if which_os==None:
            which_os=platform.system()
        if which_os=='Linux':
            line1='PREFIX=/home/mikenute/Phylolab/share/' + '\n'
            line2='/home/mikenute/Phylolab/share'
            return [line1,'$PREFIX',line2]
        else:
            line1='SET PREFIX=C:/users/michael/grad school stuff/research/phylogenetics/\n'
            line2='C:/users/michael/grad school stuff/research/phylogenetics'
            return [line1,'%PREFIX%',line2]

    def read_job_list(self):
        """
        reads list of jobs created by the pastaRunsPreparer earlier.
        """
        joblistfile=self.pf + self.cfg.get('results','job_list')
        jf=open(joblistfile,'r')
        headers=jf.readline().strip().split('\t')
        for line in jf:
            argsdict={}
            a=line.strip().split('\t')
            for i in range(len(headers)):
                argsdict[headers[i]]=a[i]
            self.job_list[a[3]]=argsdict

    def read_output_list(self):
        """reads list of output files constructed from inventory_files.py.
        """
        outputfile=self.pf + self.cfg.get('results','output_files_list')
        of=open(outputfile,'r')
        for line in of:
            args={}
            a=line.strip().split('\t')
            job=a[1].replace('.out.txt','')
            args['job']=job
            args['folder']=a[0]
            args['wall_time']=a[2]
            args['aln']=a[3]
            args['tree0']=a[4]
            args['tree1']=a[5]
            args['tree2']=a[6]
            args['tree_final']=a[7]
            self.output_list[job]=args

    def make_results_sftp_files(self):
        aln_sftp=self.pf + self.cfg.get('pasta_batch_config','results_folder')+'sftp_commands_alignment.txt'
        af=open(aln_sftp,'w')
        af.write('cd ' + self.cfg.get('pasta_batch_config','cincy_dump_folder') + '\n')
        af.write('lcd ' + self.pf +self.cfg.get('results','aln_loc') + '\n')
        for i in self.output_list.keys():
            if not os.path.exists(self.pf + self.cfg.get('results','aln_loc') + self.output_list[i]['aln']):
                if self.output_list[i]['aln']<>'':
                    af.write('get ' + self.output_list[i]['folder'] + '/' + self.output_list[i]['aln'] + '\n')
        af.close()

        tree_sftp=self.pf + self.cfg.get('pasta_batch_config','results_folder')+'sftp_commands_tree.txt'
        tf=open(tree_sftp,'w')
        tf.write('cd ' + self.cfg.get('pasta_batch_config','cincy_dump_folder') + '\n')
        tf.write('lcd ' + self.pf +self.cfg.get('results','tree_loc') + '\n')
        for i in self.output_list.keys():
            if not os.path.exists(self.pf + self.cfg.get('results','tree_loc') + self.output_list[i]['tree_final']):
                if self.output_list[i]['tree_final']<>'':
                    print self.output_list[i]['tree_final']
                    tf.write('get ' + self.output_list[i]['folder'] + '/' + self.output_list[i]['tree_final']+'\n')
        tf.close()

        bf=open(self.pf + self.cfg.get('pasta_batch_config','results_folder')+'bash_sftp_trees_from_cincy.sh','w')
        bf.write('sftp -b sftp_commands_alignment.txt -oPort=2222 mikenute@192.168.200.20 \n')
        bf.write('sftp -b sftp_commands_tree.txt -oPort=2222 mikenute@192.168.200.20 \n')


##    def read_runs_list(self,runs_list):
##        """
##        reads the list of runs according to the specifications.
##        """
##        rf=open(runs_list,'r')
##        runsflag=False
##        row=1
##        headers=None
##
##        for line in rf:
##            a=line.strip().split('\t')
##            if runsflag==False:
##                if a[0]<>'':
##                    self.prefixes[a[0]]=a[1]
##                else:
##                    runsflag=True
##                    headerRow=row+1
##            else:
##                if row==headerRow:
##                    headers=a
##                elif a[0]<>'':
##                    argsdict={}
##                    for i in range(0,len(headers)):
##                        argsdict[headers[i]]=a[i]
##                        self.run_list[a[0]]=argsdict
##            row=row+1

    def get_pasta_job_name(self,pasta_run=None,which_os=None):
        """
        finds the job name from the folder of a given pasta run, and returns None if there is no ".out.txt" file in the folder or if there are multiple ones, and prints an error message
        """
        prefs=self.get_prefix_commands(which_os)
        folder = prefs[2] + self.prefixes['results_prefix'] + pasta_run['results_folder']
##        print [i[:8] for i in os.listdir(folder)]
        jobs=[i[:-8] for i in os.listdir(folder) if i[-8:]=='.out.txt']

        if len(jobs)>1:
            print 'Job ' + pasta_run['name']+' in folder ' + self.prefixes['results_prefix'] + pasta_run['results_folder'] + ' contains results from multiple pasta runs:'
            print jobs
            return None
        elif jobs==None or len(jobs)==0:
            print 'Job ' + pasta_run['name']+' in folder ' + self.prefixes['results_prefix'] + pasta_run['results_folder'] + ' does not contain any valid console output files.'
            return None
        else:
            return jobs[0]

    def make_tree_batch_file_row(self,pasta_run=None, pasta_output=None,which_os='Windows',iteration=None):
        """
        makes a single row of a batch file that will run the tree-errors
        on all of the files in our pasta runs.
        """
        if iteration==None:
            iter_str='.tre'
            iternum='_final'
        else:
            iter_str='_temp_iteration_'+ str(iteration) + '_tree.tre'
            iternum=str(iteration)


        # get the commands
        prefs=self.get_prefix_commands(which_os)
        prefix=prefs[2]

        command=''
        if which_os=='Linux':
            command = 'perl '

        # set the files to use
        myargs={}
##        myargs['truetreeloc']=self.pf + self.prefixes['data_prefix'] + pasta_run['tree_folder'] + pasta_run['tree_name']
        myargs['truetreeloc']=prefix + pasta_run['tree_file'].replace(get_prefix('Linux'),'')
##        myargs['treeloc']=prefix + self.prefixes['results_prefix'] + pasta_run['results_folder'] + self.get_pasta_job_name(pasta_run,which_os)+iter_str
        myargs['treeloc']=prefix + self.cfg.get('results','tree_loc') + pasta_output['tree'+iternum]
##        myargs['outputloc']=prefs[2] + self.prefixes['tree_error_prefix']+'tree_error_'+pasta_run['name'] + iternum + '.txt'
        myargs['outputloc']=prefix + self.cfg.get('results','tree_error_loc') + 'tree_error_' + pasta_output['job'] + '.txt'


        # make the line
        if which_os=='Windows':
            outchar='2'
        else:
            outchar='&'

        line=command + 'CompareTree.pl -tree \"%(treeloc)s\" -versus \"%(truetreeloc)s\" ' % myargs +outchar+'>  \"%(outputloc)s\"\n' % myargs
        return line

    def make_tree_batch_file(self,batchfile_name=None, which_os='Windows', iteration=None):
        """
        assembles the batch file to compute tree errors on all the files from this run.
        """
        prefs=self.get_prefix_commands(which_os)
        fn=self.pf+self.cfg.get('pasta_batch_config','results_folder')+batchfile_name
        myf=open(fn,'w')
        if which_os=='Windows':
            myf.write('set CURDIR=%CD%\n')
##        myf.write(prefs[0])

        myf.write('cd \"' + prefs[2] + '/code/tree_comparison/\"\n' )
        for i in self.output_list.keys():
            try:
                job=self.job_list[i]
            except:
                print 'Pasta Job: ' + i + ' is not in the given job list'
                continue
            output=self.output_list[i]
            myf.write(self.make_tree_batch_file_row(job,output,which_os,iteration))
        if which_os=='Windows':
            myf.write('cd %CURDIR%\n')
        myf.close()

    def get_final_alignment_file_name(self,pasta_run,which_os):
        jobn=self.get_pasta_job_name(pasta_run,which_os)
        prefs=self.get_prefix_commands(which_os)
        folder = prefs[2] + self.prefixes['results_prefix'] + pasta_run['results_folder']
        expr=jobn+'.marker001'
        k=len(expr)
        alns=[i for i in os.listdir(folder) if i[:k]==expr]
        if len(alns)>1:
            print jobn + ' contains multiple alignment files: '
            print alns
        elif alns==None or len(alns)==0:
            print jobn + ' does not contain any files matching the alignment pattern'
        else:
            return alns[0]


    def make_alignment_batch_file_row(self,pasta_run=None, pasta_output=None):
        """
        makes a single row of a batch file that will run the tree-errors
        on all of the files in our pasta runs.
        """
##        if iteration==None:
##            #get final alignment file name
##            iternum='_final'
##        else:
##            iter_str='_temp_iteration_'+ str(iteration) + '_seq_alignment.txt'
##            iternum='_iter'+iteration

##        # get the commands
##        prefs=self.get_prefix_commands(which_os)
##        prefix=prefs[2]
##        command=''
##        if platform.system()=='Linux':
##            command = 'perl '

        # set the files to use
        myargs={}
##        myargs['truealnloc']=self.pf + self.prefixes['data_prefix'] + pasta_run['alignment_folder'] + pasta_run['alignment_name']
        myargs['truealnloc']=pasta_run['alignment_file']
##        myargs['alnloc']=prefix + self.prefixes['results_prefix'] + pasta_run['results_folder'] + self.get_final_alignment_file_name(pasta_run,which_os)
        myargs['alnloc']=self.pf + self.cfg.get('results','aln_loc') + pasta_output['aln']
##        myargs['outputloc']=prefs[2] + self.prefixes['alignment_error_prefix']+'aln_error_'+pasta_run['job'] + '.txt'
        myargs['outputloc']=self.pf + self.cfg.get('results','aln_error_loc')+'aln_error_'+pasta_output['job'] + '.txt'

        if platform.system()=='Windows':
            outchar='2'
        else:
            outchar='&'

        line='java -jar '+self.pf+'/code/misc/FastSP_1.6.0.jar -r \"%(truealnloc)s\" -e \"%(alnloc)s\" ' % myargs + outchar + '> \"%(outputloc)s\"\n' % myargs
        line = line + 'echo "'+pasta_output['job']+ ' done"\n'
        return line

    def make_alignment_batch_file(self,batchfile_name=None):
        """
        assembles the batch file to compute alignment errors on all the files from this run.
        """
##        prefs=self.get_prefix_commands(which_os)
        fn=self.pf+self.cfg.get('pasta_batch_config','results_folder')+batchfile_name
        myf=open(fn,'w')

        for i in self.output_list.keys():
            job=self.job_list[i]
            out=self.output_list[i]
            myf.write(self.make_alignment_batch_file_row(job,out))
        myf.close()

    def make_running_time_results_file(self,file_name,which_os='Windows'):
        prefs=self.get_prefix_commands(which_os)
        prefix=prefs[2]
        outfile = prefix + self.prefixes['general_repo_prefix'] + file_name
        out=open(outfile,'w')
        out.write('folder\tjob_name\trunning_time\t\n')

        rtregex=re.compile('PASTA INFO: Total time spent: (?P<time>\d+\.*\d*)s')
        for i in self.run_list.keys():
            a=self.run_list[i]
            jobn=self.get_pasta_job_name(a,which_os)
            rf=prefix + self.prefixes['results_prefix'] + a['results_folder'] + jobn + '.out.txt'
            consoleout=open(rf,'r')
            consoletext=consoleout.read()
            consoleout.close()
            data=rtregex.search(consoletext)
            myargs={}
            myargs['folder']=a['results_folder']
            myargs['job']=jobn
            try:
                myargs['time']=data.group('time')
            except:
                myargs['time']=''
            out.write('%(folder)s\t%(job)s\t%(time)s\t\n' % myargs)
        out.close()

    def assemble_all_results(self,file_name):
        rf=open(self.pf + self.cfg.get('pasta_batch_config','results_folder')+file_name,'w')
        rf.write('job\tmodel\tset\tsize\twall_time\tsplitsfound\tsplits_total\tsplits_frac\tsplts_maxlndf\tsplits_ratio\tfile_name\tSP-Score\tModeler\tSPFN\tSPFP\tCompression Ratio\tTC\tMaxLenNoGap\tNumSeq\tLenRef\tLenEst\tCells\t\n')
        str_line0='%(job)s\t%(model)s\t%(set)s\t%(size)s\t%(wall_time)s\t'
        str_line1='%(splitsfound)s\t%(total)s\t%(frac)s\t%(maxlndf)s\t%(ratio)s\t'
        str_line2='%(file_name)s\t%(sp)s\t%(modeler)s\t%(spfn)s\t%(spfp)s\t%(comp)s\t%(tc)s\t%(maxlen)s\t%(numseq)s\t%(lenref)s\t%(lenest)s\t%(cells)s\t\n'

        alnfolder=self.pf + self.cfg.get('results','aln_error_loc')
        treefolder=self.pf + self.cfg.get('results','tree_error_loc')

        for i in self.output_list.values():
            job=i['job']
            a=job.split("_")
            args={'job':job, 'model':a[1], 'set': a[2], 'size':a[3], 'wall_time': i['wall_time']}
            str_out0=str_line0 % args
            treefn='tree_error_' + job + '.txt'
            alnfn='aln_error_' + job + '.txt'
            str_out1=str_line1 % read_CompareTree_results_file(treefolder + treefn)
            str_out2=str_line2 % read_alignment_results_file(alnfolder + alnfn)
            rf.write(str_out0 + str_out1 + str_out2)

    def assemble_times_results(self,file_name=None):
        folder= self.pf + self.cfg.get('results','times_folder')
        files=os.listdir(folder)

        outf=open(self.pf + self.cfg.get('pasta_batch_config','results_folder')+'times_log_summary.txt','w')
        outf.write('file\tinitalign_mafft\tinitalign-hmmeralign\tinittree_fasttree\taln0\titer0-centroid\taln1\titer1-centroid\taln2\titer2-centroid\n')
        k=1
        n=len(files)
        for l in files:
            out=read_pasta_times_log_file(folder + l)
            str_line='%(fi)s\t%(im)s\t%(ih)s\t%(if)s\t%(a0)s\t%(t0)s\t%(a1)s\t%(t1)s\t%(a2)s\t%(t2)s\n'
            outf.write(str_line % out)
            print str(k) + ' of ' + str(n) + ': ' + l
            k+=1
        outf.close()
        pass

def read_pasta_times_log_file(filepath=None):
    cats={'initalign_mafft':[],'initalign-hmmeralign':[],'inittree_fasttree':[],
    'iter0-centroid':[],'iter1-centroid':[],'iter2-centroid':[]}
    reset=re.compile('initalign_mafft queued')
    myf=open(filepath,'r')
    linereg=re.compile('\[(?P<time>.*)\.\d*\] : (?P<cat>.*?) ')
    dt_format='%m/%d/%y %H:%M:%S'
    for l in myf:
        myreset=None
        myreset=reset.search(l)
        if myreset<>None:
            cats={'initalign_mafft':[],'initalign-hmmeralign':[],'inittree_fasttree':[],'iter0-centroid':[],'iter1-centroid':[],'iter2-centroid':[]}

        a=linereg.search(l)
##        print a.group('time')
        mydate=datetime.datetime.strptime(a.group('time'),dt_format)
        cat = ''
        cat = a.group('cat')
        if cat in cats.keys():
            cats[cat].append(mydate)
    p,f=os.path.split(filepath)
    args={'fi':f}
    try:
        args['im']=str(max(cats['initalign_mafft'])-min(cats['initalign_mafft']))
        args['ih']=str(max(cats['initalign-hmmeralign'])-min(cats['initalign-hmmeralign']))
        args['if']=str(max(cats['inittree_fasttree'])-min(cats['inittree_fasttree']))
    except:
        args['im']=''
        args['ih']=''
        args['if']=''
    try:
        args['a0']=str(min(cats['iter0-centroid'])-max(cats['inittree_fasttree']))
        args['t0']=str(max(cats['iter0-centroid'])-min(cats['iter0-centroid']))
    except:
        args['a0']=''
        args['t0']=''
    try:
        args['a1']=str(min(cats['iter1-centroid'])-max(cats['iter0-centroid']))
        args['t1']=str(max(cats['iter1-centroid'])-min(cats['iter1-centroid']))
    except:
        args['a1']=''
        args['t1']=''
    try:
        args['a2']=str(min(cats['iter2-centroid'])-max(cats['iter1-centroid']))
        args['t2']=str(max(cats['iter2-centroid'])-min(cats['iter2-centroid']))
    except:
        args['a2']=''
        args['t2']=''
    return args
##    for i in cats.keys():
##        print i + ' - ' + str(min(cats[i]))
##    print str(mydate)
##    print a.group('cat')



def read_CompareTree_results_file(filepath=None):
    myf=open(filepath,'r')
    dataregs=re.compile('Splits	Found	(?P<splitsfound>\d+)	Total	(?P<total>\d+)	Frac	(?P<frac>\d+\.*\d*)	MaxLnDf	(?P<maxlndf>\d*\.*\d*)	Ratio	(?P<ratio>\d*\.*\d*)	MaxBtDf')
    text=myf.read()
    myf.close()
    data=dataregs.search(text)
    myargs={}
    outputs=['splitsfound','total','frac','maxlndf','ratio']
    for i in outputs:
        try:
            myargs[i]=data.group(i)
        except:
            myargs[i]=''
    return myargs

def tree_results_folder_to_text(folder=None,output_file_path=None):
    outfile=open(output_file_path,'w')
    firstline='file_name\tsplitsfound\ttotal\tfrac\tmaxlndf\tratio\t\n'
    outfile.write(firstline)
    folder=folder.replace('\\','/') + '/'
    folder=folder.replace('//','/')

    for i in os.listdir(folder):
        myargs={}
        myargs=read_CompareTree_results_file(folder + i)
        myargs['name']=i
        line='%(name)s\t%(splitsfound)s\t%(total)s\t%(frac)s\t%(maxlndf)s\t%(ratio)s\t\n' % myargs
        outfile.write(line)
    outfile.close()

def read_alignment_results_file(filepath=None):
    pa, fi=os.path.split(filepath)
    myargs={}
    myargs['file_name']=fi

    # open the results file
    textfile=open(filepath,'r')
    text=textfile.read()
    textfile.close()

    # regexes to pull the results
    file_regex='SP-Score (?P<sp>\d+\.\d+[E\-\d]*)[.\n]*Modeler (?P<modeler>\d+\.\d+[E\-\d]*)[.\n]*SPFN (?P<spfn>\d+\.\d+[E\-\d]*)[.\n]*SPFP (?P<spfp>\d+\.\d+[E\-\d]*)[.\n]*Compression (?P<comp>\d+\.\d+[E\-\d]*)[.\n]*TC (?P<tc>\d+\.\d+[E\-\d]*)'
    file_regex_2='MaxLenNoGap= (?P<maxlen>\d+).*NumSeq= (?P<numseq>\d+).*LenRef= (?P<lenref>\d+).*LenEst= (?P<lenest>\d+).*Cells= (?P<cells>\d+)'
    fileregs=re.compile(file_regex)
    fileregs2=re.compile(file_regex_2)

    vals1list=['sp','modeler','spfn','spfp','comp','tc']
    vals2list=['maxlen','numseq','lenref','lenest','cells']
    vals1=fileregs.search(text)
    vals2=fileregs2.search(text)
    for i in vals1list:
        try:
            myargs[i]=vals1.group(i)
        except:
            myargs[i]=''
    for i in vals2list:
        try:
            myargs[i]=vals2.group(i)
        except:
            myargs[i]=''
    return myargs

def alignment_results_folder_to_text(folder=None, output_file_path=None):
    # quick checks on folder string
    folder=folder.replace('\\','/') + '/'
    folder=folder.replace('//','/')

    # model output line:
    line_str='%(file_name)s\t%(sp)s\t%(modeler)s\t%(spfn)s\t%(spfp)s\t%(comp)s\t%(tc)s\t%(maxlen)s\t%(numseq)s\t%(lenref)s\t%(lenest)s\t%(cells)s\t\n'

    outfile=open(output_file_path,'w')
    outfile.write('file_name\tSP-Score\tModeler\tSPFN\tSPFP\tCompression Ratio\tTC\tMaxLenNoGap\tNumSeq\tLenRef\tLenEst\tCells\t\n')

    for i in os.listdir(folder):
        myargs={}
        myargs=read_alignment_results_file(folder + i)
        outfile.write(line_str % myargs)

    outfile.close()

    pass

def get_prefix(os=None):
    if os==None:
        os=platform.system()
    prefix_general_windows='C:/Users/Michael/Grad School Stuff/Research/Phylogenetics'
    prefix_general_linux='/home/mikenute/Phylolab/share'
    if os=='Windows':
        prefix_general=prefix_general_windows
    else:
        prefix_general=prefix_general_linux
    return prefix_general

def main():
    #set prefixes
    my_prefixes={}
    # my_prefixes['tree_error_prefix']='results/indelible/error_calcs_tree/'
    # my_prefixes['alignment_error_prefix']='results/indelible/error_calcs_aln/'
    # my_prefixes['general_repo_prefix']='results/indelible/'
    # prefix='C:/Users/Michael/Grad School Stuff/Research/Phylogenetics/'
    prefix='/home/mikenute/Phylolab/share/'
    mine= groupOfCompletedPastaRuns(prefix + 'results/indelible/indelible_pasta_job_list.txt',my_prefixes)
##
##    mine.make_tree_batch_file('make_final_tree_errors.bat','Windows')
##    mine.make_tree_batch_file('make_iter0_tree_errors.bat','Windows',0)
##    mine.make_tree_batch_file('make_iter1_tree_errors.bat','Windows',1)
##    mine.make_tree_batch_file('make_iter2_tree_errors.bat','Windows',2)
    # mine.make_alignment_batch_file('make_final_alignment_erors.sh','Linux')

    folder=prefix + 'results/indelible/error_calcs_tree/'
    outputfile=prefix + 'results/indelible/tree_errors_all.txt'
    tree_results_folder_to_text(folder,outputfile)
    print 'tree results done'

    aln_folder= prefix + my_prefixes['alignment_error_prefix']
    outputfile = prefix + my_prefixes['general_repo_prefix'] + 'alignment_errors_all.txt'
    alignment_results_folder_to_text(aln_folder,outputfile)
    print 'alignment results done'

    mine.make_running_time_results_file('running_times.txt','Windows')
    print 'running times done'


if __name__ == '__main__':
##    main2()

##    print platform.system()=='Linux'
    # config_file=get_prefix()+'/results/rnasim_v2/pasta_batch_config.cfg'
    config_file=get_prefix()+'/results/rosedna_L/pasta_batch_config.cfg'
    cfg=ConfigParser.ConfigParser()
    cfg.read(config_file)
    main3(cfg)

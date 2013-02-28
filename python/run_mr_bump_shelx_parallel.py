#!/usr/bin/env python

# run mrbump anr rebuild with shelx

# python imports
import os
import re
import sys
import glob
import multiprocessing
import unittest

# our imports
import mrbump_cmd
import printTable

###################



def make_mrbump_desktop( jobid, ensemble_pdb, amopt ):
    
    path = os.getcwd()
    mr_bump = open(path+'/'+jobid+'.sub', "w")
    mr_bump.write('#!/bin/sh\n')
    jobstr = mrbump_cmd.mrbump_cmd( amopt.d, jobid=jobid, ensemble_pdb=ensemble_pdb )
    mr_bump.write(jobstr)
    mr_bump.close()
    os.system('chmod uoga=wrx '+ path + '/'+jobid+'.sub')
    os.system(path + '/'+jobid+'.sub  >'+jobid+'.log')
    os.chdir(path)
    
##End make_mrbump_desktop

##########################
def make_mrbump_desktop_domain(sigf, fp, free, jobid, local_files, mtz, seq, fixed_pdb, FIXED_INPUT, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles,  NoShelx, NoShelxCycles, MRkeys):
    path = os.getcwd()

    mr_bump = open(path+'/'+jobid+'.sub', "w")
    mr_bump.write('#!/bin/sh\n')


    mr_bump.write('mrbump HKLIN ' +mtz + ' SEQIN ' + seq+' HKLOUT ' + 'OUT.mtz  XYZOUT OUT.pdb << eof\n'  )

    mr_bump.write('LABIN ' + sigf + ' ' + fp + ' ' + free + '\n')
    mr_bump.write('JOBID '+ jobid+ '_mrbump\n')
    mr_bump.write('MRPROGRAM ' + mrbump_programs + '\n')
    mr_bump.write('LOCALFILE ' + local_files + ' CHAIN ALL RMS 1.2\n')

    if FIXED_INPUT:
        mr_bump.write('FIXED_XYZIN '+fixed_pdb+' IDENTIY 0.6 \n')
    if not FIXED_INPUT:
        mr_bump.write('LOCALFILE ' +fixed_pdb + ' CHAIN ALL RMS 1.2\n')

    mr_bump.write('SCOPSEARCH False\n')
    mr_bump.write('PQSSEARCH False\n')
    mr_bump.write('SSMSEARCH False\n')

    mr_bump.write('FAST False\n')
    mr_bump.write('DOFASTA False\n')

    mr_bump.write('MDLD False\n')
    mr_bump.write('MDLC False\n')
    mr_bump.write('MDLM False\n')
    mr_bump.write('MDLP False\n')
    mr_bump.write('MDLS False\n')
    mr_bump.write('MDLU True\n')
    mr_bump.write('UPDATE False\n')

    mr_bump.write('BUCC  '+str(Buccaneer)+'\n')
    mr_bump.write('BCYCLES  '+str(Buccaneer_cycles)+'\n')
    mr_bump.write('ARPWARP  '+str(arpwarp)+'\n')
    mr_bump.write('ACYCLES  '+str(arpwarp_cycles)+'\n')
    mr_bump.write('SHELXE  '+str(NoShelx)+'\n')
    mr_bump.write('SCYCLES  '+str(NoShelxCycles)+'\n')

    mr_bump.write('FIXSG True\n')
    mr_bump.write('PJOBS 1\n')

    mr_bump.write('CHECK False\n')
    mr_bump.write('LITE True\n')
    mr_bump.write('PICKLE False\n')
    mr_bump.write('TRYALL True\n')
    mr_bump.write('USEACORN False\n')
    mr_bump.write('USEENSEM False\n')
    mr_bump.write('DEBUG True\n')

    for k in MRkeys:
        mr_bump.write(k + '\n')

    mr_bump.write('END\n')
    mr_bump.write('eof')
    mr_bump.close()

    os.system('chmod uoga=wrx '+ path + '/'+jobid+'.sub' )
##RUN
    os.system(path + '/'+jobid+'.sub &>'+jobid+'_log')

    os.chdir(path)


###################
def make_mrbump_Cluster(sigf, fp, free, jobid, local_files, mtz, seq, mrbump_programs):
    path = os.getcwd()

    mr_bump = open(path+'/jobid.sub', "w")
    mr_bump.write('#!/bin/sh\n'
     '#$ -j y\n' +
     '#$ -cwd\n' +
     '#$ -w e\n' +
     '#$ -V\n' +
     '#$ -o '+jobid+'.log\n' +
     '#$ -N Z'+jobid+'\n\n')

    mr_bump.write('mrbump HKLIN ' + jobid.lower() + '.mtz' + ' SEQIN ' + jobid + '_.fasta HKLOUT ' + 'MRBUMP_OUT/OUT.mtz  XYZOUT /MRBUMP_OUT/OUT.pdb << eof\n'  )

    mr_bump.write('LABIN ' + sigf + ' ' + fp + ' ' + free + '\n')
    mr_bump.write('JOBID '+ jobid+ '_mrbump\n')
    mr_bump.write('MRPROGRAM ' + mrbump_programs + '\n')
    mr_bump.write('LOCALFILE ' + local_files + ' CHAIN A RMS 1.2\n')


    mr_bump.write('SCOPSEARCH False\n')
    mr_bump.write('PQSSEARCH False\n')
    mr_bump.write('SSMSEARCH False\n')

    mr_bump.write('FAST False\n')
    mr_bump.write('DOFASTA False\n')

    mr_bump.write('MDLD False\n')
    mr_bump.write('MDLC False\n')
    mr_bump.write('MDLM False\n')
    mr_bump.write('MDLP False\n')
    mr_bump.write('MDLS False\n')
    mr_bump.write('MDLU True\n')
    mr_bump.write('UPDATE False\n')


    mr_bump.write('FIXSG True\n')
    mr_bump.write('PJOBS 1\n')

    mr_bump.write('CHECK False\n')
    mr_bump.write('LITE True\n')
    mr_bump.write('PICKLE False\n')
    mr_bump.write('TRYALL True\n')
    mr_bump.write('USEACORN False\n')
    mr_bump.write('USEENSEM False\n')
    mr_bump.write('DEBUG True\n')
    mr_bump.write('END\n')
    mr_bump.write('eof')
    mr_bump.close()

    os.system('chmod uoga=wrx '+ path + '/'+ jobid +'/mr_bump_' + jobid +' >mrbump_log')
##RUN
    #os.system(path + '/'+ jobid +'/mr_bump_' + jobid)

    os.chdir(path)

##################################

def remove_hires(mtz, newpath, name):
    path = os.getcwd()


    occupancy = open(path + '/hires', "w")
    occupancy.write('#!/bin/csh -f\n')
    occupancy.write('mtzutils hklin  ' + mtz + ' hklout ' + newpath + '/' + name + '.mtz <<eof\n')
    occupancy.write('resolution 40.823 1.0\n')
    occupancy.write('END\neof')
    occupancy.close()
    os.system('chmod uoga=wrx '+ path + '/hires')
##RUN
    os.system(path + '/hires')



################################
def get_space_group(mtz):
    spacegroup = ''

    path = os.getcwd()
    os.system('mtzdmp ' + mtz + ' >mtzdmp_out')
    mtz_out = open(path + '/mtzdmp_out')
    for line in mtz_out:

        get_stat1 = re.compile('\* Space group =\s*\'(.*)\'')
        space_test = get_stat1.search(line)
        if space_test:
            space_split = re.split (get_stat1, line)
            #print space_split[1]
            spacegroup = space_split[1]

    return spacegroup

################################
def get_ensemble(fine_path, working_dir):
    truncs = ['1','2','3','4','5','6']
    rads = ['ALIGNED_rad_1','ALIGNED_rad_2','ALIGNED_rad_3','ALL']

    ensembles = []

    for trunc in truncs:

        for rad in rads:
            trunc_path = fine_path+'/THESEUS_polyala_cluster0_trunc_'+ trunc+'/'+rad

            if os.path.exists(trunc_path):

                if rad !='ALL':
                    ensemble = trunc_path+'/cluster_0_sup.pdb'
                    if os.path.exists(ensemble):
                        #print  ensemble
                        os.system('cp ' +ensemble + ' '+ working_dir+ '/trunc'+trunc+'rad_'+rad+'.pdb')
                        ensembles.append('trunc'+trunc+'rad_'+rad+'.pdb')

                if rad =='ALL':
                    ensemble = trunc_path+'/ALL_sup.pdb'
                    if os.path.exists(ensemble):
                        os.system('cp ' +ensemble + ' '+ working_dir+ '/trunc'+trunc+'rad_'+rad+'.pdb')
                        ensembles.append('trunc'+trunc+'rad_'+rad+'.pdb')

    return ensembles
###############################
def pdbcur(pdb):

    ASU = 'nan'
    curr_dir = os.getcwd()
    cur= open(curr_dir + '/pdbcur', "w")
    cur.write('#!/bin/sh\n'+
    'pdbcur xyzin  '+pdb + '<<EOF\n'+
    'SUMM \n'+
    'EOF')
    cur.close()
    os.system('chmod uoga=wrx '+curr_dir + '/pdbcur' )
    os.system(curr_dir + '/pdbcur >cur_out' )

    pattern =  re.compile('\s*Number of models =\s*(\d*)')
    cur_out =open( curr_dir + '/cur_out')
    for line in cur_out:
        if re.search(pattern, line):
                #print line
            split = re.split(pattern, line)

            ASU = split[1]
    return ASU

###############################
def  make_MRBUMP_run( pdb, name, amopt ):

    #make script
    run_dir = amopt.d['mrbump_dir']
    os.chdir(run_dir)

    make_mrbump_desktop(name, pdb, amopt)
    table = run_dir + '/search_'+name+'_mrbump/results/resultsTable.dat'
    DIE = get_table(table, run_dir, amopt.d['use_buccaneer'],  amopt.d['use_arpwarp'],  amopt.d['use_shelxe'])
    
    # Nothing below here seems to be used
    # get data for printtable
    SolutionFound = False
    #files to make:
    phaser_mtz = 'fail'
    phaser_pdb = 'fail'
    molrep_mtz = 'fail'
    molrep_pdb = 'fail'
    molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT = 'none'
    phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = 'none'
    #Phaser output:
    if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
        phaser_shelxscore ='none'

        phaser_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
        phaser_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'

    #molrep output:
    if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
        molrep_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
        molrep_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'

    os.chdir(run_dir)

    return DIE
##End make_MRBUMP_run


def make_MRBUMP_run_domain(mtz, pdb, run_dir, fasta, name, fixed_pdb, sigf, FP, free, FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles,  NoShelx, NoShelxCycles, MRkeys):

    #files to make:
    phaser_mtz = 'fail'
    phaser_pdb = 'fail'
    molrep_mtz = 'fail'
    molrep_pdb = 'fail'
    molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT = 'none'
    phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = 'none'

    #get names of flags

    os.chdir(run_dir)

    sigf = 'SIGF='+sigf
    FP = 'F='+FP
    free =  'FreeR_flag='+free

    make_mrbump_desktop_domain(sigf, FP, free, name, pdb, mtz, fasta, fixed_pdb, FIXED_INPUT, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles,  NoShelx, NoShelxCycles, MRkeys)

    if os.path.exists(run_dir + '/search_'+name+'_mrbump/results/resultsTable.dat'):

        DIE = get_table(run_dir + '/search_'+name+'_mrbump/results/resultsTable.dat', run_dir, Buccaneer, arpwarp, NoShelx)

    if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
        phaser_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
        phaser_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'

    #molrep output:
    if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
        molrep_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
        molrep_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'


    os.chdir(run_dir)
    return DIE
##End make_MRBUMP_run_domain

def get_table(table, run_dir, Buccaneer, arpwarp, use_shelx):

    DIE = False

    sys.stdout.write('\n')
    sys.stdout.write('###########################################################################################\n')
    sys.stdout.write('###########################################################################################\n')
    sys.stdout.write('##                                                                                       ##\n')
    sys.stdout.write('##                                                                                       ##\n')
    print 'Results for this Run: '


    header = ''
    for line in open(table):
        if not re.search('PHASER', line)  and not re.search('MOLREP', line):
            sys.stdout.write(line)
            header = line
        else:
            sys.stdout.write(line)

    # Set up the results table headers
    resultsTable = []
    HEADER=True

    table = []
    sys.stdout.write('\n\n  Overall Summary:\n\n')
    #print header
    for search in os.listdir(run_dir):
        if os.path.isdir(run_dir+'/'+search):
            res = run_dir+'/'+search+'/results/resultsTable.dat'
            if os.path.exists(res):
                for line in open(res):
                    if  'MR_PROGRAM' in line.upper() and 'SOLUTION' in line.upper() and HEADER:
                        HEADER=False
                    if  re.search('PHASER', line)  or re.search('MOLREP', line):
                        resultsTable.append(line.split())
                        #print line

                        split =re.split('\s*', line)

                        table.append(split)

    header  =header.split()

    Best = 'none'

    # sort results
    if use_shelx:
        y = header.index('SHELXE_CC')
    # print resultsTable
    # print y
        resultsTable.sort(key=lambda x: float(x[y]))
        resultsTable.reverse()
        best = resultsTable[0][0]
        prog = resultsTable[0][1]
        if float( resultsTable[0][y]) >=25:
            DIE = True

    # print resultsTable
    if not use_shelx:
        y = header.index('final_Rfree')

        resultsTable.sort(key=lambda x: float(x[y]))
        best = resultsTable[0][0]
        prog = resultsTable[0][1]

    n= re.sub('loc0_ALL','search', best)
    n = re.sub('UNMOD','mrbump', n)
    Best = run_dir +n+'/data/'+re.sub('_UNMOD','', best)+'/unmod/mr/'+prog.lower()+'/refine'


    resultsTable.insert(0, header)
    # for x in resultsTable:
    #   print x

    # Output the results table
    T=printTable.Table()

    out = sys.stdout
    T.pprint_table(out, resultsTable)
    sys.stdout.write("\n\n\n")
    print 'Best results so far are in : '
    print Best



    sys.stdout.write('##                                                                                       ##\n')
    sys.stdout.write('##                                                                                       ##\n')
    sys.stdout.write('###########################################################################################\n')
    sys.stdout.write('###########################################################################################\n')
    return DIE

#End get_table


############
def split(ensembles, nproc):  #split into parallel comparisons for speed

    divide= len(ensembles) / nproc
    if divide == 0 :
        divide = 1
    chunks = []
    if len(ensembles) == 0:
        print 'no ensembles '
        sys.exit()


    for i in xrange(0, len(ensembles), divide):
        a_chunk =  ensembles[i:i+divide]
        chunks.append(a_chunk)


    if len(chunks)>nproc:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    #print chunks





    return chunks

##############################
def isnumber(n):
    try :
        float(n)
        return True
    except:
        return False




############################
def run_parallel( chunk_of_ensembles, batchnum, amopt ):
    """
    #loops through each chunk
    """
    jobnum = 1
    for each_pdb in chunk_of_ensembles:
        name = re.split('/', each_pdb)
        name = name.pop()
        name = re.sub('.pdb', '', name)
        #print name
        print '=== In batch '+str(batchnum)+' running job '+str(jobnum)+' of '+str(len(chunk_of_ensembles))+'. '+str(len(chunk_of_ensembles)-jobnum)  +' left to go'
        jobnum +=1
        DIE = make_MRBUMP_run(each_pdb, name, amopt)
        if DIE and amopt.d['early_terminate']:
            print 'solution found, exiting'
            sys.exit()

################################
def run_parallel_domain(mtz, chunk_of_ensembles, run_dir, fasta,  log_name,  sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, inc,  fixed_PDB, FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles, MRkeys): #loops through each chunk

    log = open(log_name, "w")

    for each_pdb in chunk_of_ensembles:
        name = re.split('/', each_pdb)
        name = name.pop()
        name = re.sub('.pdb', '', name)
        # print name

        DIE = make_MRBUMP_run_domain(mtz, each_pdb, run_dir, fasta, name, fixed_PDB, sigf, FP, free, FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles,  NoShelx, NoShelxCycles, MRkeys)
        if DIE:
            if EarlyTerminate:
                print 'solution found, exiting'
                sys.exit()


################################
def split_into_runs( ensembles, amopt ):
    
    print 'the ensembles will be divided into '+str(amopt.d['nproc']) +' batches, each job will contain '+str(len(ensembles[0])) +' ensembles'
    
    threads = []
    batchnum=1
    for each_chunk in ensembles:
        #print '===now running job '+str(inc) +' containing  '+str(len(each_chunk)) +'ensembles. '+str(len(ensembles)-inc) +' jobs left to do==='
        thread = multiprocessing.Process(target=run_parallel, args=( each_chunk, batchnum, amopt ) )
        thread.start()
        threads.append(thread)
        batchnum+=1

    for thread in threads:
        thread.join()
        
##End split_into_runs

def split_into_runs_domains(mtz, ensembles, run_dir, domain_all_chain_fasta , NProc, sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles,  domain_all_chains_pdb, FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles, MRkeys):
    cur_dir=os.getcwd()
    inc = 1
    threads = []

    for each_chunk in ensembles:
        log_name = cur_dir+'/LOG_proc'+str(inc)

        thread = multiprocessing.Process(target=run_parallel_domain, args=(mtz, each_chunk, run_dir, domain_all_chain_fasta,  log_name,  sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, inc, domain_all_chains_pdb , FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles, MRkeys) )
        thread.start()
        threads.append(thread)
        inc+=1

    for thread in threads:
        thread.join()

class TestCases(unittest.TestCase):
    """
    Unit tests
    """

    def setUp(self):
        """
        Get paths need to think of a sensible way to do this
        """
        thisdir = os.getcwd()
        self.ampledir = os.path.abspath( thisdir+os.sep+"..")
        self.testdir = self.ampledir + os.sep + "tests"
        self.testfilesdir = self.testdir + os.sep + "testfiles"

    #def testInputString(self):
    #    """Input string for mr bump"""
        
    def testJaclyn(self):
        """
        Jaclyn's tests
        """
        ensembles=[ self.testfilesdir + os.sep + "All_atom_trunc_0.792808_rad_1.pdb",

                   ]
        #                   self.testfilesdir + os.sep + "poly_ala_trunc_17.534083_rad_2.pdb",
        #           self.testfilesdir + os.sep + "SCWRL_Reliable_sidechains_trunc_26.557637_rad_3.pdb"
        
        #for infile in glob.glob( os.path.join('/home/jaclyn/Ample_tests/toxd-example/ensembles_1', '*.pdb') ):
        #    print infile
        #    ensembles.append(infile)
        
        import ample_options
        amopt = ample_options.AmpleOptions()
        
        amopt.d['mtz'] = self.ampledir + '/examples/toxd-example/1dtx.mtz'
        amopt.d['fasta'] = self.ampledir + '/examples/toxd-example/toxd_.fasta'
        
        run_dir = os.getcwd() + os.sep + "MRBUMP"
        if not os.path.isdir( run_dir ):
            os.mkdir( run_dir )
        os.chdir(run_dir)
        
        amopt.d['mrbump_dir']=run_dir
        amopt.d['nproc'] = 1
        amopt.d['SIGF'] = 'SIGFP'
        amopt.d['F'] = 'FP'
        amopt.d['FREE'] = 'FreeR_flag'
        amopt.d['early_terminate'] = True
        amopt.d['use_shelxe'] = True
        amopt.d['shelx_cycles'] = 10
        amopt.d['old_shelx'] = False
        amopt.d['ASU'] = '0'
        amopt.d['mrbump_programs'] = 'molrep'
        amopt.d['use_buccaneer'] = False
        amopt.d['buccaneer_cycles'] = 5
        amopt.d['use_arpwarp'] = False
        amopt.d['arpwarp_cycles'] = 10
        amopt.d['mr_keys'] = []
        
        split_into_runs( [ensembles], amopt )
        
        #run_dir = '/data2/jac45/tox/toxd-example/ROSETTA_MR_1/MRBUMP/'
        #chunk = split(ensembles, nproc)
        #'split_into_runs_domains(mtz, chunk, run_dir, fasta, nproc, fixed_pdb)
        #split_into_runs(mtz,chunk  , run_dir, fasta,  nproc, sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, mrbump_programs)
        #get_table('/data2/jac45/tox/toxd-example/ROSETTA_MR_1/MRBUMP/search_SCWRL_Reliable_sidechains_trunc_31.024398_rad_3_mrbump/results/resultsTable.dat', run_dir,  Buccaneer, arpwarp, use_shelx )
        #make_MRBUMP_run_domain(mtz, pdb, run_dir, fasta, name, fixed_pdb, sigf, FP, free, FIXED_INPUT, SHELX_OLD, mrbump_programs, Buccaneer, Buccaneer_cycles, arpwarp, arpwarp_cycles,  NoShelx, NoShelxCycles)
        ###s

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/python2.6

# run mrbump anr rebuild with shelx

import signal, time
import subprocess
import os, re, sys, glob
import run_shelx


def get_flags (mtz):
  sigf = 'SIGF='
  FP = 'F='
  free = 'FreeR_flag=Unassigned'

  path = os.getcwd()
  os.system('mtzdmp ' + mtz + ' >mtzdmp_out')
  mtz_out = open(path + '/mtzdmp_out')
  flags = False
  while flags == False:
      line = mtz_out.readline()
      get_stat1 = re.compile('Column Labels')
      result_stat1 = get_stat1.search(line)
      if result_stat1:
         #print line
         next = mtz_out.readline()
         next = mtz_out.readline()
         print next
         flags = True
 
  get_stat1 = re.compile('SIGFP')
  SIGFP = get_stat1.search(next)
  if SIGFP:
     sigf = 'SIGF=SIGFP'

  get_stat1 = re.compile('SIGF')
  SIGF = get_stat1.search(next)
  if SIGF:
     sigf = 'SIGF=SIGF'


  if not SIGFP : 
    if not SIGF :
     print next
     sys.exit()

  if SIGFP and SIGF : 
     sigf = 'SIGF=SIGFP'


  get_stat1 = re.compile('\sF\W')
  F = get_stat1.search(next)
  if F:
     FP = 'F=F'

  get_stat1 = re.compile('\sFP\s')
  FP_test = get_stat1.search(next)

  if FP_test:
     FP = 'F=FP' 

  if not FP_test : 
    if not F :
     print next
     sys.exit()

  get_stat1 = re.compile('FREE')
  free_test = get_stat1.search(next)
  if free_test:
     free = 'FreeR_flag=FREE' 

  get_stat1 = re.compile('FreeR_flag')
  free_flag_test = get_stat1.search(next)
  if free_flag_test:
     free = 'FreeR_flag=FreeR_flag' 

  return next, sigf, FP, free
###################
def make_mrbump_desktop(sigf, fp, free, jobid, local_files, mtz, seq, mrbump_programs):
  path = os.getcwd()

  i_name = local_files.rstrip(',pdb')
  #os.chdir(path+ '/'+ jobid)




  mr_bump = open(path+'/'+jobid+'.sub', "w")
  mr_bump.write('#!/bin/sh\n')


  mr_bump.write('mrbump HKLIN ' +mtz + ' SEQIN ' + seq+' HKLOUT ' + 'OUT.mtz  XYZOUT OUT.pdb << eof\n'  )

  mr_bump.write('LABIN ' + sigf + ' ' + fp + ' ' + free + '\n')
  mr_bump.write('JOBID '+ jobid+ '_mrbump\n')
  mr_bump.write('MRPROGRAM ' + mrbump_programs + '\n')
  mr_bump.write('LOCALFILE ' + local_files + ' CHAIN ALL RMS 1.2\n')


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

  mr_bump.write('FIXSG True\n')
  mr_bump.write('PJOBS 1\n')

  mr_bump.write('CHECK False\n')
  mr_bump.write('LITE True\n')
  mr_bump.write('PICKLE False\n')
  mr_bump.write('TRYALL True\n')
  mr_bump.write('USEACORN False\n')
  mr_bump.write('DEBUG True\n')
  mr_bump.write('END\n')
  mr_bump.write('eof')
  mr_bump.close()

  os.system('chmod uoga=wrx '+ path + '/'+jobid+'.sub')
##RUN
  os.system(path + '/'+jobid+'.sub')

  os.chdir(path)
###################
def make_mrbump_Cluster(sigf, fp, free, jobid, local_files, mtz, seq, mrbump_programs):
  path = os.getcwd()

  i_name = local_files.rstrip(',pdb')
  #os.chdir(path+ '/'+ jobid)

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
  mr_bump.write('MDLS False\n')
  mr_bump.write('MDLU True\n')
  mr_bump.write('MDLM False\n')
  mr_bump.write('MDLP False\n')

  mr_bump.write('FIXSG True\n')
  mr_bump.write('PJOBS 1\n')

  mr_bump.write('CHECK False\n')
  mr_bump.write('LITE True\n')
  mr_bump.write('PICKLE False\n')
  mr_bump.write('TRYALL True\n')
  mr_bump.write('USEACORN False\n')
  mr_bump.write('DEBUG True\n')
  mr_bump.write('END\n')
  mr_bump.write('eof')
  mr_bump.close()

  os.system('chmod uoga=wrx '+ path + '/'+ jobid +'/mr_bump_' + jobid)
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
     print space_split[1]
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
          print  ensemble
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
       print line
       split = re.split(pattern, line)
       print split
       ASU = split[1]
   return ASU

###############################
def  make_MRBUMP_run(mtz, pdb, run_dir, fasta, name, mrbump_programs):
  
   #files to make:
   phaser_mtz = 'fail'
   phaser_pdb = 'fail'
   molrep_mtz = 'fail'
   molrep_pdb = 'fail'
   molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT = 'none'
   phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = 'none'

   #get names of flags
   next, sigf, FP, free = get_flags(mtz)  


   #make script
   #remove_hires(mtz,  curr_dir + '/' + i , i.lower()) # format and move
   
  # os.system('mkdir ' +run_dir + '/' +name)
   os.chdir(run_dir)

   make_mrbump_desktop(sigf, FP, free, name, pdb, mtz, fasta, mrbump_programs)

   
   #runshelx - beta version on refmac and phaser output:
   #Phaser output:

   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     phaser_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     phaser_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/phaser/refine/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'
     print 'got phaser'
     shelx_phaser_path = run_dir+ '/search_' +name + '_mrbump/phaser_shelx'
     os.system('mkdir ' +shelx_phaser_path)
     os.chdir(shelx_phaser_path)

     ASU = pdbcur(phaser_pdb)

     phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT, SPACE = run_shelx.RUN(phaser_mtz, phaser_pdb, ASU, fasta, shelx_phaser_path)  ####### NEEDS correct ASU, solvent content!!

   #molrep output:
   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     molrep_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     molrep_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/mr/molrep/refine/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'
     print 'got molrep'
     shelx_molrep_path = run_dir+ '/search_' +name + '_mrbump/molrep_shelx'
     os.system('mkdir ' +shelx_molrep_path)
     os.chdir(shelx_molrep_path)

     ASU = pdbcur(molrep_pdb)

     molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, SPACE = run_shelx.RUN(molrep_mtz, molrep_pdb, ASU, fasta, shelx_molrep_path)  ####### NEEDS correct ASU, solvent content!!

    
   # if the data is there, return data
   os.chdir(run_dir)
   return phaser_mtz, phaser_pdb, molrep_mtz, molrep_pdb,  molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT

     
############
if __name__ == "__main__":
  RUNNING=open('/home/jaclyn/Baker/test_cases/1AAR/TESTLOG', "w")
   
  final_ensembles=[]
  for infile in glob.glob( os.path.join('/home/jaclyn/Baker/test_cases/1AAR/ensembles', '*.pdb') ):
    final_ensembles.append(infile)
  
  mtz = '/home/jaclyn/Baker/test_cases/1AAR/1aar_unique1.mtz'
  fasta = '/home/jaclyn/Baker/test_cases/1AAR/1AAR_.fasta'
  run_dir = '/home/jaclyn/Baker/test_cases/1AAR/Mr_BUMP'
  mrbump_programs = 'molrep'
  ###
  for each_ensemble in final_ensembles:
    name = re.split('/', each_ensemble)
    name = name.pop()
    name = re.sub('.pdb', '', name)
    pdb= each_ensemble
  
  
    phaser_mtz, phaser_pdb, molrep_mtz, molrep_pdb,  molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = make_MRBUMP_run(mtz, pdb, run_dir, fasta, name, mrbump_programs)
  
    RUNNING.write(name+':\n'+
    '\nphaser done \n' 
    'phaser PDB: '+phaser_pdb+' Phaser MTZ:'+ phaser_mtz + '\n'+
    'shelx score : '+phaser_shelxscore+'\n' +
    '\nMolrep done \n' 
    'molrep PDB: '+molrep_pdb+' Phaser MTZ:'+ molrep_mtz+ '\n'+
    'shelx score : '+molrep_shelxscore+'\n'+
    '------------------------\n')
    RUNNING.flush()
















#!/usr/bin/python2.6

# run mrbump anr rebuild with shelx

import signal, time
import subprocess
import os, re, sys, glob
import run_shelx
import run_shelx_OLD
from multiprocessing import Process, Queue, JoinableQueue, Pool, Value, Array
import shutil


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
      if  re.search('Column Types', line):
         lab = mtz_out.readline()
         lab = mtz_out.readline() 

         flags = True
#  print next
#  print lab




  lab_list=re.split('\s*', lab)
  name_list=re.split('\s*', next)
  
  sigf = lab_list.index('Q')
  sigf = name_list[sigf]
  
  FP = lab_list.index('F')
  FP = name_list[FP]

  free = lab_list.index('I')
  free = name_list[free]

 # print sigf, FP, free

  return next, sigf, FP, free
###################
def make_mrbump_desktop(sigf, fp, free, jobid, local_files, mtz, seq, noASU):
  path = os.getcwd()

  i_name = local_files.rstrip(',pdb')
  #os.chdir(path+ '/'+ jobid)




  mr_bump = open(path+'/'+jobid+'.sub', "w")
  mr_bump.write('#!/bin/sh\n')


  mr_bump.write('mrbump HKLIN ' +mtz + ' SEQIN ' + seq+' HKLOUT ' + 'OUT.mtz  XYZOUT OUT.pdb << eof\n'  )

  mr_bump.write('LABIN ' + sigf + ' ' + fp + ' ' + free + '\n')
  mr_bump.write('JOBID '+ jobid+ '_mrbump\n')
  mr_bump.write('MRPROGRAM molrep phaser\n')
  mr_bump.write('LOCALFILE ' + local_files + ' CHAIN ALL RMS 0.1\n')


  mr_bump.write('SCOPSEARCH False\n')
  mr_bump.write('PQSSEARCH False\n')
  mr_bump.write('SSMSEARCH False\n')

  mr_bump.write(noASU+'\n')

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

  os.system('chmod uoga=wrx '+ path + '/'+jobid+'.sub')
##RUN
  os.system(path + '/'+jobid+'.sub  >'+jobid+'.log')
  


  os.chdir(path)
##########################
def make_mrbump_desktop_domain(sigf, fp, free, jobid, local_files, mtz, seq, fixed_pdb, FIXED_INPUT):
  path = os.getcwd()

  i_name = local_files.rstrip(',pdb')
  #os.chdir(path+ '/'+ jobid)




  mr_bump = open(path+'/'+jobid+'.sub', "w")
  mr_bump.write('#!/bin/sh\n')


  mr_bump.write('mrbump HKLIN ' +mtz + ' SEQIN ' + seq+' HKLOUT ' + 'OUT.mtz  XYZOUT OUT.pdb << eof\n'  )

  mr_bump.write('LABIN ' + sigf + ' ' + fp + ' ' + free + '\n')
  mr_bump.write('JOBID '+ jobid+ '_mrbump\n')
  mr_bump.write('MRPROGRAM molrep phaser\n')
  mr_bump.write('LOCALFILE ' + local_files + ' CHAIN ALL RMS 0.1\n')

  if FIXED_INPUT:
     mr_bump.write('FIXED_XYZIN '+fixed_pdb+' IDENTIY 0.6 \n')
  if not FIXED_INPUT:
     mr_bump.write('LOCALFILE ' +fixed_pdb + ' CHAIN ALL RMS 0.1\n')

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

  os.system('chmod uoga=wrx '+ path + '/'+jobid+'.sub')
##RUN
  os.system(path + '/'+jobid+'.sub')

  os.chdir(path)


###################
def make_mrbump_Cluster(sigf, fp, free, jobid, local_files, mtz, seq):
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
  mr_bump.write('MRPROGRAM phaser\n')
  mr_bump.write('LOCALFILE ' + local_files + ' CHAIN A RMS 1.0\n')


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
       #print split
       ASU = split[1]
   return ASU

###############################
def  make_MRBUMP_run(mtz, pdb, run_dir, fasta, name, sigf, FP, free, noASU, EarlyTerminate, NoShelx, NoShelxCycles, Resultspath, SHELX_OLD):
   SolutionFound = False 
   #files to make:
   phaser_mtz = 'fail'
   phaser_pdb = 'fail'
   molrep_mtz = 'fail'
   molrep_pdb = 'fail'
   molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT = 'none'
   phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = 'none'

   #get names of flags
   
  # next, sigf, FP, free = get_flags(mtz)  
   sigf = 'SIGF='+sigf
   FP = 'F='+FP
   free =  'FreeR_flag='+free
   #make script
   #remove_hires(mtz,  curr_dir + '/' + i , i.lower()) # format and move
 
  # os.system('mkdir ' +run_dir + '/' +name)
   os.chdir(run_dir)

   make_mrbump_desktop(sigf, FP, free, name, pdb, mtz, fasta, noASU)

   
   #runshelx - beta version on refmac and phaser output:
   #Phaser output:

   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     
     phaser_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     phaser_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'
     if NoShelx == True:
       os.system('cp '+phaser_mtz + '  '+Resultspath    )
       os.system('cp '+phaser_pdb + '  '+Resultspath    )

     if NoShelx == False:
        print '---found  phaser output, rebuilding in shelx---'
        shelx_phaser_path = run_dir+ '/search_' +name + '_mrbump/phaser_shelx'
        os.system('mkdir ' +shelx_phaser_path)
        os.chdir(shelx_phaser_path)

        ASU = pdbcur(phaser_pdb)

        if SHELX_OLD == False:
            phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT, SPACE = run_shelx.RUN(phaser_mtz, phaser_pdb, ASU, fasta, shelx_phaser_path, EarlyTerminate, NoShelxCycles, run_dir)  ####### NEEDS correct ASU, solvent content!! 
        if SHELX_OLD == True:
            phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT, SPACE = run_shelx_OLD.RUN(phaser_mtz, phaser_pdb, ASU, fasta, shelx_phaser_path, EarlyTerminate, NoShelxCycles, run_dir)     



   #molrep output:
   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     molrep_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     molrep_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'
     if NoShelx == True:  
         os.system('cp '+molrep_mtz + '  '+Resultspath    )
         os.system('cp '+molrep_pdb + '  '+Resultspath    )

     if NoShelx == False:
        print '---found molrep output, rebuilding in shelx---'
        shelx_molrep_path = run_dir+ '/search_' +name + '_mrbump/molrep_shelx'
        os.system('mkdir ' +shelx_molrep_path)
        os.chdir(shelx_molrep_path)

        ASU = pdbcur(molrep_pdb)

        molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, SPACE = run_shelx.RUN(molrep_mtz, molrep_pdb, ASU, fasta, shelx_molrep_path, EarlyTerminate, NoShelxCycles, run_dir)  ####### NEEDS correct ASU, solvent content!!

    
   # if the data is there, return data
   os.chdir(run_dir)
   return phaser_mtz, phaser_pdb, molrep_mtz, molrep_pdb,  molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT
################
def  make_MRBUMP_run_domain(mtz, pdb, run_dir, fasta, name, fixed_pdb, FIXED_INPUT):

   #files to make:
   phaser_mtz = 'fail'
   phaser_pdb = 'fail'
   molrep_mtz = 'fail'
   molrep_pdb = 'fail'
   molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT = 'none'
   phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = 'none'

   #get names of flags
   next, sigf, FP, free = get_flags(mtz)  
   os.chdir(run_dir)

   make_mrbump_desktop_domain(sigf, FP, free, name, pdb, mtz, fasta, fixed_pdb,FIXED_INPUT )


   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     phaser_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     phaser_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/phaser/refmac_phaser_loc0_ALL_'+name+'_UNMOD.pdb'
     print 'got phaser'
     shelx_phaser_path = run_dir+ '/search_' +name + '_mrbump/phaser_shelx'
     os.system('mkdir ' +shelx_phaser_path)
     os.chdir(shelx_phaser_path)

     ASU = pdbcur(phaser_pdb)

     phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT, SPACE = run_shelx.RUN(phaser_mtz, phaser_pdb, ASU, fasta, shelx_phaser_path)  ####### NEEDS correct ASU, solvent content!!

   #molrep output:
   if os.path.exists(run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz' ):
     molrep_mtz = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_HKLOUT_loc0_ALL_'+name+'_UNMOD.mtz'
     molrep_pdb = run_dir + '/search_'+name+'_mrbump/data/loc0_ALL_'+name+'/unmod/refine/molrep/refmac_molrep_loc0_ALL_'+name+'_UNMOD.pdb'
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
def split(ensembles, nproc):  #split into parallel comparisons for speed

 divide= len(ensembles) / nproc

 chunks = []

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
def run_parallel(mtz, chunk_of_ensembles, run_dir, fasta,  log_name, sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, batchname, SHELX_OLD): #loops through each chunk
  
  log = open(log_name, "w")
  inc = 1
  for each_pdb in chunk_of_ensembles:
     name = re.split('/', each_pdb)
     name = name.pop()
     name = re.sub('.pdb', '', name)
     #print name
     print '=== In batch '+str(batchname)+' running job '+str(inc)+' of '+str(len(chunk_of_ensembles))+'. '+str(len(chunk_of_ensembles)-inc)  +' left to go'
     inc +=1 
     phaser_mtz, phaser_pdb, molrep_mtz, molrep_pdb,  molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, phaser_shelxscore, phaser_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = make_MRBUMP_run(mtz, each_pdb, run_dir, fasta, name,  sigf, FP, free, noASU, EarlyTerminate, NoShelx, NoShelxCycles, Resultspath, SHELX_OLD)
     
     log.write(name+':\n'+
     '\nphaser done \n' 
     'phaser PDB: '+phaser_pdb+' Phaser MTZ:'+ phaser_mtz + '\n'+
     'shelx score : '+phaser_shelxscore+'\n' +
     '\nMolrep done \n' 
     'molrep PDB: '+molrep_pdb+' Phaser MTZ:'+ molrep_mtz+ '\n'+
     'shelx score : '+molrep_shelxscore+'\n'+
     '------------------------\n')
     log.flush()

     if isnumber(phaser_shelxscore):
       if float(phaser_shelxscore) > 25:
         #print 'Found a Solution '+ phaser_pdb
         #phaser_shelxscore =re.sub('\.','_', phaser_shelxscore)
         #phaser_refmacfreeR=re.sub('\.','_', phaser_refmacfreeR)

         os.system('cp '+ phaser_XYZOUT+'  '+ Resultspath+'/phaser'+phaser_shelxscore+'__'+phaser_refmacfreeR+'.pdb'  )
         os.system('cp '+ phaser_HKLOUT+'  '+ Resultspath+'/phaser'+phaser_shelxscore+'__'+phaser_refmacfreeR+'.mtz'  )        

         if EarlyTerminate:
           # return True
            sys.exit()
     if isnumber(molrep_shelxscore): 
      if float(molrep_shelxscore) > 25:
         #print 'Found a Solution '+ molrep_pdb
         #molrep_shelxscore =re.sub('\.','_', molrep_shelxscore)
         #molrep_refmacfreeR=re.sub('\.','_', molrep_refmacfreeR)

         os.system('cp ' +molrep_XYZOUT+'   '+ Resultspath+'/molrep'+molrep_shelxscore+'_'+molrep_refmacfreeR+'.pdb'  )
         os.system('cp ' +molrep_HKLOUT+'   '+ Resultspath+'/molrep'+molrep_shelxscore+'_'+molrep_refmacfreeR+'.mtz' )
         if EarlyTerminate:
           # return True
            sys.exit()


################################
def run_parallel_domain(mtz, chunk_of_ensembles, run_dir, fasta,  log_name,  fixed_PDB, FIXED_INPUT): #loops through each chunk
  
  log = open(log_name, "w")

  for each_pdb in chunk_of_ensembles:
     name = re.split('/', each_pdb)
     name = name.pop()
     name = re.sub('.pdb', '', name)
    # print name

     phaser_mtz, phaser_pdb, molrep_mtz, molrep_pdb,  molrep_shelxscore, molrep_refmacfreeR, molrep_HKLOUT, molrep_XYZOUT, phaser_shelxscore, molrep_refmacfreeR, phaser_HKLOUT, phaser_XYZOUT = make_MRBUMP_run_domain(mtz, each_pdb, run_dir, fasta, name, fixed_PDB, FIXED_INPUT)
   
     log.write(name+':\n'+
     '\nphaser done \n' 
     'phaser PDB: '+phaser_pdb+' Phaser MTZ:'+ phaser_mtz + '\n'+
     'shelx score : '+phaser_shelxscore+'\n' +
     '\nMolrep done \n' 
     'molrep PDB: '+molrep_pdb+' Phaser MTZ:'+ molrep_mtz+ '\n'+
     'shelx score : '+molrep_shelxscore+'\n'+
     '------------------------\n')
     log.flush()


################################
def split_into_runs(mtz, ensembles, run_dir, fasta,  nProc, sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, SHELX_OLD):
      cur_dir=os.getcwd()
      inc = 1
      threads = []
      print 'the ensembles will be divided into '+str(nProc) +' batches, each job will contain '+str(len(ensembles[0])) +' ensembles'     
      for each_chunk in ensembles:  

        #print '===now running job '+str(inc) +' containing  '+str(len(each_chunk)) +'ensembles. '+str(len(ensembles)-inc) +' jobs left to do==='
        log_name = cur_dir+'/LOG_proc'+str(inc)                        
        thread = Process(target=run_parallel, args=(mtz, each_chunk, run_dir, fasta,  log_name,  sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles, inc, SHELX_OLD))        
        thread.start() 
        threads.append(thread)
        inc+=1

      for thread in threads:

        thread.join()
#############################################
def split_into_runs_domains(mtz, ensembles, run_dir, fasta,  nProc, fixed_PDB, FIXED_INPUT):
      cur_dir=os.getcwd()
      inc = 1
      threads = []

      for each_chunk in ensembles:  
        log_name = cur_dir+'/LOG_proc'+str(inc)                        
        
        thread = Process(target=run_parallel_domain, args=(mtz, each_chunk, run_dir, fasta,  log_name, fixed_PDB, FIXED_INPUT))        
        thread.start() 
        threads.append(thread)
        inc+=1

      for thread in threads:

        thread.join()



###########################
#def START_parallel(mtz, fasta, )

if __name__ == "__main__":
 ensembles=[]
 for infile in glob.glob( os.path.join('/home/jaclyn/Ample_tests/toxd-example/ensembles_1', '*.pdb') ):
  print infile
  ensembles.append(infile)

 mtz = '/home/jaclyn/Ample_tests/toxd-example/1dtx.mtz'
 fasta = '/home/jaclyn/Ample_tests/toxd-example/toxd_.fasta'
 run_dir = '/home/jaclyn/Ample_tests/toxd-example/Ro/MRBUMP'
 nproc = 1
 #'fixed_pdb='/home/jaclyn/DOMAINS/ample/TEST/1al6/Known.pdb'
 sigf = 'SIGFP'
 FP = 'FP'
 free = 'FreeR_flag'
 EarlyTerminate = False
 Resultspath = '/home/jaclyn/Ample_tests/toxd-example/Ro'
 NoShelx = False
 NoShelxCycles = 1
 noASU = '1'

 chunk = split(ensembles, nproc)
 #'split_into_runs_domains(mtz, chunk, run_dir, fasta, nproc, fixed_pdb)
 split_into_runs(mtz,chunk  , run_dir, fasta,  nproc, sigf, FP, free, noASU, EarlyTerminate, Resultspath, NoShelx, NoShelxCycles)
###s


















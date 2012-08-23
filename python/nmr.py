#!/usr/bin/env python

#

import re
import os, glob
import sys
import subprocess
import time
import operator
import shutil, stat


def comparative_model(homolog, fasta, MR_ROSETTA, ROSETTA_DB, ALI, a9mers, a3mers, name, Nprocess):

   curdir=os.getcwd()
   
   s = open(curdir+'/remodel.sh', "w")
   s.write(MR_ROSETTA +' \\\n'+
           '-database '+ROSETTA_DB+' \\\n'+
           '-MR:mode cm \\\n'+
           '-in:file:extended_pose 1 \\\n'+
           '-in:file:fasta '+fasta+' \\\n'+
           '-in:file:alignment '+ALI+' \\\n'+
           '-in:file:template_pdb '+homolog+' \\\n'+
           '-loops:frag_sizes 9 3 1 \\\n'+
           '-loops:frag_files '+a9mers+' '+a3mers+' none \\\n'+
           '-loops:random_order \\\n'+
           '-loops:random_grow_loops_by 5 \\\n'+
           '-loops:extended \\\n'+
           '-loops:remodel quick_ccd \\\n'+
           '-loops:relax relax \\\n'+
           '-relax:default_repeats 4 \\\n'+
           '-relax:jump_move true    \\\n'+
           '-cm:aln_format grishin \\\n'+
           '-MR:max_gaplength_to_model 8 \\\n'+
           '-nstruct '+str(Nprocess)+' \\\n'+
           '-ignore_unrecognized_res \\\n'+
           '-overwrite \\\n' +
           '&')

   s.close()
   os.system('sh ' + curdir+'/remodel.sh >rosetta.log  &')
   
   
################################
def minimise(homolog, fasta, MINIMIZE,RELAX, ROSETTA_DB, ALI, a9mers, a3mers, name, Nprocess):
   print 'here'
   curdir=os.getcwd()
   #try 2 types of minimisation
   s = open(curdir+'/remodel.sh', "w")
   s.write(RELAX+' \\\n'+
       '-database '+ROSETTA_DB+' \\\n'+
       '-in:file:s '+homolog+' \\\n'+
       '-in:file:fullatom True \\\n'+
       '-ddg:out_pdb_prefix S \\\n'+
       '-out:nstruct '+ str(Nprocess))
   
   s.close()
   os.system('sh ' + curdir+'/remodel.sh')
   models=[]
   x=1
   while x< NProcess+1:
     if os.path.exists(curdir+'/S_'+name.upper()+'_000'+str(x)+'.pdb'):
      print 'got ',curdir+'/S_'+name.upper()+'_000'+str(x)+'.pdb'     
      models.append(curdir+'/S_'+name.upper()+'_000'+str(x)+'.pdb' )
     x+=1
   sys.exit()
   return  models


########################
def idealise(pdb, IDEALIZE, ROSETTA_DB  ):
   cmd = IDEALIZE+' -database '+ROSETTA_DB+'  -s  '+pdb 
   

   os.system(cmd + ' >idealise.log' )
  # p = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE,
 #                                stdout = subprocess.PIPE, stderr=subprocess.PIPE)
  # p.wait()
   curdir=os.getcwd()

   name = os.path.split(pdb)
  
   name = re.sub('.pdb','', name[-1])
   

   if os.path.exists( curdir+'/'+name+'_0001.pdb'  ):
       return curdir+'/'+name+'_0001.pdb', name+'_0001'
   else:
      print 'idealization failed'
      sys.exit()

########################
def model_with_map():
   print 'here'

########################
def Rosetta_refine():
   print 'here'
########################
def align(homolog_seq, fasta, hhsearch, name):
   print 'here'
   curdir=os.getcwd()
  
   c = open(curdir+'/cat',"w")
   f=open(fasta)
   for line in f:
      c.write(line)
      c.flush()
   t = open(homolog_seq)
   for line in t:
      c.write(line)
      c.flush()
   c.close()


   cmd = 'muscle  -in '+curdir+'/cat -out '+curdir+'/muscle'
   print cmd

   p = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.PIPE)
   p.wait()

   m = open(curdir+'/muscle')
   ali = open(curdir+'/ali', "w")
   seq = ''
   seqs=[]
   for line in m:
       if re.search('>', line):
          IN = True
          seqs.append(seq)
          seq=''
       else:
          seq+=line.rstrip('\n')
          print line
   seqs.append(seq)
   print  seqs[1]
   print seqs[2]

   ali.write('## TARGET '+name+'\n'+
             '# hhsearch\n'+
              'scores_from_program: 0 1.00\n'+
              '0 ' + seqs[1]+'\n'+
              '0 ' + seqs[2])
   ali.close
   return curdir+'/ali'
########################
def get_sequence(pdb, name):

 cur_dir=os.getcwd()
 curscript=os.path.join(cur_dir, "cur")
 cur=open(curscript, "w")  


 cur.write('#!/bin/sh\n'+
     'pdbset xyzin '+pdb+ '<<EOF\n'+
     'sequence single \n'+
     'EOF')
 cur.close()
 os.chmod(curscript, stat.S_IRWXU)
 os.system(curscript + '>cursum')

 s = open(cur_dir+'/SEQUENCE')
 o = open(cur_dir+'/'+name, "w")
 for line in s:
    if re.search('>', line):
      o.write(line)
    else:
      l = re.sub('\s*', '', line)
      o.write(l)
 o.write('\n')
 o.close()
 return cur_dir+'/'+name
#####################
def RUN_MINIMISE(homolog, NProcess):
 curdir = os.getcwd()
 os.mkdir(curdir+'/RUN')
 os.chdir(curdir+'/RUN')
 
# ROSETTA_OVER_PATH = os.environ.get("ROSETTA_PATH")
# ROSETTA_OVER_PATH = '/home/jaclyn/programs/rosetta3.1_Bundles' 
 ROSETTA_OVER_PATH = '/home/jaclyn/programs/rosetta3.3_bundles'
 if os.path.exists(ROSETTA_OVER_PATH):
   ROSETTA_PATH               =ROSETTA_OVER_PATH+'/rosetta_source/bin/AbinitioRelax.linuxgccrelease'
   ROSETTA_cluster            =ROSETTA_OVER_PATH+'/rosetta_source/bin/cluster.linuxgccrelease'
   ROSETTA_DB                 =ROSETTA_OVER_PATH+'/rosetta_database'
   Make_fragents_exe          =ROSETTA_OVER_PATH+'/rosetta_fragments/nnmake/make_fragments.pl'
   MR_ROSETTA                 =ROSETTA_OVER_PATH+'/rosetta_source/bin/mr_protocols.default.linuxgccrelease'
#   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize.linuxgccrelease'
   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize_jd2.default.linuxgccrelease'
   RELAX                      =ROSETTA_OVER_PATH+'/rosetta_source/bin/relax.linuxgccrelease'
   MINIMIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/minimize_with_cst.default.linuxgccrelease'

 #/home/jaclyn/programs/rosetta-3.2/rosetta_source/bin/mr_protocols.default.linuxgccrelease
 hhsearch = '/home/jaclyn/programs/hhpred/hhsuite-2.0.12-linux-x86_64/bin/hhsearch'

 #homolog = '/home/jaclyn/BEE/new_case/bevan_example/hom/1K3K.pdb'
 homolog, name = idealise(homolog, IDEALIZE, ROSETTA_DB )

 if os.path.exists(ROSETTA_OVER_PATH):
  fasta = '/media/Elements/160_300_set/myoglobin/fasta'
  homolog_seq = get_sequence(homolog, 'homolog.fasta')
  a9mers =''
  a3mers =''
  ALI =  align(homolog_seq, fasta, hhsearch, name)


 print name

 
 remodeled = minimise(homolog, fasta, MINIMIZE,RELAX , ROSETTA_DB , ALI, a9mers, a3mers, name, NProcess) 
 for outfile in remodeled:
    os.system('cp ' + outfile+' ' +curdir)

 os.chdir(curdir)
 os.system('rm -r ' +curdir+'/RUN')

######################
def MAFFT(homolog_seq, fasta,  name):

   curdir=os.getcwd()

   c = open(curdir+'/cat',"w")
   f=open(fasta)
   for line in f:
      c.write(line)
      c.flush()
   t = open(homolog_seq)
   for line in t:
      c.write(line)
      c.flush()
   c.close()

   cmd =  'mafft --maxiterate 1000 --localpair  --quiet  '+curdir+'/cat  >muscle'
   
  

   p = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE,
                                  stdout = subprocess.PIPE, stderr=subprocess.PIPE)
   p.wait()

   m = open(curdir+'/muscle')
   ali = open(curdir+'/ali', "w")
   seq = ''
   seqs=[]
   for line in m:
       if re.search('>', line):
          IN = True
          seqs.append(seq)
          seq=''
       else:
          seq+=line.rstrip('\n')
         
   seqs.append(seq)
   print 'using the alignment:' 
   print  seqs[1]
   print seqs[2]
   print 'if you want to use a different alignment, import using -alignment_file' 

   ali.write('## TARGET '+name+'\n'+
             '# hhsearch\n'+
              'scores_from_program: 0 1.00\n'+
              '0 ' + seqs[1]+'\n'+
              '0 ' + seqs[2])
   ali.close
   return curdir+'/ali'
######################
def CLUSTER_RUN_FORMAT_HOMS(homolog, NProcess, fasta, ROSETTA_PATH, a9mers, a3mers, NProc, homname, alignment_file, Models_dir ):
 curdir = os.getcwd()
 os.mkdir(curdir+'/RUN_'+homname)
 os.chdir(curdir+'/RUN_'+homname)
 RunDir = curdir+'/RUN_'+homname

 ROSETTA_OVER_PATH = ROSETTA_PATH
# ROSETTA_OVER_PATH = '/home/jaclyn/programs/rosetta3.3_bundles'
 if os.path.exists(ROSETTA_OVER_PATH):
   ROSETTA_PATH               =ROSETTA_OVER_PATH+'/rosetta_source/bin/AbinitioRelax.linuxgccrelease'
   ROSETTA_cluster            =ROSETTA_OVER_PATH+'/rosetta_source/bin/cluster.linuxgccrelease'
   ROSETTA_DB                 =ROSETTA_OVER_PATH+'/rosetta_database'
   Make_fragents_exe          =ROSETTA_OVER_PATH+'/rosetta_fragments/nnmake/make_fragments.pl'
   MR_ROSETTA                 =ROSETTA_OVER_PATH+'/rosetta_source/bin/mr_protocols.default.linuxgccrelease'
#   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize.linuxgccrelease'
   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize_jd2.default.linuxgccrelease'


 #/home/jaclyn/programs/rosetta-3.2/rosetta_source/bin/mr_protocols.default.linuxgccrelease

 homolog, name = idealise(homolog, IDEALIZE, ROSETTA_DB )



 homolog_seq = get_sequence(homolog, 'homolog.fasta')
 if os.path.exists(alignment_file):
     ALI = alignment_file
 if not os.path.exists(alignment_file):
     ALI =  MAFFT(homolog_seq, fasta,  name)
 return homolog,  ALI


######################
def RUN_FORMAT_HOMS(homolog, NProcess, fasta, ROSETTA_PATH, a9mers, a3mers, NProc, homname, alignment_file, Models_dir ):
 
 curdir = os.getcwd()
 os.mkdir(curdir+'/RUN_'+homname)
 os.chdir(curdir+'/RUN_'+homname)
 RunDir = curdir+'/RUN_'+homname
 
 ROSETTA_OVER_PATH = ROSETTA_PATH
# ROSETTA_OVER_PATH = '/home/jaclyn/programs/rosetta3.3_bundles'
 if os.path.exists(ROSETTA_OVER_PATH):
   ROSETTA_PATH               =ROSETTA_OVER_PATH+'/rosetta_source/bin/AbinitioRelax.linuxgccrelease'
   ROSETTA_cluster            =ROSETTA_OVER_PATH+'/rosetta_source/bin/cluster.linuxgccrelease'
   ROSETTA_DB                 =ROSETTA_OVER_PATH+'/rosetta_database'
   Make_fragents_exe          =ROSETTA_OVER_PATH+'/rosetta_fragments/nnmake/make_fragments.pl'
   MR_ROSETTA                 =ROSETTA_OVER_PATH+'/rosetta_source/bin/mr_protocols.default.linuxgccrelease'
#   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize.linuxgccrelease'
   IDEALIZE                   =ROSETTA_OVER_PATH+'/rosetta_source/bin/idealize_jd2.default.linuxgccrelease'


 #/home/jaclyn/programs/rosetta-3.2/rosetta_source/bin/mr_protocols.default.linuxgccrelease

 homolog, name = idealise(homolog, IDEALIZE, ROSETTA_DB )

 
  
 homolog_seq = get_sequence(homolog, 'homolog.fasta')
 if os.path.exists(alignment_file):
     ALI = alignment_file
 if not os.path.exists(alignment_file):   
     ALI =  MAFFT(homolog_seq, fasta,  name) 
 
 NMODELS = NProcess
 print 'splitting ', NProcess, ' jobs  onto  ',  NProc, ' processors:'
 split_jobs =  NProcess / NProc   ### split jobs between processors
 remainder =   NProcess % NProc
 jobs = [0]
 proc = 1
 while proc < NProc +1:
        jobs.insert(proc, split_jobs)
        proc +=1
 jobs[-1] = jobs[-1] + remainder


 proc = 1
 while proc < NProc +1:
      os.system('mkdir '+ RunDir + '/models_'+str(proc))
      os.chdir(RunDir + '/models_'+str(proc))



 
      comparative_model(homolog, fasta, MR_ROSETTA, ROSETTA_DB , ALI, a9mers, a3mers, name, str(jobs[proc])) 
      proc +=1

 no_models_to_make = NMODELS
 no_models_have = 0
 finished_models = 0
 while no_models_have !=  no_models_to_make:
        no_models_have = 0
        proc = 1
        while proc < NProc +1:
          list_of_files = [file for file in os.listdir(RunDir +'/models_'+str(proc) ) if file.lower().endswith('.pdb')]


          no_models_have  += len(list_of_files)
          proc+=1
        if no_models_have > finished_models:
          
           print 'Number of models made so far = ' + str(no_models_have)
          
           divisor = no_models_to_make/10
           if divisor !=0:
                if no_models_have !=0:
                  if no_models_have%divisor == 0:
                     print 'Number of models made so far = ' + str(no_models_have)

           finished_models = no_models_have
        if no_models_have ==  no_models_to_make:
           break
        time.sleep(5)
 print 'MODELLING DONE'


 proc = 1
 cat_string = ''
 while proc < NProc +1:
         
         
         
         for file in glob.glob( os.path.join(RunDir + '/models_'+str(proc), '*.pdb') ):

                 name = re.split('/', file)
                 pdbname = str(name.pop())
                 shutil.copyfile(RunDir + '/models_'+str(proc) +'/'+pdbname,Models_dir+'/' +str(proc)+'_'+pdbname )

         cat_string += RunDir + '/OUT_' + str(proc) +'  '
         #os.system('rm -r '+RunDir + '/models_'+str(proc) )
         proc+=1

 os.chdir(curdir)
 




 




if __name__ == "__main__":
 

 NProcess =4 
 homolog = '/data2/jac45/ccp4-6.3.0/share/ample/nmr/models'
 ROSETTA_PATH ='/home/jac45/rosetta'
 fasta = '/data2/jac45/tox/toxd-example/toxd_.fasta'
 a9mers = '/data2/jac45/tox/toxd-example/aat000_09_05.200_v1_3'
 a3mers = '/data2/jac45/tox/toxd-example/aat000_03_05.200_v1_3' 
 NProc = 4 
 alignment_file = ''
 for hom in os.listdir(homolog):
  pdbs =  re.split('\.', hom)
  print pdbs
  if pdbs[-1] == 'pdb':   

      name = hom  
      print hom
      RUN_FORMAT_HOMS(homolog+'/'+hom, NProcess, fasta, ROSETTA_PATH, a9mers, a3mers, NProc, name, alignment_file)
 # RUN_MINIMISE(modelspath+'/'+hom, NProcess)


#if NMR process models first


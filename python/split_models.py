#!/usr/bin/env python

#script to truncate the non secondary structure ends of a pdb

import re
import os, glob
import sys
import subprocess 
import time
import shutil
#######  EDIT 
def fix(pdb):
  cur = open(os.getcwd()+'/cur', 'w')
  cur.write('pdbset xyzin '+ pdb+' xyzout '+pdb+' <<EOF'+
             'OCCUPANCY 1 \n'+
             'CHAIN A')
  cur.close()
####### 
def check(pdb):
  i=1
  for line in open(pdb):
    if re.search('^ATOM', line):
       if   line[13:16] =='CA ':
          i+=1
  return i        
####### 
def try_NMR(pdb):
  list = [0]
  i = 2
  ISNMR = False
  out = open(os.getcwd()+'/test.pdb', 'w')
  out.write('MODEL    1\n')
  for line in  open(pdb):
  
    
    if re.search('^ATOM', line):
       out.write(line)
  
       
       if int(line[7:11]) < list[-1]:
          ISNMR = True
          out.write('ENDMDL\nMODEL    '+str(i)+'\n')
          i+=1
       list.append( int(line[7:11]) )
  out.close()
  print ISNMR  
  return ISNMR, os.getcwd()+'/test.pdb'   


####### 
def split(model, path ):
 

 SPLIT = False
 for line in open(model):
     if  re.search('^MODEL', line):
        SPLIT = True 
        modno = splitNMR(model, path)
        return modno
     
 if SPLIT == False:
    n = os.path.split(model)
    # check if run on nmr
    ISNMR, temp  = try_NMR(model)
    if ISNMR == False:
       print 'Only remodelingin one model'   
       shutil.copyfile(model, path+'/'+n[-1])  
       return 1   
    if ISNMR == True:
       model = temp
       modno = splitNMR(model, path)
       return modno
def splitNMR(model, path): 

 name = '.pdb'
 modelin= open(model)
 modno = 0
 condition= 0
 for line in modelin:
    # print line
     model_pattern = re.compile('^MODEL')
     model_result = model_pattern.match(line)
     if model_result:
       condition = 1      
     #  print 'IN'
     pdb_pattern = re.compile('ENDMDL')
     pdb_result = pdb_pattern.match(line)
     if pdb_result:           
       condition = 0
       modno = modno+1

     if condition == 1 and not model_result :
      # print 'IN'
       new_model = open(path+'/'+str(modno)+'_'+name, "a")
       #print line[21:22]
       if line[21:22] == 'A':
          #print line[21:22] 
          new_model.write(line)
          new_model.close()

 lengths = []
 for model in os.listdir(path):
  fix(path+'/'+model)
  l =  check(path+'/'+model)
  lengths.append(l)
 if len(lengths) >1:
    mina  = min(lengths, key = int)
    maxa  = max(lengths, key = int)
    if mina != maxa:
      print 'min length = ', mina, ' max length = ', maxa
      print 'All of the models need to be the same length, edit them and try again' 
      sys.exit() 
 return modno

#!/usr/bin/python2.6

#edit the sidechains to make polyala, all and reliable

import re
import os, glob
import sys




########################################## ADD ALL

def edit_sidechains(each_file, outpath):

 if os.path.exists(each_file):
 #   print 'Found ',each_file
    my_infile = open (each_file)
   
    name = re.split('/', each_file)
    pdbname = str(name.pop()) 
    my_outfile = open (outpath + 'SCWRL_Reliable_sidechains_' +pdbname, "w")
   
    os.system('cp ' + each_file +' '+outpath + 'All_atom_' +pdbname)
    for pdbline in my_infile:
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)
        if pdb_result:
         pdb_result2 = re.split(pdb_pattern, pdbline )
         #print pdb_result2
         if pdb_result2[3] == 'MET' or pdb_result2[3] == 'ASP'  or pdb_result2[3] == 'PRO'  or pdb_result2[3] == 'GLN' or pdb_result2[3] == 'LYS' or pdb_result2[3] == 'ARG' or pdb_result2[3] == 'GLU' or pdb_result2[3] == 'SER':
          if pdb_result2[2] == 'N' or pdb_result2[2] == 'CA' or pdb_result2[2] == 'C' or pdb_result2[2] == 'O' or pdb_result2[2] == 'CB':
           #print pdb_result2[2] +' ' +  pdb_result2[3]
           my_outfile.write(pdbline)
         else:
          my_outfile.write(pdbline)
        else:
          my_outfile.write(pdbline)     
   
    my_outfile.close()
    my_infile.close()
   
    my_infile = open (each_file)
    my_outfile2 = open (outpath + 'poly_ala_' +pdbname, "w")
    for pdbline in my_infile:
        #print pdbline
        pdb_pattern = re.compile('^ATOM\s*(\d*)\s*(\w*)\s*(\w*)\s*(\w)\s*(\d*)\s')
        pdb_result = pdb_pattern.match(pdbline)
        if pdb_result:
         pdb_result2 = re.split(pdb_pattern, pdbline )
        # print pdb_result2
         if pdb_result2[3] !='':# == 'MET' or pdb_result2[3] == 'ASP'  or pdb_result2[3] == 'PRO'  or pdb_result2[3] == 'GLN' or pdb_result2[3] == 'LYS' or pdb_result2[3] == 'ARG' or pdb_result2[3] == 'GLU' or pdb_result2[3] == 'SER':
          if pdb_result2[2] == 'N' or pdb_result2[2] == 'CA' or pdb_result2[2] == 'C' or pdb_result2[2] == 'O' or pdb_result2[2] == 'CB':
          # print pdb_result2[2] +' ' +  pdb_result2[3]
           my_outfile2.write(pdbline)
         else:
          my_outfile2.write(pdbline)
        else:
          my_outfile2.write(pdbline)   
   
    my_outfile2.close()
 else:
    print 'Cant Find: ',each_file
   

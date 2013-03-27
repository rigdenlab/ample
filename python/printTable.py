#!/usr/bin/env python

import locale
import os, sys, re
locale.setlocale(locale.LC_NUMERIC, "")


class Table:

   def __init__(self):   
      self.bumppath = ''
      self.cluster = ''

   def format_num(self, num):
       """Format a number according to given places.
       Adds commas, etc. Will truncate floats into ints!"""
   
       try:
           if "." in num:
             inum = float(num)
             return locale.format("%.2f", (0, inum), True)
           else:
             inum = int(num)
             return locale.format("%.*f", (0, inum), True)
   
       except (ValueError, TypeError):
           return str(num)
   
   
   def get_max_width(self, table, index):
       """Get the maximum width of the given column index"""
       return max([len(self.format_num(row[index])) for row in table])
   
   
   def pprint_table(self, out, table):
       """Prints out a table of data, padded for alignment
       @param out: Output stream (file-like object)
       @param table: The table to print. A list of lists.
       Each row must have the same number of columns. """
       col_paddings = []
   
       for i in range(len(table[0])):
           col_paddings.append(self.get_max_width(table, i))
   
       for row in table:
           # left col
           print >> out, row[0].ljust(col_paddings[0] + 1),
           # rest of the cols
           for i in range(1, len(row)):
               col = self.format_num(row[i]).rjust(col_paddings[i] + 2)
               print >> out, col,
           print >> out


   def Xmaketable(self):
       """
       DEPRECATED - we just use the resultsTable.dat files
       Go through the MRBUMP directories and parse the output files to get the results
       """
       table = []
       
       if not self.cluster:
           for run in os.listdir(self.bumppath):
               if re.search('search', run): 
                  # print run 
                   name = re.split('_', run)
                   sidechains = name[1]
                   if sidechains == 'SCWRL':
                       trunc = name[5]
                       rad = name[7]
                   else:
                       trunc = name[4]
                       rad = name[6] 
                   trunc = str(round(float(trunc)))
                   #phaser
                   phasersucc = ''
                   phaserCC = 'none'
                   phaserFreeR = 'none'
                   if os.path.exists(os.path.join(self.bumppath, run, 'phaser_shelx' ) ):
                         if os.path.exists(    os.path.join(self.bumppath, run, 'phaser_shelx', 'orig.pdb' )):
                                for line in  open(os.path.join(self.bumppath, run, 'phaser_shelx', 'orig.pdb' )):
                                           if re.search('CC', line):     
                                                 split = re.split('CC\s*=\s*(\d*\.\d*)',line)
                                                 score = split[1]
                                                 phaserCC = score  

                         if os.path.exists(    os.path.join(self.bumppath, run, 'phaser_shelx', 'refined.pdb' )):  
                                for line in  open(os.path.join(self.bumppath, run, 'phaser_shelx', 'refined.pdb' )):
                                           if re.search('REMARK\s*(\d*)\s*FREE R VALUE\s*\:', line):     
                                                 split = re.split('REMARK\s*(\d*)\s*FREE R VALUE\s*\:\s*(\d*\.\d*)',line)
                                                 score = split[2]
                                                 phaserFreeR = score  

                      
                   #print phaserCC 
                   if phaserCC != 'none': 
                      if float(phaserCC) > 25:
                         phasersucc = 'Solution'
                      else:
                         phasersucc = ''
                      table.append([phasersucc,'',trunc, sidechains, 'Phaser',phaserFreeR,float(phaserCC)])




                   #molrep
                   molrepsucc = ''
                   molrepCC = 'none'
                   molrepFreeR = 'none'
                   if os.path.exists(os.path.join(self.bumppath, run, 'molrep_shelx' ) ):
                         if os.path.exists(    os.path.join(self.bumppath, run, 'molrep_shelx', 'orig.pdb' )):
                                for line in  open(os.path.join(self.bumppath, run, 'molrep_shelx', 'orig.pdb' )):
                                           if re.search('CC', line):     
                                                 split = re.split('CC\s*=\s*(\d*\.\d*)',line)
                                                 score = split[1]
                                                 molrepCC = score  

                         if os.path.exists(    os.path.join(self.bumppath, run, 'molrep_shelx', 'refined.pdb' )):  
                                for line in  open(os.path.join(self.bumppath, run, 'molrep_shelx', 'refined.pdb' )):
                                           if re.search('REMARK\s*(\d*)\s*FREE R VALUE\s*\:', line):     
                                                 split = re.split('REMARK\s*(\d*)\s*FREE R VALUE\s*\:\s*(\d*\.\d*)',line)
                                                 score = split[2]
                                                 molrepFreeR = score   
                      
                   #print molrepCC 
                   if molrepCC  != 'none':
                      if float(molrepCC) > 25:
                         molrepsucc = 'Solution'
                      else:
                         molrepsucc = ''
                      table.append([molrepsucc,'',trunc, sidechains, 'Molrep',molrepFreeR,float(molrepCC)])



       if self.cluster:
           table = []
           for run in os.listdir(self.bumppath):
            if os.path.isdir(os.path.join(self.bumppath, run  )):
             #phaser
             #print os.path.join(self.bumppath, run  )
             for search in os.listdir(os.path.join(self.bumppath, run  )):
                   if re.search('search' , search):
                         phaser = ''
                         molrep = ''
                         phaserCC ='none'
                         molrepCC ='none'
                         phasersucc = ''
                         phaserFreeR = ''
                         molrepFreeR = ''
                         if os.path.exists(os.path.join(self.bumppath, run , search, 'phaser_shelx' )):
                            print 'got' , os.path.join(self.bumppath, run , search, 'phaser_shelx' )
                            for pdb in os.listdir(os.path.join(self.bumppath, run , search, 'phaser_shelx' )):
                                     if re.search('refmac_phaser', pdb):
                                         print pdb  

                                         name = re.split('_', pdb)
                                         print name                  
                                         sidechains = name[4]
                                         if sidechains == 'SCWRL':
                                               trunc = name[8]
                                         else:
                                               trunc = name[7]
                                         print trunc
                                         trunc = str(round(float(trunc)))
                                         if os.path.exists(os.path.join(self.bumppath, run , search, 'phaser_shelx', 'orig.pdb' )):
                                                 
                                             for line in  open(os.path.join(self.bumppath, run,search, 'phaser_shelx', 'orig.pdb' )):
                                                 if re.search('CC', line):
                                                     split = re.split('CC\s*=\s*(\d*\.\d*)',line)
                                                     score = split[1]
                                                     phaserCC = score
                                         if phaserCC != 'none':
                                            if float(phaserCC) > 25:
                                                phasersucc = 'Solution'
                                            else:
                                                phasersucc = ''
                                            table.append([phasersucc,'',trunc, sidechains, 'Phaser',phaserFreeR,float(phaserCC)])  


                                         if os.path.exists(os.path.join(self.bumppath, run , search, 'molrep_shelx', 'orig.pdb' )):
                                                 
                                             for line in  open(os.path.join(self.bumppath, run, search, 'molrep_shelx', 'orig.pdb' )):
                                                 if re.search('CC', line):
                                                     split = re.split('CC\s*=\s*(\d*\.\d*)',line)
                                                     score = split[1]
                                                     molrepCC = score

                                         if molrepCC != 'none':
                                            if float(molrepCC) > 25:
                                                molrepsucc = 'Solution'
                                            else:
                                                molrepsucc = ''
                                            table.append([molrepsucc,'',trunc, sidechains, 'Molrep',molrepFreeR,float(molrepCC)])

    


       table.sort(key=lambda x: x[6])
       table = table[::-1]
       x=1
       for i in table:
          # print i
           i[1] = x
           x+=1
       table.insert(0,["", "Ensemble No", "Truncation level", "Side chain type", "MR program", "Final Rfree", "Shelxe CC"] )
       return table
        









if __name__ == "__main__":
    T=Table()
    T.bumppath = '/home/jmht/ample-dev1/examples/toxd-example/ROSETTA_MR_3/MRBUMP/cluster_run1'
    T.cluster = True
    table = T.maketable()
    out = sys.stdout
    T.pprint_table(out, table)  
    print 'done'    

#    table = [["", "Ensemble No", "Truncation level", "Side chain type", "MR program", "Final Rfree", "Shelxe CC"],
#        ["spam", 1, 2, "All", "Phaser", 0.25, 23.51],
#        ["eggs", 2, 6, "Reliable", "Molrep", 0.45, 4.56],
#        ["lumberjacks", 3, 10, "Poly", "Phaser", 0.53, 12.12]]
#    import sys
#    out = sys.stdout
#
#    T.pprint_table(out, table)

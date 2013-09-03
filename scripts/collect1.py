#!/usr/bin/env python

import string,sys,os,shutil

MRPROG=sys.argv[1]

def getPhaserResults(logfile):

   file=open(logfile)
   #print os.path.join(os.getcwd(), logfile)

   CAPTURE=False
   solline=""
   phaserLLG=0.0
   phaserTFZ=0.0
   line=file.readline()
   while line:
      if CAPTURE:
         if "SOLU SPAC" in line:
            CAPTURE=False
         else:
            solline+=string.strip(line)+" "
      if  "Solution #1 annotation (history):" in line or  "Solution annotation (history):" in line:
         CAPTURE=True
      line=file.readline()
   file.close()

   list=string.split(solline)
   list.reverse()
   for i in list:
      if "TFZ==" in i and "*" not in i:
         phaserTFZ=float(i.replace("TFZ==",""))
         break
      if "TFZ=" in i and "TFZ==" not in i and "*" not in i:
         phaserTFZ=float(i.replace("TFZ=",""))
         break

   for i in list:
      if "LLG==" in i:
         phaserLLG=float(i.replace("LLG==",""))
         break
      if "LLG=" in i and "LLG==" not in i:
         phaserLLG=float(i.replace("LLG=",""))
         break

   return phaserLLG, phaserTFZ

def getRefmacResults(logfile):

   file=open(logfile)

   CAPTURE=False
   initRfree=1.0
   finalRfree=1.0
   initRfactor=1.0
   finalRfactor=1.0
   noResModel=0
   noChainsModel=0
   line=file.readline()
   while line:
      if "Number of residues :" in line:
         noResModel=int(string.split(line)[-1])
      if "Number of chains   :" in line:
         noChainsModel=int(string.split(line)[-1])
      if CAPTURE:
         if "R free" in line:
            initRfree=float(string.split(line)[-2])
            finalRfree=float(string.split(line)[-1])
            CAPTURE=False
         if "R factor" in line:
            initRfactor=float(string.split(line)[-2])
            finalRfactor=float(string.split(line)[-1])
      if " $TEXT:Result: $$ Final results $$" in line:
         CAPTURE=True
      line=file.readline()
   file.close()

   return initRfree, finalRfree, initRfactor, finalRfactor, noResModel, noChainsModel

def getMrbumpResults(logfile):

   file=open(logfile)

   noResTarget=0
   noChainsTarget=0
   resolution=0.0
   line=file.readline()
   while line:
      if "Number of residues:" in line:
         noResTarget=int(string.split(line)[-1])
      if "Estimated number of molecules to search for in a.s.u.:" in line:
         noChainsTarget=int(string.split(line)[-1])
      if "Resolution of collected data (angstroms):" in line:
         resolution=float(string.split(line)[-1])  
      line=file.readline()
   file.close()

   return noResTarget, noChainsTarget, resolution


class result:

   def __init__(self):
      self.modelNAME=""
      self.CCscore=""
      self.AvgChainLen=""
      self.PDBID=""

DIR=[]
quarkDIR=os.getcwd()
SOLUTIONCOUNT=0
SOL_LIST=[]
MARGINALCOUNT=0
MARG_LIST=[]
AVERAGECOUNT=0
FULLSOLCOUNT=0
TOTALGOOD=0
TOTALPROCESSED=0

for i in os.listdir(os.getcwd()):
   if os.path.isdir(i):
      DIR.append(i)

for j in DIR:
   FOUND=False
   if os.path.isdir(os.path.join(j, "cluster_1", "ROSETTA_MR_1", "MRBUMP", "cluster_1")):
      os.chdir(os.path.join(j, "cluster_1", "ROSETTA_MR_1", "MRBUMP", "cluster_1"))

      mrbDIR=[]
      for k in os.listdir(os.getcwd()):
         if os.path.isdir(k):
            mrbDIR.append(k)

      BESTCC=0.0
      BESTMODEL="None"
      BESTCYCLE=0
      BESTRESCOUNT=0
      BESTNCHAINS=0
      BESTAVERAGE=0.0
      BESTFILE="None"
      BESTPHRTIME=0.0
      BESTPHRLLG=0.0
      BESTPHRTFZ=0.0
      BESTREFMAC_I_RFREE=1.0
      BESTREFMAC_F_RFREE=1.0
      BESTREFMAC_I_RFACTOR=1.0
      BESTREFMAC_F_RFACTOR=1.0
      BEST_NORESMODEL=0
      BEST_NOCHAINSMODEL=0
      for l in mrbDIR:
         modelNAME=l.replace("search_","").replace("_mrbump","") 
         SHELDIR=os.path.join(l, "data", "loc0_ALL_"+modelNAME, "unmod", "mr", MRPROG, "build", "shelxe")         
         if os.path.isdir(SHELDIR):
            SHELLOG=os.path.join(SHELDIR, "shelxe-input.lst")
            if os.path.isfile(SHELLOG):
               if FOUND==False:
                  TOTALPROCESSED=TOTALPROCESSED+1
                  FOUND=True
               file=open(SHELLOG)
               line=file.readline()
               while line:
                  if "Best" in line:
                     BEST=float(string.split(line)[6].replace("%)",""))
                     CYCLE=int(string.split(line)[3])
                     if BEST >= BESTCC:
                        BESTCC=BEST
                        BESTMODEL=l
                        BESTCYCLE=CYCLE
                        BESTFILE=SHELLOG
                        if BEST >= 25.0:
                           TOTALGOOD+=1
                        BESTBUCCFILE=os.path.join(l, "data", "loc0_ALL_"+modelNAME, "unmod", "mr", \
                                     MRPROG, "build", "buccaneer", "buccaneer_run.log")
                        BESTREFMACFILE=os.path.join(l, "data", "loc0_ALL_"+modelNAME, "unmod", "mr", \
                                     MRPROG, "refine", "refmac_"+MRPROG+"_loc0_ALL_"+modelNAME+"_UNMOD.log")
                        BESTMRBUMPLOG=os.path.join(os.getcwd(), modelNAME + ".log")
                        if MRPROG=="phaser":
                           BESTPHASERLOG=os.path.join(l, "data", "loc0_ALL_"+modelNAME, "unmod", "mr", MRPROG, "phaser_loc0_ALL_"+modelNAME+"_UNMOD.log")
                  line=file.readline()
               file.close()
            
      #print j, ":", BESTMODEL, BESTCC, BESTCYCLE, BESTFILE
      if BESTFILE != "None":
         file=open(BESTFILE, "r")
         line=file.readline()
         while line:
            if ("Global" in line and int(string.split(line)[3])==BESTCYCLE) or BESTCYCLE == 1:
                  nline=file.readline()
                  CAPTURE=False
                  NCHAINS=0
                  while nline:
                     #print string.strip(nline)
                     if CAPTURE and ":" in nline :
                        nl=string.split(nline)
                        for i in nl:
                           if ":" in i:
                              NCHAINS+=1
                     if CAPTURE and ":" not in line:
                        CAPTURE=False
                        break
                     if "residues left after" in nline:
                        BESTRESCOUNT=string.split(nline)[0]
                        CAPTURE=True
                     nline=file.readline()
                  BESTNCHAINS=NCHAINS
                  BESTAVERAGE=float(BESTRESCOUNT)/float(BESTNCHAINS)
                  break
            line=file.readline()
         file.close()   
         PHRFILE=open(os.path.join(os.getcwd(), BESTPHASERLOG))
         phaserTime=0.0
         phaserLLG=0.0
         phaserTFZ=0.0
         for line in reversed(PHRFILE.readlines()):
            if "CPU Time" in line:
               phaserTime=float(string.split(line)[-2])
               break
         PHRFILE.close()
         BESTPHRTIME=phaserTime
         BESTPHRLLG,BESTPHRTFZ=getPhaserResults(BESTPHASERLOG)
         BESTREFMAC_I_RFREE,BESTREFMAC_F_RFREE,BESTREFMAC_I_RFACTOR,BESTREFMAC_F_RFACTOR,BEST_NORESMODEL,BEST_NOCHAINSMODEL=getRefmacResults(BESTREFMACFILE)
         BEST_NORESTARGET,BEST_NOCHAINSTARGET,RESOLUTION=getMrbumpResults(BESTMRBUMPLOG)
         
      if BESTCC >= 25.0:
         SOLUTIONCOUNT+=1
         SOL_LIST.append(j)
      if BESTCC >= 20.0 and BESTCC < 25.0 and BESTAVERAGE >= 10.0: 
         MARGINALCOUNT+=1
         MARG_LIST.append(j)
      if BESTAVERAGE >= 10.0:
         AVERAGECOUNT+=1
      if BESTCC >= 25.0 and BESTAVERAGE >= 10.0: 
         FULLSOLCOUNT+=1

      if BESTFILE != "None":
         #print "%.2lf %.2lf %d %d %s -> %s" % (BESTAVERAGE, BESTCC, BESTNCHAINS, BESTCYCLE, j, os.path.join(j, "cluster_1", "ROSETTA_MR_1", "MRBUMP", BESTFILE)) 
         #print os.path.join(os.getcwd(), BESTPHASERLOG)
         print "%.2lf %.2lf %d %d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %d %d %d %d %.2lf %s" \
            % (BESTCC, BESTAVERAGE, BESTNCHAINS, BESTCYCLE, BESTPHRTIME, BESTPHRLLG, BESTPHRTFZ, \
               BESTREFMAC_I_RFREE, BESTREFMAC_F_RFREE, BESTREFMAC_I_RFACTOR, BESTREFMAC_F_RFACTOR, \
               BEST_NORESMODEL, BEST_NOCHAINSMODEL, BEST_NORESTARGET, BEST_NOCHAINSTARGET, RESOLUTION, j) 
      #else:
      #   print j, "-> None" 

      os.chdir(quarkDIR)

print ""
print "No. of solutions with CC>=25 = %d/%d -> %.2lf success rate" % (SOLUTIONCOUNT,TOTALPROCESSED, float(SOLUTIONCOUNT)/float(TOTALPROCESSED))
print ""
print "List of successes:"
count=1
l2=[]
for i in SOL_LIST:
   if count < len(SOL_LIST):
      sys.stdout.write(i + ", ")
      count=count+1
   else:
      sys.stdout.write(i)
   l2.append(i)
sys.stdout.write("\n")
print ""
print "No. of solutions with average chain length>=10 = %d" % AVERAGECOUNT
print ""
print "No. of solutions with CC>=25 and average chain length>=10 = %d" % FULLSOLCOUNT
print ""
print "No. of solutions with CC>=25 and average chain length>=10 = %d/%d -> %.2lf success rate" % (FULLSOLCOUNT,TOTALPROCESSED, float(FULLSOLCOUNT)/float(TOTALPROCESSED))
print ""
print "No. of marginals with 20<=CC<25 and average chain length>=10 = %d" % MARGINALCOUNT
print ""
print "List of marginals:"
count=1
marg=[]
for i in MARG_LIST:
   if count < len(MARG_LIST):
      sys.stdout.write(i + ", ")
      count=count+1
   else:
      sys.stdout.write(i)
   marg.append(i)
print ""
print ""

#l1=["3DF8", "2G7O", "1Y6X", "3BJO", "1WPA", "1GVD", "3MXZ", "2F60", "1YZM", "2OUF", "3A4C", "3G2B", "3CEC", "2GKR", "1RIY", "2OXO", "3GOE", "1EZJ", "2FQ3", "3K3V", "3HRO", "3CI9", "3LAX", "2QMT", "2ZQE", "2QYW", "1UCS", "2Q2F", "2YZT", "1USM", "3JTZ", "3F2E", "2FU2", "3E21", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "3I8Z", "2O4T", "2P6V", "3NRW", "2VKL", "1TUK", "1T07", "1VBW", "1ZZK", "2QVO", "2O1K", "1OX3", "1ZVA", "3CE7", "1G6U", "2OVG", "1U84", "1GK6", "2OQQ", "2PST", "2IP6", "1USE", "2FI0", "1OKS", "1Z0P", "1VYI", "3CTR", "1YIB", "3EFG", "2O9U", "3JSR", "3KW6", "2ZQM", "2JKU", "3H36", "1R7J", "3NZL", "2RHF", "1Q8D", "3L32", "1Y0N", "3AGN", "2QSB", "1TTZ", "3HZ7"]

l1=["3DF8", "2G7O", "1Y6X", "3BJO", "1WPA", "1GVD", "3MXZ", "2F60", "1YZM", "2OUF", "3A4C", "3G2B", "3CEC", "2GKR", "1RIY", "2OXO", "3GOE", "1EZJ", "2FQ3", "3K3V", "3HRO", "3CI9", "3LAX", "2QMT", "2ZQE", "2QYW", "1UCS", "2Q2F", "2YZT", "1USM", "3JTZ", "3F2E", "2FU2", "3E21", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "3I8Z", "2O4T", "2P6V", "3NRW", "2VKL", "1TUK", "1T07", "1VBW", "1ZZK", "2QVO", "2O1K", "1OX3", "1ZVA", "3CE7", "1G6U", "2OVG", "1U84", "1GK6", "2OQQ", "2PST", "2IP6", "1USE", "2FI0", "1OKS", "1Z0P", "1VYI", "3CTR", "1YIB", "3EFG", "2O9U", "3JSR", "3KW6", "2ZQM", "2JKU", "3H36", "1R7J"]

l100=["3DF8", "3IDW", "2G7O", "1Y6X", "3BJO", "1WPA", "1GVD", "1Y0N", "3MXZ", "1GXU", "2F60", "1YZM", "1Q8D", "3FBL", "2OUF", "3A4C", "3JVL", "2RHF", "3G2B", "3CEC", "2GKR", "1RIY", "1TGR", "3BRI", "2OXO", "3FT7", "3GOE", "1EZJ", "2HDZ", "2FQ3", "3K3V", "1TUW", "3HRO", "3AGN", "3CI9", "2VC8", "3LDC", "3OSH", "3LAX", "1WHZ", "2QMT", "2ZQE", "2QYW", "2QSB", "2Q2F", "2YZT", "3JTZ", "3F2E", "2FU2", "2YVI", "3E21", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "3I8Z", "2O4T", "2HPJ", "3GHF", "2P6V", "3NRW", "2VKL", "1TUK", "1T07", "1VBW", "1ZZK", "2QVO", "2O1K", "1OX3", "1ZVA", "3CE7", "1G6U", "2OVG", "2EWT", "2NS0", "1OAP", "1U84", "1GK6", "2OQQ", "2PST", "2IP6", "3DML", "1USE", "1OKS", "1Z0P", "3BN7", "2IGP", "1VYI", "3CTR", "3NZL", "1YIB", "3EFG", "3G21", "2B8I", "3N3F", "2O9U", "3JSR", "3KW6", "2ZQM", "2JKU", "3H36", "1R7J"]

l90=["2BKF", "3DF8", "2G7O", "1Y6X", "1WPA", "1GVD", "3MXZ", "1GXU", "1FK5", "2F60", "1YZM", "1Q8D", "2OUF", "3A4C", "3JVL", "2RHF", "3G2B", "3CEC", "2GKR", "1RIY", "1TGR", "1Z96", "2OXO", "3FT7", "1EZJ", "3KKF", "2HDZ", "2FQ3", "2ES9", "3HRO", "3AGN", "3CI9", "2H8E", "3OSH", "3B64", "1WHZ", "1IQZ", "2ZQE", "2QSB", "2Q2F", "2YZT", "3JTZ", "3F2E", "2FU2", "2C60", "1VJK", "3CQ1", "2O37", "2YVI", "3E21", "2ZXY", "3L32", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "2O4T", "2HPJ", "2P6V", "2VKL", "1T07", "1ZZK", "2QVO", "2O1K", "1OX3", "1ZVA", "3CE7", "2WUJ", "1G6U", "3F14", "2OVG", "2EWT", "2RFF", "1OAP", "1U84", "1GK6", "2OQQ", "2IP6", "3FF5", "2EWH", "1USE", "2FI0", "3NBM", "2NUH", "2IGP", "1VYI", "3MSH", "3NZL", "1YIB", "3EFG", "2O9U", "3JSR", "2ZQM"]

molrepl=["1USE", "2FQ3", "1GVD", "1OKS", "3NZL", "1Y6X", "3JSR", "3HRO", "1YZM", "3GOE", "2OQQ", "1T07", "2RHF", "2FU2", "1Q8D", "2IP6", "1VBW", "3L32", "1R7J", "3BJO", "3A4C", "2QMT", "1Y0N", "3EFG", "2P5K", "3AGN", "1EZJ", "1ZVA", "1J8B", "1TUK", "2QSB", "1TTZ", "3DF8", "1OX3", "2P6V", "3MXZ", "1UCS", "3G2B", "3NRW", "3E21", "1I2T", "2OUF", "1RIY", "2O1K", "2F60", "2ZQE", "2VKL", "2GKR", "1EW4", "3CE7", "1ZZK", "2ZQM", "1GK6", "2QYW", "2QVO", "1YIB", "2YZT", "2JKU", "3HZ8", "2G7O"]

rosettaList=["1EJG", "1EW4", "1EZJ", "1FK5", "1G2R", "1G6U", "1GK6", "1GVD", "1GXU", "1I2T", "1J8B", "1OAP", "1OKS", "1OX3", "1Q8D", "1R6J", "1R7J", "1RIY", "1RW1", "1T07", "1TGR", "1U84", "1UJ8", "1USE", "1USM", "1V2Z", "1V70", "1VBW", "1VJK", "1WHZ", "1WPA", "1Y0N", "1Y6X", "1YIB", "1YU5", "1YZM", "1Z0P", "1Z96", "1ZVA", "1ZZK", "2B8I", "2C60", "2CWY", "2D3D", "2EFV", "2ES9", "2F60", "2FI0", "2FQ3", "2FU2", "2G7O", "2GKR", "2GPI", "2H9U", "2HDZ", "2HL7", "2HPJ", "2I4A", "2IGP", "2IP6", "2JKU", "2NML", "2NS0", "2NUH", "2O1K", "2O37", "2O4T", "2OQQ", "2OUF", "2OVG", "2OXO", "2P5K", "2P6V", "2PST", "2Q2F", "2QFF", "2QMT", "2QSB", "2QVO", "2QYW", "2RFF", "2RHF", "2V75", "2VC8", "2VKL", "2YZT", "2ZQE", "2ZQM", "3A4C", "3ADG", "3B64", "3BJO", "3BN0", "3BRI", "3C0F", "3CE7", "3CEC", "3CQ1", "3DF8", "3E21", "3EFG", "3F2E", "3FBL", "3FF5", "3FKC", "3FMY", "3FT7", "3G21", "3G2B", "3GOE", "3H01", "3H36", "3H8H", "3HGL", "3HRO", "3HZ7", "3IDW", "3IM3", "3JTZ", "3JVL", "3K3V", "3KW6", "3LAX", "3LBJ", "3MWZ", "3MXZ", "3NRW", "3OOU"]

oldQuark=["3EFG", "2P6V", "3CI9", "2FU2", "2RFF", "2HL7", "3F14", "3NRW", "3KW6", "2ZQM", "2O9U", "3MWZ", "2GKR", "1OKS", "2CMP", "2YZT", "2O4T", "3H01", "1T07", "2FI0", "2OVG", "2HDZ", "2QYW", "3OOU", "1J8B", "2GPI", "3BN7", "3JTZ", "1EW4", "1Y0N", "3DF8", "1VYI", "2PST", "3NZL", "3G2B", "1OX3", "2V75", "2P5K", "1Q8D", "1Y6X", "1VBW", "1USM", "2ZQE", "2FQ3", "3E21", "2QFF", "2OQQ", "3MXZ", "3HZ7", "2OUF", "2OXO", "1EZJ", "1GVD", "1G6U", "3CE7", "1GK6", "2O1K", "2IP6", "1OAP", "1I2T", "3BJO", "1R7J", "2VKL", "2QVO", "3IDW", "2RHF", "2F60", "2QSB", "1USE", "1YZM", "2G7O", "2Q2F", "1YIB", "3HRO", "1RIY", "1ZZK", "1ZVA"]

quark0_1l=["3DF8", "3IDW", "2G7O", "1Y6X", "3BJO", "1WPA", "1GVD", "1Y0N", "3MXZ", "1GXU", "2F60", "1YZM", "1Q8D", "3FBL", "2OUF", "3A4C", "3JVL", "2RHF", "3G2B", "3CEC", "2GKR", "1RIY", "1TGR", "3BRI", "2OXO", "3FT7", "3GOE", "1EZJ", "2HDZ", "2FQ3", "3K3V", "1TUW", "3HRO", "3AGN", "3CI9", "2VC8", "3LDC", "3OSH", "3LAX", "1WHZ", "2QMT", "2ZQE", "2QYW", "2QSB", "2Q2F", "2YZT", "3JTZ", "3F2E", "2FU2", "2YVI", "3E21", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "3I8Z", "2O4T", "2HPJ", "3GHF", "2P6V", "3NRW", "2VKL", "1TUK", "1T07", "1VBW", "1ZZK", "2QVO", "2O1K", "1OX3", "1ZVA", "3CE7", "1G6U", "2OVG", "2EWT", "2NS0", "1OAP", "1U84", "1GK6", "2OQQ", "2PST", "2IP6", "3DML", "1USE", "1OKS", "1Z0P", "3BN7", "2IGP", "1VYI", "3CTR", "3NZL", "1YIB", "3EFG", "3G21", "2B8I", "3N3F", "2O9U", "3JSR", "3KW6", "2ZQM", "2JKU", "3H36", "1R7J"]

quarkPH5=["2BKF", "3DF8", "3IDW", "2G7O", "1Y6X", "3BJO", "1WPA", "1GVD", "1Y0N", "3MXZ", "1GXU", "1FK5", "2F60", "1YZM", "2QFF", "1Q8D", "2OUF", "3A4C", "3JVL", "2RHF", "3G2B", "3CEC", "2GKR", "1RIY", "1TGR", "3BRI", "1Z96", "2OXO", "3FT7", "3GOE", "1EZJ", "2HDZ", "2FQ3", "2ES9", "3HRO", "3AGN", "3CI9", "2H8E", "3LDC", "3OSH", "3LAX", "1WHZ", "2QMT", "2ZQE", "2EFV", "2CWR", "2QYW", "2QSB", "2Q2F", "2YZT", "3JTZ", "3F2E", "2FU2", "2O37", "2YVI", "1YU5", "3E21", "2ZXY", "3L32", "1I2T", "1J8B", "2P5K", "1EW4", "2HL7", "3I8Z", "2O4T", "2HPJ", "3GHF", "2P6V", "3HGL", "3NRW", "2VKL", "1TUK", "1T07", "1VBW", "1ZZK", "2V75", "2QVO", "2O1K", "3HRL", "1OX3", "1ZVA", "3CE7", "2WUJ", "2OVG", "2EWT", "2NS0", "1OAP", "1U84", "1GK6", "2OQQ", "3FMY", "2PST", "2IP6", "1USE", "2FI0", "1OKS", "3IM3", "1Z0P", "2NUH", "2IGP", "1VYI", "3CTR", "3NZL", "1YIB", "3EFG", "3G21", "2B8I", "2O9U", "3JSR", "3KW6", "2ZQM", "2JKU", "3H36", "3OIZ", "1R7J"]

#l2=["3DF8", "2G7O", "3BJO", "1Y0N", "2F60", "1YZM", "2OUF", "2GKR", "1RIY", "2HDZ", "3HRO", "2ZQE", "2Q2F", "2YZT", "2FU2", "3E21", "1I2T", "1J8B", "1EW4", "2O4T", "2P6V", "2VKL", "2QVO", "2O1K", "1ZVA", "1GK6", "2OQQ", "1USE", "1YIB", "3EFG", "3JSR", "2ZQM", "1R7J"]

l1_list=[]
l1_notsolved=[]
l2_list=[]
l2_unilist=[]
l2_notRosetta=[]
l2_notoldQuark=[]
l2_molrepOnly=[]
l2_NEW=[]
l100_not=[]
l90_not=[]
l90_uni=[]

lq2_not=[]
lq2_new=[]
PH5_not=[]
PH5_uni=[]

for i in l2:
   if i not in quarkPH5:
     PH5_not.append(i)
print "Solved in this run but not in QUARK-Phaser-5min-KILL:", len(PH5_not)
print PH5_not
print ""

for i in quarkPH5:
   if i not in l2:
     PH5_uni.append(i)
print "Solved in QUARK-Phaser-5min-KILL but not yet in this run:", len(PH5_uni)
print PH5_uni
print ""

for i in l2:
   if i not in l90:
     l90_not.append(i)
print "Solved in this run but not in QUARK-Phaser-0.1-OLDKILL:", len(l90_not)
print l90_not
print ""

for i in l90:
   if i not in l2:
     l90_uni.append(i)
print "Solved in QUARK-Phaser-0.1-OLDKILL but not yet in this run:", len(l90_uni)
print l90_uni
print ""

for i in l2:
   if i not in quark0_1l:
      lq2_new.append(i)
print "Solved in this run but not in last Quark run (Phaser RMSD=1.2 + no kill)):", len(lq2_new)
print lq2_new
print ""

for i in quark0_1l:
   if i not in l2:
      lq2_not.append(i)
print "Solved in last Quark run (Phaser RMSD=1.2 + no kill)) but not in this run:", len(lq2_not)
print lq2_not
print ""

for i in l2:
   if i not in l1:
      l2_unilist.append(i)

print "Unique to this run (compared to original Quark run (Phaser RMSD=1.2)):", len(l2_unilist)
print l2_unilist
print ""

for i in l1:
   if i not in l2:
      l1_notsolved.append(i)

print "Solved in original Quark run (Phaser RMSD=1.2) but not yet solved in this run:", len(l1_notsolved)
print l1_notsolved
print ""

for i in molrepl:
   if i not in l2:
      l2_molrepOnly.append(i)

print "Molrep only solutions:", len(l2_molrepOnly)
print l2_molrepOnly
print ""

for i in l2:
   if i not in rosettaList:
      l2_notRosetta.append(i)
for i in molrepl:
   if i not in rosettaList:
      l2_notRosetta.append(i)

print "Quark solutions (phaser + molrep) not found using Rosetta:", len(l2_notRosetta)
print l2_notRosetta
print ""

for i in l2:
   if i not in oldQuark:
      l2_notoldQuark.append(i)
for i in molrepl:
   if i not in oldQuark and i not in l2:
      l2_notoldQuark.append(i)

print "Quark solutions (phaser + molrep) not found in course-grain truncation Quark run:", len(l2_notoldQuark)
print l2_notoldQuark
print ""

for i in l2:
   if i not in oldQuark and i not in rosettaList and i not in l1 and i not in molrepl:
      l2_NEW.append(i)

print "Completely new solutions:", len(l2_NEW)
print l2_NEW
print ""

#for i in `ls -d */`
#do 
#  if [ -d "${i}cluster_1/ROSETTA_MR_1/MRBUMP" ]; then
#     cd ${i}cluster_1/ROSETTA_MR_1/MRBUMP 
#     if ! [ -z `find . | grep shelxe_ | grep log | head -n 1` ]; then
#          echo $i 
#          grep Best `find . | grep shelxe_ | grep log` | gawk '{print $8 " " $0}'  | sort -n | tail -n 1 
#     fi
#     cd ../../../../ 
#   fi
#done

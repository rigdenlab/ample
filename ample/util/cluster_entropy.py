#!/usr/bin/python

#
#  divide into parallel, select one, compare to rest- measure cluster repeat untill largest
#sort and rewrite
import re
import os, glob
from multiprocessing import Process, JoinableQueue


##############
def count_files_in_dir(path):
 counter =0
 for infile in glob.glob( os.path.join(path, '*.pdb') ):
  counter=counter+1
 return counter

##################
def split(models_path, outpath, NProc, NMODELS):  #split into parallel comparisons for speed

 
 cluster_number =0
 counter=0
 divide= NMODELS / NProc #number of models in each split folder  use to divide 
 if divide <1:
  divide = 1


 divided_files=[]
 files = []

 for infile in glob.glob( os.path.join(models_path, '*.pdb') ):
  counter = counter +1
  files.append(infile)
 no_of_models =  counter




 files.sort()
 for i in xrange(0, len(files), divide):
    divided_files.append(outpath + 'divide/' + str(i) )
    #print files[i:i+divide]
    if not os.path.exists(outpath + 'divide/'):
     os.system ('mkdir ' + outpath + 'divide/')
    if not os.path.exists(outpath + 'divide/'+ str(i)):
     os.system ('mkdir ' + outpath + 'divide/' + str(i) )
    for j in files[i:i+divide]:
    
     #print 'cp ' + models_path + '/' + str(j) + ' ' + outpath + 'divide/' + str(i)
     os.system ('cp ' +str(j) + ' ' + outpath + 'divide/' + str(i) )
# print divided_files

 return divided_files


###############
def test_2_files(test1, test2, path_to_LGA, gdtout):

 string = test2
 SCORE = ''
 path = os.getcwd()
 cur_path = path

 if not os.path.exists( path + '/run_dir_'+gdtout):
  os.system ('mkdir ' + path + '/run_dir_'+gdtout)
  os.system ('cp ' +path_to_LGA + '/lga ' +path + '/run_dir_'+gdtout )
  os.system ('mkdir ' + path + '/run_dir_'+gdtout+'/MOL2') 
  os.system ('mkdir ' + path + '/run_dir_'+gdtout+'/TMP')
 
 os.chdir(path + '/run_dir_'+gdtout)
 my_LGA_run = open (path + '/run_dir_'+gdtout+'/LGA_run'+gdtout, "w")
 my_LGA_run.write('ulimit -s unlimited\n'+ path + '/run_dir_'+gdtout+'/lga run.pdb -3 -sda -atom:CA \n')
 my_LGA_run.close()
 os.system ('chmod uoga=wrx ' + path + '/run_dir_'+gdtout+'/LGA_run'+gdtout)


 test1 = open (test1)
 test2 = open (test2)
 test1_out = open(path + '/run_dir_'+gdtout+'/MOL2/run.pdb', "w")
 test1_out.write('MOLECULE 1\n')
 for line in test1:
  test1_out.write(line)
 test1_out.write('ENDMOL\n')
 test1_out.write('MOLECULE 2\n')
 for line in test2:
  test1_out.write(line)
 test1_out.write('END\n')
 test1_out.close()
#test2_out.close()



 os.system(path + '/run_dir_'+gdtout+'/LGA_run'+gdtout+ ' > RESULTS')####

 RESULTS = open(path + '/run_dir_'+gdtout+'/RESULTS')

 for line in RESULTS:
  get_stat1 = re.compile('^SUMMARY')
  result_stat1 = get_stat1.match(line)
  if result_stat1:
   # print line
    stat1_get = re.split('SUMMARY\(GDT\)\s*(\d*)\s*(\d*)\s*(\d*.\d*)\s*(\d*)\s*(\d*.\d*)\s*(\d*.\d*)', line)  
    SCORE = stat1_get[6] 
 #
 os.system('rm -r ' +path + '/run_dir_'+gdtout)
 os.chdir(cur_path)

 return SCORE
################
# compare 1 folder - thread this
#takes the split folder, top model, name, place to output 'GET',threshold for comparison (20, 25, ... 
def compare_one_split_folder(folder, random_model_to_compare, suffix, current_cluster_output, threshold, q, path_to_LGA):
     
     #Number_passed =0
     for infile in glob.glob( os.path.join(folder, '*.pdb') ): #loop through models in folder
         name= re.split('/', infile)
         name= name.pop()
         name= re.split('.pdb', name)
         name =  name[0]

         SCORE = test_2_files(random_model_to_compare, infile, path_to_LGA, suffix + '_' +name )
      #   print 'GOT ' + SCORE 
         if float(SCORE) >= float(threshold):
       #      print SCORE
             #print 'cp ' + infile + ' ' + outpath + '/' 
             os.system('mv ' + infile + ' ' + current_cluster_output + '/' )

             #Number_passed +=1
             #q.put(Number_passed)

      
########################
def SORT_clusters(LIST_of_clusters, outpath):  #sort clusters and renumber by size

  # [cluster, size]
  sorted_by_size = []
  sorted_by_size = sorted(LIST_of_clusters, key=lambda student: student[1], reverse=True)
  print sorted_by_size

  largest_size = sorted_by_size[0][1]

  sorted_cluster = 0
  for cluster in  sorted_by_size:
    os.system('mv ' + outpath +'/cluster'+str(cluster[0]) +' ' + outpath + '/sorted_cluster_' +str(sorted_cluster) )
    sorted_cluster += 1


  return largest_size, outpath + '/sorted_cluster_0'

#########################
def cluster(path_of_models, LGA, nproc, NMODELS, outpath, list_of_cores):
  #list_of_cores=[20]#, 25]#, 30, 40, 50, 60, 70, 80]
    
  ###temp dir
  if os.path.exists(outpath + '/temp'):
    os.system('rm -r ' +outpath + '/temp')
  os.system('mkdir '+ outpath + '/temp')
  for infile in glob.glob( os.path.join(path_of_models, '*.pdb') ):
    os.system('cp '+infile +' '+outpath + '/temp/')
  path_of_models = outpath + '/temp'

  random_model_to_compare = ''
  LIST_of_clusters = []
  LIST_of_cluster_sizes=[]
  used_random_models = []

  for core in list_of_cores:
    if not os.path.exists(outpath + 'cluster_' + str(core) ):
     os.system('mkdir ' +outpath + 'cluster_' + str(core)  )
  
    os.chdir(outpath + 'cluster_' + str(core) ) 
    if os.path.exists(outpath + 'split/'):
      os.system('rm -r ' +  outpath + 'split/' ) 

    if not os.path.exists(outpath + 'split/'):
      os.mkdir( outpath + 'split/')
    list_of_split_paths = split(path_of_models, outpath + 'split/', nproc, NMODELS)
    #print list_of_split_paths


    GOT_LARGEST_CLUSTER = False
    cluster = 0
    temp_path_name = path_of_models

    files = []
    for infile in glob.glob( os.path.join(path_of_models, '*.pdb') ):
       files.append(infile)

    while  GOT_LARGEST_CLUSTER == False:
      current_cluster = outpath + 'cluster_' + str(core) +'/cluster'+str(cluster)
      os.system('mkdir ' +outpath + 'cluster_' + str(core) +'/cluster'+str(cluster))




      #get a random model as a starting 'top model'

      for each_file in files:
         if not each_file in used_random_models:
          random_model_to_compare = each_file
          used_random_models.append(random_model_to_compare)
          #print used_random_models, random_model_to_compare
          break

      os.chdir(outpath)
    #run comparisons parallel against 'top model' (thread to make parallel) if similar- keep

      suffix_char_used = []
     # for letter in range(ord('a'), ord('z') + 1):
    #  print chr(letter)
     #   suffix_char.append(chr(letter))
        

      q = JoinableQueue()
      number_passed = 0

      threads = []
      counter = 0
      for split_path in list_of_split_paths:  ## Threadding to make parallel, the number of folders = number of comparisons
                                 ## as each takes little processor power, can make more than the number of cores?
                                            ##can split into 100 to run fast but with lag

        counter +=1
        suffix = str(counter)
        #compare_one_split_folder
        thread = Process(target=compare_one_split_folder, args=(split_path, random_model_to_compare, suffix, current_cluster, core, q, LGA))        
        thread.start() 
        threads.append(thread)
        #print 'number_passed ', q.get()


      for thread in threads:

        thread.join()


      

      number_in_cluster = 0
      for infile in glob.glob( os.path.join(outpath + 'cluster_' + str(core) +'/cluster'+str(cluster), '*.pdb') ): 
       number_in_cluster += 1
       name= re.split('/', infile)
       name= name.pop()
       #if  temp_path_name+'/'+name in files:
       files.remove(temp_path_name+'/'+name)

      LIST_of_clusters.append([cluster, number_in_cluster])
      LIST_of_cluster_sizes.append(int(number_in_cluster))

      print 'current cluster size= ',  number_in_cluster
     

      if float(number_in_cluster)  >= float(NMODELS/2):
   
        GOT_LARGEST_CLUSTER = True
        print 'got_largest_cluster', number_in_cluster , NMODELS, ' cluster_number ',  cluster  

        #return the largest_cluster
        largest_size, largest_cluster_path  = SORT_clusters(LIST_of_clusters, outpath + '/cluster_' + str(core) )
        return largest_size, largest_cluster_path

      else:
       cluster += 1
       #check if all done
       #if a cluster is larger than reamining files
       LIST_of_cluster_sizes.sort()
       LIST_of_cluster_sizes.reverse()

       total_files_in_all_clusters = 0
       for sizes in LIST_of_cluster_sizes:
         total_files_in_all_clusters +=sizes #total no of sorted files
       
       remaining_unsorted = NMODELS - total_files_in_all_clusters  #get number of remaining files 

       print 'CLUSTER= ', cluster, ' sorted= ', total_files_in_all_clusters, ' remaining= ', len(files), ' largest so far = ', LIST_of_cluster_sizes[0]
      # print random_model_to_compare
       if float(LIST_of_cluster_sizes[0]) >= float(len(files)):  #if there is a cluster bigger than the remaining files

         largest_size, largest_cluster_path = SORT_clusters(LIST_of_clusters, outpath + '/cluster_' + str(core) )
         #print 'RETURNING largest ', largest_size,  ' cluster ' ,    largest_cluster
         return largest_size, largest_cluster_path








############
def RUN_CLUSTER(levels, path, path_to_LGA, outpath, nproc, NMODELS):
 LOG = open(outpath +'/CLUSTERING_LOG',"w")
 largest_cluster = path
 cluster_results=[]

 index = 0
 for level in levels:

  
  using_level = [levels[index]]
  NMODELS, largest_cluster = cluster(largest_cluster, path_to_LGA, nproc, NMODELS, outpath, using_level )
  print 'DONE', largest_cluster
  print 'NEW NMODELS ' , NMODELS
  LOG.write(str(level) + ' ' + str(NMODELS) + '\n')
  if NMODELS > nproc:
    nproc = NMODELS-1
    
  cluster_results.append(NMODELS)
  index +=1
  
  if NMODELS <=2:  # stop at threshold or crash
    return cluster_results

 return cluster_results
  



####################################
#path = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/models' 
#path_to_LGA = '/home/jaclyn/Desktop/cmd/casp/LGA_package'
#outpath = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/run/'
#nproc = 1
#NMODELS = 1000#

#levels = [20,25,30,40,50,60]


#path = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/models' 


#outpath = '/home/jaclyn/Desktop/backup_scripts/NEW_WORKFLOW/MASTER_parallel/PROGRAM/Version_0.1/test/clusters/' 
#nproc = 8 
#NMODELS = 10

#cluster_results = RUN_CLUSTER(levels, path, path_to_LGA, outpath, nproc, NMODELS)

#print 'GOT', cluster_results


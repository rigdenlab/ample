/*
 * Contains code taken from:
* /Users/jmht/Downloads/openmpi-1.6.4/ompi/tools/ompi-profiler/ompi-profiler.c
 * https://svn.mcs.anl.gov/repos/darshan/branches/trac5/test/io-sample.c

 */
#include <errno.h>
#include <fcntl.h> /* for O_CREAT  etc. */
#include <libgen.h> /* for dirname, basename */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/utsname.h> /* for uname */
#include <unistd.h>
#include <sys/wait.h>

//#define MPI_CHECK(r,f) { if(r != MPI_SUCCESS) { printf("%s: (%d) failed: %d\n", f, __LINE__, r); } }

/* Global Variables */
int CHDIR=1;
int rank;
char my_lockfile[1024];
char global_lockfile[1024];


int get_lock() {
	/*
	 * Acquire the lock
	 */

	/* Create our own lockfile if it doesn't already exist */

	int jobid=100;
	int rank=1;

	/* Create hard link to the global LOCKFILE */
	int count=0;
	while ( link_lockfile( my_lockfile ) != 0 ) {
		// Sleep by rank to prevent 2 nodes getting in competition
		sleep( (rank/100)+1 );
		count++;
		if ( count > 1000 ) {
			printf("MORE THEN 1000 STEPS IN LOCK\n");
			return 1;
		}
	}

	//printf("GOT LOCK WITH COUNT %d\n",count);

	return 0;
}

int link_lockfile( char *my_lockfile ) {

	if ( link( my_lockfile, global_lockfile) == 0 ) {
		//printf("Lock was successful on link\n");
		return 0;
	}

	/* May or may not succeed, what matters is the hardlink count on our file */
    struct stat buf;
    if ( stat( my_lockfile, &buf) != 0 ) {
		char estr[256];
		sprintf( estr, "Error stating lockfile: %s\n", my_lockfile );
		perror(estr);
		exit(-1);
    }

    //printf("Size of file is: %d\n",(int)buf.st_size);
    if (buf.st_nlink == 2) {
    	//printf("SUCCESSFUL LINK COUNT\n");
    	return 0;
    } else {
    	//printf("GOT WRONG LINK COUNT\n");
    	return 1;
    }
}


void release_lock(){
	/*
	 * delete the lock files
	 */
	if ( unlink( global_lockfile ) != 0 ) {
		printf("ERROR REMOVING GLOBAL LOCKFILE\n");
	}
	return;
}


char* get_hostname(){

	/* Get the hostname */
	//int MPI_Get_processor_name( char *name, int *resultlen )
	struct utsname name;
	if ( uname( &name ) != 0 ) {
		return NULL;
	}
	return strdup(name.nodename);
}

int get_jobid() {
	/* Return the id of the job or -1 if error */
	char *tmp = getenv("LSB_JOBID");
	if ( tmp == NULL ) {
		return -1;
	}
	return atoi(tmp);
}

int get_ppn( char *hostname ){
	/* Return the number of processors per node */

	int colour=notdjb2_hash( hostname );
	//printf("Proc %d got hash %d\n", rank, colour);

	int key=0; /* follow ordering in MPI_COMM_WORLD */
	MPI_Comm hostcomm;
	if ( MPI_Comm_split(MPI_COMM_WORLD, colour, key, &hostcomm) != 0 ) {
		return -1;
	}

	int hsize,hrank;
	MPI_Comm_rank(hostcomm, &hrank);
	MPI_Comm_size(hostcomm, &hsize);
	printf("Proc %d got hrank %d size %d in hostcomm\n", rank, hrank, hsize);
	return hsize;
}


int next_job( char *filename, char *script )
{
	/*
	 * Extract the next job from the file pointed to by filename
	 * and remove that job from the file. Put the value we read from
	 * the file in script.
	 * Return 0 if we read a script, 1 if EOF or -1 for failure
	 *
	 * Initially tried allocating memory for the script in here and just returning
	 * the pointer, but got segfaults when running on > 1 processor so now we just
	 * allocate the memory in main and return an int
	 */

    char estr[1000];
    int ok = get_lock();

    /* Now assume we are safe to do our read/write */
    FILE *file = fopen( filename, "r" );
    if ( file == NULL ) {
    	sprintf( estr, "Error opening file: %s\n",filename);
    	perror(estr);
    	return -1;
    }

    // Determine the length
    struct stat buf;
    if ( fstat( fileno(file), &buf) != 0 ) {
    	sprintf( estr, "Error stating file: %s\n",filename);
    	perror(estr);
    	return -1;
    }

    //printf("Size of file is: %d\n",(int)buf.st_size);
    if (buf.st_size == 0) {
    	printf("Got empty file!\n");
    	release_lock();
    	memset(script, 0, sizeof(script)); /* Fill string with zeros */
    	return 1;
    }

    int len;
    char line[1024];
    memset(line, 0, sizeof(line)); /* Fill string with zeros */
    while (NULL != fgets(line, sizeof(line), file)) {
        len = strlen(line);

        /* printf("line len before: %d\n", len); */
        /* remove any trailing newline */
        if (line[len-1] == '\n') {
            line[len-1] = '\0';
        }

    } // End while loop
    fclose(file);

    // Copy the line into the string we will return
    //char *script = strndup( line, strlen(line) );
    strncpy( script, line, strlen(line) );

    /*
    printf("line len %d is %s\n", strlen(line), line );
    printf("script len %d is %s\n", strlen(script), script);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("proc %d got script %s\n",rank,script);
    */

    // Truncate the file by the length of the string we just read
    if ( truncate(filename, buf.st_size-len ) != 0 ) {
    	perror("Error truncating file!\n");
    }

    /* Can now release the lock */
    release_lock();

    return 0;
}

int run_job( char *script ) {
	/* Spawn a new process to run the given script */

	char estr[1000]; /* for error messages */

	/* Check the script can be found */
    if ( access( script, R_OK ) != 0 ) {
    	sprintf( estr,"Error accessing script: %s\n",script);
    	perror(estr);
    	return EXIT_FAILURE;
    }

    /*Spawn a child to run the program.*/
    pid_t pid=fork();
    if (pid==0) { /* child process */

        char * logfile;
        char * scriptpath;
    	if ( CHDIR ) {

    		/* find directory name  and script name */
    		char *script_copy1 = strdup(script);
    		char * jobdir = dirname(script_copy1);
    		/* printf("got dirname: %s\n",jobdir); */

    		char *script_copy2 = strdup(script);
    		scriptpath = basename(script_copy2);
    		/* printf("got scriptname: %s\n",scriptname); */


    		/* create log name */
    		char *suffix=".log";
    		logfile = (char *)malloc( (strlen(scriptpath) + strlen(suffix) + 1) * sizeof(char) );
    		logfile = strncpy( logfile, scriptpath, strlen(scriptpath)+1 );
    		logfile = strncat( logfile, suffix, strlen(suffix) );
    		/* printf("got logfile: %s\n",logfile); */

    		if ( chdir( jobdir ) != 0) {
    			sprintf(estr, "Error changing to directory: %s\n",jobdir);
                        perror(estr);
                        return EXIT_FAILURE ;
    		}

    	} else {
			/* create log name */
                        scriptpath = strdup( script );
			char *suffix=".log";
			logfile = (char *)malloc( (strlen(scriptpath) + strlen(suffix) + 1) * sizeof(char) );
			logfile = strncpy( logfile, scriptpath, strlen(scriptpath)+1 );
			logfile = strncat( logfile, suffix, strlen(suffix) );
			//printf("Child logfile name is: %s\n",logfile);
    	}

            /* Capture stdout and stderr into a logfile */
            int fd = open(logfile, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);

    	    dup2(fd, 1);   // make stdout go to file
    	    dup2(fd, 2);   // make stderr go to file
    	    close(fd);     // fd no longer needed - the dup'ed handles are sufficient

            char *argv[]={"sh",scriptpath, NULL};
            //char *argv[]={"sh",script, NULL};

            execv("/bin/sh",argv);

            /* only if execv fails */
            sprintf( estr, "Error executing script %s\n", script);
            perror(estr);
            return EXIT_FAILURE ;
    }
    else { /* pid!=0; parent process */
    		int status;
            if ( waitpid(pid,&status,0) == -1 ) {
            	perror("wait()");
            	return EXIT_FAILURE;
            }
            if( WIFEXITED(status) ){
               /* printf("%ld exited with return code %d\n",
                      (long)pid, WEXITSTATUS(status)); */
               if ( WEXITSTATUS(status) != 0 ) {
            	   // non-zero exit status
            	   printf("Script %s exited with code %d\n", script, WEXITSTATUS(status));
            	   return WEXITSTATUS(status);
               }
             } else {
            	 sprintf( estr, "Error executing script3: %s\n", script);
             	perror(estr);
             	return EXIT_FAILURE;
             }
    }

    return EXIT_SUCCESS;

} //End run_job


int notdjb2_hash(char *str)
{
	/*
	 * Return a hash of the given string
	 * mangled version of djb2: http://www.cse.yorku.ca/~oz/hash.html
	 *
	 */
    int hash = 5381;
    int c;
    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    if ( hash < 0 ) hash *= -1;

    return hash;
}

void exit_error( char *msg ) {
	/* Shut down everything on error */
	if ( msg != NULL ) {
		printf( msg );
	}
	printf( "%s\n" ,strerror( errno ) );
	MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	exit(EXIT_FAILURE);
}

int main(argc, argv)
	int argc;
	char *argv[];
{
	char estr[256];
	int rank, size;

	/* Get file with list of jobs */
	if ( argc != 2 ) {
		printf( "Usage is %s <jobfile>\n", argv[0] );
		exit( EXIT_FAILURE );
	}
	char *jobfile = argv[1];
	if ( access( jobfile, R_OK ) != 0 ) {
	    	sprintf( estr,"Error accessing jobfile: %s\n",jobfile);
	    	exit_error(estr);
	    }

	MPI_Init (&argc, &argv);	/* starts MPI */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &size);	/* get number of processes */

	if ( rank == 0 ) {
		printf("Taskfarm running on %d processors using jobfile: %s\n", size, jobfile);
	}

	/*
	 * Below was for checking ppn per hostname - not needed though
	char *hostname = get_hostname();
	if ( hostname == NULL ) {
		exit_error( "Error accessing hostname!\n" );
	}
	if ( get_ppn(hostname) > 1 ) {
		exit_error( "More then one processor in hostcomm!\n" );
	}
	*/

	/* Get the job id */
	int jobid = get_jobid();
	if ( jobid == -1 ) {
		exit_error("Error getting jobid\n");
	}

        /* Set global_lockfile name */
	sprintf(global_lockfile, "LOCKFILE.%d",jobid);

	/* Create our local lockfile */
	//sprintf(my_lockfile, "LOCKFILE.%d.%s",jobid, hostname);
	sprintf(my_lockfile, "LOCKFILE.%d.%d",jobid, rank);
	if ( open( my_lockfile, O_WRONLY | O_CREAT | O_EXCL, 0600 ) == -1 ) {

		sprintf( estr, "Error creating lockfile: %s\n", my_lockfile );
		exit_error( estr );
	}

	char script[1024];
	memset(script, 0, sizeof(script));
	while ( next_job( jobfile, script ) == 0 ) {
		printf("Proc %d read file: %s of length %d\n", rank, script, (int)strlen(script) );
		//script="line6";
		int ret = run_job( script );
		printf("Got ret: %d\n", ret);
                memset(script, 0, sizeof(script)); /* Clear out the string */
		//sleep(rank+5);
	}

	printf("Proc %d finished\n",rank);

	// Remove the local lockfile
	unlink( my_lockfile );

	MPI_Finalize();
	return EXIT_SUCCESS;
}

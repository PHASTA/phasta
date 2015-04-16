#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
//#define OMPI_SKIP_MPICXX 1 //Added in the CMakeList.txt file
#include <mpi.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "phastaIO.h"

inline int
cscompare( const char teststring[],
	   const char targetstring[] )
{
  char* s1 = const_cast<char*>(teststring);
  char* s2 = const_cast<char*>(targetstring);

  while( *s1 == ' ') s1++;
  while( *s2 == ' ') s2++;
  while( ( *s1 )
	 && ( *s2 )
	 && ( *s2 != '?')
	 && ( tolower( *s1 )==tolower( *s2 ) ) ) {
    s1++;
    s2++;
    while( *s1 == ' ') s1++;
    while( *s2 == ' ') s2++;
  }
  if ( !( *s1 ) || ( *s1 == '?') ) return 1;
  else return 0;
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc,&argv);

  int myrank, N_procs;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &N_procs);

  FILE * pFile;
  char target[1024], pool[256];
  char * temp, * token;
  int i, j, k, N_restart_integer, N_restart_double;
  int N_geombc_double, N_geombc_integer;
  int N_steps, N_parts, N_files;

  pFile = fopen("./IO.N2O.input","r");
  if (pFile == NULL)
    printf("Error openning\n");

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_geombc_double = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_geombc_integer = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_restart_double = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_restart_integer = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_steps = atoi(temp);

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_parts = atoi(temp);

  if(myrank==0){
    printf("numpe is %d and start is %d\n",N_parts,N_steps);
  }

  fgets( target, 1024, pFile );
  token = strtok ( target, ";" );strcpy(pool,token);
  temp = strtok ( pool, ":" );temp = strtok ( NULL, ":" );
  N_files = atoi(temp);

  double ***Dfield; int ***Ifield;
  int ***paraD, ***paraI, *expectD, *expectI;
  char **fieldNameD, **fileTypeD, **dataTypeD, **headerTypeD;
  char **fieldNameI, **fileTypeI, **dataTypeI, **headerTypeI;

  int WriteLockD[N_restart_double];
  int WriteLockI[N_restart_integer];

  int nppp = N_parts/N_procs;
  int startpart = myrank * nppp +1;
  int endpart = startpart + nppp - 1;
  char gfname[64], numTemp[128];
  int iarray[10], igeom, isize;


  if (N_parts != N_procs) {
      printf("Input error: number of parts should be equal to the number of procs!\n");
      printf("Please modify the IO.O2N.input file!\n");
      return 0;
  }



  ///////////////////// reading ///////////////////////////////

  int nppf = N_parts/N_files;
  int N_geombc = N_geombc_double + N_geombc_integer;
  int readHandle, GPID;
  char fname[255],fieldtag[255];

  int irestart;

  Dfield = new double**[N_restart_double];
  Ifield = new int**[N_restart_integer];

  paraD = new int**[N_restart_double];
  paraI = new int**[N_restart_integer];

  expectD = new int[N_restart_double];
  expectI = new int[N_restart_integer];

  fieldNameD = new char*[N_restart_double];
  fileTypeD = new char*[N_restart_double];
  dataTypeD = new char*[N_restart_double];
  headerTypeD = new char*[N_restart_double];

  fieldNameI = new char*[N_restart_integer];
  fileTypeI = new char*[N_restart_integer];
  dataTypeI = new char*[N_restart_integer];
  headerTypeI = new char*[N_restart_integer];

  if (N_restart_double>0)
    for ( i = 0; i < N_restart_double; i++ )
      {
	WriteLockD[i]=0;
	Dfield[i] = new double*[nppp];

	paraD[i] = new int*[nppp];

	fieldNameD[i] = new char[128];
	fileTypeD[i] = new char[128];
	dataTypeD[i] = new char[128];
	headerTypeD[i] = new char[128];
      }

  if (N_restart_integer>0)
    for ( i = 0; i < N_restart_integer; i++ )
      {
	WriteLockI[i]=0;
	Ifield[i] = new int*[nppp];

	paraI[i] = new int*[nppp];

	fieldNameI[i] = new char[128];
	fileTypeI[i] = new char[128];
	dataTypeI[i] = new char[128];
	headerTypeI[i] = new char[128];
      }


  ///////////////////// reading ///////////////////////////////

  int N_restart = N_restart_double + N_restart_integer;
  int readHandle1;

  bzero((void*)fname,255);
  sprintf(fname,"./%d-procs_case/restart-dat.%d.%d",N_parts,N_steps,((int)(myrank/(N_procs/N_files))+1));

  //if(myrank==0){
  //  printf("Myrank is %d - Filename is %s \n",myrank,fname);
  //}

  int nfields;
  queryphmpiio(fname, &nfields, &nppf);
  //initphmpiio(&N_restart, &nppf, &N_files,&readHandle1, "write") ;//WRONG
  initphmpiio(&nfields, &nppf, &N_files, &readHandle1, "read");
  openfile(fname, "read", &readHandle1);

  for ( i = 0; i < N_restart_double; i++ )
    {
      fgets( target, 1024, pFile );
      temp = strtok( target, ";" );
      token = strtok( temp, "," );
      strcpy( fileTypeD[i], token );
      token = strtok ( NULL, "," );
      strcpy( fieldNameD[i], token );
      token = strtok ( NULL, "," );
      strcpy( dataTypeD[i], token );
      token = strtok ( NULL, "," );
      strcpy( headerTypeD[i], token );
      token = strtok ( NULL, "," );
      strcpy( numTemp, token );
      expectD[i] = atoi (numTemp);

      for (  j = 0; j < nppp; j++  )
	{
	  paraD[i][j] = new int[expectD[i]];

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  GPID = startpart + j;
	  bzero((void*)fieldtag,255);
	  sprintf(fieldtag,"%s@%d",fieldNameD[i],GPID);

          //printf("myrank %d - filedtag %s\n",myrank,fieldtag);

	  iarray[0]=-1;
	  readheader( &readHandle1,
		       fieldtag,
		       (void*)iarray,
		       &expectD[i],
		       "double",
		       "binary" );

	  if ( iarray[0]==-1 )
	      WriteLockD[i]=1;
	  if ( WriteLockD[i]==0 )
	    {
	      for ( k = 0; k < expectD[i]; k++ )
		paraD[i][j][k] = iarray[k];

	      if ( cscompare("block",headerTypeD[i]) )
		{
		  if ( expectD[i]==1)
		    isize = paraD[i][j][0];
		  else
		    isize = paraD[i][j][0] * paraD[i][j][1];

		  Dfield[i][j] = new double[isize];
		  readdatablock( &readHandle1,
				  fieldtag,
				  (void*)Dfield[i][j],
				  &isize,
				  "double",
				  "binary" );
		}

	    }
	}
    }

  for ( i = 0; i < N_restart_integer; i++ )
    {
      fgets( target, 1024, pFile );
      temp = strtok( target, ";" );
      token = strtok( temp, "," );
      strcpy( fileTypeI[i], token );
      token = strtok ( NULL, "," );
      strcpy( fieldNameI[i], token );
      token = strtok ( NULL, "," );
      strcpy( dataTypeI[i], token );
      token = strtok ( NULL, "," );
      strcpy( headerTypeI[i], token );
      token = strtok ( NULL, "," );
      strcpy( numTemp, token );
      expectI[i] = atoi (numTemp);

      for ( j = 0; j < nppp; j++ )
	{
	  paraI[i][j] = new int[expectI[i]];

	  for ( k = 0; k < 10; k++ )
	    iarray[k]=0;

	  GPID = startpart + j;
	  bzero((void*)fieldtag,255);
	  sprintf(fieldtag,"%s@%d",fieldNameI[i],GPID);
	  iarray[0]=-1;

	  //printf("Rank %d, fieldname is %s \n",myrank,fieldtag);

	  readheader( &readHandle1,
		       fieldtag,
		       (void*)iarray,
		       &expectI[i],
		       "integer",
		       "binary" );

	  if ( iarray[0]==-1)
	      WriteLockI[i]=1;
	  if ( WriteLockI[i]==0 )
	    {
	      for ( k = 0; k < expectI[i]; k++ )
		paraI[i][j][k] = iarray[k];

	      if ( cscompare("block",headerTypeI[i]) )
		{
		  if ( expectI[i]==1)
		    isize = paraI[i][j][0];
		  else
		    isize = paraI[i][j][0] * paraI[i][j][1];

		  Ifield[i][j] = new int[isize];
		  readdatablock( &readHandle1,
				  fieldtag,
				  (void*)Ifield[i][j],
				  &isize,
				  "integer",
				  "binary" );
		}
	    }
	}

    }

  closefile(&readHandle1, "write");
  finalizephmpiio(&readHandle1);

  //////////////////////////writing////////////////////////////

  int irstou;
  int magic_number = 362436;
  int* mptr = &magic_number;
  int nitems = 1;

//MR CHANGE
  bzero((void*)fname,255);
  sprintf(fname,"./%d-procs_case-1PPP",N_parts);
  if(0<mkdir(fname,0777)) { printf("ERROR - Could not create procs_case-1PPP directory\n"); return 1; }
//MR CHANGE END

  bzero((void*)fname,255);
  sprintf(fname,"./%d-procs_case-1PPP/restart.%d.%d",N_parts,N_steps,myrank+1);
  openfile(fname,"write", &irstou);

  /* writing the top ascii header for the restart file */

  writestring( &irstou,"# PHASTA Input File Version 2.0\n");
  writestring( &irstou,
		"# format \"keyphrase : sizeofnextblock usual headers\"\n");

  bzero( (void*)fname, 255 );
  writestring( &irstou, fname );

  writestring( &irstou, fname );
  writestring( &irstou,"\n");


  isize = 1;
  nitems = 1;
  iarray[ 0 ] = 1;
  writeheader( &irstou, "byteorder magic number ",
		(void*)iarray, &nitems, &isize, "integer", "binary" );

  nitems = 1;
  writedatablock( &irstou, "byteorder magic number ",
		   (void*)mptr, &nitems, "integer", "binary" );

  for ( i = 0; i < N_restart_double; i++ )
    {
      for ( j = 0; j < nppp; j++ )
	{
	  if ( WriteLockD[i] == 0 )
	    {
	      if ( cscompare("header",headerTypeD[i]) )
		{
		  bzero( (void*)fname, 255 );
		  sprintf(fname,"%s : < 0 > %d\n", fieldNameD[i],paraD[i][j][0]);
		  writestring( &irstou, fname );
		}

	      if ( cscompare("block",headerTypeD[i]) )
		{
		  if ( expectD[i]==1 )
		    isize = paraD[i][j][0];
		  else
		    isize = paraD[i][j][0] * paraD[i][j][1];

		  for ( k = 0; k < expectD[i]; k++ )
		    iarray[k] = paraD[i][j][k];

		  if ( cscompare("header",headerTypeD[i]) )
		    isize = 0;

		  writeheader( &irstou,
				fieldNameD[i],
				(void*)iarray,
				&expectD[i],
				&isize,
				"double",
				"binary");
		  writedatablock( &irstou,
				   fieldNameD[i],
				   (void*)Dfield[i][j],
				   &isize,
				   "double",
				   "binary");
		}


              if ( cscompare("block",headerTypeD[i]) )
                delete [] Dfield[i][j];
	    }
	  delete [] paraD[i][j];
	}
    }

  for ( i = 0; i < N_restart_integer; i++ )
    {
      for ( j = 0; j < nppp; j++ )
	{

	  if ( WriteLockI[i] == 0 )
	    {

	      if ( cscompare("header",headerTypeI[i]) )
		{
		  bzero( (void*)fname, 255 );
		  sprintf(fname,"%s : < 0 > %d\n", fieldNameI[i],paraI[i][j][0]);
		  writestring( &irstou, fname );
		}

	      if ( cscompare("block",headerTypeI[i]) )
		{
		  if ( expectI[i]==1 )
		    isize = paraI[i][j][0];
		  else
		    isize = paraI[i][j][0] * paraI[i][j][1];

		  for ( k = 0; k < expectI[i]; k++ )
		    iarray[k] = paraI[i][j][k];

		  writeheader( &irstou,
				fieldNameI[i],
				(void*)iarray,
				&expectI[i],
				&isize,
				"integer",
				"binary");
		  writedatablock( &irstou,
				   fieldNameI[i],
				   (void*)Ifield[i][j],
				   &isize,
				   "integer",
				   "binary");
		}

              if ( cscompare("block",headerTypeI[i]) )
                delete [] Ifield[i][j];
	    }
	  delete [] paraI[i][j];
	}
    }


  closefile( &irstou, "write" );
  MPI_Barrier(MPI_COMM_WORLD);


  if (N_restart_double>0)
    for ( i = 0; i < N_restart_double; i++ )
      {
	delete [] Dfield[i];
	delete [] paraD[i];

	delete [] fieldNameD[i];
	delete [] fileTypeD[i];
	delete [] dataTypeD[i];
	delete [] headerTypeD[i];
      }

  if (N_restart_integer>0)
    for ( i = 0; i < N_restart_integer; i++ )
      {
	delete [] Ifield[i];
	delete [] paraI[i];

	delete [] fieldNameI[i];
	delete [] fileTypeI[i];
	delete [] dataTypeI[i];
	delete [] headerTypeI[i];
      }

  delete [] Dfield;
  delete [] Ifield;

  delete [] paraD;
  delete [] paraI;

  delete [] expectD;
  delete [] expectI;

  delete [] fieldNameD;
  delete [] fileTypeD;
  delete [] dataTypeD;
  delete [] headerTypeD;

  delete [] fieldNameI;
  delete [] fileTypeI;
  delete [] dataTypeI;
  delete [] headerTypeI;

  fclose(pFile);

  if (myrank==0)
    {
      printf("\nFinished transfer, please check data using:\n");
      printf(" grep -a ': <' filename \n\n");
    }

  MPI_Finalize();

}



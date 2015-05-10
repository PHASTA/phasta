#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include "phastaIO.h"

#include <FCMangle.h>

using namespace std;

// now we create a local "unnamed" namespace which will contain all the
// functions which will be used locally within this file.  These functions
// donot have external declerations and are not visible outside this file.

namespace {
    
    map<int,char*> LastHeaderKey;
    vector<FILE*> fileArray;
    vector<bool> byte_order;
    vector<int> header_type;
    int DataSize = 0;
    bool LastHeaderNotFound = false;
    bool Wrong_Endian = false ;
    bool Strict_Error = false ;
    bool binary_format = true;


    // the caller has the responsibility to delete the returned string 

    char*
    StringStripper( const char  istring[] ) {

		int length = strlen( istring );
		char* dest = new char [ length + 1 ];
		strcpy( dest, istring );
		dest[ length ] = '\0';

		if ( char* p = strpbrk( dest, " ") ) 
			*p = '\0';
	
		return dest;
    }


    inline int 
    cscompare( const char teststring[], 
               const char targetstring[] ) {

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
    
    inline void
    isBinary( const char iotype[] ) {

        char* fname = StringStripper( iotype );
        if ( cscompare( fname, "binary" ) ) binary_format = true;
        else binary_format = false;
        delete [] fname;

    }

  
    

    inline size_t
    typeSize( const char typestring[] ) {

        char* ts1 = StringStripper( typestring );

        if ( cscompare( "integer", ts1 ) ) {
            delete [] ts1;
            return sizeof(int);
        } else if ( cscompare( "double", ts1 ) ) { 
            delete [] ts1;
            return sizeof( double );
        } else { 
            delete [] ts1;
            fprintf(stderr,"unknown type : %s\n",ts1);
            return 0;
        }
    }
    
    int
    readHeader( FILE*       fileObject,
                const char  phrase[],
                int*        params,
                int         expect ) {
        
        char* text_header;
        char* token;
        char Line[1024];
        char junk;
        bool FOUND = false ;
        int real_length;
        int skip_size, integer_value;
        int rewind_count=0;

        if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
            rewind( fileObject );
            clearerr( fileObject );
            rewind_count++;
            fgets( Line, 1024, fileObject );
        }
        
        while( !FOUND  && ( rewind_count < 2 ) )  {
            if ( ( Line[0] != '\n' ) && ( real_length = strcspn( Line, "#" )) ){
                text_header = new char [ real_length + 1 ];
                strncpy( text_header, Line, real_length );
                text_header[ real_length ] =static_cast<char>(NULL);
                token = strtok ( text_header, ":" );
                if( cscompare( phrase , token ) ) {
                    FOUND = true ;
                    token = strtok( NULL, " ,;<>" );
                    skip_size = atoi( token );
                    int i;
                    for( i=0; i < expect && ( token = strtok( NULL," ,;<>") ); i++) {
                        params[i] = atoi( token );
                    }
                    if ( i < expect ) {
                        fprintf(stderr,"Expected # of ints not found for: %s\n",phrase );
                    }
                } else if ( cscompare(token,"byteorder magic number") ) {
                    if ( binary_format ) {
                        fread((void*)&integer_value,sizeof(int),1,fileObject);
                        fread( &junk, sizeof(char), 1 , fileObject );
                        if ( 362436 != integer_value ) Wrong_Endian = true;
                    } else{
                        fscanf(fileObject, "%d\n", &integer_value );
                    }
                } else { 
                    /* some other header, so just skip over */
                    token = strtok( NULL, " ,;<>" );
                    skip_size = atoi( token );
                    if ( binary_format) 
                        fseek( fileObject, skip_size, SEEK_CUR );
                    else 
                        for( int gama=0; gama < skip_size; gama++ ) 
                            fgets( Line, 1024, fileObject );
                }
                delete [] text_header;
            }

            if ( !FOUND ) 
                if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
                    rewind( fileObject );
                    clearerr( fileObject );
                    rewind_count++;
                    fgets( Line, 1024, fileObject );
                }
        }             
        
        if ( !FOUND ) {
#ifdef DEBUG	   
            fprintf(stderr, "Warning: Cound not find: %s\n", phrase);
#endif
            return 1;
        }
        return 0;
    }	

} // end unnamed namespace
  

// begin of publicly visible functions
void 
SwapArrayByteOrder( void* array, 
                     int   nbytes, 
                     int   nItems ) {
    /* This swaps the byte order for the array of nItems each
       of size nbytes , This will be called only locally  */
    int i,j;
    unsigned char* ucDst = (unsigned char*)array;
        
    for(i=0; i < nItems; i++) {
        for(j=0; j < (nbytes/2); j++)
            std::swap( ucDst[j] , ucDst[(nbytes - 1) - j] );
        ucDst += nbytes;
    }
}

    

void 
openfile( const char filename[],
           const char mode[],
           int*  fileDescriptor ) {
    
    FILE* file=NULL ;
    *fileDescriptor = 0;
    char* fname = StringStripper( filename );
    char* imode = StringStripper( mode );

    if ( cscompare( "read", imode ) ) file = fopen(fname, "rb" );
    else if( cscompare( "write", imode ) ) file = fopen(fname, "wb" );
    else if( cscompare( "append", imode ) ) file = fopen(fname, "ab" );

    
    if ( !file ){

        fprintf(stderr,"unable to open file : %s\n",fname ) ;

    } else {

        fileArray.push_back( file );
        byte_order.push_back( false );         
        header_type.push_back( sizeof(int) );
        *fileDescriptor = fileArray.size();
        
    }
    delete [] fname;
    delete [] imode;
}

void 
closefile( int* fileDescriptor, 
            const char mode[] ) {
    
    char* imode = StringStripper( mode );

    if( cscompare( "write", imode ) 
        || cscompare( "append", imode ) ) {
        fflush( fileArray[ *fileDescriptor - 1 ] );
    } 

    fclose( fileArray[ *fileDescriptor - 1 ] );
    delete [] imode;
}

void
readheader( int* fileDescriptor,
             const char keyphrase[],
             void* valueArray,
             int*  nItems,
             const char  datatype[],
             const char  iotype[] ) {
    
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject;
    int* valueListInt;

    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
        fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
        fprintf(stderr,"openfile_ function has to be called before \n") ;
        fprintf(stderr,"acessing the file\n ") ;
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
        return;
    }

    LastHeaderKey[ filePtr ] = const_cast< char* >( keyphrase ); 
    LastHeaderNotFound = false;

    fileObject = fileArray[ filePtr ] ;
    Wrong_Endian = byte_order[ filePtr ];

    isBinary( iotype );
    typeSize( datatype );   //redundant call, just avoid a compiler warning.

    // right now we are making the assumption that we will only write integers
    // on the header line.

    valueListInt = static_cast< int* >( valueArray );
    int ierr = readHeader( fileObject ,
                           keyphrase,
                           valueListInt,
                           *nItems ) ;

    byte_order[ filePtr ] = Wrong_Endian ;

    if ( ierr ) LastHeaderNotFound = true;
    
    return ;
}

void
readdatablock( int*  fileDescriptor,
                const char keyphrase[],
                void* valueArray,
                int*  nItems,
                const char  datatype[],
                const char  iotype[] ) {
    
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject;
    char junk;
    
    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
        fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
        fprintf(stderr,"openfile_ function has to be called before \n") ;
        fprintf(stderr,"acessing the file\n ") ;
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
        return;
    }
    
    // error check..
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.    

    if ( ! cscompare( LastHeaderKey[ filePtr ], keyphrase ) ) {
        fprintf(stderr, "Header not consistant with data block\n");
        fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ] );
        fprintf(stderr, "DataBlock: %s\n ", keyphrase ); 
        fprintf(stderr, "Please recheck read sequence \n");
        if( Strict_Error ) {
            fprintf(stderr, "fatal error: cannot continue, returning out of call\n"); 
            return;
        }
    }
    
    if ( LastHeaderNotFound ) return;

    fileObject = fileArray[ filePtr ];
    Wrong_Endian = byte_order[ filePtr ];

    size_t type_size = typeSize( datatype ); 
    int nUnits = *nItems;
    isBinary( iotype );
    
    if ( binary_format ) {
        fread( valueArray, type_size, nUnits, fileObject );
        fread( &junk, sizeof(char), 1 , fileObject );
        if ( Wrong_Endian ) SwapArrayByteOrder( valueArray, type_size, nUnits );
    } else { 

        char* ts1 = StringStripper( datatype );
        if ( cscompare( "integer", ts1 ) ) {
            for( int n=0; n < nUnits ; n++ ) 
                fscanf(fileObject, "%d\n",(int*)((int*)valueArray+n) );  
        } else if ( cscompare( "double", ts1 ) ) {
            for( int n=0; n < nUnits ; n++ ) 
                fscanf(fileObject, "%lf\n",(double*)((double*)valueArray+n) );  
        }
        delete [] ts1;
    }
    
    return;
}


void 
writeheader(const int*  fileDescriptor,
	     const char keyphrase[],
	     const void* valueArray,
	     const int* nItems,
	     const int* ndataItems,
	     const char datatype[],
	     const char iotype[]  ) {
  
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject;
    
    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
        fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
        fprintf(stderr,"openfile_ function has to be called before \n") ;
        fprintf(stderr,"acessing the file\n ") ;
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
        return;
    }
    
    LastHeaderKey[ filePtr ] = const_cast< char* >( keyphrase ); 
    DataSize = *ndataItems;
    fileObject = fileArray[ filePtr ] ;
    size_t type_size = typeSize( datatype );
    isBinary( iotype );
    header_type[ filePtr ] = type_size;
    
    int _newline = ( *ndataItems > 0 ) ? sizeof( char ) : 0;
    int size_of_nextblock = 
        ( binary_format ) ? type_size*( *ndataItems )+ _newline : *ndataItems ;
    
    fprintf( fileObject, "%s : < %d > ", keyphrase, size_of_nextblock );
    for( int i = 0; i < *nItems; i++ ) 
        fprintf(fileObject, "%d ", *((int*)((int*)valueArray+i)));
    fprintf(fileObject, "\n");

    return ;
}

void 
writedatablock( const int* fileDescriptor,
                 const char keyphrase[],
                 const void* valueArray,
                 const int* nItems,
                 const char datatype[],
                 const char iotype[] ) {
    
    int filePtr = *fileDescriptor - 1;
    
    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
        fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
        fprintf(stderr,"openfile_ function has to be called before \n") ;
        fprintf(stderr,"acessing the file\n ") ;
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
        return;
    }
    
    // error check..
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.    
    
    if ( ! cscompare( LastHeaderKey[ filePtr ], keyphrase ) ) {
        fprintf(stderr, "Header not consistant with data block\n");
        fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ] );
        fprintf(stderr, "DataBlock: %s\n ", keyphrase ); 
        fprintf(stderr, "Please recheck write sequence \n");
        if( Strict_Error ) {
            fprintf(stderr, "fatal error: cannot continue, returning out of call\n"); 
            return;
        }
    }
    
    FILE* fileObject =  fileArray[ filePtr ] ;
    size_t type_size=typeSize( datatype );
    isBinary( iotype );

    if ( header_type[filePtr] != (int)type_size ) {
        fprintf(stderr,"header and datablock differ on typeof data in the block for\n"); 
        fprintf(stderr,"keyphrase : %s\n", keyphrase);
        if( Strict_Error ) {
            fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
            return;
        }
    }

    int nUnits = *nItems;

    if ( nUnits != DataSize ) {
        fprintf(stderr,"header and datablock differ on number of data items for\n"); 
        fprintf(stderr,"keyphrase : %s\n", keyphrase);
        if( Strict_Error ) {
            fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
            return;
        }
    }
 
    if ( binary_format ) {

        fwrite( valueArray, type_size, nUnits, fileObject );
        fprintf( fileObject,"\n");
        
    } else { 
        
        char* ts1 = StringStripper( datatype );
        if ( cscompare( "integer", ts1 ) ) {
            for( int n=0; n < nUnits ; n++ ) 
                fprintf(fileObject,"%d\n",*((int*)((int*)valueArray+n)));
        } else if ( cscompare( "double", ts1 ) ) {
            for( int n=0; n < nUnits ; n++ ) 
                fprintf(fileObject,"%lf\n",*((double*)((double*)valueArray+n)));
        }
        delete [] ts1;
    }	
    return ;
}

void 
writestring( int* fileDescriptor,
              const char inString[] ) {
    
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject = fileArray[filePtr] ;
    fprintf(fileObject,"%s", inString );
    return;
}

void
Gather_Headers( int* fileDescriptor,
                vector< string >& headers ) {

    FILE* fileObject;
    char Line[1024];

    fileObject = fileArray[ (*fileDescriptor)-1 ];

    while( !feof(fileObject) ) {
        fgets( Line, 1024, fileObject);
        if ( Line[0] == '#' ) {
            headers.push_back( Line );
        } else { 
            break; 
        }
    }
    rewind( fileObject );
    clearerr( fileObject );
}
void
isWrong( void ) { (Wrong_Endian) ? fprintf(stdout,"YES\n"): fprintf(stdout,"NO\n") ; }

void 
togglestrictmode( void ) { Strict_Error = !Strict_Error; }

int
isLittleEndian( void ) { 
	// this function returns a 1 if the current running architecture is
	// LittleEndian Byte Ordered, else it returns a zero 

	union  {
		long a;
		char c[sizeof( long )];
	} endianUnion;

	endianUnion.a = 1 ;

	if ( endianUnion.c[sizeof(long)-1] != 1 ) return 1 ;
	else return 0;
}


namespace PHASTA {
const char* const PhastaIO_traits<int>::type_string = "integer";
const char* const PhastaIO_traits<double>::type_string = "double";
}

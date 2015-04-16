#include <stdio.h>
#include <stdlib.h>
#include <FCMangle.h>
void
write_hessian( double* hessian, double* gradian, int* nshg ) {

    FILE* idmap = fopen( "lihessmap.dat","r");
    int* map = ( int* )malloc( (*nshg)*sizeof( int ) );
    int i,j,k;
    int x;
   
    FILE* uhess = fopen( "uhessian.dat", "w" );
    FILE* vhess = fopen( "vhessian.dat", "w" );
    FILE* whess = fopen( "whessian.dat", "w" );


    /* FILE* ugrad = fopen( "ugradian.dat", "w" ); */
    /* FILE* vgrad = fopen( "vgradian.dat", "w" ); */
    /* FILE* wgrad = fopen( "wgradian.dat", "w" ); */

    /* int ug[3] = { 1, 4, 7 }; */
    /* int vg[3] = { 2, 5, 8 }; */
    /* int wg[3] = { 3, 6, 9 }; */
    int u[6] = { 1, 10, 19, 13, 22, 25 };
    int v[6] = { 2, 11, 20, 14, 23, 26 };
    int w[6] = { 3, 12, 21, 15, 24, 27 };

    for( x=0; x < *nshg; x++ )  fscanf(idmap,"%d\n", map+x );
    fclose( idmap );

    fprintf( uhess,"%d\n", *nshg );
    fprintf( vhess,"%d\n", *nshg );
    fprintf( whess,"%d\n", *nshg );
    /* fprintf( ugrad,"%d\n", *nshg ); */
    /* fprintf( vgrad,"%d\n", *nshg ); */
    /* fprintf( wgrad,"%d\n", *nshg ); */

    for( i=0; i< *nshg; i++ ) {

        k = map[ i ]-1;

       /* for( j=0; j<3; j++ ) { */
         /* fprintf( ugrad, "%f ", gradian[i+(ug[j]-1)*(*nshg)] ); */
         /* fprintf( vgrad, "%f ", gradian[i+(vg[j]-1)*(*nshg)] ); */
         /* fprintf( wgrad, "%f ", gradian[i+(wg[j]-1)*(*nshg)] ); */
       /* } */
         
       for( j=0; j<6; j++ ) {
         fprintf( uhess, "%f ", hessian[k+(u[j]-1)*(*nshg)] );
         fprintf( vhess, "%f ", hessian[k+(v[j]-1)*(*nshg)] );
         fprintf( whess, "%f ", hessian[k+(w[j]-1)*(*nshg)] );
        }
       fprintf( uhess, "\n" ) ;
       fprintf( vhess, "\n" ) ;
       fprintf( whess, "\n" ) ;
       /* fprintf( ugrad, "\n" ) ; */
       /* fprintf( vgrad, "\n" ) ; */
       /* fprintf( wgrad, "\n" ) ; */
    }          
  
    free( map );
    fclose( uhess );
    fclose( vhess );
    fclose( whess );
    /* fclose( ugrad ); */
    /* fclose( vgrad ); */
    /* fclose( wgrad ); */
}

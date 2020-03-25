#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include<immintrin.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( ( x ) < ( y ) ? x : y )

void LoopFive( int, int, int, double *, int, double *, int, double *, int );
void LoopFour( int, int, int, double *, int, double *, int,  double *, int );
void LoopThree( int, int, int, double *, int, double *, double *, int );
void LoopTwo( int, int, int, double *, double *, double *, int );
void LoopOne( int, int, int, double *, double *, double *, int );
void Gemm_MRxNRKernel_Packed( int, double *, double *, double *, int );
void PackBlockA_MCxKC( int, int, double *, int, double * );
void PackPanelB_KCxNC( int, int, double *, int, double * );

double* Atilde = NULL;
double* Btilde = NULL;
  
void MyGemm( int m, int n, int k, double *A, int ldA,
	     double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % NR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  LoopFive( m, n, k, A, ldA, B, ldB, C, ldC );
}

void LoopFive( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{

  int max_threads = omp_get_max_threads();
  int loadbalanced_part = (n/(NC*max_threads)) * NC * max_threads;
  int remainder = n - loadbalanced_part;

  Atilde = ( double * ) _mm_malloc( MC * KC * omp_get_max_threads() * sizeof( double ), 64 );
  Btilde = ( double * ) _mm_malloc( KC * NC * omp_get_max_threads() * sizeof( double ), 64 );

  
#pragma omp parallel for
  for ( int j=0; j< loadbalanced_part; j+=NC ) {
    int jb = NC; //min( NC, n-j );    /* Last loop may not involve a full block */
    LoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }

  int remainder_per_thread = ((remainder / max_threads) / NR) * NR;
  if (remainder_per_thread == 0) remainder_per_thread = NR;

  #pragma omp parallel for
  for (int j = loadbalanced_part; j < n; j += remainder_per_thread) {
    int jb = min (remainder_per_thread, n-j );
    LoopFour (m, jb, k , A, ldA, &beta(0, j ), ldB, &gamma(0, j ), ldC );
  }

  _mm_free(Atilde);
  _mm_free(Btilde);
}

void LoopFour( int m, int n, int k, double *A, int ldA, double *B, int ldB,
	       double *C, int ldC )
{
  //double *Btilde = ( double * ) _mm_malloc( KC * NC * sizeof( double ), 64 );
  
  for ( int p=0; p<k; p+=KC ) {
    int pb = min( KC, k-p );    /* Last loop may not involve a full block */
    PackPanelB_KCxNC( pb, n, &beta( p, 0 ), ldB, &Btilde[ NC * KC * omp_get_thread_num() ] );
    LoopThree( m, n, pb, &alpha( 0, p ), ldA, &Btilde[ NC * KC * omp_get_thread_num() ], C, ldC );
  }

  //_mm_free( Btilde); 
}

void LoopThree( int m, int n, int k, double *A, int ldA, double *Btilde, double *C, int ldC )
{
  //double *Atilde = ( double * ) _mm_malloc( MC * KC * sizeof( double ), 64 );
       
  for ( int i=0; i<m; i+=MC ) {
    int ib = min( MC, m-i );    /* Last loop may not involve a full block */
    PackBlockA_MCxKC( ib, k, &alpha( i, 0 ), ldA, &Atilde[ MC * KC * omp_get_thread_num() ] );
    LoopTwo( ib, n, k, &Atilde[ MC * KC * omp_get_thread_num() ], Btilde, &gamma( i,0 ), ldC );
  }
  
  //_mm_free( Atilde);
 
}

void LoopTwo( int m, int n, int k, double *Atilde, double *Btilde, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = min( NR, n-j );
    LoopOne( m, jb, k, Atilde, &Btilde[ j*k ], &gamma( 0,j ), ldC );
  }
}

void LoopOne( int m, int n, int k, double *Atilde, double *MicroPanelB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = min( MR, m-i );
    Gemm_MRxNRKernel_Packed( k, &Atilde[ i*k ], MicroPanelB, &gamma( i,0 ), ldC );
  }
}


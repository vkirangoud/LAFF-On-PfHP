#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#include<immintrin.h>

void Gemm_MRxNRKernel( int k, double *A, int ldA, double *B, int ldB,
		double *C, int ldC )
{
  /* Declare vector registers to hold 16x3 C and load them */
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
 
  __m256d gamma_4567_0 = _mm256_loadu_pd( &gamma( 4,0 ) );
  __m256d gamma_4567_1 = _mm256_loadu_pd( &gamma( 4,1 ) );
  __m256d gamma_4567_2 = _mm256_loadu_pd( &gamma( 4,2 ) );
 
  __m256d gamma_891011_0 = _mm256_loadu_pd( &gamma( 8,0 ) );
  __m256d gamma_891011_1 = _mm256_loadu_pd( &gamma( 8,1 ) );
  __m256d gamma_891011_2 = _mm256_loadu_pd( &gamma( 8,2 ) );

  __m256d gamma_12131415_0 = _mm256_loadu_pd( &gamma( 12,0 ) );
  __m256d gamma_12131415_1 = _mm256_loadu_pd( &gamma( 12,1 ) );
  __m256d gamma_12131415_2 = _mm256_loadu_pd( &gamma( 12,2 ) );
 
   	
  for ( int p=0; p<k; p++ ){
    /* Declare vector register for load/broadcasting beta( p,j ) */
    __m256d beta_p_j;
    
    /* Declare vector registers to hold the current column of A and load
       them with the sixteen elements of that column. */
    __m256d alpha_0to15_p    = _mm256_loadu_pd( &alpha( 0,p ) );
    __m256d alpha_4567_p     = _mm256_loadu_pd( &alpha( 4,p ) );
    __m256d alpha_891011_p   = _mm256_loadu_pd( &alpha( 8,p ) );
    __m256d alpha_12131415_p = _mm256_loadu_pd( &alpha(12, p) );

    /* Load/broadcast beta( p,0 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    
    /* update the first column of C with the current column of A times
       beta ( p,0 ) */
    gamma_0123_0     = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_4567_0     = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );
    gamma_891011_0   = _mm256_fmadd_pd( alpha_891011_p, beta_p_j, gamma_891011_0 );
    gamma_12131415_0 = _mm256_fmadd_pd( alpha_12131415_p, beta_p_j, gamma_12131415_0);
    
    /* REPEAT for second and third columns of C.  Notice that the 
       current column of A need not be reloaded. */

    /* Load/broadcast beta( p,1 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 1) );
    
    /* update the second column of C with the current column of A times
       beta ( p,1 ) */
    gamma_0123_1     = _mm256_fmadd_pd( alpha_0123_p,     beta_p_j, gamma_0123_1 );
    gamma_4567_1     = _mm256_fmadd_pd( alpha_4567_p,     beta_p_j, gamma_4567_1 );
    gamma_891011_1   = _mm256_fmadd_pd( alpha_891011_p,   beta_p_j, gamma_891011_1 );
    gamma_12131415_1 = _mm256_fmadd_pd( alpha_12131415_p, beta_p_j, gamma_12131415_1 );

    /* Load/broadcast beta( p,2 ). */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 2) );
    
    /* update the third column of C with the current column of A times
       beta ( p,2 ) */
    gamma_0123_2     = _mm256_fmadd_pd( alpha_0123_p,     beta_p_j, gamma_0123_2 );
    gamma_4567_2     = _mm256_fmadd_pd( alpha_4567_p,     beta_p_j, gamma_4567_2 );
    gamma_891011_2   = _mm256_fmadd_pd( alpha_891011_p,   beta_p_j, gamma_891011_2 );
    gamma_12131415_2 = _mm256_fmadd_pd( alpha_12131415_p, beta_p_j, gamma_12131415_2 );
  }
  
  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
 
  _mm256_storeu_pd( &gamma(4,0), gamma_4567_0 );
  _mm256_storeu_pd( &gamma(4,1), gamma_4567_1 );
  _mm256_storeu_pd( &gamma(4,2), gamma_4567_2 );
 
  _mm256_storeu_pd( &gamma(8,0), gamma_891011_0 );
  _mm256_storeu_pd( &gamma(8,1), gamma_891011_1 );
  _mm256_storeu_pd( &gamma(8,2), gamma_891011_2 );

  _mm256_storeu_pd( &gamma(12,0), gamma_12131415_0 );
  _mm256_storeu_pd( &gamma(12,1), gamma_12131415_1 );
  _mm256_storeu_pd( &gamma(12,2), gamma_12131415_2 );
 
}

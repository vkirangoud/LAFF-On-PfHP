#pragma once

#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void MyGemm(int m, int n, int k, double *A, int ldA,
	double *B, int ldB, double *C, int ldC)
{
	for (int j = 0; j<n; j++)
		for (int i = 0; i<m; i++)
			for (int p = 0; p<k; p++)
				gamma(i, j) += alpha(i, p) * beta(p, j);
}


void setMatrixZero(double* A, int rows, int cols, int rs, int cs)
{
	long i, j;

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
		{
			A[i*rs + j*cs] = 0.0; 
		}
}

// Generate a random number between min and max (inclusive)
// Assumes srand() has already been called
int genRandomNumber(int min, int max)
{
	static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
	return static_cast<int>(rand() * fraction * (max - min + 1) + min);
}

void gen_random_matrix(double* A, int rows, int cols, int rs, int cs)
{
	long i, j;
	//unsigned short xsubi[3] = { 918, 729, 123 };
	srand(static_cast<unsigned int>(time(0))); // set initial seed value to system clock

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
		{
			A[i*rs + j*cs] = (double)genRandomNumber(1, 8); // 1 + nrand48(xsubi) % 9;
		}
}

int isMatrixMatch(double* ref, double* mat, int rows, int cols, int cs)
{
	long i, j;
	bool flag = false;
	for (j = 0; j < cols; j++)
	{
		for (i = 0; i < rows; i++)
		{
			if (ref[i + j*cs] != mat[i + j*cs])
			{
				std::cout << "Mismatch at " << "( " << i << ", " << j << ")" << "\n";
				return 0;
			}
		}
	}

	return 1;
}


void writeMatrix(FILE* fin, double* A, int rows, int cols, int rs, int cs)
{
	long i, j;
	for (j = 0; j < cols; j++)
	{
		for (i = 0; i < rows; i++)
		{
			fprintf(fin, "%f ", A[i*rs + j*cs]);
		}
	}
}


void readMatrices(char* str, double** A, double** B, double** C, int* M, int* K, int* N, int* lda, int* ldb, int* ldc)
{
	FILE* fin = NULL;
	fin = fopen(str, "rt");
	if (fin == NULL)
	{
		std::cout << "Error opening the matrix file\n";
		exit(1);
	}
	fscanf(fin, "%d\t %d\t %d\t %d\t %d\t %d\n", M, K, N, lda, ldb, ldc);

	double* ap = NULL;
	double* bp = NULL;
	double* cp = NULL;

	int cs_a = *lda;
	int cs_b = *ldb;
	int cs_c = *ldc;

	int m = *M;
	int n = *N;
	int k = *K;

	int nelems_A = (m - 1) * 1 + (k - 1) * cs_a + 1; // rs_a = 1
	int nelems_B = (k - 1) * 1 + (n - 1) * cs_b + 1; // rs_b = 1
	int nelems_C = (m - 1) * 1 + (n - 1) * cs_c + 1; // rs_c = 1

	ap = (double*)malloc(sizeof(double) * nelems_A);
	if (ap == NULL) { printf("Error allocation memory A \n"); exit(1); }

	bp = (double*)malloc(sizeof(double) * nelems_B);
	if (bp == NULL) { printf("Error allocation memory B \n"); exit(1); }

	cp = (double*)malloc(sizeof(double) * nelems_C);
	if (cp == NULL) { printf("Error allocation memory C \n"); exit(1); }

	// read matrix A
	long i, j;
	for (j = 0; j < k; j++)
	{
		for (i = 0; i < m; i++)
		{
			fscanf(fin, "%lf ", &ap[i + j * (*lda)]);
		}
	}

	// read matrix B
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < k; i++)
		{
			fscanf(fin, "%lf ", &bp[i + j * (*ldb)]);
		}
	}

	// read matrix C
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			fscanf(fin, "%lf ", &cp[i + j * (*ldc)]);
		}
	}

	*A = ap;
	*B = bp;
	*C = cp;

	return;
}


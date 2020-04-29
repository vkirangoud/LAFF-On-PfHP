#include <iostream>
#include <cstdlib> // for rand() and srand()
#include <ctime> // for time()

//#include "aocl_MatricesReadWrite.h"
#include "aocl_matrices.h"

int main(int argc, char*argv[])
{
	if (argc < 8)
	{
		std::cout << " usage: m k n cs_a cs_b cs_c name of matrix \n";
		exit(1);
	}

	int m = atoi(argv[1]);
	int k = atoi(argv[2]);
	int n = atoi(argv[3]);
	int cs_a = atoi(argv[4]);
	int cs_b = atoi(argv[5]);
	int cs_c = atoi(argv[6]);

	double* ap = NULL;
	double* bp = NULL;
	double* cp = NULL;

	FILE* fp = NULL;
	fp = fopen(argv[7], "wt");
	if (fp == NULL)
	{
		std::cout << "error opening output file\n";
		exit(1);
	}

	fprintf(fp, "%d\t %d\t %d\t %d\t %d\t %d\n", m, k, n, cs_a, cs_b, cs_c);

	int nelems_A = (m - 1) * 1 + (k - 1) * cs_a + 1; // rs_a = 1
	int nelems_B = (k - 1) * 1 + (n - 1) * cs_b + 1; // rs_b = 1
	int nelems_C = (m - 1) * 1 + (n - 1) * cs_c + 1; // rs_c = 1

	ap = (double*)malloc(sizeof(double) * nelems_A);
	if (ap == NULL) { printf("Error allocation memory A \n"); exit(1); }

	bp = (double*)malloc(sizeof(double) * nelems_B);
	if (bp == NULL) { printf("Error allocation memory B \n"); exit(1); }

	cp = (double*)malloc(sizeof(double) * nelems_C);
	if (cp == NULL) { printf("Error allocation memory C \n"); exit(1); }

	gen_random_matrix(ap, m, k, 1, cs_a);

	gen_random_matrix(bp, k, n, 1, cs_b);

	gen_random_matrix(cp, m, n, 1, cs_c);

	writeMatrix(fp, ap, m, k, 1, cs_a);
	writeMatrix(fp, bp, k, n, 1, cs_b);
	writeMatrix(fp, cp, m, n, 1, cs_c);

	fclose(fp);

#if 1

	// verify Matrix read
	double* A = NULL;
	double* B = NULL;
	double* C = NULL;

	int M = 0;
	int N = 0;
	int K = 0;
	int ldA = 0;
	int ldB = 0;
	int ldC = 0;

	readMatrices(argv[7], &A, &B, &C, &M, &K, &N, &ldA, &ldB, &ldC);
	if ((M != m) || (N != n) || (K != k) || (ldA != cs_a) || (ldB != cs_b) || (ldC != cs_c))
	{
		std::cout << "Mismatch in dimensions \n";
		std::cout << "out: " << M << "\t" << K << "\t" << N << "\t" << ldA << "\t" << ldB << "\t" << ldC << "\n";
		std::cout << "in: " << m << "\t" << k << "\t" << n << "\t" << cs_a << "\t" << cs_b << "\t" << cs_c << "\n";

		free(A);
		free(B);
		free(C);
		exit(1);
	}

	// Check Matrix A
	if (isMatrixMatch(ap, A, M, K, ldA))
	{
		std::cout << "Matrix A matches \n";
	}
	else
	{
		std::cout << "Matrix A mis-matches\n";
	}

	// Check Matrix B
	if (isMatrixMatch(bp, B, K, N, ldB))
	{
		std::cout << "Matrix B matches \n";
	}
	else
	{
		std::cout << "Matrix B Mis-match\n";
	}

	// Check Matrix C
	if (isMatrixMatch(cp, C, M, N, ldC))
	{
		std::cout << "Matrix C matches \n";
	}
	else
	{
		std::cout << "Matrix C Mis-match\n";
	}

	free(A);
	free(B);
	free(C);

#endif

	free(ap);
	free(bp);
	free(cp);


	return 0;
}
#include <stdio.h>
#include <examples.h>
#include <omp.h>
#include <compiler.h>
#include <common.h>

#define N 100
#define M 100
#define P 100

void example1(void)
{
    int nthreads;
    int tid;

    printf("EXAMPLE 1\n");

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);

        /* Only master thread does this */
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }
    }

    printf("\n");
}

void example2(void)
{
    int nthreads;
    int tid;
    size_t i;

    int t1[N];
    int t2[N];
    int t3[N];

    printf("EXAMPLE 2\n");
    for (i = 0; i < ARRAY_SIZE(t1); ++i)
    {
        t1[i] = (int)i + 1;
        t2[i] = (int)i + 1;
    }

#pragma omp parallel shared(t1, t2, t3, nthreads) private(i, tid)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }

        printf("Thread %d starting...\n", tid);

        #pragma omp for
        for (i = 0; i < ARRAY_SIZE(t1); ++i)
        {
            t3[i] = t1[i] + t2[i];
            printf("Thread %d: t3[ %zu ]= %d\n", tid, i, t3[i]);
        }

    }

    printf("\n");
}

void example3(void)
{
    int t[N] = {0};
    int j;
    int x;
    size_t i;

    printf("EXAMPLE 3\n");
    j = 1;
#pragma omp parallel for firstprivate(j)
    for(i = 0; i < ARRAY_SIZE(t); ++i)
        t[i] += j * (int)i + 1;

#pragma omp parallel for firstprivate(j) ordered
    for(i = 0; i < ARRAY_SIZE(t); ++i)
    {
        t[i] += j * (int)i + 1;

#pragma omp ordered
        printf("t[%zu] = %d\n", i, t[i]);
    }

#pragma omp parallel for lastprivate(x)
    for (i = 1; i < ARRAY_SIZE(t); ++i)
        x = t[i] + t[i - 1];

    printf("X = %d\n", x);
    printf("\n");
}

void example4(void)
{
    size_t i;
    size_t j;
    size_t k;
    int temp;

    int t1[M][N];
    int t2[N][P];
    int t3[M][P];

    printf("EXAMPLE 4\n");

    for (i = 0; i < M; ++i)
        for (j = 0; j < N; ++j)
            t1[i][j] = (int)i + 1;

    for (i = 0; i < N; ++i)
        for (j = 0; j < P; ++j)
            t2[i][j] = (int)i + 1;


#pragma omp parallel for schedule(dynamic) private(temp, i, j, k)
    for (i = 0; i < M; ++i)
	    for (j = 0; j < P; ++j)
        {
            temp = 0;
            for(k = 0; k < N; k++)
                temp += t1[i][k] * t2[k][j];

			t3[i][j] = temp;
		}

    (void)t3;
    printf("\n");
}
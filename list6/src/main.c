#include <ecm.h>
#include <stdio.h>
#include <gmp.h>
#include <compiler.h>
#include <log.h>
#include <darray.h>
#include <string.h>
#include <common.h>
#include <time.h>
#include <inttypes.h>
#include <omp.h>

#define BASE 10

static int help(void);

___before_main___(1) void init(void);
___after_main___(1) void deinit(void);

___before_main___(1) void init(void)
{
    log_init(stdout, NO_LOG_TO_FILE);
}

___after_main___(1) void deinit(void)
{
    log_deinit();
}

static int help(void)
{
    (void)printf("Program to factor number\n"
                 "NEED 1 argument\n"
                 "n - number to factor\n"
                 "Output factors of n\n");

    return 0;
}

static Darray *sieve(uint32_t n)
{
    Darray *primes;
    uint32_t *sieve_array;

    const size_t sieve_size = n + 1;
    size_t i;
    size_t j;
    uint32_t temp;

    TRACE();

    primes = darray_create(DARRAY_UNSORTED, 0, sizeof(uint32_t), NULL, NULL);
    if (primes == NULL)
        ERROR("darray_create error\n", NULL);

    sieve_array = (uint32_t *)malloc(sizeof(uint32_t) * sieve_size);
    if (sieve_array == NULL)
        ERROR("malloc error\n", NULL);

    (void)memset(sieve_array, 0, sizeof(uint32_t) * sieve_size);
    sieve_array[0] = 1;
    sieve_array[1] = 1;

    for (i = 2; i < sieve_size; ++i)
        if (sieve_array[i] == 0)
        {
            temp = (uint32_t)i;
            darray_insert(primes, (const void *)&temp);
            for (j = i; j < sieve_size; j += i)
                sieve_array[j] = 1;
        }

    FREE(sieve_array);
    return primes;
}


int main(int argc, char **argv)
{
    mpz_t n;
    mpz_t factor;
    mpz_t temp;
    Darray *primes;

    uint32_t limit = 100000;
    uint32_t counter;

    if (argc < 2)
        return help();

    mpz_init(n);
    mpz_set_str(n, argv[1], BASE);

    srand((unsigned int)time(NULL));

    gmp_printf("Trying to factor %Zd\n", n);
    if (mpz_probab_prime_p(n, 10) > 0)
    {
        gmp_printf("%Zd is prime\n", n);
        mpz_clear(n);
        mpz_clear(factor);

        return 0;
    }

    primes = sieve(limit);
    if (primes == NULL)
        FATAL("Sieve error\n");

    mpz_init(temp);
    while (mpz_probab_prime_p(n, 10) == 0)
    {
        counter = 0;
        #pragma omp parallel private(factor) shared(limit, counter, n, primes, temp)
        {
            mpz_init(factor);
            mpz_set_ui(factor, 1);
            do {
                #pragma omp critical
                {
                if (counter < 10)
                    limit = 5000;
                else if (counter < 50)
                    limit = 10000;
                else if (counter < 100)
                    limit = 50000;
                else
                    limit = 100000;

                    ++counter;
                    //printf("Try %" PRIu32 " B = %" PRIu32 "\n", counter, limit);
                }
            } while (mpz_probab_prime_p(n, 10) == 0 && lenstra_ecm(n, primes, limit, factor));
            #pragma omp critical
            {
                mpz_mod(temp, n, factor);
                if (mpz_cmp_ui(factor, 1) && mpz_cmp_ui(temp, 0) == 0 && mpz_cmp(n, factor))
                {
                    mpz_div(n, n, factor);
                    gmp_printf("FACTOR: %Zd\n", factor);
                }
                mpz_clear(factor);
            }
        }
    }

    gmp_printf("FACTOR = %Zd\n", n);
    mpz_clear(n);
    mpz_clear(factor);

    return 0;
}
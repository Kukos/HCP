#include <pohling.h>
#include <stdio.h>
#include <gmp.h>
#include <compiler.h>
#include <log.h>
#include <common.h>
#include <stdlib.h>
#include <string.h>

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
    (void)printf("Program solving g^x = h (mod p)\n"
                 "NEED 3 arguments\n"
                 "g - generator\n"
                 "h - result of power\n"
                 "p - prime\n"
                 "list of f e such that (p - 1) = PRODUCT[fi^fei]\n"
                 "Output x\n");

    return 0;
}

int main(int argc, char **argv)
{
    mpz_t g;
    mpz_t h;
    mpz_t p;
    mpz_t x;

    mpz_t *factors;
    mpz_t *exponents;
    int res;
    int ret;

    size_t f_len;
    size_t i;

    /* for checking inputs */
    mpz_t ord_p;
    mpz_t temp;
    mpz_t prod;

    if (argc < 4)
        return help();

    mpz_init(g);
    mpz_init(h);
    mpz_init(p);
    mpz_init(x);

    /* set g h p */
    mpz_set_str(g, argv[1], BASE);
    mpz_set_str(h, argv[2], BASE);
    mpz_set_str(p, argv[3], BASE);

    /* load (p-1) factors */
    f_len = (size_t)(argc - 4) >> 1;
    factors = (mpz_t *)malloc(sizeof(mpz_t) * f_len);
    if (factors == NULL)
        FATAL("malloc error\n");

    exponents = (mpz_t *)malloc(sizeof(mpz_t) * f_len);
    if (exponents == NULL)
        FATAL("malloc error\n");

    mpz_init(ord_p);
    mpz_sub_ui(ord_p, p, 1);

    mpz_init(temp);

    mpz_init(prod);
    mpz_set_ui(prod, 1);

    argv += 4;
    for (i = 0; i < f_len; ++i)
    {
        mpz_init(factors[i]);
        mpz_set_str(factors[i], *(argv++), BASE);

        mpz_init(exponents[i]);
        mpz_set_str(exponents[i], *(argv++), BASE);

        mpz_powm(temp, factors[i], exponents[i], p);
        mpz_mul(prod, prod, temp);
    }

    if (mpz_cmp(ord_p, prod))
        FATAL("Ord p has incorrect factors\n");

    (void)gmp_printf("Trying find x such that %Zd^x = %Zd (mod %Zd)\n", g, h, p);
    res = pohling_discrete_log(g, h, p, factors, exponents, f_len, x);

    if (res)
        (void)printf("FAILED\n");
    else
        (void)gmp_printf("X = %Zd\n", x);

    mpz_powm(x, g, x, p);
    if (mpz_cmp(x, h) == 0)
    {
        (void)printf("SUCCESS!!!\n");
        ret = 0;
    }
    else
    {
        (void)printf("FAILED!!!\n");
        ret = 1;
    }

    mpz_clear(g);
    mpz_clear(h);
    mpz_clear(p);
    mpz_clear(x);

    for (i = 0; i < f_len; ++i)
    {
        mpz_clear(factors[i]);
        mpz_clear(exponents[i]);
    }

    FREE(exponents);
    FREE(factors);

    return ret;
}
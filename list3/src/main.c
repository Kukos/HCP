#include <pollard.h>
#include <stdio.h>
#include <gmp.h>
#include <compiler.h>
#include <log.h>

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
                 "p - strong prime such that exist q that p = 2q + 1\n"
                 "Output x\n");

    return 0;
}

int main(int argc, char **argv)
{
    mpz_t g;
    mpz_t h;
    mpz_t p;
    mpz_t x;

    int res;
    int ret;

    if (argc < 4)
        return help();

    mpz_init(g);
    mpz_init(h);
    mpz_init(p);
    mpz_init(x);

    mpz_set_str(g, argv[1], BASE);
    mpz_set_str(h, argv[2], BASE);
    mpz_set_str(p, argv[3], BASE);

    (void)gmp_printf("Trying find x such that %Zd^x = %Zd (mod %Zd)\n", g, h, p);
    res = pollard_rho_parallel_dicsrete_log(g, h, p, x);

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

    return ret;
}
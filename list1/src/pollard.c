#include <pollard.h>
#include <log.h>

#define SUBGRP 3

/*
    f(x) = { x*g mod p iff subgrp == 0
             x*h mod p iff subgrp == 1
             x^2 mod p iff subgrp == 2
           }

    PARAMS
    @IN subgrp - subgroup
    @IN x - x
    @IN g - generator
    @IN h - resul of power
    @IN p - strong prime

    RETURN
    This is a void function
*/
static void f_f(int subgrp, mpz_t x, const mpz_t g, const mpz_t h, const mpz_t p);

/*
    g(a) = { a + 1 mod q iff subgrp == 0
             a           iff subgrp == 1
             2a mod q    iff subgrp == 2
           }

    PARAMS
    @IN subgrp - subgroup
    @IN a - a
    @IN q - prime from strong prime

    RETURN
    This is a void function
*/
static void g_f(int subgrp, mpz_t a, const mpz_t q);

/*
    h(b) = { b           iff subgrp == 0
             b + 1 mod q iff subgrp == 1
             ba mod q    iff subgrp == 2
           }

    PARAMS
    @IN subgrp - subgroup
    @IN b - b
    @IN q - prime from strong prime

    RETURN
    This is a void function
*/
static void h_f(int subgrp, mpz_t b, const mpz_t q);

/*
    Single step of Pollard rho algo

    PARAMS
    @IN x - x
    @IN a - a
    @IN b - b
    @IN g - generator
    @IN h - res of power
    @IN p - strong prime
    @IN q - prime from strong prime

    RETURN
    This is a void function
*/
static void single_step(mpz_t x, mpz_t a, mpz_t b, const mpz_t g, const mpz_t h, const mpz_t p, const mpz_t q);

static void f_f(int subgrp, mpz_t x, const mpz_t g, const mpz_t h, const mpz_t p)
{
    TRACE();

    switch (subgrp)
    {
        case 0:
        {
            /* x = x * g mod p */
            mpz_mul(x, x, g);
            mpz_mod(x, x, p);
            break;
        }
        case 1:
        {
            /* x = x * h mod p */
            mpz_mul(x, x, h);
            mpz_mod(x, x, p);
            break;
        }
        case 2:
        {
            /* x = x ^ 2 mod p */
            mpz_mul(x, x, x);
            mpz_mod(x, x, p);
            break;
        }
        default:
        {
            LOG("Undefined subgrp\n");
            break;
        }
    }
}

static void g_f(int subgrp, mpz_t a, const mpz_t q)
{
    TRACE();

    switch (subgrp)
    {
        case 0:
        {
            /* a = a + 1 mod q */
            mpz_add_ui(a, a, 1);
            mpz_mod(a, a, q);
            break;
        }
        case 1:
        {
            /* a = a */
            break;
        }
        case 2:
        {
            /* a = 2a mod q */
            mpz_mul_ui(a, a, 2);
            mpz_mod(a, a, q);
            break;
        }
        default:
        {
            LOG("Undefined subgrp\n");
            break;
        }
    }
}

static void h_f(int subgrp, mpz_t b, const mpz_t q)
{
    TRACE();

    switch (subgrp)
    {
        case 0:
        {
            /* b = b*/
            break;
        }
        case 1:
        {
            /* b = b + 1 mod q */
            mpz_add_ui(b, b, 1);
            mpz_mod(b, b, q);
            break;
        }
        case 2:
        {
            /* b = 2b mod q */
            mpz_mul_ui(b, b, 2);
            mpz_mod(b, b, q);
            break;
        }
        default:
        {
            LOG("Undefined subgrp\n");
            break;
        }
    }
}

static void single_step(mpz_t x, mpz_t a, mpz_t b, const mpz_t g, const mpz_t h, const mpz_t p, const mpz_t q)
{
    int subgrp;
    mpz_t temp;

    mpz_init(temp);

    mpz_mod_ui(temp, x, SUBGRP);
    subgrp = (int)mpz_get_ui(temp);
    mpz_clear(temp);

    f_f(subgrp, x, g, h, p);
    g_f(subgrp, a, q);
    h_f(subgrp, b, q);
}


int pollard_rho_dicsrete_log(const mpz_t g, const mpz_t h, const mpz_t p, mpz_t x)
{
    mpz_t q; /* prime from strong prime */

    /* hedgehog step variables */
    mpz_t a;
    mpz_t b;

    /* rabbis step variables */
    mpz_t X;
    mpz_t A;
    mpz_t B;

    mpz_t i;

    mpz_t r;

    TRACE();

    mpz_init(q);
    mpz_init(a);
    mpz_init(b);
    mpz_init(X);
    mpz_init(A);
    mpz_init(B);

    /* p is strong prime, so exist q such that p = 2q + 1 --> q = (p - 1) / 2 */
    mpz_sub_ui(q, p, 1);
    mpz_div_ui(q, q, 2);

    /* firstly x = g*h */
    mpz_mul(x, g, h);

    /* a = 1, b = 1 */
    mpz_set_ui(a, 1);
    mpz_set_ui(b, 1);

    /* A = a, B = b, X = x */
    mpz_set(A, a);
    mpz_set(B, b);
    mpz_set(X, x);

    /* for i = 1; i < p; ++i */
    for (mpz_init(i), mpz_set_ui(i, 1); mpz_cmp(i, p) < 0; mpz_add_ui(i, i, 1))
    {
        single_step(x, a, b, g, h, p, q);
        single_step(X, A, B, g, h, p, q);
        single_step(X, A, B, g, h, p, q);

        if (mpz_cmp(X, x) == 0)
            break;
    }

    /* r = b - B */
    mpz_init(r);
    mpz_sub(r, b, B);
    if (mpz_cmp_ui(r, 0) == 0)
        ERROR("FAILURE R == 0\n", 1);

    /* x = r^-1 * (A - a) mod q */
    mpz_sub(a, A, a);
    mpz_invert(x, r, q);
    mpz_mul(x, x, a);
    mpz_mod(x, x, q);

    mpz_clear(i);
    mpz_clear(q);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(X);
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(r);

    return 0;
}
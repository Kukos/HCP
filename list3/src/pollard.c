#include <pollard.h>
#include <log.h>
#include <omp.h>
#include <time.h>
#include <stdbool.h>
#include <darray.h>
#include <common.h>
#include <stdlib.h>

#define POLLARD_TRESHOLD 40
#define POLLARD_RAND_MAX 16

typedef struct Pollard_triple
{
    mpz_t x;
    mpz_t a;
    mpz_t b;
} Pollard_triple;

/*
    Create Pollard triple

    PARAMS
    @IN x - x
    @IN a - a
    @IN b - b

    RETURN
    NULL iff failure
    Pointer to new triple iff success
*/
static Pollard_triple *pollard_triple_create(const mpz_t x, const mpz_t a, const mpz_t b);

/*
    Destroy Pollard Triple

    PARAMS
    @IN pt - pointer to Pollared tripple

    RETURN
    This is a void function
*/
static void pollard_triple_destroy(Pollard_triple *pt);

/*
    Std compare function for pollard triple

    PARAMS
    @IN pt1 - (void *)&Pollard_triple *
    @IN pt2 - (void *)&Pollard_triple *

    RETURN
    -1 iff pt1 < pt2
    1 iff pt1 > pt2
    0 iff pt1 = pt2
*/
static int pollard_triple_cmp(const Pollard_triple *pt1, const Pollard_triple *pt2);

/*
    Wrappers
*/
static int pollard_triple_cmp_wrapper(const void *a, const void *b);
static void pollard_triple_destroy_wrapper(void *p);

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

static Pollard_triple *pollard_triple_create(const mpz_t x, const mpz_t a, const mpz_t b)
{
    Pollard_triple *pt;

    pt = (Pollard_triple *)malloc(sizeof(Pollard_triple));
    if (pt == NULL)
        ERROR("Malloc error\n", NULL);

    mpz_init(pt->x);
    mpz_init(pt->a);
    mpz_init(pt->b);

    mpz_set(pt->x, x);
    mpz_set(pt->a, a);
    mpz_set(pt->b, b);

    return pt;
}

static void pollard_triple_destroy(Pollard_triple *pt)
{
    if (pt == NULL)
        return;

    mpz_clear(pt->x);
    mpz_clear(pt->a);
    mpz_clear(pt->b);

    FREE(pt);
}

static int pollard_triple_cmp(const Pollard_triple *pt1, const Pollard_triple *pt2)
{
    return mpz_cmp(pt1->x, pt2->x);
}

static void pollard_triple_destroy_wrapper(void *p)
{
    pollard_triple_destroy(*(Pollard_triple **)p);
}

static int pollard_triple_cmp_wrapper(const void *a, const void *b)
{
    Pollard_triple *pt1 = *(Pollard_triple **)a;
    Pollard_triple *pt2 = *(Pollard_triple **)b;

    return pollard_triple_cmp(pt1, pt2);
}

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


int pollard_rho_parallel_dicsrete_log(const mpz_t g, const mpz_t h, const mpz_t p, mpz_t res)
{
    mpz_t q; /* prime from strong prime */

    mpz_t a;
    mpz_t b;

    mpz_t i;

    mpz_t r;
    mpz_t x;

    mpz_t temp_g;
    mpz_t temp_h;

    gmp_randstate_t r_state;
    bool done = false;
    Pollard_triple *pt;
    Pollard_triple *pt_p;

    Darray *array;

    TRACE();

    mpz_init(q);

    /* p is strong prime, so exist q such that p = 2q + 1 --> q = (p - 1) / 2 */
    mpz_sub_ui(q, p, 1);
    mpz_div_ui(q, q, 2);

    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, (unsigned long)time(NULL));

    array = darray_create(DARRAY_SORTED, 0, sizeof(Pollard_triple *), pollard_triple_cmp_wrapper, pollard_triple_destroy_wrapper);
    if (array == NULL)
        ERROR("malloc error\n", 1);

#pragma omp parallel private(x, a, b, i, r, temp_g, temp_h, pt, pt_p) shared(g, h, p, r_state, res, done, array)
    {
        mpz_init(x);
        mpz_init(a);
        mpz_init(b);
        mpz_init(i);
        mpz_init(temp_g);
        mpz_init(temp_h);
        mpz_init(r);

        while (!done)
        {
            /* rand a and b */
            mpz_urandomb(a, r_state, POLLARD_RAND_MAX);
            mpz_urandomb(b, r_state, POLLARD_RAND_MAX);

            mpz_powm(temp_g, g, a, p);
            mpz_powm(temp_h, h, b, p);

            /* firstly x = g^a*h^b */
            mpz_mul(x, temp_g, temp_h);

            /* for i = 1; i < p; ++i */
            for (mpz_set_ui(i, 1); mpz_cmp(i, p) < 0; mpz_add_ui(i, i, 1))
            {
                single_step(x, a, b, g, h, p, q);

                /* x is distingish point */
                if (mpz_sizeinbase(x, 2) < POLLARD_TRESHOLD)
                    break;
            }

#pragma omp critical
            {
                /* check common array of pollard triples */
                pt = pollard_triple_create(x, a, b);
                if (darray_get_num_entries(array) > 0 && darray_search_first(array, (void *)&pt, (void *)&pt_p) != -1)
                {
                    LOG("Collision\n");
                    pollard_triple_destroy(pt);

                    /* r = b - B */
                    mpz_sub(r, b, pt_p->b);
                    if (mpz_cmp_ui(r, 0))
                    {
                        /* x = r^-1 * (A - a) mod q */
                        mpz_sub(a, pt_p->a, a);
                        if (mpz_invert(x, r, q)) /* inverstion exists */
                        {
                            mpz_mul(x, x, a);
                            mpz_mod(x, x, q);

                            LOG("Inversion is correct, finish work\n");

                            /* copy result and finish work on all threads */
                            mpz_set(res, x);
                            done = true;
                        }
                    }
                }
                else
                    darray_insert(array, (void *)&pt);
            }
        }

        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(r);
        mpz_clear(i);
        mpz_clear(temp_g);
        mpz_clear(temp_h);
        mpz_clear(x);
    }

    gmp_randclear(r_state);
    mpz_clear(q);

    darray_destroy_with_entries(array);

    return 0;
}
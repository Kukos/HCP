#include <pollard.h>
#include <log.h>
#include <omp.h>
#include <time.h>
#include <stdbool.h>
#include <darray.h>
#include <common.h>
#include <stdlib.h>
#include <hash.h>
#include <string.h>

#define POLLARD_TRESHOLD 40

typedef enum KANGAROO_TYPE
{
    KANGAROO_WILD,
    KANGAROO_TAME
} kangaroo_t;

typedef struct Pollard_triple
{
    kangaroo_t type;
    mpz_t dist;
    mpz_t pos; /* jump pos */
} Pollard_triple;

/*
    Calculate max jumps from formula:
    First r than (2^r - 1) / r > beta

    PARAMS
    @IN beta - beta

    RETURN
    MAX JUMPS
*/
static ___inline___ unsigned long  calculate_max_jumps(const mpz_t beta);

/*
    Create Pollard triple

    PARAMS
    @IN type - type
    @IN dist - dist
    @IN pos - pos

    RETURN
    NULL iff failure
    Pointer to new triple iff success
*/
static Pollard_triple *pollard_triple_create(kangaroo_t type, const mpz_t dist, const mpz_t pos);

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


static ___inline___ unsigned long  calculate_max_jumps(const mpz_t beta)
{
    unsigned long r = 1;
    mpz_t res;

    mpz_init(res);

    do {
        /* res = (2^r - 1) / r */
        mpz_ui_pow_ui(res, 2, r);
        mpz_sub_ui(res, res, 1);
        mpz_div_ui(res, res, r);

        ++r;
    } while (mpz_cmp(res, beta) < 0);

    mpz_clear(res);

    return r - 2;
}

static Pollard_triple *pollard_triple_create(kangaroo_t type, const mpz_t dist, const mpz_t pos)
{
    Pollard_triple *pt;

    pt = (Pollard_triple *)malloc(sizeof(Pollard_triple));
    if (pt == NULL)
        ERROR("Malloc error\n", NULL);

    mpz_init(pt->pos);
    mpz_init(pt->dist);

    pt->type = type;
    mpz_set(pt->pos, pos);
    mpz_set(pt->dist, dist);

    return pt;
}

static void pollard_triple_destroy(Pollard_triple *pt)
{
    if (pt == NULL)
        return;

    mpz_clear(pt->pos);
    mpz_clear(pt->dist);

    FREE(pt);
}

static int pollard_triple_cmp(const Pollard_triple *pt1, const Pollard_triple *pt2)
{
    return mpz_cmp(pt1->pos, pt2->pos);
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

int pollard_lambda_parallel_dicsrete_log(const mpz_t g, const mpz_t h, const mpz_t p, mpz_t res)
{
    const unsigned int nproc = (unsigned int)omp_get_max_threads();
    unsigned long r;

    unsigned long i;
    int index;
    char *str;

    mpz_t a;
    mpz_t b;
    mpz_t order_g;

    mpz_t beta;
    mpz_t v;

    mpz_t *dists;
    mpz_t *jumps;

    Darray *set;
    Pollard_triple pt;
    Pollard_triple *pt_p = &pt;
    Pollard_triple *triple;

    kangaroo_t type;
    mpz_t dist;
    mpz_t pos;
    mpz_t x;
    mpz_t step;

    bool finish = false;

    if (mpz_cmp(g, h) == 0)
    {
        mpz_set_ui(res, 1);
        return 0;
    }

    set = darray_create(DARRAY_SORTED, 0, sizeof(Pollard_triple *), pollard_triple_cmp_wrapper, pollard_triple_destroy_wrapper);
    if (set == NULL)
        ERROR("malloc error\n", 1);

    mpz_init(a);
    mpz_init(b);
    mpz_init(order_g);
    mpz_init(beta);
    mpz_init(v);

    /* order(generator(p)) = p - 1 */
    mpz_set(order_g, p);
    mpz_sub_ui(order_g, order_g, 1);

    /* range [0, order_g] */
    mpz_set_ui(a, 0);
    mpz_set(b, order_g);

    /* beta = (nproc * sqrt(b - a) / 4) */
    mpz_set(beta, b);
    mpz_sub(beta, beta, a);
    mpz_sqrt(beta, beta);
    mpz_mul_ui(beta, beta, nproc);
    mpz_div_ui(beta, beta, 4);

    /* v = beta / nproc / 2 */
    mpz_div_ui(v, beta, nproc >> 1);

    r = calculate_max_jumps(beta);

    dists = (mpz_t *)malloc(sizeof(mpz_t) * r);
    if (dists == NULL)
        ERROR("malloc error\n", 1);

    jumps = (mpz_t *)malloc(sizeof(mpz_t) * r);
    if (jumps == NULL)
        ERROR("malloc error\n", 1);

    for (i = 0; i < r; ++i)
    {
        mpz_init(dists[i]);
        mpz_init(jumps[i]);

        mpz_ui_pow_ui(dists[i], 2, i);
        mpz_powm(jumps[i], g, dists[i], p);
    }

#pragma omp parallel private(dist, pos, type, index, str, x, step) shared(a, b, g, h, p, set, jumps, dists, r, v, res, finish)
{
    if (ODD(omp_get_thread_num()))
        type = KANGAROO_WILD;
    else
        type = KANGAROO_TAME;

    mpz_init(dist);
    mpz_init(pos);

    /* start with dist = (i - 1) * v */
    mpz_set_ui(dist, (unsigned long)(omp_get_thread_num() + 2) >> 1);
    mpz_mul(dist, dist, v);

    if (type == KANGAROO_TAME)
    {
        /* start with pos = g^((a + b) / 2 + dist */
        mpz_add(pos, a, b);
        mpz_div_ui(pos, pos, 2);

        mpz_add(pos, pos, dist);
        mpz_powm(pos, g, pos, p);
    }
    else
    {
        /* start with h * g ^dist */
        mpz_powm(pos, g, dist, p);
        mpz_mul(pos, h, pos);
        mpz_mod(pos, pos, p);
    }

    mpz_init(step);
    for (mpz_set_ui(step, 0); mpz_cmp(step, order_g) < 0; mpz_add_ui(step, step, 1))
    {
        if (finish)
            break;

            str = mpz_get_str(NULL, 2, pos);
            index = (int)(hash(str, (size_t)strlen(str)) % r);

            FREE(str);

            mpz_mul(pos, pos, jumps[index]);
            mpz_mod(pos, pos, p);

            mpz_add(dist, dist, dists[index]);

        if (mpz_sizeinbase(pos, 2) < POLLARD_TRESHOLD)
        {
#pragma omp critical
            {
                if (!finish)
                {
                    triple = pollard_triple_create(type, dist, pos);
                    if (darray_get_num_entries(set) > 0 && darray_search_first(set, (void *)&triple, (void *)&pt_p) != -1 && pt_p->type != triple->type)
                    {
                        /* x = (a + b) / 2 + dTAME - dWILD */
                        mpz_init(x);
                        mpz_add(x, a, b);
                        mpz_div_ui(x, x, 2);

                        if (triple->type == KANGAROO_TAME)
                        {
                            mpz_add(x, x, triple->dist);
                            mpz_sub(x, x, pt_p->dist);
                        }
                        else
                        {
                            mpz_add(x, x, pt_p->dist);
                            mpz_sub(x, x, triple->dist);
                        }

                        mpz_set(res, x);
                        mpz_clear(x);

                        finish = true;
                    }
                    else
                        darray_insert(set, (void *)&triple);
                }
            }
        }
    }

    mpz_clear(dist);
    mpz_clear(pos);
    mpz_clear(step);
}
    mpz_mod(res, res, order_g);

    /* cleanup */
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(order_g);
    mpz_clear(beta);
    mpz_clear(v);

    for (i = 0; i < r; ++i)
    {
        mpz_clear(dists[i]);
        mpz_clear(jumps[i]);
    }

    FREE(dists);
    FREE(jumps);

    darray_destroy_with_entries(set);

    return 0;
}
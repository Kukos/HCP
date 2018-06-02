#include <ecm.h>
#include <gmp.h>
#include <log.h>
#include <stdint.h>
#include <stdlib.h>
#include <common.h>

/* 3D point */
typedef struct Point
{
    mpz_t x;
    mpz_t y;
    mpz_t z;
} Point;

/* y^2 = x^3 + ax + b in Zn */
typedef struct ECurve
{
    mpz_t a;
    mpz_t b;
    mpz_t n;
} ECurve;

/*
    Create new 3D Point

    PARAMS
    @IN n - modular

    RETURN
    NULL iff failure
    Pointer to new Point iff success
*/
static Point *point_create(void);

/*
    Destroy 3D point

    PARAMS
    @IN p - pointer to 3D point

    RETURN
    This is a void function
*/
static void point_destroy(Point *p);

/*
    Create Curve from Point in Zn

    PARAMS
    @IN p - point
    @IN n - modular

    RETURN
    NULL iff failure
    Pointer to ECurve iff success
*/
static ECurve *ecurve_create(const Point *p, const mpz_t n);

/*
    Destroy Curve

    PARAMS
    @IN ecurve - pointer to ecurve

    RETURN
    This is a void function
*/
static void ecurve_destroy(ECurve *ecurve);

/*
    Eliptic "scalar" mul, k * 3D Point p in Curve ecurve

    PARAMS
    @IN k - mult
    @IN / OUT p - point
    @IN ecurve

    RETURN
    This is a void function
*/
static void eliptic_mul(uint32_t k, Point *p, const ECurve *ecurve);

/*
    Addition in eliptic curve: inout = inout + in

    PARAMS
    @IN inout - first point and point where result will be stored
    @IN in - second point

*/
static void eliptic_add(Point *inout, const Point *in, const ECurve *ecurve);

static Point *point_create(void)
{
    Point *p;

    TRACE();

    p = (Point *)malloc(sizeof(Point));
    if (p == NULL)
        ERROR("malloc error\n", NULL);

    mpz_init(p->x);
    mpz_init(p->y);
    mpz_init(p->z);

    mpz_set_ui(p->x, 0);
    mpz_set_ui(p->y, 1);
    mpz_set_ui(p->z, 0);

    return p;
}

static Point *point_create_random(const mpz_t n)
{
    Point *p;
    gmp_randstate_t state;

    TRACE();

    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)rand());

    p = (Point *)malloc(sizeof(Point));
    if (p == NULL)
        ERROR("malloc error\n", NULL);

    mpz_init(p->x);
    mpz_init(p->y);
    mpz_init(p->z);

    mpz_urandomm(p->x, state, n);
    mpz_urandomm(p->y, state, n);
    mpz_set_ui(p->z, 1);

    gmp_randclear(state);

    return p;
}

static void point_destroy(Point *p)
{
    TRACE();

    if (p == NULL)
        return;

    mpz_clear(p->x);
    mpz_clear(p->y);
    mpz_clear(p->z);

    FREE(p);
}

static ECurve *ecurve_create(const Point *p, const mpz_t n)
{
    ECurve *ecurve;
    gmp_randstate_t state;

    mpz_t x3;
    mpz_t y2;
    mpz_t ax;

    TRACE();

    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)rand());

    ecurve = (ECurve *)malloc(sizeof(ECurve));
    if (ecurve == NULL)
        ERROR("malloc error\n", NULL);

    mpz_init(ecurve->a);
    mpz_init(ecurve->b);
    mpz_init(ecurve->n);

    mpz_init(x3);
    mpz_init(y2);
    mpz_init(ax);

    mpz_set(ecurve->n, n);
    mpz_urandomm(ecurve->a, state, n);

    /* y^2 = x^3 + ax + b --> b = y^2 - x^3 - ax */
    mpz_pow_ui(x3, p->x, 3);
    mpz_pow_ui(y2, p->y, 2);
    mpz_mul(ax, p->x, ecurve->a);
    mpz_sub(ecurve->b, y2, x3);
    mpz_sub(ecurve->b, ecurve->b, ax);
    mpz_mod(ecurve->b, ecurve->b, n);

    mpz_clear(x3);
    mpz_clear(y2);
    mpz_clear(ax);

    gmp_randclear(state);
    return ecurve;
}

static void ecurve_destroy(ECurve *ecurve)
{
    TRACE();

    if (ecurve == NULL)
        return;

    mpz_clear(ecurve->a);
    mpz_clear(ecurve->b);
    mpz_clear(ecurve->n);

    FREE(ecurve);
}

static void eliptic_add(Point *inout, const Point *in, const ECurve *ecurve)
{
    mpz_t temp;
    mpz_t temp1;
    mpz_t temp2;

    mpz_t inv;

    TRACE();

    /* if inout-> z == 0 then inout = in */
    if (mpz_cmp_ui(inout->z, 0) == 0)
    {
        mpz_set(inout->x, in->x);
        mpz_set(inout->y, in->y);
        mpz_set(inout->z, in->z);

        return;
    }

    /* if in->z == 0 then inout = inout */
    if (mpz_cmp_ui(in->z, 0) == 0)
        return;

    /* if x == x */
    if (mpz_cmp(inout->x, in->x) == 0)
    {
        /* if y + y == 0 then inout = infinity [0,1,0] */
        mpz_init(temp);
        mpz_add(temp, inout->y, in->y);
        mpz_mod(temp, temp, ecurve->n);

        if (mpz_cmp_ui(temp, 0) == 0)
        {
            mpz_set_ui(inout->x, 0);
            mpz_set_ui(inout->y, 1);
            mpz_set_ui(inout->z, 0);

            mpz_clear(temp);

            return;
        }

        /* temp1 = (3x^2 + a) mod n */
        mpz_init(temp1);
        mpz_mul(temp1, inout->x, inout->x);
        mpz_mul_ui(temp1, temp1, 3);
        mpz_add(temp1, temp1, ecurve->a);
        mpz_mod(temp1, temp1, ecurve->n);

        /* temp2 = 2y mod n */
        mpz_init(temp2);
        mpz_mul_ui(temp2, inout->y, 2);
        mpz_mod(temp2, temp2, ecurve->n);
    }
    else
    {
        /* temp1 = Y2 - Y1 mod n */
        mpz_init(temp1);
        mpz_sub(temp1, in->y, inout->y);
        mpz_mod(temp1, temp1, ecurve->n);


        /* temp2 = X2 - X1 mod n */
        mpz_init(temp2);
        mpz_sub(temp2, in->x, inout->x);
        mpz_mod(temp2, temp2, ecurve->n);
    }

    mpz_init(inv);
    /* if cannot invert, inout = non trivial factor [0,0, temp2] */
    if (mpz_invert(inv, temp2, ecurve->n) == 0)
    {
        mpz_set_ui(inout->x, 0);
        mpz_set_ui(inout->y, 0);
        mpz_set(inout->z, temp2);

        mpz_clear(temp1);
        mpz_clear(temp2);
        mpz_clear(inv);

        return;
    }

    /* temp = temp1^2 * inv^2 - X1 - X2 mod n */
    mpz_init(temp);
    mpz_mul(temp, temp1, temp1);
    mpz_mul(temp, temp, inv);
    mpz_mul(temp, temp, inv);
    mpz_sub(temp, temp, inout->x);
    mpz_sub(temp, temp, in->x);
    mpz_mod(temp, temp, ecurve->n);

    /* temp2 = temp1 * inv * (X1 - temp) - Y2 mod n */
    mpz_sub(temp2, inout->x, temp);
    mpz_mul(temp2, temp2, inv);
    mpz_mul(temp2, temp2, temp1);
    mpz_sub(temp2, temp2, inout->y);
    mpz_mod(temp2, temp2, ecurve->n);

    mpz_set(inout->x, temp);
    mpz_set(inout->y, temp2);
    mpz_set_ui(inout->z, 1);

    mpz_clear(temp);
    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(inv);
}


static void eliptic_mul(uint32_t k, Point *p, const ECurve *ecurve)
{
    Point *r;

    TRACE();

    r = point_create_random(ecurve->n);
    r = point_create();

    /* standart fast mult algorithm (k * p) */
    while (k > 0)
    {
        if (mpz_cmp_ui(p->z, 1) > 0)
        {
            point_destroy(r);
            return;
        }

        if (ODD(k))
            eliptic_add(r, p, ecurve); /* r = r + p */

        eliptic_add(p, r, ecurve); /* p = p + r */
        k >>= 1;
    }

    point_destroy(r);
}

int lenstra_ecm(const mpz_t n, Darray *primes, uint32_t limit, mpz_t factor)
{
    uint32_t prime;
    mpz_t p;
    Point *point;
    ECurve *ecurve;

    TRACE();

    point = point_create_random(n);
    if (point == NULL)
        ERROR("point_create error\n", 1);

    ecurve = ecurve_create(point, n);
    if (ecurve == NULL)
        ERROR("ecurve_create error\n", 1);

    mpz_init(p);
    for_each_data(primes, Darray, prime)
        /* p = prime; p < limit; p *= prime */
        for (mpz_set_ui(p, (unsigned long)prime); mpz_cmp_ui(p, (unsigned long)limit) < 0; mpz_mul_ui(p, p, (unsigned long)prime))
        {
            eliptic_mul(prime, point, ecurve);
            if (mpz_cmp_ui(point->z, 1) > 0) /* we have non trivial factor of n */
            {
                /* factor is gcd(n, z) */
                mpz_gcd(factor, n, point->z);

                point_destroy(point);
                ecurve_destroy(ecurve);

                return 0;
            }
        }

    point_destroy(point);
    ecurve_destroy(ecurve);

    return 1;
}
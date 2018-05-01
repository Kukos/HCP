#include <log.h>
#include <pollard.h>
#include <pohling.h>
#include <crt.h>
#include <stdlib.h>
#include <common.h>
#include <string.h>

/*
    Solve pohling subproblem

    PARAMS
    @IN g - generator
    @IN h - result
    @IN p - prime
    @IN f- factor
    @IN e - exponent
    @OUT x - x

    RETURN
    0 iff success
    Non-zero value iff failure
*/
static int solve_discrete_subproblem(const mpz_t g, const mpz_t h, const mpz_t p, const mpz_t f, const mpz_t e, mpz_t x);

/*
    Delete e == zeros in array

    PARAMS
    @IN f - factors
    @IN e - exponents
    @IN len - len of arrays

    RETURN
    LEN of new array
*/
static size_t delete_zeros(mpz_t *f, mpz_t *e, size_t len);

/*
    Calculate Order of subgroup generate by g

    PARAMS
    @IN g - generator
    @IN p - prime
    @IN factors - ord(p) factors
    @IN exponents - ord(p) exponents
    @IN / OUT ord - ord(p), after ord(g)

    RETURN
    This is a void function
*/
static void calculate_ord(const mpz_t g, const mpz_t p, mpz_t *factors, mpz_t *exponents, size_t len, mpz_t ord);

static void calculate_ord(const mpz_t g, const mpz_t p, mpz_t *factors, mpz_t *exponents, size_t len, mpz_t ord)
{

    size_t i;
    size_t j;
    size_t max;
    mpz_t temp;

    TRACE();

    mpz_init(temp);

    /* for each factors */
    for (i = 0; i < len; ++i)
    {
        /* max for each exponent */
        max = (size_t)mpz_get_ui(exponents[i]);
        for (j = 0; j < max; ++j)
        {
            /* check if g^(ord / factors[i]) == 1 */
            mpz_div(temp, ord, factors[i]);
            mpz_powm(temp, g, temp, p);
            if (mpz_cmp_ui(temp, 1) == 0) /* if yes we know that ord is not order */
            {
                mpz_sub_ui(exponents[i], exponents[i], 1);
                mpz_div(ord, ord, factors[i]);
            }
            else
                break;
        }
    }

    if (mpz_cmp_ui(ord, 1) == 0)
    {
        mpz_set(ord, factors[len - 1]);
        mpz_set_ui(exponents[len - 1], 1);
    }

    mpz_clear(temp);
}

static size_t delete_zeros(mpz_t *f, mpz_t *e, size_t len)
{
    size_t i;

    TRACE();

    for (i = 0; i < len; ++i)
        while (mpz_cmp_ui(e[i], 0) == 0)
        {
            (void)memmove(f + i, f + i + 1, (len - i - 1) * sizeof(mpz_t));
            (void)memmove(e + i, e + i + 1, (len - i - 1) * sizeof(mpz_t));
            --len;
        }

    return len;
}

static int solve_discrete_subproblem(const mpz_t g, const mpz_t h, const mpz_t p, const mpz_t f, const mpz_t e, mpz_t x)
{
    mpz_t inv;
    mpz_t less_e;

    mpz_t new_g;

    mpz_t temp1;
    mpz_t temp2;
    mpz_t temp_x;

    mpz_t temp_mod;

    mpz_t i;

    TRACE();

    mpz_init(inv);
    mpz_invert(inv, g, p);

    mpz_init(less_e);
    mpz_sub_ui(less_e, e, 1);

    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(temp_x);
    mpz_init(temp_mod);
    mpz_init(new_g);

    mpz_init(i);
    mpz_set_ui(x, 0);

    /* new_g = g^(f^(e-1)) mod p */
    mpz_powm(new_g, f, less_e, p);
    mpz_powm(new_g, g, new_g, p);

    /* x = x[0] * q^0 + x[1] * q^1 ... + x[e - 1] ^g(e - 1) */
    for (mpz_set_ui(i, 1); mpz_cmp(i, e) <= 0; mpz_add_ui(i, i, 1))
    {
        gmp_printf("\tSUB = %Zd / %Zd\n", i, e);
        mpz_powm(temp_mod, f, i, p);
        /* temp1 = (h *(g^-x))^f^(e - i) mod p */
        mpz_powm(temp1, inv, x, p);
        mpz_mul(temp1, temp1, h);

        mpz_sub(temp2, e, i);
        mpz_powm(temp2, f, temp2, p);

        mpz_powm(temp1, temp1, temp2, p);

        gmp_printf("%Zd ^x = %Zd mod %Zd\n", new_g, temp1, p);
        if (pollard_lambda_parallel_dicsrete_log(new_g, temp1, p, temp_x))
            ERROR("pollard error\n", 1);

        /* x must be in [0 ... f^e] */
        mpz_mod(temp_x, temp_x, temp_mod);

        /* X = x * f ^ (i - 1) */
        mpz_sub_ui(temp2, i, 1);
        mpz_powm(temp1, f, temp2, p);
        mpz_mul(temp_x, temp_x, temp1);
        mpz_add(x, x, temp_x);
    }

    mpz_clear(i);
    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(temp_x);
    mpz_clear(inv);
    mpz_clear(less_e);
    mpz_clear(new_g);
    mpz_clear(temp_mod);

    return 0;
}

int pohling_discrete_log(mpz_t g, mpz_t h, const mpz_t p, mpz_t *factors, mpz_t *exponents, size_t f_len, mpz_t x)
{
    mpz_t *x_array;
    mpz_t *p_array;

    mpz_t ord_p;
    mpz_t new_g;
    mpz_t new_h;

    mpz_t temp1;
    mpz_t temp2;

    mpz_t Q;

    size_t i;

    TRACE();

    x_array = (mpz_t *)malloc(sizeof(mpz_t) * f_len);
    if (x_array == NULL)
        ERROR("malloc error\n", 1);

    p_array = (mpz_t *)malloc(sizeof(mpz_t) * f_len);
    if (p_array == NULL)
        ERROR("malloc error\n", 1);

    /* ord_p = p - 1 */
    mpz_init(ord_p);
    mpz_sub_ui(ord_p, p, 1);

    calculate_ord(g, p, factors, exponents, f_len, ord_p);
    f_len = delete_zeros(factors, exponents, f_len);
    (void)gmp_printf("Order g = %Zd\n", ord_p);

    mpz_init(Q);
    mpz_powm(Q, factors[f_len - 1], exponents[f_len - 1], p);
    if (mpz_cmp(Q, ord_p))
    {
        mpz_powm(g, g, Q, p);
        mpz_powm(h, h, Q, p);
    }

    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(new_g);
    mpz_init(new_h);
    for (i = 0; i < f_len; ++i)
    {
        printf("MAIN = %zu / %zu\n", i + 1, f_len);

        /* new_g = g^(n/f^e), new_h = h^(n/f^e) */
        mpz_powm(temp1, factors[i], exponents[i], p);
        mpz_div(temp2, ord_p, temp1);

        mpz_powm(new_g, g, temp2, p);
        mpz_powm(new_h, h, temp2, p);

        mpz_init(x_array[i]);
        gmp_printf("MAIN %Zd^x = %Zd mod %Zd\n", new_g, new_h, p);
        if (solve_discrete_subproblem(new_g, new_h, p, factors[i], exponents[i], x_array[i]))
            ERROR("Cannot solve discrete log subproblem\n", 1);

        mpz_init(p_array[i]);
        mpz_set(p_array[i], temp1);
    }

    if (crt((const mpz_t *)x_array, (const mpz_t *)p_array, f_len, x))
        ERROR("Chinese remainder error\n", 1);

    mpz_mod(x, x, ord_p);

    mpz_clear(ord_p);
    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(new_g);
    mpz_clear(new_h);
    mpz_clear(Q);

    for (i = 0; i < f_len; ++i)
    {
        mpz_clear(x_array[i]);
        mpz_clear(p_array[i]);
    }

    FREE(x_array);
    FREE(p_array);

    return 0;
}
#include <gmp.h>
#include <log.h>
#include <crt.h>

int crt(const mpz_t *restrict r, const mpz_t *restrict n, size_t len, mpz_t x)
{
    mpz_t product;
    mpz_t temp;
    mpz_t inv;

    size_t i;

    TRACE();

    if (r == NULL)
        ERROR("r == NULL\n", 1);

    if (n == NULL)
        ERROR("n == NULL\n", 1);

    if (len == 0)
        ERROR("len == 0\n", 1);

    mpz_init(product);
    mpz_init(temp);
    mpz_init(inv);

    /* product = n[0] * n[1] * ... * n[len - 1] */
    mpz_set_ui(product, 1);
    for (i = 0; i < len; ++i)
        mpz_mul(product, product, n[i]);

    mpz_set_ui(x, 0);
    for (i = 0; i < len; ++i)
    {
        /* temp = Product / n[i] */
        mpz_div(temp, product, n[i]);

        /* inv = temp^-1 mod n[i] */
        mpz_invert(inv, temp, n[i]);

        /* x += temp * r[i] * inv */
        mpz_mul(temp, temp, inv);
        mpz_mul(temp, temp, r[i]);

        mpz_add(x, x, temp);
    }

    mpz_mod(x, x, product);
    mpz_clear(product);
    mpz_clear(temp);
    mpz_clear(inv);

    return 0;
}
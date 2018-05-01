#ifndef CRT_H
#define CRT_H

/*
    Chinese remainder theorem

    From list of congurenses x = r mod n
                             x = r' mod n'
                             ...

    Finding x

    Author: Michal Kukowski
    email: michalkukowski10@gmail.com

    LICENCE: GPL 3.0
*/

#include <gmp.h>

/*
    Chinese remainder theorem, note that n must be pairwise coprime

    PARAMS
    @IN r - remainders
    @IN n - modulus
    @IN len - len of arrays
    @OUT x - x

    RETURN
    0 iff success
    Non-zero value iff failure
*/
int crt(const mpz_t *restrict r, const mpz_t *restrict n, size_t len, mpz_t x);

#endif
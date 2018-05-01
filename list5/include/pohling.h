#ifndef POHLING_H
#define POHLING_H

/*
    Pohling - Hellman Algorithm to solve g^x = h mod p

    Author: Michal Kukowski
    email: michalkukowski10@gmail.com

    LICENCE: GPL3.0
*/

#include <gmp.h>

/*
    Solve discrete log problem by pohling hellman algo

    PARAMS
    @IN g - generator
    @IN h - result
    @IN p - group prime
    @IN factors - prime factors of p - 1
    @IN exponents - exponennts of factors of p - 1
    @IN f_len - factors and exponents arrays len
    @OUT x - x such that g^x = h mod p

    RETURN
    0 iff success
    Non-zero value iff failure
*/
int pohling_discrete_log(mpz_t g, mpz_t h, const mpz_t p, mpz_t *factors, mpz_t *exponents, size_t f_len, mpz_t x);

#endif
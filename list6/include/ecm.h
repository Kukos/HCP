#ifndef ECM_H
#define ECM_H

/*
    Implementation of Lenstra eliptic curve factorization

    Author: Michal Kukowski
    email: michalkukowski10@gmail.com

    LICENCE: GPL 3.0
*/

#include <gmp.h>
#include <darray.h>
#include <stdint.h>

/*
    Lenstra factorization method

    PARAMS
    @IN n - number to factor
    @IN primes - list of primes < limit
    @IN limit - max iteration
    @OUT factor - first n factor

    RETURN
    0 iff success
    Non-zero value iff failure
*/
int lenstra_ecm(const mpz_t n, Darray *primes, uint32_t limit, mpz_t factor);

#endif
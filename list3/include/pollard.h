#ifndef POLLARD_H
#define POLLARD_H

/*
    Implementation of parallel Pollard rho discrete logarithm algo

    Find x such that
    g^x = h (mod) P, where P is strong prime --> Exist Q such that P = 2Q + 1

    Author: Michal Kukowski
    email: michalkukowski10@gmail.com

    LICENCE: GPL 3.0+
*/

#include <gmp.h>

/*
    Function find X such that g^x = h (mod)p

    PARAMS
    @IN g - generator of Zp
    @IN h - result of power
    @IN p - string prime
    @OUT x - discrete log

    RETURN
    0 iff success
    Non-zero value iff failure
*/
int pollard_rho_parallel_dicsrete_log(const mpz_t g, const mpz_t h, const mpz_t p, mpz_t x);


#endif
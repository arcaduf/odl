/* Copyright 2014-2016 The ODL development group

 This file is part of ODL.

 ODL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ODL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ODL. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stddef.h>
#include <assert.h>
#include <math.h>
#include "var_lp_prox_macros.h"


/*---- Newton iteration ----*/


FLOAT IMPL(startvalue)(FLOAT val, FLOAT p, FLOAT sigma){

    FLOAT sval = val / (p * (2.0 - p) * sigma);

    if(sval <= 1)
        sval = POW(sval, p - 1.0);
    else
        sval = 1.0;

    return sval;
}


FLOAT IMPL(numerator)(FLOAT it, FLOAT p, FLOAT sigma, FLOAT val){
    return p * (p - 2) * sigma * POW(it, p - 1) + val;
}


FLOAT IMPL(denominator)(FLOAT it, FLOAT p, FLOAT sigma){
    return 1.0 + p * (p - 1) * sigma * POW(it, p - 2);
}


FLOAT IMPL(newton_iter)(FLOAT val, FLOAT sigma, FLOAT p, int niter, FLOAT start_relax){

    // Start value which guarantees convergence
    FLOAT it, startval, denom, numer;

    startval = IMPL(startvalue)(val, p, sigma);
    startval *= start_relax;

    // The iteration itself
    for(int i = 0, it = startval; i < niter; i++){
        denom = IMPL(denominator)(it, p, sigma);
        numer = IMPL(numerator)(it, p, sigma, val);
        it = numer / denom;
    }

    return it;
}


/*---- Variable Lp proximal implementation (scalar case) ----*/


FLOAT IMPL(varlp_prox_1)(FLOAT f, FLOAT s){
    return MAX(1.0 - s / ABS(f), 0.0) * f;
}


FLOAT IMPL(varlp_prox_2)(FLOAT f, FLOAT s){
    return f / (1.0 + 2.0 * s);
}


void IMPL(compute_varlp_prox_scalar)(FLOAT *f, FLOAT *p, FLOAT *out, int n, FLOAT s,
                                     int max_newton_iter){
    assert((f != NULL) && (p != NULL) & (out != NULL));

    for(int i = 0; i < n; i++){
        if(f[i] == 0.0)
            out[i] = 0.0;
        else if(p[i] <= 1.05)
            out[i] = IMPL(varlp_prox_1)(f[i], s);
        else if(p[i] >= 1.95)
            out[i] = IMPL(varlp_prox_2)(f[i], s);
        else
            out[i] = IMPL(newton_iter)(f[i], s, p[i], max_newton_iter, 0.5);
    }
}

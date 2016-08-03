/***************************************************************************\
* Copyright (c) 2016, Martin Hansen                                         *
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef CLOVER_TOOLS_H
#define CLOVER_TOOLS_H

#include "suN.h"
#include "spinor_field.h"

#if (defined(WITH_CLOVER) && defined(UPDATE_EO))
#define WITH_CLOVER_EO
#endif

#if (defined(REPR_ADJOINT) || defined(GAUGE_SON))
#define REPR_IS_REAL
#endif

#ifndef _suNfc_multiply
#define _suNfc_multiply(a,b,c) \
	_suNf_multiply(a,b,c)
#endif

#ifndef _suNfc_inverse_multiply
#define _suNfc_inverse_multiply(a,b,c) \
	_suNf_inverse_multiply(a,b,c)
#endif

double get_csw();
void compute_ldl_decomp(double);
void compute_clover_term();
void clover_la_logdet(double, double, scalar_field*);
void compute_force_logdet(double);
void clover_init(double);

#endif

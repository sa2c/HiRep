/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H

void global_sum(double *d, int n);

#include "spinor_field.h"
void complete_gf_sendrecv(suNg_field *gf);
void start_gf_sendrecv(suNg_field *gf);
void complete_sf_sendrecv(spinor_field *gf);
void start_sf_sendrecv(spinor_field *gf);

void test_spinor_field(spinor_field *p);
void complete_gf_sendrecv_flt(suNg_field_flt *gf);



#endif /* COMMUNICATIONS_H */
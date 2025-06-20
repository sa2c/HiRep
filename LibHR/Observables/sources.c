/*************************************************************************** \
 * Copyright (c) 2013 Rudy Arthur, Ari Hietanen                              *
 *                                                                           *
 *                                                                           *
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include "io.h"
#include "random.h"
#include "communications.h"
#include "ranlux.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"
#include "representation.h"

#define PI 3.141592653589793238462643383279502884197

/* Random timeslice not previously chosen */
static int random_tau()
{
  static int *slices = NULL;
  if (slices == NULL)
    slices = (int *)malloc(GLB_T * sizeof(int));
  static int counter = 0;
  int itmp, tau, i;
  double ran;

  if (counter == 0)
  {
    for (i = 0; i < GLB_T; ++i)
    {
      slices[i] = i;
    }
    counter = GLB_T;
  }
  do
  {
    ranlxd(&ran, 1);
    itmp = (int)(ran * counter);
  } while (itmp == counter);
  counter--;
  tau = slices[itmp];
  slices[itmp] = slices[counter];
  slices[counter] = tau;
  bcast_int(&tau, 1);
  return tau;
}

/***************************************************************************\

	Sources:
		point_source:
						source[spin](t,x) = \delta_{a color} \delta_{s, spin} \delta( (t,x) - (tau,0) )
		point_source_loc:
						source[spin](t,x) = \delta_{a color} \delta_{s, spin} \delta( (t,x) - (tau,0) )
		diluted_source_equal_eo:
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \delta(t - tau) \xi(x)
						source[spin](x) = \delta_{s spin} eta(t,x)  -  x even
		diluted_source_equal:
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \delta(t - tau) \xi(x)
						source[spin](t,x) = \delta_{s spin} eta(t,x)  -  x even & odd
		noise_source_equal_eo:
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \xi(t,x)
						source[spin](t,x) = \delta_{s spin} eta(t,x)  -  x even
		gauge_fixed_wall_source:
						source[spin](t,x) = \delta_{a color} \delta_{s spin} \delta(t - tau) 1
		sequential_source:
						source[spin](tf,ti,x) = \gamma_5 S(x,tf; 0,ti)
		gauge_fixed_momentum_source:
						source[spin](t,x) = \delta_{a color} \delta_{s spin} e^{ i p_\mu x_\mu }
		diluted_volume_source:
						source[spin](t,x) = \sum_{x%p == 0} \delta_{s spin} 1_color

                z2_volume_source:                  source[spin](t,x) = Z(2) x Z(2) (no dilution)
    smeared_source:
              Wuppertal smeared source, see Eq. (9) in https://arxiv.org/pdf/1602.05525.pdf.
              The involving gauge field is the original one.
    smeared_source_with_APE:
              Wuppertal smearing with APE smeared gauge field.
\***************************************************************************/
void create_point_source(spinor_field *source, int tau, int color)
{
  int beta, ix;
  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }
  if (COORD[0] == tau / T && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0)
  {
    ix = ipt(tau - zerocoord[0], 0, 0, 0);
    for (beta = 0; beta < 4; ++beta)
    {
      _FIELD_AT(&source[beta], ix)->c[beta].c[color] = 1.;
    }
  }
  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}
// creates point source for the NF color indices.
void create_full_point_source(spinor_field *source, int tau)
{
  int col, beta, idx, ix;

  for (beta = 0; beta < 4 * NF; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  if (zerocoord[0] <= tau && tau < zerocoord[0] + T && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0)
  {
    ix = ipt(tau - zerocoord[0], 0, 0, 0);
    for (col = 0; col < NF; ++col)
    {
      for (beta = 0; beta < 4; ++beta)
      {
        idx = beta + col * 4;
        _FIELD_AT(&source[idx], ix)->c[beta].c[col] = 1.;
      }
    }
  }
  for (beta = 0; beta < 4 * NF; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

void create_point_source_loc(spinor_field *source, int t, int x, int y, int z, int color)
{
  int beta, ix;
  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }
  if (zerocoord[0] <= t && t < zerocoord[0] + T && zerocoord[1] <= x && x < zerocoord[1] + X && zerocoord[2] <= y && y < zerocoord[2] + Y && zerocoord[3] <= z && z < zerocoord[3] + Z)
  {
    ix = ipt(t - zerocoord[0], x - zerocoord[1], y - zerocoord[2], z - zerocoord[3]);
    for (beta = 0; beta < 4; ++beta)
    {
      _FIELD_AT(&source[beta], ix)->c[beta].c[color] = 1.;
    }
  }
  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise
   vectors are equal in each source but placed at a different spin. Even sites only*/
int create_diluted_source_equal_eo(spinor_field *source)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;
  int tau = random_tau();
  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  //  if(COORD[0]==tau/T) {// Check that tau is in this thread.
  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  {
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((tau + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) & 1) == 0)
          {
            v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
            ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
            for (i = 1; i < 4; ++i)
            {
              v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
              *v2 = *v1;
            }
          }
        }
  }
  return tau;
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise
   vectors are equal in each source but placed at a different spin. Even sites only*/
void create_diluted_source_equal_atau_eo(spinor_field *source, int tau)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;
  //int tau = random_tau();
  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }
  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((tau + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) & 1) == 0)
          {
            v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
            ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
            for (i = 1; i < 4; ++i)
            {
              v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
              *v2 = *v1;
            }
          }
        }
  }
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise
   vectors are equal in each source but placed at a different spin. Even and Odd sites*/
int create_diluted_source_equal(spinor_field *source)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;
  int tau = random_tau();
  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
          ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
          for (i = 1; i < 4; ++i)
          {
            v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
            *v2 = *v1;
          }
        }
  }
  return tau;
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise
   vectors are equal in each source but placed at a different spin. Even and Odd sites*/
void create_diluted_source_equal_atau(spinor_field *source, int tau)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;
  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
          ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
          for (i = 1; i < 4; ++i)
          {
            v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
            *v2 = *v1;
          }
        }
  }
}

/* Creates one spinorfield  Z2xZ2 noise sources localised on time slice tau. . Even and Odd sites*/
void create_diluted_source_equal_spinorfield1(spinor_field *source, int tau)
{
  int c[4];
  suNf_vector *v1;
  spinor_field_zero_f(source);

  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          v1 = &((_FIELD_AT(source, ipt(c[0], c[1], c[2], c[3])))->c[0]);
          ranz2((double *)(v1), sizeof(suNf_spinor) / sizeof(double)); // Make new sources
        }
  }
}

/* Creates four Z2xZ2 noise sources NOT localised on time slice but spread over
   all timeslices. The noise vectors are equal in each source but placed at a different spin. Even sites only.*/
void create_noise_source_equal_eo(spinor_field *source)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;

  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((zerocoord[0] + c[0] + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) & 1) == 0)
          {
            v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
            ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
            for (i = 1; i < 4; ++i)
            {
              v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
              *v2 = *v1;
            }
          }
        }
  for (i = 0; i < 4; ++i)
  {
    start_sf_sendrecv(source + i);
    complete_sf_sendrecv(source + i);
  }
}

/* Creates four Z2xZ2 noise sources NOT localised on time slice but spread over
all timeslices. The noise vectors are equal in each source but placed at a different spin, site parity.*/
void create_noise_source_equal_oe(spinor_field *source)
{
  int c[4];
  suNf_vector *v1, *v2;
  int i;

  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((zerocoord[0] + c[0] + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) & 1) == 1)
          {
            v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0]);
            ranz2((double *)(v1), sizeof(suNf_vector) / sizeof(double)); // Make new sources
            for (i = 1; i < 4; ++i)
            {
              v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i]); //Copy previous index.
              *v2 = *v1;
            }
          }
        }
  for (i = 0; i < 4; ++i)
  {
    start_sf_sendrecv(source + i);
    complete_sf_sendrecv(source + i);
  }
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise
   vectors are equal in each source but placed at a different spin. Color dilution and Even and Odd sites */
void create_diluted_source_equal_atau_col(spinor_field *source, int tau, int col)
{
  int c[4];
  // suNf_vector *v1,*v2;
  double complex *v1;
  int i;
  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }
  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          for (i = 0; i < 4; ++i)
          {
            v1 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i].c[col]);
            ranz2((double *)(v1), sizeof(double complex) / sizeof(double)); // Make new sources
          }
        }
  }
}

/* Creates four Z2xZ2 noise sources NOT localised on time slice but spread over
   all timeslices. The noise vectors are equal in each source but placed at a different spin. Even and odd sites. Color dilution */
void create_noise_source_equal_col_dil(spinor_field *source, int col)
{
  int c[4];
  double complex *v1, *v2;
  int i;

  for (i = 0; i < 4; ++i)
  {
    spinor_field_zero_f(&source[i]);
  }

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          v1 = &((_FIELD_AT(&source[0], ipt(c[0], c[1], c[2], c[3])))->c[0].c[col]);
          ranz2((double *)(v1), sizeof(double complex) / sizeof(double)); // Make new sources
          for (i = 1; i < 4; ++i)
          {
            v2 = &((_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])))->c[i].c[col]); //Copy previous index.
            *v2 = *v1;
          }
        }
  for (i = 0; i < 4; ++i)
  {
    start_sf_sendrecv(source + i);
    complete_sf_sendrecv(source + i);
  }
}

//create a wall source at timeslice tau, all parity sites.
void create_gauge_fixed_wall_source(spinor_field *source, int tau, int color)
{
  int c[4];
  int beta;

  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  if (zerocoord[0] <= tau && tau < zerocoord[0] + T)
  { // Check that tau is in this thread.
    c[0] = tau - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          for (beta = 0; beta < 4; ++beta)
          {
            _FIELD_AT(&source[beta], ipt(c[0], c[1], c[2], c[3]))->c[beta].c[color] = 1.;
          }
        }
  }
  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

void create_sequential_source(spinor_field *source, int tf, spinor_field *prop)
{
  int c[4];
  int beta, a, ix;

  suNf_propagator sp0, sp1;

  for (beta = 0; beta < 4 * NF; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  if (zerocoord[0] <= tf && tf < zerocoord[0] + T)
  { // Check that tf is in this thread.
    c[0] = tf - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {

          ix = ipt(c[0], c[1], c[2], c[3]);
          for (a = 0; a < NF; ++a)
          {
            for (beta = 0; beta < 4; beta++)
            {
              _propagator_assign(sp0, *_FIELD_AT(&prop[a * 4 + beta], ix), a, beta);
            }
          }
          _g5_propagator(sp1, sp0);
          _propagator_transpose(sp0, sp1);

          for (a = 0; a < NF; ++a)
          {
            for (beta = 0; beta < 4; beta++)
            {
              *_FIELD_AT(&source[a * 4 + beta], ix) = sp0.c[a].c[beta];
            }
          }
        }
  }
  for (beta = 0; beta < 4 * NF; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

void create_sequential_source_stoch(spinor_field *source, int tf, spinor_field *prop)
{
  int c[4];
  int beta, ix;

  suNf_propagator sp0, sp1;

  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  if (zerocoord[0] <= tf && tf < zerocoord[0] + T)
  { // Check that tf is in this thread.
    c[0] = tf - zerocoord[0];
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {

          ix = ipt(c[0], c[1], c[2], c[3]);

          for (beta = 0; beta < 4; beta++)
          {
            _propagator_assign(sp0, *_FIELD_AT(&prop[beta], ix), 0, beta);
          }

          _g5_propagator(sp1, sp0);
          _propagator_transpose(sp0, sp1);

          for (beta = 0; beta < 4; beta++)
          {
            *_FIELD_AT(&source[beta], ix) = sp0.c[0].c[beta];
          }
        }
  }
  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

//create a e^ipx source
/*void create_gauge_fixed_momentum_source(spinor_field *source, int pt, int px, int py, int pz, int color) {
  int c[4];
  int beta;
  double pdotx;

  for (beta=0;beta<4;++beta){
  spinor_field_zero_f(&source[beta]);
  }

  for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
  pdotx = 2.*PI*( c[0]*pt/GLB_T + c[1]*px/GLB_X + c[2]*py/GLB_Y + c[3]*pz/GLB_Z );
  for (beta=0;beta<4;++beta){
  _FIELD_AT(&source[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color] = cos(pdotx) +I*sin(pdotx);
  }
  }
  for (beta=0;beta<4;++beta){
  start_sf_sendrecv(source + beta);
  complete_sf_sendrecv(source + beta);
  }
  }*/

void create_gauge_fixed_momentum_source(spinor_field *source, int pt, int px, int py, int pz, int color)
{
  int c[4];
  int beta;
  double pdotx;

  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }
  lprintf("Source", 0, "mom = (%d,%d,%d,%d)", pt, px, py, pz);

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          pdotx = 2. * PI * ((double)(c[0] + zerocoord[0]) * (double)pt / (double)GLB_T + (double)(c[1] + zerocoord[1]) * (double)px / (double)GLB_X + (double)(c[2] + zerocoord[2]) * (double)py / (double)GLB_Y + (double)(c[3] + zerocoord[3]) * (double)pz / (double)GLB_Z);
          for (beta = 0; beta < 4; ++beta)
          {
            _FIELD_AT(&source[beta], ipt(c[0], c[1], c[2], c[3]))->c[beta].c[color] = cos(pdotx) + I * sin(pdotx);
          }
        }
  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

// There don't seem to be any noise sources with momentum, so I'll just add it a posteriori using the following function (I hope it's correct!)
void add_momentum(spinor_field *out, spinor_field *in, int px, int py, int pz)
{
  int c[4];
  int beta, beta2, color;
  double pdotx;

  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&out[beta]);
  }
  lprintf("Adding momentum to the source", 0, "mom = (%d,%d,%d)", px, py, pz);

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          pdotx = 2. * PI * ((double)(c[1] + zerocoord[1]) * (double)px / (double)GLB_X + (double)(c[2] + zerocoord[2]) * (double)py / (double)GLB_Y + (double)(c[3] + zerocoord[3]) * (double)pz / (double)GLB_Z);
          for (beta = 0; beta < 4; ++beta)
            for (color = 0; color < NF; ++color)
              for (beta2 = 0; beta2 < 4; ++beta2)
              {
                _FIELD_AT(&out[beta], ipt(c[0], c[1], c[2], c[3]))->c[beta2].c[color] = (_FIELD_AT(&in[beta], ipt(c[0], c[1], c[2], c[3]))->c[beta2].c[color]) * (cos(pdotx) + I * sin(pdotx));
              }
        }

  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(out + beta);
    complete_sf_sendrecv(out + beta);
  }
}

//create a eo source
void create_diluted_volume_source(spinor_field *source, int parity_component, int mod)
{
  int c[4];
  int beta, b;

  for (beta = 0; beta < 4; ++beta)
  {
    spinor_field_zero_f(&source[beta]);
  }

  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((zerocoord[0] + c[0] + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) % mod) == parity_component)
          {
            for (beta = 0; beta < 4; ++beta)
            {
              for (b = 0; b < NF; ++b)
              {
                _FIELD_AT(&source[beta], ipt(c[0], c[1], c[2], c[3]))->c[beta].c[b] = 1.;
              }
            }
          }
        }

  for (beta = 0; beta < 4; ++beta)
  {
    start_sf_sendrecv(source + beta);
    complete_sf_sendrecv(source + beta);
  }
}

void create_z2_volume_source(spinor_field *source)
{
  z2_spinor_field(source);
}

// set to zero even or odd site of the source.
void zero_even_or_odd_site_spinorfield(spinor_field *source, int nspinor, int eo)
{
  int c[4];
  int i;
  for (c[0] = 0; c[0] < T; c[0]++)
    for (c[1] = 0; c[1] < X; c[1]++)
      for (c[2] = 0; c[2] < Y; c[2]++)
        for (c[3] = 0; c[3] < Z; c[3]++)
        {
          if (((zerocoord[0] + c[0] + zerocoord[1] + c[1] + zerocoord[2] + c[2] + zerocoord[3] + c[3]) & 1) == eo)
          {
            for (i = 0; i < nspinor; ++i)
            {
              _spinor_zero_f(*_FIELD_AT(&source[i], ipt(c[0], c[1], c[2], c[3])));
            }
          }
        }
}


/* Updates for Wuppertal smearing by HH in 2020*/

void smearing_function(spinor_field *source, int tau, double epsilon){
  /**
   * @brief Applies the Wuppertal smearing to the source field.
   * The smearing is done in the spatial directions only, basing on Eq. (9) in https://arxiv.org/pdf/1602.05525.pdf:
   * $\Phi q(x) = 1/(1+2d\epsilon)[ q(x) + \epsilon\sum_{\mu=1}^{3} U_{\mu}(x)q(x+\mu) ]$
   * $\Phi$ is the smearing function, $q(x)$ is the source field, $U_{\mu}(x)$ is the link variable in the $\mu$ direction.
   * $d = 3$ for spatial directions.
   * 
   * @param source -> the source field, $q(x)$.
   * @param epsilon -> the smearing parameter, $\epsilon$.
   * @param tau is the time slice where the source is located.
   */
    
    int ix, x, y, z, ix_up, ix_right, ix_front, ix_left, ix_back, ix_down;
    double norm_factor = 1./(1.+6.*epsilon);
    suNf_spinor spinor_OG, spinor_smeared, spinor_tmp, spinor_right, spinor_left, spinor_front, spinor_back, spinor_up, spinor_down;
    _spinor_zero_f(spinor_OG);_spinor_zero_f(spinor_tmp);_spinor_zero_f(spinor_right);_spinor_zero_f(spinor_left);
    _spinor_zero_f(spinor_front);_spinor_zero_f(spinor_back);_spinor_zero_f(spinor_up);_spinor_zero_f(spinor_down);
        
    spinor_field* smeared_source = alloc_spinor_field_f(4,&glattice);
    for (int beta=0;beta<4;++beta){
        spinor_field_zero_f(&smeared_source[beta]);
    }
    
    if(COORD[0]==tau/T){ // apply the smearing only on the time slice where the source is located.
        for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){
            
            // assign the positions: a site (tau, x, y, z) and its neighbors in the spatial directions.
            ix = ipt(tau - zerocoord[0], x,y,z);
            ix_right = iup(ix,1);
            ix_left  = idn(ix,1);
            ix_front = iup(ix,2);
            ix_back  = idn(ix,2);
            ix_up    = iup(ix,3);
            ix_down  = idn(ix,3);
    
    
    
            for (int spin=0;spin<4;spin++){
                spinor_OG    = *_FIELD_AT(&source[spin], ix); // OG: original spinor field at (x,y,z).
                spinor_right = *_FIELD_AT(&source[spin], ix_right); 
                spinor_left  = *_FIELD_AT(&source[spin], ix_left);
                spinor_front = *_FIELD_AT(&source[spin], ix_front);
                spinor_back  = *_FIELD_AT(&source[spin], ix_back);
                spinor_up    = *_FIELD_AT(&source[spin], ix_up);
                spinor_down  = *_FIELD_AT(&source[spin], ix_down);
                        
        
                _spinor_zero_f(spinor_smeared);

                // mutiply the spinor fields with link variable in the spatial directions and add to the smearing field.
                _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,1), spinor_right.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
        
                _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_left,1), spinor_left.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                    
                _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,2), spinor_front.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
        
                _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_back,2), spinor_back.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                    
                _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,3), spinor_up.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
        
                _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_down,3), spinor_down.c[spin]);
                _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                    
                _vector_mul_add_assign_f(spinor_smeared.c[spin], norm_factor, spinor_OG.c[spin]);
        
        _FIELD_AT(&smeared_source[spin], ix)->c[spin]  = spinor_smeared.c[spin];
               
        }
    }
    
    // beta is the spinor index.
    for (int beta=0;beta<4;++beta){
        spinor_field_copy_f(source + beta, smeared_source +beta);
    }
    }
    
    for (int beta=0;beta<4;++beta){
        start_sf_sendrecv(source + beta);
        complete_sf_sendrecv(source + beta);
    }
    free_spinor_field_f(smeared_source);
}

void smearing_function_with_APE(spinor_field *source, int tau, double epsilon){
  /**
   * @brief Applies the Wuppertal smearing to the source field with APE smeared gauge field,
   *        which is similar to the function `smearing_function` but changes the gauge field to the APE smeared one:
   *        `pu_gauge_f` -> `pu_gauge_APE_f`
   */
    
    int ix, x, y, z, ix_up, ix_right, ix_front, ix_left, ix_back, ix_down;
    double norm_factor = 1./(1.+6.*epsilon);
    suNf_spinor spinor_OG, spinor_smeared, spinor_tmp, spinor_right, spinor_left, spinor_front, spinor_back, spinor_up, spinor_down;
    _spinor_zero_f(spinor_OG);_spinor_zero_f(spinor_tmp);_spinor_zero_f(spinor_right);_spinor_zero_f(spinor_left);
    _spinor_zero_f(spinor_front);_spinor_zero_f(spinor_back);_spinor_zero_f(spinor_up);_spinor_zero_f(spinor_down);
    
    spinor_field* smeared_source = alloc_spinor_field_f(4, &glattice);
    for (int beta=0;beta<4;++beta){
        spinor_field_zero_f(&smeared_source[beta]);
    }
    
    represent_gauge_field_APE();
    if(COORD[0]==tau/T){
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
            
        //lprintf("SMEAR",0,"smearing function at(%d,%d,%d,%d)\n",tau,x,y,z);
        ix = ipt(tau - zerocoord[0], x, y, z);
        ix_right = iup(ix,1);
        ix_left  = idn(ix,1);
        ix_front = iup(ix,2);
        ix_back  = idn(ix,2);
        ix_up    = iup(ix,3);
        ix_down  = idn(ix,3);
            
        for (int spin=0;spin<4;spin++){
            spinor_OG    = *_FIELD_AT(&source[spin], ix);
            spinor_right = *_FIELD_AT(&source[spin], ix_right);
            spinor_left  = *_FIELD_AT(&source[spin], ix_left);
            spinor_front = *_FIELD_AT(&source[spin], ix_front);
            spinor_back  = *_FIELD_AT(&source[spin], ix_back);
            spinor_up    = *_FIELD_AT(&source[spin], ix_up);
            spinor_down  = *_FIELD_AT(&source[spin], ix_down);
                
            _spinor_zero_f(spinor_smeared);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,1), spinor_right.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_left,1), spinor_left.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,2), spinor_front.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_back,2), spinor_back.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,3), spinor_up.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_down,3), spinor_down.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _vector_mul_add_assign_f(spinor_smeared.c[spin], norm_factor, spinor_OG.c[spin]);
                
            *_FIELD_AT(&smeared_source[spin], ix)  = spinor_smeared;
        }
    }
    for (int beta=0;beta<4;++beta){
        spinor_field_copy_f(source + beta, smeared_source +beta);
    }
    }
    
    for (int beta=0;beta<4;++beta){
        start_sf_sendrecv(source + beta);
        complete_sf_sendrecv(source + beta);
    }
    free_spinor_field_f(smeared_source);
}

void create_smeared_source(spinor_field *source, int t, int x, int y, int z, int color, double epsilon, int Nsmear){
  /**
   * @brief Creates a smeared source field at a given position (t,x,y,z) with a given color by applying
   *        the Wuppertal smearing function on a point source based on Eq. (9) in https://arxiv.org/pdf/1602.05525.pdf.
   *        
   * @param t The time slice where the source is located.
   * @param x The location of the source.
   * @param y The location of the source.
   * @param z The location of the source.
   * @param Nsmear is the iteration number at the source
   * @param epsilon is the step size, $\epsilon$ in the reference.
   * 
   * @note The gauge field is the original one.
   */
    
    int beta;
    
    for (beta=0;beta<4;++beta){
        spinor_field_zero_f(&source[beta]);
    }
    
    create_point_source_loc(source, t, x, y, z, color);
    
    lprintf("SMEAR",0,"Smeared Source at (%d,%d,%d,%d) with APE smearing \n",t, x,y,z);
    lprintf("SMEAR",0,"source smearing epsilon = %f iterations: \n", epsilon);
    
    for (int n=0;n<Nsmear;n++){
        lprintf("SMEAR",0,"%d...", n+1);
        smearing_function(source, t, epsilon);
    }
}

void create_smeared_source_with_APE(spinor_field *source, int t, int x, int y, int z, int color, double epsilon, int Nsmear){
  /**
   * @brief Creates a smeared source field at a given position (t,x,y,z) with a given color by applying
   *        the Wuppertal smearing function on a point source based on Eq. (9) in https://arxiv.org/pdf/1602.05525.pdf.
   *        
   * @param t The time slice where the source is located.
   * @param x The location of the source.
   * @param y The location of the source.
   * @param z The location of the source.
   * @param Nsmear is the iteration number at the source
   * @param epsilon is the step size, $\epsilon$ in the reference.
   * 
   * @note The gauge field is the APE smeared one.
   */
    int beta;
    
    for (beta=0;beta<4;++beta){
        spinor_field_zero_f(&source[beta]);
    }
    
    create_point_source_loc(source, t, x, y, z, color);
    
    lprintf("SMEAR",0,"Smeared Source at (%d,%d,%d,%d) with APE smearing \n",t, x,y,z);
    lprintf("SMEAR",0,"source smearing epsilon = %f iterations: \n", epsilon);
    
    for (int n=0;n<Nsmear;n++){
        lprintf("SMEAR",0,"%d...", n+1);
        smearing_function_with_APE(source, t, epsilon);
    }
}

// FZ 2024: Add extra smeared sources for disconnected pieces
void smearing_function_volume(spinor_field *source, double epsilon){
    
    int ix, t, x, y, z, ix_up, ix_right, ix_front, ix_left, ix_back, ix_down;
    double norm_factor = 1./(1.+6.*epsilon);
    suNf_spinor spinor_OG, spinor_smeared, spinor_tmp, spinor_right, spinor_left, spinor_front, spinor_back, spinor_up, spinor_down;
    _spinor_zero_f(spinor_OG);_spinor_zero_f(spinor_tmp);_spinor_zero_f(spinor_right);_spinor_zero_f(spinor_left);
    _spinor_zero_f(spinor_front);_spinor_zero_f(spinor_back);_spinor_zero_f(spinor_up);_spinor_zero_f(spinor_down);
        
    spinor_field* smeared_source = alloc_spinor_field_f(4,&glattice);
    for (int beta=0;beta<4;++beta){
        spinor_field_zero_f(&smeared_source[beta]);
    }
    
    for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){
        
        ix = ipt(t, x,y,z);
        ix_right = iup(ix,1);
        ix_left  = idn(ix,1);
        ix_front = iup(ix,2);
        ix_back  = idn(ix,2);
        ix_up    = iup(ix,3);
        ix_down  = idn(ix,3);

        for (int spin=0;spin<4;spin++){
            spinor_OG    = *_FIELD_AT(&source[spin], ix);
            spinor_right = *_FIELD_AT(&source[spin], ix_right);
            spinor_left  = *_FIELD_AT(&source[spin], ix_left);
            spinor_front = *_FIELD_AT(&source[spin], ix_front);
            spinor_back  = *_FIELD_AT(&source[spin], ix_back);
            spinor_up    = *_FIELD_AT(&source[spin], ix_up);
            spinor_down  = *_FIELD_AT(&source[spin], ix_down);
                    
    
            _spinor_zero_f(spinor_smeared);
    
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,1), spinor_right.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
    
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_left,1), spinor_left.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,2), spinor_front.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
    
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_back,2), spinor_back.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix,3), spinor_up.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
    
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_f(ix_down,3), spinor_down.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
                
            _vector_mul_add_assign_f(spinor_smeared.c[spin], norm_factor, spinor_OG.c[spin]);
    
            _FIELD_AT(&smeared_source[spin], ix)->c[spin]  = spinor_smeared.c[spin];       
        }
    }
    
    for (int beta=0;beta<4;++beta){
        spinor_field_copy_f(source + beta, smeared_source +beta);
    }
    
    
    for (int beta=0;beta<4;++beta){
        start_sf_sendrecv(source + beta);
        complete_sf_sendrecv(source + beta);
    }
    free_spinor_field_f(smeared_source);
}

void smearing_function_volume_with_APE(spinor_field *source, double epsilon){
    
    int ix, t, x, y, z, ix_up, ix_right, ix_front, ix_left, ix_back, ix_down;
    double norm_factor = 1./(1.+6.*epsilon);
    suNf_spinor spinor_OG, spinor_smeared, spinor_tmp, spinor_right, spinor_left, spinor_front, spinor_back, spinor_up, spinor_down;
    _spinor_zero_f(spinor_OG);_spinor_zero_f(spinor_tmp);_spinor_zero_f(spinor_right);_spinor_zero_f(spinor_left);
    _spinor_zero_f(spinor_front);_spinor_zero_f(spinor_back);_spinor_zero_f(spinor_up);_spinor_zero_f(spinor_down);
    
    spinor_field* smeared_source = alloc_spinor_field_f(4, &glattice);
    for (int beta=0;beta<4;++beta){
        spinor_field_zero_f(&smeared_source[beta]);
    }
    
    represent_gauge_field_APE();

    for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
            
        ix = ipt(t, x, y, z);
        ix_right = iup(ix,1);
        ix_left  = idn(ix,1);
        ix_front = iup(ix,2);
        ix_back  = idn(ix,2);
        ix_up    = iup(ix,3);
        ix_down  = idn(ix,3);
            
        for (int spin=0;spin<4;spin++){
            spinor_OG    = *_FIELD_AT(&source[spin], ix);
            spinor_right = *_FIELD_AT(&source[spin], ix_right);
            spinor_left  = *_FIELD_AT(&source[spin], ix_left);
            spinor_front = *_FIELD_AT(&source[spin], ix_front);
            spinor_back  = *_FIELD_AT(&source[spin], ix_back);
            spinor_up    = *_FIELD_AT(&source[spin], ix_up);
            spinor_down  = *_FIELD_AT(&source[spin], ix_down);
                
            _spinor_zero_f(spinor_smeared);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,1), spinor_right.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_left,1), spinor_left.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,2), spinor_front.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_back,2), spinor_back.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _suNf_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix,3), spinor_up.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            _suNf_inverse_multiply(spinor_tmp.c[spin], *pu_gauge_APE_f(ix_down,3), spinor_down.c[spin]);
            _vector_mul_add_assign_f(spinor_smeared.c[spin], epsilon*norm_factor, spinor_tmp.c[spin]);
            
            _vector_mul_add_assign_f(spinor_smeared.c[spin], norm_factor, spinor_OG.c[spin]);
                
            *_FIELD_AT(&smeared_source[spin], ix)  = spinor_smeared;
        }
    }
    for (int beta=0;beta<4;++beta){
        spinor_field_copy_f(source + beta, smeared_source +beta);
    }
    
    for (int beta=0;beta<4;++beta){
        start_sf_sendrecv(source + beta);
        complete_sf_sendrecv(source + beta);
    }
    free_spinor_field_f(smeared_source);
}

void create_noise_source_equal_eo_smeared_with_APE(spinor_field *source, double epsilon, int Nsmear){
    
    int beta;
    for (beta=0;beta<4;++beta){
        spinor_field_zero_f(&source[beta]);
    }
    
    create_noise_source_equal_eo(source);
    
    lprintf("SMEAR",0,"Smeared noisy volume source with APE smearing \n");
    lprintf("SMEAR",0,"source smearing epsilon = %f iterations: \n", epsilon);
    
    for (int n=0;n<Nsmear;n++){
        lprintf("SMEAR",0,"%d...", n+1);
        smearing_function_volume_with_APE(source, epsilon);
    }
}
void create_noise_source_equal_eo_smeared(spinor_field *source, double epsilon, int Nsmear){
    
    int beta;
    for (beta=0;beta<4;++beta){
        spinor_field_zero_f(&source[beta]);
    }

    create_noise_source_equal_eo(source);
    
    lprintf("SMEAR",0,"Smeared noisy volume source \n");
    lprintf("SMEAR",0,"source smearing epsilon = %f iterations: \n", epsilon);
    
    for (int n=0;n<Nsmear;n++){
        lprintf("SMEAR",0,"%d...", n+1);
        smearing_function_volume(source, epsilon);
    }
}
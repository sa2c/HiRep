//
//  APE_smearing.c
//
//  Code added by HH on 2020.11.2 for APE smearing,
//  based on the functions written in Fortran provided by Davide Vadacchino.
//  
//

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "global.h"
#include "utils.h"
#include "suN.h"
#include "memory.h"
#include "communications.h"
#include "logger.h"
#include "complex.h"
#include "random.h"
#include "representation.h"
#include "observables.h"
#include "update.h"

void APE_smearing(double smear_val, int Nsmear){
    /**
     * APE smearing function is used to smooth the gauge field, we take the defination
     * from Eq. (16) in https://arxiv.org/abs/hep-lat/0409141v2.
     * 
     * @param smear_val: Smearing parameter $\alpha$ in the reference
     * @param Nsmear: Number of smearing iterations.
     */

    suNg staple, tr1, tr2;
    int ix, iy, iz, it;
    int mu, nu, mid, midpmu, midpnu, midmnu, midpmumnu;
    suNg v, vout, vtmp;
    
    double sm1 = 1. - smear_val;  // ( 1 - \alpha )
    double sm2 = smear_val / 6. ; // ( \alpha \ 6 )
    
    suNg_field *u_gauge_tmp = alloc_gfield(&glattice);
    u_gauge_APE = alloc_gfield(&glattice);
    
    #ifdef ALLOCATE_REPR_GAUGE_FIELD
    u_gauge_APE_f=alloc_gfield_f(&glattice);
    #endif
    
    for (it=0;it<T; it++){
        for (ix = 0; ix < X; ix++) for (iy = 0; iy < Y; iy++) for (iz = 0; iz < Z; iz++){
            for (mu = 0; mu < 4; mu++){
                mid = ipt(it,ix,iy,iz);
                *pu_gauge_APE(mid, mu) = *pu_gauge(mid, mu);
            }
        }
    }
    
    start_gf_sendrecv(u_gauge_APE);
    complete_gf_sendrecv(u_gauge_APE);
    represent_gauge_field_APE();
    
    lprintf("APE",0,"APE smearing with val=%f \n", smear_val);
    lprintf("APE",0,"N=0 <p>_s = %1.6f\n",avr_spatial_plaquette_APE());
    
    for (int N=0; N<Nsmear; N++){
        
        for (it=0;it<T; it++){
            for (ix = 0; ix < X; ix++) for (iy = 0; iy < Y; iy++) for (iz = 0; iz < Z; iz++){
                mid = ipt(it, ix, iy, iz);
                for (mu = 1; mu < 4; mu++){
                    _suNg_zero(v);
                    midpmu = iup(mid, mu);
                
                    for (int i = 0; i < 2; i++){
                        nu = (mu + i) % 3 + 1;
                        midpnu = iup(mid, nu);
                        midmnu = idn(mid, nu);
                        midpmumnu = idn(midpmu, nu);
                    
                        //Up Staple
                        _suNg_times_suNg(tr1, *pu_gauge_APE(mid, nu), *pu_gauge_APE(midpnu, mu));
                        _suNg_times_suNg_dagger(staple, tr1, *pu_gauge_APE(midpmu, nu));
                        
                        _suNg_add_assign(v, staple);

                        //Down Staple
                        _suNg_times_suNg(tr2, *pu_gauge_APE(midmnu, mu), *pu_gauge_APE(midpmumnu, nu));
                        _suNg_dagger_times_suNg(staple, *pu_gauge_APE(midmnu, nu), tr2);

                        _suNg_add_assign(v, staple);
                        }
                
                    _suNg_mul(tr1, sm1, *pu_gauge_APE(mid, mu));
                    _suNg_mul(tr2, sm2, v);
                    
                    _suNg_zero(vout);
                    _suNg_add_assign(vout,tr1);
                    _suNg_add_assign(vout,tr2);
                    
                    _suNg_mul(vtmp, 1., vout);
                    project_to_suNg(&vout); // Project to the group manifold
                    
                    #ifdef GAUGE_SPN
                    cooling_SPN(&v, &vout, &vtmp, 5); // max Re Tr(U V^+)
                    *pu_gauge_tmp(mid, mu) = v;
                    
                    #else
                    
                    *pu_gauge_tmp(mid, mu) = vout;
                    
                    #endif
                }
            }
        }
        
        for (it=0;it<T; it++){
            for (ix = 0; ix < X; ix++) for (iy = 0; iy < Y; iy++) for (iz = 0; iz < Z; iz++){
                for (mu = 0; mu < 4; mu++){
                    mid = ipt(it,ix,iy,iz);
                    *pu_gauge_APE(mid, mu) = *pu_gauge_tmp(mid, mu);
                }
            }
        }
        
        start_gf_sendrecv(u_gauge_APE);
        complete_gf_sendrecv(u_gauge_APE);
        represent_gauge_field_APE();
        
    }
    lprintf("APE",0,"N=%d <p>_s = %1.6f\n", Nsmear, avr_spatial_plaquette_APE());
    polyakov_APE();
    full_plaquette_APE();
    lprintf("APE",0,"APE smearing END \n");
    free_gfield(u_gauge_tmp);
}
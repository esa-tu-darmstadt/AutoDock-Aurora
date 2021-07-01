#ifndef _INCLUDE_LOCALSEARCH_ADADELTA_H_
#define _INCLUDE_LOCALSEARCH_ADADELTA_H_

#include "rand_gen.h"

/*
 * Ported from AutoDock-GPU v1.2
 * https://github.com/ccsb-scripps/AutoDock-GPU/blob/eed190fd08844f5a1336cf36e267d7be3d8a6b61/device/kernel_ad.cl
 */

void ls_ad(
    uint DockConst_num_of_atoms,    // FIXME: SW defines it as uchar
    
    
    float DockConst_gridsize_x_minus1, // FIXME: ADGPU defines it as int
    float DockConst_gridsize_y_minus1,
    float DockConst_gridsize_z_minus1,
    //uint DockConst_gridsize_x_times_y, // FIXME: present in ADPGU, but not in SW
    //uint DockConst_gridsize_x_times_y_times_z, // FIXME: present in ADGPU, but not in SW
    
    // ie
    const float* restrict IE_Fgrids,

    // pc
    uint DockConst_rotbondlist_length,

    // ia
    const   float*  restrict    IA_IE_atom_charges,
    const   int*    restrict    IA_IA_atom_types,
    const   int*    restrict    IA_intraE_contributors,
    const   float*  restrict    IA_reqm,
    const   float*  restrict    IA_reqm_hbond,
    const   uint*   restrict    IA_atom1_types_reqm,
    const   uint*   restrict    IA_atom2_types_reqm,
    const   float*  restrict    IA_VWpars_AC,
    const   float*  restrict    IA_VWpars_BD,
    const   float*  restrict    IA_dspars_S,
    const   float*  restrict    IA_dspars_V,
            float               DockConst_smooth,
            uint                DockConst_num_of_intraE_contributors,
            float               DockConst_grd_spacing,
            uint                DockConst_num_of_atypes,
            float               DockConst_coeff_elec,
            float               DockConst_qasp,
            float               DockConst_coeff_desolv,
    

    float in_out_genotype[][MAX_POPSIZE], // ADGPU: __global float* restrict dockpars_conformations_next,
    float* restrict in_out_energy, // ADGPU: __global float* restrict dockpars_energies_next,
    uint* restrict out_eval, // ADGPU: __global int* restrict dockpars_evals_of_new_entities,

    uint pop_size,
    uchar DockConst_num_of_genes, // ADGPU defines it as int




);

#endif
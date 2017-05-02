#include "energy_force.h"
#include "params.h"
#include "atoms.h"
#include "timer.h"

//************************************************************************
// compute_long_range_correction() function
//   - Calculates long range correction due to finite interaction cutoff.
//   - Arguments:
//       - len_jo: struct containing leonard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//       - energy_long: long-range correction to energy.
//       - force_long: long-range correction to force.
//************************************************************************
void compute_long_range_correction(const lj_params *len_jo,
                                   const misc_params *m_pars,
                                   float *energy_long, float *force_long) {

  float ulongpre =
      m_pars->float_N * 8.0 * len_jo->eps * m_pars->pi * m_pars->density;
  *energy_long = ulongpre * (len_jo->sig12 / (9.0 * len_jo->rcut9) -
                             len_jo->sig6 / (6.0 * len_jo->rcut3));

  float vlongpre = 96.0 * len_jo->eps * m_pars->pi * m_pars->density;
  *force_long = -1.0 * vlongpre * (len_jo->sig12 / (9.0 * len_jo->rcut9) -
                                   len_jo->sig6 / (6.0 * len_jo->rcut3));
}

//************************************************************************
// compute_energy_and_force() function
//   - Calculates energy and force acting on each atom.
//   - Arguments:
//       - myatoms: struct containing all atomic information.
//       - len_jo: struct containing lennard jones interaction parameters.
//       - m_pars: struct containing misc. parameters.
//************************************************************************
void compute_energy_and_force(Atoms *myatoms, const lj_params *len_jo,
                              const misc_params *m_pars) {

  timeit(1, 0);
  int atomi, atomj;
  float *restrict fx, *restrict fy, *restrict fz;
  for (atomi = 0; atomi < myatoms->N; ++atomi) {
    fx = myatoms->fx + atomi;
    fy = myatoms->fy + atomi;
    fz = myatoms->fz + atomi;
    *fx = 0.0;
    *fy = 0.0;
    *fz = 0.0;
  }
  myatoms->pot_energy = 0.0;
  myatoms->virial = 0.0;

  /* we know these don't alias due to the restriction below (atomj > atomi) */
  float *restrict xx_i, *restrict xx_j, *restrict yy_i, *restrict yy_j,
      *restrict zz_i, *restrict zz_j, *restrict fx_i, *restrict fx_j,
      *restrict fy_i, *restrict fy_j, *restrict fz_i,
      *restrict fz_j, *restrict pot_ptr = &myatoms->pot_energy,
                                *restrict virial_ptr = &myatoms->virial;
  float xxs[myatoms->N];
  float yys[myatoms->N];
  float zzs[myatoms->N];
  float pots[myatoms->N];
  float virials[myatoms->N];
#pragma omp parallel for schedule(                                             \
    guided), private(xxs, yys, zzs, pots, virials, xx_i, xx_j, yy_i, yy_j,     \
                     zz_i, zz_j, atomj, fx_i, fx_j, fy_i, fy_j, fz_i, fz_j)
  for (atomi = 0; atomi < myatoms->N; ++atomi) {
    xx_i = myatoms->xx + atomi;
    yy_i = myatoms->yy + atomi;
    zz_i = myatoms->zz + atomi;
    fx_i = myatoms->fx + atomi;
    fy_i = myatoms->fy + atomi;
    fz_i = myatoms->fz + atomi;
    memset(xxs, 0, sizeof(xxs));
    memset(yys, 0, sizeof(yys));
    memset(zzs, 0, sizeof(zzs));
    memset(pots, 0, sizeof(pots));
    memset(virials, 0, sizeof(virials));
    for (atomj = atomi + 1; atomj < myatoms->N; ++atomj) {
      xx_j = myatoms->xx + atomj;
      yy_j = myatoms->yy + atomj;
      zz_j = myatoms->zz + atomj;
      fx_j = myatoms->fx + atomj;
      fy_j = myatoms->fy + atomj;
      fz_j = myatoms->fz + atomj;
      float xxi = minimum_image(*xx_i - *xx_j, m_pars->side, m_pars->sideh);
      float yyi = minimum_image(*yy_i - *yy_j, m_pars->side, m_pars->sideh);
      float zzi = minimum_image(*zz_i - *zz_j, m_pars->side, m_pars->sideh);
      float dis2 = xxi * xxi + yyi * yyi + zzi * zzi;
      if (dis2 <= len_jo->rcut2) {
        float dis2i = 1.0 / dis2;
        float dis6i = dis2i * dis2i * dis2i;
        float dis12i = dis6i * dis6i;
        pots[atomj] = len_jo->sig12 * dis12i - len_jo->sig6 * dis6i;
        float fterm =
            dis2i * (2.0 * len_jo->sig12 * dis12i - len_jo->sig6 * dis6i);
        virials[atomj] = fterm * dis2;

        xxs[atomj] = fterm * xxi;
        yys[atomj] = fterm * yyi;
        zzs[atomj] = fterm * zzi;
      }
    }
#pragma omp critical
    for (atomj = atomi + 1; atomj < myatoms->N; ++atomj) {
      fx_j = myatoms->fx + atomj;
      fy_j = myatoms->fy + atomj;
      fz_j = myatoms->fz + atomj;
      *fx_i += xxs[atomj];
      *fy_i += yys[atomj];
      *fz_i += zzs[atomj];
      *fx_j -= xxs[atomj];
      *fy_j -= yys[atomj];
      *fz_j -= zzs[atomj];
      *pot_ptr += pots[atomj];
      *virial_ptr -= virials[atomj];
    }
  }
  for (atomi = 0; atomi < myatoms->N; ++atomi) {
    fx = myatoms->fx + atomi;
    fy = myatoms->fy + atomi;
    fz = myatoms->fz + atomi;
    *fx *= 24.0 * len_jo->eps;
    *fy *= 24.0 * len_jo->eps;
    *fz *= 24.0 * len_jo->eps;
  }
  myatoms->pot_energy *= 4.0 * len_jo->eps;
  myatoms->virial *= 24.0 * len_jo->eps;
  timeit(1, 1);
}

//**********************************************************************
// minimum_image() function
//   - Finds the nearest images of atom i and j, and returns distance.
//   - Arguments:
//       - dist: 1d distance between atoms i and j in central sim. cell.
//       - box_length: length of simulation cell.
//       - half_box_length: half of length of simulation cell.
//**********************************************************************
float minimum_image(const float dist, const float box_length,
                    const float half_box_length) {

  float min_dist = dist;
  if (dist > half_box_length)
    min_dist = dist - box_length;
  if (dist < -half_box_length)
    min_dist = dist + box_length;
  return min_dist;
}

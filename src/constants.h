#pragma once

/* Particle constants */

#define P_R  (1.0f)  // radius
#define P_H  (2  * P_R)
#define P_H2 (P_H  * P_H)
#define P_H3 (P_H  * P_H2)
#define P_H6 (P_H3 * P_H3)
#define P_H9 (P_H3 * P_H6)

#define MAX_NEIGHBORS 16


/* Math constants */

#define PI glm::pi<float>()
#define EPS 1.0E-5f  // numerical zero


/* Simulation constants */
#define SOLVER_ITERATIONS 3
#define XYBOUND 16.0f
#define ZBOUND  64.0f

#define INIT_STEP 0.64f
#define INIT_KLEVELS 4

#define PARCEL_R 3.2f
#define PARCEL_MIN_XY (PARCEL_R)
#define PARCEL_MAX_XY (XYBOUND - PARCEL_R)
#define PARCEL_Z (ZBOUND / 2)
#define PARCEL_STEP 1.0f


/* Physics constants */
#define DT 0.02f
#define G  9.81f

#define POLY6_COEFF (315.0f / (64.0f * PI * P_H9))
#define SPIKY_COEFF (45.0f  / (PI * P_H6))

#define RHO_0 1.0f
#define RHO_0_INV (1 / RHO_0)
#define RELAXATION_EPS 5.0f

#define DELTA_Q (0.1f * P_H)
#define SCORR_K 0.1f
#define SCORR_N(x) ((x) * (x) * (x) * (x))  // expands to x^n with n=4

#define XSPH_C 0.01f
#define VORTICITY_EPS 3.0f

#define SHAKE(t) (5.0f * glm::sin(2.0f * t)) 
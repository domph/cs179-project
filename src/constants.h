/* Particle constants */

#define R  1.0f  // radius
#define H  (2  * R)
#define H2 (H  * H)
#define H3 (H  * H2)
#define H6 (H3 * H3)
#define H9 (H3 * H6)

#define MAX_NEIGHBORS 64


/* Math constants */

#define PI glm::pi<float>()
#define EPS 1.0E-5f  // numerical zero


/* Simulation constants */
#define SOLVER_ITERATIONS 3


/* Physics constants */
#define DT 0.02f
#define G  9.81f

#define POLY6_COEFF (315.0f / (64.0f * PI * H9))
#define SPIKY_COEFF (45.0f  / (PI * H6))

#define RHO_0 1000.0f
#define RHO_0_INV (1 / RHO_0)
#define RELAXATION_EPS 5.0f

#define DELTA_Q (0.2f * H)
#define SCORR_K 0.1f
#define SCORR_N(x) ((x) * (x) * (x) * (x))  // expands to x^n with n=4

#define XSPH_C 0.01f
#define VORTICITY_EPS 3.0f
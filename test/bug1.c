#include <stdio.h>

#include "api.h"
#include "utils.h"

int main(void) {
  int n = 2, m = 3, ms = 2;
  double H[4] = {4.0, 1.0, 1.0, 2.0};  // Non-diagonal SPD Hessian
  double f[2] = {2.0, 2.0};            // x_unc = -H^{-1} f = [-2/7, -6/7]
  double A[2] = {1.0, 1.0};            // 1 general row: -2 <= x0 + x1 <= 2
  double bu[3] = {0.5, 0.5, 2.0};
  // x_unc[-0.286, -0.857] violates bl[1] = -0.5
  double bl[3] = {-0.5, -0.5, -2.0};
  int sense[3] = {0, 0, 0};

  DAQPProblem qp = {.n = n,
                    .m = m,
                    .ms = ms,
                    .H = H,
                    .f = f,
                    .A = A,
                    .bupper = bu,
                    .blower = bl,
                    .sense = sense,
                    .break_points = NULL,
                    .nh = 1,
                    .problem_type = 0};
  DAQPSettings settings;
  daqp_default_settings(&settings);

  // For the configuration above, the optimal solution is as follows
  const char* expected_string =
      ",   \tExpected: x = [-0.375000, -0.500000], lam[1] = -0.625000\n";

  // --- Part A: Reusing normalized Rinv with DAQP_UPDATE_unconstrained ---
  DAQPWorkspace work = {0};
  work.settings = &settings;
  setup_daqp(&qp, &work, NULL);
  double x[2] = {0}, lam[3] = {0};
  DAQPResult res = {.x = x,
                    .lam = lam,
                    .fval = 0,
                    .soft_slack = 0,
                    .exitflag = 0,
                    .iter = 0,
                    .nodes = 0,
                    .solve_time = 0,
                    .setup_time = 0};

  // True optimum at bl[1] = -0.5: 4*x0 - 0.5 + 2 = 0 => x* = [-0.375, -0.5],
  // lam[1] = -0.625
  daqp_solve(&res, &work);
  printf("Part A (initial solve) : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
         lam[1]);
  printf(expected_string);

  // Re-solve the exact same problem via a vector update with
  // DAQP_UPDATE_unconstrained:
  int vec_mask = DAQP_UPDATE_v | DAQP_UPDATE_d | DAQP_UPDATE_sense |
                 DAQP_UPDATE_unconstrained;
  daqp_update_ldp(vec_mask, &work, &qp);
  daqp_solve(&res, &work);
  printf("Part A (vector update) : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
         lam[1]);
  printf(expected_string);
  work.settings = NULL;
  free_daqp_workspace(&work);
  free_daqp_ldp(&work);

  // --- Part B: Unconstrained setup followed by constrained vector update ---
  double f_int[2] = {0.5,
                     0.5};  // x_unc = [-1/14, -3/14] strictly inside bl..bu
  qp.f = f_int;
  DAQPWorkspace work2 = {0};
  work2.settings = &settings;

  // Setup with init_mask = DAQP_UPDATE_unconstrained (early return skips M &
  // Rinv norm):
  setup_daqp_main(&qp, &work2, NULL, DAQP_UPDATE_unconstrained);
  daqp_solve(&res, &work2);

  // Now update f to [2.0, 2.0] (constrained at x* = [-0.375, -0.5]) with a
  // normal vector update:
  qp.f = f;
  daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work2, &qp);
  daqp_solve(&res, &work2);
  printf("Part B (after unc exit): x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
         lam[1]);
  printf(expected_string);
  work2.settings = NULL;
  free_daqp_workspace(&work2);
  free_daqp_ldp(&work2);

  return 0;
}

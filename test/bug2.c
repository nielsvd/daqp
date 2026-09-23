#include <stdio.h>

#include "api.h"
#include "utils.h"

int main(void) {
  int n = 2, m = 2, ms = 2;
  double f[2] = {1.0, -1.0};
  double bu[2] = {1.0, 1.0};
  double bl[2] = {-1.0, -1.0};
  int sense[2] = {0, 0};

  // Pure LP: H == NULL
  DAQPProblem qp = {.n = n,
                    .m = m,
                    .ms = ms,
                    .H = NULL,
                    .f = f,
                    .A = NULL,
                    .bupper = bu,
                    .blower = bl,
                    .sense = sense,
                    .break_points = NULL,
                    .nh = 1,
                    .problem_type = 0};
  DAQPSettings settings;
  daqp_default_settings(&settings);
  DAQPWorkspace work = {0};
  work.settings = &settings;
  setup_daqp(&qp, &work, NULL);

  // Calling daqp_update_ldp with DAQP_UPDATE_Rinv when qp.H == NULL crashes
  // unpatched:
  daqp_update_ldp(DAQP_UPDATE_Rinv | DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp);
  double x[2] = {0}, lam[2] = {0};
  DAQPResult res = {.x = x,
                    .lam = lam,
                    .fval = 0,
                    .soft_slack = 0,
                    .exitflag = 0,
                    .iter = 0,
                    .nodes = 0,
                    .solve_time = 0,
                    .setup_time = 0};
  daqp_solve(&res, &work);
  printf("LP solution: x = [%.6f, %.6f], exitflag = %d\n", x[0], x[1],
         res.exitflag);
  work.settings = NULL;
  free_daqp_workspace(&work);
  free_daqp_ldp(&work);

  return 0;
}

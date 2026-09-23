#include <stdio.h>

#include "api.h"
#include "utils.h"

int main(void) {
  // Baseline QP, bl[1] active: x*=[-0.375, -0.5] and l*=[0, -0.625]
  const char* expected_string =
      ",   \tExpected: x = [-0.375000, -0.500000], lam[1] = -0.625000\n";

  ////////////////////////////////////////////////////////////////
  /// Reusing normalized Rinv with very particular update mask ///
  ////////////////////////////////////////////////////////////////
  {
    int n = 2, m = 2, ms = 2;
    double H[4] = {4.0, 1.0, 1.0, 2.0};  // Non-diagonal SPD Hessian
    double f[2] = {2.0, 2.0};            // x_unc = -H^{-1} f = [-2/7, -6/7]
    double bu[2] = {0.5, 0.5};
    double bl[2] = {-0.5, -0.5};
    int sense[2] = {0, 0};

    // Build QP
    DAQPProblem qp = {.n = n,
                      .m = m,
                      .ms = ms,
                      .H = H,
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

    // Compare initial DAQP result with the expected solution: correct!
    daqp_solve(&res, &work);
    printf(
        "=== Problem with DAQP_UPDATE_unconstrained + DAQP_UPDATE_sense ===\n");
    printf("Initial solve : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
           lam[1]);
    printf(expected_string);

    // Keep QP exactly the same, instruct vector update with
    // DAQP_UPDATE_unconstrained and DAQP_UPDATE_sense: result is now incorrect.
    int vec_mask = DAQP_UPDATE_v | DAQP_UPDATE_d | DAQP_UPDATE_sense |
                   DAQP_UPDATE_unconstrained;
    daqp_update_ldp(vec_mask, &work, &qp);
    daqp_solve(&res, &work);
    printf("Vector update : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
           lam[1]);
    printf(expected_string);
    work.settings = NULL;

    // Clean up
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
  }

  /////////////////////////////////////////////////////////////////
  /// Unconstrained setup followed by constrained vector update ///
  /////////////////////////////////////////////////////////////////
  {
    int n = 2, m = 3, ms = 2;
    double H[4] = {4.0, 1.0, 1.0, 2.0};  // Non-diagonal SPD Hessian
    double f[2] = {2.0, 2.0};            // x_unc = -H^{-1} f = [-2/7, -6/7]
    double A[2] = {1.0, 1.0};
    double bu[3] = {0.5, 0.5, 2.0};
    double bl[3] = {-0.5, -0.5, -2.0};
    int sense[3] = {0, 0, 0};

    // We'll consider the same problem, but now with linear constraints.
    // x_unc = [-1/14, -3/14], no constraints are active in optimum
    double f_int[2] = {0.5, 0.5};

    // Build QP
    DAQPProblem qp = {.n = n,
                      .m = m,
                      .ms = ms,
                      .H = H,
                      .f = f_int,
                      .A = A,
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

    // Result struct
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

    // Setup with init_mask = DAQP_UPDATE_unconstrained
    setup_daqp_main(&qp, &work, NULL, DAQP_UPDATE_unconstrained);
    daqp_solve(&res, &work);

    printf("\n=== Problem with DAQP_UPDATE_v and update_ldp ===\n");
    printf("Initial solve : x = [%.6f, %.6f]", x[0], x[1]);
    printf(",   \tExpected: x = [%.6f, %.6f]\n", -1.0 / 14.0, -3.0 / 14.0);

    // Now update f to [2.0, 2.0] (same problem as before), with vector update
    qp.f = f;
    daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp);
    daqp_solve(&res, &work);
    printf("Vector update : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
           lam[1]);
    printf(expected_string);
    work.settings = NULL;

    // Clean up
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
  }

  /////////////////////////////////////////////////////////////////////
  /// Unconstrained AVI setup followed by constrained vector update ///
  /////////////////////////////////////////////////////////////////////
  {
    int n = 2, m = 3, ms = 2;
    double H[4] = {4.0, 1.0, 0.0,
                   2.75};  // Asymmetric matrix (problem_type = 1)
    double f[2] = {
        2.0, 2.0};  // At bl[1] = -0.5: x* = [-0.375, -0.5], lam[1] = -0.625
    double A[2] = {1.0, 1.0};
    double bu[3] = {0.5, 0.5, 2.0};
    double bl[3] = {-0.5, -0.5, -2.0};
    int sense[3] = {0, 0, 0};
    // x_unc = [-7/88, -2/11], no constraints are active in optimum
    double f_int[2] = {0.5, 0.5};
    // Build asymmetric AVI
    DAQPProblem qp = {.n = n,
                      .m = m,
                      .ms = ms,
                      .H = H,
                      .f = f_int,
                      .A = A,
                      .bupper = bu,
                      .blower = bl,
                      .sense = sense,
                      .break_points = NULL,
                      .nh = 1,
                      .problem_type = 1};
    DAQPSettings settings;
    daqp_default_settings(&settings);
    DAQPWorkspace work = {0};
    work.settings = &settings;
    // Result struct
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
    // Setup with init_mask = DAQP_UPDATE_unconstrained
    setup_daqp_main(&qp, &work, NULL, DAQP_UPDATE_unconstrained);
    daqp_solve(&res, &work);
    printf("\n=== AVI problem with DAQP_UPDATE_v and update_ldp ===\n");
    printf("Initial solve : x = [%.6f, %.6f]", x[0], x[1]);
    printf(",   \tExpected: x = [%.6f, %.6f]\n", -7.0 / 88.0, -2.0 / 11.0);
    // Now update f to [2.0, 2.0], with vector update (crashes with SIGSEGV
    // unpatched)
    qp.f = f;
    daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp);
    daqp_solve(&res, &work);
    printf("Vector update : x = [%.6f, %.6f], lam[1] = %.6f", x[0], x[1],
           lam[1]);
    printf(expected_string);
    work.settings = NULL;
    // Clean up
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
  }

  return 0;
}

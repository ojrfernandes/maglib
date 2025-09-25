#include "sode.h"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

// Test parameters for Lorenz system
typedef struct {
    double sigma;
    double rho;
    double beta;
} lorenz_params;

// Test parameters for simple harmonic oscillator
typedef struct {
    double omega;
} sho_params;

// Systems for ODE solving

int lorenz_system(double *f, double *x, double t, void *aux) {
    lorenz_params *lp = (lorenz_params *)aux;
    f[0] = lp->sigma * (x[1] - x[0]);      // dx/dt = sigma * (y - x)
    f[1] = x[0] * (lp->rho - x[2]) - x[1]; // dy/dt = x * (rho - z) - y
    f[2] = x[0] * x[1] - lp->beta * x[2];  // dz/dt = x * y - beta * z
    return 0;
}

int sho_system(double *f, double *x, double t, void *aux) {
    sho_params *sp = (sho_params *)aux;
    f[0] = x[1];                          // dx/dt = v
    f[1] = -sp->omega * sp->omega * x[0]; // dv/dt = -omega^2 * x
    return 0;
}

int linear_system(double *f, double *x, double t, void *aux) {
    f[0] = -x[0]; // dx/dt = -x
    return 0;
}

// Function to test event detection
bool zero_crossing_event(double *x, double t, void *aux) {
    return (x[0] <= 0.0); // Event when x crosses zero
}

// Test fixture for sode tests
class SodeTest : public ::testing::Test {
  protected:
    void SetUp() override {
        tolerance = 1e-6; // Default tolerance for floating point comparisons
    }
    void TearDown() override {
        // Common teardown can be done here if needed
    }

    double tolerance;

    // Additional helper functions
    bool values_are_close(double a, double b, double tol) {
        return std::fabs(a - b) < tol;
    }
};

// Test basic constructor and configuration
TEST_F(SodeTest, ConstructorAndConfiguration) {
    // Test for solver types
    sode solver_fb(SODE_RK56_FB, 3);
    sode solver_ck(SODE_RK56_CK, 3);
    sode solver_dp(SODE_RK78_DP, 3);

    // Test for stepsize configuration
    solver_fb.configure(0.01, 1e-6, 0.1);
    solver_ck.configure(0.01, 1e-6, 0.1);
    solver_dp.configure(0.01, 1e-6, 0.1);

    // Test for tolerance configuration
    solver_fb.configure(1e-8, 1e-12, 1e-12, 0.9);
    solver_ck.configure(1e-8, 1e-12, 1e-12, 0.9);
    solver_dp.configure(1e-8, 1e-12, 1e-12, 0.9);

    EXPECT_TRUE(true); // If we reach here, constructors and configurations worked
}

// Test simple linear ODE
TEST_F(SodeTest, LinearODE) {
    sode solver(SODE_RK56_CK, 1);
    solver.set_system(linear_system);
    solver.configure(1e-3, 1e-6, 0.1);
    solver.configure(1e-10, 1e-12, 1e-12, 0.9);

    double x[1] = {1.0}; // Initial condition
    double t = 0.0;      // Initial time
    double t_end = 1.0;  // End time

    solver.reset();
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        status = solver.evolve(x, &t, t_end, 0, nullptr);
    }

    // Analytical solution: x(t) = exp(-t)
    double expected = std::exp(-t_end);

    EXPECT_EQ(status, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(x[0], expected, tolerance));
    EXPECT_TRUE(values_are_close(t, t_end, tolerance));
}

// Test simple harmonic oscillator (conserves energy)
TEST_F(SodeTest, SimpleHarmonicOscillator) {
    sode solver(SODE_RK56_CK, 2);
    solver.set_system(sho_system);
    solver.configure(1e-3, 1e-6, 1e-1);
    solver.configure(1e-10, 1e-12, 1e-12, 0.9);

    sho_params params = {1.0}; // omega = 1
    double x[2] = {1.0, 0.0};  // Initial: x=1, v=0
    double t = 0.0;            // Initial time
    double t_end = 2 * M_PI;   // One full period

    // E = 0.5*(m*v^2 + k*x^2) = 0.5*(v^2 + omega^2*x^2)
    double initial_energy = 0.5 * (x[1] * x[1] + params.omega * params.omega * x[0] * x[0]);

    solver.reset();
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        status = solver.evolve(x, &t, t_end, 0, &params);
    }

    // After one period, should return to initial position
    double final_energy = 0.5 * (x[1] * x[1] + params.omega * params.omega * x[0] * x[0]);

    EXPECT_EQ(status, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(x[0], 1.0, tolerance));                    // Position should return
    EXPECT_TRUE(values_are_close(x[1], 0.0, tolerance));                    // Velocity should return
    EXPECT_TRUE(values_are_close(initial_energy, final_energy, tolerance)); // Energy conservation
}

// Test Lorenz system evolution (chaotic but should not blow up)
TEST_F(SodeTest, LorenzSystemStability) {
    sode solver(SODE_RK56_CK, 3);
    solver.set_system(lorenz_system);
    solver.configure(1e-3, 1e-6, 1e-1);
    solver.configure(1e-10, 1e-12, 1e-12, 0.9);

    lorenz_params params = {10.0, 28.0, 8.0 / 3.0}; // Parameters: sigma, rho, beta
    double x[3] = {0.1, 0.1, 0.1};                  // Initial condition
    double t = 0.0;                                 // Initial time
    double t_end = 10.0;                            // Evolve for some time

    solver.reset();
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        status = solver.evolve(x, &t, t_end, 0, &params);
    }

    EXPECT_EQ(status, SODE_SUCCESS_TIME);
    // Check that solution hasn't blown up (Lorenz attractor is bounded)
    EXPECT_LT(std::abs(x[0]), 50.0);
    EXPECT_LT(std::abs(x[1]), 50.0);
    EXPECT_LT(std::abs(x[2]), 50.0);
}

// Test event detection
TEST_F(SodeTest, EventDetection) {
    sode solver(SODE_RK56_CK, 2);
    solver.set_system(sho_system);
    solver.set_monitor(zero_crossing_event);
    solver.configure(1e-3, 1e-6, 1e-1);
    solver.configure(1e-10, 1e-12, 1e-12, 0.9);

    sho_params params = {1.0}; // omega = 1
    double x[2] = {1.0, 0.0};  // Start at x=1 (monitor = true)
    double t = 0.0;
    double t_end = 10.0;

    solver.reset();
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        status = solver.evolve(x, &t, t_end, 1, &params);
    }

    // Should stop when x crosses zero (at t = pi/2)
    EXPECT_EQ(status, SODE_SUCCESS_MONITOR);
    EXPECT_TRUE(values_are_close(x[0], 0.0, tolerance));   // Should be at zero crossing
    EXPECT_TRUE(values_are_close(t, M_PI / 2, tolerance)); // Should occur at pi/2
}

// Test different solver types give similar results
TEST_F(SodeTest, SolverConsistency) {
    lorenz_params params = {10.0, 28.0, 8.0 / 3.0};
    double t_end = 1.0;

    // Solve with different methods
    std::vector<sode_type> solvers = {SODE_RK56_FB, SODE_RK56_CK, SODE_RK78_DP};
    std::vector<std::vector<double>> results;

    for (auto solver_type : solvers) {
        sode solver(solver_type, 3);
        solver.set_system(lorenz_system);
        solver.configure(1e-3, 1e-6, 1e-1);
        solver.configure(1e-10, 1e-12, 1e-12, 0.9);

        double x[3] = {0.1, 0.1, 0.1};
        double t = 0.0;

        solver.reset();
        int status = SODE_CONTINUE_GOOD_STEP;
        while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
            status = solver.evolve(x, &t, t_end, 0, &params);
        }

        results.push_back({x[0], x[1], x[2]});
    }

    // All solvers should give similar results
    for (int i = 0; i < 3; i++) {
        EXPECT_TRUE(values_are_close(results[0][i], results[1][i], tolerance));
        EXPECT_TRUE(values_are_close(results[0][i], results[2][i], tolerance));
    }
}
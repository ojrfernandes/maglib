#include <cstring>
#include <gtest/gtest.h>

#include "../maglit/m3dc1_source.h"
#include "../maglit/maglit.h"

// Test fixture for maglit functionality
class MaglitTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::snprintf(source_path, sizeof(source_path),
                      "%s/C1.h5", TEST_DATA_DIR);
        shape_path  = std::string(TEST_DATA_DIR) + "/tcabr_first_wall.txt";
        source_vac  = new M3DC1Source(source_path, -1);
        tracer_vac  = new maglit(*source_vac);
        tracer_vac->configure(0.01, 1e-6, 0.1);
        source_resp = new M3DC1Source(source_path, 1);
        tracer_resp = new maglit(*source_resp);
        tracer_resp->configure(0.01, 1e-6, 0.1);
    }

    static void TearDownTestSuite() {
        // Delete tracers before sources (tracers hold references to sources).
        delete tracer_vac;
        delete tracer_resp;
        delete source_vac;
        delete source_resp;
        tracer_vac  = nullptr;
        tracer_resp = nullptr;
        source_vac  = nullptr;
        source_resp = nullptr;
    }

    void SetUp() override { tolerance = 1e-6; }

    double              tolerance;
    static char         source_path[256];
    static std::string  shape_path;
    static M3DC1Source *source_vac;
    static M3DC1Source *source_resp;
    static maglit      *tracer_vac;
    static maglit      *tracer_resp;

    bool values_are_close(double a, double b, double tol) {
        return std::fabs(a - b) < tol;
    }
};

// Static member definitions
M3DC1Source *MaglitTest::source_vac  = nullptr;
M3DC1Source *MaglitTest::source_resp = nullptr;
maglit      *MaglitTest::tracer_vac  = nullptr;
maglit      *MaglitTest::tracer_resp = nullptr;
char         MaglitTest::source_path[256];
std::string  MaglitTest::shape_path;

// ==================== MAGLIT INTEGRATION TESTS ====================

// Test constructor and configuration
TEST_F(MaglitTest, ConstructorAndConfiguration) {
    ASSERT_NE(tracer_vac, nullptr);
    EXPECT_TRUE(true); // If we reach here, constructors and configurations worked
}

// Test magnetic field evaluation at a specific point
TEST_F(MaglitTest, MagneticFieldEvaluation) {
    double x[3] = {0.7, 0.0, 0.0}; // R, phi, Z
    double B[3] = {0.0, 0.0, 0.0}; // B_R, B_phi, B_Z

    bool success = tracer_vac->calc_mag_field(x, B);
    ASSERT_TRUE(success);
    // Basic checks on magnetic field values
    // Compared to fusion-io python routines output
    EXPECT_TRUE(values_are_close(B[0], -0.00401528, tolerance)); // B_R
    EXPECT_TRUE(values_are_close(B[1], -0.95980778, tolerance)); // B_phi
    EXPECT_TRUE(values_are_close(B[2], -0.14235661, tolerance)); // B_Z
}

// Test poloidal flux evaluation at a specific point
TEST_F(MaglitTest, PsiFieldEvaluation) {
    double x[3] = {0.7, 0.0, 0.0}; // R, phi, Z
    double psi = 0.0;

    tracer_vac->psi_eval(x[0], x[1], x[2], &psi);
    // Basic checks on psi values
    // Compared to fusion-io python routines output
    EXPECT_TRUE(values_are_close(psi, -0.00988012, tolerance)); // psi
}

// Test normalized poloidal flux evaluation at a specific point
TEST_F(MaglitTest, PsiNFieldEvaluation) {
    double x[3] = {0.7, 0.0, 0.0}; // R, phi, Z
    double psiN = 0.0;

    tracer_vac->psin_eval(x[0], x[1], x[2], &psiN);
    // Basic checks on psiN values
    // Compared to fusion-io python routines output
    EXPECT_TRUE(values_are_close(psiN, 0.328436, tolerance)); // psiN
}

// Test inverse map functionality
TEST_F(MaglitTest, InverseMap) {
    double R = 0.7;
    double Z = 0.0;
    double phi = 0.0;
    double phi_max = M_PI / 2; // 90 degrees

    tracer_vac->reset();
    tracer_vac->inverse_map(false); // Forward map
    int status_fwd = SODE_CONTINUE_GOOD_STEP;
    while (status_fwd == SODE_CONTINUE_GOOD_STEP || status_fwd == SODE_CONTINUE_BAD_STEP) {
        status_fwd = tracer_vac->step(R, Z, phi, phi_max, 0);
    }
    EXPECT_EQ(status_fwd, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(phi, phi_max, tolerance));

    tracer_vac->reset();
    tracer_vac->inverse_map(true); // Inverse map
    phi_max += phi;                // Inverting the map does not reset phi nor makes it decrease when iterating
    int status_inv = SODE_CONTINUE_GOOD_STEP;
    while (status_inv == SODE_CONTINUE_GOOD_STEP || status_inv == SODE_CONTINUE_BAD_STEP) {
        status_inv = tracer_vac->step(R, Z, phi, phi_max, 0);
    }

    EXPECT_EQ(status_inv, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(R, 0.7, tolerance));
    EXPECT_TRUE(values_are_close(Z, 0.0, tolerance));
}

// ==================== MAGLIT MONITORING TESTS ====================

TEST_F(MaglitTest, ColliderMonitorIntegration) {
    // Initial conditions — point inside the chaotic region
    double last_R = 0;
    double last_Z = 0;
    double R = 0.5;
    double Z = -0.225;
    double phi = 0.0;
    double phi_max = 10000 * 2 * M_PI; // Large phi to ensure crossing

    tracer_resp->reset();
    tracer_resp->set_monitor(shape_path);
    tracer_resp->inverse_map(false); // Forward map

    // Integrate while checking for wall crossings (handled internally by monitor)
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        last_R = R;
        last_Z = Z;
        status = tracer_resp->step(R, Z, phi, phi_max, -1);
    }

    // If collider works, integration should stop due to monitor trigger
    EXPECT_EQ(status, SODE_SUCCESS_MONITOR);
    EXPECT_TRUE(tracer_resp->get_boundary().inside(last_R, last_Z));
    EXPECT_FALSE(tracer_resp->get_boundary().inside(R, Z));
}

TEST_F(MaglitTest, InverseColliderMonitorIntegration) {
    // Initial conditions — point inside the chaotic region
    double last_R = 0;
    double last_Z = 0;
    double R = 0.5;
    double Z = -0.225;
    double phi = 0.0;
    double phi_max = 10000 * 2 * M_PI; // Large phi to ensure crossing

    tracer_resp->reset();
    tracer_resp->set_monitor(shape_path);
    tracer_resp->inverse_map(true); // Inverse map

    // Integrate while checking for wall crossings (handled internally by monitor)
    int status = SODE_CONTINUE_GOOD_STEP;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        last_R = R;
        last_Z = Z;
        status = tracer_resp->step(R, Z, phi, phi_max, -1);
    }

    // If collider works, integration should stop due to monitor trigger
    EXPECT_EQ(status, SODE_SUCCESS_MONITOR);
    EXPECT_TRUE(tracer_resp->get_boundary().inside(last_R, last_Z));
    EXPECT_FALSE(tracer_resp->get_boundary().inside(R, Z));
    tracer_resp->inverse_map(false); // restore default
}

TEST_F(MaglitTest, OnBoundaryColliderMonitorIntegration) {
    // Initial conditions — point on the boundary
    double last_R = 0;
    double last_Z = 0;
    double R = 0.435;
    double Z = -0.235;
    double phi = 0.0;
    double phi_max = 10000 * 2 * M_PI; // Large phi to ensure crossing

    tracer_resp->reset();
    tracer_resp->set_monitor(shape_path);
    tracer_resp->inverse_map(false); // Forward map

    // Integrate while checking for wall crossings (handled internally by monitor)
    int status = SODE_CONTINUE_GOOD_STEP;
    int step_count = 0;
    while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP) {
        last_R = R;
        last_Z = Z;
        status = tracer_resp->step(R, Z, phi, phi_max, -1);
        step_count++;
    }

    // If collider works, integration should stop due to monitor trigger
    EXPECT_GT(step_count, 0);
    EXPECT_EQ(status, SODE_SUCCESS_MONITOR);
    EXPECT_TRUE(tracer_resp->get_boundary().inside(last_R, last_Z));
    EXPECT_FALSE(tracer_resp->get_boundary().inside(R, Z));
}
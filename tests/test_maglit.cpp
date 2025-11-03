#include "maglit.h"
#include "sode.h"
#include <cstring>
#include <gtest/gtest.h>

// Test fixture for maglit functionality
class MaglitTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::snprintf(source_path, sizeof(source_path),
                      "%s/C1.h5", TEST_DATA_DIR);
        tracer = new maglit(source_path, FIO_M3DC1_SOURCE, -1);
        tracer->configure(0.01, 1e-6, 0.1);
    }

    static void TearDownTestSuite() {
        delete tracer;
        tracer = nullptr;
    }

    void SetUp() override {
        tolerance = 1e-6;
    }

    double         tolerance;
    static char    source_path[256];
    static maglit *tracer;

    bool values_are_close(double a, double b, double tol) {
        return std::fabs(a - b) < tol;
    }
};

// Static member definitions
maglit *MaglitTest::tracer = nullptr;
char    MaglitTest::source_path[256];

// ==================== MAGLIT TESTS ====================

// Test constructor and configuration
TEST_F(MaglitTest, ConstructorAndConfiguration) {
    ASSERT_NE(tracer, nullptr);
    EXPECT_TRUE(true); // If we reach here, constructors and configurations worked
}

// Test magnetic field evaluation at a specific point
TEST_F(MaglitTest, MagneticFieldEvaluation) {
    double x[3] = {0.7, 0.0, 0.0}; // R, phi, Z
    double B[3] = {0.0, 0.0, 0.0}; // B_R, B_phi, B_Z

    bool success = tracer->calc_mag_field(x, B);
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

    tracer->psi_eval(x[0], x[1], x[2], &psi);
    // Basic checks on psi values
    // Compared to fusion-io python routines output
    EXPECT_TRUE(values_are_close(psi, -0.00988012, tolerance)); // psi
}

// // Test normalized poloidal flux evaluation at a specific point
// TEST_F(MaglitTest, PsiNFieldEvaluation) {
//     double x[3] = {0.7, 0.0, 0.0}; // R, phi, Z
//     double psiN = 0.0;

//     tracer->psin_eval(x[0], x[1], x[2], &psiN);
//     // Basic checks on psiN values
//     // Compared to fusion-io python routines output
//     EXPECT_TRUE(values_are_close(psiN, 0.328433, tolerance)); // psiN
// }

// Test inverse map functionality
TEST_F(MaglitTest, InverseMap) {
    double R = 0.7;
    double Z = 0.0;
    double phi = 0.0;
    double phi_max = M_PI / 2; // 90 degrees

    tracer->reset();
    tracer->alloc_hint();
    tracer->inverse_map(false); // Forward map
    int status_fwd = SODE_CONTINUE_GOOD_STEP;
    while (status_fwd == SODE_CONTINUE_GOOD_STEP || status_fwd == SODE_CONTINUE_BAD_STEP) {
        status_fwd = tracer->step(R, Z, phi, phi_max, 0);
    }
    EXPECT_EQ(status_fwd, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(phi, phi_max, tolerance));

    tracer->reset();
    tracer->inverse_map(true); // Inverse map
    phi_max += phi;            // Inverting the map does not reset phi nor makes it decrease when iterating
    int status_inv = SODE_CONTINUE_GOOD_STEP;
    while (status_inv == SODE_CONTINUE_GOOD_STEP || status_inv == SODE_CONTINUE_BAD_STEP) {
        status_inv = tracer->step(R, Z, phi, phi_max, 0);
    }
    tracer->clear_hint();

    EXPECT_EQ(status_inv, SODE_SUCCESS_TIME);
    EXPECT_TRUE(values_are_close(R, 0.7, tolerance));
    EXPECT_TRUE(values_are_close(Z, 0.0, tolerance));
}
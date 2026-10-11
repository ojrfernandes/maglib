#include <filesystem>
#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Include the mfgen headers
#include "../maglit/m3dc1_source.h"
#include "../mfgen/src/input_read.h"
#include "../mfgen/src/manifold.h"

// Turning-angle refinement threshold passed to manifold (radians); 20 deg
static constexpr double THETA_LIM = 20.0 * M_PI / 180.0;
// Magnetic axis of tests/data/C1.h5 (scalars/xmag, zmag)
static constexpr double MAG_AXIS_R = 0.63088644;
static constexpr double MAG_AXIS_Z = -0.00082328;

// Test fixture for mfgen tests
class MfgenTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::string src_path = std::string(TEST_DATA_DIR) + "/C1.h5";
        source = new M3DC1Source(src_path.c_str(), 1);
        tracer = new maglit(*source);
        tracer->configure(1e-2, 1e-6, 1e-2);
    }

    static void TearDownTestSuite() {
        delete tracer;
        delete source;
        tracer = nullptr;
        source = nullptr;
    }

    static M3DC1Source *source;
    static maglit      *tracer;
    static std::string  source_path;

    void SetUp() override {
        // Create temporary directory for test files
        test_dir = "test_input_files";
        std::filesystem::create_directories(test_dir);

        // Create valid test input file
        valid_input_file = test_dir + "/valid_input.txt";
        createValidInputFile();

        // Create input file with optional x-point coordinates
        xpoint_input_file = test_dir + "/xpoint_input.txt";
        createXPointInputFile();

        // Create invalid test files
        invalid_input_file = test_dir + "/invalid_input.txt";
        malformed_input_file = test_dir + "/malformed_input.txt";
        empty_input_file = test_dir + "/empty_input.txt";
        createInvalidInputFiles();

        // Create edge case test input file
        edge_case_input_file = test_dir + "/edge_case_input.txt";
        comment_whitespace_input_file = test_dir + "/comment_input.txt";
        createEdgeCaseInputFile();
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove_all(test_dir);
    }

    void createValidInputFile() {
        std::ofstream valid_input(valid_input_file);
        valid_input
            << "################################ MFGEN #################################\n"
            << "#\n"
            << "#=============== I/O FILES\n"
            << "        output_path = test_output.dat      # output file path and name\n"
            << "#\n"
            << "[M3DC1 SOURCE]\n"
            << "        nsources    = 1\n"
            << "        source_0    = test_source.h5\n"
            << "        timeslice_0 = 1\n"
            << "        phase_0     = 0.0\n"
            << "        amplitude_0 = 1.0\n"
            << "#\n"
            << "#=============== TRACING PARAMETERS\n"
            << "#\n"
            << "        manifold  = 0           # unstable manifold=0; stable manifold=1\n"
            << "        method = 1                      # exact manifold = 0, interpolant method = 1\n"
            << "        Phi = 0                       # Poincare section toroidal coordinate (deg)\n"
            << "#\n"
            << "#=============== MULTIPLE POINCARE SECTIONS\n"
            << "#\n"
            << "        nSections = 3              # number of Poincare sections\n"
            << "        phi_0 = 0                   # starting toroidal coordinate (deg)\n"
            << "        phi_1 = 180                 # ending toroidal coordinate (deg)\n"
            << "#\n"
            << "#=============== ADDITIONAL PARAMETERS\n"
            << "#\n"
            << "        epsilon = 1e-6            # perturbation size for manifold tracing\n"
            << "        nSegments = 5              # number of segments to trace\n"
            << "        l_lim = 0.01                # maximum segment length\n"
            << "        theta_lim = 10             # maximum angle between segment arcs (deg)\n"
            << "        h_init = 1e-2              # initial step-size for integration\n"
            << "        h_min = 1e-6               # minimum step-size for integration\n"
            << "        h_max = 1e-2               # maximum step-size for integration\n"
            << "        h_deriv = 1e-3             # step size for derivative calculation\n"
            << "        n_tol = 1e-14                 # tolerance for Newton's method\n"
            << "        max_iter = 100             # maximum iterations for Newton's method\n"
            << "        precision = 1e-8           # precision for calculations\n"
            << "        max_insertions = 50         # maximum insertions allowed\n"
            << "        verbose = 0                 # verbosity flag (0 or 1)\n";
        valid_input.close();
    }

    void createXPointInputFile() {
        std::ofstream f(xpoint_input_file);
        f   << "output_path = test_output.dat\n"
            << "[M3DC1 SOURCE]\n"
            << "nsources    = 1\n"
            << "source_0    = test_source.h5\n"
            << "timeslice_0 = 1\n"
            << "phase_0     = 0.0\n"
            << "amplitude_0 = 1.0\n"
            << "manifold  = 0\n"
            << "method = 1\n"
            << "Phi = 0\n"
            << "nSections = 1\n"
            << "phi_0 = 0\n"
            << "phi_1 = 0\n"
            << "epsilon = 1e-6\n"
            << "nSegments = 2\n"
            << "l_lim = 0.005\n"
            << "theta_lim = 20\n"
            << "h_init = 1e-2\n"
            << "h_min = 1e-6\n"
            << "h_max = 1e-2\n"
            << "h_deriv = 1e-8\n"
            << "n_tol = 1e-14\n"
            << "max_iter = 50\n"
            << "precision = 1e-14\n"
            << "max_insertions = 100\n"
            << "verbose = 0\n"
            << "R_xPoint = 0.497999\n"
            << "Z_xPoint = -0.218603\n";
        f.close();
    }

    void createInvalidInputFiles() {
        // File with missing parameters
        std::ofstream invalid_file(invalid_input_file);
        invalid_file << "# Missing required parameters\n"
                     << "[M3DC1 SOURCE]\n"
                     << "nsources = 1\n"
                     << "source_0 = test.h5\n"
                     << "# Missing output_path and numeric parameters\n";
        invalid_file.close();

        // File with malformed syntax
        std::ofstream malformed_file(malformed_input_file);
        malformed_file << "output_path test.dat  # Missing equals sign\n"
                       << "source_0 =            # Missing value\n"
                       << "num_threads = abc      # Invalid numeric value\n"
                       << "gridMin = 1.2.3        # Invalid float\n";
        malformed_file.close();

        // Empty file
        std::ofstream empty_file(empty_input_file);
        empty_file.close();
    }

    void createEdgeCaseInputFile() {
        // Edge case values — manifold=-1 and other invalid entries cause parsing to fail
        std::ofstream edge_case_input(edge_case_input_file);
        edge_case_input
            << "output_path = a\n"     // Single character
            << "[M3DC1 SOURCE]\n"
            << "nsources    = 1\n"
            << "source_0    = \n"      // Empty string — will fail validation
            << "timeslice_0 = 0\n"
            << "phase_0     = 0.0\n"
            << "amplitude_0 = 1.0\n"
            << "manifold  = -1\n"      // negative value
            << "method = 9999\n"       // very large value
            << "Phi = -360\n"          // negative angle
            << "nSections = 0\n"       // zero sections
            << "phi_0 = 360\n"         // large angle
            << "phi_1 = -360\n"        // negative angle
            << "epsilon = -1e-6\n"     // negative perturbation size
            << "nSegments = 0\n"       // zero segments
            << "l_lim = -0.01\n"       // negative maximum segment length
            << "theta_lim = -10\n"     // negative angle limit
            << "h_init = -1e-2\n"      // negative step-size
            << "h_min = -1e-6\n"       // negative minimum step-size
            << "h_max = -1e-2\n"       // negative maximum step-size
            << "h_deriv = -1e-3\n"     // negative step size for derivative calculation
            << "n_tol = -10\n"         // negative tolerance for Newton's method
            << "max_iter = -100\n"     // negative maximum iterations for Newton's method
            << "precision = -1e-8\n"   // negative precision for calculations
            << "max_insertions = -5\n" // negative maximum insertions allowed
            << "verbose = 1\n";        // verbosity flag set to 1

        edge_case_input.close();

        // Comments and whitespace variations
        std::ofstream comment_whitespace_input(comment_whitespace_input_file);
        comment_whitespace_input
            << "# This is a comment line\n"
            << "\n"
            << "   # Indented comment\n"
            << "output_path=test.dat# No spaces around equals\n"
            << "[M3DC1 SOURCE]\n"
            << "\t\tnsources    = 1\n"
            << "\t\tsource_0    = test.h5\t\t# Trailing comment with tabs\n"
            << "\t\ttimeslice_0 = 1 # Comment\n"
            << "\t\tphase_0     = 0.0\n"
            << "\t\tamplitude_0 = 1.0\n"
            << "\tmanifold\t=\t1\t\n"  // Tabs instead of spaces
            << "method = 1\n"
            << "Phi = 0 ## comment #\n"
            << "nSections = 1\n"
            << "phi_0 = 0\n"
            << "phi_1 = 1\n"
            << "epsilon = 1e-6\n"
            << "nSegments = 1\n"
            << "l_lim =0.01\n"
            << "theta_lim = 10\n"
            << "h_init = 1e-2\n"
            << "h_min=1e-6\n"
            << "h_max = 1e-2\n"
            << "h_deriv = 1e-3\n"
            << "n_tol = 10\n"
            << "max_iter = 100\n"
            << "precision = 1e-8\n"
            << "max_insertions = 5\n";
    }

    std::string test_dir;
    std::string valid_input_file;
    std::string xpoint_input_file;
    std::string invalid_input_file;
    std::string malformed_input_file;
    std::string empty_input_file;
    std::string edge_case_input_file;
    std::string comment_whitespace_input_file;
    std::string nonexistent_file = "nonexistent_file.txt";
};

// ==================== INPUT_READ TESTS ====================/

// Test: Constructor initializes reading_path correctly
TEST_F(MfgenTest, InputRead_Constructor) {
    std::string test_path = "test_path.txt";
    input_read  reader(test_path);

    // test that the constructor doesn't crash
    EXPECT_NO_THROW(input_read reader2(test_path));
}

// Test: Valid input file reading
TEST_F(MfgenTest, InputRead_ValidFile) {
    input_read reader(valid_input_file);
    bool       result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should successfully read valid input file";

    // Check string parameters
    EXPECT_EQ(reader.components[0].path, "test_source.h5");
    EXPECT_EQ(reader.output_path, "test_output.dat");

    // Check numeric parameters
    EXPECT_EQ(reader.components[0].timeslice, 1);
    EXPECT_EQ(reader.manifold, 0);
    EXPECT_EQ(reader.method, 1);
    EXPECT_DOUBLE_EQ(reader.Phi, 0.0);
    EXPECT_EQ(reader.nSections, 3);
    EXPECT_DOUBLE_EQ(reader.phi_0, 0.0);
    EXPECT_DOUBLE_EQ(reader.phi_1, 180.0);
    EXPECT_DOUBLE_EQ(reader.epsilon, 1e-6);
    EXPECT_EQ(reader.nSegments, 5);
    EXPECT_DOUBLE_EQ(reader.l_lim, 0.01);
    EXPECT_DOUBLE_EQ(reader.theta_lim, 10.0);
    EXPECT_DOUBLE_EQ(reader.h_init, 1e-2);
    EXPECT_DOUBLE_EQ(reader.h_min, 1e-6);
    EXPECT_DOUBLE_EQ(reader.h_max, 1e-2);
    EXPECT_DOUBLE_EQ(reader.h_deriv, 1e-3);
    EXPECT_DOUBLE_EQ(reader.n_tol, 1e-14);
    EXPECT_EQ(reader.max_iter, 100);
    EXPECT_DOUBLE_EQ(reader.precision, 1e-8);
    EXPECT_EQ(reader.max_insertions, 50);
    EXPECT_EQ(reader.verbose, 0);
    // Optional keys absent → default to 0 (HDF5 fallback will be used)
    EXPECT_DOUBLE_EQ(reader.R_xPoint, 0.0);
    EXPECT_DOUBLE_EQ(reader.Z_xPoint, 0.0);
}

// Test: Invalid file handling
TEST_F(MfgenTest, InputRead_InvalidFile) {
    input_read reader(invalid_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for invalid file";
}

// Test: Nonexistent file handling
TEST_F(MfgenTest, InputRead_NonexistentFile) {
    input_read reader(nonexistent_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for nonexistent file";
}

// Test: Malformed file handling
TEST_F(MfgenTest, InputRead_MalformedFile) {
    input_read reader(malformed_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for malformed file";
}

// Test: Empty file handling
TEST_F(MfgenTest, InputRead_EmptyFile) {
    input_read reader(empty_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for empty file";
}

// Test: Edge case values — file contains invalid entries (manifold=-1, method=9999,
// negative step sizes, etc.) so parsing must fail.
TEST_F(MfgenTest, InputRead_EdgeCaseValues) {
    input_read reader(edge_case_input_file);
    EXPECT_FALSE(reader.readInputFile());
}

// Test: Comments and whitespace handling
TEST_F(MfgenTest, InputRead_CommentsAndWhitespace) {
    input_read reader(comment_whitespace_input_file);
    bool       result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should handle comments and whitespace correctly";

    if (result) {
        EXPECT_EQ(reader.components[0].path, "test.h5");
        EXPECT_EQ(reader.output_path, "test.dat");
        EXPECT_EQ(reader.components[0].timeslice, 1);
        EXPECT_EQ(reader.manifold, 1);
        EXPECT_EQ(reader.method, 1);
        EXPECT_DOUBLE_EQ(reader.Phi, 0.0);
        EXPECT_EQ(reader.nSections, 1);
        EXPECT_DOUBLE_EQ(reader.phi_0, 0.0);
        EXPECT_DOUBLE_EQ(reader.phi_1, 1.0);
        EXPECT_DOUBLE_EQ(reader.epsilon, 1e-6);
        EXPECT_EQ(reader.nSegments, 1);
        EXPECT_DOUBLE_EQ(reader.l_lim, 0.01);
        EXPECT_DOUBLE_EQ(reader.theta_lim, 10.0);
        EXPECT_DOUBLE_EQ(reader.h_init, 1e-2);
        EXPECT_DOUBLE_EQ(reader.h_min, 1e-6);
        EXPECT_DOUBLE_EQ(reader.h_max, 1e-2);
        EXPECT_DOUBLE_EQ(reader.h_deriv, 1e-3);
        EXPECT_DOUBLE_EQ(reader.n_tol, 10);
        EXPECT_EQ(reader.max_iter, 100);
        EXPECT_DOUBLE_EQ(reader.precision, 1e-8);
        EXPECT_EQ(reader.max_insertions, 5);
    }
}

// Test: Multiple reads from same object
TEST_F(MfgenTest, InputRead_MultipleReads) {
    input_read reader(valid_input_file);

    // First read
    bool result1 = reader.readInputFile();
    EXPECT_TRUE(result1);

    // Store values from first read
    std::string first_source = reader.components[0].path;
    int         first_segments = reader.nSegments;

    // Second read should give same results
    bool result2 = reader.readInputFile();
    EXPECT_TRUE(result2);
    EXPECT_EQ(reader.components[0].path, first_source);
    EXPECT_EQ(reader.nSegments, first_segments);
}

// Test: HDF5 file reading (assuming test HDF5 file exists)
TEST_F(MfgenTest, InputRead_HDF5File) {
    input_read reader(valid_input_file);
    reader.nsources = 1;
    reader.components.resize(1);
    reader.components[0].path = std::string(TEST_DATA_DIR) + "/C1.h5";
    bool result = reader.readHDF5File();
    EXPECT_TRUE(result) << "Should successfully read HDF5 file";

    // Check that coordinates are within reasonable bounds
    EXPECT_NEAR(reader.R_xPoint, 0.497998952865601, 1e-6);
    EXPECT_NEAR(reader.Z_xPoint, -0.218603119254112, 1e-6);
}

// Test: R_xPoint / Z_xPoint default to 0 when absent from input file
TEST_F(MfgenTest, InputRead_XPoint_DefaultsToZero) {
    input_read reader(valid_input_file);
    bool result = reader.readInputFile();

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(reader.R_xPoint, 0.0);
    EXPECT_DOUBLE_EQ(reader.Z_xPoint, 0.0);
}

// Test: R_xPoint / Z_xPoint are parsed when present in input file
TEST_F(MfgenTest, InputRead_XPoint_ParsedFromFile) {
    input_read reader(xpoint_input_file);
    bool result = reader.readInputFile();

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(reader.R_xPoint, 0.497999);
    EXPECT_DOUBLE_EQ(reader.Z_xPoint, -0.218603);
}

// ==================== MANIFOLD TESTS ====================/

// Static member definitions
M3DC1Source *MfgenTest::source = nullptr;
maglit      *MfgenTest::tracer = nullptr;

// Test: Constructor and basic initialization
TEST_F(MfgenTest, Manifold_Constructor) {
    // Test stable manifold at phi=0
    EXPECT_NO_THROW(manifold mf(*tracer, 0, 0));

    // Test unstable manifold at phi=90
    EXPECT_NO_THROW(manifold mf(*tracer, M_PI / 2, 1));
}

// Test: Configure parameters
TEST_F(MfgenTest, Manifold_Configure) {
    manifold mf(*tracer, 0, 0);
    EXPECT_NO_THROW(mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 50));
}

// Test: set_branch accepts valid ±1 values
TEST_F(MfgenTest, Manifold_SetBranch_Valid) {
    manifold mf(*tracer, 0, 0);
    EXPECT_NO_THROW(mf.set_branch( 1,  1));
    EXPECT_NO_THROW(mf.set_branch(-1,  1));
    EXPECT_NO_THROW(mf.set_branch( 1, -1));
    EXPECT_NO_THROW(mf.set_branch(-1, -1));
}

// Test: set_branch rejects values other than ±1
TEST_F(MfgenTest, Manifold_SetBranch_InvalidValueRejected) {
    manifold mf(*tracer, 0, 0);
    EXPECT_THROW(mf.set_branch(0,  1), std::invalid_argument);
    EXPECT_THROW(mf.set_branch(1,  2), std::invalid_argument);
    EXPECT_THROW(mf.set_branch(2, -1), std::invalid_argument);
}

// Test: set_branch(-1, -1) produces a pivot point mirrored through xPoint
TEST_F(MfgenTest, Manifold_SetBranch_FullFlipMirrorspivot) {
    const double R_xp = 0.4979691771716279;
    const double Z_xp = -0.2185980054447758;
    const double eps  = 1e-6;

    manifold mf_default(*tracer, 0, 0);
    mf_default.configure(eps, 1e-8, 1e-14, 50, 1e-14, 100);
    mf_default.xPoint = {R_xp, Z_xp};
    std::vector<point> seg_default = mf_default.primarySegment(2);

    manifold mf_flipped(*tracer, 0, 0);
    mf_flipped.configure(eps, 1e-8, 1e-14, 50, 1e-14, 100);
    mf_flipped.set_branch(-1, -1);
    mf_flipped.xPoint = {R_xp, Z_xp};
    std::vector<point> seg_flipped = mf_flipped.primarySegment(2);

    // First point of default and flipped should be mirror images through xPoint.
    EXPECT_NEAR(seg_flipped[0].R, 2*R_xp - seg_default[0].R, 1e-10);
    EXPECT_NEAR(seg_flipped[0].Z, 2*Z_xp - seg_default[0].Z, 1e-10);
}

// Test: Find X-Point with valid guess
TEST_F(MfgenTest, Manifold_FindXPoint_ValidGuess) {
    manifold mf(*tracer, 0, 0);

    bool result = mf.find_xPoint(0.497999, -0.218603); // Example guess
    EXPECT_TRUE(result) << "Should find X-Point with valid guess";

    // Check that the found X-Point is within reasonable bounds
    // R: 0.4979691771716279 Z: -0.2185980054447758
    EXPECT_NEAR(mf.xPoint.R, 0.4979691771716279, 1e-6);
    EXPECT_NEAR(mf.xPoint.Z, -0.2185980054447758, 1e-6);
}

// Test: Primary segment computation
TEST_F(MfgenTest, Manifold_PrimarySegment) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    std::vector<point> segment;
    EXPECT_NO_THROW(segment = mf.primarySegment(10));

    // Check that the segment has the correct number of points
    EXPECT_EQ(segment.size(), 11);

    // Check first and last points are within reasonable bounds
    // Primary segment's first point:
    // R: 0.4979701320801332 Z: -0.2185977085445460
    // Primary segment's last point:
    // R: 0.4979855389268979 Z: -0.2185929033960707
    EXPECT_NEAR(segment.front().R, 0.4979701320801332, 1e-6);
    EXPECT_NEAR(segment.front().Z, -0.2185977085445460, 1e-6);
    EXPECT_NEAR(segment.back().R, 0.4979855389268979, 1e-6);
    EXPECT_NEAR(segment.back().Z, -0.2185929033960707, 1e-6);
}

// Test: backtracking in the refinement must not record a node twice
// (a deep segment used to contain exact-duplicate consecutive points).
TEST_F(MfgenTest, Manifold_NewSegment_NoDuplicatePoints) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;
    std::vector<point> seg = mf.primarySegment(10);

    for (int n = 1; n <= 4; ++n)
        seg = mf.newSegment(seg, 2e-2, THETA_LIM);

    ASSERT_GT(seg.size(), 100u); // deep enough that refinement backtracks
    for (size_t i = 1; i < seg.size(); ++i)
        EXPECT_FALSE(seg[i] == seg[i - 1]) << "duplicate at index " << i;
}

// Test: with no insertion budget the triples are accepted as they are
// (the segment is not truncated): one image per input node.
TEST_F(MfgenTest, Manifold_NewSegment_ZeroInsertionBudgetKeepsAllNodes) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 0);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;
    std::vector<point> primary = mf.primarySegment(10);

    std::vector<point> seg;
    EXPECT_NO_THROW(seg = mf.newSegment(primary, 1e-9, THETA_LIM));
    EXPECT_EQ(seg.size(), primary.size());
    EXPECT_NO_THROW(seg = mf.newSegment(primary, 1, 1e-9, THETA_LIM));
    EXPECT_EQ(seg.size(), primary.size());
}

// Test: a base point that is not a saddle (complex multipliers) is rejected
TEST_F(MfgenTest, Manifold_PrimarySegment_ThrowsAtNonHyperbolicPoint) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = MAG_AXIS_R; // inside the last closed surface: elliptic
    mf.xPoint.Z = MAG_AXIS_Z;
    EXPECT_THROW(mf.primarySegment(10), std::runtime_error);
}

// Test: too few points to refine is an error, not an empty result
TEST_F(MfgenTest, Manifold_NewSegment_ThrowsOnTooFewPoints) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    std::vector<point> two = {{0.4979701, -0.2185977}, {0.4979716, -0.2185972}};
    EXPECT_THROW(mf.newSegment(two, 0.005, THETA_LIM), std::invalid_argument);
    EXPECT_THROW(mf.newSegment(two, 1, 0.005, THETA_LIM), std::invalid_argument);
}

// Distance from p to the polyline through pts.
static double dist_to_polyline(const point &p, const std::vector<point> &pts) {
    double best = 1e300;
    for (size_t k = 0; k + 1 < pts.size(); ++k) {
        point a = pts[k], ab = pts[k + 1] - pts[k];
        double L2 = ab.R * ab.R + ab.Z * ab.Z;
        double u = L2 > 0 ? ((p.R - a.R) * ab.R + (p.Z - a.Z) * ab.Z) / L2 : 0.0;
        u = std::clamp(u, 0.0, 1.0);
        point c = a + ab * u;
        best = std::min(best, std::hypot(p.R - c.R, p.Z - c.Z));
    }
    return best;
}

// Test: interpolant method vs exact-map method (Ciro et al., Chaos 28, 093106,
// Sec. VI C). The exact nodes must lie within the interpolation baseline
// H_c = theta_c * l_c / 8 (Eq. 40) of the interpolant curve, and well below it
// when the curve is sampled through the arcs (sensitive to the arc sign).
TEST_F(MfgenTest, Manifold_Interpolant_MatchesExactMethodWithinBaseline) {
    const double l_lim = 5e-3;
    manifold mf0(*tracer, 0, 0), mf1(*tracer, 0, 0);
    for (manifold *m : {&mf0, &mf1}) {
        m->configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
        m->xPoint.R = 0.4979691771716279;
        m->xPoint.Z = -0.2185980054447758;
    }
    std::vector<point> primary0 = mf0.primarySegment(10);
    std::vector<point> exact = mf0.newSegment(primary0, 3, l_lim, THETA_LIM);

    std::vector<point> approx = mf1.primarySegment(10);
    for (int n = 1; n <= 3; ++n)
        approx = mf1.newSegment(approx, l_lim, THETA_LIM);

    std::vector<point> curve; // dense sampling of the interpolant arcs
    for (const auto &arc : manifold::buildInterpolants(approx))
        for (int k = 0; k < 16; ++k)
            curve.push_back(arc.evalNewPoint(k / 16.0));
    curve.push_back(approx.back());

    double max_chord = 0.0, max_arc = 0.0;
    for (const auto &x : exact) {
        max_chord = std::max(max_chord, dist_to_polyline(x, approx));
        max_arc = std::max(max_arc, dist_to_polyline(x, curve));
    }
    const double H_c = THETA_LIM * l_lim / 8.0;
    std::cout << "max dist exact->chords " << max_chord << ", ->arcs " << max_arc
              << ", H_c " << H_c << "\n";
    EXPECT_LT(max_chord, H_c);
    EXPECT_LT(max_arc, 0.05 * H_c); // 2.8e-6 measured; a wrong-sign arc gives 1.6e-5
}

// Test: a failed map evaluation is reported, not silently truncated
TEST_F(MfgenTest, Manifold_NewSegment_ThrowsOnFailedMap) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    // Points far outside the source domain: field-line tracing fails.
    std::vector<point> bad_seg = {{50.0, 50.0}, {50.1, 50.0}, {50.2, 50.0}};
    EXPECT_THROW(mf.newSegment(bad_seg, 0.005, THETA_LIM), std::runtime_error);
    EXPECT_THROW(mf.newSegment(bad_seg, 1, 0.005, THETA_LIM), std::runtime_error);
}

// Test: New segment computation from previous segment
TEST_F(MfgenTest, Manifold_NewSegment_FromPrevious) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    std::vector<point> primary_seg = {
        {0.4979701320801332, -0.218597708544546},
        {0.4979716727648097, -0.2185972280296985},
        {0.4979732134494862, -0.2185967475148509},
        {0.4979747541341626, -0.2185962670000034},
        {0.4979762948188391, -0.2185957864851559},
        {0.4979778355035156, -0.2185953059703084},
        {0.497979376188192, -0.2185948254554608},
        {0.4979809168728685, -0.2185943449406133},
        {0.497982457557545, -0.2185938644257657},
        {0.4979839982422214, -0.2185933839109182},
        {0.4979855389268979, -0.2185929033960707}};

    std::vector<point> new_seg_1, new_seg_2, new_seg_3;
    EXPECT_NO_THROW(new_seg_1 = mf.newSegment(primary_seg, 0.005, THETA_LIM));
    EXPECT_NO_THROW(new_seg_2 = mf.newSegment(new_seg_1, 0.005, THETA_LIM));
    EXPECT_NO_THROW(new_seg_3 = mf.newSegment(new_seg_2, 0.005, THETA_LIM));

    // Reference values after Stage 3 fixes (Bugs 1, 2, 5):
    //   - Bug 2: new_seg now starts empty; first 3 entries are actual map images,
    //     not raw points from the previous segment generation.
    //   - Bug 5: buildInterpolants now starts at i=0, covering the first interval.
    //   - Interpolant arcs use signed turning angles (bulge on the correct side for
    //     both turning directions): indices 9-34 moved by ~1e-5 m.
    //   Indices 4-34 are now shifted along the manifold relative to the pre-fix
    //   values, with monotonically increasing R and decreasing Z — physically correct.
    EXPECT_NEAR(new_seg_3[0].R,  0.50371942366068234, 1e-6);
    EXPECT_NEAR(new_seg_3[0].Z,  -0.21671800599570651, 1e-6);
    EXPECT_NEAR(new_seg_3[4].R,  0.51779070168116548, 1e-6);
    EXPECT_NEAR(new_seg_3[4].Z,  -0.21256465989197604, 1e-6);
    EXPECT_NEAR(new_seg_3[9].R,  0.5329413313161456, 1e-6);
    EXPECT_NEAR(new_seg_3[9].Z,  -0.20684384922194668, 1e-6);
    EXPECT_NEAR(new_seg_3[14].R, 0.544888620002773, 1e-6);
    EXPECT_NEAR(new_seg_3[14].Z, -0.19789861738018513, 1e-6);
    EXPECT_NEAR(new_seg_3[19].R, 0.5623249055389825, 1e-6);
    EXPECT_NEAR(new_seg_3[19].Z, -0.19188669144018536, 1e-6);
    EXPECT_NEAR(new_seg_3[24].R, 0.5773761419825298, 1e-6);
    EXPECT_NEAR(new_seg_3[24].Z, -0.18797022952026268, 1e-6);
    EXPECT_NEAR(new_seg_3[29].R, 0.5936666484608587,  1e-6);
    EXPECT_NEAR(new_seg_3[29].Z, -0.1813653812146642,  1e-6);
    EXPECT_NEAR(new_seg_3[34].R, 0.6134021821972827, 1e-6);
    EXPECT_NEAR(new_seg_3[34].Z, -0.17094516309055374, 1e-6);
}

// Test: New segment computation from primary segment
TEST_F(MfgenTest, Manifold_NewSegment_FromPrimary) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    std::vector<point> primary_seg = {
        {0.4979701320801332, -0.218597708544546},
        {0.4979716727648097, -0.2185972280296985},
        {0.4979732134494862, -0.2185967475148509},
        {0.4979747541341626, -0.2185962670000034},
        {0.4979762948188391, -0.2185957864851559},
        {0.4979778355035156, -0.2185953059703084},
        {0.497979376188192, -0.2185948254554608},
        {0.4979809168728685, -0.2185943449406133},
        {0.497982457557545, -0.2185938644257657},
        {0.4979839982422214, -0.2185933839109182},
        {0.4979855389268979, -0.2185929033960707}};

    std::vector<point> new_seg_3;
    EXPECT_NO_THROW(new_seg_3 = mf.newSegment(primary_seg, 3, 0.005, THETA_LIM));

    // Check a few new_seg_3 points within reasonable bounds
    // R: 0.5037194321479069 Z: -0.2167180039510279 0
    // R: 0.5177907149090401 Z: -0.2125646567853695 4
    // R: 0.5328067328393506 Z: -0.2069275324563348 9
    // R: 0.5448185439171523 Z: -0.1979494329010073 14
    // R: 0.562215480035554  Z: -0.1919076349028341 19
    // R: 0.5771532912095396 Z: -0.1880448334551187 24
    // R: 0.5936666582798458 Z: -0.1813653722927497 29
    // R: 0.6131264483728508 Z: -0.171105326569633 34
    EXPECT_NEAR(new_seg_3[0].R, 0.5037194321479069, 1e-6);
    EXPECT_NEAR(new_seg_3[0].Z, -0.2167180039510279, 1e-6);
    EXPECT_NEAR(new_seg_3[4].R, 0.5177907149090401, 1e-6);
    EXPECT_NEAR(new_seg_3[4].Z, -0.2125646567853695, 1e-6);
    EXPECT_NEAR(new_seg_3[9].R, 0.5328067328393506, 1e-6);
    EXPECT_NEAR(new_seg_3[9].Z, -0.2069275324563348, 1e-6);
    EXPECT_NEAR(new_seg_3[14].R, 0.5448185439171523, 1e-6);
    EXPECT_NEAR(new_seg_3[14].Z, -0.1979494329010073, 1e-6);
    EXPECT_NEAR(new_seg_3[19].R, 0.562215480035554, 1e-6);
    EXPECT_NEAR(new_seg_3[19].Z, -0.1919076349028341, 1e-6);
    EXPECT_NEAR(new_seg_3[24].R, 0.5771532912095396, 1e-6);
    EXPECT_NEAR(new_seg_3[24].Z, -0.1880448334551187, 1e-6);
    EXPECT_NEAR(new_seg_3[29].R, 0.5936666582798458, 1e-6);
    EXPECT_NEAR(new_seg_3[29].Z, -0.1813653722927497, 1e-6);
    EXPECT_NEAR(new_seg_3[34].R, 0.6131264483728508, 1e-6);
    EXPECT_NEAR(new_seg_3[34].Z, -0.171105326569633, 1e-6);
}

// Test: outputData is empty before any computation
TEST_F(MfgenTest, Manifold_OutputData_Empty) {
    manifold mf(*tracer, 0, 0);
    EXPECT_TRUE(mf.get_output_data().empty());
}

// Test: outputData accumulates after primarySegment and newSegment calls
TEST_F(MfgenTest, Manifold_OutputData_Accumulates) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    std::vector<point> primary = mf.primarySegment(10);
    EXPECT_EQ(mf.get_output_data().size(), 1u);
    EXPECT_EQ(mf.get_output_data()[0].size(), 11u);

    mf.newSegment(primary, 0.005, THETA_LIM); // interpolant method
    EXPECT_EQ(mf.get_output_data().size(), 2u);

    mf.newSegment(primary, 2, 0.005, THETA_LIM); // exact-map method
    EXPECT_EQ(mf.get_output_data().size(), 3u);
}

// Test: save() writes a .dat file with correct header and structure
TEST_F(MfgenTest, Manifold_Save_Dat) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;
    mf.primarySegment(10);

    std::string path = test_dir + "/manifold_test.dat";
    EXPECT_TRUE(mf.save(path));

    std::ifstream f(path);
    ASSERT_TRUE(f.is_open());

    // First line is a comment header
    std::string header;
    std::getline(f, header);
    EXPECT_EQ(header[0], '#');

    // First data line: segment index 0
    std::string first_data;
    std::getline(f, first_data);
    EXPECT_EQ(first_data[0], '0');
}

// Test: save() writes a .csv file with correct header
TEST_F(MfgenTest, Manifold_Save_Csv) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    std::vector<point> primary = mf.primarySegment(10);
    mf.newSegment(primary, 0.005, THETA_LIM);

    std::string path = test_dir + "/manifold_test.csv";
    EXPECT_TRUE(mf.save(path));

    std::ifstream f(path);
    ASSERT_TRUE(f.is_open());

    std::string header;
    std::getline(f, header);
    EXPECT_EQ(header, "seg,R,Z");

    // First data row: segment 0
    std::string row0;
    std::getline(f, row0);
    EXPECT_EQ(row0[0], '0');
    EXPECT_NE(row0.find(','), std::string::npos);

    // Skip the remaining 10 rows of seg 0 (primary has 11 pts; 1 already read)
    for (int i = 0; i < 10; ++i) std::getline(f, row0);
    std::string row1;
    std::getline(f, row1);
    EXPECT_EQ(row1[0], '1');
}

// Test: save() returns false for unsupported extension
TEST_F(MfgenTest, Manifold_Save_UnsupportedFormat) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;
    mf.primarySegment(10);

    EXPECT_FALSE(mf.save(test_dir + "/manifold_test.xyz"));
}

// Test: run() produces the correct number of segments in outputData
TEST_F(MfgenTest, Manifold_Run_OutputDataSize) {
    manifold mf(*tracer, 0, 0);
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf.xPoint.R = 0.4979691771716279;
    mf.xPoint.Z = -0.2185980054447758;

    mf.run(10, 4, 1, 0.005, THETA_LIM); // 4 total: 1 primary + 3 new (interpolant)
    EXPECT_EQ(mf.get_output_data().size(), 4u);
    EXPECT_EQ(mf.get_output_data()[0].size(), 11u); // primary: n_intervals=10 → 11 pts
}

// Test: run() (exact-map) first new segment matches newSegment called directly
TEST_F(MfgenTest, Manifold_Run_ExactMap_MatchesDirect) {
    // Reference: compute seg 1 via run(method=0)
    manifold mf1(*tracer, 0, 0);
    mf1.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf1.xPoint.R = 0.4979691771716279;
    mf1.xPoint.Z = -0.2185980054447758;
    mf1.run(10, 2, 0, 0.005, THETA_LIM); // 1 primary + 1 exact-map new

    // Direct: compute seg 1 via primarySegment + newSegment
    manifold mf2(*tracer, 0, 0);
    mf2.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100);
    mf2.xPoint.R = 0.4979691771716279;
    mf2.xPoint.Z = -0.2185980054447758;
    std::vector<point> primary = mf2.primarySegment(10);
    std::vector<point> seg1    = mf2.newSegment(primary, 1, 0.005, THETA_LIM);

    ASSERT_EQ(mf1.get_output_data()[1].size(), seg1.size());
    for (size_t i = 0; i < seg1.size(); ++i) {
        EXPECT_DOUBLE_EQ(mf1.get_output_data()[1][i].R, seg1[i].R);
        EXPECT_DOUBLE_EQ(mf1.get_output_data()[1][i].Z, seg1[i].Z);
    }
}

// ==================== INTERPOLANT GEOMETRY (no field source needed) ====================

// Points on the unit circle at constant angular step; sgn=+1 counter-clockwise, -1 clockwise.
static std::vector<point> circle_points(double sgn, double dtheta, int n) {
    std::vector<point> P; // shape=(n,2), cols=[R,Z], units=m
    for (int k = 0; k < n; ++k)
        P.push_back({std::cos(sgn * k * dtheta), std::sin(sgn * k * dtheta)});
    return P;
}

// The interpolant midpoint of an interior interval must lie on the circle for
// either turning direction (error well below the chord midpoint's sagitta ~1e-2).
TEST(ManifoldGeometry, Interpolant_Midpoint_OnCircle_BothTurningDirections) {
    const double dtheta = 0.3; // rad
    for (double sgn : {+1.0, -1.0}) {
        auto P    = circle_points(sgn, dtheta, 4);
        auto arcs = manifold::buildInterpolants(P);
        ASSERT_EQ(arcs.size(), 3u);
        point m    = arcs[1].evalNewPoint(0.5);
        point want = {std::cos(sgn * 1.5 * dtheta), std::sin(sgn * 1.5 * dtheta)};
        double err = std::hypot(m.R - want.R, m.Z - want.Z);
        EXPECT_LT(err, 1e-3) << "turning direction sgn=" << sgn;
    }
}

// Straight line: no curvature, midpoint is the chord midpoint.
TEST(ManifoldGeometry, Interpolant_Midpoint_StraightLine) {
    std::vector<point> P = {{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}};
    auto arcs = manifold::buildInterpolants(P);
    point m = arcs[1].evalNewPoint(0.5);
    EXPECT_NEAR(m.R, 1.5, 1e-14);
    EXPECT_NEAR(m.Z, 0.0, 1e-14);
}

// The last interval gets an arc too (it must be refinable).
TEST(ManifoldGeometry, Interpolant_ArcForEveryInterval) {
    auto P = circle_points(+1.0, 0.3, 5);
    EXPECT_EQ(manifold::buildInterpolants(P).size(), P.size() - 1);
}

// A near-reversal (turning angle close to pi) must not send the midpoint off to
// infinity: it stays within one chord length of the chord midpoint.
TEST(ManifoldGeometry, Interpolant_Midpoint_BoundedAtSharpTurn) {
    std::vector<point> P = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1e-6}, {1.0, 1e-6}};
    auto arcs = manifold::buildInterpolants(P);
    for (const auto &arc : arcs) {
        point m    = arc.evalNewPoint(0.5);
        point base = (arc.x0 + arc.x1) * 0.5;
        double chord = std::hypot(arc.x1.R - arc.x0.R, arc.x1.Z - arc.x0.Z);
        ASSERT_TRUE(std::isfinite(m.R) && std::isfinite(m.Z));
        EXPECT_LT(std::hypot(m.R - base.R, m.Z - base.Z), chord);
    }
}

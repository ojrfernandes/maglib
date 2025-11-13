#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Include the mfgen headers
#include "../mfgen/src/input_read.h"
#include "../mfgen/src/manifold.h"

// Test fixture for mfgen tests
class MfgenTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::string source_path = std::string(TEST_DATA_DIR) + "/C1.h5";
        tracer = new maglit(source_path.c_str(), FIO_M3DC1_SOURCE, 1);
        tracer->configure(1e-2, 1e-6, 1e-2);
    }

    static void TearDownTestSuite() {
        delete tracer;
        tracer = nullptr;
    }

    static maglit *tracer;
    static std::string source_path;

    void SetUp() override {
        // Create temporary directory for test files
        test_dir = "test_input_files";
        std::filesystem::create_directories(test_dir);

        // Create valid test input file
        valid_input_file = test_dir + "/valid_input.txt";
        createValidInputFile();

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
            << "        source_path = test_source.h5       # M3DC1 file\n"
            << "        output_path = test_output.dat      # output file path and name\n"
            << "#\n"
            << "#=============== TRACING PARAMETERS\n"
            << "#\n"
            << "        timeslice = 1           # M3DC1 timeslice \n"
            << "        manifold  = 0           # stable manifold=0; unstable manifold=1\n"
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

    void createInvalidInputFiles() {
        // File with missing parameters
        std::ofstream invalid_file(invalid_input_file);
        invalid_file << "# Missing required parameters\n"
                     << "source_path = test.h5\n"
                     << "# Missing output_path and numeric parameters\n";
        invalid_file.close();

        // File with malformed syntax
        std::ofstream malformed_file(malformed_input_file);
        malformed_file << "source_path test.h5  # Missing equals sign\n"
                       << "output_path =         # Missing value\n"
                       << "num_threads = abc    # Invalid numeric value\n"
                       << "gridMin = 1.2.3      # Invalid float\n";
        malformed_file.close();

        // Empty file
        std::ofstream empty_file(empty_input_file);
        empty_file.close();
    }

    void createEdgeCaseInputFile() {
        // Edge case values
        std::ofstream edge_case_input(edge_case_input_file);
        edge_case_input
            << "source_path = \n"      // Empty string
            << "output_path = a\n"     // Single character
            << "timeslice = 0\n"       // zero value
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
            << "\n" // Empty line
            << "   # Indented comment\n"
            << "\t\tsource_path = test.h5\t\t# Trailing comment with tabs\n"
            << "output_path=test.dat# No spaces around equals\n"
            << "timeslice = 1 # Comment\n"
            << "\tmanifold\t=\t1\t\n" // Tabs instead of spaces
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
    input_read reader(test_path);

    // test that the constructor doesn't crash
    EXPECT_NO_THROW(input_read reader2(test_path));
}

// Test: Valid input file reading
TEST_F(MfgenTest, InputRead_ValidFile) {
    input_read reader(valid_input_file);
    bool result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should successfully read valid input file";

    // Check string parameters
    EXPECT_EQ(reader.source_path, "test_source.h5");
    EXPECT_EQ(reader.output_path, "test_output.dat");

    // Check numeric parameters
    EXPECT_EQ(reader.timeslice, 1);
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
}

// Test: Invalid file handling
TEST_F(MfgenTest, InputRead_InvalidFile) {
    input_read reader(invalid_input_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for invalid file";
}

// Test: Nonexistent file handling
TEST_F(MfgenTest, InputRead_NonexistentFile) {
    input_read reader(nonexistent_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for nonexistent file";
}

// Test: Malformed file handling
TEST_F(MfgenTest, InputRead_MalformedFile) {
    input_read reader(malformed_input_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for malformed file";
}

// Test: Empty file handling
TEST_F(MfgenTest, InputRead_EmptyFile) {
    input_read reader(empty_input_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for empty file";
}

// Test: Edge case values
TEST_F(MfgenTest, InputRead_EdgeCaseValues) {
    input_read reader(edge_case_input_file);
    bool result = reader.readInputFile();

    if (result) {
        // Check string parameters
        EXPECT_EQ(reader.source_path, "");
        EXPECT_EQ(reader.output_path, "a");

        // Check numeric parameters
        EXPECT_EQ(reader.timeslice, 0);
        EXPECT_EQ(reader.manifold, -1);
        EXPECT_EQ(reader.method, 9999);
        EXPECT_DOUBLE_EQ(reader.Phi, -360.0);
        EXPECT_EQ(reader.nSections, 0);
        EXPECT_DOUBLE_EQ(reader.phi_0, 360.0);
        EXPECT_DOUBLE_EQ(reader.phi_1, -360.0);
        EXPECT_DOUBLE_EQ(reader.epsilon, -1e-6);
        EXPECT_EQ(reader.nSegments, 0);
        EXPECT_DOUBLE_EQ(reader.l_lim, -0.01);
        EXPECT_DOUBLE_EQ(reader.theta_lim, -10.0);
        EXPECT_DOUBLE_EQ(reader.h_init, -1e-2);
        EXPECT_DOUBLE_EQ(reader.h_min, -1e-6);
        EXPECT_DOUBLE_EQ(reader.h_max, -1e-2);
        EXPECT_DOUBLE_EQ(reader.h_deriv, -1e-3);
        EXPECT_DOUBLE_EQ(reader.n_tol, -10);
        EXPECT_EQ(reader.max_iter, -100);
        EXPECT_DOUBLE_EQ(reader.precision, -1e-8);
        EXPECT_EQ(reader.max_insertions, -5);
        EXPECT_EQ(reader.verbose, 1);
    }

    EXPECT_NO_THROW(reader.readInputFile());
}

// Test: Comments and whitespace handling
TEST_F(MfgenTest, InputRead_CommentsAndWhitespace) {
    input_read reader(comment_whitespace_input_file);
    bool result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should handle comments and whitespace correctly";

    if (result) {
        EXPECT_EQ(reader.source_path, "test.h5");
        EXPECT_EQ(reader.output_path, "test.dat");
        EXPECT_EQ(reader.timeslice, 1);
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
    std::string first_source = reader.source_path;
    int first_segments = reader.nSegments;

    // Second read should give same results
    bool result2 = reader.readInputFile();
    EXPECT_TRUE(result2);
    EXPECT_EQ(reader.source_path, first_source);
    EXPECT_EQ(reader.nSegments, first_segments);
}

// ==================== MANIFOLD TESTS ====================/

// Static member definitions
maglit *MfgenTest::tracer = nullptr;

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
    EXPECT_NO_THROW(mf.primarySegment(segment, 10));

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

    std::vector<point> new_seg_1;
    std::vector<point> new_seg_2;
    std::vector<point> new_seg_3;
    EXPECT_NO_THROW(mf.newSegment(primary_seg, new_seg_1, 0, 0.005, 20));
    EXPECT_NO_THROW(mf.newSegment(new_seg_1, new_seg_2, 0, 0.005, 20));
    EXPECT_NO_THROW(mf.newSegment(new_seg_2, new_seg_3, 0, 0.005, 20));

    // Check a few new_seg_3 points within reasonable bounds
    // R: 0.5012464706502761 Z: -0.2174148574827388 0
    // R: 0.5030591235329525 Z: -0.2168648161175224 4
    // R: 0.5133430437264513 Z: -0.2130464750286543 9
    // R: 0.5298324820702656 Z: -0.208634750788582  14
    // R: 0.5429352654656256 Z: -0.199387939879125  19
    // R: 0.5594125465315225 Z: -0.192432899160264  24
    // R: 0.5746497663876523 Z: -0.1888496529292251 29
    // R: 0.5901009484232285 Z: -0.1829887187256312 34
    // R: 0.6093889559868346 Z: -0.1732392616766928 39
    // R: 0.6289407305354183 Z: -0.1613738267302098 44
    // R: 0.6461994779238116 Z: -0.1496757085334106 49
    EXPECT_NEAR(new_seg_3[0].R, 0.5012464706502761, 1e-6);
    EXPECT_NEAR(new_seg_3[0].Z, -0.2174148574827388, 1e-6);
    EXPECT_NEAR(new_seg_3[4].R, 0.5030591235329525, 1e-6);
    EXPECT_NEAR(new_seg_3[4].Z, -0.2168648161175224, 1e-6);
    EXPECT_NEAR(new_seg_3[9].R, 0.5133430437264513, 1e-6);
    EXPECT_NEAR(new_seg_3[9].Z, -0.2130464750286543, 1e-6);
    EXPECT_NEAR(new_seg_3[14].R, 0.5298324820702656, 1e-6);
    EXPECT_NEAR(new_seg_3[14].Z, -0.208634750788582, 1e-6);
    EXPECT_NEAR(new_seg_3[19].R, 0.5429352654656256, 1e-6);
    EXPECT_NEAR(new_seg_3[19].Z, -0.199387939879125, 1e-6);
    EXPECT_NEAR(new_seg_3[24].R, 0.5594125465315225, 1e-6);
    EXPECT_NEAR(new_seg_3[24].Z, -0.192432899160264, 1e-6);
    EXPECT_NEAR(new_seg_3[29].R, 0.5746497663876523, 1e-6);
    EXPECT_NEAR(new_seg_3[29].Z, -0.1888496529292251, 1e-6);
    EXPECT_NEAR(new_seg_3[34].R, 0.5901009484232285, 1e-6);
    EXPECT_NEAR(new_seg_3[34].Z, -0.1829887187256312, 1e-6);
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
    EXPECT_NO_THROW(mf.newSegment(primary_seg, new_seg_3, 0, 3, 0.005, 20));

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
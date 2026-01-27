#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Include the fpgen headers
#include "../fpgen/src/footprint.h"
#include "../fpgen/src/input_read.h"

// Test fixture for fpgen tests
class FpgenTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::string source_path = std::string(TEST_DATA_DIR) + "/C1.h5";
        std::string shape_path = std::string(TEST_DATA_DIR) + "/tcabr_first_wall.txt";
        tracer = new maglit(source_path.c_str(), FIO_M3DC1_SOURCE, 1);
        tracer->configure(0.01, 1e-6, 0.1);
        tracer->set_monitor(shape_path);
    }

    static void TearDownTestSuite() {
        delete tracer;
        tracer = nullptr;
    }

    static maglit     *tracer;
    static std::string source_path;
    static std::string shape_path;

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
            << "################################ FPGEN #################################\n"
            << "#\n"
            << "#=============== I/O FILES\n"
            << "        source_path = test_source.h5       # M3DC1 file\n"
            << "        shape_path = test_shape.txt        # machine boundary shape file\n"
            << "        output_path = test_output.dat      # output file path and name\n"
            << "#\n"
            << "#=============== MAPPING PARAMETERS        \n"
            << "#       \n"
            << "        timeslice = 1           # M3DC1 timeslice \n"
            << "        manifold  = 1           # stable manifold=0; unstable manifold=1\n"
            << "\n"
            << "        grid_R1 = 0.435         # first point (R,Z) delimiting the target plate mapped surface\n"
            << "        grid_Z1 = -0.239        # first point (R,Z) delimiting the target plate mapped surface\n"
            << "\n"
            << "        grid_R2 = 0.435         # second point (R,Z) delimiting the target plate mapped surface\n"
            << "        grid_Z2 = -0.232        # second point (R,Z) delimiting the target plate mapped surface\n"
            << "\n"
            << "        nRZ = 10                # grid dimension along the (R,Z) plane\n"
            << "        nPhi = 20               # grid dimension along toroidal direction (Phi)\n"
            << "#\n"
            << "#=============== ADDITIONAL PARAMETERS\n"
            << "#\n"
            << "        num_threads = 1         # number of threads for OpenMP execution\n"
            << "        max_turns = 10000       # maximum toroidal turns for field line integration\n"
            << "        h_init = 1e-2           # initial step-size for integration\n"
            << "        h_min = 1e-6            # minimum step-size for integration\n"
            << "        h_max = 1e-2            # maximum step-size for integration\n";
        valid_input.close();
    }

    void createInvalidInputFiles() {
        // File with missing parameters
        std::ofstream invalid_input(invalid_input_file);
        invalid_input << "# Missing required parameters\n"
                      << "source_path = test.h5\n"
                      << "# Missing shape_path, output_path, and numeric parameters\n";
        invalid_input.close();

        // File with malformed syntax
        std::ofstream malformed_input(malformed_input_file);
        malformed_input << "source_path test.h5  # Missing equals sign\n"
                        << "shape_path =         # Missing value\n"
                        << "num_threads = abc    # Invalid numeric value\n"
                        << "gridMin = 1.2.3      # Invalid float\n";
        malformed_input.close();

        // Empty file
        std::ofstream empty_input(empty_input_file);
        empty_input.close();
    }

    void createEdgeCaseInputFile() {
        // Edge case values
        std::ofstream edge_case_input(edge_case_input_file);
        edge_case_input
            << "source_path = \n" // Empty string
            << "shape_path = a\n" // Single character
            << "output_path = very_long_filename_with_many_characters_to_test_string_handling.dat\n"
            << "timeslice = 999\n" // Large value
            << "manifold = 0\n"    // Lower bound
            << "grid_R1 = -1.0\n"  // Negative float
            << "grid_Z1 = 1e6\n"   // Very large float
            << "grid_R2 = 0.0\n"   // Zero value
            << "grid_Z2 = -1e6\n"  // Very large negative float
            << "nRZ = 1\n"         // Minimum grid size
            << "nPhi = 1000000\n"  // Very large grid size
            << "num_threads = 0\n" // Edge case for threads
            << "max_turns = 1\n"   // Minimum turns
            << "h_init = 1e-10\n"  // Very small step size
            << "h_min = 0.0\n"     // Zero step size
            << "h_max = 1e10\n";   // Very large step size

        edge_case_input.close();

        // Comments and whitespace variations
        std::ofstream comment_whitespace_input(comment_whitespace_input_file);
        comment_whitespace_input
            << "# This is a comment line\n"
            << "\n" // Empty line
            << "   # Indented comment\n"
            << "\t\tsource_path = test.h5\t\t# Trailing comment with tabs\n"
            << "    shape_path    =    test.txt    # Lots of spaces\n"
            << "output_path=test.dat# No spaces around equals\n"
            << "timeslice = 1 # Comment\n"
            << "\tmanifold\t=\t1\t\n" // Tabs instead of spaces
            << "grid_R1 = 0.435\n"
            << "grid_Z1 = -0.239\n"
            << "grid_R2 = 0.435\n"
            << "grid_Z2 = -0.232\n"
            << "nRZ = 10\n"
            << "nPhi = 20\n"
            << "num_threads = 8\n" // Trailing spaces
            << "max_turns = 10000\n"
            << "h_init = 1e-2\n"
            << "h_min = 1e-6\n"
            << "h_max = 1e-2\n";

        comment_whitespace_input.close();
    }

    std::string test_dir;
    std::string valid_input_file;
    std::string invalid_input_file;
    std::string malformed_input_file;
    std::string empty_input_file;
    std::string edge_case_input_file;
    std::string comment_whitespace_input_file;
    std::string nonexistent_input_file = "nonexistent_file.txt";
};

// ==================== INPUT_READ TESTS ====================/

// Test: Constructor initializes reading_path correctly
TEST_F(FpgenTest, InputRead_Constructor) {
    std::string test_path = "test_path.txt";
    input_read  reader(test_path);

    // test that the constructor doesn't crash
    EXPECT_NO_THROW(input_read reader2(test_path));
}

// Test: Valid input file parsing
TEST_F(FpgenTest, InputRead_ValidFile) {
    input_read reader(valid_input_file);
    bool       result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should successfully read valid input file";

    // Check string parameters
    EXPECT_EQ(reader.source_path, "test_source.h5");
    EXPECT_EQ(reader.shape_path, "test_shape.txt");
    EXPECT_EQ(reader.output_path, "test_output.dat");

    // Check numeric parameters
    EXPECT_EQ(reader.timeslice, 1);
    EXPECT_EQ(reader.manifold, 1);
    EXPECT_DOUBLE_EQ(reader.grid_R1, 0.435);
    EXPECT_DOUBLE_EQ(reader.grid_Z1, -0.239);
    EXPECT_DOUBLE_EQ(reader.grid_R2, 0.435);
    EXPECT_DOUBLE_EQ(reader.grid_Z2, -0.232);
    EXPECT_EQ(reader.nRZ, 10);
    EXPECT_EQ(reader.nPhi, 20);
    EXPECT_EQ(reader.num_threads, 1);
    EXPECT_EQ(reader.max_turns, 10000);
    EXPECT_DOUBLE_EQ(reader.h_init, 1e-2);
    EXPECT_DOUBLE_EQ(reader.h_min, 1e-6);
    EXPECT_DOUBLE_EQ(reader.h_max, 1e-2);
}

// Test: Invalid file handling
TEST_F(FpgenTest, InputRead_InvalidFile) {
    input_read reader(invalid_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for invalid file";
}

// Test: Nonexistent file handling
TEST_F(FpgenTest, InputRead_NonexistentFile) {
    input_read reader(nonexistent_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for nonexistent file";
}

// Test: Malformed file handling
TEST_F(FpgenTest, InputRead_MalformedSyntax) {
    input_read reader(malformed_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for malformed input file";
}

// Test: Empty file handling
TEST_F(FpgenTest, InputRead_EmptyFile) {
    input_read reader(empty_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for empty file";
}

// Test: Edge case values
TEST_F(FpgenTest, InputRead_EdgeCaseValues) {
    input_read reader(edge_case_input_file);
    bool       result = reader.readInputFile();

    if (result) {
        EXPECT_EQ(reader.source_path, "");
        EXPECT_EQ(reader.shape_path, "a");
        EXPECT_EQ(reader.output_path,
                  "very_long_filename_with_many_characters_to_test_string_handling.dat");
        EXPECT_EQ(reader.timeslice, 999);
        EXPECT_EQ(reader.manifold, 0);
        EXPECT_DOUBLE_EQ(reader.grid_R1, -1.0);
        EXPECT_DOUBLE_EQ(reader.grid_Z1, 1e6);
        EXPECT_DOUBLE_EQ(reader.grid_R2, 0.0);
        EXPECT_DOUBLE_EQ(reader.grid_Z2, -1e6);
        EXPECT_EQ(reader.nRZ, 1);
        EXPECT_EQ(reader.nPhi, 1000000);
        EXPECT_EQ(reader.num_threads, 0);
        EXPECT_EQ(reader.max_turns, 1);
        EXPECT_DOUBLE_EQ(reader.h_init, 1e-10);
        EXPECT_DOUBLE_EQ(reader.h_min, 0.0);
        EXPECT_DOUBLE_EQ(reader.h_max, 1e10);
    }

    EXPECT_NO_THROW(reader.readInputFile());
}

// Test: Comments and whitespace handling
TEST_F(FpgenTest, InputRead_CommentsAndWhitespace) {
    input_read reader(comment_whitespace_input_file);
    bool       result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should handle comments and whitespace correctly";

    if (result) {
        EXPECT_EQ(reader.source_path, "test.h5");
        EXPECT_EQ(reader.shape_path, "test.txt");
        EXPECT_EQ(reader.output_path, "test.dat");
        EXPECT_EQ(reader.timeslice, 1);
        EXPECT_EQ(reader.manifold, 1);
        EXPECT_EQ(reader.num_threads, 8);
    }
}

// Test: Multiple reads from same object
TEST_F(FpgenTest, InputRead_MultipleReads) {
    input_read reader(valid_input_file);

    // First read
    bool result1 = reader.readInputFile();
    EXPECT_TRUE(result1);

    // Store values from first read
    std::string first_source = reader.source_path;
    int         first_threads = reader.num_threads;

    // Second read should give same results
    bool result2 = reader.readInputFile();
    EXPECT_TRUE(result2);
    EXPECT_EQ(reader.source_path, first_source);
    EXPECT_EQ(reader.num_threads, first_threads);
}

// ==================== FOOTPRINT TESTS ====================/

// Static member definitions
maglit *FpgenTest::tracer = nullptr;

// Test: Constructor and basic initialization
TEST_F(FpgenTest, Footprint_Constructor) {
    // Test floor plate
    EXPECT_NO_THROW(footprint fp(0, 0.435, -0.239, 0.435, -0.232, 10, 20, 1000));

    // Test wall plate
    EXPECT_NO_THROW(footprint fp(1, 0.435, -0.239, 0.435, -0.232, 15, 30, 1000));
}

// Test: Output data structure initialization
TEST_F(FpgenTest, Footprint_OutputDataInitialization) {
    int       nRZ = 10;
    int       nPhi = 20;
    footprint fp(0, 0.5, 0.6, 0.4, 0.5, nRZ, nPhi, 1000);

    // Check that outputData is properly sized
    EXPECT_EQ(fp.outputData.size(), nRZ * nPhi);

    // Check that each row has 6 columns
    if (fp.outputData.size() > 0) {
        EXPECT_EQ(fp.outputData[0].size(), 6);
    }

    // Check all rows are properly initialized
    for (const auto &row : fp.outputData) {
        EXPECT_EQ(row.size(), 6);
    }
}

// Test: Simple 2 x 2 tcabr wall grid run limited to 50 turns
TEST_F(FpgenTest, Footprint_RunGrid_Simple2x2_Wall) {
    int       manifold = 0; // stable manifold
    double    grid_R1 = 0.435;
    double    grid_Z1 = -0.239;
    double    grid_R2 = 0.435;
    int       nR = 2;
    int       nZ = 2;
    int       max_turns = 50;
    double    grid_Z2 = -0.232;
    footprint fp(manifold, grid_R1, grid_Z1, grid_R2, grid_Z2, nR, nZ, max_turns);

    // Run the grid
    EXPECT_NO_THROW(fp.runGrid(*tracer));

    // Check that outputData is populated
    EXPECT_EQ(fp.outputData.size(), nR * nZ);

    // Check that each row has 6 columns
    if (fp.outputData.size() > 0) {
        EXPECT_EQ(fp.outputData[0].size(), 6);
    }

    // Check that data values (within 3 digit precision) are
    // 0.435000 -0.239000 0.000000 4.699479 0.995964 1
    // 0.435000 -0.232000 0.000000 10.802152 1.001697 3
    // 0.435000 -0.239000 3.141593 4.852588 0.995996 2
    // 0.435000 -0.232000 3.141593 10.739208 1.009234 3
    EXPECT_NEAR(fp.outputData[0][0], 0.435000, 1e-3);
    EXPECT_NEAR(fp.outputData[1][0], 0.435000, 1e-3);
    EXPECT_NEAR(fp.outputData[2][0], 0.435000, 1e-3);
    EXPECT_NEAR(fp.outputData[3][0], 0.435000, 1e-3);

    EXPECT_NEAR(fp.outputData[0][1], -0.239000, 1e-3);
    EXPECT_NEAR(fp.outputData[1][1], -0.232000, 1e-3);
    EXPECT_NEAR(fp.outputData[2][1], -0.239000, 1e-3);
    EXPECT_NEAR(fp.outputData[3][1], -0.232000, 1e-3);

    EXPECT_NEAR(fp.outputData[0][2], 0.000000, 1e-3);
    EXPECT_NEAR(fp.outputData[1][2], 0.000000, 1e-3);
    EXPECT_NEAR(fp.outputData[2][2], 3.141593, 1e-3);
    EXPECT_NEAR(fp.outputData[3][2], 3.141593, 1e-3);

    EXPECT_NEAR(fp.outputData[0][3], 4.699479, 1e-3);
    EXPECT_NEAR(fp.outputData[1][3], 10.802152, 1e-3);
    EXPECT_NEAR(fp.outputData[2][3], 4.852588, 1e-3);
    EXPECT_NEAR(fp.outputData[3][3], 10.739208, 1e-3);

    EXPECT_NEAR(fp.outputData[0][4], 0.995964, 1e-3);
    EXPECT_NEAR(fp.outputData[1][4], 1.001697, 1e-3);
    EXPECT_NEAR(fp.outputData[2][4], 0.995996, 1e-3);
    EXPECT_NEAR(fp.outputData[3][4], 1.009234, 1e-3);

    EXPECT_EQ(fp.outputData[0][5], 1);
    EXPECT_EQ(fp.outputData[1][5], 3);
    EXPECT_EQ(fp.outputData[2][5], 1);
    EXPECT_EQ(fp.outputData[3][5], 3);
}

// Test: Simple 2 x 2 tcabr floor grid run limited to 50 turns
TEST_F(FpgenTest, Footprint_RunGrid_Simple2x2_Floor) {
    int       manifold = 1; // unstable manifold
    double    grid_R1 = 0.51;
    double    grid_Z1 = -0.24;
    double    grid_R2 = 0.54;
    double    grid_Z2 = -0.24;
    int       nR = 2;
    int       nZ = 2;
    int       max_turns = 50;
    footprint fp(manifold, grid_R1, grid_Z1, grid_R2, grid_Z2, nR, nZ, max_turns);

    // Run the grid
    EXPECT_NO_THROW(fp.runGrid(*tracer));

    // Check that outputData is populated
    EXPECT_EQ(fp.outputData.size(), nR * nZ);

    // Check that each row has 6 columns
    if (fp.outputData.size() > 0) {
        EXPECT_EQ(fp.outputData[0].size(), 6);
    }

    // Check that data values (within 3 digit precision) are
    // 0.51 -0.24 0 10.5302 0.974279 3
    // 0.54 -0.24 0 8.5844 1.05251 2
    // 0.51 -0.24 3.14159 47.3727 0.882983 13
    // 0.54 -0.24 3.14159 8.80513 1.0541 2
    EXPECT_NEAR(fp.outputData[0][0], 0.51, 1e-3);
    EXPECT_NEAR(fp.outputData[1][0], 0.54, 1e-3);
    EXPECT_NEAR(fp.outputData[2][0], 0.51, 1e-3);
    EXPECT_NEAR(fp.outputData[3][0], 0.54, 1e-3);

    EXPECT_NEAR(fp.outputData[0][1], -0.24, 1e-3);
    EXPECT_NEAR(fp.outputData[1][1], -0.24, 1e-3);
    EXPECT_NEAR(fp.outputData[2][1], -0.24, 1e-3);
    EXPECT_NEAR(fp.outputData[3][1], -0.24, 1e-3);

    EXPECT_NEAR(fp.outputData[0][2], 0.000000, 1e-3);
    EXPECT_NEAR(fp.outputData[1][2], 0.000000, 1e-3);
    EXPECT_NEAR(fp.outputData[2][2], 3.141593, 1e-3);
    EXPECT_NEAR(fp.outputData[3][2], 3.141593, 1e-3);

    EXPECT_NEAR(fp.outputData[0][3], 10.5302, 1e-3);
    EXPECT_NEAR(fp.outputData[1][3], 8.5844, 1e-3);
    EXPECT_NEAR(fp.outputData[2][3], 47.3727, 1e-3);
    EXPECT_NEAR(fp.outputData[3][3], 8.80513, 1e-3);

    EXPECT_NEAR(fp.outputData[0][4], 0.974279, 1e-3);
    EXPECT_NEAR(fp.outputData[1][4], 1.05251, 1e-3);
    EXPECT_NEAR(fp.outputData[2][4], 0.882983, 1e-3);
    EXPECT_NEAR(fp.outputData[3][4], 1.0541, 1e-3);

    EXPECT_EQ(fp.outputData[0][5], 3);
    EXPECT_EQ(fp.outputData[1][5], 2);
    EXPECT_EQ(fp.outputData[2][5], 13);
    EXPECT_EQ(fp.outputData[3][5], 2);
}

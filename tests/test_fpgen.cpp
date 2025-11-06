#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Include the fpgen headers
#include "../fpgen/src/footprint.h"
#include "../fpgen/src/input_read.h"

// Test fixture for maglit functionality
class FpgenTest : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        std::string source_path = std::string(TEST_DATA_DIR) + "/C1.h5";
        tracer = new maglit(source_path.c_str(), FIO_M3DC1_SOURCE, 1);
        tracer->configure(0.01, 1e-6, 0.1);
    }

    static void TearDownTestSuite() {
        delete tracer;
        tracer = nullptr;
    }

    static maglit     *tracer;
    static std::string source_path;

    void SetUp() override {
        // Create temporary directory for test files
        test_dir = "test_input_files";
        std::filesystem::create_directories(test_dir);

        // Create valid test input file
        valid_input_file = test_dir + "/valid_input.txt";
        createValidInputFile(valid_input_file);

        // Create invalid test files
        invalid_input_file = test_dir + "/invalid_input.txt";
        malformed_input_file = test_dir + "/malformed_input.txt";
        empty_input_file = test_dir + "/empty_input.txt";

        createInvalidInputFiles();
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove_all(test_dir);
    }

    void createValidInputFile(const std::string &filename) {
        std::ofstream file(filename);
        file << "################################ FPGEN #################################\n"
             << "#\n"
             << "#=============== I/O FILES\n"
             << "        source_path = test_source.h5       # M3DC1 file\n"
             << "        shape_path = test_shape.txt        # machine boundary shape file\n"
             << "        output_path = test_output.dat      # output file path and name\n"
             << "#\n"
             << "#=============== MAPPING PARAMETERS        \n"
             << "#       \n"
             << "        timeslice = 1           # M3DC1 timeslice \n"
             << "        manifold  = 1           # unstable manifold=0; stable manifold=1\n"
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
        file.close();
    }

    void createInvalidInputFiles() {
        // File with missing parameters
        std::ofstream invalid_file(invalid_input_file);
        invalid_file << "# Missing required parameters\n"
                     << "source_path = test.h5\n"
                     << "# Missing shape_path, output_path, and numeric parameters\n";
        invalid_file.close();

        // File with malformed syntax
        std::ofstream malformed_file(malformed_input_file);
        malformed_file << "source_path test.h5  # Missing equals sign\n"
                       << "shape_path =         # Missing value\n"
                       << "num_threads = abc    # Invalid numeric value\n"
                       << "gridMin = 1.2.3      # Invalid float\n";
        malformed_file.close();

        // Empty file
        std::ofstream empty_file(empty_input_file);
        empty_file.close();
    }

    std::string test_dir;
    std::string valid_input_file;
    std::string invalid_input_file;
    std::string malformed_input_file;
    std::string empty_input_file;
    std::string nonexistent_file = "nonexistent_file.txt";
};

// ==================== INPUT_READ TESTS ====================/

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

// Test: Constructor initializes reading_path correctly
TEST_F(FpgenTest, InputRead_Constructor) {
    std::string test_path = "test_path.txt";
    input_read  reader(test_path);

    // test that the constructor doesn't crash
    EXPECT_NO_THROW(input_read reader2(test_path));
}

// Test: Nonexistent file handling
TEST_F(FpgenTest, InputRead_NonexistentFile) {
    input_read reader(nonexistent_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for nonexistent file";
}

// Test: Empty file handling
TEST_F(FpgenTest, InputRead_EmptyFile) {
    input_read reader(empty_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for empty file";
}

// Test: File with malformed syntax
TEST_F(FpgenTest, InputRead_MalformedSyntax) {
    input_read reader(malformed_input_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for malformed input file";
}

// Test: Edge case values
TEST_F(FpgenTest, InputRead_EdgeCaseValues) {
    std::string   edge_case_file = test_dir + "/edge_case.txt";
    std::ofstream file(edge_case_file);
    file << "source_path = \n" // Empty string
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

    file.close();

    input_read reader(edge_case_file);
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
    std::string   whitespace_file = test_dir + "/whitespace_test.txt";
    std::ofstream file(whitespace_file);
    file << "# This is a comment line\n"
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

    file.close();

    input_read reader(whitespace_file);
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

// Test: Parameter validation
TEST_F(FpgenTest, InputRead_ParameterValidation) {
    // Test logical constraints
    std::string   constraint_file = test_dir + "/constraint_test.txt";
    std::ofstream file(constraint_file);
    file << "source_path = test.h5\n"
         << "shape_path = test.txt\n"
         << "output_path = test.dat\n"
         << "timeslice = 1\n"
         << "manifold = 1\n"
         << "grid_R1 = 0.435\n"
         << "grid_Z1 = -0.239\n"
         << "grid_R2 = 0.435\n"
         << "grid_Z2 = -0.232\n"
         << "nRZ = 10\n"
         << "nPhi = 20\n"
         << "num_threads = 8\n"
         << "max_turns = 10000\n"
         << "h_init = 1e-2\n"
         << "h_min = 1e-6\n"
         << "h_max = 1e-8\n"; // h_max < h_min should be invalid
    file.close();

    input_read reader(constraint_file);
    bool       result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for invalid parameter constraints";
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
maglit     *FpgenTest::tracer = nullptr;
std::string FpgenTest::source_path;

// Test: Constructor and basic initialization
TEST_F(FpgenTest, Footprint_Constructor) {
    // Test floor plate (plate = 0)
    EXPECT_NO_THROW(footprint fp(0, 0.435, -0.239, 0.435, -0.232, 10, 20));

    // Test wall plate (plate = 1)
    EXPECT_NO_THROW(footprint fp(1, 0.435, -0.239, 0.435, -0.232, 15, 30));
}

// Test: Output data structure initialization
TEST_F(FpgenTest, Footprint_OutputDataInitialization) {
    int       nRZ = 10;
    int       nPhi = 20;
    footprint fp(0, 0.5, 0.6, 0.4, 0.5, nRZ, nPhi);

    // Check that outputData is properly sized
    EXPECT_EQ(fp.outputData.size(), nRZ * nPhi);

    // Check that each row has 5 columns
    if (fp.outputData.size() > 0) {
        EXPECT_EQ(fp.outputData[0].size(), 5);
    }

    // Check all rows are properly initialized
    for (const auto &row : fp.outputData) {
        EXPECT_EQ(row.size(), 5);
    }
}
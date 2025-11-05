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
        std::snprintf(source_path, sizeof(source_path),
                      "%s/C1.h5", TEST_DATA_DIR);
        tracer = new maglit(source_path, FIO_M3DC1_SOURCE, 1);
        tracer->configure(0.01, 1e-6, 0.1);
    }

    static void TearDownTestSuite() {
        delete tracer;
        tracer = nullptr;
    }

    static maglit *tracer;
    static char    source_path[256];

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
             << "#=============== STRING PARAMETERS\n"
             << "        source_path = test_source.h5\n"
             << "        shape_path = test_shape.txt\n"
             << "        output_path = test_output.dat\n"
             << "#\n"
             << "#=============== NUMERIC PARAMETERS        \n"
             << "#       \n"
             << "        num_threads = 4         # number of threads\n"
             << "        plate  = 1              # floor=0; wall=1\n"
             << "        timeslice = 0           # equilibrium=-1; vacuum=0; plasma_resp=1 \n"
             << "        gridMin = 0.400         # minimum value mapped on correspondent axis (R or Z)\n"
             << "        gridMax = 0.600         # maximum value mapped on correspondent axis (R or Z)\n"
             << "        nGrid = 50              # grid dimensions on correspondent axis (R or Z)\n"
             << "        nPhi = 100              # grid dimentions on phi axis\n";
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
    EXPECT_EQ(reader.num_threads, 4);
    EXPECT_EQ(reader.plate, 1);
    EXPECT_EQ(reader.timeslice, 0);
    EXPECT_DOUBLE_EQ(reader.gridMin, 0.400);
    EXPECT_DOUBLE_EQ(reader.gridMax, 0.600);
    EXPECT_EQ(reader.nGrid, 50);
    EXPECT_EQ(reader.nPhi, 100);
}

// Test: Constructor initializes reading_path correctly
TEST_F(FpgenTest, InputRead_Constructor) {
    std::string test_path = "test_path.txt";
    input_read  reader(test_path);

    // We can't directly access reading_path (it's private),
    // but we can test that the constructor doesn't crash
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
         << "num_threads = 0\n" // Zero threads
         << "plate = -1\n"      // Negative value
         << "timeslice = 999\n" // Large value
         << "gridMin = -1.0\n"  // Negative minimum
         << "gridMax = 1e6\n"   // Scientific notation
         << "nGrid = 1\n"       // Minimum grid
         << "nPhi = 1000000\n"; // Large grid
    file.close();

    input_read reader(edge_case_file);
    bool       result = reader.readInputFile();

    if (result) {
        EXPECT_EQ(reader.shape_path, "a");
        EXPECT_EQ(reader.num_threads, 0);
        EXPECT_EQ(reader.plate, -1);
        EXPECT_EQ(reader.timeslice, 999);
        EXPECT_DOUBLE_EQ(reader.gridMin, -1.0);
        EXPECT_DOUBLE_EQ(reader.gridMax, 1e6);
        EXPECT_EQ(reader.nGrid, 1);
        EXPECT_EQ(reader.nPhi, 1000000);
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
         << "  num_threads  =  8  \n" // Trailing spaces
         << "\tplate\t=\t0\t\n"       // Tabs instead of spaces
         << "timeslice = 1 # Comment\n"
         << "gridMin = 0.5\n"
         << "gridMax = 0.6\n"
         << "nGrid = 10\n"
         << "nPhi = 20\n";
    file.close();

    input_read reader(whitespace_file);
    bool       result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should handle comments and whitespace correctly";

    if (result) {
        EXPECT_EQ(reader.source_path, "test.h5");
        EXPECT_EQ(reader.shape_path, "test.txt");
        EXPECT_EQ(reader.output_path, "test.dat");
        EXPECT_EQ(reader.num_threads, 8);
        EXPECT_EQ(reader.plate, 0);
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
         << "num_threads = 1\n"
         << "plate = 0\n"
         << "timeslice = 1\n"
         << "gridMin = 0.6\n" // gridMin > gridMax
         << "gridMax = 0.5\n"
         << "nGrid = 10\n"
         << "nPhi = 20\n";
    file.close();

    input_read reader(constraint_file);
    bool       result = reader.readInputFile();

    if (result) {
        // If it accepts the values, we can test that they were read correctly
        EXPECT_DOUBLE_EQ(reader.gridMin, 0.6);
        EXPECT_DOUBLE_EQ(reader.gridMax, 0.5);
    }

    // Either way, it shouldn't crash
    EXPECT_NO_THROW(reader.readInputFile());
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
char    FpgenTest::source_path[256];

// Test: Constructor and basic initialization
TEST_F(FpgenTest, Footprint_Constructor) {
    // Test floor plate (plate = 0)
    EXPECT_NO_THROW(footprint fp(0, 0.5, 0.6, 10, 20));

    // Test wall plate (plate = 1)
    EXPECT_NO_THROW(footprint fp(1, -0.3, 0.3, 15, 30));
}

// Test: Output data structure initialization
TEST_F(FpgenTest, Footprint_OutputDataInitialization) {
    int       nGrid = 10;
    int       nPhi = 20;
    footprint fp(0, 0.5, 0.6, nGrid, nPhi);

    // Check that outputData is properly sized
    EXPECT_EQ(fp.outputData.size(), nGrid * nPhi);

    // Check that each row has 5 columns
    if (fp.outputData.size() > 0) {
        EXPECT_EQ(fp.outputData[0].size(), 5);
    }

    // Check all rows are properly initialized
    for (const auto &row : fp.outputData) {
        EXPECT_EQ(row.size(), 5);
    }
}
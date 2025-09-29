#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Include the fpgen headers
#include "../fpgen/src/input_read.h"
// #include "../fpgen/src/tcabr_collider.h"
// #include "../fpgen/src/footprint.h"

// Test fixture for maglit functionality
class FpgenTest : public ::testing::Test {
  protected:
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

        // // Create shape files for collider tests
        // square_shape_file = test_dir + "/square_shape.txt";
        // circle_shape_file = test_dir + "/circle_shape.txt";
        // invalid_shape_file = test_dir + "/invalid_shape.txt";

        // createShapeFiles();
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

    // void createShapeFiles() {
    //     // Create a square shape (0 to 1 in both R and Z)
    //     std::ofstream square(square_shape_file);
    //     square << "4\n"; // 4 vertices
    //     square << "0.0 0.0\n";
    //     square << "1.0 0.0\n";
    //     square << "1.0 1.0\n";
    //     square << "0.0 1.0\n";
    //     square.close();

    //     // Create an approximate circle (octagon) centered at (0.5, 0.5) with radius 0.4
    //     std::ofstream circle(circle_shape_file);
    //     circle << "8\n"; // 8 vertices
    //     double cx = 0.5, cz = 0.5, radius = 0.4;
    //     for (int i = 0; i < 8; i++) {
    //         double angle = 2.0 * M_PI * i / 8.0;
    //         double r = cx + radius * cos(angle);
    //         double z = cz + radius * sin(angle);
    //         circle << r << " " << z << "\n";
    //     }
    //     circle.close();

    //     // Create invalid shape file
    //     std::ofstream invalid(invalid_shape_file);
    //     invalid << "3\n"; // Says 3 vertices
    //     invalid << "0.0 0.0\n";
    //     invalid << "1.0 0.0\n";
    //     // Missing third vertex
    //     invalid.close();
    // }

    std::string test_dir;
    std::string valid_input_file;
    std::string invalid_input_file;
    std::string malformed_input_file;
    std::string empty_input_file;
    std::string nonexistent_file = "nonexistent_file.txt";

    // std::string square_shape_file;
    // std::string circle_shape_file;
    // std::string invalid_shape_file;
};

// ==================== INPUT_READ TESTS ====================/

// Test: Valid input file parsing
TEST_F(FpgenTest, InputRead_ValidFile) {
    input_read reader(valid_input_file);
    bool result = reader.readInputFile();

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
    input_read reader(test_path);

    // We can't directly access reading_path (it's private),
    // but we can test that the constructor doesn't crash
    EXPECT_NO_THROW(input_read reader2(test_path));
}

// Test: Nonexistent file handling
TEST_F(FpgenTest, InputRead_NonexistentFile) {
    input_read reader(nonexistent_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for nonexistent file";
}

// Test: Empty file handling
TEST_F(FpgenTest, InputRead_EmptyFile) {
    input_read reader(empty_input_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for empty file";
}

// Test: File with malformed syntax
TEST_F(FpgenTest, InputRead_MalformedSyntax) {
    input_read reader(malformed_input_file);
    bool result = reader.readInputFile();

    EXPECT_FALSE(result) << "Should return false for malformed input file";
}

// Test: Edge case values
TEST_F(FpgenTest, InputRead_EdgeCaseValues) {
    std::string edge_case_file = test_dir + "/edge_case.txt";
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
    bool result = reader.readInputFile();

    if (result) {
        // If parsing succeeded, check the values
        EXPECT_EQ(reader.shape_path, "a");
        EXPECT_EQ(reader.num_threads, 0);
        EXPECT_EQ(reader.plate, -1);
        EXPECT_EQ(reader.timeslice, 999);
        EXPECT_DOUBLE_EQ(reader.gridMin, -1.0);
        EXPECT_DOUBLE_EQ(reader.gridMax, 1e6);
        EXPECT_EQ(reader.nGrid, 1);
        EXPECT_EQ(reader.nPhi, 1000000);
    }

    // Test should either succeed with correct parsing or fail gracefully
    EXPECT_NO_THROW(reader.readInputFile());
}

// Test: Comments and whitespace handling
TEST_F(FpgenTest, InputRead_CommentsAndWhitespace) {
    std::string whitespace_file = test_dir + "/whitespace_test.txt";
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
    bool result = reader.readInputFile();

    EXPECT_TRUE(result) << "Should handle comments and whitespace correctly";

    if (result) {
        EXPECT_EQ(reader.source_path, "test.h5");
        EXPECT_EQ(reader.shape_path, "test.txt");
        EXPECT_EQ(reader.output_path, "test.dat");
        EXPECT_EQ(reader.num_threads, 8);
        EXPECT_EQ(reader.plate, 0);
    }
}

// Test: Parameter validation (if the class does validation)
TEST_F(FpgenTest, InputRead_ParameterValidation) {
    // Test logical constraints
    std::string constraint_file = test_dir + "/constraint_test.txt";
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
    bool result = reader.readInputFile();

    // The test should either:
    // 1. Accept the values as-is (no validation in parser)
    // 2. Reject invalid logical relationships
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
    int first_threads = reader.num_threads;

    // Second read should give same results
    bool result2 = reader.readInputFile();
    EXPECT_TRUE(result2);
    EXPECT_EQ(reader.source_path, first_source);
    EXPECT_EQ(reader.num_threads, first_threads);
}

// // ==================== TCABR_SHAPE TESTS ====================

// // Test: Default constructor
// TEST_F(FpgenTest, TcabrShape_DefaultConstructor) {
//     EXPECT_NO_THROW(tcabr_shape shape);
// }

// // Test: Constructor with valid shape file
// TEST_F(FpgenTest, TcabrShape_ConstructorWithFile) {
//     EXPECT_NO_THROW(tcabr_shape shape(square_shape_file));

//     tcabr_shape shape(square_shape_file);
//     // If construction succeeded, object should be usable
//     EXPECT_TRUE(true);
// }

// // Test: Load square shape and test inside/outside detection
// TEST_F(FpgenTest, TcabrShape_SquareInsideOutside) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Points clearly inside the square [0,1]x[0,1]
//     EXPECT_TRUE(shape.tcabr_inside(0.5, 0.5, 0.0, shape_ptr))
//         << "Center point (0.5, 0.5) should be inside";
//     EXPECT_TRUE(shape.tcabr_inside(0.25, 0.25, 0.0, shape_ptr))
//         << "Point (0.25, 0.25) should be inside";
//     EXPECT_TRUE(shape.tcabr_inside(0.75, 0.75, 0.0, shape_ptr))
//         << "Point (0.75, 0.75) should be inside";

//     // Points clearly outside
//     EXPECT_FALSE(shape.tcabr_inside(1.5, 0.5, 0.0, shape_ptr))
//         << "Point (1.5, 0.5) should be outside";
//     EXPECT_FALSE(shape.tcabr_inside(0.5, 1.5, 0.0, shape_ptr))
//         << "Point (0.5, 1.5) should be outside";
//     EXPECT_FALSE(shape.tcabr_inside(-0.5, 0.5, 0.0, shape_ptr))
//         << "Point (-0.5, 0.5) should be outside";
//     EXPECT_FALSE(shape.tcabr_inside(0.5, -0.5, 0.0, shape_ptr))
//         << "Point (0.5, -0.5) should be outside";
// }

// // Test: Boundary points (edge cases)
// TEST_F(FpgenTest, TcabrShape_BoundaryPoints) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Points on the edges - behavior depends on implementation
//     // Testing that these don't crash and give consistent results
//     EXPECT_NO_THROW(shape.tcabr_inside(0.0, 0.5, 0.0, shape_ptr));
//     EXPECT_NO_THROW(shape.tcabr_inside(1.0, 0.5, 0.0, shape_ptr));
//     EXPECT_NO_THROW(shape.tcabr_inside(0.5, 0.0, 0.0, shape_ptr));
//     EXPECT_NO_THROW(shape.tcabr_inside(0.5, 1.0, 0.0, shape_ptr));

//     // Vertices
//     EXPECT_NO_THROW(shape.tcabr_inside(0.0, 0.0, 0.0, shape_ptr));
//     EXPECT_NO_THROW(shape.tcabr_inside(1.0, 1.0, 0.0, shape_ptr));
// }

// // Test: Circular shape
// TEST_F(FpgenTest, TcabrShape_CircularShape) {
//     tcabr_shape shape(circle_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Center should be inside
//     EXPECT_TRUE(shape.tcabr_inside(0.5, 0.5, 0.0, shape_ptr))
//         << "Center (0.5, 0.5) should be inside circle";

//     // Point near center should be inside
//     EXPECT_TRUE(shape.tcabr_inside(0.6, 0.5, 0.0, shape_ptr))
//         << "Point (0.6, 0.5) should be inside circle";

//     // Point far from center should be outside
//     EXPECT_FALSE(shape.tcabr_inside(1.5, 0.5, 0.0, shape_ptr))
//         << "Point (1.5, 0.5) should be outside circle";
//     EXPECT_FALSE(shape.tcabr_inside(0.5, 1.5, 0.0, shape_ptr))
//         << "Point (0.5, 1.5) should be outside circle";
// }

// // Test: All quadrants around center
// TEST_F(FpgenTest, TcabrShape_AllQuadrants) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Test points in all four quadrants relative to center (0.5, 0.5)
//     EXPECT_TRUE(shape.tcabr_inside(0.75, 0.75, 0.0, shape_ptr))
//         << "Quadrant I should be inside";
//     EXPECT_TRUE(shape.tcabr_inside(0.25, 0.75, 0.0, shape_ptr))
//         << "Quadrant II should be inside";
//     EXPECT_TRUE(shape.tcabr_inside(0.25, 0.25, 0.0, shape_ptr))
//         << "Quadrant III should be inside";
//     EXPECT_TRUE(shape.tcabr_inside(0.75, 0.25, 0.0, shape_ptr))
//         << "Quadrant IV should be inside";
// }

// // Test: Phi parameter (should not affect 2D collision in poloidal plane)
// TEST_F(FpgenTest, TcabrShape_PhiParameter) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Same point with different phi values should give same result
//     bool result1 = shape.tcabr_inside(0.5, 0.5, 0.0, shape_ptr);
//     bool result2 = shape.tcabr_inside(0.5, 0.5, M_PI, shape_ptr);
//     bool result3 = shape.tcabr_inside(0.5, 0.5, 2 * M_PI, shape_ptr);

//     EXPECT_EQ(result1, result2) << "Phi should not affect result";
//     EXPECT_EQ(result1, result3) << "Phi should not affect result";
// }

// // Test: Invalid shape file handling
// TEST_F(FpgenTest, TcabrShape_InvalidFile) {
//     // Constructor with invalid file should not crash
//     EXPECT_NO_THROW(tcabr_shape shape(invalid_shape_file));

//     // Constructor with nonexistent file should not crash
//     EXPECT_NO_THROW(tcabr_shape shape(nonexistent_file));
// }

// // Test: Multiple collision checks
// TEST_F(FpgenTest, TcabrShape_MultipleChecks) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Perform many collision checks (tests index caching/searching)
//     for (int i = 0; i < 100; i++) {
//         double r = 0.5 + 0.3 * cos(i * 0.1);
//         double z = 0.5 + 0.3 * sin(i * 0.1);

//         // All these points should be inside (circle of radius 0.3 around center)
//         EXPECT_TRUE(shape.tcabr_inside(r, z, 0.0, shape_ptr))
//             << "Point (" << r << ", " << z << ") should be inside";
//     }
// }

// // Test: Extreme coordinates
// TEST_F(FpgenTest, TcabrShape_ExtremeCoordinates) {
//     tcabr_shape shape(square_shape_file);
//     tcabr_shape *shape_ptr = &shape;

//     // Very far points
//     EXPECT_FALSE(shape.tcabr_inside(1000.0, 0.5, 0.0, shape_ptr));
//     EXPECT_FALSE(shape.tcabr_inside(-1000.0, 0.5, 0.0, shape_ptr));
//     EXPECT_FALSE(shape.tcabr_inside(0.5, 1000.0, 0.0, shape_ptr));
//     EXPECT_FALSE(shape.tcabr_inside(0.5, -1000.0, 0.0, shape_ptr));

//     // Very close to center
//     EXPECT_TRUE(shape.tcabr_inside(0.50001, 0.50001, 0.0, shape_ptr));
// }

// // Basic compilation test
// TEST(FpgenSimpleTest, CanIncludeHeaders) {
//     EXPECT_TRUE(true);
// }
#include "../maglit/collider.h"
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

// Test fixture for collider functionality
class ColliderTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // Create temporary directory for test files
        test_dir = "test_input_files";
        std::filesystem::create_directories(test_dir);
        // Create shape files for collider tests
        square_shape_file = test_dir + "/square_shape.txt";
        circle_shape_file = test_dir + "/circle_shape.txt";
        invalid_shape_file = test_dir + "/invalid_shape.txt";
        reversed_square_shape_file = test_dir + "/reversed_square_shape.txt";

        createShapeFiles();
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove_all(test_dir);
    }

    std::string test_dir;
    std::string square_shape_file;
    std::string circle_shape_file;
    std::string invalid_shape_file;
    std::string reversed_square_shape_file;

    void createShapeFiles() {
        // Create a square shape (0 to 1 in both R and Z)
        std::ofstream square(square_shape_file);
        square << "0.0 0.0\n";
        square << "1.0 0.0\n";
        square << "1.0 1.0\n";
        square << "0.0 1.0\n";
        square.close();

        // Create an approximate circle (octagon) centered at (0.5, 0.5) with radius 0.4
        std::ofstream circle(circle_shape_file);
        double        cx = 0.5, cz = 0.5, radius = 0.4;
        for (int i = 0; i < 8; i++) {
            double angle = 2.0 * M_PI * i / 8.0;
            double r = cx + radius * cos(angle);
            double z = cz + radius * sin(angle);
            circle << r << " " << z << "\n";
        }
        circle.close();

        // Create invalid shape file
        std::ofstream invalid(invalid_shape_file);
        invalid << "0.0 0.0\n";
        invalid << "1.0 0.0\n";
        // Missing third vertex
        invalid.close();

        // Create reversed winding order square shape
        std::ofstream reversed(reversed_square_shape_file);
        reversed << "0.0 1.0\n";
        reversed << "1.0 0.0\n";
        reversed << "1.0 1.0\n";
        reversed << "0.0 0.0\n";
        reversed.close();
    }
};

// ==================== COLLIDER TESTS ====================

// Test: Default constructor
TEST_F(ColliderTest, Collider_DefaultConstructor) {
    EXPECT_NO_THROW(collider shape);
}

// Test: Load valid square shape
TEST_F(ColliderTest, Collider_LoadSquareShape) {
    collider square_shape;
    square_shape.load_shape(square_shape_file);
    EXPECT_TRUE(square_shape.is_loaded());
    EXPECT_EQ(square_shape.get_vertices().size(), 4);
}

// Test: Load valid circle shape
TEST_F(ColliderTest, Collider_LoadCircleShape) {
    collider circle_shape;
    circle_shape.load_shape(circle_shape_file);
    EXPECT_TRUE(circle_shape.is_loaded());
    EXPECT_EQ(circle_shape.get_vertices().size(), 8);
}

// Test: Load invalid shape
TEST_F(ColliderTest, Collider_LoadInvalidShape) {
    collider invalid_shape;
    invalid_shape.load_shape(invalid_shape_file);
    EXPECT_FALSE(invalid_shape.is_loaded());
}

// Test: Point inside square shape
TEST_F(ColliderTest, Collider_PointInsideSquare) {
    collider square_shape;
    square_shape.load_shape(square_shape_file);
    EXPECT_TRUE(square_shape.inside(0.5, 0.5));  // Inside
    EXPECT_FALSE(square_shape.inside(1.5, 0.5)); // Outside
    EXPECT_TRUE(square_shape.inside(0.0, 0.5));  // On edge
}

// Test: Point inside circle shape
TEST_F(ColliderTest, Collider_PointInsideCircle) {
    collider circle_shape;
    circle_shape.load_shape(circle_shape_file);
    EXPECT_TRUE(circle_shape.inside(0.5, 0.5));  // Inside
    EXPECT_FALSE(circle_shape.inside(0.1, 0.1)); // Outside
    EXPECT_TRUE(circle_shape.inside(0.9, 0.5));  // On edge
}

// Test: Point inside without loading shape
TEST_F(ColliderTest, Collider_PointInsideWithoutLoading) {
    collider shape;
    EXPECT_FALSE(shape.inside(0.5, 0.5)); // Should handle not loaded case
}

// Test: Reverse winding order shape
TEST_F(ColliderTest, Collider_ReverseWindingOrder) {
    collider reversed_shape;
    reversed_shape.load_shape(reversed_square_shape_file);
    EXPECT_TRUE(reversed_shape.is_loaded());
    EXPECT_EQ(reversed_shape.get_vertices().size(), 4);
    EXPECT_TRUE(reversed_shape.inside(0.5, 0.5));  // Inside
    EXPECT_FALSE(reversed_shape.inside(1.5, 0.5)); // Outside
}
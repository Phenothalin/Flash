#include <gtest/gtest.h>
#include "element_matrix_builder.hpp"
#include <cmath>

using namespace randflash;

// Test 1: Parse simple formulas
TEST(ElementMatrixBuilderTest, ParseSimpleFormulas) {
    // CH4 (methane)
    auto ch4 = parseFormula("CH4", "CH4");
    EXPECT_EQ(ch4.name, "CH4");
    EXPECT_EQ(ch4.atoms.size(), 2);
    EXPECT_EQ(ch4.atoms["C"], 1);
    EXPECT_EQ(ch4.atoms["H"], 4);

    // C2H6 (ethane)
    auto c2h6 = parseFormula("C2H6", "C2H6");
    EXPECT_EQ(c2h6.atoms["C"], 2);
    EXPECT_EQ(c2h6.atoms["H"], 6);

    // H2O (water)
    auto h2o = parseFormula("H2O", "H2O");
    EXPECT_EQ(h2o.atoms["H"], 2);
    EXPECT_EQ(h2o.atoms["O"], 1);
}

// Test 2: Parse complex formulas
TEST(ElementMatrixBuilderTest, ParseComplexFormulas) {
    // C10H22 (decane)
    auto c10h22 = parseFormula("C10H22", "C10H22");
    EXPECT_EQ(c10h22.atoms["C"], 10);
    EXPECT_EQ(c10h22.atoms["H"], 22);

    // CO2 (carbon dioxide)
    auto co2 = parseFormula("CO2", "CO2");
    EXPECT_EQ(co2.atoms["C"], 1);
    EXPECT_EQ(co2.atoms["O"], 2);
}

// Test 3: Build element matrix for hydrocarbons
TEST(ElementMatrixBuilderTest, BuildHydrocarbonMatrix) {
    std::vector<SpeciesFormula> species = {
        parseFormula("CH4", "CH4"),    // C1
        parseFormula("C2H6", "C2H6"),  // C2
        parseFormula("C3H8", "C3H8")   // C3
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Check dimensions
    EXPECT_EQ(A.size(), 2);      // 2 elements (C, H)
    EXPECT_EQ(A[0].size(), 3);   // 3 species

    // Check element names (sorted alphabetically)
    EXPECT_EQ(elementNames.size(), 2);
    EXPECT_EQ(elementNames[0], "C");
    EXPECT_EQ(elementNames[1], "H");

    // Check matrix values
    //      CH4  C2H6  C3H8
    // C  [  1    2     3  ]
    // H  [  4    6     8  ]
    EXPECT_DOUBLE_EQ(A[0][0], 1.0);  // C in CH4
    EXPECT_DOUBLE_EQ(A[0][1], 2.0);  // C in C2H6
    EXPECT_DOUBLE_EQ(A[0][2], 3.0);  // C in C3H8
    EXPECT_DOUBLE_EQ(A[1][0], 4.0);  // H in CH4
    EXPECT_DOUBLE_EQ(A[1][1], 6.0);  // H in C2H6
    EXPECT_DOUBLE_EQ(A[1][2], 8.0);  // H in C3H8
}

// Test 4: Compute element moles
TEST(ElementMatrixBuilderTest, ComputeElementMoles) {
    std::vector<SpeciesFormula> species = {
        parseFormula("CH4", "CH4"),
        parseFormula("C2H6", "C2H6"),
        parseFormula("C3H8", "C3H8")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Species mole fractions: z = [0.5, 0.3, 0.2]
    std::vector<double> z = {0.5, 0.3, 0.2};

    auto b = computeElementMoles(A, z);

    // Expected element moles:
    // b_C = 0.5*1 + 0.3*2 + 0.2*3 = 0.5 + 0.6 + 0.6 = 1.7
    // b_H = 0.5*4 + 0.3*6 + 0.2*8 = 2.0 + 1.8 + 1.6 = 5.4
    EXPECT_EQ(b.size(), 2);
    EXPECT_NEAR(b[0], 1.7, 1e-12);  // Carbon
    EXPECT_NEAR(b[1], 5.4, 1e-12);  // Hydrogen
}

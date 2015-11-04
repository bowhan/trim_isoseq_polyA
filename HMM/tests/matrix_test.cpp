#include <iostream>
#include <cmath>
#include "matrix.hpp"
#include "gmock/gmock.h"

namespace
{
    class MatrixTest : public ::testing::Test
    {
    protected:
        MatrixTest() : mtrix_(3, 8)
        {
            mtrix_ = 2.3;
        }
    public:
        Matrix<double> mtrix_;
    };
    
    TEST_F(MatrixTest, MoveConstructor)
    {
        Matrix<double> mtrix2 = std::move(mtrix_);
        EXPECT_EQ(mtrix2.col(), 8);
        EXPECT_EQ(mtrix2.row(), 3);
        EXPECT_EQ(mtrix2.size(), 24);
        EXPECT_DOUBLE_EQ(mtrix2(0,0), 2.3);
    }
    
    TEST_F(MatrixTest, MoveAssignmentOperator)
    {
        Matrix<double> mtrix2;
        mtrix2 = std::move(mtrix_);
        EXPECT_EQ(mtrix2.col(), 8);
        EXPECT_EQ(mtrix2.row(), 3);
        EXPECT_EQ(mtrix2.size(), 24);
        EXPECT_DOUBLE_EQ(mtrix2(0,0), 2.3);
    }
    
    TEST_F(MatrixTest, Column)
    {
        EXPECT_EQ(mtrix_.col(), 8);
    }
    
    TEST_F(MatrixTest, Row)
    {
        EXPECT_EQ(mtrix_.row(), 3);
    }
    
    TEST_F(MatrixTest, Size)
    {
        EXPECT_EQ(mtrix_.size(), 24);
    }
    
    TEST_F(MatrixTest, ValueRetrieve)
    {
        EXPECT_DOUBLE_EQ(mtrix_(0,0), 2.3);
        EXPECT_DOUBLE_EQ(mtrix_(1,4), 2.3);
        EXPECT_DOUBLE_EQ(mtrix_(2,7), 2.3);
    }
    
    TEST_F(MatrixTest, ValueAssignment)
    {
        mtrix_(0,2) = 1.24;
        mtrix_(2,1) = 33.9;
        mtrix_(2,7) = .321;
        EXPECT_DOUBLE_EQ(mtrix_(0,2), 1.24);
        EXPECT_DOUBLE_EQ(mtrix_(2,1), 33.9);
        EXPECT_DOUBLE_EQ(mtrix_(2,7), .321);
    }
    
    TEST_F(MatrixTest, ValueAssignment2)
    {
        mtrix_ = {
            1, 2, 3, 4, 5, 6, 7, 8,
            9,10,11,12,13,14,15,16,
            17,18,19,20,21,22,23,24
        };
        EXPECT_DOUBLE_EQ(mtrix_(0,2), 3);
        EXPECT_DOUBLE_EQ(mtrix_(2,1), 18);
        EXPECT_DOUBLE_EQ(mtrix_(2,7), 24);
    }
    
    TEST_F(MatrixTest, Add)
    {
        mtrix_(0,3) = 1.2;
        mtrix_(2,5) = 2.3;
        mtrix_(1,7) = 23.4;
        mtrix_ += 3.0;
        EXPECT_DOUBLE_EQ(mtrix_(0,3), 1.2 + 3.);
        EXPECT_DOUBLE_EQ(mtrix_(2,5), 2.3 + 3.);
        EXPECT_DOUBLE_EQ(mtrix_(1,7), 23.4+ 3.);
    }
    
    TEST_F(MatrixTest, Substract)
    {
        mtrix_(0,3) = 1.2;
        mtrix_(2,5) = 2.3;
        mtrix_(1,7) = 23.4;
        mtrix_ -= 3.0;
        EXPECT_DOUBLE_EQ(mtrix_(0,3), 1.2 - 3.);
        EXPECT_DOUBLE_EQ(mtrix_(2,5), 2.3 - 3.);
        EXPECT_DOUBLE_EQ(mtrix_(1,7), 23.4- 3.);
    }
    
    TEST_F(MatrixTest, Multiplication)
    {
        const double m = 1.2112;
        mtrix_(0,3) = 1.2;
        mtrix_(2,5) = 2.3;
        mtrix_(1,7) = 23.4;
        mtrix_ *= m;
        EXPECT_DOUBLE_EQ(mtrix_(0,3), 1.2 * m);
        EXPECT_DOUBLE_EQ(mtrix_(2,5), 2.3 * m);
        EXPECT_DOUBLE_EQ(mtrix_(1,7), 23.4* m);
    }
    
    TEST_F(MatrixTest, Division)
    {
        const double m = 3.2222;
        mtrix_(0,3) = 1.2;
        mtrix_(2,5) = 2.3;
        mtrix_(1,7) = 23.4;
        mtrix_ /= m;
        EXPECT_DOUBLE_EQ(mtrix_(0,3), 1.2 / m);
        EXPECT_DOUBLE_EQ(mtrix_(2,5), 2.3 / m);
        EXPECT_DOUBLE_EQ(mtrix_(1,7), 23.4/ m);
    }
    
    TEST_F(MatrixTest, Log2)
    {
        mtrix_(0,3) = 1.2;
        mtrix_(2,5) = 2.3;
        mtrix_(1,7) = 23.4;
        mtrix_.log();
        EXPECT_DOUBLE_EQ(mtrix_(0,3), std::log2(1.2));
        EXPECT_DOUBLE_EQ(mtrix_(2,5), std::log2(2.3));
        EXPECT_DOUBLE_EQ(mtrix_(1,7), std::log2(23.4));
    }
    
    TEST_F(MatrixTest, ColumnSum)
    {
        mtrix_ = {
            1, 2, 3, 4, 5, 6, 7, 8,
            9,10,11,12,13,14,15,16,
            17,18,19,20,21,22,23,24
        };
        EXPECT_DOUBLE_EQ(mtrix_.colSum(0), 27);
        EXPECT_DOUBLE_EQ(mtrix_.colSum(1), 30);
        EXPECT_DOUBLE_EQ(mtrix_.colSum(2), 33);
        EXPECT_DOUBLE_EQ(mtrix_.colSum(7), 48);
    }
    
    TEST_F(MatrixTest, RowSum)
    {
        mtrix_ = {
            1, 2, 3, 4, 5, 6, 7, 8,
            9,10,11,12,13,14,15,16,
            17,18,19,20,21,22,23,24
        };
        EXPECT_DOUBLE_EQ(mtrix_.rowSum(0), 36);
        EXPECT_DOUBLE_EQ(mtrix_.rowSum(1), 100);
        EXPECT_DOUBLE_EQ(mtrix_.rowSum(2), 164);
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}


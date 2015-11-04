
#include <string>
#include "polyA_hmm_model.hpp"
#include "hmm_utilities.h"
#include "gmock/gmock.h"

using namespace std;
namespace
{
    class PolyAHmmModeTest : public ::testing::Test
    {
    protected:
        PolyAHmmModeTest() : hmm()
        {
            hmm.initialProb(0) = 0.5;
            hmm.initialProb(1) = 0.5;

            hmm.transProb(PolyAHmmMode<string>::States::POLYA, PolyAHmmMode<string>::States::POLYA) = 0.7;
            hmm.transProb(PolyAHmmMode<string>::States::POLYA, PolyAHmmMode<string>::States::NONPOLYA) = 0.3;
            hmm.transProb(PolyAHmmMode<string>::States::NONPOLYA, PolyAHmmMode<string>::States::POLYA) = 0.0;
            hmm.transProb(PolyAHmmMode<string>::States::NONPOLYA, PolyAHmmMode<string>::States::NONPOLYA) = 1.0;
            
            hmm.emitProb(PolyAHmmMode<string>::States::POLYA, to_idx['A']) = 0.96;
            hmm.emitProb(PolyAHmmMode<string>::States::POLYA, to_idx['C']) = 0.01;
            hmm.emitProb(PolyAHmmMode<string>::States::POLYA, to_idx['G']) = 0.01;
            hmm.emitProb(PolyAHmmMode<string>::States::POLYA, to_idx['T']) = 0.01;
            hmm.emitProb(PolyAHmmMode<string>::States::NONPOLYA, to_idx['A']) = 0.3;
            hmm.emitProb(PolyAHmmMode<string>::States::NONPOLYA, to_idx['C']) = 0.2;
            hmm.emitProb(PolyAHmmMode<string>::States::NONPOLYA, to_idx['G']) = 0.2;
            hmm.emitProb(PolyAHmmMode<string>::States::NONPOLYA, to_idx['T']) = 0.3;
        }
    public:
        PolyAHmmMode<string> hmm;
    };
    
    TEST_F(PolyAHmmModeTest, Constructor)
    {
        EXPECT_DOUBLE_EQ(hmm.states(),  2);
        EXPECT_DOUBLE_EQ(hmm.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm.initialProb(1), 0.5);
        EXPECT_DOUBLE_EQ(hmm.transProb(0,0), 0.7);
        EXPECT_DOUBLE_EQ(hmm.transProb(0,1), 0.3);
        EXPECT_DOUBLE_EQ(hmm.transProb(1,0), 0);
        EXPECT_DOUBLE_EQ(hmm.transProb(1,1), 1.0);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0,0), 0.96);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0,1), 0.01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0,2), 0.01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0,3), 0.01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1,0), 0.3);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1,1), 0.2);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1,2), 0.2);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1,3), 0.3);
    }
    
    TEST_F(PolyAHmmModeTest, MoveConstructor)
    {
        PolyAHmmMode<string> hmm2 = std::move(hmm);
        EXPECT_DOUBLE_EQ(hmm2.states(),  2);
        EXPECT_DOUBLE_EQ(hmm2.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(1), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0,0), 0.7);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0,1), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1,0), 0);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1,1), 1.0);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,0), 0.96);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,1), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,2), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,3), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,0), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,1), 0.2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,2), 0.2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,3), 0.3);
    }
    
    TEST_F(PolyAHmmModeTest, MoveAssignmentOperator)
    {
        PolyAHmmMode<string> hmm2; hmm2 = std::move(hmm);
        EXPECT_DOUBLE_EQ(hmm2.states(),  2);
        EXPECT_DOUBLE_EQ(hmm2.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(1), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0,0), 0.7);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0,1), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1,0), 0);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1,1), 1.0);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,0), 0.96);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,1), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,2), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0,3), 0.01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,0), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,1), 0.2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,2), 0.2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1,3), 0.3);
    }
    
    TEST_F(PolyAHmmModeTest, ForwardAlgorithm)
    {
        const Matrix<double>& fwd = hmm.calculateForward("GAAC");
        Matrix<double> answer (2,4);
        answer = {-7.643856, -8.217323, -8.790790, -15.94922,
                -3.321928, -5.037414, -6.727395, -8.94932};
        for(size_t i = 0; i < 2; ++i)
        {
            for(size_t j = 0; j < 4; ++j)
            {
                EXPECT_NEAR(answer(i,j), fwd(i,j), 1e-5);
            }
        }
    }
    
    TEST_F(PolyAHmmModeTest, BackwardAlgorithm)
    {
        const Matrix<double>& bak = hmm.calculateBackward("GAAC");
        Matrix<double> answer (2,4);
        answer = {-4.388291, -3.987955, -3.899695, 0,
                -5.795859, -4.058894, -2.321928, 0};
        for(size_t i = 0; i < 2; ++i)
        {
            for(size_t j = 0; j < 4; ++j)
            {
                EXPECT_NEAR(answer(i,j), bak(i,j), 1e-5);
            }
        }
    }
    
    TEST_F(PolyAHmmModeTest, VirtabiAlgorithm)
    {
        hmm.calculateVirtabi("AAAC");
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

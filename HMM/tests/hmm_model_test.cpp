#include <string>
#include "hmm_model.hpp"
#include "gmock/gmock.h"

using namespace std;
namespace
{
    class HmmModeBaseTest : public ::testing::Test
    {
    protected:
        HmmModeBaseTest(): hmm(2,4)
        {
            hmm.initialProb(0) = .5;
            hmm.initialProb(1) = .5;
            
            hmm.transProb(0, 0) = .7;
            hmm.transProb(0, 1, .3);
            hmm.transProb(1, 0, .3);
            hmm.transProb(1, 1, .7);
            
            hmm.emitProb(0, 0, .96);
            hmm.emitProb(0, 1, .01);
            hmm.emitProb(0, 2, .01);
            hmm.emitProb(0, 3) = .01;
            
            hmm.emitProb(1, 0, .3);
            hmm.emitProb(1, 1, .2);
            hmm.emitProb(1, 2, .3);
            hmm.emitProb(1, 3) = .2;
        }
        
    public:
        HmmModeBase hmm;
    };
    
    TEST_F(HmmModeBaseTest, Constructor)
    {
        EXPECT_EQ(hmm.states(), 2);
        EXPECT_EQ(hmm.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm.initialProb(1), 0.5);
        
        EXPECT_DOUBLE_EQ(hmm.transProb(0, 0), 0.7);
        EXPECT_DOUBLE_EQ(hmm.transProb(0, 1), 0.3);
        EXPECT_DOUBLE_EQ(hmm.transProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm.transProb(1, 1), .7);
        
        EXPECT_DOUBLE_EQ(hmm.emitProb(0, 0), .96);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0, 1), .01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0, 2), .01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(0, 3), .01);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1, 1), .2);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1, 2), .3);
        EXPECT_DOUBLE_EQ(hmm.emitProb(1, 3), .2);
    }
    
    TEST_F(HmmModeBaseTest, MoveConstructor)
    {
        HmmModeBase hmm2 = std::move(hmm);
        EXPECT_EQ(hmm2.states(), 2);
        EXPECT_EQ(hmm2.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(1), 0.5);
        
        EXPECT_DOUBLE_EQ(hmm2.transProb(0, 0), 0.7);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0, 1), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1, 1), .7);
        
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 0), .96);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 1), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 2), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 3), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 1), .2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 2), .3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 3), .2);
    }
    
    TEST_F(HmmModeBaseTest, MoveAssignmentOperator)
    {
        HmmModeBase hmm2; hmm2 = std::move(hmm);
        EXPECT_EQ(hmm2.states(), 2);
        EXPECT_EQ(hmm2.symbols(), 4);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(0), 0.5);
        EXPECT_DOUBLE_EQ(hmm2.initialProb(1), 0.5);
        
        EXPECT_DOUBLE_EQ(hmm2.transProb(0, 0), 0.7);
        EXPECT_DOUBLE_EQ(hmm2.transProb(0, 1), 0.3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm2.transProb(1, 1), .7);
        
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 0), .96);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 1), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 2), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(0, 3), .01);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 0), .3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 1), .2);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 2), .3);
        EXPECT_DOUBLE_EQ(hmm2.emitProb(1, 3), .2);
    }
    
}

int main(int argc, char** argv)
{
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

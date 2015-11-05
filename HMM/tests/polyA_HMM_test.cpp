
#include <vector>
#include <list>
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
            hmm.initialProb(PolyAHmmMode::States::POLYA)    = 0.5;
            hmm.initialProb(PolyAHmmMode::States::NONPOLYA) = 0.5;

            hmm.transProb(PolyAHmmMode::States::POLYA, PolyAHmmMode::States::POLYA)       = 0.7;
            hmm.transProb(PolyAHmmMode::States::POLYA, PolyAHmmMode::States::NONPOLYA)    = 0.3;
            hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::POLYA)    = 0.0;
            hmm.transProb(PolyAHmmMode::States::NONPOLYA, PolyAHmmMode::States::NONPOLYA) = 1.0;
            
            hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['A'])    = 0.96;
            hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['C'])    = 0.01;
            hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['G'])    = 0.01;
            hmm.emitProb(PolyAHmmMode::States::POLYA, to_idx['T'])    = 0.01;
            hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['A']) = 0.3;
            hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['C']) = 0.2;
            hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['G']) = 0.2;
            hmm.emitProb(PolyAHmmMode::States::NONPOLYA, to_idx['T']) = 0.3;
        }
    public:
        PolyAHmmMode hmm;
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
        PolyAHmmMode hmm2 = std::move(hmm);
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
        PolyAHmmMode hmm2; hmm2 = std::move(hmm);
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
    
    TEST_F(PolyAHmmModeTest, ForwardAlgorithmString)
    {
		const Matrix<double>& fwd = hmm.calculateForward(string{"GAAC"});
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
	
	TEST_F(PolyAHmmModeTest, ForwardAlgorithmCString)
	{ // string literal
		const Matrix<double>& fwd = hmm.calculateForward("AAAAAACAGTCGACGAAAAA");
		Matrix<double> answer (2,20);
		answer = {-1.058894,-1.632361,-2.205827,-2.779294,-3.352761,-3.926228,-11.084657,-11.658124,-18.816554,-25.974983,-33.133412,-40.291842,-40.865309,-48.023738,-55.182167,-55.755634,-56.329101,-56.902568,-57.476035,-58.049502,
			-2.736966,-3.503078,-4.171618,-4.789639,-5.383449,-5.966092,-7.128636,-8.837980,-11.099890,-12.834799,-15.156679,-17.478606,-19.215571,-21.537499,-23.859427,-25.596393,-27.333359,-29.070324,-30.807290,-32.544255
		};
		for(size_t i = 0; i < 2; ++i)
		{
			for(size_t j = 0; j < 20; ++j)
			{
				EXPECT_NEAR(answer(i,j), fwd(i,j), 1e-3);
			}
		}
	}
	
	TEST_F(PolyAHmmModeTest, ForwardAlgorithmCaseInsensistiveString)
	{ // caseInsensitiveString
		const Matrix<double>& fwd = hmm.calculateForward(caseInsensitiveString{"GAAC"});
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

	TEST_F(PolyAHmmModeTest, ForwardAlgorithmIterator)
	{ // iterator
		vector<char> str = {'G', 'A', 'A', 'C'};
		const Matrix<double>& fwd = hmm.calculateForward(str.begin(), str.size());
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
	
    TEST_F(PolyAHmmModeTest, BackwardAlgorithmCString)
	{ // string literal
		const Matrix<double>& bak = hmm.calculateBackward("AAAAAACAGTCGACGAAAAA");
		Matrix<double> answer (2,20);
		answer = {-31.4992583,-30.9392593,-30.3964246,-29.8940298,-29.4937411,-29.4056497,-25.5110163,-25.4545512,-23.1488508,-21.3901001,-18.9374239,-14.8182778,-14.6651073,-9.5510579,-2.5597362,-1.9924837,-1.4330348,-0.8914715,-0.3921371,0.0000000,-36.512121,-34.775156,-33.038190,-31.301225,-29.564259,-27.827293,-25.505365,-23.768400,-21.446472,-19.709506,-17.387578,-15.065650,
			-13.328684,-11.006756,-8.684828,-6.947862,-5.210897,-3.473931,-1.736966,0.000000};
		for(size_t i = 0; i < 2; ++i)
		{
			for(size_t j = 0; j < 20; ++j)
			{
				EXPECT_NEAR(answer(i,j), bak(i,j), 1e-5);
			}
		}
	}

	TEST_F(PolyAHmmModeTest, BackwardAlgorithmString)
	{ // std::string
		const Matrix<double>& bak = hmm.calculateBackward(string{"GAAC"});
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
	
	TEST_F(PolyAHmmModeTest, BackwardAlgorithmCaseInsensitiveString)
	{ // caseInsensitiveString
		const Matrix<double>& bak = hmm.calculateBackward(caseInsensitiveString{"GAAC"});
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
	
	TEST_F(PolyAHmmModeTest, BackwardAlgorithmIterator)
	{
		list<char> str {'G', 'A', 'A', 'C'};
		/*
		 */
		auto iter = str.begin();
		std::advance(iter, str.size() - 1);
		
		const Matrix<double>& bak = hmm.calculateBackward(iter, str.size());
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

    TEST_F(PolyAHmmModeTest, VirtabiAlgorithm1)
    {
		const Matrix<int>& path = hmm.calculateVirtabi("AAAAAAC");
		Matrix<int> answer(1,7); answer = {0,0,0,0,0,0,1};
		for (int i=0; i < 7; ++i)
		{
			EXPECT_DOUBLE_EQ(answer[i], path[i]);
		}
    }
	
	TEST_F(PolyAHmmModeTest, VirtabiAlgorithm2)
	{
		const Matrix<int>& path = hmm.calculateVirtabi("AAAAAACAGTCGACGAAAAA");
		Matrix<int> answer(1,path.size()); answer = {0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
		for (int i=0; i < 20; ++i)
		{
			EXPECT_DOUBLE_EQ(answer[i], path[i]);
		}
	}
	
}

int main(int argc, char** argv)
{
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

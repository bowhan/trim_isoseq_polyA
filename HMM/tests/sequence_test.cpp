#include <string>
#include "gmock/gmock.h"
#include "sequence.hpp"

using namespace std;
namespace
{
    class SequenceTest : public ::testing::Test
    {
    protected:
        SequenceTest() : seq("GGATCGATCcatcga")
        { }
    public:
        Sequence<> seq;
    };
    
    TEST_F(SequenceTest, Constructor1)
    {
        Sequence<> seq2 {"GGATCGATCCATCGA"};
        fprintf(stderr, "%s\n", seq.seq_.c_str());
        fprintf(stderr, "%s\n", seq2.seq_.c_str());
        EXPECT_TRUE(seq == seq2);
    }
    
    TEST_F(SequenceTest, Reverse)
    {
        Sequence<> seq2 {"agctacCTAGCTAGG"};
        seq.reverse();
        EXPECT_TRUE(seq == seq2);
    }
    
    TEST_F(SequenceTest, ReverseCopy)
    {
        Sequence<> seq2 {"agctacCTAGCTAGG"};
        auto seq3 = seq.reverse_copy();
        EXPECT_TRUE(seq3 == seq2);
    }
    
    TEST_F(SequenceTest, Complement1)
    {
        Sequence<> seq2 {"CCTAGCTAGGTAGCT"};
        seq.complement();
        EXPECT_TRUE(seq == seq2);
    }

    TEST_F(SequenceTest, Complement2)
    {
        Sequence<> seq2 {"CCTAGCTAGGTAGC"};
        Sequence<> seq3 {"GGATCGATCcatcg"};
        seq2.complement();
        EXPECT_TRUE(seq3 == seq2);
    }
    
    TEST_F(SequenceTest, ComplementCopy)
    {
        auto seq2 = seq.complement_copy();
        Sequence<> seq3 {"CCTAGCTAGGTAGCT"};
        EXPECT_TRUE(seq3 == seq2);
    }
    
    TEST_F(SequenceTest, ReverseComplement)
    {
        Sequence<> seq2 {"TCGATGGATCGATCC"};
        seq.reverse_complement();
        EXPECT_TRUE(seq == seq2);
    }
    
    TEST_F(SequenceTest, ReverseComplementCopy)
    {
        Sequence<> seq2 {"TCGATGGATCGATCC"};
        auto seq3 = seq.reverse_complement_copy();
        EXPECT_TRUE(seq3 == seq2);
    }
    
    TEST_F(SequenceTest, Size)
    {
        EXPECT_TRUE(seq.size() == seq.seq_.size());
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
}

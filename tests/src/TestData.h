#ifndef TESTDATA_H
#define TESTDATA_H

#include <string>

namespace tests {

    const std::string Source_Dir   = std::string("/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/ext/pi/trim_isoseq_polyA/tests/src/");
    const std::string Bin_Dir      = std::string("/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/ext/pi/trim_isoseq_polyA/tests/bin/");
    const std::string Data_Dir     = std::string("/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/ext/pi/trim_isoseq_polyA/tests/data/");
    const std::string Out_Dir      = std::string("/home/UNIXHOME/yli/git/depot/software/smrtanalysis/bioinformatics/ext/pi/trim_isoseq_polyA/tests/out/");


    const std::string testReader_Fasta      = Data_Dir + "testReader.fa";
    const std::string polyA_Fasta           = Data_Dir + "polyA.fa";
    const std::string polyA_Fastq           = Data_Dir + "polyA.fq";
    const std::string polyA_train_Fasta     = Data_Dir + "polyA_train.fa";
    const std::string non_polyA_train_Fasta = Data_Dir + "non_polyA_train.fa";
    const std::string hmm_Default_Out       = Out_Dir + "HMM_default.txt";
} // namespace tests

#endif // TESTDATA_H

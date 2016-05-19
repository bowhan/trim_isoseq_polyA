# Copyright (c) 2015, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.

# Author: Bo Han

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# this is an R implementation of polyA trimmer using its HMM package            #
# it is EXTREMELY slow; but the purpose is to validate our C++ implementation   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 
library(HMM)
library(seqinr)

# arg parse
argv = commandArgs (TRUE)
input_fa = argv[1]

# build HMM model
States = c("polyA", "non-polyA")
Symbols = c("A", "C", "G", "T")
# init prob
polyA_start_pro = 0.5
startProbs= c(polyA_start_pro, 1 - polyA_start_pro)
# trans prob
polyA_to_non_polyA_rate = 0.7
transProbs = matrix(data=c(polyA_to_non_polyA_rate, 1.0 - polyA_to_non_polyA_rate, 0.0, 1.0), nrow=length(States), ncol=length(States), byrow=T)
# emit prob
polyA_emission = c(0.96, 0.01, 0.01, 0.01)
non_polyA_emission = c(0.3, 0.2, 0.2, 0.3)
emissionProbs = matrix(data=c(polyA_emission, non_polyA_emission), nrow=length(States), ncol=length(Symbols), byrow=T)
# model
HMM = initHMM(States, Symbols, startProbs=startProbs, transProbs=transProbs, emissionProbs=emissionProbs)

# aux functions
reverseStr = function(a) 
{ 
    paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="") 
}

trimPloyA = function(HMM, seq)
{
    states = viterbi(HMM,unlist(strsplit(reverseStr(seq), split="")))
    return (states)
}

findFirstNonPolyA = function(s)
{
    for(i in 1:length(s)) 
    { 
        if(s[i]=="non-polyA") return (i); 
    }
    return -1;
}

# process fasta
fh = read.fasta(file=input_fa, as.string=T, forceDNAtolower=F)
output_fa_file = paste(input_fa, "trimmed_by_R", sep=".")
trimedfa = lapply(fh, 
    function(fa) 
    {
        states = trimPloyA(HMM, toupper(fa))
        polyA_len = findFirstNonPolyA(states) - 1
        return (substr(fa, 1, nchar(fa) - polyA_len)) # TODO: need to avoid the case where all sequences are polyA
    })
write.fasta(trimedfa, names(fh), output_fa_file, open="w")

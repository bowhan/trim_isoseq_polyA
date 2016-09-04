#include <assert.h>
#include "flnc.hpp"

namespace {
/* distance in the header section of flnc file from "C" in "_CCS" to the first digit after "fiveend=" */
const size_t k_header_distance1 = 57;
}

void adjustHeader(std::string& s, size_t polyalen) {
    // fl.trimmed.fasta file only:
    // format: <movie_name>/<ZMW name>/<start>_<end>_<CCS> strand=[+|-];fiveseen=1;polyAseen=1;threeseen=1;fiveend=[\d+];polyAend=[\d+];threeend=[\d+];primer=1;chimra=NA

    // example: (- strand), length of CCS: 1564
    //          >m160403_065056_42175_c100993391270000001823222007191686_s1_p0/13/1533_53_CCS strand=-;fiveseen=1;polyAseen=1;threeseen=1;fiveend=31;polyAend=1511;threeend=1535;primer=1;chimera=NA
    //                                                                            | start; 1533 == 1564 - fiveseen(31)
    //                                                                                | end <- update this one; 53 = 1564 - polyAend (1511)

    // example: (+ strand), length of CCS: 1534 (not used)
    //          >m160403_065056_42175_c100993391270000001823222007191686_s1_p0/9/30_1487_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=1487;threeend=1514;primer=1;chimera=NA
    //                                                                           | start; fiveend (30)
    //                                                                              | end <- update this one; 1487 = polyAend
    auto l = s.cbegin();
    while (*l++ != '/'); // l is now point to the char pass the 1st '\'       l
    while (*l++ != '/'); // l is now point to the char pass the 2nd '\'         l
    auto r = l; //                                                              r
    std::string newstr{s.cbegin(), l}; // >m16................................../
    newstr.reserve(s.size() + 3);
    while (*r++ != '_'); //                                                         r
    // this will likely throw exception if the format is not Iso-Seq specific, trailing _ is find for stoi
    int start = std::stoi(std::string{l, r});
    l = r; //                                                                       l
    while (*r++ != '_'); //                                                              r
    assert(*r == 'C' && "invalid flnc file!");
    int end = std::stoi(std::string{l, r}); // this will likely throw exception if the format is not Iso-Seq specific
    if (start < end) { // + strand, update end
        end -= polyalen;
    } else { // - strand, update start
        end += polyalen;
    }
    newstr += std::to_string(start) + '_' + std::to_string(end) + '_';
    newstr.append(s.c_str() + (r - s.cbegin()), k_header_distance1
        - 1); // append "CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend="
    r += k_header_distance1; // jump from "C" in "_CCS" to the first digit after in "fiveend="
    l = r - 1; // l at "=" of "fiveend=""
    while (*r++ != '='); // r is the first digit of polyAseen=
    newstr.append(l, r); // append the value of fiveend
    l = r; // l is not at the first digit of the original polyAseen
    while (*++r != ';'); // r is at ';'
    int polyAend = std::stoi(std::string{l, r});
    polyAend -= polyalen; // polyAend further away from threeend
    newstr.append(std::to_string(polyAend));
    newstr.append(r, s.cend());
    s.swap(newstr);
}


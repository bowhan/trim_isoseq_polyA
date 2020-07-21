# trim_isoseq_polyA
This is a program to trim the polyA tails of DNA sequences in a FastQ format, using HMM.

## Prerequisite
- C++11
- Boost
- pthread
- ~~Zlib~~ 
- ~~BZib2~~ <br>
Zlib and BZlib2 are needed if you need to support gzip/bzip2 compressed input files.


## Install

```bash
mkdir build && \
cd build && \
cmake ../ -DBOOST_ROOT=PATH_TO_YOUR_BOOST_ROOT_DIR -DCMAKE_BUILD_TYPE=Release && \
make
```
    Executable `trim_isoseq_polyA` will be saved in directory trim_isoseq_polyA/bin.
    Please replace PATH_TO_YOUR_BOOST_ROOT_DIR with your own boost root directory.
    e.g., -DBOOST_ROOT=~/mylib/boost/boost_1_60_0/

#### With native support for gzip/bzip2 compressed input
```bash
mkdir build && \
cd build && \
cmake ../ -DBOOST_ROOT=PATH_TO_YOUR_BOOST_ROOT_DIR -DCMAKE_BUILD_TYPE=Release \
    -DSUPPORT_COMPRESSED_INPUT=ON && \
make
```

#### Build with unit tests
```bash
mkdir buildwtest && \
cd buildwtest && \
cmake ../ -DBOOST_ROOT=PATH_TO_YOUR_BOOST_ROOT_DIR -DTrimIsoseqPolyA_build_tests=ON && \
make && make test
```

## Usage
To process Iso-Seq `classfy` output
```bash
trim_isoseq_polyA -i isoseq.flnc.fq -t 8 > isoseq.flnc.atrim.fq 2> isoseq.flnc.atrim.log
```

To process generic fasta files
```bash
trim_isoseq_polyA -i input.fq -t 8 -G > input.atrim.fq 2> input.atrim.log
```
`input.atrim.fq` file contain the fasta entries with polyA trimmed, based on a default HMM model trained with PacBio data.

`input.atrim.log` is a tab file with length of polyA been trimmed.

To visualize polyA (colored red when visualized by `cat`)
```bash
trim_isoseq_polyA -i isoseq.flnc.fq -t 8 -c 2>/dev/null
```

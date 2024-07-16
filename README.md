# isONclust3
A rust implementation of a novel de novo clustering algorithm.
isONclust3 is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene family. Output is a tsv file with each read assigned to a cluster-ID and a folder 'fastq' containing one fastq file per cluster generated. Detailed information is available in [paper](https://link.springer.com/chapter/10.1007/978-3-030-17083-7_14).

# Table of contents
1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Output](#output)
4. [Running isONform](#Running)
    1. [Running a test](#runtest)
5. [Contact](#contact)
6. [Credits](#credits)
## Installation <a name="installation"></a>
At the moment building from source is the only option to install the tool. For this please install rust via
`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` (for macOS and Linux or other Unix-based OS). For windows please follow the instructions on the following site 
`cd isONclust3` 
`cargo build --release` ( Compile the current package, the executable is located in target/release)

## Introduction <a name="introduction"></a>




## Running isONform <a name="Running"></a>

To run the algorithm:<br />

## Contact <a name="contact"></a>
If you encounter any problems, please raise an issue on the issues page, you can also contact the developer of this repository via:
alexander.petri[at]math.su.se


## Credits <a name="credits"></a>

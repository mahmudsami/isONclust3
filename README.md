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
## Installation Guide <a name="installationguide"></a>
At the moment building from source is the only option to install the tool. This requires users to install the Rust programming language onto their system.
## Installing Rust <a name="installingrust"></a>
You can install rust via<br />

`curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` (for macOS and Linux or other Unix-based OS). For Windows please follow the instructions on the following site: https://forge.rust-lang.org/infra/other-installation-methods.html .<br />
## Installation <a name="installation"></a>
After cloning the repository via `git clone https://github.com/aljpetri/isONclust3.git` use the following two commands to compile the code: <br />
`cd isONclust3` <br />
`cargo build --release` ( Compile the current package, the executable is then located in target/release) <br />

## Introduction <a name="introduction"></a>




## Running isONclust3 <a name="Running"></a>

To run the algorithm:<br />

## Contact <a name="contact"></a>
If you encounter any problems, please raise an issue on the issues page, you can also contact the developer of this repository via:
alexander.petri[at]math.su.se


## Credits <a name="credits"></a>

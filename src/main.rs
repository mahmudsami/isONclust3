//#![allow(warnings)]

use std::fs::File;
use std::collections::{HashMap, HashSet, VecDeque};
use rayon::prelude::*;
extern crate rayon;
use std::time::Instant;
use std::time::Duration;
//use crate::generate_sorted_fastq_new_version::{filter_minimizers_by_quality, Minimizer,get_kmer_minimizers};
//use clap::{arg, command, Command};
use clap::Parser;

pub mod file_actions;
mod clustering;
//mod generate_sorted_fastq_for_cluster;
mod generate_sorted_fastq_new_version;
mod generate_sorted_fastq_for_cluster;
mod gff_handling;
use std::path::{PathBuf, Path};
//mod isONclust;
mod structs;
use crate::structs::{FastaRecord, FastqRecord_isoncl_init};
use std::thread;

extern crate clap;
mod write_output;
use memory_stats::memory_stats;
use rustc_hash::{FxHashMap,FxHashSet};

use std::io::Read;

use bio::io::fasta;
use bio::io::fasta::FastaRead;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use std::convert::TryFrom;
use crate::clustering::post_clustering_new;


fn compute_d() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0_f64).powf(exponent).min(0.79433);
    }
    d
}

fn expected_number_errornous_kmers(quality_string: &str, k: usize, d:&[f64;128]) -> f64 {
    //computes the expeced number of errornous kmers for a read by analysing the quality entry
    let prob_error: Vec<f64> = quality_string.chars().map(|char_| d[char_ as u8 as usize]).collect();
    let mut sum_of_expectations = 0.0;
    let mut qurrent_prob_no_error = 1.0;
    let mut window: VecDeque<f64> = VecDeque::with_capacity(k);

    for (i, &p_e) in prob_error.iter().enumerate() {
        qurrent_prob_no_error *= 1.0 - p_e;

        if i >= k {
            let p_to_leave = window.pop_front().unwrap();
            qurrent_prob_no_error /= p_to_leave;
        }

        sum_of_expectations += qurrent_prob_no_error;

        if i >= k - 1 {
            window.push_back(1.0 - p_e);
        }
    }
    println!("Quality string: {}",(quality_string.len() - k + 1)as f64);
    println!("SoE: {}",sum_of_expectations);
    (quality_string.len() - k + 1) as f64 - sum_of_expectations
}


fn calculate_error_rate(qual: &str, d_no_min: &[f64; 128]) -> f64 {
    let mut poisson_mean = 0.0;
    let mut total_count = 0;

    for char_ in qual.chars().collect::<HashSet<_>>() {
        let count = qual.chars().filter(|&c| c == char_).count();
        let index = char_ as usize;
        poisson_mean += count as f64 * d_no_min[index];
        total_count += count;
    }

    poisson_mean / total_count as f64

}




fn get_sorted_entries(mini_map_filtered: FxHashMap<i32, Vec<structs::Minimizer_hashed>>)->Vec<(i32, Vec<structs::Minimizer_hashed>)>{
    // Sort by the length of vectors in descending order
    let mut sorted_entries: Vec<(i32, Vec<structs::Minimizer_hashed>)> = mini_map_filtered
        .into_iter()
        .collect();
    sorted_entries.sort_by(|a, b| {
            b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0))
        });


    sorted_entries
}


fn filter_fastq_records(mut fastq_records:Vec<FastqRecord_isoncl_init>,d_no_min:[f64;128], q_threshold:f64,k:usize,d :[f64;128])->Vec<FastqRecord_isoncl_init>{
    fastq_records.par_iter_mut().for_each(|fastq_record| {
    //calculate the error rate and store it in vector errors
        if fastq_record.sequence.len() > k {
            fastq_record.set_error_rate(calculate_error_rate(fastq_record.get_quality(), &d_no_min));
            let exp_errors_in_kmers = expected_number_errornous_kmers(fastq_record.get_quality(), k, &d);
            let p_no_error_in_kmers = 1.0 - exp_errors_in_kmers / (fastq_record.get_sequence().len() - k + 1) as f64;
            //calculate the final score and add it to fastq_record (we have a dedicated field for that that was initialized with 0.0)
            fastq_record.set_score(p_no_error_in_kmers * ((fastq_record.get_sequence().len() - k + 1) as f64))
        }
    });
    //filter out records that have a too high error rate
    fastq_records.retain(|record| 10.0_f64*-record.get_err_rate().log(10.0_f64) > q_threshold);
    println!("Nr of records that passed the filtering: {}",fastq_records.len());
    fastq_records
}


fn parse_cli(k:usize ,w:usize,s:usize,t:usize,quality_threshold:f64, cli:Cli){

}


fn convert_cl_id(v: usize) -> Option<i32> {
    i32::try_from(v).ok()
}

#[derive(Parser,Debug)]
#[command(name = "isONclust3")]
#[command(author = "Alexander J. Petri <alexjpetri@gmail.com>")]
#[command(version = "0.0.1")]
#[command(about = "Clustering of long-read sequencing data into gene families", long_about = "isONclust is a tool for clustering either PacBio Iso-Seq reads, or Oxford Nanopore reads into clusters, where each cluster represents all reads that came from a gene." )]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[arg(long, short, help="Path to consensus fastq file(s)")]
    fastq: String,
    #[arg(long, short,help="Path to initial clusters (stored in fastq format)")]
    init_cl: Option<String>,
    #[arg(short,  help="Kmer length")]
    k: Option<usize>,
    #[arg(short, help=" window size")]
    w: Option<usize>,
    #[arg(short, help=" syncmer length")]
    s: Option<usize>,
    #[arg(short, help=" minimum syncmer position")]
    t: Option<usize>,
    #[arg(long, short, help="Path to outfolder")]
    outfolder: String,
    #[arg(long,short,default_value_t= 1, help="Minimum number of reads for cluster")]
    n: usize,
    #[arg(long, short,help="Path to gff3 file (optional parameter)")]
    gff: Option<String>,
    #[arg(long,help="we do not want to use canonical seeds")]
    noncanonical: Option<bool>,
    #[arg(long,help="Run mode of isONclust (pacbio or ont")]
    mode: String,
    //TODO:add argument telling us whether we want to use syncmers instead of kmers, maybe also add argument determining whether we want to use canonical_minimizers
    #[arg(long,help="seeding approach we choose")]
    seeding: Option<String>,
    #[arg(long,help="quality threshold used for the data (standard: 0.9) ")]
    quality_threshold: Option<f64>,
    #[arg(long,help="print additional information")]
    verbose: Option<bool>,
}

fn main() {

    //INITIALIZATION
    let cli = Cli::parse();
    println!("n: {:?}", cli.n);
    println!("outfolder {:?}", cli.outfolder);

    let mode = cli.mode;
    let n = cli.n;
    let mut k;
    let mut w;
    let mut s;
    let mut t;
    let mut quality_threshold;

    //right now we only have two modes( custom settings for variables k, w, s, and t: 'ont' for reads with  3% error rate or more and 'pacbio' for reads with less than 3% error rate)
    if mode=="ont"{
        k = 13;
        w = 22;//->22, standard: 20
        quality_threshold = 0.9_f64.powi(k as i32);
        s = 9;
        t = 2;
    }
    else if mode == "pacbio"{
        k = 15;
        w = 50;
        quality_threshold = 0.97_f64.powi(k as i32);
        s = 9;
        t = 3;
    }
    else {
        if cli.quality_threshold.is_some(){
            let qt=cli.quality_threshold.unwrap();
            if cli.k.is_some(){
                k = cli.k.unwrap();
            }
            else{
                k=0;
            }
            w=0;
            t=0;
            s=0;
            quality_threshold = qt.powi(k as i32);

        }
        else { panic!("Please set the quality_threshold") }
    }
    let verbo = cli.verbose;
    let mut verbose = false;
    if let Some(verb) = verbo{
        verbose = true;
    }
    let noncan = cli.noncanonical;
    let mut noncanonical_bool= false;
    if let Some(noncanonical)= noncan{
        noncanonical_bool = true;
    }
    let seeding_input = cli.seeding.as_deref();
    let mut seeding= "minimizer";
    if let Some(seed) = seeding_input {
        seeding = seed;
    }
    if seeding =="syncmer"{
        if cli.s.is_some(){
            s = cli.s.unwrap();
            if k-s+1%2!=0{
                panic!("Please set k and s so that (k-s)+1 yields an odd number")
            }
        }
        if cli.t.is_some(){
            t = cli.t.unwrap();
            if (k-s)/2!=t{
                panic!("Please set k,s and t to fulfill (k-s)/2=t")
            }
        }
    }
    else if seeding =="minimizer" {
        if cli.w.is_some(){
            w = cli.w.unwrap();
            if w<k{
                panic!("Please set w greater than k")
            }
        }
    }
    if cli.k.is_some(){
        k = cli.k.unwrap();
    }
    println!("k: {:?}", k);
    println!("w: {:?}", w);
    println!("s: {:?}", s);
    println!("t: {:?}", t);
    //let k = cli.k;
    let outfolder = cli.outfolder;
    let gff_path = cli.gff.as_deref();
    //makes the read  identifiable and gives us the possibility to only use ids during the clustering step
    let mut id_map = FxHashMap::default();
    let mut clusters: FxHashMap<i32, Vec<i32>> = FxHashMap::default();


    let mut annotation_based= false;

    println!("Using {}s as seeds",seeding);
    let now1 = Instant::now();
    {//main scope (holds all the data structures that we can delete when the clustering is done
        //holds the mapping of which minimizer belongs to what clusters
        let mut shared_seed_info: FxHashMap<i32,i32>=FxHashMap::default();
        let mut cluster_map: FxHashMap<u64, Vec<i32>> = FxHashMap::default();
        let initial_clustering_path = cli.init_cl.as_deref();
        if gff_path.is_some(){
            gff_handling::gff_based_clustering(gff_path, initial_clustering_path, &mut clusters, &mut cluster_map, k, w,seeding,s,t, noncanonical_bool);
            println!("{} s used for parsing the annotation information", now1.elapsed().as_secs());
            print!("{:?}",clusters);
            annotation_based = true;
        }
            //let initial_clustering_path = &cli.init_cl.unwrap_or_else(||{"".to_string()});

        //let noncanonical= cli.noncanonical.as_deref();
        //let mut noncanonical=false;
        //if noncanonical.is_some(){
        //    noncanonical=true;
        //}


        //let mut cl_name_map = FxHashMap::default();
        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


        //GENERATION OF ANNOTATION-BASED CLUSTERS


        //cl_id is used to appoint a cluster id to a cluster
        let mut cl_id = 0;

        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


        //SORTING STEP

        let q_threshold = 7.0;
        //path to the sorted file
        let filename = outfolder.clone() + "/clustering/sorted.fastq";
        //count the number of reads that were too short to be clustered
        let mut skipped_cter = 0;
        //d_no_min contains a translation for chars into quality values
        let d_no_min = generate_sorted_fastq_new_version::compute_d_no_min();
        println!("{}", filename);
        let now2 = Instant::now();
        generate_sorted_fastq_for_cluster::sort_fastq_for_cluster(k, q_threshold, &cli.fastq, &outfolder, &quality_threshold, w, seeding,s,t,noncanonical_bool);
        println!("{} s for sorting the fastq file", now2.elapsed().as_secs());
        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


        let now3 = Instant::now();
        //CLUSTERING STEP
        println!("{:?}", clusters);
        {//Clustering scope ( we define a scope so that variables die that we do not use later on)
            //the read id stores an internal id to represent our read
            let mut read_id = 0;
            //this gives the percentage of high_confidence seeds that the read has to share with a cluster to be added to it
            let mut min_shared_minis = 0.8;

            /*if annotation_based{
                min_shared_minis=0.6;
            }*/
            println!("{}", min_shared_minis);
            //parse the file:
            let mut reader = fastq::Reader::from_file(Path::new(&filename)).expect("We expect the file to exist");
            for record in reader.records().into_iter(){
                let seq_rec = record.expect("invalid record");
                let header_new = seq_rec.id();
                if verbose{
                    println!("ID: {}",header_new);
                }

                let sequence = seq_rec.seq();
                let quality = seq_rec.qual();
                //add the read id and the real header to id_map
                id_map.insert(read_id, header_new.to_string());
                if sequence.len() > k {
                    let mut this_minimizers = vec![];
                    let mut filtered_minis = vec![];

                    if seeding == "minimizer"{
                        if noncanonical_bool{
                            generate_sorted_fastq_new_version::get_kmer_minimizers_hashed(&sequence.clone(), k, w, &mut this_minimizers);
                        }
                        else{
                            generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(&sequence.clone(), k, w, &mut this_minimizers);
                        }

                    }
                    else if seeding =="syncmer"{
                        generate_sorted_fastq_new_version::syncmers_canonical(&sequence.clone(), k, s,t , &mut this_minimizers);
                    }

                    generate_sorted_fastq_new_version::filter_seeds_by_quality(&this_minimizers,  quality, k, d_no_min, &mut filtered_minis, &quality_threshold,verbose);
                    // perform the clustering step
                    clustering::cluster(&filtered_minis, min_shared_minis, &this_minimizers, &mut clusters, &mut cluster_map, read_id, &mut cl_id, &mut shared_seed_info);
                    read_id += 1;
                }
                else {
                    skipped_cter += 1
                }
            }
            println!("Finished clustering");
            println!("{} reads used for clustering",read_id);
            println!("Skipped {} reads due to being too short", skipped_cter);

            println!("{} s for reading the sorted fastq file and clustering of the reads", now3.elapsed().as_secs());

            println!("Starting post-clustering to refine clusters");

            clustering::post_clustering_new(&mut clusters,&mut cluster_map,min_shared_minis);
        }


        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }
        //clustering::post_clustering(&mut clusters, &mut cluster_map);
    }


    //FILE OUTPUT STEP

    let now4 = Instant::now();
    write_output::write_output(outfolder, &clusters, cli.fastq, id_map, n );
    println!("{} s for file output", now4.elapsed().as_secs());
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    println!("{} overall runtime", now1.elapsed().as_secs());
}


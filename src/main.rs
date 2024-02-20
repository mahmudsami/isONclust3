//#![allow(warnings)]

use std::fs::File;
use std::collections::{HashMap, HashSet, VecDeque};
use rayon::prelude::*;
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
mod gtf_handling;
use std::path::{PathBuf, Path};
//mod isONclust;
mod structs;
use crate::structs::{FastaRecord, FastqRecord_isoncl_init};
use std::thread;
extern crate rayon;
extern crate clap;
mod write_output;
use memory_stats::memory_stats;
use rustc_hash::{FxHashMap,FxHashSet};

use std::io::Read;
use bio::io::fastq;
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use bio::io::fastq::FastqRead;





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


/*fn main() {
    let args: Vec<String> = env::args().collect();
    //let fasta_path = &args[1];
    let fastq_path = &args[1];
    let mut initial_clustering_records=vec![];
    if args.len()>1{
        let initial_clustering_path = &args[2];
        let initial_clustering_file = File::open(initial_clustering_path).unwrap();
        initial_clustering_records = file_actions::parse_fasta(initial_clustering_path).unwrap();

    }

    let k=9;
    let window_size = 20;

    //let fasta_file = File::open(fasta_path).unwrap();
    //let fasta_records = FileActions::parse_fasta(fasta_file).unwrap();
    let fastq_file = File::open(fastq_path).unwrap();

    let fastq_records = file_actions::parse_fastq(fastq_file).unwrap();

    //test_minimizer_gens();
    //let mut mini_map: HashMap<i32, Vec<generate_sorted_fastq_new_version::Minimizer>> = HashMap::with_capacity(fastq_records.len());
    let mut mini_map_filtered: HashMap<i32, Vec<generate_sorted_fastq_new_version::Minimizer>> = HashMap::with_capacity(fastq_records.len());

    //let mut idmap: HashMap<&str, i32> = HashMap::with_capacity(fastq_records.len());
    //let mut iterator=fastq_records.iter();
    for fastq_record in fastq_records{
        let this_minimizers = generate_sorted_fastq_new_version::get_kmer_minimizers(&*fastq_record.get_sequence(), k, window_size);
        //mini_map.insert(fastq_record.internal_id, this_minimizers.clone());
        let filtered_minis = generate_sorted_fastq_new_version::filter_minimizers_by_quality(this_minimizers,&fastq_record.sequence, &fastq_record.quality,window_size,k);
        mini_map_filtered.insert(fastq_record.internal_id,filtered_minis);
    }
    //sorted_entries: a Vec<(i32,Vec<Minimizer)> sorted by the number of significant minimizers: First read has the most significant minimizers->least amount of significant minimizers
    let sorted_entries = get_sorted_entries(mini_map_filtered);

    if initial_clustering_records.len() > 0{
        let init_cluster_map= clustering::get_initial_clustering(initial_clustering_records,k,window_size);
        clustering::cluster_from_initial(sorted_entries, init_cluster_map);
    }
    else{
        clustering::cluster_sorted_entries(sorted_entries);
    }
    let d_no_min =generate_sorted_fastq_new_version::compute_d_no_min();
    println!("{:?}",d_no_min)


    //generate_sorted_fastq_for_cluster::sort_fastq_for_cluster(15,7.0,"/home/alexanderpetri/Rust/100k_sample.fastq")
    //let seq=dna!("ACTGACTACACAT");
    //let dna_sequence = "ATCGA";
    //let reversed_complement = generate_sorted_fastq_for_cluster::reverse_complement(dna_sequence);
    //println!("Reverse complement: {}", reversed_complement);

    /*for rec in fasta_records {
        println!("{}", rec);
    }
    for rec in fastq_records {
        println!("{}",rec);
    }*/






    /*
    THREAD test field: here I test how to use threads in Rust to actually make it possible to run several threads in parallel

     */
    /*let handle=thread::spawn(|| {
        for i in 1..10 {
            println!("hi number {} from the spawned thread!", i);
            thread::sleep(Duration::from_millis(1));
        }
    });
    handle.join().unwrap();
    for i in 1..5 {
        println!("hi number {} from the main thread!", i);
        thread::sleep(Duration::from_millis(1));
    }
    let threads: Vec<_> = (0..500)
        .map(|i| {
            thread::spawn(move || {
                println!("Thread #{} started!",i);
                thread::sleep(Duration::from_millis(5000));
                println!("Thread #{} finished!",i);
            })
        })
        .collect();

    for handle in threads {
        handle.join().unwrap();
    }*/

}
*/

/*fn main() {
    let matches = command!() // requires `cargo` feature
        .propagate_version(true)
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("index_kmers")
                .arg(
                    arg!([FASTQ])
                        .value_parser(value_parser!(String))
                        .long("fastq")
                        .help("Path to consensus fastq file(s)")
                )
                .arg(
                    arg!([K])
                        .value_parser(value_parser!(u32))
                        .short('k')
                        .default_value("15")
                        .help("kmer size")
                )
                .arg(
                    arg!([T])
                        .value_parser(value_parser!(u32))
                        .short('t')
                        .default_value("8")
                        .help("kmindex number of threads")
                )
        )
        .get_matches();
        match matches.subcommand() {
        Some(("index_kmers", sub_matches)) => {
        index_kmers(&sub_matches)
        }
    _ => unreachable!("Exhausted list of subcommands and subcommand_required prevents `None`"),
    }

}*/

use std::convert::TryFrom;
use crate::gtf_handling::resolve_gtf;

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
    /*#[arg(short,  help="Kmer length")]
    k: Option<usize>,
    #[arg(short, help=" window size")]
    w: Option<usize>,*/
    #[arg(long, short, help="Path to outfolder")]
    outfolder: String,
    #[arg(long,short,default_value_t= 1, help="Minimum number of reads for cluster")]
    n: usize,
    #[arg(long, short,help="Path to gtf file (optional parameter)")]
    gtf: Option<String>,
    #[arg(long,help="Path to gtf file (optional parameter)")]
    noncanonical: Option<bool>,
    #[arg(long,help="Run mode of isONclust (pacbio or ont")]
    mode: String,
    //TODO:add argument telling us whether we want to use syncmers instead of kmers, maybe also add argument determining whether we want to use canonical_minimizers

}

fn main() {


    //INITIALIZATION


    let cli = Cli::parse();
    //println!("k: {:?}", cli.k);
    //println!("w: {:?}", cli.w);
    println!("n: {:?}", cli.n);
    println!("outfolder {:?}", cli.outfolder);

    let now1 = Instant::now();
    let mode = cli.mode;
    let mut k;
    let mut w;
    let mut quality_threshold;

    if mode=="ont"{
        k = 13;
        w = 20;
        quality_threshold = 0.9_f64.powi(k as i32);
    }
    else{
        k = 15;
        w = 50;
        quality_threshold = 0.97_f64.powi(k as i32);
    }


    println!("k: {:?}", k);
    println!("w: {:?}", w);
    //let k = cli.k;
    let window_size = w;
    let w = window_size - k;
    let outfolder = cli.outfolder;
    //makes the read  identifiable and gives us the possiblility to only use ids during the clustering step
    let mut id_map = FxHashMap::default();
    let mut clusters: FxHashMap<i32, Vec<i32>> = FxHashMap::default();
    let gtf_path = cli.gtf.as_deref();

    {//main scope (holds all the data structures that we can delete when the clustering is done

        let initial_clustering_path = cli.init_cl.as_deref();
        //resolve_gtf(gtf_path,initial_clustering_path,outfolder.clone());
        //let initial_clustering_path = &cli.init_cl.unwrap_or_else(||{"".to_string()});

        //let noncanonical= cli.noncanonical.as_deref();
        //let mut noncanonical=false;
        //if noncanonical.is_some(){
        //    noncanonical=true;
        //}

        //holds the mapping of which minimizer belongs to what clusters
        let mut cluster_map: FxHashMap<u64, Vec<i32>> = FxHashMap::default();
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
        //initial clusters were given (annotation based clustering)-> we generate the initial clusters
        if let Some(clustering_path) = initial_clustering_path {
            if initial_clustering_path.is_some() {
                //parse the fasta file containing the information
                //let mut reader = parse_fastx_file(&clustering_path).expect("valid path/file");
                //iterate over each read int the fasta file
                //while let Some(record) = reader.next(){
                let mut reader = fasta::Reader::from_file(Path::new(&clustering_path)).expect("We expect the file to exist");
                //let mut reader = parse_fastx_file(&filename).expect("valid path/file");
                for record in reader.records().into_iter(){
                    //retreive the current record
                    let seq_rec = record.expect("invalid record");
                    //in the next lines we make sure that we have a proper header and store it as string
                    let header_str = seq_rec.id();
                    //let header_str = match std::str::from_utf8(header) {
                    //    Ok(v) => v,
                    //    Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
                    //};
                    let header_new = header_str.to_string();
                    //retreive the sequence of the read
                    let sequence = seq_rec.seq();
                    //make the cluster identifiable for later use
                    //cl_name_map.insert(identifier, header_new);
                    //we make sure that the sequence is not smaller than k
                    if sequence.len() > k {
                        //this_minimizers stores the minimizers generated for the current cluster read
                        let mut this_minimizers = vec![];
                        //generate the minimizers for this gene family
                        generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(&sequence.clone(), k, window_size, &mut this_minimizers);
                        //fill up the cluster_map, a FxHashMap having the Minimizer hash as key and a vector of  clusters the minimizer was found in as value
                        clustering::generate_initial_cluster_map(&this_minimizers, &mut cluster_map,cl_id);
                        //we add an empty vector to clusters for each of the clusters we found in the annotation file
                        let id_vec= vec![];
                        clusters.insert(cl_id,id_vec);
                        //increase the cl_id
                        cl_id += 1;
                    }

                }
            }

        }
        println!("{} s used for parsing the initial clustering file", now1.elapsed().as_secs());
        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


        //SORTING STEP

        let q_threshold = 7.0;
        //path to the sorted file
        let filename = outfolder.clone() + "/sorted.fastq";
        //count the number of reads that were too short to be clustered
        let mut skipped_cter = 0;
        //d_no_min contains a translation for chars into quality values
        let d_no_min = generate_sorted_fastq_new_version::compute_d_no_min();
        println!("{}", filename);
        let now2 = Instant::now();
        generate_sorted_fastq_for_cluster::sort_fastq_for_cluster(k, q_threshold, &cli.fastq, &outfolder, &quality_threshold, window_size);
        println!("{} s for sorting the fastq file", now2.elapsed().as_secs());
        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


        let now3 = Instant::now();
        //CLUSTERING STEP

        {//Clustering scope ( we define a scope so that variables die that we do not use later on)
            //the read id stores an internal id to represent our read
            let mut read_id = 0;
            //this gives the percentage of high_confidence seeds that the read has to share with a cluster to be added to it
            let min_shared_minis = 0.8;
            //parse the file and do for each read in it:

            let mut reader = fastq::Reader::from_file(Path::new(&filename)).expect("We expect the file to exist");
            //let mut reader = parse_fastx_file(&filename).expect("valid path/file");
            for record in reader.records().into_iter(){
            //while let Some(record) = reader.next() {
                let seq_rec = record.expect("invalid record");
                let header_new = seq_rec.id();
                //let header_str = match std::str::from_utf8(header) {
                //    Ok(v) => v,
                //    Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
                //};
                //let header_new = header_str.to_string();
                let sequence = seq_rec.seq();

                let quality = seq_rec.qual();//.expect("We also should have a quality");
                //add the read id and the real header to id_map
                id_map.insert(read_id, header_new.to_string());
                if sequence.len() > k {
                    let mut this_minimizers = vec![];
                    let mut filtered_minis = vec![];
                    generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(&sequence.clone(), k, window_size, &mut this_minimizers);
                    generate_sorted_fastq_new_version::filter_minimizers_by_quality(&this_minimizers,  quality, k, d_no_min, &mut filtered_minis, &quality_threshold);
                    //perform the clustering step
                    clustering::cluster(&filtered_minis, min_shared_minis, &this_minimizers, &mut clusters, &mut cluster_map, read_id, &mut cl_id);
                    read_id += 1;
                }
                else {
                    skipped_cter += 1
                }

            }
            println!("Finished clustering");
            println!("{} reads used for clustering",read_id);
            println!("Skipped {} reads due to being too short", skipped_cter);
        }



        println!("{} s for reading the sorted fastq file and clustering of the reads", now3.elapsed().as_secs());
        if let Some(usage) = memory_stats() {
            println!("Current physical memory usage: {}", usage.physical_mem);
            println!("Current virtual memory usage: {}", usage.virtual_mem);
        } else {
            println!("Couldn't get the current memory usage :(");
        }


    }


    //FILE OUTPUT STEP


    let now4 = Instant::now();
    write_output::write_output(outfolder, &clusters, cli.fastq, id_map);
    println!("{} s for file output", now4.elapsed().as_secs());
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    println!("{} overall runtime", now1.elapsed().as_secs());
}


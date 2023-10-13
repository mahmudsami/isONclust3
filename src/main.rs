use std::time::Duration;
use std::time::Instant;
use std::fs::File;
use std::collections::HashMap;
use rayon::prelude::*;
//use crate::generate_sorted_fastq_new_version::{filter_minimizers_by_quality, Minimizer,get_kmer_minimizers};
//use clap::{arg, command, Command};
use clap::Parser;
pub mod file_actions;
mod clustering;
mod generate_sorted_fastq_for_cluster;
mod generate_sorted_fastq_new_version;
use std::path::PathBuf;
mod isONclust;
mod structs;
use crate::structs::FastaRecord;
mod write_output;

fn compare_minimizer_gens(){
    //let input = "ATGCTAGCATGCTAGCATGCTAGC";
    let input ="GCTTCAGCTCGCCGCCGCGCTTCGGCCGCGGCCGAGCGGGCGGGAAGGGATCTCTCTAGGAGCTGACACTCGAACCTTCACGACCCATTCGGATTTTT";
    //let input = "GCTTCAGCTCGCCGCCGCGCTTCGGCCGCGGCCGAGCGGGCGGGAAGGGATCTCTCTAGGAGCTGACACTCGAACCTTCACGACCCATTCGGATTTTTCCAGGACTCAGGAGTGGCACTGGGAAGAAGGGGACCGCTTCTGCAATTGGCCTCGACACTGGCTGCCAAGAAGACCTGTCGCCTTTGTTTTTAAGTCTCCAGAAATGGAAGAAGAAGCGAAGTCAGTTGAAGTCACGAGAAATCAGCGAGGCATTTGAAGGCGCTTCCTGAAAACCTTTATTCCCTGGCACGTTGTCTGTTTCAACAGCCCTGCCCCTCCTCGGAGCCTGCTTGTGGAATTCTTCCCCTTCGGGTGTGTGGTGGCATTCCCGCCACGTCCAATGTGGACTCCAAAGATTTGTCCCTTCTTTGTCATTTGAAATGAAATGATGGCCACCGCGTCGGACTGGTCTGTCTGAGGGAGATGGTGACAAGCTCAAGGCCTGCGAGGAAAACAAAGAAACCGCTTGGAACCAATGGATACCATTTTGTTAAGCAAGTTAAAGAGGAGGACCTGCTTTTGAAGCTGGATTATGTACAGGTGACCGAATTATAAAAGTCATGGAGAAAGTGTTATTGGCAAAACCTATTCCCAAGTAATTGCTTTAATTCAAACAGTGATACAACATTGGAACTTAGTGTTATGCCAAAAGATGAAGACATTCTCCAAGTGGCATATTCTCAAGATGCCTACCTGAAAGGCAACGAAGCTTATAGCGGCAATGCCCGCAATATACCTGAACCTCCACCAATCTGCTATCCTGGCTGCCATCTGCCCCATCAGCCATGGCACAGCCAGTTGAAATATCTCCTCCTGACTCATCATTGAGCAAACAGCAAACCAGTACACCAGTACTGACACAACCTGGTAGGGCCTATAGAATGGAAATACAAATGCCTCCATCACCAACAGATGTTGCAAATTCAAACACAGCAGTGTGTGTTTGCAATGAAAGTGTAAGGACTGTCATTGTGCCTTCTGAGAAGGTTGTAGATTTGTTATCCAATAGAAACAACCATACAGGTCCTTCACATAGAACTGAAGAAGTGAGGTATGGCGTGAGTGAGAGCAGACCTCTTTAAAAACAGTGTCAAGAACCACATCACCACCATTATCAATTCCCACCACTCATCTAATTCATCAGCCTGCAGGCTCCAGATCACTGGAACCTTCTGGAATTTTACTTAAGTCTGGAAATTACAGTGGACATTCTGATGGAATCTCAAGCAGCAGATCTCAAGCTGTGGAGGCTCCCTCTGTATCTGTTAATCACTATTCGCCAAATTCCCATCAGCACATAGACTGGAAAAACTATAAAACTTACAAAGAGTATATTGATAACAGACGATTGCACATAGGTTGTCGGACAATACAAGAAAGATTAGATAGTTTAAGAGCAGCATCTCAAAGCACGACAGATTATAACCAGTCGTCCCCAACCGCACTACTTTGCAGGGACGACGTCGAAGCACCTCTCATGATCGAGTGCCCAGTCTGTCCAGATACGGCAACGCAGTGTGTCCCAAGAAAGACTGGAAGATTCTGTGCTAATGAAGTATTGTCCAAGAAGTGCATCTCAAGGAGCACTGACGTCTCCATCTGTTAGTTTTAGTAATCATAGAACTCGTTCATGGGATTATATTGAGGGACAGGATGAAACCTTAGAAAATGTCAATTCTGGAACCCAATACCTGATTCCAATGGAGAGAAAAAACAGACTTACAAGTGGAGTGGGTTTACTGAACAGGATGATAGACGAGGTATTTGTGAAAGACCTAGGCAGCAAGAAATTCATAAATCTTTTCGAGGTTCCAATTTTACTGTGGCTCCAAGCGTTGTTAATTCTGATAACAGGCGAATGAGTGGTAGAGGAGTGGGATCTGTGTCGCAGTTTAAAAAAATTCCACCAGATCTAAAAACATTGCAGTCAAACAGAAATTTTCAGACTACTTGTGGAATGTCACTGCCTCGGGTATTTCACAAGACAGGTCACCTCTTGTGAAAGTCCGAAGTAATTCTCTGAAAGCTCCTTCCACGCATGTCACAAAACATCATTTAGCCAGAAATCATTTGTTTCTATCAAAGACCAAAGACCAGTAAATCACTTGCATCAGAACAGTCTGTTGAATCAGCAGACATGGGTAAGGACTGACAGTGCCCCCGATCAGCAAGTGGAGACTGGGAAATCCCCCTCTTTATCTGGAGCCTCTGCCAAGCCTGCCCTCAGTCGAGTGAAAACGCTGGTACTTCAGATTTAGAACTACTGTCAGTCAAAGGAATCAAGATTTAAGTTTACAAGAGGCTGAAACTGAGCAATCAGATACTTTAGATAATAAAGAAGCTGTCATCCTAAGGGAAAAAACCTCCATCTGGACGCCAGACACCGCAGCCTTTAAGGCATCAGTCTTACATCTTGGCAGTAAATGACCAGGAGACCGGGTCAGACACTACCTGCTGGCTGCCCAATGAAGCACGTCGAGAGGTCCACATAAAAAGAATGGAGGAAAAAAAGCCTCGAGTACCAGTCCGCCTGGCGATTCTTTGGCTTCCATCCCATTTATAGATGAACCAACTAGCCCTAGCATTGATCATGATATTGCACATATCCCTGCCTCTGCTGTTATATCAGCCTCTACCTCTCAGGTCCCCTCCATAGCAACAGTTCCTCCTTGCCTCACAACTTCAGCTCCATTAATTCGCCGTCAGCTCTCACATGACCACGAATCTGTTGGCCCTCCTAGCCTGGATGCTCAGCCCAACTCAAAGACAGAAAGATCAAAATCATATGATGAGGGTCTGGATGATTACAGAGAAGATGCAAAATTGTCCTTTAAGCACGTATCTAGTCTGAAGGGAATCAAGATCGCAGACAGCCAAAAGTCATCAGAAGACTCTGGGTCCAGAAAAGATTCTTCCTCAGAGGTCTTCAGTGATGCTGCCAAGGAAGGGTGGCTTCATTTCCGACCCCTTGTCACCGATAAGGGCAAGCGAGTTGGTGGAAGTATTCGGCCATGAAACAGATGTATGTTGTCCTTCGGGGTCATTCACTTTACCTGTACAAAGATAAAAGAGAGCAGACGACTCCGTCTGAGGAAGAGCAGCCCATCAGTGTTAATGCTTGCTTGATAGACATCTCTTACAGTGAGACCAAGAGGAAAAATGTGTTTCGACTCACCACGTCCGACTGTGAATGCCTGTTTCAGCTGAAGACAGAGATGATATGCTAGCTTGGATCAAGACGATCCAGGAGAGCAGCAACCTAAACGAAGAGGACACTGGAGTCACTAACAGGGATCTAATTAGTCGAAGAATAAAAGAATACAACAATCTGATGAGCAAAGCAGAACAGTTGCCAAAAACACCTCGCCAGAGTCTCAGCATCAGGCAAACTTTGCTTGGTGCTAAATCAGAGCCAAAGACTCAAAGCCCACACTCTCCGAAGGAAGAGTCGGAAAGGAAACTTCTCAGTAAAGATGATACCAGTCCCCCAAAAGACAAAGGCACATGGAGAAAAGGCATTCCAAGTATCATGAGAAAGACATTGAGAAAAAGCCAACTGCTACAGGAACTTTCGGCGTCCGACTAGATGACTGCCACACAGCTCATACTAATCGGGTAGGAGACACCAGCCCACCGGATTTCCTCAAGTCATTAATTCAGTTTGATTATAATGGGATCTGAGTATATAGCCTTTCAGACAACAGAGACTGTCAATCTGTTTTCCTTTTCATAGCAGACCTTAGTTTCCCTAAGTATTTCTGACCTGTAATGCTTTTTTTCTAATGTATTGACTCCGAAATTGAATAATTATTCAACCATTAGGAATTTGTCTGAAATGGTCTGTAAGAATTCGTCTTCAAATCGTGTGTTTCTGCTTTTAGCATTTTAAAACAAGGATGAGAACAAGGACCTGATTAAGAAAGTATAATACATGTAAATGTATGTGCATATATTTAAAACACCCATGCTAGGATTCGTATCATGAATTGATATTACTCTCTCCTCCTCTTTTCTTGACTCTTTATTCAGCGGAACGGTTTCCACATTCTCTTTACCTGTTGTTACTTTCCTTTTTAGAATGTTATCTCTCTGGGGGTCCCTCTCTTTCTTTCAGTCTCTTTGTTCTCTTCTTACTGTCTCAAATATCTCTTGCAGGTTTTATCTATTCCCAGATCTATAATGTCTGTATTCTTTATATTACATAACAGCAATTTATGTTTCCTTTTGGATGAAAGGGGCATTTTTGCTGTGATATACAACTACTAAACATTTGTGCTGCATTAAATAAAATGTTTATTTCTTGTTTTGCAGTATATTCCATTAATAGTTGACATATGTTGCAAATTAGTTGAAGAAAGAGGTCTTGAATATACAGGTATTTATAGAGTTCCTGGAAATAATGCAGCCATCTCAAGTATGCAAGAAGAACTCAACAAGGGAATGGCTGATATTGATATACAAGATGATGTAAGTTTGCTTATTTGAATAACTGATCCTTTTTTTTCCCCCTTTGGTGGAACTGGATATGAAAATACTATATTTCTATTTTGTACACCTAGAAATGGCGAGATTTGAATGTGATAAGCAGTTTACTAAAATCCTTCTTCAGAAAACTCCCTGAGCCTCTCTTCACAAATGGTGAGTTAACTCCCGTGAAAGATACCAGATGAAATGGTTCAAGTTACAATTGGAAAACCATAGAATATTAAATTGAACAAGGACAAAATGTCTTGGTCAAAGGAAGCAGGAGGCTTTATTTAGTTCTCCCTCATATATTAAAAACATAAAAGTTATTTAGGTTCAGTGAAATATGTATTATCTTGTTCATACATGTTTGCAGCCAAATTAATATTATTGTTTCCAGTCAAATGGCTTCTGACATCATAAGGTACAATGAGTAGAAAAATGTTTTTTCAGTCAAGCCAAGTATTTAAGTGCCATGTAAGACACTTTACTAGATGTTGTTGGAGGGTGGGGGAGGAATGGACAGGCTGTGTAAATATAGGATGGTTCCTGCAGGGAACTTACATTTTAGAAATTTACATTTCAAAACTATTAGAACTCTTTCCTATGTTGTGGTTTAAAATTGATCCTTAAAGGGAAACTAGATTACAATTCACCTAGGCAATTTAAGATAAATTTGTTTTATTGTTCTTGTCGTCAGATAAATATGCTGATTTTATTGAAGCCAATCGTAAAGAAGATCCTCTAGATCGTCTGAAAACATTAAAAAGACTACCTCTGAATCCATTGGAAAATAATCTGTTTCAGATTCACGATTTGCCTGAACATCATTATGAAACACTTAAGTTCCTTTCAGCTCATCTGAAGACAGTGGCAGAAAATTCAGAAAAAAATAAGATGGAACCAAGAAACCTAGCAATAGTGTTTGGTCCCACCCTTGTTCGAACATCAGAAGACAACATGACCCACATGGTCACCCACATGCCTGACCAGTACAAGATTGTAGAAACGCTCATCCAGCACGTAAGTTCATGCTTCATCTCTGTAAGTTCATGCTTCATCTCTGCAAGAAAATTCATTTTCACGGGCAGCCGCTATTTGGTCAGACAGAGGTTGGTTGTAATATCTTTCCCTAAGAAAATTACCCTTTTTCTATACTATCCTCCAAGATCTCAGTATATTTAAAGGCTTGTTTAATTTCTCTTCTCACCTTTTTTAAATTTTAGCATGACTGGTTTTTCACAGAAGAAGGTGCTGAAGAGCCTCTTGTAAGTATTGTCTGTCAGCATTTGTACTCAGTAGTTTCATTTGGGACACAAGATCCTGGCGCGCCAACAATGTCAACTTCCTGCTTTGCTATTGTCCGCTTTGCGCCCCGGAAGCAGGTCTCTAGCTCAGTTTCTGATCTAATGAATATACAAGTACGGTGCTTATTTGCAGTTAATTATGCGTAAACAATGTAAACACATTAAGCACAGGTGCTTTACTTAGTAATGTGTCTTTTACTAAATGTTGTTTAAATTTTAAAGAAAGTATCTTTGAAGAGATGATGAAGAGATCTTTGTTTATTAACCATAGACAACAGTGCAGGAGGAAAGCACAGTAGACTCCCAGCCAGTGCCAAACATAGATCATTTACTCACCAACATTGGAAGGACAGGAGTCTCCCCAGGAGATGTATCAGGTAACTCTTGCCTGGGCAGTCTTGACCTCTAGGCTGAAAGCTACATTCTTACAATATGAAAGAGCAAATCAATCATAGTGAATTGTTACTTATGATTTTTTCTTTTTCTCTTCCTTTCTACTAATAAAAAACCTTAATCTTCCTTTTCCAATTGATCTGAATTATCTTCCAGTGACTCTAATGCCATTGGCAGTGCATGCAGTTGTAATGTTTGTCAGAATGGAACTATGAGGGAAATGAGCATGCAAACACTGTTTTGTTCTTTTTTTTTTTTTTTTTTTTTAATTCACAGTTTCCATGTAATAGTCATCACAAGGTAAAAACAGATATATGCTACTGAGTATGTGGTTGCTGTGATGGCAGAAAGTATGATCTGGGCTGTATAAGTAGTTCTTATTCTTGAGCCAGGAGTTCTGTTTAAAAGCAGCAAGTTAAATCTGTTGAACTGGTAATGGGCACCATTTTCATAAATACATTTTCATAAATACAGTGTATTATTTGTCACAGTAAAAAGCCAGCCTATTTGACATTTATTTTCTTCTGTTTCCTGTTGCCACCCTAATGCTTGGTTAAGTATTAGTTCTTCTAGGCTCTTTTGATGAAAACTTTCTCCTGCCAAACTGTTGCACTAATCTAGCATGTTGTTTTACACAATGTTCTCCCTTTGTCTCTAGTTCTTGTTTTAGTCCTAAAGGGTTTTTTCAGCTCTTACTGACACAGGCAAATCTATACTGATGTTGAACAGGCGAAGTGGGAGGAGTCTATCACACGGTAATGTCACTAGTGGATTCTGTGATGTCCCTGATGGACAGCTGGAACTGTAGGAAGAAGGAGGACTGTTGTCCACACTGTAGCTGCCTTGTCTGCAGAAACCCCCGATTCTTGTAAATGCAAACATAGTTGATCTGTGGTGTAAAGTTTCTCATAAGATTATTGTAAACAGATCAATACTACTTTCTTAAAAGTTTTTCTGCTGTTTCTTGTTATTTCTGAAGCAGACCCTGAGAGAAATAGCATGGTATGGTCTACCAGTTCTTAAAGTCATAAAGCATTTAGCCATGTGTTAAAAGGGGCGTGGCATGGGGCTCAGCTCATTAATTCTGAAAACCAATTTCCCTAATGTTGAAAGATGAAACTTCAAAATTGGTATTATCCTGAATATTGCTAGCCCCTTCCTACCCTCTTTATAATTGTATAATCTGAGTTTTGCCAGACAGTAAAATTATGGTCTCCATTTTTCCATTTCCATAACAGCATTTGTTAATAGATAATTGCAAAACTGGTAGCAATTTAACCAAACTAGTCAAAAGGATCCCTCAACTGCACTGCTATGGAACTGGGGTCTGAAGCACAAAGATAAACATTATTAATGTAAAACCAAAGTGACCTGGAAGCAAGATTATCTTAATGCCTCTTTTGCATCCATACGCATACAGCATAATGCTCAGAAAATTTAAGTCTGACAAAAATCAACATAAACATGCCAGTATGAGTGTCTCAATTGTAGCTAGGGCGGCCAACACTCTTCAAAGTTCTACATTGTGTATCAGTCTCAGATACAGGATGTGGTACATTAAAAAAGAAAATTATTACACTGTCCTGAAGCTTCTTTATTGCTAGATTTAATATTTAGAAAGTCTAATAAACACGAGTGTTTCTGTTTTTGTTTAGCACAGGGCTCCTTTGCCATTTCTATCATTGAGCTGGTTCTGCCTACAATTCACAGCATCCTTTCATGGCGGTGGCTTTCTGTGGTGCAGCTTGCATCTTTGTGAGCTCTGACTATCCACTGTGTTCTATTGTTTTGAGTTTTGCTTTCACATTTTGTTCTGTGCCTGTCTTTCCATAAGTTTTGCATTGTTTGTGGAAAGCAGTAAAGATTGGCTAATACTGACCATGTTTTATTTAATTTTTCACAGATTCAGCTACTAGTGACTCAACAAAATCTAAGGCTCAGGTGATCCTCTGACCTCGGCCTCCCAAGTAGCCAGGACCACAGGGTTCTTGGGGATCTGGAAAGGATCAGTATAGCAGGGAACTGCTTGTGTCCTCCATCTTTGCAGCTGCTAGTCGCAAGAGGAAGAAGCCGAAAGAAAAAGCACAGCCTAGCAGCTCAGAAGATGAACTGGACAATGTATTTTTTAAGAAAGAAAATGTGGAACAGTGTCACAATGATACTAAAGAGGAGTCCAAAAAAGAAAGTGAGACACTGGGCAGAAAACAGAAGATCATCATTGCCAAAGAAAACAGCACTAGGAAAGACCCCAGCACGAAAAAGATGAAAAGATATCACTAGGAAAAGAGAGCACGCCTTCTGAAAACCCTCACCACCACACAACTCAAAACACAACAAGTCACCAACTCTCAGCTGTCGCTTTGCCATCCTGAAAGAGAGCCCCAGGTCACTTCTGGCACAGAAGTCCTCCCACCTTGAAGAGACAGGCTCTGACTCTGGCACTTTGCTCAGCACGTCTTCCCAGGCCTCCCTGGCAAGGTTTTCCATGAAGAAATCAACCAGTCCAGAAACGAAACATAGCGAGTTTTTGGCCAACGTCAGCACCATCACCTCAGATTATTCCACCACATCGTCTGCTACATACTTGACTAGCCTGGACTCCAGTCGACTGAGCCCTGAGGTGCAATCCGTGGCAGAGAGCAAGGGGGACGAGGCAGATGACGAGAGAAGCGAACTCATCAGTGAAGGGCGGCCTGTGGAAACCGACAGCGAGAGCGAGTTTCCCGTGTTCCCCACAGCCTTGACTTCAGAGAGGCTTTTCCGAGGAAAACTGCAAGAAGTGACTAAGAGCAGCCGGAGAAATTCTGAAGGAAGTGAATTAAGTTGCACGAGGGAAGTTTAACATCAAGTTTAGATAGCCGGAGACAGCTCTTCAGTTCCCATAAACTCATCGAATGTGATACTCTTTCCAGGAAAAAATCAGCTAGATTCAAGTCAGATAGTGGAAGTCTAGGAGATGCCAAGAATGAGAAAGAAGCACCTTCGTTAACTAAAGTGTTTGATGTTATGAAAAAAGGAAAGTCAACTGGGAGTTTACTGACACCCACCAGAGGCGAATCCGAAAAACAGGAACCCACATGGAAAACGAAAATAGCAGATCGGTTAAAACTGAGACCCAGAGCCCCTGCGGATGACATGTTTGGAGTAGGGAATCACAAAGTGAATGCCGAGACTGCTAAAAGGAAAAGCATCCGGCGCAGACATACACTAGGAGGGCACAGAGATGCTACCGAAATCAGCGTTTTGAATTTTTGGAAAGTGCATGAGCAGAGCGGGGAGAGAGAATCTGAACTTTCAGCTGTAAACCGGTTAAAACCAAAATGCTCAGCCCAGGACCTTTCCATCTCAGACTGGCTGGCCAGGGAACGCCTACGCACCAGTACCTCTGACCTTAGCAGAGGAGAAATCGGAGATCCCCAGACAGAGAACCCAAGCACACGAGAAATAGCCACGACCGACACACCTTTGTCTCTTCATTGCAACACAGGCAGTTCTTCCAGCACCTTGGCTTCAACAAACAGGCCCCTTCTTTCCATACCACCACAGTCACCTGACCAAATAAACGGAGAAAGCTTCCAGAACGTGAGCAAAAATGCTAGTTCTGCAGCGAATGCCCAACCTCATAAACTGTCTGAAACCCAGGCAGTAAAGCAGAGTTTCATCCCTGTCTTTAAACTGGGGGTATGTCCACTCTAGCAAGTAAAAAAACTACTGTTACACGTTCCAGTAACTCTGTCAATATTTTCTTGTATCAGAATTGTTATTATGCAGCCTTCATTTGGGCTGGTTTCATCATTTTGCACTGTGAAATAGCTTTACAGTGCATTACTACAGCCAGAAGAACATAATATATATATATATTTAAAAATATATCGGATAGTTGTATACAAATGAGCAAGGTATTTGTTGCAACTTACTACATAGCATATACCCAAAATCACTGAAGAAAATCGCTGGCATCAGTGTGCAGCAAATTTGTTCTTTTGGTTTCATCACTAACAAAGTGCCTCATCATAAAAATACAGTTGGTTTTTAGGGTGCCATATTGTTAAAATTAGATAACTTACTTACATTGAATAAACGAATGCGTTTTATTGGTAACAGATATCATTACATTTACCAGTTTTAACACAGGTGGATACAGAACTTCCATTCTTTAGTCATTCCAGGTGGATCTGAGTTTTATATTCAAACTTTTAATACAGTTTTTGAGTTTTGTGTGACTTGAATTTTTAATCTTTCTGTAAATACGTAACTTAAATGAACATATTAAATGTGTATCTTTCTTCAGATACCAGATTTGATATAATGTTGTAACATAGGTGTGTAGATAGTGGATCCTGGATGGAACTGGCTTCTTTATCGAGAAGAATATAATTCTGCATGAGGACTTAATGAATCCAAACCTGTGTCATGCCTGTGTGCATACCCAATTAAACACTGGAAATAAAAATTGTTTTGGC";

    let window_size = 5;
    let k= 3;
    let now = Instant::now();
    let minimizers = generate_sorted_fastq_new_version::get_kmer_minimizers_efficient(input, k, window_size);
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    //println!("Generated Minimizers: {:?}", minimizers);
    let now = Instant::now();
    let minimizers_ineff = generate_sorted_fastq_new_version::get_kmer_minimizers(input, k, window_size);
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    //println!("Generated Minimizers ineff: {:?}", minimizers_ineff);
    let a = generate_sorted_fastq_new_version::Minimizer { sequence: "AGC".to_string(), position: 24 };
    let b = generate_sorted_fastq_new_version::Minimizer { sequence: "AGC".to_string(), position: 24 };
    assert_eq!(minimizers_ineff, minimizers)

}
fn get_sorted_entries(mini_map_filtered: HashMap<i32, Vec<generate_sorted_fastq_new_version::Minimizer>>)->Vec<(i32, Vec<generate_sorted_fastq_new_version::Minimizer>)>{
    // Sort by the length of vectors in descending order
    let mut sorted_entries: Vec<(i32, Vec<generate_sorted_fastq_new_version::Minimizer>)> = mini_map_filtered
        .into_iter()
        .collect();

    sorted_entries.sort_by_key(|(_, v)| std::cmp::Reverse(v.len()));

    sorted_entries
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
    #[arg(short, default_value_t = 13, help="Kmer length")]
    k: usize,
    #[arg(short, default_value_t = 20)]
    w: usize,
    #[arg(long, short, help="Path to outfolder")]
    outfolder: String,
    #[arg(long,short,default_value_t= 1, help="Minimum number of reads for cluster")]
    n: usize,

}

fn main() {
    let cli = Cli::parse();
    println!("k: {:?}", cli.k);
    println!("t: {:?}", cli.w);
    println!("n: {:?}",cli.n);
    println!("outfolder {:?}",cli.outfolder);
    //
    // READ the files (the initial_clusters_file as well as the fastq file containing the reads)
    //
    let fastq_file = File::open(cli.fastq).unwrap();
    //let initial_clustering_path = &cli.init_cl.unwrap_or_else(||{"".to_string()});
    let  initial_clustering_path =cli.init_cl.as_deref();
    let mut initial_clustering_records=vec![];
    let mut init_clust_rec_both_dir=vec![];
    if initial_clustering_path.is_some(){
        let clustering_path = initial_clustering_path.unwrap();
        initial_clustering_records = file_actions::parse_fasta(clustering_path).unwrap();
        init_clust_rec_both_dir = clustering::add_rev_comp_seqs_annotation(initial_clustering_records);
    }
    let k = cli.k;
    let window_size = cli.w;
    let outfolder = cli.outfolder;
    let fastq_records= file_actions::parse_fastq_old(fastq_file).unwrap();
    //let (fastq_records,id_map) = file_actions::parse_fastq(fastq_file);
    //
    //Generate the minimizers for the initial clusters
    //
    let mut mini_map_filtered: HashMap<i32, Vec<generate_sorted_fastq_new_version::Minimizer>> = HashMap::with_capacity(fastq_records.len());
    //
    // Generate minimizers for the fastq file and filter by significance
    //
    let mut int_id_cter= 0;
    let mut id_map=HashMap::new();
    for fastq_record in &fastq_records{
        id_map.insert(int_id_cter,(*fastq_record.header.clone()).to_string());
        //TODO: add hashmap holding the internal_id and read_id
        if fastq_record.sequence.len()>k{
            let this_minimizers = generate_sorted_fastq_new_version::get_kmer_minimizers(&fastq_record.sequence, k, window_size);
            //mini_map.insert(fastq_record.internal_id, this_minimizers.clone());
            let filtered_minis = generate_sorted_fastq_new_version::filter_minimizers_by_quality(this_minimizers,&fastq_record.sequence, &fastq_record.quality,window_size,k);
            mini_map_filtered.insert(int_id_cter,filtered_minis);
            //println!("{} : {} ",int_id_cter, fastq_record.header);
            int_id_cter += 1;
        }
        else {
            println!("Read too short- skipped {}",fastq_record.header)
        }

    }

    //sorted_entries: a Vec<(i32,Vec<Minimizer)> sorted by the number of significant minimizers: First read has the most significant minimizers->least amount of significant minimizers
    let sorted_entries = get_sorted_entries(mini_map_filtered);
    //
    //Perform the clustering
    //
    let mut clusters:HashMap<i32,Vec<i32>> = HashMap::new();
    if init_clust_rec_both_dir.len() > 0{
        let init_cluster_map= clustering::get_initial_clustering(init_clust_rec_both_dir,k,window_size);
        //println!("{:?}",init_cluster_map);
        clusters = clustering::cluster_from_initial(sorted_entries, init_cluster_map);
        println!("{:?}",clusters);
    }
    else{
        //min_shared_minis: The minimum amount of minimizers shared with the cluster to assign the read to the cluster
        let min_shared_minis= 10;
        clusters = clustering::cluster_sorted_entries(sorted_entries, min_shared_minis);
        //println!("{:?}",clusters);
        //TODO: would it make sense to add a post_clustering? i.e. find the overlap between all clusters and merge if > min_shared_minis
    }
    write_output::write_output(outfolder, clusters,fastq_records, id_map);

}

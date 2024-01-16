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
use std::path::PathBuf;
mod isONclust;
mod structs;
use crate::structs::{FastaRecord, FastqRecord_isoncl_init};
use std::thread;
use crate::clustering::calculate_hash;
extern crate rayon;
extern crate clap;
mod write_output;
use memory_stats::memory_stats;

fn compute_d() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0 as f64).powf(exponent).min(0.79433);
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

    let error_rate = poisson_mean / total_count as f64;
    error_rate
}



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
    let a = structs::Minimizer { sequence: "AGC".to_string(), position: 24 };
    let b = structs::Minimizer { sequence: "AGC".to_string(), position: 24 };
    assert_eq!(minimizers_ineff, minimizers)

}
fn get_sorted_entries(mini_map_filtered: HashMap<i32, Vec<structs::Minimizer_hashed>>)->Vec<(i32, Vec<structs::Minimizer_hashed>)>{
    // Sort by the length of vectors in descending order
    let mut sorted_entries: Vec<(i32, Vec<structs::Minimizer_hashed>)> = mini_map_filtered
        .into_iter()
        .collect();

    sorted_entries.sort_by_key(|(_, v)| std::cmp::Reverse(v.len()));

    sorted_entries
}


fn filter_fastq_records(mut fastq_records:Vec<FastqRecord_isoncl_init>,d_no_min:[f64;128], q_threshold:f64,k:usize,d :[f64;128])->Vec<FastqRecord_isoncl_init>{
    fastq_records.par_iter_mut().for_each(|fastq_record| {
    //calculate the error rate and store it in vector errors
        if fastq_record.sequence.len() > k {
            fastq_record.set_error_rate(calculate_error_rate(&fastq_record.get_quality(), &d_no_min));
            let exp_errors_in_kmers = expected_number_errornous_kmers(&fastq_record.get_quality(), k, &d);
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
    #[arg(long, short,help="Path to gtf file (optional parameter)")]
    gtf: Option<String>,
    #[arg(long, short,help="Path to gtf file (optional parameter)")]
    noncanonical: Option<bool>,
    //TODO:add argument telling us whether we want to use syncmers instead of kmers, maybe also add argument determining whether we want to use canonical_minimizers
    
}

fn main() {
    let cli = Cli::parse();
    println!("k: {:?}", cli.k);
    println!("w: {:?}", cli.w);
    println!("n: {:?}",cli.n);
    println!("outfolder {:?}",cli.outfolder);
    //
    // READ the files (the initial_clusters_file as well as the fastq file containing the reads)
    //
    let now1 = Instant::now();
    let fastq_file = File::open(cli.fastq).unwrap();
    println!("{} s used for file input", now1.elapsed().as_secs());
    //let initial_clustering_path = &cli.init_cl.unwrap_or_else(||{"".to_string()});
    let  initial_clustering_path =cli.init_cl.as_deref();
    //let noncanonical= cli.noncanonical.as_deref();
    //let mut noncanonical=false;
    //if noncanonical.is_some(){
    //    noncanonical=true;
    //}

    let gtf_path = cli.gtf.as_deref();
    //let mut initial_clustering_records: Vec<_>=vec![];
    let mut init_clust_rec_both_dir=vec![];
    if initial_clustering_path.is_some(){
        let clustering_path = initial_clustering_path.unwrap();
        init_clust_rec_both_dir = file_actions::parse_fasta(clustering_path).unwrap();
        //init_clust_rec_both_dir = clustering::add_rev_comp_seqs_annotation(initial_clustering_records);
    }
    let mut gtf_entries=vec![];
    if gtf_path.is_some(){
        let gtf_path_u = gtf_path.unwrap();
        gtf_entries = file_actions::parse_gtf(gtf_path_u).unwrap();
        println!("gtf file parsed")
    }

    /*for gtf_e in gtf_entries{
        println!("{}",gtf_e)
    }*/
    let quality_threshold=7.0;
    let k = cli.k;
    let window_size = cli.w;
    let w = window_size - k;
    let outfolder = cli.outfolder;
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    let fastq_records= file_actions::parse_fastq(fastq_file).unwrap();
    //let (fastq_records,id_map) = file_actions::parse_fastq(fastq_file);
    //
    //Generate the minimizers for the initial clusters
    //
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    let mut mini_map_filtered: HashMap<i32, Vec<structs::Minimizer_hashed>> = HashMap::with_capacity(fastq_records.len());
    let mut mini_map_unfiltered: HashMap<i32, Vec<structs::Minimizer_hashed>> = HashMap::with_capacity(fastq_records.len());
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    //
    // Generate minimizers for the fastq file and filter by significance
    //
    let mut int_id_cter= 0;
    let mut id_map=HashMap::new();
    //count the number of reads that were too short to be clustered
    let mut skipped_cter=0;
    //d_no_min contains a translation for chars into quality values
    let d_no_min=generate_sorted_fastq_new_version::compute_d_no_min();
    let d =compute_d();
    //TODO: only use the actual kmers position nothing surrounding
    let mini_range_len = 2 * (w - 1) + k-1;
    let mini_range_len = k;
    println!("Mini range len: {}",mini_range_len);
    //quality_threshold gives at what point minimizers are too low quality to be used in our algo
    //let filtered_records=filter_fastq_records(fastq_records,d_no_min,quality_threshold,k,d);
    println!("Parsed the files");
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    let now2 = Instant::now();
    println!("{} s before minimizergen", now1.elapsed().as_secs());
    for fastq_record in &fastq_records{
        //println!("int id {}",int_id_cter);
        id_map.insert(int_id_cter,(*fastq_record.header.clone()).to_string());
        if fastq_record.sequence.len() > mini_range_len{
            //let this_minimizers = generate_sorted_fastq_new_version::get_kmer_minimizers(&fastq_record.sequence, k, window_size);
            //let this_minimizers = generate_sorted_fastq_new_version::get_kmer_syncmers(&fastq_record.sequence, k,5,-1);
            let this_minimizers = generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(&fastq_record.sequence, k, window_size);
            mini_map_unfiltered.insert(int_id_cter,this_minimizers.clone());
            //mini_map.insert(fastq_record.internal_id, this_minimizers.clone());
            let filtered_minis = generate_sorted_fastq_new_version::filter_minimizers_by_quality(this_minimizers.clone(),&fastq_record.sequence, &fastq_record.quality,w,k,d_no_min);
            mini_map_filtered.insert(int_id_cter, filtered_minis);
            //println!("{} : {} ",int_id_cter, fastq_record.header);
            int_id_cter += 1;
        }
        else {
            skipped_cter+=1
            //println!("Read too short- skipped {}",fastq_record.header)
        }
    }
    println!("{} s for minimizer gen and filtering of minis", now2.elapsed().as_secs());
    println!("Skipped {} reads due to being too short",skipped_cter);
    if let Some(usage) = memory_stats() {
        println!("Current physical memory usage: {}", usage.physical_mem);
        println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    //sorted_entries: a Vec<(i32,Vec<Minimizer)>, sorted by the number of significant minimizers: First read has the most significant minimizers->least amount of significant minimizers
    let sorted_sign_minis = get_sorted_entries(mini_map_filtered);
    //
    //Perform the clustering
    //
    let mut clusters:HashMap<i32,Vec<i32>> = HashMap::new();
    let now3 = Instant::now();
    //annotation based clustering
    if init_clust_rec_both_dir.len() > 0{
        let init_cluster_map= clustering::get_initial_clustering(init_clust_rec_both_dir,k,window_size);
        //println!("{:?}",init_cluster_map);
        //TODO fix cluster_from initial to work with mini hashs and not the strings
        //clusters = clustering::cluster_from_initial(sorted_sign_minis, init_cluster_map);
        //println!("{:?}",clusters);
    }
    //de novo clustering
    else{
        //min_shared_minis: The minimum amount of minimizers shared with the cluster to assign the read to the cluster
        let min_shared_minis= 0.8;
        //the clustering step
        clusters = clustering::cluster_de_novo(sorted_sign_minis, min_shared_minis, mini_map_unfiltered);
        //println!("{:?}",clusters);
        //TODO: would it make sense to add a post_clustering? i.e. find the overlap between all clusters and merge if > min_shared_minis
    }
    println!("{} s for denovo clustering", now3.elapsed().as_secs());

    if let Some(usage) = memory_stats() {
    println!("Current physical memory usage: {}", usage.physical_mem);
    println!("Current virtual memory usage: {}", usage.virtual_mem);
    } else {
        println!("Couldn't get the current memory usage :(");
    }
    println!("Finished clustering");
    let now4 = Instant::now();
    write_output::write_output(outfolder, clusters,fastq_records, id_map);
    println!("{} s for file output", now4.elapsed().as_secs());
    println!("{} overall runtime", now1.elapsed().as_secs());
}

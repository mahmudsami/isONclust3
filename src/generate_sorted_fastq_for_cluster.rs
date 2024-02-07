use std::path::Path;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use crate::file_actions;
use rayon::prelude::*;
use std::time::Instant;
use std::borrow::Borrow;
use std::collections::VecDeque;
use crate::clustering::reverse_complement;
use crate::structs::FastqRecord_isoncl_init;
use crate::write_output;


//https://doc.rust-lang.org/std/primitive.char.html#method.decode_utf16  for parsing of quality values
fn compress_sequence(seq: &str) -> String {
    //compresses the sequence seq by keeping only the first character of each consecutive group of equal characters. The resulting compressed sequence is stored in the variable seq_hpol_comp.

    let mut seq_hpol_comp = String::new();
    let mut last_char: Option<char> = None;

    for ch in seq.chars() {
        if last_char.is_none() || last_char.unwrap() != ch {
            seq_hpol_comp.push(ch);
        }
        last_char = Some(ch);
    }

    seq_hpol_comp
}


fn expected_number_errornous_kmers(quality_string: &str, k: usize, d:&[f64;128]) -> f64 {
    //computes the expected number of errornous kmers for a read by analysing the quality entry
    let prob_error: Vec<f64> = quality_string.chars().map(|char_| d[char_ as u8 as usize]).collect();
    let mut sum_of_expectations = 0.0;
    let mut current_prob_no_error = 1.0;
    //holds k nucleotides to represent a kmer
    let mut window: VecDeque<f64> = VecDeque::with_capacity(k);
    //iterates over the quality values
    for (i, &p_e) in prob_error.iter().enumerate() {
        //the probability that we do not have an error is mulitplied by the probability that we did not have an error at position p_e
        current_prob_no_error *= 1.0 - p_e;

        if i >= k {
            let p_to_leave = window.pop_front().unwrap();
            current_prob_no_error /= p_to_leave;
        }

        sum_of_expectations += current_prob_no_error;

        if i >= k - 1 {
            window.push_back(1.0 - p_e);
        }
    }
    //
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

//D_no_min = {chr(i) : 10**( - (ord(chr(i)) - 33)/10.0 )  for i in range(128)}
fn compute_d_no_min() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0 as f64).powf(exponent);
    }
    d
}



fn analyse_fastq_and_sort(k:usize, q_threshold:f64, in_file_path:&str)->Vec<FastqRecord_isoncl_init>{
    /*
    Reads, filters and sorts reads from a fastq file so that we are left with reads having a reasonable quality score, that are sorted by score
     */
    //read the fastq file and store the result in fastq_records (a vector of FastqRecord_isoncl_init)
    let fastq_file = File::open(in_file_path).unwrap();
    let mut fastq_records = file_actions::parse_fastq(fastq_file).unwrap();
    println!("{} reads recorded",fastq_records.len());
    //filter fastq_records: We only keep reads having a sequence length>2*k and that do not have a shorter compression than k
    fastq_records.retain(|record| record.get_sequence().len() >= 2*k && compress_sequence(&*record.get_sequence()).len() >= k );
    println!("{} reads accepted",fastq_records.len());
    //compute d_no_min and d, two arrays that we use for the calculations
    let d_no_min=compute_d_no_min();
    let d =compute_d();
    
    //iterate over all fastq_records and calculate error_rate as well as score
    fastq_records.par_iter_mut().for_each(|fastq_record| {
        //calculate the error rate and store it in vector errors
        fastq_record.set_error_rate( calculate_error_rate(&fastq_record.get_quality(), &d_no_min));
        let exp_errors_in_kmers = expected_number_errornous_kmers(&fastq_record.get_quality(), k, &d);
        let p_no_error_in_kmers = 1.0 - exp_errors_in_kmers / (fastq_record.get_sequence().len() - k +1) as f64;
        //calculate the final score and add it to fastq_record (we have a dedicated field for that that was initialized with 0.0)
        fastq_record.set_score( p_no_error_in_kmers  * ((fastq_record.get_sequence().len() - k +1) as f64))
    });
    //filter out records that have a too high error rate
    fastq_records.retain(|record| 10.0_f64*-record.get_err_rate().log(10.0_f64) > q_threshold);

    println!("{} reads accepted",fastq_records.len());
    //sort the vector fastq_records by scores
    fastq_records.par_sort_by(|a, b| b.get_score().partial_cmp(&a.get_score()).unwrap());
    //fastq_records.reverse();
    fastq_records
}





fn print_statistics(fastq_records:&Vec<FastqRecord_isoncl_init>){
    /*
    Prints the statistics for the resulting file TODO: add median and mean
     */
    let min_e = fastq_records[0].get_err_rate();
    let max_e = fastq_records[fastq_records.len()-1].get_err_rate();
    println!("Lowest read error rate: {}",min_e);
    println!("Highest read error rate: {}",max_e);
    //logfile.write("Median read error rate:{0}\n".format(median_e))
    //logfile.write("Mean read error rate:{0}\n".format(mean_e))
}


pub(crate) fn sort_fastq_for_cluster(k:usize, q_threshold:f64, in_file_path:&str,outfolder: &String ) {

    let now = Instant::now();
    let fastq_records = analyse_fastq_and_sort(k, q_threshold, in_file_path);
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    write_output::write_ordered_fastq(&fastq_records,outfolder);
    //print_statistics(fastq_records.borrow());
}
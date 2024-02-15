use std::path::Path;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use crate::file_actions;
use rayon::prelude::*;
use std::time::Instant;
use std::borrow::Borrow;
use std::collections::VecDeque;
use crate::clustering;
use crate::structs::FastqRecord_isoncl_init;
use crate::write_output;
use crate::generate_sorted_fastq_new_version;
use crate::write_output::path_exists;
use std::fs;


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
/// Generates positional minimizers from an input string.
/// A positional minimizer is the lexicographically smallest substring of a given window size
/// as the window slides through the input string.
///
/// # Arguments
///
/// * `input` - The input string to generate minimizers from.
/// * `window_size` - The size of the sliding window for generating minimizers.
/// * `k` - The length of k-mers to use for generating minimizers.
///
/// # Returns
///
/// A vector containing `Minimizer` structs, each containing the lexicographically smallest
///substring and its starting position in the input string.
pub fn get_kmer_minimizers<'a>(seq: &'a str, k_size: usize, w_size: usize, mut minimizers:  Vec<usize> )  {
    let mut w= 0;
    if w_size > k_size{
        w = w_size - k_size;
    }
    else {
        w = 1;
    }
    let mut window_kmers: VecDeque<(&'a str, usize)> = VecDeque::with_capacity(w + 1);
    if w+ k_size < seq.len() + 1{
        for i in 0..w {
            window_kmers.push_back((&seq[i..i + k_size], i));
        }
    }
    // Initialize the window_kmers deque
    else if seq.len() + 1 >= k_size{
        let short_w = seq.len() + 1 - k_size;
        for i in 0..short_w {
            window_kmers.push_back((&seq[i..i + k_size], i));
        }
    }

    //store the final positional minimizers in a vector
    //let mut minimizers = vec![];
    // Find the initial minimizer (minimizer of initial window)
    let (curr_min, min_pos) = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
    //add the initial minimizer to the vector
    let mini =  *min_pos;
    minimizers.push(mini);
    //we always store the previous minimizer to compare to the newly found one
    let mut prev_minimizer =(*min_pos,*curr_min);
    //iterate further over the sequence and generate the minimizers thereof
    for (i, new_kmer) in seq[w..].as_bytes().windows(k_size).enumerate() {
        let new_kmer_pos = i  + w;
        let new_kmer_str = std::str::from_utf8(new_kmer).unwrap();
        // updating window
        window_kmers.pop_front().unwrap();
        window_kmers.push_back((new_kmer_str, new_kmer_pos));
        // Find the new minimizer
        let (curr_min, min_pos) = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
        //make sure that the minimal string is a new minimizer not just the previously found one
        if  *min_pos != prev_minimizer.0{ //&& *curr_min != prev_minimizer.1 {
            //add the minimizer into the vector and store the minimizer as previously detected minimizer
            let mini =  *min_pos ;
            //println!("minimizer {:?}",mini);
            minimizers.push(mini);
            prev_minimizer=(*min_pos, *curr_min);
        }
    }

}
pub fn is_significant(quality_interval: &str, d_no_min:[f64;128])->bool{
    let mut significance_indicator= false;
    let mut qualities :Vec<f64> = vec![];
    let mut quality = 1.0;
    let mut index;
    let mut q_value;
    let mut probability_error;
    //for each character in quality string:
    for (i, c) in quality_interval.chars().enumerate() {
        index = c as usize;
        //q_value gives the PHRED quality score: i.e. '+' gives us 0.1
        q_value = d_no_min[index];
        //here we get the base call accuracy
        probability_error= 1.0 - q_value;
        //TODO: if we have a position having a worse quality char than '+' maybe we should not let this minimizer be significant
        //if probability_error <0.9{
        //    significance_indicator = false;}
        qualities.push(probability_error);
        quality *= probability_error
    }

    let quality_threshold=0.9_f64.powi(quality_interval.len() as i32);
    //TODO: let quality be dependent on length of quality_interval (e.g. 1*E-len)
    if quality > quality_threshold {
        significance_indicator = true;
    }
    significance_indicator
}

pub fn get_canonical_kmer_minimizers_hashed(seq: &str, k_size: usize, w_size: usize, this_minimizers: &mut Vec<usize>)  {
    //make sure that we have suitable values for k_size and w_size (w_size should be larger)
    let mut w= 0;
    if w_size > k_size{
        w = w_size - k_size;
    }
    //k_size was chosen larger than w_size. To not fail we use every k-mer as minimizer (maybe have an error message?)
    else {
        w = 1;
    }
    //let mut rc_vec=VecDeque::with_capacity(w+1);
    let mut window_kmers: VecDeque<(u64, usize)> = VecDeque::with_capacity(w + 1);
    let mut k_mer_str;
    let mut rc_string;
    let mut forward_hash;
    let mut reverse_hash;
    //we can only get a minimizer if the sequence is longer than w + k_size - 1 (else we do not even cover one full window)
    if w + k_size < seq.len() + 1{
        for i in 0..w {
            k_mer_str = &seq [i..i + k_size];
            rc_string = clustering::reverse_complement(k_mer_str).clone();

            //generate the hashes of the kmers
            forward_hash = generate_sorted_fastq_new_version::calculate_hash(&k_mer_str.to_string());
            reverse_hash = generate_sorted_fastq_new_version::calculate_hash(&rc_string.to_string());
            //we now want to find the canonical minimizer: we only push the smaller k-mer of k_mer_str and rc_String into the window
            if forward_hash <= reverse_hash {
                window_kmers.push_back((forward_hash, i));
            }
            else{
                window_kmers.push_back((reverse_hash, i))
            }

        }
    }
    //println!("kmers in window: {:?}", window_kmers);
    //store the final positional minimizers in a vector
    if !window_kmers.is_empty(){
        // Find the initial minimizer (minimizer of initial window)
        let mut binding=window_kmers.clone();
        let (curr_min, min_pos) = binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap();
        //add the initial minimizer to the vector
        let mut mini = *min_pos ;
        this_minimizers.push(mini.clone());
        //we always store the previous minimizer to compare to the newly found one
        let mut prev_minimizer = mini;
        let mut new_kmer_pos;
        let mut new_kmer_str;
        let mut rc_string;
        let mut forward_hash;
        let mut reverse_hash;
        //iterate further over the sequence and generate the minimizers thereof
        for (i, new_kmer) in seq[w..].as_bytes().windows(k_size).enumerate() {
            new_kmer_pos = i  + w;
            new_kmer_str = std::str::from_utf8(new_kmer).unwrap();
            rc_string = clustering::reverse_complement(new_kmer_str).clone();
            // updating  by removing first kmer from window
            window_kmers.pop_front().unwrap();
            forward_hash=generate_sorted_fastq_new_version::calculate_hash(&new_kmer_str.to_string());
            reverse_hash=generate_sorted_fastq_new_version::calculate_hash(&rc_string.to_string());
            if reverse_hash > forward_hash{
                window_kmers.push_back((forward_hash, new_kmer_pos));
            }
            else {
                window_kmers.push_back((reverse_hash,new_kmer_pos))
            }
            // Find the new minimizer, we need a ds that was cloned from window_kmers to abide ownership rules in rust
            binding = window_kmers.clone();
            let (curr_min, min_pos) = *binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap();
            //make sure that the minimal string is a new minimizer not just the previously found one
            if  min_pos != prev_minimizer{ //&& *curr_min != prev_minimizer.1 {
                //add the minimizer into the vector and store the minimizer as previously detected minimizer
                mini = min_pos ;
                //println!("minimizer {:?}",mini);
                this_minimizers.push(mini.clone());
                prev_minimizer = mini.clone();
            }
            rc_string.clear();
        }
    }
}
pub fn filter_minimizers_by_quality(this_minimizers: &Vec<usize>,fastq_sequence: &str, fastq_quality:&str, w: usize, k: usize, d_no_min:[f64;128], minimizers_filtered: &mut Vec<usize>) {
    //let mut minimizers_filtered = vec![];
    let minimizer_range = w - 1;
    let mut skipped_cter= 0;
    let mut minimizer_pos;
    let mut minimizer_range_start;
    let mut minimizer_range_end;
    let mut significant;
    //println!("Number of minimizers: {}",this_minimizers.len());
    for mini in this_minimizers{
        //println!("{:?}",mini);
        minimizer_pos = *mini;
        minimizer_range_start = 0;
        //set the start of the minimizer_range that we want to inspect
        if minimizer_pos > minimizer_range {
            //minimizer_range_start = minimizer_pos - minimizer_range;
            minimizer_range_start = minimizer_pos;
        }

        minimizer_range_end= fastq_sequence.len();
        if minimizer_pos + minimizer_range + k < minimizer_range_end{
            //minimizer_range_end = minimizer_pos + minimizer_range + k ;
            minimizer_range_end= minimizer_pos + k;
        }
        let qualitiy_interval = &fastq_quality[minimizer_range_start..minimizer_range_end - 1];
        //println!("Quality_interval len {}",qualitiy_interval.len());
        significant= is_significant(qualitiy_interval, d_no_min);
        if significant{
            minimizers_filtered.push(mini.clone())
        }
        else{
            skipped_cter += 1;
        }
    }
    //println!("{} minimizers filtered out due to bad quality", skipped_cter);
    //println!("Length after filter: {}",minimizers_filtered.len());
    //minimizers_filtered
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

    let window_size=20;
    let w=7;
    //iterate over all fastq_records and calculate error_rate as well as score
    fastq_records.par_iter_mut().for_each(|fastq_record| {
        let mut this_minimizers=vec![];
        let mut filtered_minis = vec![];
        //let window_size=
        get_canonical_kmer_minimizers_hashed(fastq_record.get_sequence(), k, window_size,&mut this_minimizers);
        filter_minimizers_by_quality(&this_minimizers,fastq_record.get_sequence(), &fastq_record.get_quality(),w,k,d_no_min,&mut filtered_minis);
        //calculate the error rate and store it in vector errors
        fastq_record.set_error_rate( calculate_error_rate(&fastq_record.get_quality(), &d_no_min));
        //calculate the final score and add it to fastq_record (we have a dedicated field for that that was initialized with 0.0)
        fastq_record.set_score(filtered_minis.len() as f64)
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
    if !path_exists(&outfolder){
        fs::create_dir(outfolder.clone()).expect("We should be able to create the directory");
    }
    write_output::write_ordered_fastq(&fastq_records,outfolder);
    //print_statistics(fastq_records.borrow());
}
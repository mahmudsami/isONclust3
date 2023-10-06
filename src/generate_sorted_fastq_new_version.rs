use std::collections::VecDeque;
use std::ops::Index;
use rayon::prelude::*;
//use crate::file_actions::FastqRecord_isoncl_init;
use std::cmp::max;
use crate::structs::{FastqRecord_isoncl_init, FastaRecord};
//fn get_positional_minimizers(&seq:String,k:u32,w:u32)->(str,u32){
//    let window: VecDeque<u32> = VecDeque::new();
//    OK(mini_seq,mini_pos)
//}


/// Represents a minimizer along with its starting position in the input string.
#[derive(Debug, PartialEq,Clone)]
pub struct Minimizer {
    pub sequence: String,
    pub position: usize,
}
/*
/// Computes the Probability of incorrect base call for the quality scores we receive from the fastq format
/// #Returns:
///            d[i]: Probability of incorrect base call for ith character
///
*/

pub fn compute_d_no_min() -> [f64; 128] {
    let mut d = [0.0; 128];

    for i in 0..128 {
        let chr_i = i as u8 as char;
        let ord_i = chr_i as i8;
        let exponent = -(ord_i - 33) as f64 / 10.0;
        d[i] = (10.0 as f64).powf(exponent);
    }
    d
}
/*
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
///substring and its starting position in the input string.*/
pub fn get_kmer_minimizers<'a>(seq: &'a str, k_size: usize, w_size: usize) -> Vec<Minimizer> {
    let mut w=0;
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
    else{
        if seq.len()+1>=k_size{
            let short_w = seq.len() + 1 - k_size;
            for i in 0..short_w {
                window_kmers.push_back((&seq[i..i + k_size], i));
            }
        }

    }
    //store the final positional minimizers in a vector
    let mut minimizers = vec![];
    // Find the initial minimizer (minimizer of initial window)
    let (curr_min, min_pos) = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
    //add the initial minimizer to the vector
    let mini =Minimizer {sequence: curr_min.to_string(),position: *min_pos };
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
        if  *min_pos !=prev_minimizer.0{ //&& *curr_min != prev_minimizer.1 {
            //add the minimizer into the vector and store the minimizer as previously detected minimizer
            let mini =Minimizer {sequence: curr_min.to_string(),position: *min_pos };
            //println!("minimizer {:?}",mini);
            minimizers.push(mini);
            prev_minimizer=(*min_pos, *curr_min);
        }
    }
    minimizers
}

//calculates the average of  a list of f64s and returns it as f64
fn average(numbers: &[f64]) -> f64 {
    numbers.iter().sum::<f64>()/ numbers.len() as f64
}

//used to detect significant minimizers
pub fn is_significant(quality_interval: &str)->bool{
    let mut significance_indicator= false;
    let mut qualities :Vec<f64> = vec![];
    let mut quality = 1.0;
    //d_no_min contains a translation for chars into quality values
    let d_no_min=compute_d_no_min();
    //for each character in quality string:
    for (i, c) in quality_interval.chars().enumerate() {
        let index = c as usize;
        let q_value = d_no_min[index];
        let probability_error= 1.0 - q_value;
        qualities.push(probability_error);
        quality = quality * probability_error;
    }

    if quality < 0.001{
        significance_indicator = true;
    }
    significance_indicator
}


//filter out minimizers for which the quality of the minimizer_impact range is too bad
pub fn filter_minimizers_by_quality(this_minimizers: Vec<Minimizer>,fastq_sequence: &str, fastq_quality:&str, w: usize, k: usize)-> Vec<Minimizer>{
    let mut minimizers_filtered = vec![];
    let minimizer_range = w - 1;
    //println!("Length of minimizers: {}",this_minimizers.len());
    for mini in this_minimizers{
        //println!("{:?}",mini);
        let minimizer_pos= mini.position;
        let mut minimizer_range_start= 0;
        if minimizer_pos>minimizer_range{
            minimizer_range_start = minimizer_pos - minimizer_range;
        }
        let mut minimizer_range_end = fastq_sequence.len();
        if minimizer_pos+minimizer_range + k < minimizer_range_end{
            minimizer_range_end = minimizer_pos + minimizer_range + k ;
        }
        let qualitiy_interval= &fastq_quality[minimizer_range_start..minimizer_range_end - 1];
        let significant= is_significant(&qualitiy_interval);
        if significant{
            minimizers_filtered.push(mini.clone())
        }

    }
    //println!("Length after filter: {}",minimizers_filtered.len());
    minimizers_filtered
}


//test implementation currently deprecated
pub fn get_kmer_minimizers_efficient<'a>(seq: &'a str, k_size: usize, w_size: usize) -> Vec<Minimizer> {
    let w = w_size - k_size;
    let mut window_kmers: VecDeque<(&'a str, usize)> = VecDeque::with_capacity(w + 1);
    // Initialize the window_kmers deque
    for i in 0..w {
        window_kmers.push_back((&seq[i..i + k_size], i));
    }
    //store the final positional minimizers in a vector
    let mut minimizers = vec![];
    // Find the initial minimizer (minimizer of initial window)
    let kmer_option = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
    let (curr_min,min_pos) = *kmer_option;
    let mini = Minimizer {sequence: curr_min.to_string(),position: min_pos };
    //add the initial minimizer to the vector
    minimizers.push(mini);

    //throw out all elements that are before the prev_minimizer including prev_minimizer itself
    let mut last_mini = minimizers.last().expect("We already have at least one minimizer in minimizers");
    let mut last_mini_pos = last_mini.position;
    let mut pos = window_kmers.front().expect("The window should contain elements");
    let mut stored_value= pos.1;
    while stored_value <= last_mini_pos{
        window_kmers.pop_front();
        if window_kmers.is_empty(){
            stored_value = last_mini_pos+1;
        }
        else {
            pos = window_kmers.front().expect("The window should contain elements still");
            stored_value = pos.1;
        }
    }
    //we always store the previous minimizer to compare to the newly found one
    for (i, new_kmer) in seq[w..].as_bytes().windows(k_size).enumerate() {
        let new_kmer_pos = i  + w;
        let new_kmer_str = std::str::from_utf8(new_kmer).unwrap();

        if new_kmer_str < &last_mini.sequence{
            let mini = Minimizer {sequence: new_kmer_str.to_string(),position: new_kmer_pos };
            //println!("Minimizer found: {:?}",mini);
            //add the minimizer into the vector
            minimizers.push(mini.clone());
            last_mini = minimizers.last().expect("We already have at least one minimizer in minimizers");
            continue
        }

        // updating window
        window_kmers.push_back((new_kmer_str, new_kmer_pos));
        if window_kmers.len() == w {
            // Find the new minimizer
            let (curr_min, min_pos) = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
            let mini = Minimizer {sequence: curr_min.to_string(),position: *min_pos };
            //println!("Minimizer found: {:?}",mini);
            //add the minimizer into the vector
            minimizers.push(mini.clone());
            last_mini = minimizers.last().expect("We already have at least one minimizer in minimizers");
            last_mini_pos = last_mini.position;
            pos=window_kmers.front().expect("The window should contain elements");
            stored_value = pos.1;
            while stored_value <= last_mini_pos{
                window_kmers.pop_front();
                if window_kmers.is_empty(){
                    stored_value = last_mini_pos+1;
                }
                else{
                    pos = window_kmers.front().expect("The window should contain elements still");
                    stored_value = pos.1;
                }
            }
        }
    }
    if !window_kmers.is_empty(){
        let (curr_min,min_pos) = window_kmers.iter().min_by_key(|&&(kmer, _)| kmer).unwrap();
        let mini = Minimizer {sequence: curr_min.to_string(),position: *min_pos };
        //add the initial minimizer to the vector
        minimizers.push(mini);
    }


    minimizers
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_minimizers_0() {
        let input = "ATGCTAGCATGCTAGCATGCTAGC";
        let window_size = 8;
        let k = 3;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "ATG".to_string(), position: 0 },
            Minimizer { sequence: "AGC".to_string(), position: 5 },
            Minimizer { sequence: "ATG".to_string(), position: 8 },
            Minimizer { sequence: "AGC".to_string(), position: 13 },
            Minimizer { sequence: "ATG".to_string(), position: 16 },
            Minimizer { sequence: "AGC".to_string(), position: 21 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_minimizers_1() {
        let input = "CAATTTAAGGCCCGGG";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AATTT".to_string(), position: 1 },
            Minimizer { sequence: "AAGGC".to_string(), position: 6 },
            Minimizer { sequence: "AGGCC".to_string(), position: 7 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }

    #[test]
    fn test_kmer_minimizers_2() {
        let input = "CAAAGTAAGGCCCTCC";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AAAGT".to_string(), position: 1 },
            Minimizer { sequence: "AAGGC".to_string(), position: 6 },
            Minimizer { sequence: "AGGCC".to_string(), position: 7 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_minimizers_3() {
        let input = "CAATGA";
        let window_size = 10;
        let k = 5;
        let actual_minimizers = get_kmer_minimizers(input, k, window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AATGA".to_string(), position: 1 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_average_0(){
        let mut input=vec![];
        input.push(0.5);
        input.push(0.75);
        input.push(0.25);
        let average_res=average(&*input);
        assert_eq!(average_res,0.5);
    }
    #[test]
    fn test_average_1(){
        let mut input=vec![];
        input.push(1.0);
        input.push(2.0);
        input.push(3.0);
        let average_res=average(&*input);
        assert_eq!(average_res,2.0);
    }
}
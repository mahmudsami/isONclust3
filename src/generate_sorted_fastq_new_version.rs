use std::collections::VecDeque;
use std::ops::Index;
use rayon::prelude::*;
//use crate::file_actions::FastqRecord_isoncl_init;
use std::cmp::max;
use crate::structs::{FastqRecord_isoncl_init, FastaRecord};
use crate::clustering::reverse_complement;
use std::borrow::Borrow;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
//fn get_positional_minimizers(&seq:String,k:u32,w:u32)->(str,u32){
//    let window: VecDeque<u32> = VecDeque::new();
//    OK(mini_seq,mini_pos)
//}


/// Represents a minimizer along with its starting position in the input string.
/// TODO: rename to indexer or similar
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
pub fn get_kmer_minimizers<'a>(seq: &'a str, k_size: usize, w_size: usize) -> Vec<Minimizer> {
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


/// Generates positional canonical minimizers from an input string.
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
///
pub fn get_canonical_kmer_minimizers(seq: &str, k_size: usize, w_size: usize) -> Vec<Minimizer> {
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
    let mut window_kmers: VecDeque<(String, usize)> = VecDeque::with_capacity(w + 1);

    //we can only get a minimizer if the sequence is longer than w + k_size - 1 (else we do not even cover one full window)
    if w + k_size < seq.len() + 1{
        for i in 0..w {
            let k_mer_str = &seq [i..i + k_size];
            let rc_string = reverse_complement(k_mer_str).clone();

            //we now want to find the canonical minimizer: we only push the smaller k-mer of k_mer_str and rc_String into the window
            if k_mer_str <= rc_string.as_str() {
                window_kmers.push_back(((*k_mer_str).to_string(), i));
            }
            else{
                window_kmers.push_back((rc_string, i))
            }

        }
    }
    //println!("kmers in window: {:?}", window_kmers);
    //store the final positional minimizers in a vector
    let mut minimizers = vec![];
    if !window_kmers.is_empty(){
        // Find the initial minimizer (minimizer of initial window)
        let mut binding=window_kmers.clone();
        let (curr_min, min_pos) = binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap();
        //add the initial minimizer to the vector
        let mini =Minimizer {sequence: curr_min.to_string(),position: *min_pos };
        minimizers.push(mini.clone());
        //we always store the previous minimizer to compare to the newly found one
        let mut prev_minimizer = mini;
        //iterate further over the sequence and generate the minimizers thereof
        for (i, new_kmer) in seq[w..].as_bytes().windows(k_size).enumerate() {
            let new_kmer_pos = i  + w;
            let new_kmer_str = std::str::from_utf8(new_kmer).unwrap();
            let rc_string = reverse_complement(new_kmer_str).clone();
            // updating  by removing first kmer from window
            window_kmers.pop_front().unwrap();
            if rc_string > new_kmer_str.to_string(){
                window_kmers.push_back(((*new_kmer_str).to_string(), new_kmer_pos));
            }
            else {
                window_kmers.push_back((rc_string,new_kmer_pos))
            }

            // Find the new minimizer, we need a ds that was cloned from window_kmers to abide ownership rules in rust
            binding = window_kmers.clone();
            let (curr_min, min_pos) = binding.iter().min_by_key(|&(kmer, _)| kmer).unwrap().clone();
            //make sure that the minimal string is a new minimizer not just the previously found one
            if  min_pos !=prev_minimizer.position{ //&& *curr_min != prev_minimizer.1 {
                //add the minimizer into the vector and store the minimizer as previously detected minimizer
                let mini =Minimizer {sequence: curr_min.to_string(),position: min_pos };
                //println!("minimizer {:?}",mini);
                minimizers.push(mini.clone());
                prev_minimizer = mini.clone();
            }
        }
    }
    minimizers
}



//calculates the average of  a list of f64s and returns it as f64
fn average(numbers: &[f64]) -> f64 {
    numbers.iter().sum::<f64>()/ numbers.len() as f64
}

///Used to detect significant minimizers by checking the read qualities and estimating the overall quality of the area of the read
/// Input: quality_interval: the quality values of the area we want to check
/// Output: significance_indicator: a bool stating whether the minimizer is significant( true: yes, false: no)
///
pub fn is_significant(quality_interval: &str, d_no_min:[f64;128])->bool{
    let mut significance_indicator= false;
    let mut qualities :Vec<f64> = vec![];
    let mut quality = 1.0;
    //for each character in quality string:
    for (i, c) in quality_interval.chars().enumerate() {
        let index = c as usize;
        //q_value gives the PHRED quality score: i.e. '+' gives us 0.1
        let q_value = d_no_min[index];
        //here we get the base call accuracy
        let probability_error= 1.0 - q_value;
        //TODO: if we have a position having a worse quality char than '+' maybe we should not let this minimizer be significant
        //if probability_error <0.9{
        //    significance_indicator = false;}
        qualities.push(probability_error);
        quality = quality * probability_error;
    }

    let quality_threshold=0.9_f64.powi(quality_interval.len() as i32);
    //TODO: let quality be dependent on length of quality_interval (e.g. 1*E-len)
    if quality > quality_threshold {
        significance_indicator = true;
    }
    significance_indicator
}


//filter out minimizers for which the quality of the minimizer_impact range is too bad
pub fn filter_minimizers_by_quality(this_minimizers: Vec<Minimizer>,fastq_sequence: &str, fastq_quality:&str, w: usize, k: usize, d_no_min:[f64;128])-> Vec<Minimizer>{
    let mut minimizers_filtered = vec![];
    let minimizer_range = w - 1;
    let mut skipped_cter=0;
    //println!("Number of minimizers: {}",this_minimizers.len());
    for mini in this_minimizers{
        //println!("{:?}",mini);
        let minimizer_pos= mini.position;
        let mut minimizer_range_start= 0;
        //set the start of the minimizer_range that we want to inspect
        if minimizer_pos > minimizer_range{
            //minimizer_range_start = minimizer_pos - minimizer_range;
            minimizer_range_start=minimizer_pos;
        }

        let mut minimizer_range_end = fastq_sequence.len();
        if minimizer_pos + minimizer_range + k < minimizer_range_end{
            //minimizer_range_end = minimizer_pos + minimizer_range + k ;
            minimizer_range_end=minimizer_pos+k;
        }
        let qualitiy_interval= &fastq_quality[minimizer_range_start..minimizer_range_end - 1];
        //println!("Quality_interval len {}",qualitiy_interval.len());
        let significant= is_significant(&qualitiy_interval, d_no_min);
        if significant{
            minimizers_filtered.push(mini.clone())
        }
        else{
            skipped_cter= skipped_cter+1;
        }
    }
    //println!("{} minimizers filtered out due to bad quality", skipped_cter);
    //println!("Length after filter: {}",minimizers_filtered.len());
    minimizers_filtered
}




///Method used to generate syncmers from reads
/// INPUT:  seq: a string reference to the original read sequence
///         k_size: The size of the k_mer used
///         s_size: The size of s
///         t: The size of parameter t
///OUtput:  syncmers: A vector storing all syncmers (we use the minimizer struct to store them as essentially the same infos)
pub(crate) fn get_kmer_syncmers(seq: &str, k_size: usize, s_size: usize, t: isize) -> Vec<Minimizer> {
    let w = k_size - s_size;

    // get t, the position of s-mer
    // t is chosen to be in the middle of k-mer for the chosen syncmer
    let mut t = t;
    if t < 0 {
        let diff = k_size as isize - s_size as isize;
        t = if diff % 2 == 0 {
            diff / 2
        } else {
            (diff + 1) / 2
        };
        t -= 1;
    }

    let mut syncmers = Vec::new();
    // get list of all s-mers in the first k-mer

    let mut kmer_smers: VecDeque<&str> = (0..=w).map(|i| &seq[i..i + s_size]).collect();
    for i in 0..seq.len() - k_size {
        // add a new syncmer to the list if its smallest s-mer is at place t
        if kmer_smers.iter().position(|&x| x == *kmer_smers.iter().min().unwrap()) == Some(t as usize) {
            syncmers.push(Minimizer {sequence: (&seq[i..i + k_size]).parse().unwrap(), position: i});
        }
        // move the window one step to the right by popping the leftmost
        // s-mer and adding one to the right
        kmer_smers.pop_front();
        kmer_smers.push_back(&seq[i + k_size - s_size + 1..i + k_size + 1]);
    }

    syncmers
}


/*fn syncmers_canonical(seq: &str, k: usize, s: usize, t: usize) -> Vec<Minimizer> {
    let seq_rc = reverse_complement(seq);
    let mut hasher = DefaultHasher::new();
    let mut window_smers_fw: VecDeque<u64> = (0..k - s + 1).map(|i| seq[i..i + s].hash()).collect();
    let mut window_smers_rc: VecDeque<u64> = (0..k - s + 1).map(|i| seq_rc[i..i + s].hash()).collect();

    let curr_min_fw = *window_smers_fw.iter().min().unwrap();
    let curr_min_rc = *window_smers_rc.iter().min().unwrap();

    let pos_min_fw = window_smers_fw.iter().position(|&x| x == curr_min_fw).unwrap();
    let pos_min_rc = window_smers_rc.iter().position(|&x| x == curr_min_rc).unwrap();

    let (pos_min, seq_tmp) = if curr_min_fw < curr_min_rc {
        (pos_min_fw, &seq[0..k])
    } else {
        (pos_min_rc, &seq_rc[0..k])
    };

    let mut syncmers:Vec<Minimizer> = vec![];
    if pos_min == t {
        syncmers.push(Minimizer{sequence:seq_tmp.to_string(), position: 0});
    }

    for i in k - s + 1..seq.len() - s {
        seq[i..i + s].hash(&mut hasher);
        let new_smer_fw = hasher.finish();
        seq_rc[i..i + s].hash(&mut hasher);
        let new_smer_rc=hasher.finish();
        // Updating windows
        let _ = window_smers_fw.pop_front();
        window_smers_fw.push_back(new_smer_fw);
        let _ = window_smers_rc.pop_front();
        window_smers_rc.push_back(new_smer_rc);

        let curr_min_fw = *window_smers_fw.iter().min().unwrap();
        let curr_min_rc = *window_smers_rc.iter().min().unwrap();

        let pos_min_fw = window_smers_fw.iter().position(|&x| x == curr_min_fw).unwrap();
        let pos_min_rc = window_smers_rc.iter().position(|&x| x == curr_min_rc).unwrap();

        let (pos_min, seq_tmp) = if curr_min_fw < curr_min_rc {
            (pos_min_fw, &seq[i - (k - s)..i - (k - s) + k])
        } else {
            (pos_min_rc, &seq_rc[i - (k - s)..i - (k - s) + k])
        };

        if pos_min == t {
            let kmer = seq_tmp.to_string();
            syncmers.push(Minimizer {sequence:kmer, position:i - (k - s)});
        }
    }

    syncmers
}*/



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
    fn test_canonical_minis_1(){
        let input ="GGGTAACTTTTCA";
        let window_size=12;
        let k =6;
        let actual_minimizers=get_canonical_kmer_minimizers(input,k,window_size);
        println!("Generated Minimizers: {:?}", actual_minimizers);
        let expected_minimizers = vec![
            Minimizer { sequence: "AAAAGT".to_string(), position: 5 },
        ];
        assert_eq!(actual_minimizers, expected_minimizers);
    }
    #[test]
    fn test_kmer_syncmers(){
        let input ="CATTCAGGAATC";
        let k=5;
        let s=2;
        let t=2;
        let retreived_syncmers=get_kmer_syncmers(input,k,s,t);
        let expected_syncmers=vec![
            Minimizer { sequence: "TCAGG".to_string(), position: 3 },
            Minimizer { sequence: "GGAAT".to_string(),position: 6},

        ];
        assert_eq!(retreived_syncmers, expected_syncmers);
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
use std::path::Path;
use std::collections::HashSet;
use std::process::exit;
use std::thread::sleep;
use rayon::prelude::*;
use std::time::Instant;
use std::collections::VecDeque;
use crate::structs::{FastqRecord_isoncl_init, Minimizer_hashed};
use crate::write_output;
use crate::seeding_and_filtering_seeds;
use crate::write_output::path_exists;
//use crate::bio_rust_file_read;
use std::fs;
use bio::io::fastq;
use rustc_hash::FxHashMap;
use minimizer_iter::MinimizerBuilder;


//https://doc.rust-lang.org/std/primitive.char.html#method.decode_utf16  for parsing of quality values
fn compress_sequence(seq: &[u8]) -> String {
    //compresses the sequence seq by keeping only the first character of each consecutive group of equal characters. The resulting compressed sequence is stored in the variable seq_hpol_comp.

    let mut seq_hpol_comp = String::new();
    let mut last_char: Option<&u8> = None;
    for ch in seq{
        if last_char.is_none() || last_char.unwrap()!= ch{
            seq_hpol_comp.push(*ch as char);
        }
        last_char= Some(ch);
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
        //the probability that we do not have an error is multiplied by the probability that we did not have an error at position p_e
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



fn calculate_error_rate(qual: &[u8], d_no_min: &[f64; 128]) -> f64 {
    let mut counts = vec![0; 128];
    let mut total_count = 0;
    let mut poisson_mean = 0.0;

    for &char_byte in qual.iter() {
        counts[char_byte as usize] += 1;
    }

    for (idx, cnt) in counts.iter().enumerate()  {
        poisson_mean += *cnt as f64 * d_no_min[idx];
        total_count += *cnt;
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


fn analyse_fastq_and_sort(k:usize, q_threshold:f64, in_file_path:&str, quality_threshold: &f64, window_size: usize, score_vec: &mut Vec<(i32,usize)>, id_map: &mut FxHashMap<i32,String>, seeding: &str, s: usize, t: usize, noncanonical_bool: bool){
    /*
    Reads, filters and sorts reads from a fastq file so that we are left with reads having a reasonable quality score, that are sorted by score
     */
    let d_no_min= compute_d_no_min();
    //read_id holds the internal id we appoint to a read
    let mut read_id = 0;
    //generate a Reader object that parses the fastq-file (taken from rust-bio)
    let mut reader = fastq::Reader::from_file(Path::new(&in_file_path)).expect("We expect the file to exist");
    //make sure that we have suitable values for k_size and w_size (w_size should be larger)
    let mut w;
    if window_size > k{
        w = window_size - k+ 1; // the minimizer generator will panic if w is even to break ties
        if w % 2 == 0 {
            w += 1;
        }
    }
    //k_size was chosen larger than w_size. To not fail we use every k-mer as minimizer (maybe have an error message?)
    else {
        w = 1;
    }
    //iterate over the records
    for record in reader.records().into_iter(){
        let seq_rec = record.expect("invalid record");
        let header_new = seq_rec.id();
        let sequence = seq_rec.seq();
        let quality = seq_rec.qual();//.expect("We also should have a quality");
        //add the read id and the real header to id_map
        if sequence.len() >= 2*k && compress_sequence(sequence).len() >= k{
            let mut this_minimizers = vec![];
            let mut filtered_minis = vec![];
            if seeding == "minimizer" {
                if noncanonical_bool{
                    let min_iter = MinimizerBuilder::<u64, _>::new()
                        .minimizer_size(k)
                        .width((w) as u16)
                        .iter(sequence);
                    for (minimizer, position) in min_iter {
                        let mini = Minimizer_hashed { sequence: minimizer, position: position };
                        this_minimizers.push(mini);
                    }
                }
                else {
                    let min_iter = MinimizerBuilder::<u64, _>::new()
                        .canonical()
                        .minimizer_size(k)
                        .width((w) as u16)
                        .iter(sequence);
                    for (minimizer, position, _) in min_iter {
                        let mini = Minimizer_hashed {sequence: minimizer,position: position };
                        this_minimizers.push(mini);
                    }

                    //println!("minimizers NEW len: {:?}", this_minimizers.len());

                    //generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(sequence, k, window_size, &mut this_minimizers);
                    //println!("minimizers OLD len: {:?}", &this_minimizers.len());
                }
            }
            else if seeding =="syncmer"{
                if noncanonical_bool{
                    seeding_and_filtering_seeds::get_kmer_syncmers(sequence, k, s, t, &mut this_minimizers);
                }
                else {
                    seeding_and_filtering_seeds::syncmers_canonical(sequence, k, s, t, &mut this_minimizers);
                }
            }
            else{
                println!("No seeding");
            }
            seeding_and_filtering_seeds::filter_seeds_by_quality(&this_minimizers, quality, k, d_no_min, &mut filtered_minis, quality_threshold, false);
            let error_rate= calculate_error_rate(quality, &d_no_min);
            if 10.0_f64 * - error_rate.log(10.0_f64) > q_threshold{
                id_map.insert(read_id, header_new.to_string());
                let score_tuple=(read_id, filtered_minis.len());
                score_vec.push(score_tuple);
                read_id += 1;
            }
        }
    }
    //sort the score_vec by the number of high-confidence seeds (most to least)
    score_vec.par_sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());//TODO: replace by par_sort_by
    println!("{} reads accepted",score_vec.len());
    //println!("{:?}",score_vec.pop());
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


pub(crate) fn sort_fastq_for_cluster(k:usize, q_threshold:f64, in_file_path:&str, outfolder: &String, quality_threshold:&f64, window_size: usize, seeding: &str, s: usize, t: usize, noncanonical_bool: bool, memory_restriction: bool, p: usize, c: usize){ 
    println!("Sorting the fastq_file");
    println!("Memory restriction: {}", memory_restriction);
    let now = Instant::now();
    //holds the internal ids and scores as tuples to be able to sort properly
    let mut score_vec=vec![];
    //holds the internal read id
    let mut id_map=FxHashMap::default();
    //the main step of the sort_fastq_for_cluster step: Gets the number of high-confidence seeds for each read and writes them into score_vec
    analyse_fastq_and_sort(k, q_threshold, in_file_path,quality_threshold,window_size,&mut score_vec, &mut id_map, seeding,s,t, noncanonical_bool);
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    if !path_exists(outfolder){
        fs::create_dir(outfolder.clone()).expect("We should be able to create the directory");
    }
    //write a fastq-file that contains the reordered reads
    if memory_restriction{
        println!("Memory restriction is active");
        let num_part = p;
        let num_chunks = c;
        
        let total_size = score_vec.len();
        let step = total_size / num_chunks;
        let rem = total_size - (total_size / num_chunks) * (num_chunks-1);

        let (id_maps, score_vecs) = partition_id_map(& id_map, & score_vec, num_part);
        
        for i in 0..num_part{
            write_output::write_ordered_fastq_offset(&score_vecs[i], outfolder, &id_maps[i], in_file_path, num_chunks, i==0, step, rem);
        }
    }
    else{
    write_output::write_ordered_fastq(&score_vec, outfolder,&id_map,in_file_path);
    }
    println!("Wrote the sorted fastq file");
    sleep(std::time::Duration::from_secs(10));
    //print_statistics(fastq_records.borrow());
}

fn partition_id_map( id_map: & FxHashMap<i32,String>, score_vec: & Vec<(i32,usize)>, num_part: usize) -> (Vec<FxHashMap<i32,String>>, Vec<Vec<(i32,usize)>>){
    let mut id_maps = vec![FxHashMap::default(); num_part];
    let mut score_vecs = vec![vec![]; num_part];
    let n = score_vec.len();
    let p = (n + (num_part - n % num_part )) / num_part;
    for (i, (id, score))  in score_vec.iter().enumerate(){
        let part = i / p;
        id_maps[part].insert(*id, id_map.get(id).unwrap().clone());
        score_vecs[part].push((*id, *score));
    }
    (id_maps, score_vecs)
}
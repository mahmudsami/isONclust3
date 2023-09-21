use std::collections::HashSet;
use crate::file_actions::{FastqRecord_isoncl_init, FastaRecord};
use std::collections::HashMap;
use crate::generate_sorted_fastq_new_version::{filter_minimizers_by_quality, Minimizer, get_kmer_minimizers};
use crate::generate_sorted_fastq_new_version;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::num::ParseIntError;

enum Cluster<T, U> {
    read_ids(HashSet<T>),
    mini_seqs(HashSet<U>),
}


//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster_sorted_entries(sorted_entries: Vec<(i32,Vec<Minimizer>)>) -> HashMap<i32,Vec<i32>>{
    let mut clusters: HashMap<i32,Vec<i32>>=HashMap::new();
    let mut init_cluster_map: HashMap<u64, Vec<i32>> = HashMap::new();
    for entry in &sorted_entries{
        let id = entry.0;
        let sign_minis= &entry.1;

        //if we already have at least one cluster: compare the new read to the cluster(s)
        if clusters.len() > 0{

        }
        //we do not yet have a cluster and therefore need to fill the first read into the first
        else{
            for minimizer in sign_minis{
                //we get the minimizer sequence from the minimizer object
                let mini_seq= &minimizer.sequence;
                //calculate the hash of the minimizer sequence
                let mini_hash= calculate_hash(&mini_seq);
                init_cluster_map
                    .entry(mini_hash)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != id);
                // Check if id was retained (not a duplicate) and push it if needed
                let vec = init_cluster_map.get_mut(&mini_hash).unwrap();
                if !vec.contains(&id) {
                    vec.push(id);
                }
            }
        }
        //println!("{:?}",init_cluster_map)
        //println!("{}",sorted_entry.0);
        println!("{:?}",clusters);
    }
clusters
}


//clustering method for the case that we have an annotation to compare our reads against
pub(crate) fn cluster_from_initial(sorted_entries: Vec<(i32,Vec<Minimizer>)>, initial_clusters: HashMap<u64, Vec<i32>>) ->HashMap<i32, Vec<i32>>{
    //the hashmap containing the clusters we found TODO: verify:(key: cluster_id, value: vector of internal read ids
    let mut clusters: HashMap<i32, Vec<i32>> = HashMap::new();
    for sorted_entry in sorted_entries{
        let mut this_clusters= vec![];
        for mini in sorted_entry.1{
            let this_mini_hash= calculate_hash(&mini.sequence);
            if initial_clusters.contains_key(&this_mini_hash){
                let value= initial_clusters.get(&this_mini_hash).unwrap();
                //this_clusters.append(&mut value.clone());
                if value.len()>0{
                    this_clusters.push(value);
                }
                else{
                    //TODO: align the sequence against each reference and find best match
                    //Maybe add an offset as command line argument to isONclust, alternative: second iteration with shorter minimizers
                }

            }
        }
        println!("{}",sorted_entry.0);
        println!("{:?}",this_clusters);

        //clusters.add(this_clusters)
    }

    clusters
    }
pub(crate) fn add_backward_seqs(initial_clustering_records:Vec<FastaRecord>) ->Vec<FastaRecord>{
    let mut both_dir_records =vec![];
    let mut init_len=initial_clustering_records.len();
    for record in initial_clustering_records{
        init_len= init_len+1;
        let reversed: String = record.sequence.chars().rev().collect();
        let head = record.header.clone();

        both_dir_records.push(record.clone());
        let rev_record= FastaRecord{ sequence: reversed, header: init_len.to_string()};
        println!("{}, {}",head,rev_record);
        both_dir_records.push(rev_record)
    }

    both_dir_records
}

//if we have an annotation file we are interested to take the cluster number from the header of the respective read
fn get_id_from_header(head: String) -> i32 {
    // Attempt to parse the extracted number as i32
    let number_part: String = head.chars().filter(|c| c.is_digit(10)).collect();
    let parsed_number: Result<i32, _> = number_part.parse();

    match parsed_number {
        Ok(parsed) => {
            // Successfully parsed the number
            parsed
        }
        Err(_) => {
            println!("Failed to parse the extracted number as i32.");
            1 // Default value when parsing fails
        }
    }
}


//takes a string and hashes it via DefaultHasher. Used to improve search for minimizers in the data
fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}


pub(crate) fn get_initial_clustering(initial_clustering_records: Vec<FastaRecord>,k: usize,window_size: usize)->HashMap<u64, Vec<i32>>{
    let mut init_cluster_map: HashMap<u64, Vec<i32>> = HashMap::new();
    //iterate over the records of initial_clustering_records
    for record in initial_clustering_records{
        //get the header of record
        let head= record.header;
        //transform the header (expected e.g. SIRV1 to a number (i32), here 1
        let id= get_id_from_header(head);
        //println!("{}",id);
        //generate the minimizers for each record and store them in this_minimizers
        let this_minimizers= generate_sorted_fastq_new_version::get_kmer_minimizers(&*record.sequence, k, window_size);
        //let sub_minis= &this_minimizers[0..100];
        //println!("{:?}",sub_minis);
        //now iterate over all minimizers that we generated for this record
        for minimizer in this_minimizers{
            //we get the minimizer sequence from the minimizer object
            let mini_seq= minimizer.sequence;
            //calculate the hash of the minimizer sequence
            let mini_hash= calculate_hash(&mini_seq);
            //if the hash exists in the map: add the value (a new vector) or if the hash exists: extend the vector making up the value
            init_cluster_map
                    .entry(mini_hash)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != id);

            // Check if id was retained (not a duplicate) and push it if needed
            let vec = init_cluster_map.get_mut(&mini_hash).unwrap();
            //vec.push(id);
            if !vec.contains(&id) {
                vec.push(id);
            }

        }
    }
    //println!("{:?}",init_cluster_map);
    init_cluster_map
}
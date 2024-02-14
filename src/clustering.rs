use crate::structs:: FastaRecord;
use std::collections::HashMap;
use crate::structs::Minimizer_hashed;
use crate::{generate_sorted_fastq_new_version, structs};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use rustc_hash::FxHashMap;
use std::borrow::Cow;

pub(crate) fn reverse_complement(dna: &str) -> String {
    let reverse_complement: String = dna.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        }).clone()
        .collect();
    reverse_complement
}


fn calculate_shared_perc(nr_sign_minis:usize,value:i32)->f64{
    value as f64/nr_sign_minis as f64
}

fn detect_whether_shared(min_shared_minis:f64, shared_mini_infos: &FxHashMap<i32,i32>, minimizers: &Vec<Minimizer_hashed>) -> (bool, i32) {
    let mut most_shared= 0.0;
    let mut shared= false;
    let mut most_shared_cluster= 0;
    let mut nr_minis:usize;
    let mut shared_perc:f64;
    for (key, value) in shared_mini_infos {
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
        nr_minis = minimizers.len();
        shared_perc = calculate_shared_perc(nr_minis,*value);
        //println!("shared percentage between read and cluster {} : {}",key, shared_perc);
        if shared_perc > min_shared_minis && shared_perc > most_shared{
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared{
                shared = true;
            }
        }
    }
    (shared,most_shared_cluster)
}
//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster_de_novo(sign_minis: &Vec<Minimizer_hashed>,min_shared_minis:f64,minimizers: &Vec<Minimizer_hashed>,clusters:&mut FxHashMap<i32,Vec<i32>>,cluster_map: &mut FxHashMap<u64, Vec<i32>>, id: i32,  cl_id: &mut i32 ){
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of cluster ids as a value
    //we only need cl_id if cluster 0 already exists so we start with '1'
    //let mut cl_id= 1;
    let shared_perc_mini= min_shared_minis / 2.0_f64;
    //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
    let mut shared_mini_infos= FxHashMap::default();
    let mut shared_mini_infos_norm= FxHashMap::default();
    //TODO: This can be heavily improved if I add a field, high_confidence, to the seed object (here minimizer) We then can only pass over the minis with high_confidence=false
    //entry represents a read in our data
    //if we already have at least one cluster: compare the new read to the cluster(s)
    if !(clusters.is_empty()){
        //if sign_minis.len() > min_shared_minis as usize {
        for minimizer in sign_minis {
            //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
            if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                //iterate over belongs_to to update the counts of shared minimizers for each cluster
                for &belong_cluster in belongs_to {
                    //if the minimizer is already appointed to the cluster
                    if let Some(new_val) = shared_mini_infos.get(&belong_cluster) {
                        //increase the counter of minimizers shared with cluster
                        shared_mini_infos.insert(belong_cluster, *new_val + 1);
                    } else {
                        //add a new counter for the cluster element
                        shared_mini_infos.insert(belong_cluster, 1);
                    }
                }
            }
        }
        let mut most_shared_cluster= 0;
        let mut shared= false;
            //key: cluster_id, value: count of shared minimizers
            //TODO: sort wrt length of value and then only look at the first (best fit)
        (shared,most_shared_cluster) = detect_whether_shared(min_shared_minis, &shared_mini_infos, sign_minis);
        if !shared{
            for minimizer in minimizers{
                //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
                if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    for &belong_cluster in belongs_to {
                        //if the minimizer is already appointed to the cluster
                        if let Some(new_val) = shared_mini_infos_norm.get(&belong_cluster) {
                            //increase the counter of minimizers shared with cluster
                            shared_mini_infos_norm.insert(belong_cluster, *new_val + 1);
                        } else {
                            //add a new counter for the cluster element
                            shared_mini_infos_norm.insert(belong_cluster, 0);
                        }
                    }
                }
            }
            (shared,most_shared_cluster) = detect_whether_shared(shared_perc_mini, &shared_mini_infos_norm, &minimizers);
        }
        //if we have a cluster that we share enough minimizers with
        if shared {
            //add the read id to read_list
            let read_list = clusters.get_mut(&most_shared_cluster).unwrap();
            if !read_list.contains(&id) {
                read_list.push(id);
            }
            for sign_mini in sign_minis {
                cluster_map
                    .entry(sign_mini.sequence)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != most_shared_cluster);
                // Check if id was retained (not a duplicate) and push it if needed
                let vect = cluster_map.get_mut(&sign_mini.sequence).unwrap();
                if !vect.contains(&most_shared_cluster) {
                    vect.push(most_shared_cluster);
                }
            }
        }
        //we did not find a cluster that we could put the read into-> generate a new cluster
        else{
            for sign_mini in sign_minis {
                cluster_map
                    .entry(sign_mini.sequence)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != *cl_id);
                // Check if id was retained (not a duplicate) and push it if needed
                let vect = cluster_map.get_mut(&sign_mini.sequence).unwrap();
                if !vect.contains(&cl_id) {
                    vect.push(*cl_id);
                }
            }
            let id_vec=vec![id];
            clusters.insert(*cl_id,id_vec);
            *cl_id += 1;
        }
    }
    //we do not yet have a cluster and therefore need to fill the first read into the first
    else{
        //fill_first_cluster(&mut clusters, id, sign_minis, cluster_map);
        let init_id = 0;
        for minimizer in sign_minis {
            //fill cluster_map with the minimizers that we found in the first read
            cluster_map
                .entry(minimizer.sequence)
                .or_insert_with(Vec::new)
                .retain(|&existing_id| existing_id != init_id);
            // Check if id was retained (not a duplicate) and push it if needed
            let vect = cluster_map.get_mut(&minimizer.sequence).unwrap();
            if !vect.contains(&init_id) {
                vect.push(init_id);
            }
        }
        let id_vec=vec![id];
        clusters.insert(init_id,id_vec);
    }
}






pub(crate) fn add_rev_comp_seqs_annotation(initial_clustering_records:Vec<FastaRecord>) ->Vec<FastaRecord>{
    //Adds the reverse_compliment of the annotation
    //INPUT:   initial_clustering_records: The annotation
    //
    // OUTPUT:     both_dir_records: the records in both directions
    let mut both_dir_records =vec![];
    let mut init_len=initial_clustering_records.len();
    for record in initial_clustering_records{
        init_len += 1;
        //let reversed: String = record.sequence.chars().rev().collect();
        let reversed=reverse_complement(&record.sequence);
        let head = record.header.clone();
        both_dir_records.push(record.clone());
        let rev_record= FastaRecord{ sequence: reversed, header: init_len.to_string()};
        //println!("{}, {}",head,rev_record);
        both_dir_records.push(rev_record)
    }

    both_dir_records
}


pub(crate) fn add_rev_comp_seqs(initial_clustering_records:Vec<FastaRecord>) ->Vec<FastaRecord>{
    let mut both_dir_records =vec![];
    let mut init_len=initial_clustering_records.len();
    //iterate over the records in initial_clustering_records
    for record in initial_clustering_records{
        init_len += 1;
        //let reversed: String = record.sequence.chars().rev().collect();
        let reversed=reverse_complement(&record.sequence);
        let head = record.header.clone();
        //add the forward direction to both_dir_records
        both_dir_records.push(record.clone());
        //generate a new record
        let rev_record= FastaRecord{ sequence: reversed, header: init_len.to_string()};
        println!("{}, {}",head,rev_record);
        //add the reversed reads to both_dir_records
        both_dir_records.push(rev_record)
    }
    both_dir_records
}



//if we have an annotation file we are interested to take the cluster number from the header of the respective read
fn get_id_from_header(head: String) -> i32 {
    // Attempt to parse the extracted number as i32
    let number_part: String = head.chars().filter(|c| c.is_ascii_digit()).collect();
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


//takes an object T and hashes it via DefaultHasher. Used to improve search for minimizers in the data
pub fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
pub(crate) fn generate_initial_cluster_map(this_minimizers: &Vec<Minimizer_hashed>, init_cluster_map: &mut FxHashMap<u64, Vec<i32>>,identifier: i32){
    for minimizer in this_minimizers{
        init_cluster_map
            .entry(minimizer.sequence)
            .or_insert_with(Vec::new)
            .retain(|&existing_id| existing_id != identifier);
        // Check if id was retained (not a duplicate) and push it if needed
        let vec = init_cluster_map.get_mut(&minimizer.sequence).unwrap();
        //vec.push(id);
        if !vec.contains(&identifier) {
            vec.push(identifier);
        }
    }
}
//generates a mapping of minimizers to the reference genes (i.e. how many minimizers are shared with which cluster (stored in init_cluster_map))
/*pub(crate) fn generate_initial_clusters(initial_clustering_records: Vec<FastaRecord>,k: usize,window_size: usize,init_cluster_map:&mut FxHashMap<u64, Vec<i32>>){
    //holds the respective

    //let mut head;
    let mut id;
    let mut this_minimizers:Vec<Minimizer_hashed> =vec![];
    //let mut vec= vec![];
    //iterate over the records of initial_clustering_records
    for record in initial_clustering_records{
        //get the header of record
        head = record.header;
        //transform the header (expected e.g. SIRV1 to a number (i32), here 1
        id = get_id_from_header(head);
        //println!("{}",id);

        //generate the minimizers for each record and store them in this_minimizers
        generate_sorted_fastq_new_version::canonical_minimizers_hashed_old(&record.sequence, k, window_size, &mut this_minimizers);
        //let sub_minis= &this_minimizers[0..100];
        //println!("{:?}",sub_minis);
        //now iterate over all minimizers that we generated for this record
        for minimizer in &this_minimizers{
            //we get the minimizer sequence from the minimizer object
            //let mini_seq= minimizer.sequence;
            //calculate the hash of the minimizer sequence
            //if the hash exists in the map: add the value (a new vector) or if the hash exists: extend the vector making up the value
            init_cluster_map
                    .entry(minimizer.sequence)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != id);
            // Check if id was retained (not a duplicate) and push it if needed
            let vec = init_cluster_map.get_mut(&minimizer.sequence).unwrap();
            //vec.push(id);
            if !vec.contains(&id) {
                vec.push(id);
            }
        }
        this_minimizers.clear();
    }
    //println!("{:?}",init_cluster_map);
    init_cluster_map
}*/



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_reverse_complement() {
        let rev_comp = reverse_complement("GGGGATCATCAGGGCTA");
        assert_eq!(rev_comp,"TAGCCCTGATGATCCCC");
        let rev_comp2 = reverse_complement("ATCGA");
        assert_eq!(rev_comp2,"TCGAT");
    }

}
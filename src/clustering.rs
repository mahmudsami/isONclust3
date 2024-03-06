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
        .rev()//TODO: test whether into_par_iter works here
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
    for (key, value) in shared_mini_infos {//TODO: test whether into_par_iter works here
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
        nr_minis = minimizers.len();
        shared_perc = calculate_shared_perc(nr_minis,*value);
        //println!("shared percentage between read and cluster {} : {}",key, shared_perc);
        if shared_perc > min_shared_minis && shared_perc > most_shared && *value >= 3 {
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared {
                shared = true;
            }
        }
    }
    (shared,most_shared_cluster)
}
//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster(sign_minis: &Vec<Minimizer_hashed>,min_shared_minis:f64,minimizers: &Vec<Minimizer_hashed>,clusters:&mut FxHashMap<i32,Vec<i32>>,cluster_map: &mut FxHashMap<u64, Vec<i32>>, id: i32,  cl_id: &mut i32 ){
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
        for minimizer in sign_minis {//TODO: test whether into_par_iter works here
            //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
            if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                //iterate over belongs_to to update the counts of shared minimizers for each cluster
                for &belong_cluster in belongs_to {//TODO: test whether into_par_iter works here
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
        (shared,most_shared_cluster) = detect_whether_shared(min_shared_minis, &shared_mini_infos, sign_minis);
        if !shared{
            for minimizer in minimizers{//TODO: test whether into_par_iter works here
                //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
                if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    for &belong_cluster in belongs_to {//TODO: test whether into_par_iter works here
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
                cluster_map //TODO: test whether into_par_iter works here
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
                cluster_map //TODO: test whether into_par_iter works here
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
        for minimizer in sign_minis {//TODO: test whether into_par_iter works here
            //fill cluster_map with the minimizers that we found in the first read
            cluster_map //TODO: test whether into_par_iter works here
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


pub(crate) fn generate_initial_cluster_map(this_minimizers: &Vec<Minimizer_hashed>, init_cluster_map: &mut FxHashMap<u64, Vec<i32>>,identifier: i32){
    for minimizer in this_minimizers{//TODO: test whether into_par_iter works here
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
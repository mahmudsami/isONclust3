use crate::structs:: FastaRecord;
use std::collections::HashMap;
use crate::structs::Minimizer_hashed;
use crate::{generate_sorted_fastq_new_version, structs};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use rustc_hash::{FxHashMap, FxHashSet};
use std::borrow::{Cow, Borrow};
use bio::alignment::sparse::HashMapFx;


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

fn detect_whether_shared(min_shared_minis:f64, shared_seed_infos: &FxHashMap<i32,i32>, minimizers: &Vec<Minimizer_hashed>) -> (bool, i32) {
    let mut most_shared= 0.0;
    let mut shared= false;
    let mut most_shared_cluster= 0;
    let mut nr_minis: usize;
    let mut shared_perc: f64;
    for (key, nr_shared) in shared_seed_infos {//TODO: test whether into_par_iter works here
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
        nr_minis = minimizers.len();
        shared_perc = calculate_shared_perc(nr_minis, *nr_shared);
        //println!("shared percentage between read and cluster {} : {}",key, shared_perc);
        if shared_perc > min_shared_minis && shared_perc > most_shared{//} && *nr_shared >=0 {
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared {
                shared = true;
            }
        }
    }
    (shared,most_shared_cluster)
}

//TODO:
//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster(sign_minis: &Vec<Minimizer_hashed>, min_shared_minis:f64, minimizers: &Vec<Minimizer_hashed>, clusters:&mut FxHashMap<i32,Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, id: i32,  cl_id: &mut i32 ){
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of cluster ids as a value
    //we only need cl_id if cluster 0 already exists so we start with '1'
    //let mut cl_id= 1;
    let shared_perc_mini= min_shared_minis / 2.0_f64;
    //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
    let mut shared_seed_infos= FxHashMap::default();
    let mut shared_seed_infos_norm= FxHashMap::default();
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
                    if let Some(new_val) = shared_seed_infos.get(&belong_cluster) {
                        //increase the counter of minimizers shared with cluster
                        shared_seed_infos.insert(belong_cluster, *new_val + 1);
                    } else {
                        //add a new counter for the cluster element
                        shared_seed_infos.insert(belong_cluster, 1);
                    }
                }
            }
        }
        let mut most_shared_cluster= 0;
        let mut shared= false;
        //key: cluster_id, value: count of shared minimizers
        (shared,most_shared_cluster) = detect_whether_shared(min_shared_minis, &shared_seed_infos, sign_minis);
        if !shared{
            for minimizer in minimizers{//TODO: test whether into_par_iter works here
                //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
                if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    for &belong_cluster in belongs_to {//TODO: test whether into_par_iter works here
                        //if the minimizer is already appointed to the cluster
                        if let Some(new_val) = shared_seed_infos_norm.get(&belong_cluster) {
                            //increase the counter of minimizers shared with cluster
                            shared_seed_infos_norm.insert(belong_cluster, *new_val + 1);
                        } else {
                            //add a new counter for the cluster element
                            shared_seed_infos_norm.insert(belong_cluster, 0);
                        }
                    }
                }
            }
            (shared,most_shared_cluster) = detect_whether_shared(shared_perc_mini, &shared_seed_infos_norm, &minimizers);
        }
        //if we have a cluster that we share enough minimizers with
        if shared {
            //add the read id to read_list
            let read_list = clusters.get_mut(&most_shared_cluster).unwrap();
            if !read_list.contains(&id) {
                read_list.push(id);
            }
            //the following for-loop updates cluster_map
            for sign_mini in sign_minis {
                cluster_map //TODO: test whether into_par_iter works here
                    .entry(sign_mini.sequence)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != most_shared_cluster);
                // Check if id was retained (not a duplicate) and push it if needed
                //add the new cluster_id to the cluster_map should it not have been in there
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
//TODO: we could make shared_seed_infos available to already have an idea how many seeds a cluster has. Then we would need to update the ds. We then use cluster map and the update from the clustering and should have all we need
pub(crate) fn post_clustering_new(clusters: &mut FxHashMap<i32,Vec<i32>>,clusters_map: &mut FxHashMap<u64, Vec<i32>>){
    //cluster_seeds holds the cluster id and a vector of seed hashes (the ones contained in the cluster)
    let cluster_seeds: FxHashMap<i32,FxHashSet<u64>> = FxHashMap::default();
    //STEP1: get a list of minimizers for each cluster (stored in FXHashMap<cl_id,Vec<minimizer_hash>>)
    //iterate over the clusters_map to find out the overlaps between clusters (the two measures to calculate whether the clusters should be merged)
    //TODO: sort in ascending order so that we start with the lowest number of ids to the highest (Sort by vec_of_ids.len())
    //TODO: use hashset instead of vector to store the seed_hashes->easy intersection!
    for (seed_hash, vec_of_ids) in clusters_map.into_iter() {
        //iterate over the ids that we have stored in the value of each minimizer
        for (i, &id) in vec_of_ids.iter().enumerate() {
            if let Some(vec_of_ids) = cluster_seeds.get(&id) {
                if !vec_of_ids.contains(seed_hash){
                    vec_of_ids.push(*seed_hash);
                }

            }
            else{
                let mini_vec=vec![*mini];
                cluster_seeds.insert(*cl_id,mini_vec);

            }
        }
    }
}

pub(crate) fn post_clustering(clusters:&mut FxHashMap<i32,Vec<i32>>,clusters_map: &mut FxHashMap<u64, Vec<i32>>){
    let min_shared_perc=0.8;
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of cluster ids as a value

    //cluster_seeds_hash_map contains the cl_id as key and the count of seeds as a value
    let mut cluster_seeds_hash_map:FxHashMap<i32,Vec<u64>>= FxHashMap::default();
    //cluster_seed_shared_map contains the two cluster_ids (smaller, higher) as a key and the count of common seeds as a value
    let mut cluster_seeds_shared_map: FxHashMap<String,i32> = FxHashMap::default();
    let mut small_id= 0;
    let mut high_id= 0;
    //iterate over the clusters_map to find out the overlaps between clusters and how many seeds each cluster has (the two measures to calculate whether the clusters should be merged)
    //TODO: sort in ascending order so that we start with the lowest number of ids to the highest (Sort by vec_of_ids.len())
    for (mini, vec_of_ids) in clusters_map.into_iter(){
        //iterate over the ids that we have stored in the value of each minimizer
        for (i, &id) in vec_of_ids.iter().enumerate() {
        //for id in vec_of_ids{
            //the cluster is already in the cluster_seeds_hash_map ->increase count by one
            cluster_seeds_hash_map.entry(id).or_insert_with(Vec::new).push(*mini);
            //we also need to get the counts of overlaps( this is what we do here)
            for &id2 in vec_of_ids.iter().skip(i + 1) {
            //for id2 in vec_of_ids[id..]{
                if id < id2{
                    small_id = id;
                    high_id = id2;
                }
                else if id2 < id{
                    small_id = id2;
                    high_id = id;
                }
                let full_id = small_id.to_string()+","+ &*high_id.to_string();
                if let Some(shared_count) = cluster_seeds_shared_map.get(full_id.as_str()) {
                    //increase the counter of minimizers shared with cluster
                    cluster_seeds_shared_map.insert(full_id, *shared_count + 1);
                }
                else {
                    //add a new counter for the cluster element
                    cluster_seeds_shared_map.insert(full_id, 1);
                }
            }
        }
    }
    //clustered_bool_map holds the cl_id and a bool to indicate whether the cluster has been merged into another cluster
    let  mut clustered_bool_map:FxHashMap<i32,bool> = FxHashMap::default();
    //merged_into is a hashmap containing the cl_id of the subcluster as a key and the overcluster as value
    let  mut merged_into:FxHashMap<i32,i32> = FxHashMap::default();
    for (ids,shared_nr) in cluster_seeds_shared_map{
        //retreive the ids from the appended ids_string
        let mut ids_s: Vec<_>= ids.split(",").collect();
        let first_id:i32 = ids_s[0].parse::<i32>().unwrap();
        let second_id:i32 = ids_s[1].parse::<i32>().unwrap();
        //retreive the counts of the overall seeds for the clusters
        let first_count: usize = cluster_seeds_hash_map.get(&first_id).unwrap().len();
        let snd_count: usize = cluster_seeds_hash_map.get(&second_id).unwrap().len();
        //calculate the rate of shared clusters from each total number of seeds
        let first_perc = shared_nr as f64/first_count as f64;
        let snd_perc = shared_nr as f64/snd_count as f64;
        let mut small_cl_id= 0;
        let mut large_cl_id = 0;
        let mut merge_cl_small:Vec<u64> = vec![];
        //the first cluster has more shared minis with the second cluster (is smaller)
        if first_perc >= snd_perc{
            small_cl_id = second_id;
            large_cl_id = first_id;
            merged_into.insert(first_id,second_id);
            merge_cl_small = cluster_seeds_hash_map.get(&second_id).unwrap().to_vec();
            clustered_bool_map.insert(small_cl_id,true)
        }
        //the second cluster has more shared minis with the first cluster (is smaller)
        else{
            small_cl_id = first_id;
            large_cl_id = second_id;
            merged_into.insert(second_id,first_id);
            merge_cl_small = cluster_seeds_hash_map.get(&first_id).unwrap().to_vec();
        }
        //TODO: add this properly into the code ( adds the significant seeds that were not in the larger cluster but in the smaller one
        //here we want to add the significant seeds from the small cluster also into the large cluster
        for add_seed in merge_cl_small {
            cluster_seeds_hash_map //TODO: test whether into_par_iter works here
                .entry(add_seed.sequence)
                .or_insert_with(Vec::new)
                .retain(|&existing_id| existing_id != most_shared_cluster);
        }
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
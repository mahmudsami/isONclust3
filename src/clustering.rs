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
    let mut most_shared_cluster= -1;
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

//shared_seed_infos: hashmap that holds read_id->nr shared minimizers with clusters->not updated when cluster changes!
//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster(sign_minis: &Vec<Minimizer_hashed>, min_shared_minis:f64, minimizers: &Vec<Minimizer_hashed>, clusters:&mut FxHashMap<i32,Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, id: i32,  cl_id: &mut i32,shared_seed_infos:&mut FxHashMap<i32,i32> ){
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of ids as a value
    //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
    //let mut shared_seed_infos= FxHashMap::default();
    let shared_perc_mini = min_shared_minis / 2.0_f64;
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
            if !sign_minis.is_empty(){
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
            else{
                let id_vec=vec![id];
                clusters.insert(*cl_id,id_vec);
                *cl_id += 1;
            }
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

//takes clusters_map as input and generates cl_set_map: a Hashmap containing the cluster id as key and a hashset of seedhashes as value.
fn generate_post_clustering_ds(cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>,clusters_map: &mut FxHashMap<u64, Vec<i32>>){
    //TODO: count overlaps between clusters to sort and identify all cluster_combinations
    for (mini, vec_of_ids) in clusters_map.into_iter() {
        //iterate over the ids that we have stored in the value of each minimizer
        for id in vec_of_ids.iter() {
            //the cluster is already in the cluster_seeds_hash_map ->add the seed hash to the hashset, otherwise add new hashset with the seed hash
            if cl_set_map.contains_key(id){
                cl_set_map.get_mut(id).unwrap().insert(*mini);
            }
            else{
                let mut this_set:FxHashSet<u64>=FxHashSet::default();
                this_set.insert(*mini);
                cl_set_map.insert(*id,this_set);
            }
        }
    }
}


//helper function for the post_clustering step: Updates the 'clusters' and 'clusters_map' data structures
fn update_clusters(clusters: &mut FxHashMap<i32,Vec<i32>>, clusters_map: &mut FxHashMap<u64,Vec<i32>>, smallHS:&FxHashSet<u64>, large_cluster_id:&i32, small_cluster_id:&i32){
    //get the infos of clusters that belong to the two clusters we want to merge
    let binding = clusters.clone();
    let small_cl_info= binding.get(small_cluster_id).unwrap();
    let large_cl_info= clusters.get_mut(large_cluster_id).unwrap();
    //add the reads of the small cluster into the large cluster
    large_cl_info.extend(small_cl_info);
    //also add the hashes of the small cluster into the large cluster
    for seed_hash in smallHS{
        let mut cl_vec=clusters_map.get_mut(&seed_hash).unwrap();
        cl_vec.push(*large_cluster_id);
    }
}


//merges the clusters for cl_set_map and calls update_clusters
fn merge_clusters(clusters: &mut FxHashMap<i32,Vec<i32>>, clusters_map: &mut FxHashMap<u64,Vec<i32>>, cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, large_cluster_id:&i32, small_cluster_id:&i32){
    let binding=cl_set_map.clone();
    let mut largeHs: &mut FxHashSet<u64> = cl_set_map.get_mut(large_cluster_id).unwrap();
    let mut smallHs: &FxHashSet<u64> = binding.get(small_cluster_id).unwrap();
    //TODO: fix type here
    let mut intersect: FxHashSet<u64> = largeHs.intersection(&smallHs).cloned().collect();
    largeHs = &mut intersect;
    //largeHs = largeHs.intersection(&smallHs).cloned().collect();
    update_clusters(clusters,clusters_map,smallHs,large_cluster_id,small_cluster_id);
}


pub(crate) fn post_clustering_new(clusters: &mut FxHashMap<i32,Vec<i32>>,clusters_map: &mut FxHashMap<u64, Vec<i32>>, min_shared_minis:f64){
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    let mut cl_set_map: FxHashMap<i32,FxHashSet<u64>> = FxHashMap::default();
    //TODO update mergedHS and updatedHS and use them to our advantage
    let mut mergedHS: FxHashSet<i32>=FxHashSet::default();
    let mut updatedHS: FxHashSet<i32>=FxHashSet::default();

    generate_post_clustering_ds(&mut cl_set_map, clusters_map);
    //clusters holds the cluster_id as key and a vector of read_ids as value
    //cluster_map holds the hash of a seed as key and a vector of cluster_ids as value
    //shared_seed_infos holds
    //cluster_seeds holds the cluster id and a hashset of seed hashes (the ones contained in the cluster)
    //STEP1: get a list of minimizers for each cluster (stored in FXHashMap<cl_id,Vec<minimizer_hash>>)
    //iterate over the clusters_map to find out the overlaps between clusters (the two measures to calculate whether the clusters should be merged)
    //TODO: sort in ascending order so that we start with the lowest number of ids to the highest (Sort by vec_of_ids.len())
    //TODO: use hashset instead of vector to store the seed_hashes->easy intersection!
    let cl_binding= cl_set_map.clone();
    // Convert the HashMap into a Vec of references to its entries
    let keys: Vec<&i32> = cl_binding.keys().collect();
    // Iterate over each pair of entries without repeating comparisons
    for i in 0..keys.len() {
        //get key and value of the first entry
        let key1 = keys[i];
        let value1 = cl_binding.get(&key1).unwrap();
        //iterate over second entry but do not repeat already compared entries
        for j in i + 1..keys.len() {
            //get key and value of the second entry
            let key2 = keys[j];
            if mergedHS.contains(key2){
                //retrieve the value to key2 (the seed_hashes for the cluster)
                let value2 = cl_binding.get(&key2).unwrap();
                //get the number of seed_hashes of the intersection of the two clusters
                let intersect_len= value1.intersection(&value2).collect::<Vec<_>>().len();
                //compute the rates of how much of the seed_hashes of the clusters are shared
                let shared1 = intersect_len as f64/ value1.len() as f64;
                let shared2 = intersect_len as f64/ value2.len() as f64;
                //cases covered: shared1 but not shared 2, shared 2 but not shared 1, shared1 and shared2,
                if shared1 > min_shared_minis{
                    if shared2> min_shared_minis{
                        //shared1 and shared2
                        if shared1 >= shared2{
                            //TODO: merge shared1 into shared2
                            merge_clusters(clusters, clusters_map, &mut cl_set_map, key2, &key1);
                            mergedHS.insert(*key1);
                            updatedHS.insert(*key2);
                        }
                        else{
                            //TODO: merge shared2 into shared1
                            merge_clusters(clusters, clusters_map, &mut cl_set_map, &key1, key2);
                            mergedHS.insert(*key2);
                            updatedHS.insert(*key1);
                        }
                    }
                    else{
                        //shared1 but not shared2
                        //TODO: merge shared 1 into shared2
                        merge_clusters(clusters, clusters_map, &mut cl_set_map, key2, &key1);
                        mergedHS.insert(*key1);
                        updatedHS.insert(*key2);
                    }
                }
                else if shared2 > min_shared_minis{
                    //shared2 but not shared1
                    //TODO: merge shared2 into shared1
                    merge_clusters(clusters, clusters_map, &mut cl_set_map, &key1, key2);
                    mergedHS.insert(*key2);
                    updatedHS.insert(*key1);
                }
                else{
                    //Do not merge anything
                }
            }

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
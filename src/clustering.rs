use crate::structs::Minimizer_hashed;

use rustc_hash::{FxHashMap, FxHashSet};
//use rayon::prelude::*;
use crate::{Cluster_ID_Map, Seed_Map};

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


fn calculate_shared_perc(nr_sign_minis: usize,value: i32)->f64{
    value as f64 / nr_sign_minis as f64
}


fn detect_whether_shared(min_shared_minis:f64, shared_seed_infos: &FxHashMap<i32,i32>, minimizers: &Vec<Minimizer_hashed>) -> (bool, i32) {
    let mut most_shared= 0.0;
    let mut shared= false;
    let mut most_shared_cluster= -1;
    let mut nr_minis=minimizers.len();
    let mut shared_perc: f64;
    for (key, nr_shared) in shared_seed_infos {//TODO: test whether into_par_iter works here
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
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
pub(crate) fn cluster(sign_minis: &Vec<Minimizer_hashed>, min_shared_minis:f64, minimizers: &Vec<Minimizer_hashed>, clusters:&mut Cluster_ID_Map, cluster_map: &mut Seed_Map, id: i32,  cl_id: &mut i32 ){
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of ids as a value
    //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
    //let mut shared_seed_infos= FxHashMap::default();
    let shared_perc_mini = min_shared_minis / 2.0_f64;
    //let mut shared_seed_infos_norm= FxHashMap::default();
    let mut shared_seed_infos_norm_vec: Vec<i32> = vec![0; clusters.len()];

    //TODO: This can be heavily improved if I add a field, high_confidence, to the seed object (here minimizer) We then can only pass over the minis with high_confidence=false
    //entry represents a read in our data
    //if we already have at least one cluster: compare the new read to the cluster(s)
    if !(clusters.is_empty()){
        //if sign_minis.len() > min_shared_minis as usize {
/*        for minimizer in sign_minis {//TODO: test whether into_par_iter works here
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
        }*/
        let mut most_shared_cluster= -1;
        let mut shared= false;
        //key: cluster_id, value: count of shared minimizers
        // (shared,most_shared_cluster) = detect_whether_shared(min_shared_minis, &shared_seed_infos, sign_minis);
        //println!("{} shared when using filtered", shared);
        if !shared{
            for minimizer in minimizers{//TODO: test whether into_par_iter works here
                //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
                if let Some(belongs_to) = cluster_map.get(&minimizer.sequence) {
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    for &belong_cluster in belongs_to {//TODO: test whether into_par_iter works here
                        //iterate over belongs_to to update the counts of shared minimizers for each cluster
                        // println!("belong_cluster {} belong_cluster as usize {}", belong_cluster, belong_cluster as usize );
                        shared_seed_infos_norm_vec[belong_cluster as usize] += 1;
                        //if the minimizer is already appointed to the cluster
                        // if let Some(new_val) = shared_seed_infos_norm.get(&belong_cluster) {
                        //     //increase the counter of minimizers shared with cluster
                        //     shared_seed_infos_norm.insert(belong_cluster, *new_val + 1);
                        // } else {
                        //     //add a new counter for the cluster element
                        //     shared_seed_infos_norm.insert(belong_cluster, 1);
                        // }
                    }
                }
            }

            //println!("shared_seed_infos_norm_vec {:?}", shared_seed_infos_norm_vec );

            if let Some((max_cluster_id, max_shared)) = shared_seed_infos_norm_vec.iter().enumerate().max_by_key(|&(_, value)| value) {
                let nr_minis= minimizers.len();
                let mut shared_perc: f64;
                //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
                shared_perc = calculate_shared_perc(nr_minis, *max_shared);
                //println!("max_cluster_id: {}, max_shared: {}, shared_perc: {}, nr_minis: {}",max_cluster_id,max_shared,shared_perc,nr_minis );
                // println!("shared percentage between read and cluster : {} min required: {}", shared_perc, shared_perc_mini);
                if shared_perc > shared_perc_mini {
                    shared = true;
                    most_shared_cluster= max_cluster_id as i32;
                }
            }

            //(shared,most_shared_cluster) = detect_whether_shared(shared_perc_mini, &shared_seed_infos_norm, &minimizers);
            //println!("{} shared when using ALL", shared);

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
                    if !vect.contains(cl_id) {
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
fn generate_post_clustering_ds(cl_set_map: &mut FxHashMap<i32, Vec<u64>>, clusters_map: &mut Seed_Map){
    //TODO: count overlaps between clusters to sort and identify all cluster_combinations
    //let mut clusters_that_overlap: FxHashMap<i32,FxHashSet<i32>> = FxHashMap::default();
    //println!("Clusters_map len {}",clusters_map.len());
    for (mini, vec_of_ids) in clusters_map {
        //iterate over the ids that we have stored in the value of each minimizer
        for i in 0..vec_of_ids.len() {
            let id = vec_of_ids[i];
            //the cluster is already in the cluster_seeds_hash_map ->add the seed hash to the hashset, otherwise add new hashset with the seed hash
            if cl_set_map.contains_key(&id){
                cl_set_map.get_mut(&id).unwrap().push(*mini);
            }
            else{
                let mut this_set: Vec<u64> = vec![*mini];
                //this_set.push(*mini);
                cl_set_map.insert(id,this_set);
            }
            for j in i + 1..vec_of_ids.len() {
                let id2 = vec_of_ids[j];
                let mut small_id;
                let mut large_id;
                if id < id2{
                    small_id = id;
                    large_id = id2;
                }
                else{
                    small_id = id2;
                    large_id = id;
                }
            }
        }
    }
}



//helper function for the post_clustering step: Updates the 'clusters' and 'clusters_map' data structures
fn update_clusters(clusters: &mut Cluster_ID_Map, clusters_map: &mut Seed_Map, small_hs: &Vec<u64>, large_cluster_id: &i32, small_cluster_id:&i32){

    //get the infos of clusters that belong to the two clusters we want to merge
    let binding = clusters.clone();
    let small_cl_info= binding.get(small_cluster_id).unwrap();
    let large_cl_info= clusters.get_mut(large_cluster_id).unwrap();
    //add the reads of the small cluster into the large cluster
    large_cl_info.extend(small_cl_info);
    clusters.remove_entry(small_cluster_id);
    //also add the hashes of the small cluster into the large cluster
    for seed_hash in small_hs{
        let mut cl_vec= clusters_map.get_mut(seed_hash).unwrap();
        //TODO:if higher cluster id is not in the vector just replace the small cl id by the large cl id
        if !cl_vec.contains(large_cluster_id){
            cl_vec.push(*large_cluster_id);
        }
        //delete small_cluster_id from the seed hashes so we do not have any indication of the cluster in the ds
        cl_vec.retain(|x| *x != *small_cluster_id);
    }
}


//merges the clusters for cl_set_map and calls update_clusters
fn merge_clusters(clusters: &mut Cluster_ID_Map, clusters_map: &mut Seed_Map, cl_set_map: & FxHashMap<i32,Vec<u64>>, large_cluster_id:&i32, small_cluster_id:&i32){
    println!("Merging cl {} into cluster {}",small_cluster_id,large_cluster_id);
    //let binding= cl_set_map.clone();
    let mut small_hs: &Vec<u64> = cl_set_map.get(small_cluster_id).unwrap();
    update_clusters(clusters, clusters_map, small_hs, large_cluster_id, small_cluster_id);
}


fn detect_overlaps(nr_clusters: usize, cl_set_map: & FxHashMap<i32,Vec<u64>>, cluster_map: &mut Seed_Map, merge_into: &mut Vec<(i32,i32)>, min_shared_minis: f64, small_hs: &mut FxHashSet<i32>){
    //shared_seed_infos_vec: a vector
    let mut shared_seed_infos_vec: Vec<i32> = vec![0; nr_clusters];
    //println!("ALL KEYS: {:?}", cl_set_map.keys());

    for (cl_id,hashes) in cl_set_map.iter(){
        //println!("cl_id: {} nr minimizers: {}",cl_id, hashes.len());
        //iterate over the hashes for each cl_id
        for hash in hashes.iter() {
            if let Some(belongs_to) = cluster_map.get(hash) {
                //iterate over belongs_to to update the counts of shared minimizers for each cluster
                for &belong_cluster in belongs_to {
                    if belong_cluster != *cl_id {
                        shared_seed_infos_vec[belong_cluster as usize] += 1;
                    }
                }
            }
        }
        //println!("Post cluster vec: {:?}", shared_seed_infos_vec);

        if let Some((max_cluster_id, max_shared)) = shared_seed_infos_vec.iter().enumerate().max_by_key(|&(_, value)| value) {
            let nr_minis= hashes.len();
            let mut shared_perc: f64;
            let most_shared_cluster_id= max_cluster_id as i32;
            //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
            shared_perc = calculate_shared_perc(nr_minis, *max_shared);
            //println!("shared percentage between cluster: {} and max_cluster_id {} : {}. Min required: {}. Max shared {}",cl_id, max_cluster_id, shared_perc, min_shared_minis, max_shared);
            if shared_perc > min_shared_minis {
                //println!("ENTERING MERGE");
                //println!("Nr elems_this {}, other {}", nr_elems_this,nr_elems_other);
                //we have a new best cluster as soon as
                if  nr_minis < cl_set_map.get(&most_shared_cluster_id).unwrap().len() {
                    merge_into.push((*cl_id, most_shared_cluster_id));
                    small_hs.insert(*cl_id);
                }
                else{
                    merge_into.push(( most_shared_cluster_id,*cl_id));
                    small_hs.insert(most_shared_cluster_id);
                }
            }
        }
        // clear count vector for next cluster
        for item in &mut shared_seed_infos_vec { *item = 0; }
        //println!("Post cluster vec CLEARING: {:?}", shared_seed_infos_vec);
    }
}


/*fn detect_whether_shared_other(this_id: i32,min_shared_minis:f64, shared_seed_infos: &FxHashMap<i32,i32>, minimizers:  Vec<u64>) -> (bool, i32) {
    let mut most_shared= 0.0;
    let mut shared= false;
    let mut most_shared_cluster= -1;
    let mut nr_minis= minimizers.len();
    let mut shared_perc: f64;
    for (key, nr_shared) in shared_seed_infos {//TODO: test whether into_par_iter works here
        //we have more shared minis with the cluster than our threshold and this is the cluster we share the most minimizers with
        shared_perc = calculate_shared_perc(nr_minis, *nr_shared);
        //println!("shared percentage between read and cluster {} : {}",key, shared_perc);
        if shared_perc > min_shared_minis && shared_perc > most_shared && this_id != *key {
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared {
                shared = true;
            }
        }
    }
    (shared,most_shared_cluster)
}*/


fn merge_clusters_from_merge_into(merge_into: &mut Vec<(i32,i32)>, clusters_map: &mut  Seed_Map, clusters: &mut Cluster_ID_Map, cl_set_map: &mut FxHashMap<i32,Vec<u64>>, small_hs: &FxHashSet<i32>){
    //println!("Merge_into_len: {}",merge_into.len());
    let mut not_mergeable_cter= 0;
    for (id , value) in merge_into{
        let mut large_id = value;
        //we might already have deleted large_id from clusters during this iteration
        if clusters.contains_key(large_id) {
            //idea here: we merge the ids into larger clusters first, smaller clusters are still bound to merge into the new cluster later
            if !small_hs.contains(large_id){
                merge_clusters( clusters, clusters_map, cl_set_map,large_id,id)
            }
            else{
                not_mergeable_cter += 1;
            }
        }
    }
}


pub(crate) fn post_clustering(clusters: &mut Cluster_ID_Map, cluster_map: &mut Seed_Map, min_shared_minis:f64){
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    let mut cl_set_map: FxHashMap<i32,Vec<u64>> = FxHashMap::default();
    //merge_into is a vector of a tuple(cl_id1,cl_id2)
    let mut merge_into: Vec<(i32,i32)> = vec![];
    //small_hs is a HashSet storing all cluster ids that were merged into other clusters during this iteration
    let mut small_hs: FxHashSet<i32> = FxHashSet::default();
    //used to have do-while structure
    let mut first_iter= true;
    let nr_clusters= clusters.len();
    //continue merging as long as we still find clusters that we may merge
    while !merge_into.is_empty() || first_iter {
        //let nr_cl_cter = cl_set_map.len();
        //let mut hash_cter = 0;
        //println!("Post Cluster iter");
        println!("Merge into: {}", merge_into.len());
        println!("Nr clusters: {}", clusters.len());
        //clear merge_into as this is the indicator how often we attempt to merge further (the while loop depends on it)
        merge_into.clear();
        small_hs.clear();
        //set first_iter to be false to not stay in a infinity loop
        first_iter = false;
        //merge_into contains the information about which clusters to merge into which
        //generate the data structure giving us merge infos
        generate_post_clustering_ds(&mut cl_set_map,  cluster_map);
        //for (cl_id, hash_set) in cl_set_map.iter(){
        //    hash_cter += hash_set.len();
        //}
        //println!("Nr clusters: {}, nr hashes {}",nr_cl_cter, hash_cter);
        detect_overlaps(nr_clusters, &cl_set_map, cluster_map, &mut merge_into, min_shared_minis, &mut small_hs);
        //println!("# Elements in small_hs {}",small_hs.len());
        //merges the clusters
        merge_clusters_from_merge_into(&mut merge_into, cluster_map, clusters, &mut cl_set_map, &small_hs);
        cl_set_map.clear();
    }
}


pub(crate) fn generate_initial_cluster_map(this_minimizers: &Vec<Minimizer_hashed>, init_cluster_map: &mut Seed_Map,identifier: i32){
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
use crate::structs:: FastaRecord;
use std::collections::{HashMap, HashSet};
use crate::structs::Minimizer_hashed;
use crate::{generate_sorted_fastq_new_version, structs};

use rustc_hash::{FxHashMap, FxHashSet};
use std::borrow::{Cow, Borrow};
use bio::alignment::sparse::HashMapFx;
use rayon::prelude::*;

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
pub(crate) fn cluster(sign_minis: &Vec<Minimizer_hashed>, min_shared_minis:f64, minimizers: &Vec<Minimizer_hashed>, clusters:&mut FxHashMap<i32,Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, id: i32,  cl_id: &mut i32 ){
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of ids as a value
    //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
    let mut shared_seed_infos= FxHashMap::default();
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
                            shared_seed_infos_norm.insert(belong_cluster, 1);
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
    /*else{
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
            let id_vec=vec![id];
            clusters.insert(*cl_id,id_vec);
            *cl_id += 1;
        }
    }*/
}

//takes clusters_map as input and generates cl_set_map: a Hashmap containing the cluster id as key and a hashset of seedhashes as value.
fn generate_post_clustering_ds_other(cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, clusters_map: &mut FxHashMap<u64, Vec<i32>>){
    //TODO: count overlaps between clusters to sort and identify all cluster_combinations
    let mut clusters_that_overlap: FxHashMap<i32,FxHashSet<i32>>=FxHashMap::default();
    //println!("Clusters_map len {}",clusters_map.len());
    for (mini, vec_of_ids) in clusters_map.into_iter() {
        //iterate over the ids that we have stored in the value of each minimizer
        for i in 0..vec_of_ids.len() {
            let id = vec_of_ids[i];
            //the cluster is already in the cluster_seeds_hash_map ->add the seed hash to the hashset, otherwise add new hashset with the seed hash
            if cl_set_map.contains_key(&id){
                cl_set_map.get_mut(&id).unwrap().insert(*mini);
            }
            else{
                let mut this_set: FxHashSet<u64> = FxHashSet::default();
                this_set.insert(*mini);
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
                if clusters_that_overlap.contains_key(&small_id){
                    let other_ids= clusters_that_overlap.get_mut(&small_id).unwrap();
                    other_ids.insert(large_id);
                }
                else{
                    let mut other_ids=FxHashSet::default();
                    other_ids.insert(large_id);
                    clusters_that_overlap.insert(small_id,other_ids);
                }
            }
        }
    }
}


//takes clusters_map as input and generates cl_set_map: a Hashmap containing the cluster id as key and a hashset of seedhashes as value.
fn generate_post_clustering_ds(cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, cl_overlaps: &mut FxHashMap<i32,FxHashMap<i32,i32>>, clusters_map: &mut FxHashMap<u64, Vec<i32>>){
    //println!("Clusters_map len {}",clusters_map.len());
    for (mini, vec_of_ids) in clusters_map.into_iter() {
        //iterate over the ids that we have stored in the value of each minimizer
        for i in 0..vec_of_ids.len() {
            let id = vec_of_ids[i];
            //cl_set_map contains the cl_id as key and the set of seed hashes as value
            //the cluster is already in the cluster_seeds_hash_map ->add the seed hash to the hashset, otherwise add new hashset with the seed hash
            if cl_set_map.contains_key(&id){
                cl_set_map.get_mut(&id).unwrap().insert(*mini);
            }
            else{
                let mut this_set: FxHashSet<u64> = FxHashSet::default();
                this_set.insert(*mini);
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
                if cl_overlaps.contains_key(&small_id){
                    let id_overlaps= cl_overlaps.get_mut(&small_id).unwrap();
                    if id_overlaps.contains_key(&large_id){
                        let mut overlap = id_overlaps.get_mut(&large_id).unwrap();
                        *overlap += 1;
                    }
                    else {
                        id_overlaps.insert(large_id,1);
                    }
                }
                else{
                    let mut id_overlaps= FxHashMap::default();
                    id_overlaps.insert(large_id,1);
                    cl_overlaps.insert(small_id,id_overlaps);
                }
            }
        }
    }
}


//helper function for the post_clustering step: Updates the 'clusters' and 'clusters_map' data structures
fn update_clusters(clusters: &mut FxHashMap<i32,Vec<i32>>, clusters_map: &mut FxHashMap<u64,Vec<i32>>, small_hs:&FxHashSet<u64>, large_cluster_id:&i32, small_cluster_id:&i32){

    //get the infos of clusters that belong to the two clusters we want to merge
    let binding = clusters.clone();
    let small_cl_info= binding.get(small_cluster_id).unwrap();
    let large_cl_info= clusters.get_mut(large_cluster_id).unwrap();
    //add the reads of the small cluster into the large cluster
    large_cl_info.extend(small_cl_info);
    clusters.remove_entry(small_cluster_id);
    //also add the hashes of the small cluster into the large cluster
    for seed_hash in small_hs{
        let mut cl_vec= clusters_map.get_mut(&seed_hash).unwrap();
        //TODO:if higher cluster id is not in the vector just replace the small cl id by the large cl id
        if !cl_vec.contains(large_cluster_id){
            cl_vec.push(*large_cluster_id);
        }
        //TODO:delete the
        //delete small_cluster_id from the seed hashes so we do not have any indication of the cluster in the ds
        cl_vec.retain(|x| *x != *small_cluster_id);
    }
}


//merges the clusters for cl_set_map and calls update_clusters
fn merge_clusters(clusters: &mut FxHashMap<i32,Vec<i32>>, clusters_map: &mut FxHashMap<u64,Vec<i32>>, cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, large_cluster_id:&i32, small_cluster_id:&i32){
    println!("Merging cl {} into cluster {}",small_cluster_id,large_cluster_id);
    let binding= cl_set_map.clone();
    //let mut large_hs: &mut FxHashSet<u64> = cl_set_map.get_mut(large_cluster_id).unwrap();
    let mut small_hs: &FxHashSet<u64> = binding.get(small_cluster_id).unwrap();
    //let mut intersect: FxHashSet<u64> = large_hs.intersection(&small_hs).cloned().collect();
    //large_hs = &mut intersect;
    //largeHs = largeHs.intersection(&smallHs).cloned().collect();
    update_clusters(clusters, clusters_map, small_hs, large_cluster_id, small_cluster_id);
}



fn fill_merge_into(cl_overlaps: &mut FxHashMap<i32,FxHashMap<i32,i32>>,merge_into: &mut FxHashMap<i32,(i32,f64)>, min_shared_minis: f64, cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>,small_hs: &mut FxHashSet<i32>){
    //iterate over cl_overlaps to retreive the first id and
    for (first_id,overlap_hm) in cl_overlaps {
        for (second_id, overlap) in overlap_hm {
            let nr_seeds1 = cl_set_map.get(&first_id).unwrap().len();
            let nr_seeds2 = cl_set_map.get(&second_id).unwrap().len();
            //compute the rates of how much of the seed_hashes of the clusters are shared
            let shared1 = *overlap as f64/ nr_seeds1 as f64;
            let shared2 = *overlap as f64/ nr_seeds2 as f64;
            //println!("Shareds {},{}: {}, {}",first_id,second_id,shared1,shared2);
            //cases covered: shared1 but not shared 2, shared 2 but not shared 1, shared1 and shared2,
            //shared1 and shared2
            if shared1 > min_shared_minis{
                if shared2 > min_shared_minis{
                    //shared1>shared2: we want to merge cluster1 into cluster2
                    if shared1 >= shared2 {
                        if merge_into.contains_key(&first_id) {
                            let mut value = merge_into.get_mut(&first_id).unwrap();
                            let mut other_id= value.0;
                            let mut shared = value.1;
                            if shared1 > shared {
                                other_id = *second_id;
                                shared = shared1;
                            }
                        }
                        else {
                            small_hs.insert(*first_id);
                            merge_into.insert(*first_id, (*second_id, shared1));
                        }
                    }
                    else{
                        //TODO: merge shared1 into shared2
                        if merge_into.contains_key(&second_id) {
                            let mut value  = merge_into.get_mut(&second_id).unwrap();
                            let mut other_id= value.0;
                            let mut shared = value.1;
                            if shared2 > shared {
                                other_id = *first_id;
                                shared = shared2;
                            }
                        }
                        else{
                            small_hs.insert(*second_id);
                            merge_into.insert(*second_id,(*first_id,shared2));
                        }
                    }
                }
                else{
                    //shared1 but not shared2
                    //TODO: merge shared 1 into shared 2
                    //merge_clusters(clusters, clusters_map, &mut cl_set_map, key2, &key1);
                    if merge_into.contains_key(&first_id) {
                        let mut value = merge_into.get_mut(&first_id).unwrap();
                        let mut other_id= value.0;
                        let mut shared = value.1;
                        if shared1 > shared {
                            other_id = *second_id;
                            shared = shared1;
                        }
                    }
                    else {
                        small_hs.insert(*first_id);
                        merge_into.insert(*first_id, (*second_id, shared1));
                    }
                }
            }
            else if shared2 > min_shared_minis{
                //shared2 but not shared1
                //TODO: merge shared2 into shared1
                if merge_into.contains_key(&second_id) {
                    let mut value = merge_into.get_mut(&second_id).unwrap();
                    let mut other_id= value.0;
                    let mut shared = value.1;
                    if shared2 > shared {
                        other_id = *first_id;
                        shared = shared2;

                    }
                }
                else{
                    small_hs.insert(*second_id);
                    merge_into.insert(*second_id,(*first_id,shared2));
                }
            }
        }
    }
}



fn detect_overlaps(cl_set_map: FxHashMap<i32,FxHashSet<u64>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, merge_into: &mut Vec<(i32,i32)>,min_shared_minis: f64, small_hs: &mut FxHashSet<i32>){
    for (cl_id,hashes) in &cl_set_map{
        let mut shared_seed_infos= FxHashMap::default();
        for hash in hashes.iter() {
            if let Some(belongs_to) = cluster_map.get(&hash) {
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
        (shared,most_shared_cluster) = detect_whether_shared_other(*cl_id,min_shared_minis, &shared_seed_infos, hashes.clone());
        if shared{

            //println!("Nr elems_this {}, other {}", nr_elems_this,nr_elems_other);
            if  hashes.len() < cl_set_map.get(&most_shared_cluster).unwrap().len() {
                merge_into.push((*cl_id, most_shared_cluster));
                small_hs.insert(*cl_id);
            }
            else{
                merge_into.push(( most_shared_cluster,*cl_id));
                small_hs.insert(most_shared_cluster);
            }

        }
    }
}


fn detect_whether_shared_other(this_id: i32,min_shared_minis:f64, shared_seed_infos: &FxHashMap<i32,i32>, minimizers:  FxHashSet<u64>) -> (bool, i32) {
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
        if shared_perc > min_shared_minis && shared_perc > most_shared && !(this_id==*key){//} && *nr_shared >=0 {
            most_shared = shared_perc;
            most_shared_cluster = *key;
            if !shared {
                shared = true;
            }
        }
    }
    (shared,most_shared_cluster)
}


fn merge_clusters_from_merge_into_new(merge_into: &mut Vec<(i32,i32)>, clusters_map: &mut  FxHashMap<u64, Vec<i32>>, clusters: &mut FxHashMap<i32,Vec<i32>>, cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, small_hs: FxHashSet<i32>){
    //println!("Merge_into_len: {}",merge_into.len());
    let mut not_mergeable_cter=0;
    for (id , value) in merge_into{
        let mut large_id = value;
        //we might already have deleted large_id from clusters during this iteration
        if clusters.contains_key(&large_id) {
            //idea here: we merge the ids into larger clusters first, smaller clusters are still bound to merge into the new cluster later
            if !small_hs.contains(&large_id){
                merge_clusters( clusters, clusters_map, cl_set_map,&large_id,id)
            }
            else{
                not_mergeable_cter += 1;
            }
        }
    }
    //println!("Not mergeable {}",not_mergeable_cter);
}

fn merge_clusters_from_merge_into(merge_into: &mut FxHashMap<i32,(i32,f64)>, clusters_map: &mut  FxHashMap<u64, Vec<i32>>, clusters: &mut FxHashMap<i32,Vec<i32>>, cl_set_map: &mut FxHashMap<i32,FxHashSet<u64>>, small_hs: FxHashSet<i32>){

    for (id , value) in merge_into{
        let mut large_id = value.0;
        //we might already have deleted large_id from clusters during this iteration
        if clusters.contains_key(&large_id) {
            //idea here: we merge the ids into larger clusters first, smaller clusters are still bound to merge into the new cluster later
            if !small_hs.contains(&large_id){
                merge_clusters( clusters, clusters_map, cl_set_map,&large_id,id)
            }


        }
    }

}

pub(crate) fn post_clustering(clusters: &mut FxHashMap<i32,Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, min_shared_minis:f64){
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    let mut cl_set_map: FxHashMap<i32,FxHashSet<u64>> = FxHashMap::default();
    //TODO; use the new ds instead of cl_overlaps to hopefully reduce RAM usage significantly
    //let mut cl_overlaps: FxHashMap<i32,Vec<i32,i32>> = FxHashMap::default();
    let mut merge_into: Vec<(i32,i32)> = vec![];
    //small_hs is a HashSet storing all cluster ids that were merged into other clusters during this iteration
    let mut small_hs: FxHashSet<i32> = FxHashSet::default();
    //used to have do-while structure
    let mut first_iter= true;
    //continue merging as long as we still find clusters that we may merge
    while (merge_into.len() > 0) || first_iter {
        //println!("Post Cluster iter");
        println!("Merge into : {}", merge_into.len());
        println!("Nr clusters: {}", clusters.len());
        //clear merge_into as this is the indicator how often we attempt to merge further (the while loop depends on it)
        merge_into.clear();
        cl_set_map.clear();
        small_hs.clear();
        //set first_iter to be false to not stay in a infinity loop
        first_iter = false;
        //generate the data structure giving us merge infos
        generate_post_clustering_ds_other(&mut cl_set_map,  cluster_map);
        //print!("pc_ds_generated\n");
        //println!("cl_overlaps_len {}",cl_overlaps.len());
        detect_overlaps(cl_set_map.clone(), cluster_map, &mut merge_into, min_shared_minis, &mut small_hs);
        //println!("# Elements in small_hs {}",small_hs.len());
        //merge_into contains the information about which clusters to merge into which
        //fill_merge_into(&mut cl_overlaps, &mut merge_into, min_shared_minis, &mut cl_set_map, &mut small_hs);
        //println!("merge into filled");
        //cl_overlaps.clear();
        //println!("{:?}",merge_into);
        //merges the clusters
        merge_clusters_from_merge_into_new(&mut merge_into, cluster_map, clusters, &mut cl_set_map, small_hs.clone());
        //it_id += 1;
    }
}
pub(crate) fn post_clustering_new(clusters: &mut FxHashMap<i32,Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>, min_shared_minis:f64){
    //cl_set_map is a hashmap with cl_id -> Hashset of seed hashes
    let mut cl_set_map: FxHashMap<i32,FxHashSet<u64>> = FxHashMap::default();
    //TODO; use the new ds instead of cl_overlaps to hopefully reduce RAM usage significantly
    let mut cl_overlaps: FxHashMap<i32,FxHashMap<i32,i32>> = FxHashMap::default();
    let mut merge_into: FxHashMap<i32,(i32,f64)> = FxHashMap::default();
    //small_hs is a HashSet storing all cluster ids that were merged into other clusters during this iteration
    let mut small_hs: FxHashSet<i32>= FxHashSet::default();
    //used to have do-while structure
    let mut first_iter= true;
    //continue merging as long as we still find clusters that we may merge
    while (merge_into.len() > 0) || first_iter {
        //println!("Post Cluster iter");
        //println!("Merge into : {}", merge_into.len());
        //println!("Nr clusters: {}", clusters.len());
        //clear merge_into as this is the indicator how often we attempt to merge further (the while loop depends on it)
        merge_into.clear();
        cl_set_map.clear();
        small_hs.clear();
        //set first_iter to be false to not stay in a infinity loop
        first_iter = false;
        //generate the data strucutre giving us merge infos
        //fills cl_set_map: cl:id->hashset<seed_hashes> and cl_overlaps: FxHashMap<cl_id1,FxHashMap<cl_id2,overlap>>
        generate_post_clustering_ds(&mut cl_set_map, &mut cl_overlaps, cluster_map);
        //print!("pc_ds_generated\n");
        //println!("cl_overlaps_len {}",cl_overlaps.len());
        //merge_into contains the information about which clusters to merge into which
        fill_merge_into(&mut cl_overlaps, &mut merge_into, min_shared_minis, &mut cl_set_map, &mut small_hs);
        //println!("merge into filled");
        cl_overlaps.clear();
        //println!("{:?}",merge_into);
        //merges the clusters
        merge_clusters_from_merge_into(&mut merge_into, cluster_map, clusters, &mut cl_set_map, small_hs.clone());
        //it_id += 1;
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
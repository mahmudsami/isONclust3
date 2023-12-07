use crate::structs:: FastaRecord;
use std::collections::HashMap;
use crate::generate_sorted_fastq_new_version::Minimizer;
use crate::generate_sorted_fastq_new_version;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::borrow::Borrow;


fn fill_first_cluster(clusters: &mut HashMap<i32,Vec<i32>>, id: i32, sign_minis:&Vec<Minimizer>, mut cluster_map:  HashMap<u64, Vec<i32>>){
    let init_id = 0;
    for minimizer in sign_minis {
        //we get the minimizer sequence from the minimizer object
        let mini_seq = &minimizer.sequence;
        //calculate the hash of the minimizer sequence
        let mini_hash = calculate_hash(&mini_seq);
        //fill cluster_map with the minimizers that we found in the first read
        cluster_map
            .entry(mini_hash)
            .or_insert_with(Vec::new)
            .retain(|&existing_id| existing_id != init_id);
        // Check if id was retained (not a duplicate) and push it if needed
        let vect = cluster_map.get_mut(&mini_hash).unwrap();
        if !vect.contains(&init_id) {
            vect.push(init_id);
        }
    }
    let id_vec=vec![id];
    clusters.insert(init_id,id_vec);
}



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



//clustering method for the case that we do not have any annotation to compare the reads against
pub(crate) fn cluster_de_novo(sorted_entries: Vec<(i32,Vec<Minimizer>)>,min_shared_minis:i32) -> HashMap<i32,Vec<i32>>{
    //clusters contains the main result we are interested in: it will contain the cluster id as key and the read_ids of reads from the cluster as value
    let mut clusters: HashMap<i32,Vec<i32>>=HashMap::new();
    //cluster_map contains a hashmap in which we have a hash_value for a minimizer as key and a vector of read ids as a value
    let mut cluster_map: HashMap<u64, Vec<i32>> = HashMap::new();
    let mut cl_id=1;

    //entry represents a read in our data
    for entry in &sorted_entries{
        //shared mini_infos contains the cluster (key) as well as the number of minimizers appointed to it (value)
        let mut shared_mini_infos=HashMap::new();
        println!("Clusters len {}",clusters.len());
        let id = entry.0;
        let sign_minis= &entry.1;
        //println!("id: {}",id);

        //if we already have at least one cluster: compare the new read to the cluster(s)
        if clusters.len() > 0{
            let mut mini_hashs_vec=vec![];
            for minimizer in sign_minis {
                //we get the minimizer sequence from the minimizer object
                let mini_seq = &minimizer.sequence;
                //calculate the hash of the minimizer sequence
                let mini_hash = calculate_hash(&mini_seq);

                //if we find the minimizer hash in cluster_map: store the clusters in belongs_to
                if let Some(belongs_to) = cluster_map.get(&mini_hash){
                    //iterate over belongs_to to update the counts of shared minimizers for each cluster
                    for &belong_cluster in belongs_to {
                        //if the cluster is already appointed to the read
                        if let Some(new_val) = shared_mini_infos.get(&belong_cluster){
                            //increase the counter of minimizers shared with cluster
                            shared_mini_infos.insert(belong_cluster, *new_val+1);
                        }
                        else{
                            //add a new counter for the cluster element
                            shared_mini_infos.insert(belong_cluster, 0);
                        }
                    }
                }
                mini_hashs_vec.push(mini_hash);
            }
            let mut shared= false;
            let mut most_shared=0;
            let mut most_shared_cluster=0;
            //key: cluster_id, value: count of shared minimizers
            for (key, value) in shared_mini_infos{
                //we have more shared minis with the cluster than our threshold
                if value > min_shared_minis{
                    println!("{} for old clusters" ,cl_id);
                    //add the new minimizers to the cluster
                    for mini_hash in &mini_hashs_vec {
                        cluster_map
                            .entry(*mini_hash)
                            .or_insert_with(Vec::new)
                            .retain(|&existing_id| existing_id != key);
                        // Check if id was retained (not a duplicate) and push it if needed
                        let vect = cluster_map.get_mut(&mini_hash).unwrap();
                        if !vect.contains(&key) {
                            vect.push(key);
                        }
                    }
                    // add the new read_id to the cluster
                    if value > most_shared{
                        most_shared_cluster = key;
                    }
                    shared = true;
                }
                //println!("id: {} , shared_with_init_cl: {}",key, value);
            }
            //if we have a cluster that we share enough minimizers with
            if shared{
                //add the read id to read_list
                let read_list= clusters.get_mut(&most_shared_cluster).unwrap();
                if !read_list.contains(&id){
                    read_list.push(id);
                }
            }
            //we did not find a cluster that we could put the read into-> generate a new cluster
            else{
                for mini_hash in mini_hashs_vec{
                    cluster_map
                        .entry(mini_hash)
                        .or_insert_with(Vec::new)
                        .retain(|&existing_id| existing_id != cl_id);
                    // Check if id was retained (not a duplicate) and push it if needed
                    let vect = cluster_map.get_mut(&mini_hash).unwrap();
                    if !vect.contains(&cl_id) {
                        vect.push(cl_id);
                    }
                }
                let id_vec=vec![id];
                clusters.insert(cl_id,id_vec);
                cl_id+=1;
            }
        }
        //we do not yet have a cluster and therefore need to fill the first read into the first
        else{
            //fill_first_cluster(&mut clusters, id, sign_minis, cluster_map);
            let init_id = 0;
            for minimizer in sign_minis {
                //we get the minimizer sequence from the minimizer object
                let mini_seq = &minimizer.sequence;
                //calculate the hash of the minimizer sequence
                let mini_hash = calculate_hash(&mini_seq);
                //fill cluster_map with the minimizers that we found in the first read
                cluster_map
                    .entry(mini_hash)
                    .or_insert_with(Vec::new)
                    .retain(|&existing_id| existing_id != init_id);
                // Check if id was retained (not a duplicate) and push it if needed
                let vect = cluster_map.get_mut(&mini_hash).unwrap();
                if !vect.contains(&init_id) {
                    vect.push(init_id);
                }
            }
            let id_vec=vec![id];
            clusters.insert(init_id,id_vec);
        }
    }
clusters
}



fn get_final_cl_init(this_clusters: Vec<&Vec<i32>>) -> i32 {
    //let mut final_cl=vec![];
    let mut clustering: HashMap<i32,f32> = HashMap::new();
    for tcl in this_clusters{
        if tcl.len()==1{
            let tcl_elem= tcl.first().unwrap();
            println!("tcl elem {}",tcl_elem);
            if clustering.contains_key(tcl_elem){
                let mut elem=clustering.get_mut(tcl_elem).unwrap();
                *elem = *elem + 1.0;
            }
            else{
                clustering.insert(*tcl_elem,1.0);
            }
        }
    }
    let mut max_val:f32 = 0.0;
    let mut max_clust=0;
    let mut max_eq=0;
    for (key,val) in clustering{
        if val>max_val{
            max_val=val;
            max_clust=key;
            max_eq=0;
        }
        if val==max_val{
            max_eq=key;
        }
    }
    if !max_eq==0{
        println!("Not a singular cluster: {},{}", max_clust,max_eq)
    }
    if max_val<10.0{
        println!("low support {}",max_val);
    }

    return max_clust
}



//clustering method for the case that we have an annotation to compare our reads against
pub(crate) fn cluster_from_initial(sorted_entries: Vec<(i32,Vec<Minimizer>)>, initial_clusters: HashMap<u64, Vec<i32>>) ->HashMap<i32, Vec<i32>>{
    //the hashmap containing the clusters we found
    let mut clusters: HashMap<i32, Vec<i32>> = HashMap::new();
    for sorted_entry in sorted_entries{
        let mut this_clusters= vec![];
        for mini in sorted_entry.1{
            let this_mini_hash= calculate_hash(&mini.sequence);
            if initial_clusters.contains_key(&this_mini_hash){
                let value= initial_clusters.get(&this_mini_hash).unwrap();
                //this_clusters.append(&mut value.clone());
                //we append the read to this
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
        let final_cl = get_final_cl_init(this_clusters);
        if clusters.contains_key(&final_cl){
            let mut id_vec=clusters.get_mut(&final_cl).unwrap();
            id_vec.push(sorted_entry.0)
        }
        else {
            let mut id_vec=vec![sorted_entry.0];
            clusters.insert(final_cl,id_vec);
        }

    }

    clusters
    }



pub(crate) fn add_rev_comp_seqs_annotation(initial_clustering_records:Vec<FastaRecord>) ->Vec<FastaRecord>{
    ///Adds the reverse_compliment of the annotation
    ///INPUT:   initial_clustering_records: The annotation
    ///
    /// OUTPUT:     both_dir_records: the records in both directions
    let mut both_dir_records =vec![];
    let mut init_len=initial_clustering_records.len();
    for record in initial_clustering_records{
        init_len= init_len+1;
        //let reversed: String = record.sequence.chars().rev().collect();
        let reversed=reverse_complement(&*record.sequence);
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
        init_len= init_len+1;
        //let reversed: String = record.sequence.chars().rev().collect();
        let reversed=reverse_complement(&*record.sequence);
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
        let this_minimizers= generate_sorted_fastq_new_version::get_canonical_kmer_minimizers(&*record.sequence, k, window_size);
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
use crate::generate_sorted_fastq_new_version::{ Minimizer};
use std::collections::HashSet;


enum Cluster<T, U> {
    read_ids(HashSet<T>),
    mini_seqs(HashSet<U>),
}

pub(crate) fn cluster_sorted_entries(sorted_entries: Vec<(i32,Vec<Minimizer>)>) ->Vec<Vec<i32>>{
    let mut clusters= vec![];
    for entry in sorted_entries{
        let id = entry.0;
        let sign_minis=entry.1;
        /*if entry==sorted_entries.0{
            //let mut cluster0=Cluster { i32,String };
            for sign_mini in sign_minis{
                //cluster0.mini_seqs.insert(sign_mini.sequence)
            }
        }*/


    }
    clusters
}
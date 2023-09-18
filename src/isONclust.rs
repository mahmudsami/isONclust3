use clap::{ArgMatches};
use std::fs;
use std::process::Command as RunCommand;


pub fn index_kmers(sub_matches: &ArgMatches) {
   // let (kmtricks_path, kmindex_path) = get_km_paths(sub_matches);
    //let inkmers = sub_matches.get_one::<String>("IN_KMERS").map(|s| s.clone()).unwrap();
    //let outkmers = sub_matches.get_one::<String>("OUT_INDEX").map(|s| s.clone()).unwrap();

    let k = sub_matches.get_one::<u32>("K").map(|s| s.clone()).unwrap();
    println!("K: {}",k);
    let t = sub_matches.get_one::<u32>("T").map(|s| s.clone()).unwrap();
    println!("T: {}",t);
}
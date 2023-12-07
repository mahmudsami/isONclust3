use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufReader, BufRead, Error};
use std::path::Path;
use std::collections::HashMap;
use std::fs;
use crate::structs::{FastqRecord, FastqRecord_isoncl_init};
use std::collections::hash_map::RandomState;



fn write_final_clusters_tsv(outfolder: String, clusters: HashMap<i32,Vec<i32>>, id_map:HashMap<i32,String>)->HashMap<String,i32>{
    let mut header_cluster_map=HashMap::new();
    let file_path = PathBuf::from(outfolder).join("final_clusters.tsv");
    let mut f = File::create(file_path).expect("unable to create file");
    println!("{} different clusters identified",clusters.len());
    for (cl_id, r_int_ids) in clusters.into_iter(){
        println!("cl_id {}, nr_reads {}",cl_id,r_int_ids.len());
        for r_int_id in r_int_ids{
            let read_id = id_map.get(&r_int_id).unwrap().to_string();
            //write!(f ,"{}\t{}\n", cl_id, read_id);
            header_cluster_map.insert(read_id,cl_id);
        }
    }
    header_cluster_map
}



fn create_final_ds(header_cluster_map: HashMap<String,i32>, fastq_vec: Vec<FastqRecord_isoncl_init>)->HashMap<i32,Vec<FastqRecord_isoncl_init>>{
    let mut cluster_map= HashMap::new();
    for read in fastq_vec{
        let id =read.header.clone();
        //println!("id {}",id);
        if header_cluster_map.contains_key(&id) {
            let cluster_id = header_cluster_map.get(&id).unwrap();
            if cluster_map.contains_key(cluster_id) {
                let mut id_vec: &mut Vec<FastqRecord_isoncl_init> = cluster_map.get_mut(cluster_id).unwrap();
                id_vec.push(read)
            } else {
                let mut id_vec = vec![read];
                cluster_map.insert(*cluster_id, id_vec);
            }
        }
    }
    cluster_map
}



fn write_fastq_files(outfolder: &Path, cluster_map: HashMap<i32, Vec<FastqRecord_isoncl_init>>){
   //Writes the fastq files using the data structure cluster_map HashMap<i32, Vec<FastqRecord_isoncl_init>>

    for (cl_id, records) in cluster_map.into_iter(){
        let filename = cl_id.to_string()+".fastq";
        let file_path = PathBuf::from(outfolder).join(filename);
        let mut f = File::create(file_path).expect("unable to create file");
        for record in records{
            write!(f ,"{} \n {} \n + \n {} \n", record.header, record.sequence,record.quality);
        }
            //let read_id=id_map.get(&r_int_id);
        //let index = test.iter().position(|&r| r == "two").unwrap();
            //


    }
}



pub fn path_exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}
fn convert_infos_for_writing(id_map:HashMap<i32,String>, clusters:HashMap<i32,Vec<i32>>, fastq_vec:&Vec<FastqRecord_isoncl_init>){



}


pub(crate) fn write_output(outfolder:String,clusters:HashMap<i32,Vec<i32>>,fastq_vec:Vec<FastqRecord_isoncl_init>, id_map:HashMap<i32,String>){
    if !path_exists(&*outfolder){
        fs::create_dir(outfolder.clone());
    }
    let fastq_path=Path::new(&outfolder).join("fastq_files");
    if !fastq_path.exists(){
        fs::create_dir(fastq_path.clone());
    }
    convert_infos_for_writing(id_map.clone(), clusters.clone(), &fastq_vec);
    let header_cluster_map=write_final_clusters_tsv(outfolder,clusters.clone(),id_map.clone());
    let cluster_hashmap_fastq_record = create_final_ds(header_cluster_map, fastq_vec);
    write_fastq_files(&*fastq_path, cluster_hashmap_fastq_record);
}
use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufReader, BufRead, Error};
use std::path::Path;
use std::collections::HashMap;
use std::fs;
use crate::structs::{FastqRecord, FastqRecord_isoncl_init};
use std::collections::hash_map::RandomState;



fn write_final_clusters_tsv(outfolder: String, clusters: HashMap<i32,Vec<i32>>, id_map:HashMap<i32,String>){
    let file_path = PathBuf::from(outfolder).join("final_clusters.tsv");
    let mut f = File::create(file_path).expect("unable to create file");
    for (cl_id, r_int_ids) in clusters.into_iter(){
        for r_int_id in r_int_ids{
            let read_id=id_map.get(&r_int_id);
            write!(f ,"{}\t{}\n", cl_id, read_id.unwrap());
        }
    }
}



/*fn write_fastq_files(outfolder: &Path, clusters: HashMap<i32,Vec<i32>>, fastq_vec: Vec<FastqRecord_isoncl_init>){
    for (cl_id, r_int_ids) in clusters.into_iter(){
        let filename = cl_id.to_string()+".fastq";
        let file_path = PathBuf::from(outfolder).join(filename);
        let mut f = File::create(file_path).expect("unable to create file");

            //let read_id=id_map.get(&r_int_id);
        //let index = test.iter().position(|&r| r == "two").unwrap();
            //write!(f ,"{}\t{}\n", cl_id, read_id.unwrap());


    }
}*/



pub fn path_exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}
fn convert_infos_for_writing(id_map:HashMap<i32,String>,clusters:HashMap<i32,Vec<i32>>,fastq_vec:Vec<FastqRecord_isoncl_init>){
    for read in fastq_vec{
        
    }

}


pub(crate) fn write_output(outfolder:String,clusters:HashMap<i32,Vec<i32>>,fastq_vec:Vec<FastqRecord_isoncl_init>, id_map:HashMap<i32,String>){
    if !path_exists(&*outfolder){
        fs::create_dir(outfolder.clone());
    }
    let fastq_path=Path::new(&outfolder).join("fastq_files");
    if !fastq_path.exists(){
        fs::create_dir(fastq_path.clone());
    }
    convert_infos_for_writing(id_map.clone(),clusters.clone(),fastq_vec);
    write_final_clusters_tsv(outfolder,clusters.clone(),id_map.clone());
    //write_fastq_files(&*fastq_path, clusters, fastq_vec);
}
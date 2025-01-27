use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::path::Path;
use std::fs;
use crate::structs::FastqRecord_isoncl_init;
use rustc_hash::FxHashMap;
use crate::{Cluster_ID_Map, file_actions};
use rayon::{prelude::*, vec};


pub(crate) fn write_ordered_fastq(score_vec: &Vec<(i32,usize)>, outfolder: &String,id_map: &FxHashMap<i32,String>,fastq: &str){
    //writes the fastq file
    let _ = fs::create_dir_all(PathBuf::from(outfolder).join("clustering"));
    let fastq_file = File::open(fastq).unwrap();
    let mut fastq_records= FxHashMap::default();
    file_actions::parse_fastq_hashmap(fastq_file,&mut fastq_records);
    let f = File::create(outfolder.to_owned()+"/clustering/sorted.fastq").expect("Unable to create file");
    let mut buf_write = BufWriter::new(&f);
    //for record in fastq_records {
    for score_tup in score_vec.iter(){
        let this_header = id_map.get(&score_tup.0).unwrap();
        let record= fastq_records.get(this_header).unwrap();
        if &record.header == this_header{
            write!(buf_write, "@{}\n{}\n+\n{}\n", record.get_header(), record.get_sequence(),record.get_quality()).expect("Could not write file");
        }
    }
    buf_write.flush().expect("Failed to flush the buffer");
}

pub(crate) fn write_ordered_fastq_offset(score_vec: &Vec<(i32,usize)>, outfolder: &String,id_map: &FxHashMap<i32,String>,fastq: &str, num_chunks: usize, create: bool, step: usize, rem: usize ){
    //writes the fastq file
    if create {
        let _ = fs::create_dir_all(PathBuf::from(outfolder).join("clustering"));
    }
    //let mut fastq_records= FxHashMap::default();
    let f = if create {
        File::create(outfolder.to_owned()+"/clustering/sorted.fastq").expect("Unable to create file")
    } else {
        std::fs::OpenOptions::new().append(true).open(outfolder.to_owned()+"/clustering/sorted.fastq").expect("Unable to create file")
    };
    let mut global_records: Vec<Vec<FastqRecord_isoncl_init>> = vec![Vec::new();num_chunks];
    let mut max = vec![0;num_chunks];
    let mut current = vec![0;num_chunks];

    for i in 0..num_chunks{
        let start = i*step;
        let size = if i == num_chunks-1 {rem} else {step};
        let mut records = FxHashMap::default();
        file_actions::parse_fastq_hashmap_offset(fastq,&mut records,start,size);
        for score_tup in score_vec.iter(){
            let this_header = id_map.get(&score_tup.0).unwrap();
            if records.get(this_header).is_some(){
                let record:&FastqRecord_isoncl_init= records.get(this_header).unwrap();
                global_records[i].push((*record).clone());
                max[i] += 1;
            }
        }
    }

    let mut buf_write = BufWriter::new(&f);
    for (read_id,_) in score_vec.iter(){
        for i in 0..num_chunks{
            if current[i] < max[i]{
                if global_records[i][current[i]].header == *id_map.get(&read_id).unwrap(){
                    write!(buf_write, "@{}\n{}\n+\n{}\n", global_records[i][current[i]].get_header(), global_records[i][current[i]].get_sequence(),global_records[i][current[i]].get_quality()).expect("Could not write file");
                    current[i] += 1;
                    break;
                }
            }
        }
    }
    buf_write.flush().expect("Failed to flush the buffer");
    

}

fn write_final_clusters_tsv(outfolder: &Path, clusters: Cluster_ID_Map, id_map:FxHashMap<i32,String>, header_cluster_map:&mut FxHashMap<String,i32>){
    let file_path = PathBuf::from(outfolder).join("final_clusters.tsv");

    let f = File::create(file_path).expect("unable to create file");
    let mut buf_write = BufWriter::new(&f);
    let mut nr_reads= 0;
    println!("{} different clusters identified",clusters.len());
    //let nr_clusters=clusters.len();
    for (cl_id, r_int_ids) in clusters.into_iter(){
        //println!("cl_id {}, nr_reads {:?}",cl_id,nr_reads);
        for r_int_id in r_int_ids{
            let read_id = id_map.get(&r_int_id).unwrap();
            nr_reads += 1;
            let _ = writeln!(buf_write ,"{}\t{}", cl_id, read_id);
            header_cluster_map.insert(read_id.to_string(),cl_id);
        }
    }
    // Flush the buffer to ensure all data is written to the underlying file
    buf_write.flush().expect("Failed to flush the buffer");
    //println!("{} different clusters identified",nr_clusters);
    //println!("HCLM {:?}",header_cluster_map);
    println!("{} reads added to tsv file", nr_reads);
}


//TODO: this is the current RAM bottleneck: we read the whole file to then have the reads when we write the output
//Outline: sort the fastq file by cluster and then write the entries from the sorted fastq file to not having to read the full file
fn create_final_ds(header_cluster_map: FxHashMap<String,i32>, fastq: String, cluster_map: &mut FxHashMap<i32,Vec<FastqRecord_isoncl_init>>){
    let fastq_file = File::open(fastq).unwrap();
    let mut fastq_vec= vec![];
    //parse the fastq file to store the data in fastq_vec
    file_actions::parse_fastq(fastq_file,&mut fastq_vec);
    //iterate over fastq_vec and add the reads to cluster_map
    for read in fastq_vec{
        let id = read.header.clone();
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
}

fn create_final_ds_offset(header_cluster_map: &FxHashMap<String,i32>, fastq: String, cluster_map: &mut FxHashMap<i32,Vec<FastqRecord_isoncl_init>>, num_chunks: usize , step: usize, rem: usize){
    for i in 0..num_chunks{
        let start = i*step;
        let size = if i == num_chunks-1 {rem} else {step};
        let mut fastq_vec= vec![];
        //parse the fastq file to store the data in fastq_vec
        file_actions::parse_fastq_offset(fastq.clone(),&mut fastq_vec,start,size);
        //iterate over fastq_vec and add the reads to cluster_map
        for read in fastq_vec{
            let id = read.header.clone();
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
    }

}


fn write_fastq_files(outfolder: &Path, cluster_map: FxHashMap<i32, Vec<FastqRecord_isoncl_init>>, n: usize){
    let mut new_cl_id = 0;
    let mut read_cter= 0;
    //fs::create_dir_all(PathBuf::from(outfolder).join("fastq_files"));
    let fastq_outfolder= PathBuf::from(outfolder);
    //Writes the fastq files using the data structure cluster_map HashMap<i32, Vec<FastqRecord_isoncl_init>>
    for (cl_id, records) in cluster_map.into_iter(){
        if records.len() >= n { //only write the records if we have n or more reads supporting the cluster
            let filename = new_cl_id.to_string()+".fastq";
            let file_path = fastq_outfolder.join(filename);
            let f = File::create(file_path).expect("unable to create file");
            let mut buf_write = BufWriter::new(&f);
            for record in records{
                write!(buf_write ,"@{}\n{}\n+\n{}\n", record.header, record.sequence,record.quality).expect("We should be able to write the entries");
                read_cter += 1;
            }
            buf_write.flush().expect("Failed to flush the buffer");
            new_cl_id += 1; //this is the new cl_id as we skip some on the way
        }


        //println!("cl id for writing: {}, {}",cl_id,read_cter);
    }
    println!("{} reads written",read_cter);
}



pub fn path_exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}


fn write_fastq_files_offset(outfolder: &Path, cluster_map: &FxHashMap<i32, Vec<FastqRecord_isoncl_init>>, n: usize, num_part: usize, part: usize){
    let mut new_cl_id = part;
    let mut read_cter= 0;
    //fs::create_dir_all(PathBuf::from(outfolder).join("fastq_files"));
    let fastq_outfolder= PathBuf::from(outfolder);
    //Writes the fastq files using the data structure cluster_map HashMap<i32, Vec<FastqRecord_isoncl_init>>
    for (cl_id, records) in cluster_map.into_iter(){
        if records.len() >= n { //only write the records if we have n or more reads supporting the cluster
            let filename = new_cl_id.to_string()+".fastq";
            let file_path = fastq_outfolder.join(filename);
            let f = File::create(file_path).expect("unable to create file");
            let mut buf_write = BufWriter::new(&f);
            for record in records{
                write!(buf_write ,"@{}\n{}\n+\n{}\n", record.header, record.sequence,record.quality).expect("We should be able to write the entries");
                read_cter += 1;
            }
            buf_write.flush().expect("Failed to flush the buffer");
            new_cl_id += num_part; //this is the new cl_id as we skip some on the way
        }


        //println!("cl id for writing: {}, {}",cl_id,read_cter);
    }
    println!("{} reads written",read_cter);
}

pub(crate) fn write_output(outfolder: String, clusters: &Cluster_ID_Map, fastq: String, id_map: FxHashMap<i32,String>, n: usize, no_fastq: bool, memory_restrict: bool, p: usize, c: usize){

    if !path_exists(&outfolder){
        let _ = fs::create_dir(outfolder.clone()).expect("We should be able to create the directory");
    }
    //clustering_path: the outfolder of isONclust3
    let clustering_path= Path::new(&outfolder).join("clustering");
    if !clustering_path.exists(){
        let _ = fs::create_dir(clustering_path.clone());
    }
    //the subfolder of clustering in which we write the fastq_files ( as done in isONclust1)
    let fastq_path= clustering_path.join("fastq_files");
    if !fastq_path.exists(){
        let _ = fs::create_dir(fastq_path.clone());
    }
    let mut cluster_hashmap_fastq_record= FxHashMap::default();
    let mut header_cluster_map= FxHashMap::default();
    write_final_clusters_tsv(&clustering_path, clusters.clone(), id_map.clone(), &mut  header_cluster_map);
    //no_fastq: true -> we do not want to write the fastq files
    if !no_fastq{
        if memory_restrict {
            let num_part = p;
            let num_chunks = c;
            let total_size = id_map.len();
            let step = total_size / num_chunks;
            let rem = total_size - (total_size / num_chunks) * (num_chunks-1);
            let header_cluster_maps = partition_cluster_map(&header_cluster_map,num_part);
            for i in 0..num_part{
                create_final_ds_offset(&header_cluster_maps[i], fastq.clone(), &mut cluster_hashmap_fastq_record, num_chunks, step, rem);
                write_fastq_files_offset(&fastq_path, &cluster_hashmap_fastq_record, n, num_part, i);
                cluster_hashmap_fastq_record.clear();
            }
        }else{
            //create a data structure that we use to generate the proper fastq files
            create_final_ds(header_cluster_map, fastq,&mut cluster_hashmap_fastq_record);
            //println!("Cluster_hashmap: {}",cluster_hashmap_fastq_record.len());
            println!("Writing the fastq files");
            write_fastq_files(&fastq_path, cluster_hashmap_fastq_record, n);
        }
    }
}

fn partition_cluster_map(header_cluster_map: &FxHashMap<String,i32>, num_part: usize) -> Vec<FxHashMap<String,i32>>{
    let mut header_cluster_maps = vec![FxHashMap::default();num_part];
    for (header,cl_id) in header_cluster_map.iter(){
        header_cluster_maps[(cl_id%num_part as i32) as usize].insert(header.clone(),*cl_id);
    }
    header_cluster_maps
}
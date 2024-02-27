use crate::structs::{GtfEntry, FastaRecord, Coord_obj};
use std::fs::File;
use std::io::{BufReader, BufRead, Read};
use std::str::FromStr;
use rustc_hash::{FxHashMap, FxHasher, FxHashSet};
use bio::io::gff;
use bio::io::gff::{GffType, Record};
use bio::io::gff::GffType::GFF3;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::path::Path;
use bio::io::fasta;
use bio::io::fasta::FastaRead;
extern crate rayon;
use crate::generate_sorted_fastq_new_version;
use crate::clustering;
use std::time::Instant;

//TODO: add overlap detection
//TODO: remove multiple occasions of minimizers in the same gene if the exons overlap
//TODO: possibly use the gene id for cluster_identification

fn detect_overlaps(gene_map: &FxHashMap<i32,Vec<Coord_obj>>, this_gene_id:&i32, this_coords: &Vec<Coord_obj>, overlap_ctr: &mut i32){
    for (gene_id_other,coords) in gene_map{
        if gene_id_other > this_gene_id{
            for coord in coords{
                for this_coord in this_coords{
                    if coord.overlaps_with(this_coord){
                        *overlap_ctr +=1;
                    }
                }
            }
        }
    }
}



fn parse_fasta_and_gen_clusters(fasta_path: Option<&str>, coords: FxHashMap<String,FxHashMap<i32,Vec<Coord_obj>>>,clusters: &mut FxHashMap<i32, Vec<i32>>, init_cluster_map: &mut FxHashMap<u64, Vec<i32>>,k: usize,w: usize){
    println!("parse_fasta");
    let path=fasta_path.unwrap();
    let mut reader = fasta::Reader::from_file(Path::new(path)).expect("We expect the file to exist");
    //let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    reader.records().into_iter().for_each(|record|{

        let mut internal_id=0;
        let mut record_minis=vec![];
        //let mut record_seq: String =String::new();
        //retreive the current record
        //println!("It over records");
        let seq_rec = record.expect("invalid record");
        let sequence = std::str::from_utf8(seq_rec.seq()).unwrap().to_uppercase();
        let mut overlap_ctr= 0;
        //let local_seq= std::str::from_utf8(sequence).expect("The genomic sequence should be utf8").to_string();
        //in the next lines we make sure that we have a proper header and store it as string
        let id = seq_rec.id().to_string().split(' ').collect::<Vec<_>>()[0].to_string();
        println!("Now to the coords_ds");
        if let Some(gene_map) = coords.get(id.as_str()) {
            //iterate over the genes in the gene_map
            for (gene_id, exon_coords) in gene_map { //TODO:test whether into_par_iter works here
                let mut coords_in_gene=FxHashSet::default();
                for exon_coord in exon_coords{
                    if ! coords_in_gene.contains(exon_coord){
                        let exon_seq = &sequence[exon_coord.startpos as usize..exon_coord.endpos as usize].to_string();
                        let mut exon_minis=vec![];
                        generate_sorted_fastq_new_version::get_canonical_kmer_minimizers_hashed(exon_seq.as_bytes(),k,w,&mut exon_minis);
                        record_minis.append( &mut exon_minis);
                        coords_in_gene.insert(exon_coord);
                    }
                }
                detect_overlaps(gene_map,gene_id,exon_coords, &mut overlap_ctr);
                //println!("Record_seq {}: {}",id,record_minis);
                clustering::generate_initial_cluster_map(&record_minis, init_cluster_map, *gene_id);
                let id_vec= vec![];
                clusters.insert(*gene_id,id_vec);
                //println!("{:?}",init_cluster_map);
                internal_id += 1;
            }
        }
        println!("{} overlaps between genes (their exons) ",overlap_ctr);
    });
}


fn parse_gtf_and_collect_coords(gtf_path: Option<&str>, coords:&mut FxHashMap<String,FxHashMap<i32,Vec<Coord_obj>>>){
    let reader = gff::Reader::from_file(gtf_path.unwrap(),GFF3);
    let mut gene_id=0;
    let mut coords_in_gene=FxHashSet::default();
    let mut true_gene=false;
    for record in reader.expect("The reader should find records").records() {
        let mut rec = record.ok().expect("Error reading record.");
        //we have a new gene
        if rec.feature_type() == "gene" {//|| rec.feature_type() == "pseudogene"{
            true_gene = true;
            //see if we are in a new chromosome/scaffold
            if !coords.contains_key(rec.seqname()){
                //we are in a new chromosome/scaffold
                //reset the gene_id
                gene_id = 0;
                //
                let sname= rec.seqname().to_string();
                coords.insert(sname,FxHashMap::default());

            }
            //increase the gene_id by 1
            gene_id += 1;
            coords_in_gene.clear();
        }
        else if rec.feature_type()=="exon"{
            //we only are interested in exons from true genes
            if true_gene {
                //println!("{} {} {}",rec.seqname(),rec.feature_type(),gene_id);
                if let Some(gene_map) = coords.get_mut(rec.seqname()) {
                    if let Some(coord_vec) = gene_map.get_mut(&gene_id) {
                        let coord_o = Coord_obj { startpos: *rec.start(), endpos: *rec.end() };
                        if !coords_in_gene.contains(&coord_o) {
                            coord_vec.push(coord_o.clone());
                            coords_in_gene.insert(coord_o);
                        }
                    } else {
                        let coord_o = Coord_obj { startpos: *rec.start(), endpos: *rec.end() };
                        if !coords_in_gene.contains(&coord_o) {
                            let mut coord_vec = vec![];
                            coord_vec.push(coord_o.clone());
                            gene_map.insert(gene_id, coord_vec);
                            coords_in_gene.insert(coord_o);
                        }
                    }
                }
            }
        }
        //we skip any exons found in pseudogenes for now
        else if rec.feature_type()=="pseudogene"{
            true_gene = false;
        }
    }
}


pub(crate) fn resolve_gff(gff_path: Option<&str>, fasta_path: Option<&str>,clusters: &mut FxHashMap<i32, Vec<i32>>, cluster_map: &mut FxHashMap<u64, Vec<i32>>,k:usize ,w:usize) {
    println!("Resolving GFF file ");
    let now1 = Instant::now();
    let mut coords=FxHashMap::default();//: HashMap<K, HashMap<i32, Vec<Coord_obj>, BuildHasherDefault<FxHasher>>, BuildHasherDefault<FxHasher>> = FxHashMap::default();
    parse_gtf_and_collect_coords(gff_path, &mut coords);
    println!("{} s used for parsing the gff file", now1.elapsed().as_secs());
    println!("First step done");
    let now2 = Instant::now();
    parse_fasta_and_gen_clusters(fasta_path,coords, clusters, cluster_map,k, w);
    println!("Generated {} initial clusters from the reference",clusters.len());
    println!("{} s used for parsing the fasta file", now2.elapsed().as_secs());
    //detectOverlaps(coords);
    //for coord in &coords{
    //    println!("id: {}",coord.0);
    //    for coord_e in coord.1{
    //        println!(" {}", coord_e)
    //    }
    //}
    println!("GTF resolved");
}



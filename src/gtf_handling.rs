use crate::structs::{GtfEntry, FastaRecord, Coord_obj};
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::str::FromStr;
use rustc_hash::FxHashMap;
use bio::io::gff;
use bio::io::gff::{GffType, Record};
use bio::io::gff::GffType::GFF3;

fn gtf_parse_strand(entry:String) -> bool {
    let return_bool:bool;
    let char_vec:Vec<char>=entry.chars().collect();
    if char_vec.len()!=1{
        panic!("Expected String to be of length 1, but got {}",entry);

    }
    else {
        let char=char_vec[0];
        if char=='+'{
            return_bool=true
        }
        else if char=='-' {
            return_bool=false
        }
        else {
            panic!("Expected '+' or '-' in this field but got {}",char)
        }
    }
    return_bool
}


fn gtf_parse_frame(entry:String) -> i8 {
    entry.parse::<i8>().unwrap()
}

fn gtf_parse_score(entry:&str) ->f64 {
    let score:f64=f64::from_str(entry).unwrap();
    /*if entry=="."{
        score=0.0;
    }
    else{
        score;
    }*/
    score
}

fn get_coords(rec: Record, mut coords: &mut FxHashMap<String,Vec<Coord_obj>>){
    if rec.feature_type() == "exon"{
        if let Some(x) = coords.get_mut(rec.seqname()) {
            x.push(Coord_obj{startpos: *rec.start(),endpos: *rec.end() });
        }
        else{
            let mut coord_vec=vec![];
            coord_vec.push(Coord_obj{startpos: *rec.start(),endpos: *rec.end() });
            coords.insert(rec.seqname().to_string(), coord_vec);
        }
    }
}
/*fn get_coords( gtf_entries: Vec<GtfEntry>,mut coords: &mut FxHashMap<String,Vec<Coord_obj>>){
    for entry in gtf_entries{
        if let Some(x) = coords.get_mut(entry.seqname.as_str()) {
            x.push(Coord_obj{startpos: entry.start,endpos: entry.end });
        }
        else{
            let mut coord_vec=vec![];
            coord_vec.push(Coord_obj{startpos: entry.start,endpos: entry.end });
            coords.insert(entry.seqname,coord_vec);
        }
    }
}*/


fn get_reads(ref_rec_head: String, ref_rec_seq: String, coords: FxHashMap<String,Vec<Coord_obj>>,new_seqs: Vec<FastaRecord>){

}


fn lineToGTF(line:&str) -> GtfEntry {
    let bound_line=line.to_string();
    let splitline: Vec<&str> = bound_line.split('\t').collect();
    let seqname=splitline[0];
    let source=splitline[1];
    let feature=splitline[2];
    let start:usize =splitline[3].trim()
        .parse()
        .expect("Wanted a number");
    let end:usize =splitline[4].trim()
        .parse()
        .expect("Wanted a number");
    //println!("seqname: {},float: {}",splitline[0],splitline[5]);
    //let score: f64  = gtf_parse_score(splitline[5]);
    let strand: bool = gtf_parse_strand(splitline[6].to_string());
    //let frame: i8 = gtf_parse_frame(splitline[7].to_string());
    //let attribute = splitline[8];
    let gtf_entry = GtfEntry{
        seqname: seqname.to_string(),
        source: source.to_string(),
        feature: feature.to_string(),
        start,
        end,
        //score,
        strand,
        //frame,
        //attribute: attribute.to_string()
    };
    gtf_entry
}


fn parse_gtf(file_path:&str, mut gtfRecords: &mut Vec<GtfEntry>){
    let file = File::open(file_path).expect("The gtf file was not found");
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if !line.as_ref().unwrap().starts_with("#"){
            let gtf_entry = lineToGTF(line.unwrap().as_str());
            gtfRecords.push(gtf_entry)
        }
    }
}


pub(crate) fn resolve_gtf_old(gtf_path: Option<&str>, fastq_path: Option<&str>, outfolder: String){
    let mut gtfRecords =vec![];
    let reference_records: Vec<FastaRecord>=vec![];
    //let mut gtf_entries = vec![];
    let mut coords:FxHashMap<String,Vec<Coord_obj>>=FxHashMap::default();
    if let Some(gtf_path_u) = gtf_path {
        parse_gtf(gtf_path_u, &mut gtfRecords);
        println!("gtf file parsed")
    }
    //get_coords(gtfRecords,&mut coords);
    for coord in &coords{
        for coord_e in coord.1{
            println!("{}: {}", coord.0, coord_e)
        }
    }

    //reference_records.iter().for_each(|reference_record| {
    //    get_coords(*reference_record.header,reference_record.sequence, gtf_coords);
    //});
}
pub(crate) fn resolve_gtf(gtf_path: Option<&str>, fastq_path: Option<&str>, outfolder: String) {
    println!("Resolving GFF file ");
    let mut reader = gff::Reader::from_file(gtf_path.unwrap(),GFF3);
    let mut coords: FxHashMap<String,Vec<Coord_obj>> = FxHashMap::default();
    for record in reader.expect("The reader should find records").records() {
        let rec = record.ok().expect("Error reading record.");
        get_coords(rec, &mut coords);

        //Coord_obj{startpos: *rec.start(),endpos: *rec.end() }
    }
    for coord in &coords{
        println!("id: {}",coord.0);
        for coord_e in coord.1{
            println!(" {}", coord_e)
        }
    }
    println!("GTF resolved")
}



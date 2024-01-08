use std::fs::{File};
use std::io::{BufReader, BufRead};
use std::error::Error;


use crate::structs;
use crate::structs::{FastqRecord_isoncl_init, GtfEntry};
use crate::structs::FastaRecord;
use rayon::iter::split;
use std::str::FromStr;


pub(crate) fn parse_fasta(file_path: &str) -> Result<Vec<FastaRecord>, std::io::Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut current_header = String::new();
    let mut current_sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // New header line
            if !current_header.is_empty() {
                // Store the previous sequence
                records.push(FastaRecord {
                    header: current_header.clone(),
                    sequence: current_sequence.concat(),
                });
                current_sequence.clear();
            }
            current_header = line[1..].trim().to_string(); // Remove '>'
        } else {
            // Sequence line
            current_sequence.push(line);
        }
    }

    // Store the last sequence
    if !current_header.is_empty() {
        records.push(FastaRecord {
            header: current_header.clone(),
            sequence: current_sequence.concat(),
        });
    }

    Ok(records)
}

fn shorten_header(header:&str)-> &str{
    let header_parts: Vec<_>=header.split(" ").collect();
    let header_new = header_parts[0];

    header_new
}


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
    let frame_ind = entry.parse::<i8>().unwrap();
    frame_ind
}
fn gtf_parse_score(entry:&str) ->f64 {
    let score:f64;
    if entry=="."{
        score=0.0;
    }
    else{
        score=f64::from_str(entry).unwrap();
    }
    score
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
    let score: f64  = gtf_parse_score(splitline[5]);
    let strand: bool = gtf_parse_strand(splitline[6].to_string());
    let frame: i8 = gtf_parse_frame(splitline[7].to_string());
    let attribute = splitline[8];
    let gtf_entry=structs::GtfEntry{
        seqname: seqname.to_string(),
        source: source.to_string(),
        feature: feature.to_string(),
        start,
        end,
        score,
        strand,
        frame,
        attribute: attribute.to_string()
    };
    gtf_entry
}


pub(crate) fn parse_gtf(file_path:&str) -> Result<Vec<GtfEntry>, std::io::Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut gtfRecords =vec![];
    for line in reader.lines() {
        let gtf_entry=lineToGTF(line.unwrap().as_str());
        gtfRecords.push(gtf_entry)
    }
    Ok(gtfRecords)
}




pub(crate) fn parse_fastq(file: File) -> Result<Vec<structs::FastqRecord_isoncl_init>, Box<dyn Error>> {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    let mut records=vec![];
    let mut id_int=0;
    loop {
        let mut header = String::new();
        let header_read = reader.read_line(&mut header)?;
        if header_read == 0 {
            break;
        }

        header = header.trim().to_owned();
        header = (&header[1..]).to_string();
        header = shorten_header(&header).parse().unwrap();
        let mut sequence = String::new();
        let sequence_read = reader.read_line(&mut sequence)?;
        if sequence_read == 0 {
            break;
        }
        sequence = sequence.trim().to_owned();

        let mut quality_header = String::new();
        let quality_header_read = reader.read_line(&mut quality_header)?;
        if quality_header_read == 0 {
            break;
        }
        quality_header = quality_header.trim().to_owned();

        let mut quality = String::new();
        let quality_read = reader.read_line(&mut quality)?;
        if quality_read == 0 {
            break;
        }
        quality = quality.trim().to_owned();
        let score=0.0_f64;
        let error_rate=0.0_f64;
        let internal_id=id_int;
        records.push(FastqRecord_isoncl_init { header, sequence,/* quality_header,*/ quality, score,error_rate});
        id_int += 1;
    }
    Ok(records)
}

/*pub(crate) fn parse_fastq(file: File) -> (Vec<FastqRecord_isoncl_init>, HashMap<i32,String>) {
    let file = File::open(file)?;
    let mut id_map=HashMap::new();
    let reader = BufReader::new(file);
    let mut records=vec![];
    let mut line_cter=0;
    let mut fastq_record:FastqRecord;
    let mut header;
    let mut sequence;
    let mut quality_header;
    let mut quality;
    for line in reader.lines().advance_by() {
        if line_cter % 4==0{

            header= line.unwrap();

        }
        else if  line_cter %4 ==1{
            sequence =line.unwrap();
        }
        else if line_cter%4 ==2{
            quality_header=line.unwrap();
        }
        else {
            quality =line.unwrap();
        }
        println!("{}", line?);
        line_cter+=1;
    }

    return (records,id_map)
}
pub(crate) fn parse_fastq(file: File) -> (Vec<FastqRecord_isoncl_init>, HashMap<i32,String>) {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    let mut records=vec![];
    let mut id_int=0;
    let mut id_map=HashMap::new();
    loop {
        let mut header = String::new();
        let header_read = reader.read_line(&mut header);
        if header_read == 0 {
            break;
        }

        header = header.trim().to_owned();
        header = (&header[1..]).to_string();
        let mut sequence = String::new();
        let sequence_read = reader.read_line(&mut sequence);
        if sequence_read == 0 {
            break;
        }
        sequence = sequence.trim().to_owned();

        let mut quality_header = String::new();
        let quality_header_read = reader.read_line(&mut quality_header);
        if quality_header_read == 0 {
            break;
        }
        quality_header = quality_header.trim().to_owned();

        let mut quality = String::new();
        let quality_read = reader.read_line(&mut quality);
        if quality_read == 0 {
            break;
        }
        quality = quality.trim().to_owned();
        let score=0.0_f64;
        let error_rate=0.0_f64;
        let internal_id=id_int;
        id_map.insert(id_int, header);
        records.push(FastqRecord_isoncl_init { header,internal_id, sequence,/* quality_header,*/ quality, score,error_rate});
        id_int += 1;
    }
    return(records, id_map)
}*/


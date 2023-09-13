use std::fs::File;
use std::io::{BufReader, BufRead};
use std::error::Error;
use std::fmt;
use std::collections::HashMap;



pub(crate) struct FastaRecord {
    //a struct used to store fasta records
    pub header: String,
    pub sequence: String,
}


impl fmt::Display for FastaRecord {
    // enables displaying the fasta record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}


pub(crate) struct FastqRecord {
    //a struct used to store fastq records
    header: String,
    sequence: String,
    quality_header: String,
    quality: String,
}
impl fmt::Display for FastqRecord {
    // enables displaying the fastq record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}


pub struct FastqRecord_isoncl_init {
    //a struct used to store fastq records
    pub header: String,
    pub internal_id: i32,
    pub sequence: String,
    //pub(crate) quality_header: String,
    pub quality: String,
    pub score: f64,
    pub error_rate:f64,
}


impl FastqRecord_isoncl_init{
    pub fn get_header(&self)->&str{
        &self.header
    }
    pub fn get_int_id(&self)->&i32{
        &self.internal_id
    }
    pub fn get_sequence(&self)->&str{
        &self.sequence
    }
    pub fn get_quality(&self)->&str{
        &self.quality
    }
    pub fn get_score(&self)->&f64{
        &self.score
    }
    pub fn get_err_rate(&self)->&f64{
        &self.error_rate
    }
    pub fn set_error_rate(&mut self, new_error_rate: f64){
        self.error_rate = new_error_rate
    }
    pub fn set_score(&mut self, new_score: f64) {
        self.score = new_score
    }
}


impl fmt::Display for FastqRecord_isoncl_init {
    // enables displaying the fastq record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}", self.header, self.sequence)
    }
}


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


pub(crate) fn parse_fastq(file: File) -> Result<Vec<FastqRecord_isoncl_init>, Box<dyn Error>> {
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
        records.push(FastqRecord_isoncl_init { header,internal_id, sequence,/* quality_header,*/ quality, score,error_rate});
        id_int += 1;
    }
    Ok(records)
}


use std::fs::File;
use std::io::{BufReader, BufRead};
use std::error::Error;
use std::fmt;



pub(crate) struct FastaRecord {
    //a struct used to store fasta records
    header: String,
    sequence: String,
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



pub(crate) struct FastqRecord_isoncl_init {
    //a struct used to store fastq records
    pub(crate) header: String,
    pub(crate) sequence: String,
    //pub(crate) quality_header: String,
    pub(crate) quality: String,
    pub(crate) score: f64,
    pub(crate) error_rate:f64,
}

impl fmt::Display for FastqRecord {
    // enables displaying the fastq record
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}\n{}\n{}\n{}", self.header, self.sequence,self.quality_header,self.quality)
    }
}

pub(crate) fn parse_fasta(file: File) -> Result<Vec<FastaRecord>, Box<dyn Error>> {
    //Parses a fasta file, returns a vector of FastaRecords
    let reader = BufReader::new(file);
    let mut records = vec![];
    let mut lines = reader.lines();

    while let Some(line) = lines.next() {
        let line = line?;
        if line.starts_with(">") {
            let header = line[1..].to_owned();
            let mut sequence = String::new();
            while let Some(line) = lines.next() {
                let line = line?;
                if line.starts_with(">") {
                    break;
                }
                sequence.push_str(&line);
            }
            records.push(FastaRecord { header, sequence });
        }
    }

    Ok(records)
}


pub(crate) fn parse_fastq(file: File) -> Result<Vec<FastqRecord_isoncl_init>, Box<dyn Error>> {
    //Parses a fastq file, returns a vector of FastqRecords
    let mut reader = BufReader::new(file);
    let mut records=vec![];

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
        records.push(FastqRecord_isoncl_init { header, sequence,/* quality_header,*/ quality, score,error_rate});
    }
    Ok(records)
}


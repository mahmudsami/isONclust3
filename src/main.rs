use std::time::Duration;
mod file_actions;
mod generate_sorted_fastq_for_cluster;

fn main() {
    /*let args: Vec<String> = env::args().collect();
    let fasta_path = &args[1];
    let fastq_path = &args[2];
    let fasta_file = File::open(fasta_path).unwrap();
    let fasta_records = FileActions::parse_fasta(fasta_file).unwrap();
    let fastq_file = File::open(fastq_path).unwrap();
    let fastq_records = FileActions::parse_fastq(fastq_file).unwrap();
*/
    generate_sorted_fastq_for_cluster::sort_fastq_for_cluster(15,7.0,"/home/alexanderpetri/Rust/100k_sample.fastq")
    //let seq=dna!("ACTGACTACACAT");
    //let dna_sequence = "ATCGA";
    //let reversed_complement = generate_sorted_fastq_for_cluster::reverse_complement(dna_sequence);
    //println!("Reverse complement: {}", reversed_complement);

    /*for rec in fasta_records {
        println!("{}", rec);
    }
    for rec in fastq_records {
        println!("{}",rec);
    }*/






    /*
    THREAD test field: here I test how to use threads in Rust to actually make it possible to run several threads in parallel

     */
    /*let handle=thread::spawn(|| {
        for i in 1..10 {
            println!("hi number {} from the spawned thread!", i);
            thread::sleep(Duration::from_millis(1));
        }
    });
    handle.join().unwrap();
    for i in 1..5 {
        println!("hi number {} from the main thread!", i);
        thread::sleep(Duration::from_millis(1));
    }
    let threads: Vec<_> = (0..500)
        .map(|i| {
            thread::spawn(move || {
                println!("Thread #{} started!",i);
                thread::sleep(Duration::from_millis(5000));
                println!("Thread #{} finished!",i);
            })
        })
        .collect();

    for handle in threads {
        handle.join().unwrap();
    }*/

}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_reverse_complement() {
        let rev_comp = generate_sorted_fastq_for_cluster::reverse_complement("GGGGATCATCAGGGCTA");
        assert_eq!(rev_comp,"TAGCCCTGATGATCCCC");
        let rev_comp2 = generate_sorted_fastq_for_cluster::reverse_complement("ATCGA");
        assert_eq!(rev_comp,"TCGAT");
    }
    
}
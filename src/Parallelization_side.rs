

fn split_in_batches(n_threads: usize, data: &[u8])-> Vec<Vec<u8>>{
    let each_len= data.len() / n_threads + if data.len() % n_threads == 0 {0} else {1};
    let mut out = vec![Vec::with_capacity(each_len); n_threads];
    for(i,d) in data.iter().copied().enumerate(){
        out[i % n_threads].push(d);
    }
    out
}



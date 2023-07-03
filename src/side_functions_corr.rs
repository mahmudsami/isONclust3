


fn solve_WIS(intervals_sorted_by_finish: Vector<i32>)->Vector<i32>{
    let mut p = vec![];
    fill_p2(p,intervals_sorted_by_finish);
    let eps=0.0001f64;
    let mut v: Vec<Option<i32>> = vec![None];
    for (start, stop, w, _) in all_intervals_sorted_by_finish.iter() {
        let value = (w - 1) * (stop - start + epsilon);
        v.push(Some(value));
    }
    let opt= vec![];
    for j in 1..len(intervals_sorted_by_finish)+1{
        opt.append(std::cmp::max(v[j]+opt[p[j]],opt[j-1]));
    }
    let mut opt_indices = vec![];
    while j>=0{
        if j==0{
            break;
        }
        if v[j] + opt[p[j]] > opt[j-1]{
            opt_indicies.append(j - 1);
        }
        else{
            j-=1;
        }
    }
    Ok(opt_indices);
}


fn main(){
    let mut all_intervals =vec![];

}
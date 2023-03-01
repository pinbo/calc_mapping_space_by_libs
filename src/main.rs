// v0.1.0: 1st version
// to do: add lib names and add number to the 1st element of vector
// to make dynamic smaller app: cargo rustc --release -- -C prefer-dynamic
use std::io::{self, BufRead, BufReader};
use std::collections::HashMap;
use clap::{App, Arg};

fn main() {
    // parse options
    let matches = App::new("Calc_mapping_space")
        .version("0.1.0")
        .author("Junli Zhang <zhjl86@gmail.com>")
        .about("Calculate mapping space and %GC from MAPS output parsed_Kronos_mpileup.txt.gz")
        .arg(Arg::from_usage("-m, --max_missing=[NUMBER] 'maxium number of mising libs'"))
        .arg(Arg::from_usage("-c, --min_cov=[NUMBER]              'minimum coverage at positions'"))
        .arg(Arg::from_usage("-n, --max_depth=[NUMBER] 'maxium depth to summarize'"))
        .get_matches();
    // get the values
    let min_cov: u32 = matches.value_of("min_cov").unwrap_or("1").parse().expect("Please give a number of minimum coverage"); // minimum coverage at a position
    println!("minimum coverage is {}", min_cov);
    let max_missing: usize = matches.value_of("max_missing").unwrap_or("4").parse().expect("Please give a number of maximum missing libs");
    println!("Maximum missing libs allowed is {}", max_missing);
    let max_depth: u32 = matches.value_of("max_depth").unwrap_or("10").parse().expect("Please give a number of maximum missing libs");
    println!("Maximum depth to report is {}", max_depth);


    let input = io::stdin();
    let mut reader = BufReader::new(input.lock());
    let mut first_line = String::new();
    let _ = reader.read_line(&mut first_line);
    let ncol = first_line.split("\t").count();
    let nsample = (ncol - 3) / 4;
    println!("First line has {} columns and {} samples.", ncol, nsample);
    let zero_vec = vec![0; nsample];
    let min_lib_count = nsample - max_missing;
    let mut good_lines = 0;
    let mut map = HashMap::new(); // count A, T, G, C
    let mut map_libs = HashMap::new(); // key is coverage, value is a vector of occurance for each lib
    let mut map_gc = HashMap::new();
    let nts = "ATGCatgc";
    let gc = "GCgc";
    // get the lib names
    let mut lib_names = vec![""; nsample];
    let mut n = 0;
    let mut target = 4;
    for word in first_line.split("\t") {
        n += 1;
        if n == target {
            lib_names[(target - 4) / 4] = word; //.to_string();
            target += 4;
        }
    }

    for line in reader.lines() {
        let mut ngood = 0;
        let mut n = 0;
        let mut target = 7;
        let mut nt = String::from("A"); // initial value for column 2
        let mut sample = 0; // nth sample
        for word in line.unwrap().split("\t") {
            n += 1;
            if n == 3 {
                nt = word.to_string();
            }
            if n == target {
                target += 4;
                if word != "." {
                    let nn: u32 = word.parse().unwrap();
                    let kk = if nn < max_depth {nn} else {max_depth};
                    let count = map_libs.entry(kk).or_insert(zero_vec.clone());
                    count[sample] += 1;
                    let count2 = map_gc.entry(kk).or_insert(zero_vec.clone());
                    if gc.contains(&nt) {count2[sample] += 1;}
                    if nn >= min_cov {
                        ngood += 1;
                    }
                }
                sample += 1; // nth sample
            } else {
                continue;
            }
        }
        if ngood >= min_lib_count {
            good_lines += 1;
            let count = map.entry(nt).or_insert(0);
            *count += 1;
        }
    }
    println!("The input has {} good lines.", good_lines);
    let mut total = 0;
    let mut gc_count = 0;
    for (key, value) in &map {
        println!("{}\t{}", &key, value);
        if nts.contains(key) {total += value;}
        if gc.contains(key) {gc_count += value;}
    }
    println!("GC% is {:.2}%", (gc_count as f64 / total as f64) * 100.0);

    // organize individual counts
    println!("Cov\t{:?}", lib_names);

    println!("\nCoverage counts are below:");
    for (key, value) in &map_libs {
        println!("{}\t{:?}", key, value);
    }
    println!("\nGC counts are below");
    for (key, value) in &map_gc {
        println!("{}\t{:?}", key, value);
    }

}



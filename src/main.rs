// v0.1.0: 1st version
// v0.2.0: add more options same as beta-maps1-v2.py
// to do: add lib names and add number to the 1st element of vector
// to make dynamic smaller app: cargo rustc --release -- -C prefer-dynamic
use std::io::{self, BufRead, BufReader};
use std::collections::HashMap;
use clap::{App, Arg};

fn main() {
    // parse options
    let matches = App::new("calc_mapping_space_by_libs")
        .version("0.2.0")
        .author("Junli Zhang <zhjl86@gmail.com>")
        .about("Calculate mapping space and %GC from MAPS output parsed_Kronos_mpileup.txt.gz")
        .arg(Arg::from_usage("-c, --min_cov=[NUMBER]              'minimum coverage at positions (1)'"))
        .arg(Arg::from_usage("-n, --max_depth=[NUMBER] 'maxium depth to summarize in the output (10)'"))
        .arg(Arg::from_usage("-l, --min_libs=[NUMBER] 'minimum number of libraries with at least 1 coverage to be considered a valid position (15)'"))
        .arg(Arg::from_usage("-C, --max_cov=[NUMBER] 'maximum position coverage (sum of all libs) to be considered as a vaild position (10000)'"))
        .get_matches();
    // get the values
    let min_cov: u32 = matches.value_of("min_cov").unwrap_or("1").parse().expect("Please give a number of minimum coverage"); // minimum coverage at a position
    println!("minimum coverage is {}", min_cov);
    let max_depth: u32 = matches.value_of("max_depth").unwrap_or("10").parse().expect("Please give a number of maximum depth to summarize");
    println!("Maximum depth to report is {}", max_depth);
    let min_libs: usize = matches.value_of("min_libs").unwrap_or("15").parse().expect("Please give a number of minimum libs to consider");
    println!("Minimum libraries covered at least once is {}", min_libs);
    let max_cov: u32 = matches.value_of("max_cov").unwrap_or("10000").parse().expect("Please give a number of maximum coverage for all libs for a valid position");
    println!("Maximum coverage for this position in all libs is {}", max_cov);


    let input = io::stdin();
    let mut reader = BufReader::new(input.lock());
    let mut first_line = String::new();
    let _ = reader.read_line(&mut first_line);
    let ncol = first_line.split("\t").count();
    let nsample = (ncol - 3) / 4;
    println!("First line has {} columns and {} samples.", ncol, nsample);
    let zero_vec = vec![0; nsample];
    let max_missing = nsample - min_libs;
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
    // tmp store for each line
    let mut vec0 = vec![0; nsample];
    let mut gc0 = vec![0; nsample];

    for line in reader.lines() {
        let mut nmiss = 0;
        let mut ngood = 0;
        let mut n = 0;
        let mut target = 7;
        let mut nt = String::from("A"); // initial value for column 2
        let mut sample = 0; // nth sample
        let mut pos_coverage = 0; // calculate the total coverage for this position (line)
        for word in line.unwrap().split("\t") {
            n += 1;
            if n == 3 {
                nt = word.to_string();
            }
            if n == target {
                target += 4;
                if word == "." {
                    nmiss += 1;
                    if nmiss > max_missing {break}
                    vec0[sample] = 0;
                    gc0[sample] = 0;
                } else {
                    let nn: u32 = word.parse().unwrap();
                    pos_coverage += nn;
                    if pos_coverage > max_cov {break}
                    // tmp store all the numbers
                    vec0[sample] = nn;
                    gc0[sample] = if gc.contains(&nt) {1} else {0};
                    if nn >= min_cov {
                        ngood += 1;
                    }
                }
                sample += 1; // nth sample
            } else {
                continue;
            }
        }
        // check whether this position is valid
        if nmiss <= max_missing && pos_coverage <= max_cov {
            for i in 0..nsample {
                let dep = vec0[i];
                if dep > 0 {
                    let kk = if dep < max_depth {dep} else {max_depth};
                    let count = map_libs.entry(kk).or_insert(zero_vec.clone());
                    count[i] += 1;
                    let count2 = map_gc.entry(kk).or_insert(zero_vec.clone());
                    count2[i] += gc0[i];
                }
            }
        }

        // get good
        if ngood >= min_libs {
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



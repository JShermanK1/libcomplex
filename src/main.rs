
use clap::{Arg};
use polars::prelude::*;
use rayon::current_num_threads;
use std::fs::File;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::io::{BufReader, 
              BufWriter,
              Write,
              BufRead};

fn cli() -> clap::Command {
clap::Command::new("downsample")
            .args(&[
                Arg::new("input")
                    .long("input")
                    .short('i')
                    .required(true)
                    .value_name("FILE"),
                Arg::new("output")
                    .long("output")
                    .short('o')
                    .required(true)
                    .value_name("FILE"),
            ])
}

fn main() -> Result<()> {

    println!("starting with {} threads", current_num_threads());    

    let matches = cli().get_matches();

    let input = matches.get_one::<String>("input")
                        .unwrap();
    let output = matches.get_one::<String>("output")
                        .unwrap();

    let mut no_m;
    {

        let mut bed = CsvReader::from_path(PathBuf::from(input))?
                            .infer_schema(Some(5))
                            .with_projection(Some(vec![0, 1, 3, 5, 8, 9]))
                            .with_delimiter('\t' as u8)
                            .has_header(false)
                            .finish()?;
        println!("{}", bed);
        bed.set_column_names(&["chrom1", "start1", "chrom2", "end2", "strand1", "strand2"])?;

        no_m = bed.filter(&bed["chrom1"].not_equal("chrM"))?
                    .sort(&["chrom1", "start1"], vec![false, false])?;
    }
    let total = no_m.height() as f64;

    no_m.with_row_count_mut("index", Some(1));

    let distinct_frame = no_m.distinct_stable(Some(&["chrom1".to_string(), 
                                                     "start1".to_string(),
                                                     "chrom2".to_string(),
                                                     "end2".to_string(),
                                                     "strand1".to_string(),
                                                     "strand2".to_string()]), DistinctKeepStrategy::Last)?;
    let distinct = distinct_frame.height() as f64;

    let idx = distinct_frame.column("index")?;
    let shift = idx.shift(-1);
    println!("idx:{}", idx.head(Some(10)));
    println!("shift:{}", shift.head(Some(10)));
    let mut counts = &shift - idx;
    counts = counts.shift(1);
    counts = counts.u32()?
                    .fill_null_with_values(idx.u32()?.get(0).unwrap())?
                    .into_series();
    let counts_table = counts.value_counts()?;
    println!("{}", counts_table);
    
    let unique = counts_table["counts"].u32()?.get(0).unwrap() as f64;
    let pairs = counts_table["counts"].u32()?.get(1).unwrap_or(0) as f64;

    
    let nrf = distinct / total;
    let pbc1 = unique / distinct;
    let pbc2 = unique / pairs;
    
    let mut output = File::create(PathBuf::from(output)).expect("Error creating outputfile");

    output.write(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}", total, distinct, unique, pairs, nrf, pbc1, pbc2).as_bytes())?;
                                
    Ok(())
}
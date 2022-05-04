extern crate clap;
extern crate env_logger;
extern crate exitcode;

use clap::{crate_version, Arg, Command};
use log::{error, info};
use mimalloc::MiMalloc;
use std::fs::File;
use std::io;

use msbwt2::bwt_converter::{convert_to_vec, save_bwt_numpy};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    //initialize logging for our benefit later
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    //this is the CLI block, params that get populated appear before
    let mut in_fn: String = "stdin".to_string();

    let matches = Command::new("msbwt2 BWT Converter")
        .version(crate_version!())
        .author("J. Matthew Holt <jholt@hudsonalpha.org>")
        .about("msbwt2 BWT Converter - this will convert an external BWT to our expected representation")
        .arg(Arg::new("in_fn")
            .short('i')
            .long("--input")
            .takes_value(true)
            .help("The raw uncompressed BWT (default: stdin)"))
        .arg(Arg::new("COMP_MSBWT.NPY")
            .help("The location to store the compressed BWT")
            .required(true)
            .index(1))
        .get_matches();

    //pull out required values
    let bwt_fn: String = matches.value_of("COMP_MSBWT.NPY").unwrap().to_string();

    //optional values
    in_fn = matches.value_of_t("in_fn").unwrap_or(in_fn);

    info!("Input parameters (required):");
    info!("\tInput BWT: \"{}\"", in_fn);
    let input_reader: Box<dyn io::Read> = if &in_fn == "stdin" {
        Box::new(io::stdin())
    } else {
        Box::new(match File::open(&in_fn) {
            Ok(fp) => fp,
            Err(e) => {
                error!("Failed to open BWT file: {:?}", e);
                std::process::exit(exitcode::NOINPUT);
            }
        })
    };

    info!("\tOutput BWT: \"{}\"", bwt_fn);
    match File::create(&bwt_fn) {
        Ok(file) => file,
        Err(e) => {
            error!("Failed to create output BWT file: {:?}", e);
            std::process::exit(exitcode::NOINPUT);
        }
    };

    //this is where the work happens
    let comp_bwt = convert_to_vec(input_reader);
    match save_bwt_numpy(&comp_bwt[..], &bwt_fn) {
        Ok(_) => {}
        Err(e) => {
            error!("Error saving BWT to file: {:?}", bwt_fn);
            error!("Error: {:?}", e);
            std::process::exit(exitcode::IOERR);
        }
    };

    info!("RLE-BWT conversion complete.");
}
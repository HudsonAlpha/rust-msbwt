
extern crate clap;
extern crate env_logger;
extern crate exitcode;
extern crate log;

use clap::{App, Arg, value_t, values_t};
use log::{error, info};
use mimalloc::MiMalloc;
use std::fs::File;

use msbwt2::bwt_converter::save_bwt_runs_numpy;
use msbwt2::dynamic_bwt::{create_from_fastx,DynamicBWT};
use msbwt2::string_util::INT_TO_STRING;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");

fn main() {
    //initialize logging for our benefit later
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let matches = App::new("msbwt2 BWT Builder")
        .version(VERSION.unwrap_or("?"))
        .author("J. Matthew Holt <jholt@hudsonalpha.org>")
        .about("msbwt2 BWT Builder - will construct a BWT from one or more FASTX files")
        .arg(
            Arg::with_name("out_bwt")
            .short("o")
            .long("--out-bwt")
            .takes_value(true)
            .help("The output BWT (default: stdout)")
        )
        .arg(
            Arg::with_name("FASTX")
                .help("The FASTQ/A file(s) to load into the BWT, gzip accepted")
                .required(true)
                .multiple(true)
                .index(1)
        )
        .get_matches();

    let fastx_fns: Vec<String> = values_t!(matches.values_of("FASTX"), String).unwrap_or_else(|_| vec![]);
    let out_fn: String = value_t!(matches.value_of("out_bwt"), String).unwrap_or_else(|_| "stdout".to_string());
    //TODO: make this a command line option
    let unsorted_strings: bool = value_t!(matches.value_of("unsorted"), bool).unwrap_or_else(|_| false);
    let sorted_strings: bool = !unsorted_strings;
    
    info!("Input parameters (required):");
    info!("\tFASTX: {:?}", fastx_fns);
    info!("\tout_bwt: {:?}", out_fn);
    info!("Optional Parameters:");
    info!("\tsort order: {}", match sorted_strings {
        true => "lexicographical",
        false => "chronological"
    });

    //check inputs
    for fastx_fn in fastx_fns.iter() {
        match File::open(fastx_fn) {
            Ok(_) => {}
            Err(e) => {
                error!("Failed to open FASTX file: {:?}", fastx_fn);
                error!("Error: {:?}", e);
                std::process::exit(exitcode::NOINPUT);
            }
        };
    }

    //check outputs
    if out_fn != "stdout" {
        match File::create(&out_fn) {
            Ok(file) => file,
            Err(e) => {
                error!("Failed to create output BWT file: {:?}", out_fn);
                error!("Error: {:?}", e);
                std::process::exit(exitcode::CANTCREAT);
            }
        };
    }

    //create the BWT
    let bwt: DynamicBWT = match create_from_fastx(&fastx_fns, sorted_strings) {
        Ok(result) => result,
        Err(e) => {
            error!("Error while parsing FASTX files: {:?}", fastx_fns);
            error!("Error: {:?}", e);
            std::process::exit(exitcode::IOERR);
        }
    };
    
    //this is what we should put in the file
    if out_fn == "stdout" {
        for sym in bwt.iter() {
            print!("{}", INT_TO_STRING[sym as usize] as char);
        }
        println!();
    }
    else {
        info!("Saving results to file: {:?}", out_fn);
        match save_bwt_runs_numpy(bwt.run_iter(), &out_fn) {
            Ok(_) => {},
            Err(e) => {
                error!("Error saving BWT to file: {:?}", out_fn);
                error!("Error: {:?}", e);
                std::process::exit(exitcode::IOERR);
            }
        };
    }

    info!("Processes successfully finished.")
}
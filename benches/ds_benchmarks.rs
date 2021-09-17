
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use msbwt2::bwt_converter::save_bwt_runs_numpy;
use msbwt2::dynamic_bwt::DynamicBWT;
use msbwt2::msbwt_core::BWT;
use msbwt2::rle_bplus_tree::RLEBPlusTree;
use msbwt2::rle_bwt::RleBWT;
use msbwt2::run_block_av_flat::RLEBlock;
use msbwt2::string_util;
use msbwt2::wavelet_tree::WaveletTree;
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

fn get_random_insert(length: usize) -> (Vec<u64>, Vec<u8>) {
    //this is how to provide a constant "random" set of inserts to play with
    let mut rng = StdRng::seed_from_u64(0);

    //get random symbols AND positions into the array
    let mut inserted: Vec<u8> = vec![];
    let mut positions: Vec<u64> = vec![];
    for i in 0..length {
        let symbol: u8 = rng.gen_range(0, 6);
        let position: u64 = rng.gen_range(0, i+1) as u64;
        inserted.push(symbol);
        positions.push(position);    
    }
    (positions, inserted)
}

fn get_biased_insert(length: usize) -> (Vec<u64>, Vec<u8>) {
    //this is how to provide a constant "random" set of inserts to play with
    let mut rng = StdRng::seed_from_u64(0);

    //get random symbols AND positions into the array
    let mut inserted: Vec<u8> = vec![];
    let mut positions: Vec<u64> = vec![];
    for i in 0..length {
        let symbol: u8 = 0;//rng.gen_range(0, 6);
        let position: u64 = rng.gen_range(0, i+1) as u64;
        inserted.push(symbol);
        positions.push(position);    
    }
    (positions, inserted)
}

// should have 4 length blocks + 3 single-character runs
fn get_constant_insert() -> (Vec<u64>, Vec<u8>) {
    let mut inserted: Vec<u8> = vec![0; 1000];
    let mut positions: Vec<u64> = vec![0; 1000];

    inserted.push(1);
    positions.push(500);
    inserted.push(1);
    positions.push(250);
    inserted.push(1);
    positions.push(750);

    (positions, inserted)
}

fn create_fixed_bwt() -> String {
    let out_fn: String = "./test_data/example_output_bench_001.npy".to_string();
    let mut bwt: DynamicBWT = Default::default();
    for _ in 0..10000 {
        bwt.insert_string("AAAAACCCCCTTTTTGGGGGACGTACGTTGCATGCA", true);
    }
    save_bwt_runs_numpy(bwt.run_iter(), &out_fn).unwrap();
    out_fn
}

pub fn bench_rle_bplus_tree(c: &mut Criterion) {
    let (positions, symbols) = get_random_insert(100000);
    
    c.bench_function("rle_bplus_tree_100k_random", |b| b.iter(|| {
        let mut tree: RLEBPlusTree = Default::default();
        for (position, symbol) in positions.iter().zip(symbols.iter()) {
            black_box(tree.insert_and_count(*position, *symbol));
        }
    }));

    c.bench_function("wavelet_tree_100k_random", |b| b.iter(|| {
        let mut tree: WaveletTree = Default::default();
        for (position, symbol) in positions.iter().zip(symbols.iter()) {
            black_box(tree.insert_and_count(*position, *symbol));
        }
    }));

    let (positions, symbols) = get_biased_insert(1000000);
    c.bench_function("rle_bplus_tree_1M_biased", |b| b.iter(|| {
        let mut tree: RLEBPlusTree = Default::default();
        for (position, symbol) in positions.iter().zip(symbols.iter()) {
            black_box(tree.insert_and_count(*position, *symbol));
        }
    }));

    c.bench_function("wavelet_tree_1M_biased", |b| b.iter(|| {
        let mut tree: WaveletTree = Default::default();
        for (position, symbol) in positions.iter().zip(symbols.iter()) {
            black_box(tree.insert_and_count(*position, *symbol));
        }
    }));
}

pub fn bench_rle_block(c: &mut Criterion) {
    let (positions, symbols) = get_constant_insert();
    c.bench_function("run_block_av_flat_insert_and_split", |b| b.iter(|| {
        let mut block: RLEBlock = Default::default();
        for (position, symbol) in positions.iter().zip(symbols.iter()) {
            black_box(block.insert_and_count(*position, *symbol));
        }
        
        for i in 0..1000 {
            black_box(block.count(i, 0));
        }
        black_box(block.split())
    }));
}

pub fn bench_bwt_queries(c: &mut Criterion) {
    let bwt_fn: String = create_fixed_bwt();
    let mut rle_bwt: RleBWT = Default::default();
    rle_bwt.load_numpy_file(&bwt_fn).unwrap();
    let mut dyn_bwt: DynamicBWT = Default::default();
    dyn_bwt.load_numpy_file(&bwt_fn).unwrap();

    let query1 = string_util::convert_stoi(&"ACGT");
    let query2 = string_util::convert_stoi(&"AACC");

    c.bench_function("rle_bwt_count_kmer", |b| b.iter(|| {
        black_box(rle_bwt.count_kmer(&query1));
        black_box(rle_bwt.count_kmer(&query2));
    }));

    c.bench_function("dyn_bwt_count_kmer", |b| b.iter(|| {
        black_box(dyn_bwt.count_kmer(&query1));
        black_box(dyn_bwt.count_kmer(&query2));
    }));
}

criterion_group!(benches, bench_rle_bplus_tree, bench_rle_block, bench_bwt_queries);
criterion_main!(benches);
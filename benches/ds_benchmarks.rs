
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use msbwt2::rle_bplus_tree::RLEBPlusTree;
use msbwt2::run_block_av_flat::RLEBlock;
use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;
//use rand_core::SeedableRng;

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

pub fn bench_rle_bplus_tree(c: &mut Criterion) {
    let (positions, symbols) = get_random_insert(10000);
    
    c.bench_function("rle_bplus_tree_10k_random", |b| b.iter(|| {
        let mut tree: RLEBPlusTree = Default::default();
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

criterion_group!(benches, bench_rle_bplus_tree, bench_rle_block);
criterion_main!(benches);
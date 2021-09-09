[![Crates.io](https://img.shields.io/crates/d/msbwt2.svg)](https://crates.io/crates/msbwt2)
[![Crates.io](https://img.shields.io/crates/v/msbwt2.svg)](https://crates.io/crates/msbwt2)
[![Crates.io](https://img.shields.io/crates/l/msbwt2.svg)](https://crates.io/crates/msbwt2)
[![Build status](https://github.com/HudsonAlpha/rust-msbwt/actions/workflows/quickstart-ci.yml/badge.svg)](https://github.com/HudsonAlpha/rust-msbwt/actions)

# msbwt2
The intent of crate is to provide Rust functionality for querying a Multi-String BWT (MSBWT), and is mostly based on the same methodology used by the original [msbwt](https://github.com/holtjma/msbwt).

NOTE: This is very much a work-in-progress and currently only being updated as a side project during spare time.
If you have any feature requests, feel free to submit a new issue on GitHub.
Here is a current list of planned additions:

1. Incorporate the high-memory BWT implementation from `fmlrc2`
2. Add a built-in BWT construction tool (instead of relying on external tools like `ropebwt2`)
3. Add some more query functionality

## Installation
All installation options assume you have installed [Rust](https://www.rust-lang.org) along with the `cargo` crate manager for Rust.

### From Cargo
```bash
cargo install msbwt2
msbwt2-convert -h
```

### From GitHub
```bash 
git clone https://github.com/HudsonAlpha/rust-msbwt.git
cd rust-msbwt
#testing optional
cargo test --release
cargo build --release
./target/release/msbwt2-convert -h
```

## Usage
### MSBWT Building
The Multi-String Burrows Wheeler Transform (MSBWT or BWT) must be built prior to performing any queries. 
Currently, there is no built in builder, but we expect to have one included soon.
For now, the [original instructions](https://github.com/holtjma/fmlrc#building-the-short-read-bwt) can be used.

Given a FASTQ file of reads (`reads.fq.gz`), you can also use the following command from this crate to create a BWT at `comp_msbwt.npy`.  
Note that this command requires the [ropebwt2](https://github.com/lh3/ropebwt2) executable to be installed:
```
gunzip -c reads.fq.gz | \
    awk 'NR % 4 == 2' | \
    sort | \
    tr NT TN | \
    ropebwt2 -LR | \
    tr NT TN | \
    msbwt2-convert comp_msbwt.npy
```

#### Optional Construction Speedup
If you are **only** using the BWT for k-mer queries, then the `sort` can be removed from the above command. 
This will reduce construction time significantly, but loses the read recovery property of the BWT.

### Queries
The general use case of the library is k-mer queries, which can be performed as follows:
```rust
use msbwt2::msbwt_core::BWT;
use msbwt2::rle_bwt::RleBWT;
use msbwt2::string_util;
let mut bwt = RleBWT::new();
let filename: String = "test_data/two_string.npy".to_string();
bwt.load_numpy_file(&filename);
assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACGT")), 1);
```

## Reference
`msbwt2` does not currently have a pre-print or paper. If you use `msbwt2`, please cite the one of the `msbwt` papers:

[Holt, James, and Leonard McMillan. "Merging of multi-string BWTs with applications." Bioinformatics 30.24 (2014): 3524-3531.](https://doi.org/10.1093/bioinformatics/btu584)

[Holt, James, and Leonard McMillan. "Constructing Burrows-Wheeler transforms of large string collections via merging." Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics. 2014.](https://doi.org/10.1145/2649387.2649431)

## License
Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution
Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
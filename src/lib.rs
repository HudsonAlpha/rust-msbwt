/*!
# MSBWT v2
This library provides access to a rust-based implementation of Multi-String BWT (MSBWT) queries.
Currently, the BWT is assumed to have been generated externally (typically with a tool like ropebwt2) and stored in the same numpy format as expected by the original `msbwt` tool.

## Example
```
use msbwt2::msbwt_core::BWT;
use msbwt2::rle_bwt::RleBWT;
use msbwt2::string_util;
let mut bwt = RleBWT::new();
let filename: String = "test_data/two_string.npy".to_string();
bwt.load_numpy_file(&filename);
assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACGT")), 1);
```
*/

/// Contains the function for reformating a BWT string into the expected run-length format or numpy file
pub mod bwt_converter;
/// Contains helper functions related to BWT construction, primarily for testing purposes
pub mod bwt_util;
/// Includes the trait for a multi-string BWT
pub mod msbwt_core;
/// This is the classic RLE implementation from the original msbwt package
pub mod rle_bwt;
/// Contains inline functions for converting between strings and integer formats
pub mod string_util;

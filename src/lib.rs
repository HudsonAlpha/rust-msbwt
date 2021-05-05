
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

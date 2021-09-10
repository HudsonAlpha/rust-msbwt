
extern crate log;

use log::info;
use std::convert::TryInto;
use std::io::prelude::*;
use std::fs;

use crate::msbwt_core::*;

use crate::run_block_av_flat::VC_LEN;
use crate::rle_bplus_tree::RLEBPlusTree;
use crate::string_util::convert_stoi;

const INITIAL_QUERY: usize = 10;
const COST_FACTOR: f64 = 0.000001;

pub struct DynamicBWT {
    tree_bwt: RLEBPlusTree,
    symbol_counts: [u64; VC_LEN],
    start_index: [u64; VC_LEN],
    string_count: u64,
    total_count: u64,
    sort_query_len: f64,
    short_circuits: [usize; 3],
}

impl Default for DynamicBWT {
    fn default() -> Self {
        DynamicBWT {
            string_count: 0,
            total_count: 0,
            tree_bwt: Default::default(),
            symbol_counts: [0; VC_LEN],
            start_index: [0; VC_LEN],
            sort_query_len: INITIAL_QUERY as f64,
            short_circuits: [0; 3]
        }
    }
}

impl BWT for DynamicBWT {
    /// Initializes the BWT from a compressed BWT vector.
    /// # Arguments
    /// * `bwt` - the run-length encoded BWT stored in a Vec<u8> 
    /// # Examples
    /// ```rust
    /// use msbwt2::msbwt_core::BWT;
    /// use msbwt2::dynamic_bwt::DynamicBWT;
    /// use msbwt2::bwt_converter::convert_to_vec;
    /// //strings "ACGT" and "CCGG"
    /// let seq = "TG$$CAGCCG";
    /// let vec = convert_to_vec(seq.as_bytes());
    /// let mut bwt = DynamicBWT::new();
    /// bwt.load_vector(vec);
    /// ```
    fn load_vector(&mut self, bwt: Vec<u8>) {
        info!("Initializing BWT with {:?} compressed values...", bwt.len());

        //reset all of these 
        self.tree_bwt = Default::default();
        self.symbol_counts = [0; VC_LEN];
        self.sort_query_len = INITIAL_QUERY as f64;
        self.short_circuits = [0; 3];

        //general strategy here is to build up the count of a character and then insert them
        let mut prev_char: u8 = 255;
        let mut current_char: u8;
        let mut power_multiple: u64 = 1;
        let mut current_count: u64;

        //go through each compressed block in the RLE encoded vector to calculate total character counts
        let mut current_index: u64 = 0;
        for value in bwt {
            current_char = value & MASK;
            if current_char == prev_char {
                power_multiple *= NUM_POWER as u64;
            }
            else {
                power_multiple = 1;
            }
            prev_char = current_char;
            current_count = (value >> LETTER_BITS) as u64 * power_multiple;

            for _ in 0..current_count {
                let _dummy = self.tree_bwt.insert_and_count(current_index, current_char);
                current_index += 1;
            }

            //add this to the total counts
            self.symbol_counts[current_char as usize] += current_count;
        }

        self.string_count = self.symbol_counts[0];
        self.total_count = self.symbol_counts.iter().sum();
        self.start_index = self.symbol_counts.iter().scan(0, |sum, &count| {
            *sum += count;
            Some(*sum - count)
        }).collect::<Vec<u64>>().try_into().unwrap();
        info!("Loaded BWT with symbol counts: {:?}", self.symbol_counts);
        info!("Finished BWT initialization.")
    }

    /// Initializes the BWT from the numpy file format for compressed BWTs
    /// # Arguments
    /// * `filename` - the name of the file to load into memory
    /// # Examples
    /// ```rust
    /// use msbwt2::msbwt_core::BWT;
    /// use msbwt2::dynamic_bwt::DynamicBWT;
    /// use msbwt2::string_util;
    /// let mut bwt = DynamicBWT::new();
    /// let filename: String = "test_data/two_string.npy".to_string();
    /// bwt.load_numpy_file(&filename);
    /// assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACGT")), 1);
    /// ```
    fn load_numpy_file(&mut self, filename: &str) -> std::io::Result<()> {
        //read the numpy header: http://docs.scipy.org/doc/numpy-1.10.1/neps/npy-format.html
        //get the initial file size
        let file_metadata: fs::Metadata = fs::metadata(&filename)?;
        let full_file_size: u64 = file_metadata.len();

        //read the initial fixed header
        let mut file = fs::File::open(&filename)?;
        let mut init_header: Vec<u8> = vec![0; 10];
        let read_count: usize = file.read(&mut init_header[..])?;
        if read_count != 10 {
            panic!("Could not read initial 10 bytes of header for file {:?}", filename);
        }

        //read the dynamic header
        let header_len: usize = init_header[8] as usize + 256 * init_header[9] as usize;
        let mut skip_bytes: usize = 10+header_len;
        if skip_bytes % 16 != 0 {
            skip_bytes = ((skip_bytes / 16)+1)*16;
        }
        let mut skip_header: Vec<u8> = vec![0; skip_bytes-10];
        match file.read_exact(&mut skip_header[..]) {
            Ok(()) => {},
            Err(e) => {
                return Err(
                    std::io::Error::new(
                        e.kind(),
                        format!("Could not read bytes 10-{:?} of header for file {:?}, root-error {:?}", skip_bytes, filename, e)
                    )
                );
            }
        }
        
        //parse the header string for the expected length, requires a lot of manipulation of the string because of numpy header styling
        let header_string = String::from_utf8(skip_header).unwrap()
            .replace("\'", "\"")
            .replace("False", "false")
            .replace("(", "[")
            .replace(")", "]")
            .replace(", }", "}")
            .replace(", ]", "]")
            .replace(",]", "]");
        let header_dict: serde_json::Value = serde_json::from_str(&header_string)
            .unwrap_or_else(|_| panic!("Error while parsing header string: {:?}", header_string));
        let expected_length: u64 = header_dict["shape"][0].as_u64().unwrap();
        
        //check that the disk size matches our expectation
        let bwt_disk_size: u64 = full_file_size - skip_bytes as u64;
        if expected_length != bwt_disk_size {
            return Err(
                std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!("Header indicates shape of {:?}, but remaining file size is {:?}", expected_length, bwt_disk_size)
                )
            );
        }

        //finally read in everything else
        info!("Loading BWT with {:?} compressed values from disk...", bwt_disk_size);
        let mut bwt_data: Vec<u8> = Vec::<u8>::with_capacity(bwt_disk_size as usize);
        let read_count: usize = file.read_to_end(&mut bwt_data)?;
        if read_count as u64 != bwt_disk_size {
            return Err(
                std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!("Only read {:?} of {:?} bytes of BWT body for file {:?}", read_count, bwt_disk_size, filename)
                )
            );
        }
        
        //we loaded the file into memory, now just do the load from vec
        self.load_vector(bwt_data);

        Ok(())
    }

    /// Returns the total number of occurences of a given symbol
    /// # Arguments
    /// * `symbol` - the symbol in integer form
    /// # Examples
    /// ```rust
    /// # use msbwt2::msbwt_core::BWT;
    /// # use msbwt2::dynamic_bwt::DynamicBWT;
    /// # use msbwt2::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let vec = convert_to_vec(seq.as_bytes());
    /// # let mut bwt = DynamicBWT::new();
    /// # bwt.load_vector(vec);
    /// let string_count = bwt.get_symbol_count(0);
    /// assert_eq!(string_count, 2);
    /// ```
    #[inline]
    fn get_symbol_count(&self, symbol: u8) -> u64 {
        self.symbol_counts[symbol as usize] as u64
    }

    /// This will return the total number of symbols contained by the BWT
    /// # Examples
    /// ```rust
    /// # use msbwt2::msbwt_core::BWT;
    /// # use msbwt2::dynamic_bwt::DynamicBWT;
    /// # use msbwt2::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let vec = convert_to_vec(seq.as_bytes());
    /// # let mut bwt = DynamicBWT::new();
    /// # bwt.load_vector(vec);
    /// let total_size = bwt.get_total_size();
    /// assert_eq!(total_size, 10);
    /// ```
    #[inline]
    fn get_total_size(&self) -> u64 {
        self.total_count as u64
    }

    /// Performs a range constraint on a BWT range. This implicitly represents prepending a character `sym` to a k-mer
    /// represented by `input_range` to create a new range representing a (k+1)-mer.
    /// # Arguments
    /// * `sym` - the symbol to pre-pend in integer form
    /// * `input_range` - the range to pre-pend to
    /// # Safety
    /// This function is unsafe because there are no guarantees that the symbol or bounds will be checked by the implementing structure.
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange {
        BWTRange {
            l: (self.start_index[sym as usize]+self.tree_bwt.count(input_range.l, sym)),
            h: (self.start_index[sym as usize]+self.tree_bwt.count(input_range.h, sym))
        }
    }
}

impl DynamicBWT {
    /// Allocation function for the BWT, look at `load_vector(...)` for initialization.
    /// # Examples
    /// ```rust
    /// use msbwt2::dynamic_bwt::DynamicBWT;
    /// let mut bwt = DynamicBWT::new();
    /// ```
    pub fn new() -> Self {
        Default::default()
    }

    #[inline]
    pub fn get_symbol_counts(&self) -> [u64; VC_LEN] {
        self.symbol_counts
    }

    #[inline]
    pub fn get_height(&self) -> usize {
        self.tree_bwt.get_height()
    }

    #[inline]
    pub fn get_node_count(&self) -> usize {
        self.tree_bwt.get_node_count()
    }
    
    /// This is the main function for adding strings to the DynamicBWT.
    /// # Arguments
    /// * val - the string to add
    /// * sorted - if true, this will add the string to it's sorted position, otherwise the end
    /// # Examples
    /// ```rust
    /// use msbwt2::dynamic_bwt::DynamicBWT;
    /// let data: String = "ACGNT".to_string();
    /// let bwt: Vec<u8> = vec![5, 0, 1, 2, 3, 4];
    /// let mut ubwt: DynamicBWT = Default::default();
    /// ubwt.insert_string(&data, false);
    /// assert_eq!(ubwt.to_vec(), bwt);
    /// ```
    #[inline]
    pub fn insert_string(&mut self, val: &str, sorted: bool) {
        let int_form: Vec<u8> = convert_stoi(val);

        //initial position is the total number of string
        let mut next_insert: u64;
        
        if sorted {
            let mut start_index: u64 = 0;
            next_insert = self.total_count;

            //attempt a short circuit
            let query_len = std::cmp::min(self.sort_query_len as usize, int_form.len());
            for pred_symbol in int_form[..query_len].iter().rev() {
                start_index = self.tree_bwt.count(start_index, *pred_symbol)+self.start_index[*pred_symbol as usize];
                next_insert = self.tree_bwt.count(next_insert, *pred_symbol)+self.start_index[*pred_symbol as usize];
            }
            start_index = self.tree_bwt.count(start_index, 0);
            next_insert = self.tree_bwt.count(next_insert, 0);

            if start_index != next_insert {
                let original_ni: u64 = next_insert;

                //short circuit failed
                for pred_symbol in int_form.iter().rev() {
                    next_insert = self.tree_bwt.count(next_insert, *pred_symbol)+self.start_index[*pred_symbol as usize];
                }
                next_insert = self.tree_bwt.count(next_insert, 0);
                
                if original_ni == next_insert {
                    //the full search did nothing this is a copy sequence
                    //cutting down the search will save on the initial k-mer query
                    self.sort_query_len -= 2.0 * COST_FACTOR * query_len as f64;
                    self.short_circuits[1] += 1;
                } else {
                    //the full search did make it more specific, so we should increase the query short circuit size
                    //this could save at most the length of the full sequence
                    self.sort_query_len += COST_FACTOR * int_form.len() as f64;
                    self.short_circuits[2] += 1;
                }
            } else {
                //short circuit success, making it smaller will save 2 bases queries, and we also have downweighted this
                //self.sort_query_len -= 0.2 * COST_FACTOR;
                self.short_circuits[0] += 1;
            }
        } else {
            next_insert = self.string_count;
        }
        
        //go through the characters in reverse
        let mut symbol: u8 = 0; // $
        for pred_symbol in int_form.iter().rev() {
            next_insert = self.tree_bwt.insert_and_count(next_insert, *pred_symbol);
            self.symbol_counts[*pred_symbol as usize] += 1;
            for i in (symbol+1) as usize..VC_LEN {
                self.start_index[i] += 1;
            }

            //after any adjustments add in the new offset for the symbol we're currently at
            next_insert += self.start_index[*pred_symbol as usize];
            symbol = *pred_symbol;
        }

        //one final insert for the $
        self.tree_bwt.insert_and_count(next_insert, 0);
        self.symbol_counts[0] += 1;
        for i in (symbol+1) as usize..VC_LEN {
            self.start_index[i] += 1;
        }

        self.total_count += (int_form.len()+1) as u64;
        self.string_count += 1;

        if self.string_count % 10000 == 0 {
            println!("{} {:?} {}", self.string_count, self.short_circuits, self.sort_query_len);
            self.short_circuits = [0; 3];
        }
    }

    /// This will return the data in a plain Vector format, with one symbol per index.
    /// # Examples
    /// ```rust
    /// use msbwt2::dynamic_bwt::DynamicBWT;
    /// let data: String = "ACGNT".to_string();
    /// let bwt: Vec<u8> = vec![5, 0, 1, 2, 3, 4];
    /// let mut ubwt: DynamicBWT = Default::default();
    /// ubwt.insert_string(&data, false);
    /// assert_eq!(ubwt.to_vec(), bwt);
    /// ```
    #[inline]
    pub fn to_vec(&self) -> Vec<u8> {
        self.tree_bwt.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bwt_converter::*;
    use crate::bwt_util::naive_bwt;
    use crate::string_util;
    use tempfile::{Builder, NamedTempFile};

    #[test]
    fn test_init() {
        let ubwt: DynamicBWT = Default::default();
        assert_eq!(ubwt.to_vec(), Vec::<u8>::new());
    }

    #[test]
    fn test_single_string() {
        let data: String = "ACGNT".to_string();
        let bwt: Vec<u8> = vec![5, 0, 1, 2, 3, 4];
        let mut ubwt: DynamicBWT = Default::default();
        ubwt.insert_string(&data, false);
        assert_eq!(ubwt.to_vec(), bwt);
    }

    #[test]
    fn test_multi_string_unsorted() {
        let mut data: Vec<&str> = vec!["CCGT", "ACG", "N"];
        let bwt = convert_stoi(&"GTN$$ACCC$G");
        
        //to get an identical result to sorted rope, we have to sort
        data.sort();
        
        //insert the strings in the now sorted order
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, false);
        }
        assert_eq!(ubwt.to_vec(), bwt);
    }

    #[test]
    fn test_multi_string_sorted() {
        let data: Vec<&str> = vec!["ACG", "N", "CCGT", "N", "ACG", "ACG", "CCGT", "N"];
        let bwt = string_util::convert_stoi(&naive_bwt(&data));
        
        //insert the strings in the now sorted order
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), bwt);
    }

    #[test]
    fn test_multi_length() {
        //getting bigger in order
        let data: Vec<&str> = vec!["A", "AA", "AAA", "AAAA", "AAAAA"];
        let bwt = string_util::convert_stoi(&naive_bwt(&data));
        
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), bwt);

        //getting smaller in order
        let data: Vec<&str> = vec!["AAAAA", "AAAA", "AAA", "AA", "A"];
        let bwt = string_util::convert_stoi(&naive_bwt(&data));
        
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), bwt);
    }

    #[test]
    fn test_sampled_bwt() {
        let genome: String = "ACCGTGTTGCCGTAGTGAAAAGTGACGACGTGAGATGGCCAAAGTGGGTCTCTGTG".to_string();
        let read_length: usize = 20;
        let coverage: usize = 32;//make sure we get some runs
        //let read_length: usize = 5;
        //let coverage: usize = 1;
        let mut data: Vec<&str> = vec![];
        for s in 0..genome.len()-read_length {//-43 {
            for _ in 0..coverage {
                data.push(&genome[s..s+read_length]);
            }
        }

        //let data: Vec<&str> = vec!["ACCGT", "TGCCGC"];
        //let data: Vec<&str> = vec!["ACCGT", "CCGTG"];
        //println!("data: {:?}", data);

        //use this function to make a bwt the naive (i.e. slow) way
        let naive = string_util::convert_stoi(&naive_bwt(&data));
        
        //insert the strings in the now sorted order
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), naive);
    }

    #[test]
    fn test_load_dynamicbwt_from_vec() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        
        //stream and compress the BWT
        //let bwt_stream = stream_bwt_from_fastqs(&fastq_filenames).unwrap();
        let bwt_stream = naive_bwt(&data);
        let compressed_bwt = convert_to_vec(bwt_stream.as_bytes());
        
        let mut bwt = DynamicBWT::new();
        bwt.load_vector(compressed_bwt);

        let expected_totals = vec![3, 1, 3, 2, 1, 1];
        for i in 0..6 {
            //make sure the total counts are correct
            assert_eq!(bwt.get_symbol_count(i as u8), expected_totals[i]);
        }
    }

    #[test]
    fn test_load_dynamicbwt_from_npy() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        
        //stream and compress the BWT
        //let bwt_stream = stream_bwt_from_fastqs(&fastq_filenames).unwrap();
        let bwt_stream = naive_bwt(&data);
        let compressed_bwt = convert_to_vec(bwt_stream.as_bytes());
        
        //save the output to a temporary numpy file
        let bwt_file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".npy").tempfile().unwrap();
        let filename: String = bwt_file.path().to_str().unwrap().to_string();
        save_bwt_numpy(&compressed_bwt[..], &filename).unwrap();
        
        //load it back in and verify counts
        let mut bwt = DynamicBWT::new();
        bwt.load_numpy_file(&filename).unwrap();

        let expected_totals = vec![3, 1, 3, 2, 1, 1];
        for i in 0..6 {
            //make sure the total counts are correct
            assert_eq!(bwt.get_symbol_count(i as u8), expected_totals[i]);
        }
    }

    #[test]
    fn test_constrain_range() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        
        //stream and compress the BWT
        let bwt_stream = naive_bwt(&data);
        assert_eq!(bwt_stream, "GTN$$ACCC$G");
        let bwt_int_form = string_util::convert_stoi(&bwt_stream);
        let compressed_bwt = convert_to_vec(bwt_stream.as_bytes());
        //[G, T, N, 2$, A, 3C, $, G]
        assert_eq!(compressed_bwt.len(), 8);
        
        //load it back in and verify counts
        let mut bwt: DynamicBWT = Default::default();
        bwt.load_vector(compressed_bwt.clone());

        let initial_range = BWTRange {
            l: 0,
            h: bwt_stream.len() as u64
        };
        
        //this is verifying that all single-symbol queries get the start/end range
        for sym in 0..VC_LEN {
            let new_range = unsafe {
                bwt.constrain_range(sym as u8, &initial_range)
            };
            assert_eq!(new_range, BWTRange{
                l: bwt.start_index[sym] as u64, 
                h: (bwt.start_index[sym]+bwt.symbol_counts[sym]) as u64
            });
        }

        //now lets verify that we get all ascending symbols
        for sym in 0..VC_LEN {
            let mut sym_count = 0;
            for ind in 0..(bwt_stream.len()+1) {
                //test from 0 to the current index
                let initial_range = BWTRange {
                    l: 0,
                    h: ind as u64
                };

                let new_range = unsafe {
                    bwt.constrain_range(sym as u8, &initial_range)
                };
                assert_eq!(new_range, BWTRange {
                    l: bwt.start_index[sym] as u64,
                    h: (bwt.start_index[sym]+sym_count) as u64
                });

                //test from the current index to the high point
                let initial_range = BWTRange {
                    l: ind as u64,
                    h: bwt_stream.len() as u64
                };

                let new_range = unsafe {
                    bwt.constrain_range(sym as u8, &initial_range)
                };
                assert_eq!(new_range, BWTRange {
                    l: (bwt.start_index[sym]+sym_count) as u64,
                    h: (bwt.start_index[sym]+bwt.symbol_counts[sym]) as u64
                });

                //check if we need to adjust our expected values at all
                if ind < bwt_stream.len() && bwt_int_form[ind] == sym as u8 {
                    sym_count += 1;
                }
            }
        }
    }

    #[test]
    fn test_count_kmer() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let data: Vec<&str> = vec!["CCGTACGTA", "GGTACAGTA", "ACGACGACG"];
        
        //stream and compress the BWT
        let bwt_stream = naive_bwt(&data);
        let compressed_bwt = convert_to_vec(bwt_stream.as_bytes());
        
        //load it back in and verify counts
        let mut bwt: DynamicBWT = Default::default();
        bwt.load_vector(compressed_bwt.clone());

        //simple sanity checks, make sure our single-character symbols matches the total count
        for c in 0..VC_LEN as u8 {
            let test_seq: Vec<u8> = vec![c];
            assert_eq!(bwt.get_symbol_count(c), bwt.count_kmer(&test_seq));
        }
        
        //check that each string shows up once
        for seq in data.iter() {
            let test_seq = string_util::convert_stoi(seq);
            assert_eq!(bwt.count_kmer(&test_seq), 1);
        }

        //now lets check some semi-arbitrary substrings
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACG")), 4);
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"CC")), 1);
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"TAC")), 2);
    }

    #[test]
    fn test_load_and_add() {
        //strings - "CCGT\nACG\nN"
        //build the BWT
        let mut data: Vec<&str> = vec!["CCGTACGTA", "GGTACAGTA", "ACGACGACG"];
        
        //stream and compress the BWT
        //let bwt_stream = stream_bwt_from_fastqs(&fastq_filenames).unwrap();
        let bwt_stream = naive_bwt(&data);
        let compressed_bwt = convert_to_vec(bwt_stream.as_bytes());
        
        //load it back in and verify counts
        let mut bwt: DynamicBWT = Default::default();
        bwt.load_vector(compressed_bwt.clone());

        // now lets add a new string
        let new_string = "AAGTCATAT";
        bwt.insert_string(&new_string, true);
        data.push(new_string);

        //simple sanity checks, make sure our single-character symbols matches the total count
        for c in 0..VC_LEN as u8 {
            let test_seq: Vec<u8> = vec![c];
            assert_eq!(bwt.get_symbol_count(c), bwt.count_kmer(&test_seq));
        }
        
        //check that each string shows up once
        for seq in data.iter() {
            let test_seq = string_util::convert_stoi(seq);
            assert_eq!(bwt.count_kmer(&test_seq), 1);
        }

        //now lets check some semi-arbitrary substrings
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACG")), 4);
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"CC")), 1);
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"TAC")), 2);

        //these should have changed with the new string
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"AA")), 1);
        assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"GT")), 5);
    }
}
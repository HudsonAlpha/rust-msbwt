
// Not all of these are directly tied to the msbwt "core", but it's better to keep them together IMO
/// The number of characters in our alphabet
pub const VC_LEN: usize = 6;      //$ A C G N T
/// The number of bits for storing the character in a byte
pub const LETTER_BITS: usize = 3; //defined
/// The number of bit for storing quantity in a byte
pub const NUMBER_BITS: usize = 5; //8-letterBits
/// Multiplier for multi-byte runs
pub const NUM_POWER: usize = 32;  //2**numberBits
/// Contains the character mask
pub const MASK: u8 = 0x07;        //255 >> numberBits
/// Contains the right-shifted number mask
pub const COUNT_MASK: u8 = 0x1F;

/// Basic struct for containing a range in a BWT.
/// Only contains fields `l` and `h`, representing a range [l, h).
#[derive(Clone,Copy,Default,Debug,Eq,PartialEq)]
pub struct BWTRange {
    /// the lower bound, inclusive
    pub l: u64,
    /// the upper bound, exclusive
    pub h: u64
}

/// This is the primary functionality that will be shared between BWT types.
/// The main functions are loading data and querying for k-mers.
pub trait BWT {
    /// Initializes the BWT from a compressed BWT vector.
    /// # Arguments
    /// * `bwt` - the run-length encoded BWT stored in a Vec<u8> 
    /// # Examples
    /// ```rust
    /// use msbwt2::msbwt_core::BWT;
    /// use msbwt2::rle_bwt::RleBWT;
    /// use msbwt2::bwt_converter::convert_to_vec;
    /// //strings "ACGT" and "CCGG"
    /// let seq = "TG$$CAGCCG";
    /// let vec = convert_to_vec(seq.as_bytes());
    /// let mut bwt = RleBWT::new();
    /// bwt.load_vector(vec);
    /// ```
    fn load_vector(&mut self, bwt: Vec<u8>);

    /// Initializes the BWT from the numpy file format for compressed BWTs
    /// # Arguments
    /// * `filename` - the name of the file to load into memory
    /// # Examples
    /// ```rust
    /// use msbwt2::msbwt_core::BWT;
    /// use msbwt2::rle_bwt::RleBWT;
    /// use msbwt2::string_util;
    /// let mut bwt = RleBWT::new();
    /// let filename: String = "test_data/two_string.npy".to_string();
    /// bwt.load_numpy_file(&filename);
    /// assert_eq!(bwt.count_kmer(&string_util::convert_stoi(&"ACGT")), 1);
    /// ```
    fn load_numpy_file(&mut self, filename: &str) -> std::io::Result<()>;

    /// Returns the total number of occurences of a given symbol
    /// # Arguments
    /// * `symbol` - the symbol in integer form
    /// # Examples
    /// ```rust
    /// # use msbwt2::msbwt_core::BWT;
    /// # use msbwt2::rle_bwt::RleBWT;
    /// # use msbwt2::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let vec = convert_to_vec(seq.as_bytes());
    /// # let mut bwt = RleBWT::new();
    /// # bwt.load_vector(vec);
    /// let string_count = bwt.get_symbol_count(0);
    /// assert_eq!(string_count, 2);
    /// ```
    fn get_symbol_count(&self, symbol: u8) -> u64;

    /// This will return the total number of symbols contained by the BWT
    /// # Examples
    /// ```rust
    /// # use msbwt2::msbwt_core::BWT;
    /// # use msbwt2::rle_bwt::RleBWT;
    /// # use msbwt2::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let vec = convert_to_vec(seq.as_bytes());
    /// # let mut bwt = RleBWT::new();
    /// # bwt.load_vector(vec);
    /// let total_size = bwt.get_total_size();
    /// assert_eq!(total_size, 10);
    /// ```
    fn get_total_size(&self) -> u64;

    /// Performs a range constraint on a BWT range. This implicitly represents prepending a character `sym` to a k-mer
    /// represented by `input_range` to create a new range representing a (k+1)-mer.
    /// # Arguments
    /// * `sym` - the symbol to pre-pend in integer form
    /// * `input_range` - the range to pre-pend to
    /// # Safety
    /// This function is unsafe because there are no guarantees that the symbol or bounds will be checked by the implementing structure.
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange;

    /// Returns the total number of occurrences of a given k-mer in the BWT.
    /// # Arguments
    /// * `kmer` - the integer-encoded kmer sequence to count
    /// # Examples
    /// ```rust
    /// # use msbwt2::msbwt_core::BWT;
    /// # use msbwt2::rle_bwt::RleBWT;
    /// # use msbwt2::bwt_converter::convert_to_vec;
    /// # use msbwt2::string_util;
    /// # let seq = "TG$$CAGCCG";
    /// # let vec = convert_to_vec(seq.as_bytes());
    /// # let mut bwt = RleBWT::new();
    /// # bwt.load_vector(vec);
    /// //strings "ACGT" and "CCGG"
    /// let kmer_count = bwt.count_kmer(&vec![1, 2, 3, 5]); //ACGT
    /// assert_eq!(kmer_count, 1);
    /// let kmer_count = bwt.count_kmer(&vec![2, 3]); //CG
    /// assert_eq!(kmer_count, 2);
    /// //if you want to go directly from a string, use stoi to help
    /// let kmer_count = bwt.count_kmer(&string_util::convert_stoi(&"ACGT"));
    /// assert_eq!(kmer_count, 1);
    /// ```
    #[inline]
    fn count_kmer(&self, kmer: &[u8]) -> u64 {
        //init to everything
        assert!(kmer.iter().all(|&v| v < VC_LEN as u8));
        let mut ret: BWTRange = BWTRange {
            l: 0,
            h: self.get_total_size()
        };
        
        /*
        //check for cache entry
        let cut_kmer: &[u8];
        if kmer.len() >= self.cache_k {
            ret = self.kmer_cache[self.get_cache_index(&kmer[kmer.len()-self.cache_k..])];
            cut_kmer = &kmer[..kmer.len()-self.cache_k];
        } else {
            cut_kmer = kmer;
            ret = BWTRange {
                l: 0,
                h: self.total_size
            };
        }*/
        
        //go through what remains in reverse
        //for c in cut_kmer.iter().rev() {
        for c in kmer.iter().rev() {
            if ret.h == ret.l {
                return 0;
            }
            unsafe {
                ret = self.constrain_range(*c, &ret);
            }
        }

        //return the delta
        ret.h-ret.l
    }
}
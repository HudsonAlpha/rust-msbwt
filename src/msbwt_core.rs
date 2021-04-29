
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

pub trait BWT {
    /// Initializes the BWT from a compressed BWT vector.
    /// # Arguments
    /// * `bwt` - the run-length encoded BWT stored in a Vec<u8> 
    /// # Examples
    /// ```rust
    /// use std::io::Cursor;
    /// use msbwt::msbwt_core::BWT;
    /// use msbwt::rle_bwt::RleBWT;
    /// use msbwt::bwt_converter::convert_to_vec;
    /// //strings "ACGT" and "CCGG"
    /// let seq = "TG$$CAGCCG";
    /// let seq = Cursor::new(seq);
    /// let vec = convert_to_vec(seq);
    /// let mut bwt = RleBWT::new();
    /// bwt.load_vector(vec);
    /// ```
    fn load_vector(&mut self, bwt: Vec<u8>);

    /// Initializes the BWT from the numpy file format for compressed BWTs
    /// # Arguments
    /// * `filename` - the name of the file to load into memory
    /// # Examples
    /// ```ignore
    /// use msbwt::rle_bwt::RleBWT;
    /// let mut bwt = RleBWT::new();
    /// let filename: String = "/path/to/my/file/comp_msbwt.npy".to_string();
    /// bwt.load_numpy_file(&filename);
    /// ```
    fn load_numpy_file(&mut self, filename: &str) -> std::io::Result<()>;

    /// Returns the total number of occurences of a given symbol
    /// # Arguments
    /// * `symbol` - the symbol in integer form
    /// # Examples
    /// ```rust
    /// # use std::io::Cursor;
    /// # use msbwt::msbwt_core::BWT;
    /// # use msbwt::rle_bwt::RleBWT;
    /// # use msbwt::bwt_converter::convert_to_vec;
    /// # let seq = "TG$$CAGCCG";
    /// # let seq = Cursor::new(seq);
    /// # let vec = convert_to_vec(seq);
    /// # let mut bwt = RleBWT::new();
    /// # bwt.load_vector(vec);
    /// let string_count = bwt.get_total_counts(0);
    /// assert_eq!(string_count, 2);
    /// ```
    fn get_total_counts(&self, symbol: u8) -> u64;

    /// Performs a range constraint on a BWT range. This implicitly represents prepending a character `sym` to a k-mer
    /// represented by `input_range` to create a new range representing a (k+1)-mer.
    /// # Arguments
    /// * `sym` - the symbol to pre-pend in integer form
    /// * `input_range` - the range to pre-pend to
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange;
}
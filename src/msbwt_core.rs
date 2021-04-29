
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
    /// Performs a range constraint on a BWT range. This implicitly represents prepending a character `sym` to a k-mer
    /// represented by `input_range` to create a new range representing a (k+1)-mer.
    /// # Arguments
    /// * `sym` - the symbol to pre-pend in integer form
    /// * `input_range` - the range to pre-pend to
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange;
}
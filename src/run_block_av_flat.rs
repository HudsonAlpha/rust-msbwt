
use arrayvec::ArrayVec;
use likely_stable::{likely,unlikely};

/// This is the number of bytes used for data storage.
/// When it reaches this size, the block should be split to avoid issues.
pub const MAX_BLOCK_SIZE: usize = 63; // best for current B+
const CAPACITY_BUFFER: usize = MAX_BLOCK_SIZE+2;

/// The number of characters in our alphabet
pub const VC_LEN: usize = 6;      //$ A C G N T
/// the number of bits that are encoding the symbol
const SYMBOL_BITS: usize = 3;
/// Contains the right-shifted number mask
const SYMBOL_MASK: u16 = 0x07;
/// Stores the set bits for a half full block, e.g. (2**12) << 3
const HALF_FULL: u16 = 0x8000;
/// Stores the set bits for a single count block e.g. 1 << 3
const SINGLE_COUNT: u16 = 0x0008;

/// A run-length encoded block of data implemented with an ArrayVec
#[derive(Clone,Debug,PartialEq)]
pub struct RLEBlock {
    runs: ArrayVec<u16, CAPACITY_BUFFER>,
    symbol_counts: [u64; VC_LEN],
    values_contained: u64
}

impl Default for RLEBlock {
    /// Default function
    #[inline]
    fn default() -> Self {
        let runs = ArrayVec::<u16, CAPACITY_BUFFER>::new();
        Self {
            runs,
            symbol_counts: [0; VC_LEN],
            values_contained: 0
        }
    }
}
/*
#[inline] 
fn calc_run_bytes(count: u64) -> usize {
    if count < (1_u64 << 13) {
        1
    } else {
        2
    }
}
*/
#[inline]
fn encode_run(symbol: u8, count: u16) -> u16 {
    //assert!(count < (1_u16 << 13));
    //1 bit for length, 3 bits for symbol, leaves 16-4 = 12 for length
    (symbol as u16) | (count << SYMBOL_BITS) as u16
}

#[inline]
fn decode_run(run: u16) -> (u8, u16) {
    (
        (run & SYMBOL_MASK) as u8, 
        (run >> SYMBOL_BITS)
    )
}

impl RLEBlock {
    /// Returns the total number of symbols encoded by this RLEBlock
    #[inline]
    pub fn get_values_contained(&self) -> u64 {
        self.values_contained
    }

    /// Returns the symbol counts encoded by this RLEBlock
    #[inline]
    pub fn get_symbol_counts(&self) -> [u64; VC_LEN] {
        self.symbol_counts
    }

    /// Returns the current usage of the RLEBlock's underlying storage.
    /// If this reaches the `MAX_BLOCK_SIZE`, then the block should be split.
    #[inline]
    pub fn block_len(&self) -> usize {
        self.runs.len()
    }

    /// Performs a rank/count operation.
    /// For a given data index, this will count all occurences of symbol `value` up to that index.
    /// # Arguments
    /// * `position` - the index to count to 
    /// * `symbol` - the symbol to count
    /// # Example
    /// ```rust
    /// use msbwt2::run_block_av_flat::RLEBlock;
    /// let mut block: RLEBlock = Default::default();
    /// //"data" is inserted in-order
    /// let data: Vec<u8> =            vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
    /// //"expected_count" is the number of occurences of the inserted symbol BEFORE this one
    /// let expected_count: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
    /// for (i, v) in data.iter().enumerate() {
    ///     assert_eq!(block.count(i as u64, *v), expected_count[i]);
    ///     assert_eq!(block.insert_and_count(i as u64, *v), expected_count[i]);
    /// }
    /// assert_eq!(block.to_vec(), data);
    /// ```
    #[inline]
    pub fn count(&self, position: u64, symbol: u8) -> u64 {
        //make sure we're inserting valid symbols into a valid position
        assert!(symbol < VC_LEN as u8);
        assert!(position <= self.values_contained);
        
        //its more efficient to count all symbols and return the one we care about than do checks inside the loop
        let mut pos_end: u64 = 0;
        let mut total_counts: [u64; VC_LEN] = [0; VC_LEN];
        
        //now iterate
        let mut i: usize = 0;
        let mut run: (u8, u16) = (0, 0);
        while pos_end < position {
            //get the symbol and the added counts
            run = decode_run(self.runs[i]);

            //add the count to the total so far
            total_counts[run.0 as usize] += run.1 as u64;
            pos_end += run.1 as u64;
            i += 1;
        }
        //we potentially went past the target, so subtract those back out
        total_counts[run.0 as usize] -= pos_end - position;
        
        //return the accumulated counts for the symbol we care about
        total_counts[symbol as usize]
    }

    /// Performs a rank/count operation while also inserting that symbol at the provided index.
    /// For a given data index, this will count all occurences of symbol `value` up to that index, and then insert an addition `value` at that index.
    /// # Arguments
    /// * `position` - the index to count to 
    /// * `symbol` - the symbol to count
    /// # Example
    /// ```rust
    /// use msbwt2::run_block_av_flat::RLEBlock;
    /// let mut block: RLEBlock = Default::default();
    /// //"data" is inserted in-order
    /// let data: Vec<u8> =            vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
    /// //"expected_count" is the number of occurences of the inserted symbol BEFORE this one
    /// let expected_count: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
    /// for (i, v) in data.iter().enumerate() {
    ///     assert_eq!(block.count(i as u64, *v), expected_count[i]);
    ///     assert_eq!(block.insert_and_count(i as u64, *v), expected_count[i]);
    /// }
    /// assert_eq!(block.to_vec(), data);
    /// ```
    #[inline]
    pub fn insert_and_count(&mut self, position: u64, symbol: u8) -> u64 {
        //make sure we're inserting valid symbols into a valid position
        assert!(symbol < VC_LEN as u8);
        assert!(position <= self.values_contained);

        //adding a value will always increment these counts
        self.symbol_counts[symbol as usize] += 1;
        self.values_contained += 1;

        //we need to find the run we're going into
        let mut pos_end: u64 = 0;
        let mut total_counts: [u64; VC_LEN+1] = [0; VC_LEN+1];
        
        //now iterate
        let mut i: usize = 0;
        let mut run: (u8, u16) = (VC_LEN as u8, 0);
        while pos_end < position {
            //get the symbol and the added counts
            run = decode_run(self.runs[i]);

            //add the count to the total so far
            total_counts[run.0 as usize] += run.1 as u64;
            pos_end += run.1 as u64;
            i += 1;
        }

        //we potentially went past the target, so subtract those back out
        total_counts[run.0 as usize] -= pos_end - position;

        //first, find the end of the run we're currently in
        if run.0 == symbol {
            //the run we're looking at is the same as what we're inserting
            self.increment_run(i-1);
        } else if position < pos_end {
            //our insert is in the middle of the run and different, we need to split the run
            let counts_after: u16 = (pos_end - position) as u16;
            let counts_before: u16 = run.1 - counts_after;
            self.runs[i-1] = encode_run(run.0, counts_before);

            //let new_length = 1+calc_run_bytes(counts_before)+calc_run_bytes(counts_after) - run.2;
            let new_length = 2;

            //this method only shifts the values once, so more efficient than double insert
            self.runs.try_extend_from_slice(&vec![0; new_length]).unwrap();
            for j in (i+new_length..self.runs.len()).rev() {
                self.runs[j] = self.runs[j-new_length];
            }
            self.runs[i] = encode_run(symbol, 1);
            self.runs[i+1] = encode_run(run.0, counts_after);
        } else if i < self.runs.len() && (self.runs[i] & SYMBOL_MASK) as u8 == symbol {
            //right on a run boundary AND the next run matches our symbol, so increment the next run
            self.increment_run(i);
        } else {
            //right on a position boundary (or at the very end) and the next run does not match, so do an insert
            //self.runs.insert(i, symbol | SINGLE_COUNT);
            self.runs.insert(i, (symbol as u16) | SINGLE_COUNT);
            //encode_run(&mut self.runs[i..], symbol, 1);
        }

        //return the accumulated counts for the symbol we care about
        total_counts[symbol as usize]
    }

    #[inline]
    fn increment_run(&mut self, run_index: usize) {
        /*
        let curr_run = decode_run(&self.runs[run_index..]);
        let curr_bytes = curr_run.2;
        let next_bytes = calc_run_bytes(curr_run.1+1);
        
        for _ in 0..next_bytes-curr_bytes {
            self.runs.insert(run_index+1, 0);
        }

        encode_run(&mut self.runs[run_index..], curr_run.0, curr_run.1+1);
        */
        //increment
        let overflow = self.runs[run_index].overflowing_add(SINGLE_COUNT);
        self.runs[run_index] = overflow.0;
        if unlikely(overflow.1) {
            println!("split");
            //if it overflows, split into two halfs at the same position
            self.runs[run_index] |= HALF_FULL;
            self.runs.insert(run_index+1, self.runs[run_index]);
        }
    }
    
    #[inline]
    pub fn split(&mut self) -> RLEBlock {
        //break right in the middle
        let breakpoint: usize = self.runs.len() / 2;
        let mut new_block_counts: [u64; VC_LEN] = [0; VC_LEN];
        let mut new_block_size: u64 = 0;
        
        //drain from the breakpoint, counting the symbols as we go
        let new_runs_data: ArrayVec<u16, CAPACITY_BUFFER> = self.runs.drain(breakpoint..)
            .map(
                |v| {
                    let run = decode_run(v);
                    new_block_counts[run.0 as usize] += run.1 as u64;
                    v
                }
            ).collect();

        //count up the total also
        new_block_size += new_block_counts.iter().sum::<u64>();

        //build the block
        let ret: RLEBlock = RLEBlock {
            runs: new_runs_data,
            symbol_counts: new_block_counts,
            values_contained: new_block_size
        };

        //now update internals before returning
        for (i, &nbc) in new_block_counts.iter().enumerate() {
            self.symbol_counts[i] -= nbc;
        }
        self.values_contained -= new_block_size;
        ret
        
        /*
        //first figure out what this block will change to
        let mut breakpoint = 0;
        let mut block_counts: [u64; VC_LEN] = [0; VC_LEN];
        while breakpoint < self.runs.len() / 2 {
            let run = decode_run(&self.runs[breakpoint..]);
            block_counts[run.0 as usize] += run.1;
            breakpoint += 1;
        }
        let block_size: u64 = block_counts.iter().sum();

        //now figure out the new block
        let new_block_size: u64 = self.values_contained - block_size;
        let mut new_block_counts: [u64; VC_LEN] = [0; VC_LEN];
        for (i, &bc) in block_counts.iter().enumerate() {
            new_block_counts[i] = self.symbol_counts[i] - bc;
        }
        let new_runs_data: ArrayVec<u16, CAPACITY_BUFFER> = self.runs.drain(breakpoint..).collect();
        //let mut new_runs_data: ArrayVec<u8, CAPACITY_BUFFER> = Default::default();
        //new_runs_data.try_extend_from_slice(&self.runs[breakpoint..]).unwrap();
        //self.runs.truncate(breakpoint);

        //build the block
        let ret: RLEBlock = RLEBlock {
            runs: new_runs_data,
            symbol_counts: new_block_counts,
            values_contained: new_block_size
        };

        //finally update this block
        self.symbol_counts = block_counts;
        self.values_contained = block_size;

        //now return
        ret
        */
    }

    /// Returns a vector representation of the contained data, one symbol per entry
    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = Vec::<u8>::with_capacity(self.values_contained as usize);
        
        //now iterate
        let mut i = 0;
        while i < self.runs.len() {
            //get the symbol and the added counts
            let run = decode_run(self.runs[i]);
            
            //add the run symbol to the vec
            for _ in 0..run.1 {
                ret.push(run.0 as u8);
            }
            i += 1;
        }
        
        ret
    }

    /*
    /// Returns the raw data iterator, one (u8, u8) run block per entry
    pub fn raw_iter(&self) -> std::slice::Iter<(u8, u8)> {
        self.runs.iter().map(
            |&b| {
                let d = decode_run(&[b]);
                (b.0, b.1 as u8)
            }
        )
    }
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init() {
        let block: RLEBlock = Default::default();
        assert_eq!(block.runs.to_vec(), Vec::<u16>::new());
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), vec![]);
        for symbol in 0..VC_LEN {
            assert_eq!(0, block.count(0, symbol as u8));
        }
    }
    
    
    #[test]
    fn test_insert() {
        let mut block: RLEBlock = Default::default();
        assert_eq!(block.count(0, 0), 0);
        assert_eq!(block.insert_and_count(0, 0), 0);

        for _ in 0..256 {
            assert_eq!(block.count(1, 1), 0);
            assert_eq!(block.insert_and_count(1, 1), 0);
        }
        
        for _ in 0..256 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        
        for _ in 0..255 {
            assert_eq!(block.count(0, 2), 0);
            assert_eq!(block.insert_and_count(0, 2), 0);
        }
        
        //let correct_data = vec![(2, 255), (0, 129), (0, 128), (1, 128), (1, 128)];
        //assert_eq!(block.runs.to_vec(), correct_data);
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }
    
    #[test]
    fn test_insert_splitting() {
        //middle split
        let mut block: RLEBlock = Default::default();
        for _ in 0..255 {
            assert_eq!(block.count(0, 0), 0);
            assert_eq!(block.insert_and_count(0, 0), 0);
        }
        assert_eq!(block.count(128, 1), 0);
        assert_eq!(block.insert_and_count(128, 1), 0);

        //let correct_data = vec![(0, 128), (1, 1), (0, 127)];
        let mut correct_data = vec![0; 128];
        correct_data.extend_from_slice(&vec![1]);
        correct_data.extend_from_slice(&vec![0; 127]);
        assert_eq!(block.to_vec(), correct_data);
    }
    
    #[test]
    fn test_increment_run() {
        //simple multiple block extension
        let mut block: RLEBlock = Default::default();
        block.insert_and_count(0, 0);
        for _ in 0..256 {
            block.increment_run(0);
        }

        let correct_data = vec![0; 257];//vec![(0, 129), (0, 128)];
        assert_eq!(block.to_vec(), correct_data);
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }

    /*
    #[test]
    fn test_block_splits() {
        
        //middle split
        let mut block: RLEBlock = Default::default();
        for _ in 0..512 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(256, 1);
        let new_block: RLEBlock = block.split();
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(0, 128), (0, 128)]).unwrap();
        assert_eq!(block, RLEBlock {
            runs,
            symbol_counts: [256, 0, 0, 0, 0, 0],
            values_contained: 256
        });
        
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(1, 1), (0, 128), (0, 128)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            runs,
            symbol_counts: [256, 1, 0, 0, 0, 0],
            values_contained: 257
        });

        //early split
        let mut block: RLEBlock = Default::default();
        for _ in 0..256 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(1, 1);
        let new_block: RLEBlock = block.split();
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(0, 1), (1, 1)]).unwrap();
        assert_eq!(block, RLEBlock {
            runs,
            symbol_counts: [1, 1, 0, 0, 0, 0],
            values_contained: 2
        });
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(0, 127), (0, 128)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            runs,
            symbol_counts: [255, 0, 0, 0, 0, 0],
            values_contained: 255
        });
        
        //late split
        let mut block: RLEBlock = Default::default();
        for _ in 0..256 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(255, 1);
        let new_block: RLEBlock = block.split();
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(0, 128), (0, 127)]).unwrap();
        assert_eq!(block, RLEBlock {
            runs,
            symbol_counts: [255, 0, 0, 0, 0, 0],
            values_contained: 255
        });
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(1, 1), (0, 1)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            runs,
            symbol_counts: [1, 1, 0, 0, 0, 0],
            values_contained: 2
        });
    }
    */

    #[test]
    fn test_to_vec() {
        //big run test
        let mut block: RLEBlock = Default::default();
        for _ in 0..1024 {
            block.insert_and_count(0, 1);
        }
        assert_eq!(block.to_vec(), vec![1; 1024]);

        //semi-random collection test
        let mut block: RLEBlock = Default::default();
        let data: Vec<u8> = vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
        let expected_count: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
        for (i, v) in data.iter().enumerate() {
            assert_eq!(block.count(i as u64, *v), expected_count[i]);
            assert_eq!(block.insert_and_count(i as u64, *v), expected_count[i]);
        }
        assert_eq!(block.to_vec(), data);
    }
    
    #[test]
    fn test_block_failure_001() {
        let mut block: RLEBlock = Default::default();
        let mut data: Vec<u8> = vec![];
        let inserted: Vec<u8> = vec![4, 4, 4, 5, 2, 4, 4];
        let positions: Vec<u64> = vec![0, 0, 1, 1, 4, 1, 5];
        for i in 0..inserted.len() {
            data.insert(positions[i] as usize, inserted[i]);
            let mut expected_count = 0;
            for j in 0..positions[i] {
                if data[j as usize] == inserted[i] {
                    expected_count += 1;
                }
            }

            let pre_count = block.count(positions[i], inserted[i]);
            assert_eq!(pre_count, expected_count);

            let count = block.insert_and_count(positions[i], inserted[i]);
            assert_eq!(block.to_vec(), data);
            assert_eq!(count, expected_count);
        }
    }
}
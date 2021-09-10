
use arrayvec::ArrayVec;

pub const MAX_BLOCK_SIZE: usize = 63; // best for current B+
const CAPACITY_BUFFER: usize = MAX_BLOCK_SIZE+2;

/// The number of characters in our alphabet
pub const VC_LEN: usize = 6;      //$ A C G N T
/// The number of bit for storing quantity in a byte
pub const NUMBER_BITS: usize = 8; //8-letterBits
/// Contains the right-shifted number mask
pub const COUNT_MASK: u8 = 0xFF;
pub const HALF_FULL: u8 = 128;
/// Multiplier for multi-byte runs
pub const NUM_POWER: usize = 256;  //2**numberBits

//TODO: if this is useful, change the name
#[derive(Clone,Debug,PartialEq)]
pub struct RLEBlock {
    runs: ArrayVec<(u8, u8), CAPACITY_BUFFER>,
    symbol_counts: [u64; VC_LEN],
    values_contained: u64
}

impl Default for RLEBlock {
    /// Default function
    #[inline]
    fn default() -> Self {
        let runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        Self {
            runs,
            symbol_counts: [0; VC_LEN],
            values_contained: 0
        }
    }
}

impl RLEBlock {
    #[inline]
    pub fn get_values_contained(&self) -> u64 {
        self.values_contained
    }

    #[inline]
    pub fn get_symbol_counts(&self) -> [u64; VC_LEN] {
        self.symbol_counts
    }

    #[inline]
    pub fn block_len(&self) -> usize {
        self.runs.len()
    }

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
        let mut sym: u8 = 0;
        let mut count: u8;
        while pos_end < position {
            sym = self.runs[i].0;
            count = self.runs[i].1;
            total_counts[sym as usize] += count as u64;
            pos_end += count as u64;
            i += 1;
        }
        //we potentially went past the target, so subtract those back out
        total_counts[sym as usize] -= pos_end - position;
        
        //return the accumulated counts for the symbol we care about
        total_counts[symbol as usize]
    }

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
        let mut sym: u8 = VC_LEN as u8;
        let mut count: u8 = 0;
        while pos_end < position {
            sym = self.runs[i].0;
            count = self.runs[i].1;
            total_counts[sym as usize] += count as u64;
            pos_end += count as u64;
            i += 1;
        }

        //we potentially went past the target, so subtract those back out
        total_counts[sym as usize] -= pos_end - position;

        if sym == symbol {
            //this can only happen if i != 0, so increment i-1
            self.increment_run(i-1);
        } else if position < pos_end {
            //our insert is in the middle of the run and different, we need to split the run
            let counts_after: u8 = (pos_end - position) as u8;
            let counts_before: u8 = count - counts_after;

            self.runs[i-1].1 = counts_before;

            //the lines below are doing this two insert operation
            //self.runs.insert(i, (symbol, 1));
            //self.runs.insert(i+1, (sym, counts_after));

            //this method only shifts the values once, so more efficient than double insert
            self.runs.push((0, 0));
            self.runs.push((0, 0));
            for j in (i+2..self.runs.len()).rev() {
                self.runs[j] = self.runs[j-2];
            }
            self.runs[i] = (symbol, 1);
            self.runs[i+1] = (sym, counts_after);
        } else if i < self.runs.len() && self.runs[i].0 == symbol {
            //right on a position boundary AND the next run matches our symbol, so increment
            self.increment_run(i);
        } else {
            //right on a position boundary (or at the very end) and the next run does not match, so do an insert
            self.runs.insert(i, (symbol, 1));
        }
        
        //return the accumulated counts
        total_counts[symbol as usize]
    }

    #[inline]
    fn increment_run(&mut self, block_index: usize) {
        //increment
        self.runs[block_index].1 = self.runs[block_index].1.wrapping_add(1);
        if self.runs[block_index].1 == 0 {
            //if it overflows, split into two halfs at the same position
            self.runs[block_index].1 = HALF_FULL;
            self.runs.insert(block_index+1, self.runs[block_index]);
        }
    }
    
    #[inline]
    pub fn split(&mut self) -> RLEBlock {
        //break right in the middle
        let breakpoint: usize = self.runs.len() / 2;
        let mut new_block_counts: [u64; VC_LEN] = [0; VC_LEN];
        let mut new_block_size: u64 = 0;
        
        //drain from the breakpoint, counting the symbols as we go
        let new_runs_data: ArrayVec<(u8, u8), CAPACITY_BUFFER> = self.runs.drain(breakpoint..)
            .map(
                |(sym, count)| {
                    new_block_counts[sym as usize] += count as u64;
                    (sym, count)
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
    }

    #[inline]
    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = vec![0; self.values_contained as usize];
        let mut ret_index: usize = 0;
        for &(curr_sym, count) in self.runs.iter() {
            ret[ret_index..ret_index+count as usize].fill(curr_sym);
            ret_index += count as usize;
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init() {
        let block: RLEBlock = Default::default();
        assert_eq!(block.runs.to_vec(), vec![]);
    }
    
    #[test]
    fn test_insert() {
        let mut block: RLEBlock = Default::default();
        assert_eq!(block.insert_and_count(0, 0), 0);
        //assert_eq!(block.run_symbols, vec![0]);
        //assert_eq!(block.run_counts, vec![1]);

        for _ in 0..256 {
            assert_eq!(block.insert_and_count(1, 1), 0);
        }
        //assert_eq!(block.run_symbols, vec![0, 1, 1]);
        //assert_eq!(block.run_counts, vec![1, 0, 1]);
        
        for _ in 0..256 {
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        //assert_eq!(block.run_symbols, vec![0, 0, 1, 1]);
        //assert_eq!(block.run_counts, vec![1, 1, 0, 1]);
        
        for _ in 0..255 {
            assert_eq!(block.insert_and_count(0, 2), 0);
        }
        //assert_eq!(block.run_symbols, vec![2, 0, 0, 1, 1]);
        //assert_eq!(block.run_counts, vec![255, 1, 1, 0, 1]);
        //assert_eq!(block.runs.to_vec(), vec![(2, 255), (0, 1), (0, 1), (1, 0), (1,1)]);
        assert_eq!(block.runs.to_vec(), vec![(2, 255), (0, 129), (0, 128), (1, 128), (1, 128)])
    }
    
    #[test]
    fn test_insert_splitting() {
        //middle split
        let mut block: RLEBlock = Default::default();
        for _ in 0..255 {
            assert_eq!(block.insert_and_count(0, 0), 0);
        }
        assert_eq!(block.insert_and_count(128, 1), 0);
        //assert_eq!(block.run_symbols, vec![0, 0, 1, 0, 0]);
        //assert_eq!(block.run_counts, vec![0, 128, 1, 0, 128]);
        assert_eq!(block.runs.to_vec(), vec![(0, 128), (1, 1), (0, 127)]);
    }
    
    #[test]
    fn test_increment_run() {
        //simple multiple block extension
        let mut block: RLEBlock = Default::default();
        block.insert_and_count(0, 0);
        for _ in 0..256 {
            block.increment_run(0);
        }
        //assert_eq!(block.run_symbols, vec![0, 0, 0]);
        //assert_eq!(block.run_counts, vec![0, 0, 1]);
        assert_eq!(block.runs.to_vec(), vec![(0, 129), (0, 128)]);
    }

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
            //run_symbols: vec![0, 0],
            //run_counts: vec![0, 128],
            runs,
            symbol_counts: [256, 0, 0, 0, 0, 0],
            values_contained: 256
        });
        
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(1, 1), (0, 128), (0, 128)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            //run_symbols: vec![1, 0, 0],
            //run_counts: vec![1, 0, 128],
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
            //run_symbols: vec![0, 1],
            //run_counts: vec![1, 1],
            runs,
            symbol_counts: [1, 1, 0, 0, 0, 0],
            values_contained: 2
        });
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(0, 127), (0, 128)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            //run_symbols: vec![0, 0],
            //run_counts: vec![0xFF, 0xFF],
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
            //run_symbols: vec![0, 0],
            //run_counts: vec![0xFF, 0xFF],
            runs,
            symbol_counts: [255, 0, 0, 0, 0, 0],
            values_contained: 255
        });
        let mut runs = ArrayVec::<(u8, u8), CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![(1, 1), (0, 1)]).unwrap();
        assert_eq!(new_block, RLEBlock {
            //run_symbols: vec![1, 0],
            //run_counts: vec![1, 1],
            runs,
            symbol_counts: [1, 1, 0, 0, 0, 0],
            values_contained: 2
        });
        
    }
    
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

            let count = block.insert_and_count(positions[i], inserted[i]);
            assert_eq!(block.to_vec(), data);
            assert_eq!(count, expected_count);
        }
    }
}
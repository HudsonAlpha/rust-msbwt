use arrayvec::ArrayVec;

/// This is the number of bytes used for data storage.
/// When it reaches this size, the block should be split to avoid issues.
pub const MAX_BLOCK_SIZE: usize = 60; // best for current B+
const CAPACITY_BUFFER: usize = MAX_BLOCK_SIZE+4;

/// The number of characters in our alphabet
pub const VC_LEN: usize = 2;      //0 or 1
/// The number of bits for storing the symbol
const SYMBOL_BITS: usize = 1;
/// The number of bit for storing quantity in a byte
const NUMBER_BITS: usize = 7; //8-letterBits
/// Contains the right-shifted number mask
const COUNT_MASK: u8 = 0x7F;

const SYMBOL_MASK: u8 = 0x1;
/// Multiplier for multi-byte runs
//pub const NUM_POWER: usize = 128;  //2**numberBits
const SINGLE_COUNT: u8 = 0x2;

const HALF_FULL: u8 = 0x80; // this is 64 shifted left 1

/// A run-length encoded block of data implemented with an ArrayVec
#[derive(Clone,Debug,PartialEq)]
pub struct BinaryBlock {
    runs: ArrayVec<u8, CAPACITY_BUFFER>,
    symbol_counts: [u64; VC_LEN],
    values_contained: u64
}

impl Default for BinaryBlock {
    /// Default function
    #[inline]
    fn default() -> Self {
        let runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        Self {
            runs,
            symbol_counts: [0; VC_LEN],
            values_contained: 0
        }
    }
}

#[inline] 
fn calc_run_bytes(count: u64) -> usize {
    if count < 64 {
        1
    } else if count < 2048 {
        2
    } else if count == 2048 {
        4
    } else {
        panic!("no impl > 2048");
    }
}

#[inline]
fn encode_run(run: &mut [u8], symbol: u8, count: u64) -> usize {
    //run[0] = (symbol & SYMBOL_MASK) | (count << SYMBOL_BITS) as u8;
    //1
    if count < 64 {
        //1 bit for symbol, 1 for control, so 6 for length
        run[0] = symbol | (count << SYMBOL_BITS) as u8;
        1
    } else if count < 2048 {
        //1 bit for symbol, 2*2 for control, so 11 for length
        run[0] = 0xC0 | symbol | (count >> 6 << SYMBOL_BITS) as u8;
        run[1] = 0x80 | (count & 0x3F) as u8;
        2
    }
    else if count == 2048 {
        println!("split");
        //1 bit for symbol, 2*2 for control, so 11 for length
        run[0] = 0xC0 | symbol | (1024 >> 6 << SYMBOL_BITS) as u8;
        run[1] = 0x80 | (1024 & 0x3F) as u8;
        run[2] = run[0];
        run[3] = run[1];
        4
    } else {
        panic!("no impl > 2048");
    }
}

#[inline]
fn decode_run(run: &[u8]) -> (u8, u64, usize) {
    if run[0] & 0x80 == 0 {
        (run[0] & SYMBOL_MASK, (run[0] >> SYMBOL_BITS) as u64, 1)
    } else {
        (run[0] & SYMBOL_MASK, ((run[0] & 0x3E) as u64) << 5 | (run[1] & 0x3F) as u64, 2)
    }
    //(run[0] & SYMBOL_MASK, (run[0] >> SYMBOL_BITS) as u64, 1)
}

impl BinaryBlock {
    #[inline]
    pub fn get_symbol_counts(&self) -> [u64; VC_LEN] {
        self.symbol_counts
    }

    #[inline]
    pub fn get_values_contained(&self) -> u64 {
        self.values_contained
    }

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
        let mut total_counts: [u64; VC_LEN+1] = [0; VC_LEN+1];
        
        //now iterate
        let mut i: usize = 0;
        let mut run: (u8, u64, usize) = (VC_LEN as u8, 0, 0);
        while pos_end < position {
            //get the symbol and the added counts
            run = decode_run(&self.runs[i..]);

            //add the count to the total so far
            total_counts[run.0 as usize] += run.1;
            pos_end += run.1;
            i += run.2 as usize;
        }
        //we potentially went past the target, so subtract those back out
        total_counts[run.0 as usize] -= pos_end - position;
        
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
        
        //its more efficient to count all symbols and return the one we care about than do checks inside the loop
        let mut pos_start: u64 = 0;
        let mut pos_end: u64 = 0;
        let mut total_counts: [u64; VC_LEN+1] = [0; VC_LEN+1];
        
        //now iterate
        let mut i: usize = 0;
        let mut run: (u8, u64, usize) = (VC_LEN as u8, 0, 0);
        while pos_end < position {
            //get the symbol and the added counts
            run = decode_run(&self.runs[i..]);

            //add the count to the total so far
            total_counts[run.0 as usize] += run.1;
            pos_end += run.1;
            i += run.2 as usize;
        }

        //we potentially went past the target, so subtract those back out
        total_counts[run.0 as usize] -= pos_end - position;

        //first, find the end of the run we're currently in
        if run.0 == symbol {
            //the run we're looking at is the same as what we're inserting
            self.increment_run(i-run.2);
        } else if position < pos_end {
            //our insert is in the middle of the run and different, we need to split the run
            let counts_after: u64 = pos_end - position;
            let counts_before: u64 = run.1 - counts_after;
            encode_run(&mut self.runs[i-run.2..], run.0, counts_before);

            let new_length = 1+calc_run_bytes(counts_before)+calc_run_bytes(counts_after) - run.2;

            //this method only shifts the values once, so more efficient than double insert
            self.runs.try_extend_from_slice(&vec![0; new_length]).unwrap();
            for j in (i+new_length..self.runs.len()).rev() {
                self.runs[j] = self.runs[j-new_length];
            }
            encode_run(&mut self.runs[i-run.2+calc_run_bytes(counts_before)..], symbol, 1);
            encode_run(&mut self.runs[i-run.2+calc_run_bytes(counts_before)+1..], run.0, counts_after);
        } else if i < self.runs.len() && (self.runs[i] & SYMBOL_MASK) == symbol {
            //right on a run boundary AND the next run matches our symbol, so increment the next run
            self.increment_run(i);
        } else {
            //right on a position boundary (or at the very end) and the next run does not match, so do an insert
            //self.runs.insert(i, symbol | SINGLE_COUNT);
            self.runs.insert(i, 0);
            encode_run(&mut self.runs[i..], symbol, 1);
        }

        //return the accumulated counts for the symbol we care about
        total_counts[symbol as usize]
    }

    #[inline]
    fn increment_run(&mut self, run_index: usize) {
        /*
        self.runs[run_index] = self.runs[run_index].wrapping_add(2);
        if self.runs[run_index] & 0xFE == 0 {
            self.runs[run_index] += HALF_FULL;
            self.runs.insert(run_index+1, self.runs[run_index]);
        }
        */
        let curr_run = decode_run(&self.runs[run_index..]);
        let curr_bytes = curr_run.2;
        let next_bytes = calc_run_bytes(curr_run.1+1);
        
        for _ in 0..next_bytes-curr_bytes {
            self.runs.insert(run_index+1, 0);
        }

        encode_run(&mut self.runs[run_index..], curr_run.0, curr_run.1+1);
    }

    #[inline]
    pub fn split(&mut self) -> BinaryBlock {
        //break right in the middle
        let mut breakpoint: usize = self.runs.len() / 2;
        while (self.runs[breakpoint] & 0xC0) == 0x80 {
            breakpoint -= 1;
        }
        let mut new_block_counts: [u64; VC_LEN] = [0; VC_LEN];
        let mut new_block_size: u64 = 0;
        
        //drain from the breakpoint, counting the symbols as we go
        let new_runs_data: ArrayVec<u8, CAPACITY_BUFFER> = self.runs.drain(breakpoint..).collect();
        let mut i = 0;
        while i < new_runs_data.len() {
            let run = decode_run(&new_runs_data[i..]);
            new_block_counts[run.0 as usize] += run.1;
            i += run.2 as usize;
        }

        //count up the total also
        new_block_size += new_block_counts.iter().sum::<u64>();

        //build the block
        let ret: BinaryBlock = BinaryBlock {
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

    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = Vec::<u8>::with_capacity(self.values_contained as usize);
        
        //now iterate
        //let mut sym: u8;
        //let mut count: u8;
        let mut i = 0;
        //for &run_byte in self.runs.iter() {
        while i < self.runs.len() {
            //get the symbol and the added counts
            //sym = run_byte & SYMBOL_MASK;
            //count = run_byte >> SYMBOL_BITS;
            let run = decode_run(&self.runs[i..]);
            
            //add the run symbol to the vec
            for _ in 0..run.1 {
                ret.push(run.0 as u8);
            }
            i += run.2 as usize;
        }
        
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init() {
        let block: BinaryBlock = Default::default();
        assert_eq!(block.runs.to_vec(), Vec::<u8>::new());
        assert_eq!(block.symbol_counts, [0, 0]);
        assert_eq!(block.values_contained, 0);
        
        //TODO add these?
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<u8>>(), vec![]);
        for symbol in 0..VC_LEN {
            assert_eq!(0, block.count(0, symbol as u8));
        }
    }
    
    
    #[test]
    fn test_insert() {
        let mut block: BinaryBlock = Default::default();

        //insert a single 0
        assert_eq!(block.count(0, 0), 0);
        assert_eq!(block.insert_and_count(0, 0), 0);
        let correct_data = vec![2];
        assert_eq!(block.runs.to_vec(), correct_data);

        //insert 256 1s at index 1
        for _ in 0..31 {
            assert_eq!(block.count(1, 1), 0);
            assert_eq!(block.insert_and_count(1, 1), 0);
        }
        let correct_data = vec![2, 0x1 | (31 << 1)];
        assert_eq!(block.runs.to_vec(), correct_data);

        //insert 256 0s at index 1
        for _ in 0..31 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        let correct_data = vec![0xC0 | 0, 0x80 | 32, 0x1 | (31 << 1)];
        assert_eq!(block.runs.to_vec(), correct_data);
        
        for _ in 0..32 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        let correct_data = vec![0xC0 | (64 >> 6) << 1, 0x80 | 0, 0x1 | (31 << 1)];
        assert_eq!(block.runs.to_vec(), correct_data);
        

        //insert 255 0s at index 1
        for _ in 0..255 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        //let correct_data = vec![0, 8, 1, 5];
        //assert_eq!(block.runs.to_vec(), correct_data);
        
        //TODO
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }
    
    #[test]
    fn test_insert_splitting() {
        //middle split
        let mut block: BinaryBlock = Default::default();
        for _ in 0..255 {
            assert_eq!(block.count(0, 0), 0);
            assert_eq!(block.insert_and_count(0, 0), 0);
        }
        assert_eq!(block.count(128, 1), 0);
        assert_eq!(block.insert_and_count(128, 1), 0);
        // we should have 128 0s, 1 1, and then 127 0s

        //let correct_data = vec![0, 2, 3, 0xFE];
        //assert_eq!(block.runs.to_vec(), correct_data);
        //TODO
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }
    
    #[test]
    fn test_increment_run() {
        //simple multiple block extension
        let mut block: BinaryBlock = Default::default();
        block.insert_and_count(0, 0);
        for _ in 0..127 {
            block.increment_run(0);
        }
        let correct_data = vec![0xC0 | (128 >> 6 << 1), 0x80 | (128 & 0x3F)];
        assert_eq!(block.runs.to_vec(), correct_data);

        //now break and increment that run
        block.insert_and_count(1, 1);
        let correct_data = vec![2, 3, 0xC0 | (127 >> 6 << 1), 0x80 | (127 & 0x3F)];
        assert_eq!(block.runs.to_vec(), correct_data);

        for _ in 0..127 {
            block.increment_run(1);
        }
        //one 0, 128 1s, 127 0s
        let correct_data = vec![2,
            0xC0 | (128 >> 6 << 1) | 0x1, 0x80 | (128 & 0x3F),
            0xC0 | (127 >> 6 << 1), 0x80 | (127 & 0x3F)
        ];
        assert_eq!(block.runs.to_vec(), correct_data);

        //TODO
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }

    
    #[test]
    fn test_block_splits() {
        //middle split
        let mut block: BinaryBlock = Default::default();
        for _ in 0..512 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(256, 1);

        let new_block: BinaryBlock = block.split();
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![0xC0 | (256_u64 >> 6 << 1) as u8, 0x80 | (256_u64 & 0x3F) as u8]).unwrap();

        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [256, 0],
            values_contained: 256
        });
        /*
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![3, 128, 128, 128, 128]).unwrap();
        assert_eq!(new_block, BinaryBlock {
            runs,
            symbol_counts: [256, 1],
            values_contained: 257
        });

        //early split
        let mut block: BinaryBlock = Default::default();
        for _ in 0..256 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(1, 1);
        let new_block: BinaryBlock = block.split();
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![2, 3, 126]).unwrap();
        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [64, 1],
            values_contained: 65
        });
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![128, 128, 128]).unwrap();
        assert_eq!(new_block, BinaryBlock {
            runs,
            symbol_counts: [192, 0],
            values_contained: 192
        });
        
        //late split
        let mut block: BinaryBlock = Default::default();
        for _ in 0..256 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(255, 1);
        let new_block: BinaryBlock = block.split();
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![128, 128, 128]).unwrap();
        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [192, 0],
            values_contained: 192
        });
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![126, 3, 2]).unwrap();
        assert_eq!(new_block, BinaryBlock {
            runs,
            symbol_counts: [64, 1],
            values_contained: 65
        });
        */
    }
    
    #[test]
    fn test_to_vec() {
        //big run test
        let mut block: BinaryBlock = Default::default();
        for _ in 0..1024 {
            block.insert_and_count(0, 1);
        }
        assert_eq!(block.to_vec(), vec![1; 1024]);

        //semi-random collection test
        let mut block: BinaryBlock = Default::default();
        let data: Vec<u8> =            vec![0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0];
        let expected_count: Vec<u64> = vec![0, 0, 1, 2, 1, 2, 3, 3, 4, 4, 5, 6, 5];
        for (i, v) in data.iter().enumerate() {
            assert_eq!(block.count(i as u64, *v), expected_count[i]);
            assert_eq!(block.insert_and_count(i as u64, *v), expected_count[i]);
        }
        assert_eq!(block.to_vec(), data);
    }
}
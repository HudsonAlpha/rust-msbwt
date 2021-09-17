use arrayvec::ArrayVec;

/// This is the number of bytes used for data storage.
/// When it reaches this size, the block should be split to avoid issues.
pub const MAX_BLOCK_SIZE: usize = 63; // best for current B+
const CAPACITY_BUFFER: usize = MAX_BLOCK_SIZE+2;

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

impl BinaryBlock {
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
        let mut prev_sym: u8 = VC_LEN as u8;
        let mut sym: u8 = VC_LEN as u8;
        let mut curr_shift = 0;
        let mut count: u64;
        while pos_end < position {
            //get the symbol and the added counts
            sym = self.runs[i] & SYMBOL_MASK;
            if sym == prev_sym {
                count = ((self.runs[i] as u64) >> SYMBOL_BITS) << curr_shift;
                curr_shift += NUMBER_BITS;
            } else {
                count = (self.runs[i] as u64) >> SYMBOL_BITS;
                curr_shift = NUMBER_BITS;
            }

            //add the count to the total so far
            total_counts[sym as usize] += count as u64;
            pos_end += count as u64;
            prev_sym = sym;
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
        
        //its more efficient to count all symbols and return the one we care about than do checks inside the loop
        let mut pos_start: u64 = 0;
        let mut pos_end: u64 = 0;
        let mut total_counts: [u64; VC_LEN+1] = [0; VC_LEN+1];
        
        //now iterate
        let mut i: usize = 0;
        let mut prev_sym: u8 = VC_LEN as u8;
        let mut sym: u8 = VC_LEN as u8;
        let mut curr_shift = 0;
        let mut count: u64;
        let mut run_start = 0;
        while pos_end < position {
            //get the symbol and the added counts
            sym = self.runs[i] & SYMBOL_MASK;
            if sym == prev_sym {
                count = ((self.runs[i] as u64) >> SYMBOL_BITS) << curr_shift;
                curr_shift += NUMBER_BITS;
            } else {
                count = (self.runs[i] as u64) >> SYMBOL_BITS;
                curr_shift = NUMBER_BITS;
                run_start = i;
                pos_start = pos_end;
            }

            //add the count to the total so far
            total_counts[sym as usize] += count as u64;
            pos_end += count as u64;
            prev_sym = sym;
            i += 1;
        }

        //we potentially went past the target, so subtract those back out
        total_counts[sym as usize] -= pos_end - position;

        //first, find the end of the run we're currently in
        while i < self.runs.len() && (self.runs[i] & 0x1) == sym {
            pos_end += ((self.runs[i] as u64) >> SYMBOL_BITS) << curr_shift;
            curr_shift += NUMBER_BITS;
            i += 1;
        }
        let run_end = i;
        
        if sym == symbol {
            //this can only happen if i != 0, so increment i-1
            self.increment_run(run_start);
        } else if position < pos_end {
            //our insert is in the middle of the run and different, we need to split the run
            let mut counts_after: u64 = pos_end - position;
            let mut counts_before: u64 = position - pos_start;

            //build up the before run, the insert single character, and the after run
            let mut splice_in: Vec<u8> = Vec::<u8>::with_capacity(10);
            while counts_before > 0 {
                let new_byte = sym | (((counts_before & COUNT_MASK as u64) as u8) << 1);
                splice_in.push(new_byte);
                counts_before >>= NUMBER_BITS;
            }
            splice_in.push(symbol | 0x2);
            while counts_after > 0 {
                let new_byte = sym | (((counts_after & COUNT_MASK as u64) as u8) << 1);
                splice_in.push(new_byte);
                counts_after >>= NUMBER_BITS;
            }

            //calculate the new size of the array relative to the old size
            let byte_shift = splice_in.len() - (run_end - run_start);

            //add buffer at the end
            for _ in 0..byte_shift {
                self.runs.push(0);
            }
            //now shift values into the buffer
            for j in (run_end+byte_shift..self.runs.len()).rev() {
                self.runs[j] = self.runs[j-byte_shift];
            }

            //now splice in the new values in their space
            for (index, &value) in splice_in.iter().enumerate() {
                self.runs[run_start+index] = value;
            }

        } else if i < self.runs.len() {
            //right on a position boundary AND the next run must match our symbol, so increment
            self.increment_run(i);
        } else {
            //right on a position boundary (or at the very end) and the next run does not match, so do an insert
            self.runs.insert(i, symbol | 0x2);
        }

        //return the accumulated counts for the symbol we care about
        total_counts[symbol as usize]
    }

    #[inline]
    fn increment_run(&mut self, run_index: usize) {
        let mut current_index: usize = run_index;
        let symbol: u8 = self.runs[run_index] & 0x1;

        //zero out while we have to
        let mut to_increment: bool = true;
        while current_index < self.runs.len() && 
            (self.runs[current_index] & SYMBOL_MASK) == symbol &&
            to_increment
        {
            let overflow_result = self.runs[current_index].overflowing_add(SINGLE_COUNT);
            self.runs[current_index] = overflow_result.0;
            to_increment = overflow_result.1;
            current_index += 1;
        }

        //we need to extend the run into the next byte
        if to_increment {
            self.runs.insert(current_index, symbol+SINGLE_COUNT);
        }
    }

    #[inline]
    pub fn split(&mut self) -> BinaryBlock {
        let midpoint: usize = self.runs.len() / 2;
        let mut closest_down: usize = midpoint;
        while closest_down > 0 && (self.runs[closest_down-1] & SYMBOL_MASK) == (self.runs[closest_down] & SYMBOL_MASK) {
            closest_down -= 1;
        }

        let mut closest_up: usize = midpoint;
        while closest_up < self.runs.len() && (self.runs[closest_up-1] & SYMBOL_MASK) == (self.runs[closest_up] & SYMBOL_MASK) {
            closest_up += 1;
        }

        if closest_down == 0 && closest_up == self.runs.len() {
            panic!("single run, block cannot be split");
        }

        let breakpoint: usize = if midpoint - closest_down < closest_up - midpoint {
            //break using down 
            closest_down
        } else {
            //break using up
            closest_up
        };

        let mut new_block_counts: [u64; VC_LEN] = [0; VC_LEN];
        let mut new_block_size = 0;
        let mut curr_sym: u8 = 0;
        let mut curr_count: u64 = 0;
        let mut curr_shift: usize = 0;
        for c in self.runs[breakpoint..].iter() {
            if curr_sym == (c & SYMBOL_MASK) {
                //same symbol in run, use multiplier and increase it
                curr_count += ((c >> SYMBOL_BITS) as u64) << curr_shift;
                curr_shift += NUMBER_BITS;
            } else {
                new_block_counts[curr_sym as usize] += curr_count;
                new_block_size += curr_count;

                //reset to new symbol
                curr_sym = c & SYMBOL_MASK;
                curr_count = (c >> SYMBOL_BITS) as u64;
                curr_shift = NUMBER_BITS;
            }
        }
        new_block_counts[curr_sym as usize] += curr_count;
        new_block_size += curr_count;

        //let mut new_block_data: Vec<u8> = Vec::<u8>::with_capacity(CAPACITY_BUFFER);
        let mut new_runs_data = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        new_runs_data.try_extend_from_slice(&self.runs[breakpoint..]).unwrap();
        self.runs.truncate(breakpoint);

        //build the block
        let ret: BinaryBlock = BinaryBlock {
            runs: new_runs_data,
            symbol_counts: new_block_counts,
            values_contained: new_block_size
        };

        //now update internals before returning
        for i in 0..VC_LEN {
            self.symbol_counts[i] -= new_block_counts[i];
        }
        self.values_contained -= new_block_size;
        ret
    }

    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = Vec::<u8>::with_capacity(self.values_contained as usize);
        
        //now iterate
        let mut prev_sym: u8 = VC_LEN as u8;
        let mut sym: u8;
        let mut curr_shift = 0;
        let mut count: u64;
        for &run_byte in self.runs.iter() {
            //get the symbol and the added counts
            sym = run_byte & SYMBOL_MASK;
            if sym == prev_sym {
                count = ((run_byte as u64) >> SYMBOL_BITS) << curr_shift;
                curr_shift += NUMBER_BITS;
            } else {
                count = (run_byte as u64) >> SYMBOL_BITS;
                curr_shift = NUMBER_BITS;
            }

            //add the run symbol to the vec
            prev_sym = sym;
            for _ in 0..count {
                ret.push(sym);
            }
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
        for _ in 0..256 {
            assert_eq!(block.count(1, 1), 0);
            assert_eq!(block.insert_and_count(1, 1), 0);
        }
        let correct_data = vec![2, 1, 5];
        assert_eq!(block.runs.to_vec(), correct_data);

        //insert 256 0s at index 1
        for _ in 0..256 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        let correct_data = vec![2, 4, 1, 5];
        assert_eq!(block.runs.to_vec(), correct_data);

        //insert 255 0s at index 1
        for _ in 0..255 {
            assert_eq!(block.count(1, 0), 1);
            assert_eq!(block.insert_and_count(1, 0), 1);
        }
        let correct_data = vec![0, 8, 1, 5];
        assert_eq!(block.runs.to_vec(), correct_data);
        
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

        let correct_data = vec![0, 2, 3, 0xFE];
        assert_eq!(block.runs.to_vec(), correct_data);
        //TODO
        //assert_eq!(block.raw_iter().cloned().collect::<Vec<(u8, u8)>>(), correct_data);
    }
    
    #[test]
    fn test_increment_run() {
        //simple multiple block extension
        let mut block: BinaryBlock = Default::default();
        block.insert_and_count(0, 0);
        for _ in 0..128_usize.pow(3) {
            block.increment_run(0);
        }
        let correct_data = vec![2, 0, 0, 2];
        assert_eq!(block.runs.to_vec(), correct_data);

        //now break and increment that run
        block.insert_and_count(1, 1);
        let correct_data = vec![2, 3, 0, 0, 0, 2];
        assert_eq!(block.runs.to_vec(), correct_data);

        for _ in 0..128_usize.pow(3) {
            block.increment_run(1);
        }
        let correct_data = vec![2, 3, 1, 1, 3, 0, 0, 0, 2];
        assert_eq!(block.runs.to_vec(), correct_data);

        //one final increment on the original big run just to make sure things are cool
        block.increment_run(5);
        let correct_data = vec![2, 3, 1, 1, 3, 2, 0, 0, 2];
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
        runs.try_extend_from_slice(&vec![0, 4]).unwrap();

        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [256, 0],
            values_contained: 256
        });
        
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![3, 0, 4]).unwrap();
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
        runs.try_extend_from_slice(&vec![2, 3]).unwrap();
        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [1, 1],
            values_contained: 2
        });
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![0xFE, 2]).unwrap();
        assert_eq!(new_block, BinaryBlock {
            runs,
            symbol_counts: [255, 0],
            values_contained: 255
        });
        
        //late split
        let mut block: BinaryBlock = Default::default();
        for _ in 0..256 {
            block.insert_and_count(0, 0);
        }
        block.insert_and_count(255, 1);
        let new_block: BinaryBlock = block.split();
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![0xFE, 2]).unwrap();
        assert_eq!(block, BinaryBlock {
            runs,
            symbol_counts: [255, 0],
            values_contained: 255
        });
        let mut runs = ArrayVec::<u8, CAPACITY_BUFFER>::new();
        runs.try_extend_from_slice(&vec![3, 2]).unwrap();
        assert_eq!(new_block, BinaryBlock {
            runs,
            symbol_counts: [1, 1],
            values_contained: 2
        });
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
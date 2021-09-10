
use crate::run_block_av_flat::VC_LEN;
use crate::rle_bplus_tree::RLEBPlusTree;
use crate::string_util::convert_stoi;

const INITIAL_QUERY: usize = 10;
const COST_FACTOR: f64 = 0.000001;

pub struct DynamicBWT {
    string_count: usize,
    symbol_count: usize,
    tree_bwt: RLEBPlusTree,
    //tree_bwt: RLEBinaryTree,
    //tree_bwt: RLESkipList,
    total_counts: [usize; VC_LEN],
    offsets: [usize; VC_LEN],
    sort_query_len: f64,
    short_circuits: [usize; 3],
}

impl Default for DynamicBWT {
    fn default() -> Self {
        DynamicBWT {
            string_count: 0,
            symbol_count: 0,
            tree_bwt: Default::default(),
            total_counts: [0; VC_LEN],
            offsets: [0; VC_LEN],
            sort_query_len: INITIAL_QUERY as f64,
            short_circuits: [0; 3]
        }
    }
}

impl DynamicBWT {
    #[inline]
    pub fn get_total_counts(&self) -> [usize; VC_LEN] {
        self.total_counts
    }

    #[inline]
    pub fn get_height(&self) -> usize {
        self.tree_bwt.get_height()
    }

    #[inline]
    pub fn get_node_count(&self) -> usize {
        self.tree_bwt.get_node_count()
    }
    ///*
    //this is a much simpler version, but is also slower
    #[inline]
    pub fn insert_string(&mut self, val: &str, sorted: bool) {
        let int_form: Vec<u8> = convert_stoi(val);

        //initial position is the total number of string
        let mut next_insert;
        
        if sorted {
            let mut start_index = 0;
            next_insert = self.symbol_count;

            //attempt a short circuit
            let query_len = std::cmp::min(self.sort_query_len as usize, int_form.len());
            for pred_symbol in int_form[..query_len].iter().rev() {
                start_index = self.tree_bwt.count(start_index, *pred_symbol)+self.offsets[*pred_symbol as usize];
                next_insert = self.tree_bwt.count(next_insert, *pred_symbol)+self.offsets[*pred_symbol as usize];
            }
            start_index = self.tree_bwt.count(start_index, 0);
            next_insert = self.tree_bwt.count(next_insert, 0);

            if start_index != next_insert {
                let original_ni: usize = next_insert;

                //short circuit failed
                for pred_symbol in int_form.iter().rev() {
                    next_insert = self.tree_bwt.count(next_insert, *pred_symbol)+self.offsets[*pred_symbol as usize];
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
            self.total_counts[*pred_symbol as usize] += 1;
            for i in (symbol+1) as usize..VC_LEN {
                self.offsets[i] += 1;
            }

            //after any adjustments add in the new offset for the symbol we're currently at
            next_insert += self.offsets[*pred_symbol as usize];
            symbol = *pred_symbol;
        }

        //one final insert for the $
        self.tree_bwt.insert_and_count(next_insert, 0);
        self.total_counts[0] += 1;
        for i in (symbol+1) as usize..VC_LEN {
            self.offsets[i] += 1;
        }

        self.symbol_count += int_form.len()+1;
        self.string_count += 1;

        if self.string_count % 10000 == 0 {
            println!("{} {:?} {}", self.string_count, self.short_circuits, self.sort_query_len);
            self.short_circuits = [0; 3];
        }
    }

    #[inline]
    pub fn to_vec(&self) -> Vec<u8> {
        self.tree_bwt.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::string_util::STRING_TO_INT;

    fn naive_bwt(inputs: &Vec<&str>) -> Vec<u8> {
        let mut rotations: Vec<String> = vec![];
        for s in inputs.iter() {
            let dollar_string = s.to_string()+&"$".to_string();
            for l in 0..dollar_string.len() {
                rotations.push(dollar_string[l..].to_string()+&dollar_string[..l]);
            }
        }
        rotations.sort();
        let mut ret: Vec<u8> = Vec::<u8>::with_capacity(rotations.len());
        for r in rotations.iter() {
            ret.push(STRING_TO_INT[r.as_bytes()[r.len()-1] as usize]);
        }
        ret
    }

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
        let bwt = naive_bwt(&data);
        
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
        let bwt = naive_bwt(&data);

        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), bwt);

        //getting smaller in order
        let data: Vec<&str> = vec!["AAAAA", "AAAA", "AAA", "AA", "A"];
        let bwt = naive_bwt(&data);

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
        let naive = naive_bwt(&data);

        //insert the strings in the now sorted order
        let mut ubwt: DynamicBWT = Default::default();
        for s in data.iter() {
            ubwt.insert_string(s, true);
        }
        assert_eq!(ubwt.to_vec(), naive);
    }
}
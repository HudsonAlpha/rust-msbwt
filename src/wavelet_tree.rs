
use crate::binary_bplus_tree::BinaryBPlusTree;

const VC_LEN: usize = 6;

const BINARY_TRANSLATIONS: [[u8; 3]; VC_LEN] = [
    [0, 0, 0], //$
    [0, 1, 0], //A
    [0, 1, 1], //C
    [1, 0, 2], //G
    [0, 0, 1], //N
    [1, 1, 2]  //T
];

const BINARY_TO_SYMBOL: [u8; 8] = [
    0, // $ = 000
    4, // N = 001
    1, //A = 010
    2, //C = 011
    //the last bit isn't used for these, so we set both
    3, //G = 10
    3, //G = 10 
    5, //T = 11
    5  //T = 11
];

const TREE_ORDER: [[usize; 3]; VC_LEN] = [
    [0, 1, 3], //$
    [0, 1, 4], //A
    [0, 1, 4], //C
    [0, 2, VC_LEN], //G
    [0, 1, 3], //N
    [0, 2, VC_LEN]  //T
];

/// TODO
pub struct WaveletTree {
    /// The actual B+ tree structure
    //         0 
    //       /    \
    //      1     2
    //     / \   /  \
    //     3  4  G   T
    //    / \ / \
    //    $ N A C
    wavelets: [BinaryBPlusTree; VC_LEN-1]
}

impl Default for WaveletTree {
    /// TODO
    fn default() -> Self {
        WaveletTree {
            wavelets: Default::default()
        }
    }
}

impl WaveletTree {
    pub fn get_height(&self) -> usize {
        self.wavelets.iter().map(|w| w.get_height()).max().unwrap()
    }

    pub fn get_node_count(&self) -> usize {
        self.wavelets.iter().map(|w| w.get_node_count()).sum()
    }

    pub fn count(&self, index: u64, value: u8) -> u64 {
        let query_order = TREE_ORDER[value as usize];
        let bit_order = BINARY_TRANSLATIONS[value as usize];
        let mut curr_index = index;

        for (i, &qo) in query_order.iter().enumerate() {
            if qo == VC_LEN {
                //we're done, this is a sentinel
            } else {
                curr_index = self.wavelets[qo].count(curr_index, bit_order[i]);
            }
        }
        curr_index
    }

    pub fn insert_and_count(&mut self, index: u64, value: u8) -> u64 {
        let query_order = TREE_ORDER[value as usize];
        let bit_order = BINARY_TRANSLATIONS[value as usize];
        let mut curr_index = index;

        for (i, &qo) in query_order.iter().enumerate() {
            if qo == VC_LEN {
                //we're done, this is a sentinel
            } else {
                curr_index = self.wavelets[qo].insert_and_count(curr_index, bit_order[i]);
            }
        }
        curr_index
    }

    //TODO this is pretty inefficient because of the vec building
    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = Default::default();
        let mut vec_forms: Vec<Vec<u8>> = Default::default();
        for wl in self.wavelets.iter() {
            vec_forms.push(wl.to_vec());
        }

        //get iterators over all the wavelets
        let mut iter_forms: Vec<_> = Default::default();
        for wl in vec_forms.iter() {
            iter_forms.push(wl.iter());
        }

        //now go through all in b1
        while let Some(&b1) = iter_forms[0].next() {
            //b2 depends on b1
            let b2: u8 = *(iter_forms[1+b1 as usize].next().unwrap());
            // b3 only occurs if b1 == 0, but then depends on b2
            let b3 = if b1 == 0 {
                *(iter_forms[3+b2 as usize].next().unwrap())
            } else {
                0
            };

            //assemble the full binary lookup, and then translate to our original symbol
            let binary_sym = ((b1 << 2) + (b2 << 1) + b3) as usize;
            let symbol = BINARY_TO_SYMBOL[binary_sym];
            ret.push(symbol);
        }

        ret
    }

    pub fn run_iter(&self) -> WaveletTreeRunIterator<'_> {
        //let next_child_index = self.next_child[0];
        //let current_block_iter = self.data_arena[0].raw_iter();
        WaveletTreeRunIterator {
            tree: self,
            //next_child_index,
            //current_block_iter,
            //next_sym: 0,
            //next_count: 0
        }
    }
}

impl<'a> IntoIterator for &'a WaveletTree {
    type Item = u8;
    type IntoIter = WaveletTreeIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        let vec_iter = self.to_vec().into_iter();
        WaveletTreeIterator {
            tree: self,
            vec_iter
        }
    }
}

//TODO: make this more efficient
pub struct WaveletTreeIterator<'a> {
    tree: &'a WaveletTree,
    vec_iter: std::vec::IntoIter<u8>
}

impl<'a> Iterator for WaveletTreeIterator<'a> {
    type Item = u8;
    /// Will return the next symbol contained by the compressed B+ tree data
    fn next(&mut self) -> Option<u8> {
        match self.vec_iter.next() {
            //current block has data, so return it
            Some(x) => Some(x),
            None => {
                None
            }
        }
    }
}

pub struct WaveletTreeRunIterator<'a> {
    tree: &'a WaveletTree
}

impl<'a> Iterator for WaveletTreeRunIterator<'a> {
    type Item = (u8, u64);
    /// Will return the next run contained by the compressed B+ tree data
    fn next(&mut self) -> Option<(u8, u64)> {
        None
    }
}

#[cfg(test)]
mod tests {
    extern crate rand;
    use super::*;
    use rand::Rng;
    
    #[test]
    fn test_init() {
        let wavelet: WaveletTree = Default::default();
        assert_eq!(wavelet.to_vec(), Vec::<u8>::new());
    }

    #[test]
    fn test_simple_inserts() {
        let mut tree: WaveletTree = Default::default();
        let data: Vec<u8> =             vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
        let expected_counts: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
        for (i, v) in data.iter().enumerate() {
            let pre_count = tree.count(i as u64, *v);
            assert_eq!(pre_count, expected_counts[i]);

            let count = tree.insert_and_count(i as u64, *v);
            //println!("{} {:?}", i, tree.to_vec());
            assert_eq!(count, expected_counts[i]);
        }
        assert_eq!(tree.to_vec(), data);
        assert_eq!(tree.into_iter().collect::<Vec<u8>>(), data);
        /*
        //check the pairs since this is in-order
        let runs: Vec<(u8, u64)> = vec![(0, 1), (1, 3), (2, 1), (0, 1), (2, 1), (3, 1), (4, 1), (1, 3), (0, 1)];
        let full_pairing = tree.run_iter().collect::<Vec<(u8, u64)>>();
        assert_eq!(full_pairing, runs);
        */
    }

    #[test]
    fn test_10krandom_insert() {
        let mut tree: WaveletTree = Default::default();
        let mut data: Vec<u8> = vec![];
        let mut rng = rand::thread_rng();
        let mut inserted: Vec<u8> = vec![];
        let mut positions: Vec<usize> = vec![];
        for _ in 0..10000 {
            let symbol: u8 = rng.gen_range(0, 6);
            let position: usize = rng.gen_range(0, data.len()+1);
            inserted.push(symbol);
            positions.push(position);
            data.insert(position, symbol);
            println!("RANDOM_DATA: {:?}", data);
            println!("inserted: {:?}", inserted);
            println!("positions: {:?}", positions);
            
            let pre_count = tree.count(position as u64, symbol);
            let post_count = tree.insert_and_count(position as u64, symbol);
            assert_eq!(pre_count, post_count);
            
            assert_eq!(tree.to_vec(), data);
            
            //TODO
            assert_eq!(tree.into_iter().collect::<Vec<u8>>(), data);
        }
    }
}
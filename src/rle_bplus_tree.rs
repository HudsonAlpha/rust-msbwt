
extern crate arrayvec;

use arrayvec::ArrayVec;

use crate::run_block_av_flat::{VC_LEN, MAX_BLOCK_SIZE, RLEBlock}; // current best

const MAX_NODE_SIZE: usize = 64; // must be a power of 2
const NODE_MIDPOINT: usize = MAX_NODE_SIZE / 2;

/// A B+ tree structure that has special functionality to encode runs of the same symbol in a defined alphabet
#[derive(Clone)]
pub struct RLEBPlusTree {
    /// the current height of the tree
    height: usize,
    /// the collection of nodes in the tree, root is always at index 0
    nodes: Vec<RLEBPlusNode>,
    /// the collection of data blocks referenced by the tree nodes, these are the actual symbols
    data_arena: Vec<RLEBlock>,
    /// for the RLEBlocks, this says which block is next for an in-order traversal
    next_child: Vec<usize>
}

/// A node in the B+ tree
#[derive(Clone,Debug)]
struct RLEBPlusNode {
    /// stores if this is a leaf node or internal node
    is_leaf: bool,
    /// the parent node ID
    parent: usize,
    /// stores the total symbol counts for that child node and all its children
    total_counts: ArrayVec<u64, MAX_NODE_SIZE>,
    /// stores the individual symbol counts for that child node and all its children
    total_symbols: ArrayVec<[u64; VC_LEN], MAX_NODE_SIZE>,
    /// the indices of any child nodes; if it's a leaf this refers to RLEBlocks, otherwise other RLEBPlusNodes
    children: ArrayVec<usize, MAX_NODE_SIZE>
}

impl RLEBPlusNode {
    #[inline]
    fn count(&self, relative_index: u64, symbol: u8) -> (usize, u64, u64) {
        let mut total_count: u64 = 0;
        let bs_index: usize = self.total_counts.iter().position(|&v| {
            total_count += v;
            total_count >= relative_index
        }).unwrap();
        total_count -= self.total_counts[bs_index];

        let symbol_count: u64 = self.total_symbols[..bs_index].iter().fold(0, |acc, v| acc + v[symbol as usize]);
        (bs_index, symbol_count, total_count)
    }

    //for now, this is a clone except for two minor increments on the found node
    #[inline]
    fn insert_and_count(&mut self, relative_index: u64, symbol: u8) -> (usize, u64, u64) {
        let mut total_count: u64 = 0;
        let bs_index: usize = self.total_counts.iter().position(|&v| {
            total_count += v;
            total_count >= relative_index
        }).unwrap();
        total_count -= self.total_counts[bs_index];

        let symbol_count: u64 = self.total_symbols[..bs_index].iter().fold(0, |acc, v| acc + v[symbol as usize]);

        //increment total counts and symbols for this node
        self.total_counts[bs_index] += 1;
        self.total_symbols[bs_index][symbol as usize] += 1;

        (bs_index, symbol_count, total_count)
    }

    #[inline]
    fn build_indices(&mut self) {
        //place-holder
    }
}

impl Default for RLEBPlusTree {
    #[inline]
    fn default() -> Self {
        let children: Vec<RLEBlock> = vec![Default::default()];
        let mut child_indices: ArrayVec<usize, MAX_NODE_SIZE> = Default::default();
        child_indices.push(0);
        let mut total_counts = ArrayVec::<u64, MAX_NODE_SIZE>::new();
        total_counts.push(0);
        let mut total_symbols = ArrayVec::<[u64; VC_LEN], MAX_NODE_SIZE>::new();
        total_symbols.push([0; VC_LEN]);
        let root: RLEBPlusNode = RLEBPlusNode {
            is_leaf: true,
            parent: 0,
            total_counts,
            total_symbols,
            children: child_indices
        };

        RLEBPlusTree {
            height: 0,
            nodes: vec![root],
            data_arena: children,
            next_child: vec![0]
        }
    }
}

impl<'a> IntoIterator for &'a RLEBPlusTree {
    type Item = u8;
    type IntoIter = RLEBPlusTreeIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        let next_child_index = self.next_child[0];
        let current_block_iter = self.data_arena[0].to_vec().into_iter();
        RLEBPlusTreeIterator {
            tree: self,
            next_child_index,
            current_block_iter
        }
    }
}

impl RLEBPlusTree {
    /// Returns the current height of the B+ tree
    #[inline]
    pub fn get_height(&self) -> usize {
        self.height
    }

    /// Returns the total number of data nodes (i.e. leaf nodes) in the B+ tree
    #[inline]
    pub fn get_node_count(&self) -> usize {
        self.data_arena.len()
    }

    /// Performs a rank/count operation.
    /// For a given data index, this will count all occurences of symbol `value` up to that index.
    /// # Arguments
    /// * `index` - the index to count to 
    /// * `value` - the symbol to count
    /// # Examples
    /// ```rust
    /// use msbwt2::rle_bplus_tree::RLEBPlusTree;
    /// let mut tree: RLEBPlusTree = Default::default();
    /// let data: Vec<u8> =             vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
    /// let expected_counts: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
    /// for (i, v) in data.iter().enumerate() {
    ///     let pre_count = tree.count(i as u64, *v);
    ///     assert_eq!(pre_count, expected_counts[i]);
    ///     let count = tree.insert_and_count(i as u64, *v);
    ///     println!("{} {:?}", i, tree.to_vec());
    ///     assert_eq!(count, expected_counts[i]);
    /// }
    /// ```
    #[inline]
    pub fn count(&self, index: u64, value: u8) -> u64 {
        //start with root node 0
        let mut current_node_index: usize = 0;
        let mut current_node: &RLEBPlusNode = &self.nodes[current_node_index];
        let mut relative_index: u64 = index;
        let mut total_count: u64 = 0;
            
        //iterate downwards until we find a leaf node point to run blocks
        while !current_node.is_leaf {
            //linear search tended to be faster in practice
            let (bs_index, sc, tc) = current_node.count(relative_index, value);

            //add symbol count to the running total
            total_count += sc;
            //adjust our relative index by how much we went through
            relative_index -= tc;

            //go down to the indexing node
            current_node_index = current_node.children[bs_index];
            current_node = &self.nodes[current_node_index];
        }

        //we're in a leaf, still need to search for the data node
        let (bs_index, sc, tc) = current_node.count(relative_index, value);
        
        //add symbol count to the running total
        total_count += sc;
        //adjust our relative index by how much we went through
        relative_index -= tc;

        //get the data node from the search index
        let arena_index: usize = current_node.children[bs_index];
        
        //add in the final counts and return
        total_count += self.data_arena[arena_index].count(relative_index, value);
        total_count
    }

    /// Performs a rank/count operation while also inserting that symbol at the provided index.
    /// For a given data index, this will count all occurences of symbol `value` up to that index, and then insert an addition `value` at that index.
    /// # Arguments
    /// * `index` - the index to count to 
    /// * `value` - the symbol to count
    /// # Examples
    /// ```rust
    /// use msbwt2::rle_bplus_tree::RLEBPlusTree;
    /// let mut tree: RLEBPlusTree = Default::default();
    /// let data: Vec<u8> =             vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
    /// let expected_counts: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
    /// for (i, v) in data.iter().enumerate() {
    ///     let pre_count = tree.count(i as u64, *v);
    ///     assert_eq!(pre_count, expected_counts[i]);
    ///     let count = tree.insert_and_count(i as u64, *v);
    ///     println!("{} {:?}", i, tree.to_vec());
    ///     assert_eq!(count, expected_counts[i]);
    /// }
    /// ```
    #[inline]
    pub fn insert_and_count(&mut self, index: u64, value: u8) -> u64{
        //start with root node 0
        let mut current_node_index: usize = 0;
        let mut current_node = &mut self.nodes[current_node_index];
        let mut relative_index: u64 = index;
        let mut total_count: u64 = 0;

        //iterate downwards until we find a leaf node point to run blocks
        while !current_node.is_leaf {
            //linear search tended to be faster in practice
            let (bs_index, sc, tc) = current_node.insert_and_count(relative_index, value);

            //add symbol count to the running total
            total_count += sc;
            //adjust our relative index by how much we went through
            relative_index -= tc;

            //update current node
            current_node_index = current_node.children[bs_index];
            current_node = &mut self.nodes[current_node_index];
        }

        //we're in a leaf, get the leaf index and then the arena index for the data node
        let (bs_index, sc, tc) = current_node.insert_and_count(relative_index, value);
        
        //add symbol count to the running total
        total_count += sc;
        //adjust our relative index by how much we went through
        relative_index -= tc;
        
        //increment the total count by what comes back from the block insertion
        let arena_index: usize = current_node.children[bs_index];
        total_count += self.data_arena[arena_index].insert_and_count(relative_index, value);
        
        //handle any internal splitting that needs to occur due to nodes filling up
        self.split_internal_nodes(arena_index, current_node_index, bs_index);

        total_count
    }

    /// This is a helper function for splitting internal nodes to maintain the structure.
    /// It's currently only used in the insert function, but it's useful to have it separated out for readability.
    #[inline]
    fn split_internal_nodes(&mut self, arena_index: usize, node_index: usize, bs_index: usize) {
        //check if the block is big enough to get split
        let mut current_node_index = node_index;
        let current_node = &mut self.nodes[current_node_index];
        if self.data_arena[arena_index].block_len() >= MAX_BLOCK_SIZE {
            //first we need to actually split the block
            let new_block: RLEBlock = self.data_arena[arena_index].split();
            let new_arena_index = self.data_arena.len();
            self.data_arena.push(new_block);

            //correct the next child options
            self.next_child.push(self.next_child[arena_index]);
            self.next_child[arena_index] = new_arena_index;

            //insert into the before counts
            current_node.total_counts[bs_index] = self.data_arena[arena_index].get_values_contained();
            current_node.total_counts.insert(bs_index+1, self.data_arena[new_arena_index].get_values_contained());
            
            //insert into the before symbols also
            current_node.total_symbols[bs_index] = self.data_arena[arena_index].get_symbol_counts();
            current_node.total_symbols.insert(bs_index+1, self.data_arena[new_arena_index].get_symbol_counts());
            
            //insert the child into the list
            current_node.children.insert(bs_index+1, new_arena_index);

            current_node.build_indices();

            //recursively push as needed into internal nodes
            while self.nodes[current_node_index].children.len() == MAX_NODE_SIZE {
                assert!(self.nodes[current_node_index].total_counts.len() == MAX_NODE_SIZE);
                
                let current_node = &mut self.nodes[current_node_index];
                
                //changes to left child mostly occur in-place by removing the right side data
                let rtc = current_node.total_counts.drain(NODE_MIDPOINT..).collect();
                let rts = current_node.total_symbols.drain(NODE_MIDPOINT..).collect();
                let rtch = current_node.children.drain(NODE_MIDPOINT..).collect();
                
                //get the total counts for the new left node
                let left_total_count: u64 = current_node.total_counts.iter().sum();
                let mut left_total_symbols: [u64; VC_LEN] = [0; VC_LEN];
                for ts in current_node.total_symbols.iter() {
                    for i in 0..VC_LEN {
                        left_total_symbols[i] += ts[i];
                    }
                }

                //right child is create from the removed data
                let mut right_child = RLEBPlusNode {
                    is_leaf: current_node.is_leaf,
                    parent: current_node.parent,
                    total_counts: rtc,
                    total_symbols: rts,
                    children: rtch
                };

                //all node changes have been made to current & right at this point
                current_node.build_indices();
                right_child.build_indices();

                //get the total counts for the new right node
                let right_total_count: u64 = right_child.total_counts.iter().sum();
                let mut right_total_symbols: [u64; VC_LEN] = [0; VC_LEN];
                for ts in right_child.total_symbols.iter() {
                    for i in 0..VC_LEN {
                        right_total_symbols[i] += ts[i];
                    }
                }
                
                if current_node_index == 0 {
                    //special things because it's root:
                    //we increase the height and we actually create a full new left child because root needs to be a new node
                    self.height += 1;
                    let left_child = current_node.clone();

                    //parent currently does not need to be set because root's parent is 0 (e.g. also root)

                    //push the two new children
                    let left_child_index = self.nodes.len();
                    let right_child_index = left_child_index+1;
                    if !left_child.is_leaf {
                        for &cindex in left_child.children.iter() {
                            self.nodes[cindex].parent = left_child_index;
                        }
                        for &cindex in right_child.children.iter() {
                            self.nodes[cindex].parent = right_child_index;
                        }
                    }
                    self.nodes.push(left_child);
                    self.nodes.push(right_child);
                    
                    //update root with new info
                    self.nodes[0].is_leaf = false;
                    
                    self.nodes[0].total_counts[0] = left_total_count;
                    self.nodes[0].total_counts[1] = right_total_count;
                    self.nodes[0].total_counts.truncate(2);
                    
                    self.nodes[0].total_symbols[0] = left_total_symbols;
                    self.nodes[0].total_symbols[1] = right_total_symbols;
                    self.nodes[0].total_symbols.truncate(2);
                    
                    self.nodes[0].children[0] = left_child_index;
                    self.nodes[0].children[1] = right_child_index;
                    self.nodes[0].children.truncate(2);
                    
                    //root node changes are complete
                    self.nodes[0].build_indices();

                } else {
                    //this is not root, so left child is completely done at this point; we need to insert right child though
                    //push the midpoint before counts into the parent and add the new (right) child
                    let right_child_index = self.nodes.len();
                    if !right_child.is_leaf {
                        for cindex in right_child.children.iter() {
                            self.nodes[*cindex].parent = right_child_index;
                        }
                    }
                    
                    //calculate the midpoint information
                    let parent_index = right_child.parent;
                    let parent_node = &mut self.nodes[parent_index];
                    let child_index = parent_node.children.iter().position(|&cindex| cindex == current_node_index).unwrap();
                    parent_node.children.insert(child_index+1, right_child_index);
                    
                    //update the before symbols & counts
                    parent_node.total_counts[child_index] = left_total_count;
                    parent_node.total_counts.insert(child_index+1, right_total_count);
                    
                    parent_node.total_symbols[child_index] = left_total_symbols;
                    parent_node.total_symbols.insert(child_index+1, right_total_symbols);

                    //parent node changes are complete
                    parent_node.build_indices();

                    //finally, push the right node as a new node
                    self.nodes.push(right_child);

                    //now go up a level to make sure it isn't too big
                    current_node_index = parent_index;
                }
            }
        }
    }

    /// This will convert the stored data into a plain Vec<u8> for easy manipulation.
    /// This is really most useful for debugging and testing.
    /// # Example
    /// ```rust
    /// use msbwt2::rle_bplus_tree::RLEBPlusTree;
    /// let mut tree: RLEBPlusTree = Default::default();
    /// let data: Vec<u8> =             vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
    /// let expected_counts: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
    /// for (i, v) in data.iter().enumerate() {
    ///     let _count = tree.insert_and_count(i as u64, *v);
    /// }
    /// assert_eq!(tree.to_vec(), data);
    /// ```
    pub fn to_vec(&self) -> Vec<u8> {
        let mut ret: Vec<u8> = vec![];
        let mut level_vec: Vec<usize>;
        let mut next_level_vec: Vec<usize> = vec![0];

        while !self.nodes[next_level_vec[0]].is_leaf {
            level_vec = next_level_vec;
            next_level_vec = vec![];

            for lv in level_vec.iter() {
                next_level_vec.extend(&self.nodes[*lv].children);
            }
        }

        for lv in next_level_vec.iter() {
            for block_index in self.nodes[*lv].children.iter() {
                ret.extend_from_slice(&self.data_arena[*block_index].to_vec());
            }
        }
        ret

        /*
        TODO: functionally, the below statement can replace the entire above code, but I would guess it
        is less efficient, do we care? I don't think we plan to use to_vec in practice anywhere, just for
        testing currently
        benchmarks have the current to_vec about 40 usecs vs 100 usecs in the iter method
        */
        //self.into_iter().collect()
    }

    /// This will return a run iterator over the data stored in the tree.
    /// # Example
    /// ```rust
    /// use msbwt2::rle_bplus_tree::RLEBPlusTree;
    /// let mut tree: RLEBPlusTree = Default::default();
    /// let data: Vec<u8> =  vec![0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0];
    /// for (i, v) in data.iter().enumerate() {
    ///     let _count = tree.insert_and_count(i as u64, *v);
    /// }
    /// let expected_runs: Vec<(u8, u64)> = vec![(0, 1), (1, 3), (2, 6), (0, 1)];
    /// assert_eq!(tree.run_iter().collect::<Vec<(u8, u64)>>(), expected_runs);
    /// ```
    pub fn run_iter(&self) -> RLEBPlusTreeRunIterator<'_> {
        let next_child_index = self.next_child[0];
        let current_block_iter = Box::new(self.data_arena[0].run_iter());
        RLEBPlusTreeRunIterator {
            tree: self,
            next_child_index,
            current_block_iter,
            next_sym: 0,
            next_count: 0
        }
    }
}

/// An iterator over the data nodes of the B+ tree that will perform an in-order traversal of the contained characters.
/// Any runs are automatically decompressed into single symbols.
pub struct RLEBPlusTreeIterator<'a> {
    tree: &'a RLEBPlusTree,
    next_child_index: usize,
    current_block_iter: std::vec::IntoIter<u8>
}

impl<'a> Iterator for RLEBPlusTreeIterator<'a> {
    type Item = u8;
    /// Will return the next symbol contained by the compressed B+ tree data
    fn next(&mut self) -> Option<u8> {
        match self.current_block_iter.next() {
            //current block has data, so return it
            Some(x) => Some(x),
            None => {
                //current block does not have data
                if self.next_child_index == 0 {
                    //no more blocks left (we've circled around), so no more data
                    None
                } else {
                    //get the next block iterator and update the next child index before returning a value
                    self.current_block_iter = self.tree.data_arena[self.next_child_index].to_vec().into_iter();
                    self.next_child_index = self.tree.next_child[self.next_child_index];

                    //this *should* always work given our current setup, there shouldn't be any empty
                    //blocks unless the entire tree is empty
                    self.current_block_iter.next()
                }
            }
        }
    }
}

/// An iterator over the data nodes of the B+ tree that will perform an in-order traversal of the contained characters.
/// Each run is returned as a tuple (symbol (u8), count (u64)) instead of an individual symbol at a time.
pub struct RLEBPlusTreeRunIterator<'a> {
    tree: &'a RLEBPlusTree,
    next_child_index: usize,
    current_block_iter: Box<dyn Iterator<Item=(u8, u64)> + 'a>, //std::slice::Iter<'a, (u8, u64)>,
    next_sym: u8,
    next_count: u64
}

impl<'a> Iterator for RLEBPlusTreeRunIterator<'a> {
    type Item = (u8, u64);
    /// Will return the next run contained by the compressed B+ tree data
    fn next(&mut self) -> Option<(u8, u64)> {
        //None
        loop {
            let next_pair = match self.current_block_iter.next() {
                //current block has data, so return it
                Some(x) => Some(x),
                None => {
                    //current block does not have data
                    if self.next_child_index == 0 {
                        //no more blocks left (we've circled around), so no more data
                        None
                    } else {
                        //get the next block iterator and update the next child index before returning a value
                        self.current_block_iter = Box::new(self.tree.data_arena[self.next_child_index].run_iter());
                        self.next_child_index = self.tree.next_child[self.next_child_index];

                        //this *should* always work given our current setup, there shouldn't be any empty
                        //blocks unless the entire tree is empty
                        self.current_block_iter.next()
                    }
                }
            };

            match next_pair {
                Some(pair_values) => {
                    //check if this run matches the ongoing run
                    if pair_values.0 == self.next_sym {
                        //part of the same run, increment and loop back
                        self.next_count += pair_values.1 as u64;
                    } else {
                        //not part of the same run, store the new run values and return the current one
                        let ret_value = (self.next_sym, self.next_count);
                        self.next_sym = pair_values.0;
                        self.next_count = pair_values.1 as u64;
                        if ret_value.1 > 0 {
                            return Some(ret_value);
                        }
                    }
                },
                None => {
                    //we're at the end of the iterators
                    if self.next_count > 0 {
                        //one last pair to send
                        let ret_value = (self.next_sym, self.next_count);
                        self.next_count = 0;
                        return Some(ret_value);
                    } else {
                        //nothing remains
                        return None;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate rand;
    use super::*;
    use rand::Rng;

    #[test]
    fn test_init() {
        let tree: RLEBPlusTree = Default::default();
        assert_eq!(tree.to_vec(), Vec::<u8>::new());
        assert_eq!(tree.into_iter().collect::<Vec<u8>>(), Vec::<u8>::new());
        let full_pairing = tree.run_iter().collect::<Vec<(u8, u64)>>();
        assert_eq!(full_pairing, vec![]);
    }

    #[test]
    fn test_simple_inserts() {
        let mut tree: RLEBPlusTree = Default::default();
        let data: Vec<u8> =               vec![0, 1, 1, 1, 2, 0, 2, 3, 4, 1, 1, 1, 0];
        let expected_counts: Vec<u64> = vec![0, 0, 1, 2, 0, 1, 1, 0, 0, 3, 4, 5, 2];
        for (i, v) in data.iter().enumerate() {
            let pre_count = tree.count(i as u64, *v);
            assert_eq!(pre_count, expected_counts[i]);

            let count = tree.insert_and_count(i as u64, *v);
            println!("{} {:?}", i, tree.to_vec());
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
    fn test_failure_case_001() {
        let mut tree: RLEBPlusTree = Default::default();
        let data: Vec<u8> =               vec![3, 0, 3, 0, 5, 2, 5, 2, 3, 1, 2];
        let expected_counts: Vec<u64> = vec![0, 0, 1, 1, 0, 0, 1, 1, 2, 0, 2];
        for (i, v) in data.iter().enumerate() {
            let pre_count = tree.count(i as u64, *v);
            assert_eq!(pre_count, expected_counts[i]);

            let count = tree.insert_and_count(i as u64, *v);
            println!("{} {:?}", i, tree.to_vec());
            assert_eq!(tree.to_vec(), data[..i+1].to_vec());
            assert_eq!(count, expected_counts[i]);
        }
        assert_eq!(tree.to_vec(), data);
        
        /*
        //check the pairs since this is in-order
        let runs: Vec<(u8, u64)> = vec![(3, 1), (0, 1), (3, 1), (0, 1), (5, 1), (2, 1), (5, 1), (2, 1), (3, 1), (1, 1), (2, 1)];
        let full_pairing = tree.run_iter().collect::<Vec<(u8, u64)>>();
        assert_eq!(full_pairing, runs);
        */
    }

    #[test]
    fn test_failure_case_002() {
        /*
        //this one was cause by used binary search on the children array, silly me
        [2, 3, 3, 1, 3, 2, 4, 2, 3, 4, 5, 1, 2,    5]
        [2, 3, 3, 1, 3, 2, 4, 2, 3, 4, 5, 1, 2, 1, 5]
        */
        let mut tree: RLEBPlusTree = Default::default();
        let mut data: Vec<u8> = vec![];
        let inserted: Vec<u8> = vec![4, 5, 5, 0, 0, 0, 2, 3, 0, 1, 5, 3, 3, 5];
        let positions: Vec<usize> = vec![0, 0, 2, 1, 1, 1, 5, 0, 5, 6, 3, 9, 0, 10];
        for i in 0..inserted.len() {
            data.insert(positions[i], inserted[i]);
            let mut expected_count = 0;
            for j in 0..positions[i] {
                if data[j] == inserted[i] {
                    expected_count += 1;
                }
            }

            let pre_count = tree.count(positions[i] as u64, inserted[i]);
            assert_eq!(pre_count, expected_count);

            let count = tree.insert_and_count(positions[i] as u64, inserted[i]);
            println!("{} {:?}", i, tree.to_vec());
            assert_eq!(tree.to_vec(), data);
            assert_eq!(count, expected_count);
        }
    }

    #[test]
    fn test_failure_case_003() {
        /*
        caught by when I added an addition loop to fix a problem that actually wasn't the solution,
        instead causing this problem
        */
        let mut tree: RLEBPlusTree = Default::default();
        let mut data: Vec<u8> = vec![];
        let inserted: Vec<u8> = vec![2, 3, 1, 1, 1, 4, 3, 0, 2, 2, 1, 4, 3];
        let positions: Vec<usize> = vec![0, 1, 0, 3, 2, 5, 1, 7, 3, 7, 2, 2, 9];
        for i in 0..inserted.len() {
            data.insert(positions[i], inserted[i]);
            let mut expected_count = 0;
            for j in 0..positions[i] {
                if data[j] == inserted[i] {
                    expected_count += 1;
                }
            }

            let pre_count = tree.count(positions[i] as u64, inserted[i]);
            assert_eq!(pre_count, expected_count);
            
            let count = tree.insert_and_count(positions[i] as u64, inserted[i]);
            println!("{} {:?}", i, tree.to_vec());
            assert_eq!(tree.to_vec(), data);
            assert_eq!(count, expected_count);
        }
    }

    #[test]
    fn test_10krandom_insert() {
        let mut tree: RLEBPlusTree = Default::default();
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
            assert_eq!(tree.into_iter().collect::<Vec<u8>>(), data);
        }
    }

    #[test]
    fn test_run_iter() {
        //this test is needed to make sure we are testing the run iterations across data blocks
        let mut tree: RLEBPlusTree = Default::default();
        let total_symbols = MAX_BLOCK_SIZE*256*MAX_NODE_SIZE;
        for _ in 0..total_symbols {
            let _post_count = tree.insert_and_count(0, 0);
        }

        //we expect 1 really big run of 0s
        let full_pairing = tree.run_iter().collect::<Vec<(u8, u64)>>();
        assert_eq!(full_pairing, vec![(0, total_symbols as u64)]);
        assert!(tree.get_node_count() > 1);

        //now let's break the run somewhere and just make sure things are still fine
        let break_point: u64 = (MAX_BLOCK_SIZE*256) as u64;
        tree.insert_and_count(break_point, 1);
        let full_pairing = tree.run_iter().collect::<Vec<(u8, u64)>>();
        assert_eq!(full_pairing, vec![(0, break_point), (1, 1), (0, total_symbols as u64-break_point)]);
    }
}
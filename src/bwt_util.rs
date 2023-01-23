use bitvec::prelude::*;
use itertools::Itertools;
use std::clone::Clone;
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

/// implementation of the mergeIter function as described in Holt & McMillan 2014
/// Takes the current state of the interleave bitvector, the two burrows wheeler
/// transform slices, and a mapping from character to number of lexigraphically
/// smaller characters in the BWTs provided and does a round of interleave
/// resolution
pub fn pairwise_merge_iter<T: Ord + Hash + Clone + Debug>(
    current_interleave: &BitVec<u64, Msb0>,
    bwt0: &[T],
    bwt1: &[T],
    offsets: &HashMap<T, usize>,
) -> BitVec<u64, Msb0> {
    // initial conditions
    let mut next_interleave = bitvec![u64, Msb0; 0; current_interleave.len()];
    let mut current_pos0 = 0usize;
    let mut current_pos1 = 0usize;

    let mut temp_index = offsets.clone();

    // iterate through each bit
    for b in current_interleave.iter() {
        let c = if *b {
            let v = bwt0[current_pos0].clone();
            current_pos0 += 1;
            v
        } else {
            let v = bwt1[current_pos1].clone();
            current_pos1 += 1;
            v
        };
        // copy b into the next position for c
        let next_position = temp_index
            .get_mut(&c)
            .expect("Unknown character encountered in merging BWTs");
        next_interleave.set(*next_position, *b);
        // update the tempIndex to match the FM-index
        *next_position += 1;
    }

    next_interleave
}

/// Determines the offset hashmap for a provided slice of BWTs.
/// The offset hashmap is the mapping from a character to the number of
/// lexigraphically lesser characters present in the combined string provided.
fn generate_offset_hashmap<T: Ord + Hash + Clone + Debug>(bwts: &[&[T]]) -> HashMap<T, usize> {
    let mut num_occurrences: HashMap<T, usize> = HashMap::new();
    for &bwt in bwts.iter() {
        for c in bwt {
            if let Some(i) = num_occurrences.get_mut(c) {
                *i += 1;
            } else {
                num_occurrences.insert(c.clone(), 1);
            }
        }
    }
    let ordered_chars = num_occurrences.keys().sorted().cloned().collect::<Vec<T>>();
    let mut total = 0usize;
    let mut offset_map = HashMap::with_capacity(num_occurrences.len());
    for c in ordered_chars {
        offset_map.insert(c.clone(), total);
        total += num_occurrences[&c];
    }
    offset_map
}

/// Merges two Burrows-Wheeler transformed slices into a new, merged Burrows-Wheeler
/// slice.
pub fn pairwise_bwt_merge<T: Ord + Hash + Clone + Debug>(bwt0: &[T], bwt1: &[T]) -> Vec<T> {
    let total_len = bwt0.len() + bwt1.len();
    let initial_offsets = generate_offset_hashmap(&[bwt0, bwt1]);

    let mut interleave = bitvec![u64, Msb0; 0; total_len];
    // initialize ret array to 0s followed by 1s
    let mut final_interleave = bitvec![u64, Msb0; 0; total_len];
    for i in bwt0.len()..total_len {
        // TODO - figure out if there's a way to do this faster than iterating through each element
        // and setting them one by one. (BitVec should have option to do a bitwise OR)
        final_interleave.set(i, true)
    }
    while interleave != final_interleave {
        // copy the old interleaving and re-iterate
        interleave = final_interleave.clone();
        final_interleave = pairwise_merge_iter(&interleave, bwt0, bwt1, &initial_offsets);
    }

    let mut return_val = Vec::with_capacity(total_len);

    let mut current_pos0 = 0usize;
    let mut current_pos1 = 0usize;

    for val in final_interleave.into_iter() {
        if val {
            return_val.push(bwt0[current_pos0].clone());
            current_pos0 += 1;
        } else {
            return_val.push(bwt1[current_pos1].clone());
            current_pos1 += 1;
        }
    }
    return_val
}

/// This function will take a collection of strings and naively calculate the MSBWT for those strings.
/// This process is can be very slow, so this is really only useful for small datasets and testing simple use cases to verify correctness.
/// # Arguments
/// * `inputs` - the collection of strings to get converted into a MSBWT
/// # Examples
/// ```rust
/// use msbwt2::bwt_util::naive_bwt;
/// let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
/// let bwt_stream = naive_bwt(&data);
/// assert_eq!(bwt_stream, "GTN$$ACCC$G");
/// ```
pub fn naive_bwt(inputs: &[&str]) -> String {
    let mut rotations: Vec<String> = vec![];
    for s in inputs.iter() {
        let dollar_string = s.to_string() + "$";
        for l in 0..dollar_string.len() {
            rotations.push(
                //we have to loop the string twice in the event they are not all equal lengths to break
                dollar_string[l..].to_string() + &dollar_string + &dollar_string[..l],
            );
        }
    }
    rotations.sort();
    let mut ret: String = String::with_capacity(rotations.len());
    for r in rotations.iter() {
        ret.push(r.as_bytes()[r.len() - 1] as char);
    }
    ret
}

/*
// we might need this later
pub fn write_strings_to_fqgz(data: Vec<&str>) -> NamedTempFile {
    let file: NamedTempFile = Builder::new().prefix("temp_data_").suffix(".fq.gz").tempfile().unwrap();
    let mut gz = GzBuilder::new().write(file, Compression::default());
    let mut i: usize = 0;
    for s in data {
        writeln!(gz, "@seq_{}\n{}\n+\n{}", i, s, "F".repeat(s.len())).unwrap();
        i += 1;
    }

    //have to keep the file handle or everything blows up
    gz.finish().unwrap()
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        //build the BWT and make sure it's right
        let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
        let bwt_stream = naive_bwt(&data);
        assert_eq!(bwt_stream, "GTN$$ACCC$G");
    }

    #[test]
    fn test_diff_len() {
        //build the BWT and make sure it's right
        let data: Vec<&str> = vec!["A", "AA", "AAA"];
        /*
        $A$A
        $AA$AA
        $AAA$AAA
        A$A$
        A$AA$A
        A$AAA$AA
        AA$AA$
        AA$AAA$A
        AAA$AAA$
        */
        let bwt_stream = naive_bwt(&data);
        assert_eq!(bwt_stream, "AAA$AA$A$");
    }

    #[test]
    fn test_cycle_breaker() {
        //build the BWT and make sure it's right
        let data: Vec<&str> = vec!["ACA", "CA"];
        /*
        //this test case break if you don't loop the sequence twice
        $ACA
        $CA
        A$AC
        A$C
        ACA$
        CA$A
        CA$
        */
        let bwt_stream = naive_bwt(&data);
        assert_eq!(bwt_stream, "AACC$A$");
    }

    #[test]
    fn merging_paper_example_works() {
        let bwt0 = "AC$CA".as_bytes();
        let bwt1 = "AAAC$".as_bytes();
        let merged_bwt = pairwise_bwt_merge(bwt0, bwt1);
        assert_eq!(merged_bwt, "AACAAC$C$A".as_bytes())
    }
}

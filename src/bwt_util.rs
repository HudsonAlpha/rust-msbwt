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
///
/// # Arguments
/// * `current_interleave` - The bitvector indicating which bin the corresponding character in the
///                          merge currently lies.
/// * `bwt0`               - The first burrows wheeler transform
/// * `bwt1`               - The second burrows wheeler transform
/// * `offsets`            - A mapping of a key of type ``T`` to the number of elements of ``T``
///                          in the bwt that are smaller than the key.
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
        let c: &T = if *b {
            let v = &bwt0[current_pos0];
            current_pos0 += 1;
            v
        } else {
            let v = &bwt1[current_pos1];
            current_pos1 += 1;
            v
        };
        // copy b into the next position for c
        let next_position = temp_index
            .get_mut(c)
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
///
/// # Arguments
/// * `bwts` - the collection of strings that will be used in producing a MSBWT.
/// # Examples
/// ```rust
/// use msbwt2::bwt_util::generate_offset_hashmap;
/// use std::collections::HashMap;
/// let data: Vec<&[u8]> = vec!["CCGT".as_bytes(), "ACG".as_bytes()];
///
/// let offsets: HashMap<u8, usize> = generate_offset_hashmap(&data);
///
/// assert_eq!(offsets, vec![(b'A', 0), (b'C', 1), (b'G', 4), (b'T', 6)].into_iter().collect::<HashMap<u8, usize>>() );
/// ```
pub fn generate_offset_hashmap<T: Ord + Hash + Clone + Debug>(bwts: &[&[T]]) -> HashMap<T, usize> {
    let mut num_occurrences: HashMap<&T, usize> = HashMap::new();
    for &bwt in bwts.iter() {
        for c in bwt {
            num_occurrences.entry(c).and_modify(|counter| *counter += 1).or_insert(1);
        }
    }
    let ordered_chars = num_occurrences.keys().sorted().map(|&c| c.clone()).collect::<Vec<T>>();
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
/// # Arguments
/// * `bwt0` - the first of the BWT values to merge pairwise
/// * `bwt1` - the second of the BWT values to merge pairwise
/// # Examples
/// ```rust
/// use msbwt2::bwt_util::naive_bwt;
/// use msbwt2::bwt_util::pairwise_bwt_merge;
///
/// let data: Vec<&str> = vec!["CCGT", "ACG"];
///
/// let bwt_stream = naive_bwt(&data);
///
/// let bwts: Vec<String> = data.into_iter().map(|s| {naive_bwt(&[s])}).collect::<Vec<_>>();
/// let pairwise_stream = pairwise_bwt_merge(bwts[0].as_bytes(), bwts[1].as_bytes());
/// assert_eq!(bwt_stream.as_bytes(), pairwise_stream);
/// ```
pub fn pairwise_bwt_merge<T: Ord + Hash + Clone + Debug>(bwt0: &[T], bwt1: &[T]) -> Vec<T> {
    let total_len = bwt0.len() + bwt1.len();
    let initial_offsets = generate_offset_hashmap(&[bwt0, bwt1]);

    let mut interleave = bitvec![u64, Msb0; 0; total_len];
    // initialize ret array to 0s followed by 1s
    let mut final_interleave = bitvec![u64, Msb0; 0; total_len];
    for i in 0..bwt0.len() {
        // TODO - figure out if there's a way to do this faster than iterating through each element
        // and setting them one by one. (BitVec should have option to do a bitwise OR)
        final_interleave.set(i, true)
    }
    while interleave != final_interleave {
        // copy the old interleaving and re-iterate
        interleave = final_interleave;
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

    fn test_recursive_merging_result_vs_naive_result(data: Vec<&str>) {
        // run the naive burrow wheeler transform
        let bwt_stream = naive_bwt(&data);

        // run the merging burrows wheeler transforms
        let mut bwts: Vec<String> = data.iter().map(|s| naive_bwt(&[s])).collect::<Vec<_>>();
        let mut current_bwt = bwts.pop().unwrap().as_bytes().to_vec();
        while let Some(next_bwt) = bwts.pop() {
            current_bwt = pairwise_bwt_merge(&current_bwt, next_bwt.as_bytes());
        }

        // make sure they're identical
        assert_eq!(bwt_stream.as_bytes(), current_bwt);
    }

    #[test]
    fn merging_paper_example_works() {
        let data: Vec<&str> = vec!["ACCA", "CAAA"];
        test_recursive_merging_result_vs_naive_result(data)
    }

    #[test]
    fn merging_samples_of_different_size_example_works() {
        let data: Vec<&str> = vec!["ACCA", "CA"];
        test_recursive_merging_result_vs_naive_result(data)
    }

    #[test]
    fn merging_samples_of_high_similarity_works() {
        // getting bigger in order
        let data1: Vec<&str> = vec!["A", "AA", "AAA", "AAAA", "AAAAA"];
        test_recursive_merging_result_vs_naive_result(data1);

        // getting smaller in order
        let data2: Vec<&str> = vec!["AAAAA", "AAAA", "AAA", "AA", "A"];
        test_recursive_merging_result_vs_naive_result(data2);
    }
}

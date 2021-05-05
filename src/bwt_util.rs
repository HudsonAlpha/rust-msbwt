
/// This function will take a collection of strings and naively calculate the MSBWT for those strings.
/// This process is can be very slow, so this is really only useful for small datasets and testing simple use cases to verify correctness.
/// # Arguments
/// * `inputs` - the collection of strings to get converted into a MSBWT
/// # Examples
/// ```rust
/// use msbwt::bwt_util::naive_bwt;
/// let data: Vec<&str> = vec!["CCGT", "N", "ACG"];
/// let bwt_stream = naive_bwt(&data);
/// assert_eq!(bwt_stream, "GTN$$ACCC$G");
/// ```
pub fn naive_bwt(inputs: &Vec<&str>) -> String {
    let mut rotations: Vec<String> = vec![];
    for s in inputs.iter() {
        let dollar_string = s.to_string()+&"$".to_string();
        for l in 0..dollar_string.len() {
            rotations.push(
                //we have to loop the string twice in the event they are not all equal lengths to break
                dollar_string[l..].to_string()+
                &dollar_string+
                &dollar_string[..l]
            );
        }
    }
    rotations.sort();
    let mut ret: String = String::with_capacity(rotations.len());
    for r in rotations.iter() {
        ret.push(r.as_bytes()[r.len()-1] as char);
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
}
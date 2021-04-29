
extern crate log;
extern crate serde_json;

use log::info;
use std::io::prelude::*;
use std::fs;

use crate::msbwt_core::*;

pub struct RleBWT {
    bwt: Vec<u8>
}

impl Default for RleBWT {
    fn default() -> Self {
        Self {
            bwt: vec![]
        }
    }
}

impl BWT for RleBWT {
    fn load_numpy_file(&mut self, filename: &str) -> std::io::Result<()> {
        //read the numpy header: http://docs.scipy.org/doc/numpy-1.10.1/neps/npy-format.html
        //get the initial file size
        let file_metadata: fs::Metadata = fs::metadata(&filename)?;
        let full_file_size: u64 = file_metadata.len();

        //read the initial fixed header
        let mut file = fs::File::open(&filename)?;
        let mut init_header: Vec<u8> = vec![0; 10];
        let read_count: usize = file.read(&mut init_header[..])?;
        if read_count != 10 {
            panic!("Could not read initial 10 bytes of header for file {:?}", filename);
        }

        //read the dynamic header
        let header_len: usize = init_header[8] as usize + 256 * init_header[9] as usize;
        let mut skip_bytes: usize = 10+header_len;
        if skip_bytes % 16 != 0 {
            skip_bytes = ((skip_bytes / 16)+1)*16;
        }
        let mut skip_header: Vec<u8> = vec![0; skip_bytes-10];
        match file.read_exact(&mut skip_header[..]) {
            Ok(()) => {},
            Err(e) => {
                return Err(
                    std::io::Error::new(
                        e.kind(),
                        format!("Could not read bytes 10-{:?} of header for file {:?}, root-error {:?}", skip_bytes, filename, e)
                    )
                );
            }
        }
        
        //parse the header string for the expected length, requires a lot of manipulation of the string because of numpy header styling
        let header_string = String::from_utf8(skip_header).unwrap()
            .replace("\'", "\"")
            .replace("False", "false")
            .replace("(", "[")
            .replace(")", "]")
            .replace(", }", "}")
            .replace(", ]", "]")
            .replace(",]", "]");
        let header_dict: serde_json::Value = serde_json::from_str(&header_string)
            .unwrap_or_else(|_| panic!("Error while parsing header string: {:?}", header_string));
        let expected_length: u64 = header_dict["shape"][0].as_u64().unwrap();
        
        //check that the disk size matches our expectation
        let bwt_disk_size: u64 = full_file_size - skip_bytes as u64;
        if expected_length != bwt_disk_size {
            return Err(
                std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!("Header indicates shape of {:?}, but remaining file size is {:?}", expected_length, bwt_disk_size)
                )
            );
        }

        //finally read in everything else
        self.bwt = Vec::<u8>::with_capacity(bwt_disk_size as usize);//vec![0; bwt_disk_size as usize];
        let read_count: usize = file.read_to_end(&mut self.bwt)?;
        if read_count as u64 != bwt_disk_size {
            return Err(
                std::io::Error::new(
                    std::io::ErrorKind::UnexpectedEof,
                    format!("Only read {:?} of {:?} bytes of BWT body for file {:?}", read_count, bwt_disk_size, filename)
                )
            );
        }
        info!("Loading BWT with {:?} compressed values", bwt_disk_size);

        //we loaded the file into memory, standard init now
        self.standard_init();

        Ok(())
    }

    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange {
        BWTRange {
            l: 0,
            h: 0
        }
    }
}

impl RleBWT {
    fn standard_init(&mut self) {
        //we will call this function when the bwt is fully loaded into memory
        //TODO: the actual implementation
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_constrain_range() {
        let rle_bwt: RleBWT = Default::default();
        let initial_range = BWTRange {
            l: 0,
            h: 0
        };
        let new_range = unsafe {
            rle_bwt.constrain_range(0, &initial_range)
        };

        assert_eq!(new_range, BWTRange{l: 0, h:0});
    }
}
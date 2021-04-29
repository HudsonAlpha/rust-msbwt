
use crate::msbwt_core::*;

pub struct RleMsbwt {

}

impl Default for RleMsbwt {
    fn default() -> Self {
        Self {
            
        }
    }
}

impl Msbwt for RleMsbwt {
    unsafe fn constrain_range(&self, sym: u8, input_range: &BWTRange) -> BWTRange {
        BWTRange {
            l: 0,
            h: 0
        }
    }
}

impl Msbwt {
    
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_constrain_range() {
        let rle_bwt: RleMsbwt = Default::default();
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
// this file is copied from thrombe/fracGen look there for history

pub fn xyrange(pow: f64, x: f64, y: f64) -> (f64, f64, f64, f64) {
    let nudge = 2.0_f64.powf(pow-1.0);
    (x-nudge, x+nudge, y+nudge, y-nudge)
}

#[derive(Clone)]
pub struct MapRange {
    scale: f64,
    toff: f64,
    foff: f64,
}
impl MapRange {
    /// gives a struct. use struct.map(num) to map num
    pub fn new(fs: f64, fe: f64, ts: f64, te: f64) -> MapRange {
        MapRange {
            scale: ((te-ts)/(fe-fs)),
            toff: (te+ts)/2.0,
            foff: (fe+fs)/2.0,
            }
    }
    
    #[inline(always)]
    pub fn map(&self, num: f64) -> f64 {
        (num-self.foff)*self.scale + self.toff
    }
    
}

#[derive(Clone)]
pub struct ChopRange { // generic version of thihs would be nice
    start: f64,
    end: f64,
}
impl ChopRange {
    pub fn new (start: f64, end: f64) -> ChopRange {
        ChopRange {
            start,
            end,
        }
    }

    #[inline(always)]
    pub fn map(&self, num: f64) -> f64 {
        if num < self.start {self.start} else if num > self.end {self.end} else {num}
    }
}

#[inline(always)]
pub fn min(a: f64, b: f64) -> f64 {
    if a < b {a} else {b}
}

#[inline(always)]
pub fn max(a: f64, b: f64) -> f64 {
    if a > b {a} else {b}
}

#[inline(always)]
pub fn abs(x: f64) -> f64 {
    if x > 0.0 {x} else {0.0}
}
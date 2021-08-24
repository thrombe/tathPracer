// this file is copied from thrombe/fracGen look there for history

use std::io::Write; // for progress indicator

pub struct ProgressIndicator {
    oneperc: usize,
    current: usize,
    ended: bool,
}

impl ProgressIndicator {
    /// total is the total amount of work done by the worker
    pub fn new(total: usize) -> ProgressIndicator {
        let mut pi = ProgressIndicator {
            oneperc: total/100,
            current: 0,
            ended: false,
        };
        if total == 0 {pi.ended = true}
        pi
    }

    /// prog is the current work done by the worker, should end at total
    #[inline(always)]
    pub fn indicate(&mut self, mut prog: usize) {
        if self.ended {return}
        if prog % self.oneperc == 0 {
            prog = prog/self.oneperc;
            if self.current != prog {
                self.current = prog;
                if prog != 0 {print!("\x08\x08\x08")}
                match prog {
                    0..=9 => print!("0{}%", prog),
                    10..=99 => print!("{}%", prog),
                    _ => {
                        println!("");
                        self.ended = true;
                        return            
                    },
                }
                std::io::stdout().flush().unwrap(); // needs std::io::Write
            }
        }
    }
}
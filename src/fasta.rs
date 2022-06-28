// This module is for reading and writing FASTA format files

use std::io::BufRead;
use std::str::FromStr;
use crate::DnaSequence;

pub struct SimpleFastaParser {
}

impl SimpleFastaParser {
    pub fn parse<R: BufRead>(handle: &mut R) -> Result<Vec<(String, String)>, anyhow::Error> {
        let mut result: Vec<(String, String)> = vec![];

        let mut lines: Vec<String> = vec![];
        let mut title: Option<String> = None;

        for line in handle.lines() {
            match line {
                Ok(l) => {
                    if let Some(ch) = l.chars().nth(0) {
                        if ch == '>' {
                            if let Some(existing_title) = title {
                                result.push((existing_title, lines.join("")));
                                lines.clear();
                            }

                            let mut c = l.chars();
                            c.next();
                            title = Some(c.as_str().trim().to_string());
                        } else {
                            match title {
                                None => {
                                    // There is no title, skip this content!
                                }
                                Some(_) => {
                                    lines.push(l.replace(" ", "").replace("\r", ""));
                                }
                            }

                        }
                    }
                }
                _ => {}
            }
        }

        if let Some(t) = title {
            result.push((t, lines.join("")));
        }

        return Ok(result);
    }
}

pub struct DNAFastaParser {

}

impl DNAFastaParser {
    pub fn parse<R: BufRead>(handle: &mut R) -> Result<Vec<(String, DnaSequence)>, anyhow::Error> {
        SimpleFastaParser::parse(handle)?.into_iter().map(|(title, sequence)| -> Result<(String, DnaSequence), _> {
            let dna = DnaSequence::from_str(&*sequence)?;
            Ok((title, dna))
        }).collect()
    }
}


#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use std::str::FromStr;
    use crate::DnaSequence;
    use crate::fasta::{DNAFastaParser, SimpleFastaParser};

    #[test]
    fn test_empty_fasta() {
        let r = SimpleFastaParser::parse(&mut BufReader::new("".as_bytes())).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_empty_fasta_with_newlines() {
        let r = SimpleFastaParser::parse(&mut BufReader::new("\n\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_empty_fasta_with_many_newlines() {
        let r = SimpleFastaParser::parse(&mut BufReader::new("  \n\n \n  \r\r \r\n  \r \n".as_bytes())).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_fasta_with_no_content() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "".to_string()));
    }

    #[test]
    fn test_fasta_with_empty_content() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus\n\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "".to_string()));
    }

    #[test]
    fn test_fasta_with_stuff_no_header() {
        let r = SimpleFastaParser::parse(&mut BufReader::new("// this is a file comment\n@author is foo\n\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_fasta_with_stuff_before_header() {
        let r = SimpleFastaParser::parse(&mut BufReader::new("// this is a file comment\n@author is foo\n\n>Virus\n\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "".to_string()));
    }

    #[test]
    fn test_fasta_with_single_line_content() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus\nCAAAGT\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "CAAAGT".to_string()));
    }

    #[test]
    fn test_fasta_with_multi_line_content() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus\nAAAA\nCCCC\nGGGG\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "AAAACCCCGGGG".to_string()));
    }

    #[test]
    fn test_fasta_mutiple_contents() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus1\nAAAA\n>Virus2\nCCCC\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], ("Virus1".to_string(), "AAAA".to_string()));
        assert_eq!(r[1], ("Virus2".to_string(), "CCCC".to_string()));
    }

    #[test]
    fn test_fasta_mutiple_contents_multiline() {
        let r = SimpleFastaParser::parse(&mut BufReader::new(">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], ("Virus1".to_string(), "AAAAAAAA".to_string()));
        assert_eq!(r[1], ("Virus2".to_string(), "CCCCCCCC".to_string()));
    }

    #[test]
    fn test_dna_fasta() {
        let r: Vec<(String, DnaSequence)> = DNAFastaParser::parse(&mut BufReader::new(">Virus1\nAAAA".as_bytes())).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus1".to_string(), DnaSequence::from_str("AAAA").unwrap()));
    }

    #[test]
    fn test_dna_fasta_multiple() {
        let r: Vec<(String, DnaSequence)> = DNAFastaParser::parse(&mut BufReader::new(">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n".as_bytes())).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], ("Virus1".to_string(), DnaSequence::from_str("AAAAAAAA").unwrap()));
        assert_eq!(r[1], ("Virus2".to_string(), DnaSequence::from_str("CCCCCCCC").unwrap()));
    }

    #[test]
    fn test_dna_invalid_dna() {
        let r = DNAFastaParser::parse(&mut BufReader::new(">Virus1\nAAAAelephant".as_bytes()));
        assert!(matches!(r, Err(_)));
    }

    #[test]
    fn test_dna_invalid_dna_multiple() {
        let r = DNAFastaParser::parse(&mut BufReader::new(">Virus1\nAAAA\n>Virus2\nAAAAelephant".as_bytes()));
        assert!(matches!(r, Err(_)));
    }
}
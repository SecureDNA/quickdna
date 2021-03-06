//! This module is for reading and writing FASTA format files

use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

use thiserror::Error;

use crate::{DnaSequence, TranslationError};

/// Note that the parser is buffered automatically, so you should not wrap handle in a buffered reader like io::BufReader.
pub struct SimpleFastaParser;

impl SimpleFastaParser {
    pub fn parse<R: Read>(handle: R) -> Result<Vec<(String, String)>, std::io::Error> {
        let mut result: Vec<(String, String)> = vec![];

        let mut lines: Vec<String> = vec![];
        let mut title: Option<String> = None;

        let handle = BufReader::new(handle);

        for line in handle.lines() {
            let line = line?;

            if let Some(ch) = line.chars().next() {
                if ch == '>' {
                    if let Some(existing_title) = title {
                        result.push((existing_title, lines.join("")));
                        lines.clear();
                    }

                    let mut c = line.chars();
                    c.next();
                    title = Some(c.as_str().trim().to_string());
                } else {
                    match title {
                        None => {
                            // There is no title, skip this content!
                        }
                        Some(_) => {
                            lines.push(line.replace(' ', "").replace('\r', ""));
                        }
                    }
                }
            }
        }

        if let Some(t) = title {
            result.push((t, lines.join("")));
        }

        Ok(result)
    }
}

/// Note that the parser is buffered automatically, so you should not wrap handle in a buffered reader like io::BufReader.
pub struct DnaFastaParser;

#[derive(Error, Debug)]
pub enum DnaFastaParseError {
    #[error("could not translate FASTA contents")]
    TranslationError(#[from] TranslationError),
    #[error("could not read FASTA")]
    IOError(#[from] std::io::Error),
}

impl DnaFastaParser {
    pub fn parse<R: Read>(handle: R) -> Result<Vec<(String, DnaSequence)>, DnaFastaParseError> {
        SimpleFastaParser::parse(handle)?
            .into_iter()
            .map(|(title, sequence)| -> Result<(String, DnaSequence), _> {
                let dna = DnaSequence::from_str(&*sequence)?;
                Ok((title, dna))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_fasta() {
        let r = SimpleFastaParser::parse("".as_bytes()).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_empty_fasta_with_newlines() {
        let r = SimpleFastaParser::parse("\n\n".as_bytes()).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_empty_fasta_with_many_newlines() {
        let r = SimpleFastaParser::parse("  \n\n \n  \r\r \r\n  \r \n".as_bytes()).unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_fasta_with_no_content() {
        let r = SimpleFastaParser::parse(">Virus\n".as_bytes()).unwrap();
        assert_eq!(r.len(), 1);
        assert_eq!(r[0], ("Virus".to_string(), "".to_string()));
    }

    #[test]
    fn test_fasta_with_empty_content() {
        let r = SimpleFastaParser::parse(">Virus\n\n".as_bytes()).unwrap();
        assert_eq!(r, vec![("Virus".to_string(), "".to_string())]);
    }

    #[test]
    fn test_fasta_with_stuff_no_header() {
        let r =
            SimpleFastaParser::parse("// this is a file comment\n@author is foo\n\n".as_bytes())
                .unwrap();
        assert_eq!(r.len(), 0);
    }

    #[test]
    fn test_fasta_with_stuff_before_header() {
        let r = SimpleFastaParser::parse(
            "// this is a file comment\n@author is foo\n\n>Virus\n\n".as_bytes(),
        )
        .unwrap();
        assert_eq!(r, vec![("Virus".to_string(), "".to_string())]);
    }

    #[test]
    fn test_fasta_with_single_line_content() {
        let r = SimpleFastaParser::parse(">Virus\nCAAAGT\n".as_bytes()).unwrap();
        assert_eq!(r, vec![("Virus".to_string(), "CAAAGT".to_string())]);
    }

    #[test]
    fn test_fasta_with_multi_line_content() {
        let r = SimpleFastaParser::parse(">Virus\nAAAA\nCCCC\nGGGG\n".as_bytes()).unwrap();
        assert_eq!(r, vec![("Virus".to_string(), "AAAACCCCGGGG".to_string())]);
    }

    #[test]
    fn test_fasta_mutiple_contents() {
        let r = SimpleFastaParser::parse(">Virus1\nAAAA\n>Virus2\nCCCC\n".as_bytes()).unwrap();
        assert_eq!(
            r,
            vec![
                ("Virus1".to_string(), "AAAA".to_string()),
                ("Virus2".to_string(), "CCCC".to_string())
            ]
        );
    }

    #[test]
    fn test_fasta_mutiple_contents_multiline() {
        let r = SimpleFastaParser::parse(">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n".as_bytes())
            .unwrap();
        assert_eq!(
            r,
            vec![
                ("Virus1".to_string(), "AAAAAAAA".to_string()),
                ("Virus2".to_string(), "CCCCCCCC".to_string())
            ]
        );
    }

    #[test]
    fn test_dna_fasta() {
        let r: Vec<(String, DnaSequence)> =
            DnaFastaParser::parse(">Virus1\nAAAA".as_bytes()).unwrap();

        assert_eq!(
            r,
            vec![("Virus1".to_string(), DnaSequence::from_str("AAAA").unwrap())]
        );
    }

    #[test]
    fn test_dna_fasta_multiple() {
        let r: Vec<(String, DnaSequence)> =
            DnaFastaParser::parse(">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n".as_bytes()).unwrap();
        assert_eq!(
            r,
            vec![
                (
                    "Virus1".to_string(),
                    DnaSequence::from_str("AAAAAAAA").unwrap()
                ),
                (
                    "Virus2".to_string(),
                    DnaSequence::from_str("CCCCCCCC").unwrap()
                )
            ]
        );
    }

    #[test]
    fn test_dna_invalid_dna() {
        let r = DnaFastaParser::parse(">Virus1\nAAAAelephant".as_bytes());
        assert!(matches!(r, Err(_)));
    }

    #[test]
    fn test_dna_invalid_dna_multiple() {
        let r = DnaFastaParser::parse(">Virus1\nAAAA\n>Virus2\nAAAAelephant".as_bytes());
        assert!(matches!(r, Err(_)));
    }
}

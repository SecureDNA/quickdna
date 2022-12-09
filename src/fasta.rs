//! This module is for reading and writing FASTA format files

use std::fmt::Display;
use std::io::{self, BufRead};
use std::str::FromStr;

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord<T> {
    /// The header of this record, without the leading character (usually '>')
    /// Depending on the content and parser settings, this header may be empty, and may contain newlines.
    pub header: String,
    /// The contents of this record
    pub contents: T,
    /// The starting and ending line numbers of this record, start inclusive, end exclusive, 1-indexed.
    /// The record header is included in this range.
    pub line_range: (usize, usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaFile<T> {
    /// The records parsed from the file.
    pub records: Vec<FastaRecord<T>>,
}

impl FastaRecord<String> {
    /// Try to parse a string-containing record into a type that impls [`core::str::FromStr`].
    pub fn parse<T: FromStr>(self) -> Result<FastaRecord<T>, T::Err> {
        let Self {
            header,
            contents,
            line_range,
        } = self;
        let contents = contents.parse()?;
        Ok(FastaRecord {
            header,
            contents,
            line_range,
        })
    }
}

impl<T: ToString> Display for FastaRecord<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if !self.header.is_empty() {
            writeln!(f, ">{}", self.header.replace('\n', "\n>"))?;
        }
        writeln!(f, "{}", self.contents.to_string())
    }
}

impl<T: ToString> Display for FastaFile<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for record in &self.records {
            write!(f, "{}", record)?;
        }
        Ok(())
    }
}

/// Settings for a fasta parser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastaParseSettings {
    /// If this flag is true, then successive headers in a FASTA file will be concatenated, instead of
    /// generating empty records.
    ///
    /// ```rust
    /// use quickdna::{FastaParser, FastaParseSettings, FastaRecord};
    ///
    /// let nonconcat = FastaParser::<String>::new(
    ///     FastaParseSettings::new().concatenate_headers(false)
    /// ).parse_str(">a\n>b\n...").unwrap();
    /// assert_eq!(
    ///     nonconcat.records,
    ///     vec![
    ///         FastaRecord {
    ///             header: "a".to_string(),
    ///             contents: "".to_string(),
    ///             line_range: (1, 2),
    ///         },
    ///         FastaRecord {
    ///             header: "b".to_string(),
    ///             contents: "...".to_string(),
    ///             line_range: (2, 4),
    ///         }
    ///     ],
    /// );
    ///
    /// let concat = FastaParser::<String>::new(
    ///     FastaParseSettings::new().concatenate_headers(true)
    /// ).parse_str(">a\n>b\n...").unwrap();
    /// assert_eq!(
    ///     concat.records,
    ///     vec![
    ///         FastaRecord {
    ///             header: "a\nb".to_string(),
    ///             contents: "...".to_string(),
    ///             line_range: (1, 4),
    ///         }
    ///     ],
    /// );
    /// ```
    concatenate_headers: bool,

    /// If this flag is true, content before the first header will be treated as a file comment and ignored. Otherwise,
    /// it will be parsed as a sequence with an empty header.
    ///
    /// ```rust
    /// use quickdna::{FastaParser, FastaParseSettings, FastaRecord};
    ///
    /// let comment = FastaParser::<String>::new(
    ///     FastaParseSettings::new().allow_preceding_comment(true)
    /// ).parse_str("comment\n>a\nsequence").unwrap();
    /// assert_eq!(
    ///     comment.records,
    ///     vec![
    ///         FastaRecord {
    ///             header: "a".to_string(),
    ///             contents: "sequence".to_string(),
    ///             line_range: (2, 4),
    ///         }
    ///     ]
    /// );
    ///
    /// let comment = FastaParser::<String>::new(
    ///     FastaParseSettings::new().allow_preceding_comment(false)
    /// ).parse_str("comment\n>a\nsequence").unwrap();
    /// assert_eq!(
    ///     comment.records,
    ///     vec![
    ///         FastaRecord {
    ///             header: "".to_string(),
    ///             contents: "comment".to_string(),
    ///             line_range: (1, 2),
    ///         },
    ///         FastaRecord {
    ///             header: "a".to_string(),
    ///             contents: "sequence".to_string(),
    ///             line_range: (2, 4),
    ///         }
    ///     ]
    /// );
    /// ```
    allow_preceding_comment: bool,
}

// "Builder-lite" pattern: https://matklad.github.io/2022/05/29/builder-lite.html
impl FastaParseSettings {
    /// Initializes settings "strictly" (See [`Self::strict()`]).
    pub fn new() -> Self {
        Self::strict()
    }

    /// Initializes settings "strictly":
    /// * [`Self::concatenate_headers`]: false
    /// * [`Self::allow_preceding_comment`]: false
    pub fn strict() -> Self {
        Self {
            concatenate_headers: false,
            allow_preceding_comment: false,
        }
    }

    /// Initializes settings "laxly":
    /// * [`Self::concatenate_headers`]: true
    /// * [`Self::allow_preceding_comment`]: true
    pub fn lax() -> Self {
        Self {
            concatenate_headers: true,
            allow_preceding_comment: true,
        }
    }

    /// Changes the setting for [`Self::concatenate_headers`]
    pub fn concatenate_headers(mut self, concatenate_headers: bool) -> Self {
        self.concatenate_headers = concatenate_headers;
        self
    }

    /// Changes the setting for [`Self::allow_preceding_comment`]
    pub fn allow_preceding_comment(mut self, allow_preceding_comment: bool) -> Self {
        self.allow_preceding_comment = allow_preceding_comment;
        self
    }
}

impl Default for FastaParseSettings {
    fn default() -> Self {
        Self::new()
    }
}

enum ParserState {
    StartOfFile {
        /// only non-empty if settings.allow_preceding_comment is false
        contents: String,
    },
    InHeader {
        start_line_number: usize,
        header: String,
    },
    InRecord {
        start_line_number: usize,
        header: String,
        contents: String,
    },
}

impl ParserState {
    /// Pump the state machine and maybe emit a record, if this line completes a record, along with
    /// a new state
    fn advance_line(
        self,
        settings: &FastaParseSettings,
        line: &str,
        line_number: usize,
    ) -> (Self, Option<FastaRecord<String>>) {
        let new_header = try_parse_header(line);
        match (self, new_header) {
            // start of file, and we have a header line => start a new record,
            // maybe emiting the preceding content as a headerless record depending
            // on parse settings
            (ParserState::StartOfFile { contents }, Some(new_header)) => {
                // don't emit if the settings don't want it, or if the record would just be whitespace
                let record = if settings.allow_preceding_comment || contents.trim().is_empty() {
                    None
                } else {
                    Some(FastaRecord {
                        header: "".to_string(),
                        contents,
                        line_range: (1, line_number),
                    })
                };

                (
                    Self::InHeader {
                        start_line_number: line_number,
                        header: new_header.to_string(),
                    },
                    record,
                )
            }
            // start of file, no header line yet => maybe store this line in content,
            // if we need to emit it later
            (ParserState::StartOfFile { mut contents }, None) => {
                if !settings.allow_preceding_comment {
                    // only bother to keep contents updated if we need to emit it
                    contents.push_str(line);
                }
                (Self::StartOfFile { contents }, None)
            }

            // In header, and we have a new header => either concatenate the headers
            // or emit an empty record, depending on settings
            (
                ParserState::InHeader {
                    start_line_number,
                    mut header,
                },
                Some(new_header),
            ) => {
                if settings.concatenate_headers {
                    header.push('\n');
                    header.push_str(new_header);
                    (
                        Self::InHeader {
                            start_line_number,
                            header,
                        },
                        None,
                    )
                } else {
                    (
                        Self::InHeader {
                            start_line_number: line_number,
                            header: new_header.to_string(),
                        },
                        Some(FastaRecord {
                            header,
                            contents: "".to_string(),
                            line_range: (start_line_number, line_number),
                        }),
                    )
                }
            }
            // in header and we don't have a new header => start of record content
            (
                ParserState::InHeader {
                    start_line_number,
                    header,
                },
                None,
            ) => (
                Self::InRecord {
                    start_line_number,
                    header,
                    contents: line.to_string(),
                },
                None,
            ),

            // in record and we have a new header => start of a new header
            (
                ParserState::InRecord {
                    start_line_number,
                    header,
                    contents,
                },
                Some(new_header),
            ) => (
                Self::InHeader {
                    start_line_number: line_number,
                    header: new_header.to_string(),
                },
                Some(FastaRecord {
                    header,
                    contents,
                    line_range: (start_line_number, line_number),
                }),
            ),
            // in record and we don't have a new header => continue record
            (
                ParserState::InRecord {
                    start_line_number,
                    header,
                    mut contents,
                },
                None,
            ) => {
                if line.is_empty() {
                    // don't push an empty line to a record at the end of a file with a trailing newline
                    // (if this isn't the EOF, the line numbers will continue and pushing an empty line would
                    // have been a no-op anyways)
                    (
                        Self::InRecord {
                            start_line_number,
                            header,
                            contents,
                        },
                        None,
                    )
                } else {
                    contents.push_str(line);
                    (
                        Self::InRecord {
                            start_line_number,
                            header,
                            contents,
                        },
                        None,
                    )
                }
            }
        }
    }

    /// At the end of the file, we want to pump the state machine one more time to maybe
    /// spit out a trailing record. We don't return a new state because after EOF we don't care.
    ///
    /// `eof_line_number` here should be 1+ the last line number in the file, so that if something
    /// e.g., spans a 5-line file, the end of its range can be 6 because the end is exclusive.
    fn advance_eof(
        self,
        settings: &FastaParseSettings,
        eof_line_number: usize,
    ) -> Option<FastaRecord<String>> {
        match self {
            // start of file => maybe emit contents as a headerless record, depending on settings
            ParserState::StartOfFile { contents } => {
                // again, don't emit if settings don't want it or the record would just be whitespace
                if settings.allow_preceding_comment || contents.trim().is_empty() {
                    None
                } else {
                    Some(FastaRecord {
                        header: "".to_string(),
                        contents,
                        line_range: (1, eof_line_number),
                    })
                }
            }

            // in header => emit the header as empty record
            ParserState::InHeader {
                start_line_number,
                header,
            } => Some(FastaRecord {
                header,
                contents: "".to_string(),
                line_range: (start_line_number, eof_line_number),
            }),

            // in record => emit as final record
            ParserState::InRecord {
                start_line_number,
                header,
                contents,
            } => Some(FastaRecord {
                header,
                contents,
                line_range: (start_line_number, eof_line_number),
            }),
        }
    }
}

#[derive(Debug, Error)]
pub enum FastaParseError<ParseError> {
    #[error("error reading from reader: {0}")]
    IOError(#[from] io::Error),
    #[error("error parsing record: {0}")]
    ParseError(#[source] ParseError), // can't use #[from] due to generic impl clash
}

pub struct FastaParser<T: FromStr> {
    settings: FastaParseSettings,
    _marker: std::marker::PhantomData<T>,
}

impl<T: FromStr> FastaParser<T> {
    /// Construct a new FastaParser with the given [`FastaParseSettings`]
    pub fn new(settings: FastaParseSettings) -> Self {
        Self {
            settings,
            _marker: Default::default(),
        }
    }

    /// Construct a new FastaParser with strict settings (see [`FastaParseSettings::strict()`])
    pub fn strict() -> Self {
        Self::new(FastaParseSettings::strict())
    }

    /// Construct a new FastaParser with lax settings (see [`FastaParseSettings::lax()`])
    pub fn lax() -> Self {
        Self::new(FastaParseSettings::lax())
    }

    pub fn parse<R: BufRead>(&self, handle: R) -> Result<FastaFile<T>, FastaParseError<T::Err>> {
        let mut records: Vec<FastaRecord<T>> = vec![];
        let mut state = ParserState::StartOfFile {
            contents: "".to_string(),
        };

        let mut line_number = 0;
        for (idx, line) in handle.lines().enumerate() {
            let line = line?;
            line_number = idx + 1;

            let (new_state, record) = state.advance_line(&self.settings, &line, line_number);
            state = new_state;
            if let Some(record) = record {
                records.push(record.parse().map_err(FastaParseError::ParseError)?);
            }
        }

        if let Some(record) = state.advance_eof(&self.settings, line_number + 1) {
            records.push(record.parse().map_err(FastaParseError::ParseError)?);
        }

        Ok(FastaFile { records })
    }

    pub fn parse_str(&self, s: &str) -> Result<FastaFile<T>, FastaParseError<T::Err>> {
        self.parse(s.as_bytes())
    }
}

/// Try to parse a FASTA header (prefixed with > or ;), returning the line without the prefix char.
fn try_parse_header(line: &str) -> Option<&str> {
    let head = line.chars().next();
    // ; is semi-obsolete alternative header char
    if head == Some('>') || head == Some(';') {
        Some(&line[1..]) // we know it's an ASCII char so this slice is panic-safe
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::{DnaSequence, Nucleotide, NucleotideAmbiguous, ProteinSequence, TranslationError};
    use std::time::Duration;

    macro_rules! assert_parse {
        ($testcase:expr, $parser:expr, $expected:expr) => {
            assert_eq!(
                $parser.parse_str($testcase).unwrap().records,
                $expected,
                "settings = {:?}",
                $parser.settings
            );
        };
    }

    macro_rules! assert_parse_err {
        ($testcase:expr, $parser:expr, $match:pat) => {
            // we use matches! since io::Error doesn't have a PartialEq impl
            assert!(
                matches!($parser.parse_str($testcase).unwrap_err(), $match),
                "settings = {:?}",
                $parser.settings
            )
        };
    }

    /// Assert a parse matches the expected value with any combination of parser settings
    fn assert_parse_with_all_settings(s: &str, expected: Vec<FastaRecord<String>>) {
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: true,
                allow_preceding_comment: true
            }),
            expected
        );
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: true,
                allow_preceding_comment: false
            }),
            expected
        );
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: false,
                allow_preceding_comment: true
            }),
            expected
        );
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: false,
                allow_preceding_comment: false
            }),
            expected
        );
    }

    /// Assert a parse matches the expected value with the given value of concatenate_headers,
    /// and all values for other settings
    fn assert_parse_with_concatenate_headers(
        s: &str,
        concatenate_headers: bool,
        expected: Vec<FastaRecord<String>>,
    ) {
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers,
                allow_preceding_comment: true
            }),
            expected
        );
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers,
                allow_preceding_comment: false
            }),
            expected
        );
    }

    /// Assert a parse matches the expected value with the given value of allow_preceding_comment,
    /// and all values for other settings
    fn assert_parse_with_allow_preceding_comment(
        s: &str,
        allow_preceding_comment: bool,
        expected: Vec<FastaRecord<String>>,
    ) {
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: true,
                allow_preceding_comment
            }),
            expected
        );
        assert_parse!(
            s,
            FastaParser::<String>::new(FastaParseSettings {
                concatenate_headers: false,
                allow_preceding_comment
            }),
            expected
        );
    }

    /// Helper to panic if a closure doesn't complete within a specified Duration.
    /// Author: @shepmaster, https://github.com/rust-lang/rfcs/issues/2798#issuecomment-552949300
    fn panic_after<T: Send + 'static, F: FnOnce() -> T + Send + 'static>(d: Duration, f: F) -> T {
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        let handle = std::thread::spawn(move || {
            let val = f();
            done_tx.send(()).expect("Unable to send completion signal");
            val
        });

        match done_rx.recv_timeout(d) {
            Ok(_) => handle.join().expect("Thread panicked"),
            Err(_) => panic!("Thread took too long"),
        }
    }

    #[test]
    fn test_empty_fasta() {
        assert_parse_with_all_settings("", vec![]);
    }

    #[test]
    fn test_empty_fasta_with_newlines() {
        assert_parse_with_all_settings("\n\n", vec![]);
    }

    #[test]
    fn test_empty_fasta_with_many_newlines() {
        assert_parse_with_all_settings("  \n\n \n  \r\r \r\n  \r \n", vec![]);
    }

    #[test]
    fn test_fasta_with_no_content() {
        assert_parse_with_all_settings(
            ">Virus\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (1, 2),
            }],
        );

        assert_parse_with_all_settings(
            ";Virus\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (1, 2),
            }],
        );
    }

    #[test]
    fn test_fasta_with_empty_content() {
        assert_parse_with_all_settings(
            ">Virus\n\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (1, 3),
            }],
        );
        assert_parse_with_all_settings(
            ";Virus\n\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (1, 3),
            }],
        );
    }

    #[test]
    fn test_comment_fasta_with_allow_preceding_comment() {
        assert_parse_with_allow_preceding_comment(
            "this is a file comment\n@author is foo",
            true,
            vec![],
        );
    }

    #[test]
    fn test_comment_fasta_with_no_allow_preceding_comment() {
        assert_parse_with_allow_preceding_comment(
            "this is a file comment\n@author is foo",
            false,
            vec![FastaRecord {
                header: "".to_string(),
                contents: "this is a file comment@author is foo".to_string(),
                line_range: (1, 3),
            }],
        );
    }

    #[test]
    fn test_comment_and_record_fasta_with_allow_preceding_comment() {
        assert_parse_with_allow_preceding_comment(
            "this is a file comment\n@author is foo\n\n>Virus\n\n",
            true,
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (4, 6),
            }],
        );
    }

    #[test]
    fn test_comment_and_record_fasta_with_no_allow_preceding_comment() {
        assert_parse_with_allow_preceding_comment(
            "this is a file comment\n@author is foo\n\n>Virus\n\n",
            false,
            vec![
                FastaRecord {
                    header: "".to_string(),
                    contents: "this is a file comment@author is foo".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus".to_string(),
                    contents: "".to_string(),
                    line_range: (4, 6),
                },
            ],
        );
    }

    #[test]
    fn test_no_allow_preceding_comment_doesnt_emit_whitespace_record() {
        assert_parse_with_all_settings(
            "   \t\n\r\t   \n\n>Virus\n\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "".to_string(),
                line_range: (4, 6),
            }],
        );
    }

    #[test]
    fn test_fasta_with_single_line_content() {
        assert_parse_with_all_settings(
            ">Virus\nCAAAGT\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "CAAAGT".to_string(),
                line_range: (1, 3),
            }],
        );
        assert_parse_with_all_settings(
            ";Virus\nCAAAGT\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "CAAAGT".to_string(),
                line_range: (1, 3),
            }],
        );
    }

    #[test]
    fn test_fasta_eof() {
        assert_parse_with_all_settings(
            ">Virus\nCAAAGT",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "CAAAGT".to_string(),
                line_range: (1, 3),
            }],
        );
        assert_parse_with_all_settings(
            ";Virus\nCAAAGT",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "CAAAGT".to_string(),
                line_range: (1, 3),
            }],
        );
    }

    #[test]
    fn test_fasta_with_multi_line_content() {
        assert_parse_with_all_settings(
            ">Virus\nAAAA\nCCCC\nGGGG\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "AAAACCCCGGGG".to_string(),
                line_range: (1, 5),
            }],
        );
        assert_parse_with_all_settings(
            ";Virus\nAAAA\nCCCC\nGGGG\n",
            vec![FastaRecord {
                header: "Virus".to_string(),
                contents: "AAAACCCCGGGG".to_string(),
                line_range: (1, 5),
            }],
        );
    }

    #[test]
    fn test_fasta_mutiple_contents() {
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\n>Virus2\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAA".to_string(),
                    line_range: (1, 3),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCC".to_string(),
                    line_range: (3, 5),
                },
            ],
        );
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\n;Virus2\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAA".to_string(),
                    line_range: (1, 3),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCC".to_string(),
                    line_range: (3, 5),
                },
            ],
        );
        assert_parse_with_all_settings(
            ";Virus1\nAAAA\n>Virus2\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAA".to_string(),
                    line_range: (1, 3),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCC".to_string(),
                    line_range: (3, 5),
                },
            ],
        );
        assert_parse_with_all_settings(
            ";Virus1\nAAAA\n;Virus2\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAA".to_string(),
                    line_range: (1, 3),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCC".to_string(),
                    line_range: (3, 5),
                },
            ],
        );
    }

    #[test]
    fn test_fasta_mutiple_contents_multiline() {
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
            ],
        );
        assert_parse_with_all_settings(
            ";Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
            ],
        );
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\nAAAA\n;Virus2\nCCCC\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
            ],
        );
        assert_parse_with_all_settings(
            ";Virus1\nAAAA\nAAAA\n;Virus2\nCCCC\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
            ],
        );
    }

    #[test]
    fn test_fasta_triple() {
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n>Virus3\nCCCC\nCCCC\n",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
                FastaRecord {
                    header: "Virus3".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (7, 10),
                },
            ],
        );
    }

    #[test]
    fn test_fasta_trailing_header_only() {
        assert_parse_with_all_settings(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n>Virus3",
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".to_string(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".to_string(),
                    line_range: (4, 7),
                },
                FastaRecord {
                    header: "Virus3".to_string(),
                    contents: "".to_string(),
                    line_range: (7, 8),
                },
            ],
        );
    }

    #[test]
    fn test_fasta_header_concat() {
        assert_parse_with_concatenate_headers(
            ">a\n>b\ntest",
            true,
            vec![FastaRecord {
                header: "a\nb".to_string(),
                contents: "test".to_string(),
                line_range: (1, 4),
            }],
        );
    }

    #[test]
    fn test_fasta_header_noconcat() {
        assert_parse_with_concatenate_headers(
            ">a\n>b\ntest",
            false,
            vec![
                FastaRecord {
                    header: "a".to_string(),
                    contents: "".to_string(),
                    line_range: (1, 2),
                },
                FastaRecord {
                    header: "b".to_string(),
                    contents: "test".to_string(),
                    line_range: (2, 4),
                },
            ],
        );
    }

    #[test]
    fn test_fasta_header_concat_no_content() {
        assert_parse_with_concatenate_headers(
            ">a\n>b",
            true,
            vec![FastaRecord {
                header: "a\nb".to_string(),
                contents: "".to_string(),
                line_range: (1, 3),
            }],
        );
    }

    #[test]
    fn test_fasta_unicode() {
        // U+2009 is THIN-SPACE, which is White_Space=yes so should be ignored even if allow_comment is false
        // Note that some "spacey" unicode characters, like U+200B "ZERO WIDTH SPACE", are White_Space=no
        assert_parse_with_all_settings("\u{2009}\r\n>ὦ Ᾰ̓θηνᾶ, Heizölrückstoßabdämpfungを持つ!\nPchnąć w tę łódź jeża lub ośm skrzyń fig", vec![FastaRecord {
            header: "ὦ Ᾰ̓θηνᾶ, Heizölrückstoßabdämpfungを持つ!".to_string(),
            contents: "Pchnąć w tę łódź jeża lub ośm skrzyń fig".to_string(),
            line_range: (2, 4),
        }])
    }

    #[test]
    fn test_fasta_ignores_alternative_linebreaks() {
        use std::fmt::Write;

        let bad_linebreaks = [
            '\x0B',     // LINE TABULATION  (\v)
            '\x0C',     // FORM FEED (\f)
            '\x0D',     // CARRIAGE RETURN (\r)
            '\u{0085}', // NEXT LINE
            '\u{2028}', // LINE SEPARATOR
            '\u{2029}', // PARAGRAPH SEPARATOR
        ];
        let mut test_case = "contents-should-be-single-line:".to_string();
        for lb in bad_linebreaks {
            write!(
                test_case,
                " {}>U+{:x} should not be a header",
                lb, lb as u32
            )
            .unwrap();
        }

        assert_parse_with_allow_preceding_comment(
            &test_case.clone(),
            false,
            vec![FastaRecord {
                header: "".to_string(),
                contents: test_case,
                line_range: (1, 2),
            }],
        )
    }

    #[test]
    fn test_fasta_handles_windows_newline() {
        // the Rust BufRead trait's `lines()` method should handle \r\n lines by removing the \r as well,
        // so this just checks that we're not doing anything to mess that up

        assert_parse_with_all_settings(
            "\r\n>i love compatability\r\nwindows is awesome\r\n",
            vec![FastaRecord {
                header: "i love compatability".to_string(),
                contents: "windows is awesome".to_string(),
                line_range: (2, 4),
            }],
        )
    }

    #[test]
    fn test_concat_headers_not_quadratic() {
        // if this test gets slow, we probably accidentally introduced quadratic behavior somewhere

        let mut test_case = String::new();
        for _ in 0..10_000 {
            test_case.push_str(">header\n");
        }
        let header = test_case.trim().replace('>', "");

        panic_after(Duration::from_millis(200), move || {
            assert_parse_with_concatenate_headers(
                &test_case,
                true,
                vec![FastaRecord {
                    header,
                    contents: "".to_string(),
                    line_range: (1, 10_001),
                }],
            )
        });
    }

    #[test]
    fn test_comment_not_quadratic() {
        // if this test gets slow, we probably accidentally introduced quadratic behavior somewhere

        let mut test_case = String::new();
        for _ in 0..10_000 {
            test_case.push_str("comment\n");
        }
        let contents = test_case.replace('\n', "");

        panic_after(Duration::from_millis(200), move || {
            assert_parse_with_allow_preceding_comment(
                &test_case,
                false,
                vec![FastaRecord {
                    header: "".to_string(),
                    contents,
                    line_range: (1, 10_001),
                }],
            )
        });
    }

    #[test]
    fn test_contents_not_quadratic() {
        // if this test gets slow, we probably accidentally introduced quadratic behavior somewhere

        let mut contents = String::new();
        for _ in 0..10_000 {
            contents.push_str("contents\n");
        }
        let test_case = format!(">header\n{contents}");
        let contents = contents.replace('\n', "");

        panic_after(Duration::from_millis(200), move || {
            assert_parse_with_all_settings(
                &test_case,
                vec![FastaRecord {
                    header: "header".to_string(),
                    contents,
                    line_range: (1, 10_002),
                }],
            )
        });
    }

    #[test]
    fn test_dna_fasta() {
        assert_parse!(
            ">Virus1\nAAAA",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            vec![FastaRecord {
                header: "Virus1".to_string(),
                contents: "AAAA".parse().unwrap(),
                line_range: (1, 3),
            }]
        );
    }

    #[test]
    fn test_dna_fasta_ambiguous() {
        assert_parse!(
            ">Virus1\nABCD",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            vec![FastaRecord {
                header: "Virus1".to_string(),
                contents: "ABCD".parse().unwrap(),
                line_range: (1, 3),
            }]
        );
    }

    #[test]
    fn test_dna_fasta_strict() {
        // Strict as in "when parsing Nucleotide, disallow ambiguity codes".
        // No relation to the FastaParser `strict()` call.
        assert_parse_err!(
            ">Virus1\nABCD",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            FastaParseError::ParseError(TranslationError::UnexpectedAmbiguousNucleotide('B'))
        );
    }

    #[test]
    fn test_dna_fasta_multiple() {
        assert_parse!(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".parse().unwrap(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".parse().unwrap(),
                    line_range: (4, 7),
                },
            ]
        );
        assert_parse!(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nRRRR\n",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".parse().unwrap(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCRRRR".parse().unwrap(),
                    line_range: (4, 7),
                },
            ]
        );
    }

    #[test]
    fn test_dna_dna_whitespace() {
        assert_parse!(
            ">Virus1\nAAAA \n",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            vec![FastaRecord {
                header: "Virus1".to_string(),
                contents: "AAAA".parse().unwrap(),
                line_range: (1, 3),
            }]
        );

        assert_parse!(
            ">Virus1\n  AAAA\tBCD \t\n",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            vec![FastaRecord {
                header: "Virus1".to_string(),
                contents: "AAAABCD".parse().unwrap(),
                line_range: (1, 3),
            }]
        );
    }

    #[test]
    fn test_dna_invalid_dna() {
        assert_parse_err!(
            ">Virus1\nAAAelephant",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            FastaParseError::ParseError(TranslationError::BadNucleotide('e'))
        );
        assert_parse_err!(
            ">Virus1\nAAAelephant",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            FastaParseError::ParseError(TranslationError::BadNucleotide('e'))
        );
    }

    #[test]
    fn test_dna_invalid_dna_multiple() {
        assert_parse_err!(
            ">Virus1\nAAAA\n>Virus2\nAAAAelephant",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            FastaParseError::ParseError(TranslationError::BadNucleotide('e'))
        );
        assert_parse_err!(
            ">Virus1\nAAAA\n>Virus2\nAAAAelephant",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            FastaParseError::ParseError(TranslationError::BadNucleotide('e'))
        );
    }

    #[test]
    fn test_dna_invalid_dna_unicode() {
        assert_parse_err!(
            ">Virus1\nAAčCCG\n",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            FastaParseError::ParseError(TranslationError::NonAsciiByte(196))
        );
        assert_parse_err!(
            ">Virus1\nAAčCCG\n",
            FastaParser::<DnaSequence<NucleotideAmbiguous>>::strict(),
            FastaParseError::ParseError(TranslationError::NonAsciiByte(196))
        );
    }

    #[test]
    fn test_protein_fasta() {
        assert_parse!(
            ">Virus1\nAAAA",
            FastaParser::<ProteinSequence>::strict(),
            vec![FastaRecord {
                header: "Virus1".to_string(),
                contents: "AAAA".parse().unwrap(),
                line_range: (1, 3),
            }]
        );
    }

    #[test]
    fn test_protein_fasta_multiple() {
        assert_parse!(
            ">Virus1\nAAAA\nAAAA\n>Virus2\nCCCC\nCCCC\n",
            FastaParser::<ProteinSequence>::strict(),
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".parse().unwrap(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus2".to_string(),
                    contents: "CCCCCCCC".parse().unwrap(),
                    line_range: (4, 7),
                },
            ]
        );
    }

    #[test]
    fn test_duplicate_header_lines() {
        assert_parse!(
            ">Virus1\nAAAA\nAAAA\n>Virus1\nCCCC\nCCCC\n",
            FastaParser::<DnaSequence<Nucleotide>>::strict(),
            vec![
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "AAAAAAAA".parse().unwrap(),
                    line_range: (1, 4),
                },
                FastaRecord {
                    header: "Virus1".to_string(),
                    contents: "CCCCCCCC".parse().unwrap(),
                    line_range: (4, 7),
                },
            ]
        );
    }

    #[test]
    fn test_to_string() {
        let parser = FastaParser::<DnaSequence<Nucleotide>>::lax();
        let string = ">Virus1\nAC\nT\n>Empty\n\n>Virus2\n>with many\n>comment lines\nC  AT";
        let parsed = parser.parse_str(string).unwrap();

        // Test: if we to_string the parsed file and parse it again, we should
        // get the same records again, ignoring line_range.
        let restrung = parsed.to_string();
        let reparsed = parser.parse_str(&restrung).unwrap();

        assert_eq!(parsed.records.len(), 3);
        assert_eq!(reparsed.records.len(), 3);

        // Compare all records, but ignore line_range:
        for i in 0..3 {
            assert_eq!(parsed.records[i].header, reparsed.records[i].header);
            assert_eq!(parsed.records[i].contents, reparsed.records[i].contents);
        }
    }

    // TODO: when we add validation for ProteinSequence, add tests for that here
}

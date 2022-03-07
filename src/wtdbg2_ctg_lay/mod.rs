use std::cmp::Ordering;
use std::num::NonZeroU64;
use std::str::FromStr;

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Wtdbg2CtgLayLineWithContext {
    pub line: Wtdbg2CtgLayLine,
    pub context: LineContext,
}

#[derive(Clone, Eq, PartialEq, Debug)]
pub enum Wtdbg2CtgLayLine {
    Contig {
        name: String,
        node_count: u64,
        length: u64,
    },

    Edge {
        offset: u64,
        from_node: String,
        /// True for forwards (+), false for backwards (-).
        from_direction: bool,
        to_node: String,
        /// True for forwards (+), false for backwards (-).
        to_direction: bool,
    },

    Alignment {
        read_id: Vec<u8>,
        /// True for forwards (+), false for backwards (-).
        direction: bool,
        offset: usize,
        length: usize,
        original_length: usize,
    },
}

impl FromStr for Wtdbg2CtgLayLine {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.chars().next() {
            Some('>') => {
                let mut columns = s[1..].split('\t');
                let name = columns.next().ok_or(())?.to_owned();
                let node_count = columns.next().ok_or(())?[6..].parse().map_err(|_| ())?;
                let length = columns.next().ok_or(())?[4..].parse().map_err(|_| ())?;
                Ok(Self::Contig {
                    name,
                    node_count,
                    length,
                })
            }
            Some('E') => {
                let mut columns = s[1..].split('\t');
                columns.next().ok_or(())?;
                let offset = columns.next().ok_or(())?.parse().map_err(|_| ())?;
                let from_node = columns.next().ok_or(())?.to_owned();
                let from_direction = match columns.next().ok_or(())? {
                    "+" => true,
                    "-" => false,
                    _ => return Err(()),
                };
                let to_node = columns.next().ok_or(())?.to_owned();
                let to_direction = match columns.next().ok_or(())? {
                    "+" => true,
                    "-" => false,
                    _ => return Err(()),
                };
                Ok(Self::Edge {
                    offset,
                    from_node,
                    from_direction,
                    to_node,
                    to_direction,
                })
            }
            Some('S') => {
                let mut columns = s[1..].split('\t');
                columns.next().ok_or(())?;
                let read_id = columns.next().ok_or(())?.as_bytes().to_owned();
                let direction = match columns.next().ok_or(())? {
                    "+" => true,
                    "-" => false,
                    _ => return Err(()),
                };
                let offset = columns.next().ok_or(())?.parse().map_err(|_| ())?;
                let length = columns.next().ok_or(())?.parse().map_err(|_| ())?;
                Ok(Self::Alignment {
                    read_id,
                    direction,
                    offset,
                    length,
                    original_length: length,
                })
            }
            _ => Err(()),
        }
    }
}

impl PartialOrd for LineContext {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for LineContext {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.contig_index.cmp(&other.contig_index) {
            Ordering::Equal => {
                assert_eq!(
                    self.previous_contig_edge_count,
                    other.previous_contig_edge_count
                );
                match self.edge_index.cmp(&other.edge_index) {
                    Ordering::Equal => {
                        assert_eq!(
                            self.previous_edge_alignment_count,
                            other.previous_edge_alignment_count
                        );
                        self.alignment_index.cmp(&other.alignment_index)
                    }
                    ordering => ordering,
                }
            }
            ordering => ordering,
        }
    }
}

impl ToString for Wtdbg2CtgLayLine {
    fn to_string(&self) -> String {
        match self {
            Wtdbg2CtgLayLine::Contig {
                name,
                node_count,
                length,
            } => format!(">{name}\tnodes={node_count}\tlen={length}"),
            Wtdbg2CtgLayLine::Edge {
                offset,
                from_node,
                from_direction,
                to_node,
                to_direction,
            } => {
                let from_direction = if *from_direction { "+" } else { "-" };
                let to_direction = if *to_direction { "+" } else { "-" };
                format!("E\t{offset}\t{from_node}\t{from_direction}\t{to_node}\t{to_direction}")
            }
            Wtdbg2CtgLayLine::Alignment {
                read_id,
                direction,
                offset,
                length,
                ..
            } => {
                let read_id = String::from_utf8(read_id.clone()).unwrap();
                let direction = if *direction { "+" } else { "-" };
                format!("S\t{read_id}\t{direction}\t{offset}\t{length}\t")
            }
        }
    }
}

#[derive(Clone, Eq, PartialEq, Debug, Default)]
pub struct LineContext {
    pub contig_index: u64,
    pub edge_index: u64,
    pub alignment_index: u64,
    pub previous_contig_edge_count: Option<NonZeroU64>,
    pub previous_edge_alignment_count: Option<NonZeroU64>,
}

impl LineContext {
    pub fn directly_precedes(&self, other: &Self) -> bool {
        if self.contig_index != other.contig_index {
            if self.contig_index != other.contig_index - 1 {
                return false;
            }

            self.edge_index == other.previous_contig_edge_count.unwrap().get()
                && self.alignment_index == other.previous_edge_alignment_count.unwrap().get()
        } else if self.edge_index != other.edge_index {
            if self.edge_index != other.edge_index - 1 {
                return false;
            }

            self.alignment_index == other.previous_edge_alignment_count.unwrap().get()
        } else if self.alignment_index != other.alignment_index {
            if self.alignment_index != other.alignment_index - 1 {
                return false;
            }

            true
        } else {
            false
        }
    }
}

use std::path::PathBuf;
use log::{info, LevelFilter};
use simplelog::{ColorChoice, TerminalMode, TermLogger};
use clap::Parser;
use crate::fasta_sequence_index::FastaSequenceIndex;

mod fasta_sequence_index;

#[derive(Parser)]
struct Configuration {
    /// The input file. Must be in wtdbg2's .ctg.lay format.
    #[clap(long, parse(from_os_str))]
    input: PathBuf,

    /// The output file. Must be in wtdbg2's .ctg.lay format.
    #[clap(long, parse(from_os_str))]
    output: PathBuf,

    /// A fasta file containing the normal reads.
    #[clap(long, parse(from_os_str))]
    normal_reads: PathBuf,

    /// A fasta file containing the homopolymer-compressed reads.
    #[clap(long, parse(from_os_str))]
    hoco_reads: PathBuf,

    /// The size of the buffers between threads.
    #[clap(long, default_value = "32768")]
    buffer_size: usize,

    /// The number of compute threads to use for decompression.
    /// Note that the input and output threads are not counted under this number.
    #[clap(long, default_value = "1")]
    threads: usize,
}

fn initialise_logging() {
    TermLogger::init(
        LevelFilter::Debug,
        Default::default(),
        TerminalMode::Stderr,
        ColorChoice::Auto,
    )
        .unwrap();
    info!("Logging initialised successfully")
}

fn main() {
    let configuration = Configuration::parse();
    initialise_logging();

    let mut normal_sequence_index_path = configuration.output.clone().into_os_string();
    normal_sequence_index_path.push(".normal_index");
    let normal_sequence_index_path = normal_sequence_index_path;
    let mut hoco_sequence_index_path = configuration.output.clone().into_os_string();
    hoco_sequence_index_path.push(".hoco_index");
    let hoco_sequence_index_path = hoco_sequence_index_path;

    info!("Building reads sequence indices...");
    let (normal_sequence_index, hoco_sequence_index) = crossbeam::scope(|scope| {
        let normal_sequence_index = scope.builder().name("normal_index_builder_thread".to_string()).spawn(|scope| {
            FastaSequenceIndex::build_parallel(configuration.normal_reads, normal_sequence_index_path, scope, configuration.buffer_size)
        }).unwrap();
        let hoco_sequence_index = scope.builder().name("hoco_index_builder_thread".to_string()).spawn(|scope| {
            FastaSequenceIndex::build_parallel(configuration.hoco_reads, hoco_sequence_index_path, scope, configuration.buffer_size)
        }).unwrap();

        (normal_sequence_index.join().unwrap(), hoco_sequence_index.join().unwrap())
    }).unwrap();
    info!("Built reads sequence indices");

    crossbeam::scope(|scope| {
        // TODO
    }).unwrap();
}

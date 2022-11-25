use crate::decompress::decompress;
use crate::fasta_sequence_index::FastaSequenceIndex;
use crate::wtdbg2_ctg_lay::{LineContext, Wtdbg2CtgLayLine, Wtdbg2CtgLayLineWithContext};
use clap::Parser;
use crossbeam::channel;
use log::{info, trace, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use std::collections::BTreeMap;
use std::fs;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::str::FromStr;

mod decompress;
mod fasta_sequence_index;
mod wtdbg2_ctg_lay;

#[derive(Parser, Clone)]
struct Configuration {
    /// The input file. Must be in wtdbg2's .ctg.lay format.
    #[clap(long, parse(from_os_str))]
    input: PathBuf,

    /// The output file. Must be in wtdbg2's .ctg.lay format.
    #[clap(long, parse(from_os_str))]
    output: PathBuf,

    /// A fasta file containing the normal (uncompressed) reads.
    #[clap(long, parse(from_os_str))]
    normal_reads: PathBuf,

    /// The size of the queues between threads.
    #[clap(long, default_value = "32768")]
    queue_size: usize,

    /// The size of the I/O buffers in bytes.
    #[clap(long, default_value = "67108864")]
    io_buffer_size: usize,

    /// The number of compute threads to use for decompression.
    /// Note that the input and output threads are not counted under this number.
    #[clap(long, default_value = "1")]
    compute_threads: usize,

    /// The level of log messages to be produced.
    #[clap(long, default_value = "Info")]
    log_level: LevelFilter,
}

fn initialise_logging(log_level: &LevelFilter) {
    TermLogger::init(
        *log_level,
        Default::default(),
        TerminalMode::Stderr,
        ColorChoice::Auto,
    )
    .unwrap();
    info!("Logging initialised successfully")
}

fn main() {
    let configuration = Configuration::parse();
    initialise_logging(&configuration.log_level);

    let mut normal_sequence_index_path = configuration.output.clone().into_os_string();
    normal_sequence_index_path.push(".normal_index");
    let normal_sequence_index_path = normal_sequence_index_path;
    let mut contig_tmp_path = configuration.output.clone().into_os_string();
    contig_tmp_path.push(".current_contig");
    let contig_tmp_path = contig_tmp_path;
    // Create/open the files here already to abort early if it cannot be created.
    let tmp_file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(&contig_tmp_path)
        .unwrap();
    let input_file = File::open(&configuration.input).unwrap();
    let output_file = File::create(&configuration.output).unwrap();

    info!("Building reads sequence indices...");
    // parallel builds seem to be a little faster on my laptop.
    let mut normal_sequence_index = {
        let configuration = configuration.clone();
        crossbeam::scope(|scope| {
            let normal_sequence_index = scope
                .builder()
                .name("normal_index_builder_thread".to_string())
                .spawn(|scope| {
                    FastaSequenceIndex::build_parallel(
                        configuration.normal_reads,
                        normal_sequence_index_path,
                        scope,
                        configuration.queue_size,
                        configuration.io_buffer_size,
                    )
                    //FastaSequenceIndex::build(configuration.normal_reads, normal_sequence_index_path, configuration.io_buffer_size)
                })
                .unwrap();

            normal_sequence_index.join().unwrap()
        })
        .unwrap()
    };
    info!("Built read sequence indices");

    info!("Decompressing...");
    crossbeam::scope(|scope| {
        // Read input file.
        let (input_sender, input_receiver) = channel::bounded(configuration.queue_size);
        scope
            .builder()
            .name("input_reader".to_string())
            .spawn(move |_| {
                let mut line_number = 0;
                for line in
                    BufReader::with_capacity(configuration.io_buffer_size, input_file).lines()
                {
                    let line = line.unwrap();
                    line_number += 1;
                    trace!("Read line {line_number}");
                    input_sender.send(line).unwrap();
                }
            })
            .unwrap();

        // Parse input lines.
        let (alignment_sender, alignment_receiver) = channel::bounded(configuration.queue_size);
        let (decompressed_alignment_sender, decompressed_alignment_receiver) =
            channel::bounded(configuration.queue_size);
        {
            let decompressed_alignment_sender = decompressed_alignment_sender.clone();
            scope
                .builder()
                .name("input_parser".to_string())
                .spawn(move |_| {
                    let mut context = LineContext::default();

                    while let Ok(input) = input_receiver.recv() {
                        let line = Wtdbg2CtgLayLine::from_str(&input)
                            .unwrap_or_else(|_| panic!("Could not parse: {input}"));
                        match line {
                            Wtdbg2CtgLayLine::Contig { .. } => {
                                trace!("Parsed contig line {input}");
                                if context.contig_index != -1 {
                                    context.previous_contig_edge_count = context.edge_index + 1;
                                    context.previous_edge_alignment_count =
                                        context.alignment_index + 1;
                                }
                                context.contig_index += 1;
                                context.edge_index = -1;
                                context.alignment_index = -1;
                                decompressed_alignment_sender
                                    .send((
                                        Wtdbg2CtgLayLineWithContext {
                                            line,
                                            context: context.clone(),
                                        },
                                        None,
                                    ))
                                    .unwrap();
                            }
                            Wtdbg2CtgLayLine::Edge { .. } => {
                                trace!("Parsed edge line {input}");
                                assert!(context.contig_index >= 0);
                                if context.edge_index != -1 {
                                    context.previous_edge_alignment_count =
                                        context.alignment_index + 1;
                                }
                                context.edge_index += 1;
                                context.alignment_index = -1;
                                decompressed_alignment_sender
                                    .send((
                                        Wtdbg2CtgLayLineWithContext {
                                            line,
                                            context: context.clone(),
                                        },
                                        None,
                                    ))
                                    .unwrap();
                            }
                            Wtdbg2CtgLayLine::Alignment { .. } => {
                                trace!("Parsed alignment line");
                                assert!(context.contig_index >= 0);
                                assert!(context.edge_index >= 0);
                                context.alignment_index += 1;
                                alignment_sender
                                    .send(Wtdbg2CtgLayLineWithContext {
                                        line,
                                        context: context.clone(),
                                    })
                                    .unwrap();
                            }
                        }
                    }
                })
                .unwrap();
        }

        // Decorate alignments with read sequences.
        let (decorated_alignment_sender, decorated_alignment_receiver) =
            channel::bounded(configuration.queue_size);
        scope
            .builder()
            .name("read_sequence_reader".to_string())
            .spawn(move |_| {
                while let Ok(line_with_context) = alignment_receiver.recv() {
                    let read_id = match &line_with_context.line {
                        Wtdbg2CtgLayLine::Alignment { read_id, .. } => read_id,
                        _ => unreachable!("Not an alignment: {line_with_context:?}"),
                    };
                    let mut sequence = Vec::new();
                    let read_id_string = String::from_utf8(read_id.clone()).unwrap();
                    trace!("Reading read {read_id_string}");
                    normal_sequence_index.get_sequence(read_id, &mut sequence);
                    decorated_alignment_sender
                        .send((line_with_context, sequence))
                        .unwrap();
                }
            })
            .unwrap();

        // Decompress.
        for thread_index in 0..configuration.compute_threads {
            let decorated_alignment_receiver = decorated_alignment_receiver.clone();
            let decompressed_alignment_sender = decompressed_alignment_sender.clone();
            scope
                .builder()
                .name(format!("decompressor_{thread_index}"))
                .spawn(move |_| {
                    while let Ok((
                        Wtdbg2CtgLayLineWithContext {
                            line:
                                Wtdbg2CtgLayLine::Alignment {
                                    read_id,
                                    direction,
                                    offset,
                                    length,
                                    original_length,
                                },
                            context,
                        },
                        sequence,
                    )) = decorated_alignment_receiver.recv()
                    {
                        trace!("Decompressing {context:?}");
                        let limit = offset + length;
                        let (shifted_offset, shifted_limit) = decompress(offset, limit, &sequence);
                        let shifted_length = shifted_limit - shifted_offset;
                        let shifted_sequence = &sequence[shifted_offset..shifted_limit];
                        decompressed_alignment_sender
                            .send((
                                Wtdbg2CtgLayLineWithContext {
                                    line: Wtdbg2CtgLayLine::Alignment {
                                        read_id,
                                        direction,
                                        offset: shifted_offset,
                                        length: shifted_length,
                                        original_length,
                                    },
                                    context,
                                },
                                Some(if direction {
                                    shifted_sequence.to_owned()
                                } else {
                                    reverse_complement(shifted_sequence.iter().cloned())
                                }),
                            ))
                            .unwrap();
                    }
                })
                .unwrap();
        }

        // Decompression with multiple threads will shuffle the alignments a bit, so we need to put them back into order.
        let (output_sender, output_receiver) = channel::bounded(configuration.queue_size);
        scope
            .builder()
            .name("sorter".to_owned())
            .spawn(move |_| {
                let mut current_context = LineContext::default();
                let mut sorted_lines = BTreeMap::new();
                let mut alignment_count = 0;
                let mut original_alignment_length_sum = 0;
                let mut shifted_alignment_length_sum = 0;
                let mut original_previous_offset = 0;
                let mut shifted_previous_offset = 0;

                while let Ok((Wtdbg2CtgLayLineWithContext { line, context }, shifted_sequence)) =
                    decompressed_alignment_receiver.recv()
                {
                    trace!("Received {context:?}");
                    assert!(sorted_lines
                        .insert(context, (line, shifted_sequence))
                        .is_none());

                    while let Some(context) = sorted_lines.keys().next().cloned() {
                        trace!(
                            "Last context is {current_context:?}, and next known is {context:?}"
                        );
                        if current_context.directly_precedes(&context) {
                            let (mut line, shifted_sequence) =
                                sorted_lines.remove(&context).unwrap();
                            match &mut line {
                                Wtdbg2CtgLayLine::Contig { .. } => {
                                    alignment_count = 0;
                                    original_alignment_length_sum = 0;
                                    shifted_alignment_length_sum = 0;
                                    original_previous_offset = 0;
                                    shifted_previous_offset = 0;
                                    assert!(shifted_sequence.is_none());
                                    output_sender.send((line, None)).unwrap()
                                }
                                Wtdbg2CtgLayLine::Edge { offset, .. } => {
                                    let original_offset = *offset;
                                    *offset = shifted_previous_offset
                                        + ((*offset - original_previous_offset) as f64
                                            * shifted_alignment_length_sum as f64
                                            / original_alignment_length_sum as f64)
                                            .round()
                                            as u64;
                                    alignment_count = 0;
                                    original_alignment_length_sum = 0;
                                    shifted_alignment_length_sum = 0;
                                    original_previous_offset = original_offset;
                                    shifted_previous_offset = *offset;
                                    assert!(shifted_sequence.is_none());
                                    output_sender.send((line, None)).unwrap()
                                }
                                Wtdbg2CtgLayLine::Alignment {
                                    length,
                                    original_length,
                                    ..
                                } => {
                                    alignment_count += 1;
                                    original_alignment_length_sum += *original_length;
                                    shifted_alignment_length_sum += *length;
                                    let estimated_length = (shifted_alignment_length_sum as f64
                                        / alignment_count as f64)
                                        .round()
                                        as u64;
                                    assert!(shifted_sequence.is_some());
                                    output_sender
                                        .send((
                                            line,
                                            shifted_sequence.map(|shifted_sequence| {
                                                (shifted_sequence, estimated_length)
                                            }),
                                        ))
                                        .unwrap();
                                }
                            }
                            current_context = context;
                        } else {
                            break;
                        }
                    }
                }
            })
            .unwrap();

        // Write output.
        scope
            .builder()
            .name("output_writer".to_owned())
            .spawn(move |_| {
                let mut output_writer =
                    BufWriter::with_capacity(configuration.io_buffer_size, output_file);
                let mut tmp_writer =
                    BufWriter::with_capacity(configuration.io_buffer_size, tmp_file);
                let mut append_file_buffer = vec![0; configuration.io_buffer_size];
                let mut current_offset = 0;
                let mut current_last_edge_length = 0;
                let mut current_contig_line = None;
                while let Ok((mut line, sequence_and_length)) = output_receiver.recv() {
                    trace!("Writing line {line:?}");
                    match &mut line {
                        Wtdbg2CtgLayLine::Contig { .. } => {
                            tmp_writer = finalise_contig_line(
                                &mut current_contig_line,
                                &mut current_offset,
                                &mut current_last_edge_length,
                                &mut output_writer,
                                tmp_writer,
                                &mut append_file_buffer,
                                &configuration,
                            );

                            current_contig_line = Some(line);
                        }
                        Wtdbg2CtgLayLine::Edge { offset, .. } => {
                            current_offset = *offset;
                            tmp_writer.write_all(line.to_string().as_bytes()).unwrap();
                            tmp_writer.write_all(&[b'\n']).unwrap();
                        }
                        Wtdbg2CtgLayLine::Alignment { .. } => {
                            let (sequence, edge_length) = sequence_and_length.unwrap();
                            current_last_edge_length = edge_length;

                            tmp_writer.write_all(line.to_string().as_bytes()).unwrap();
                            tmp_writer.write_all(&sequence).unwrap();
                            tmp_writer.write_all(&[b'\n']).unwrap();
                        }
                    }
                }

                finalise_contig_line(
                    &mut current_contig_line,
                    &mut current_offset,
                    &mut current_last_edge_length,
                    &mut output_writer,
                    tmp_writer,
                    &mut append_file_buffer,
                    &configuration,
                );
            })
            .unwrap();
    })
    .unwrap();

    // Remove tmp file as it is not needed anymore.
    fs::remove_file(&contig_tmp_path).unwrap();

    info!("Done");
}

fn finalise_contig_line<OutputWriter: Write>(
    current_contig_line: &mut Option<Wtdbg2CtgLayLine>,
    current_offset: &mut u64,
    current_last_edge_length: &mut u64,
    output_writer: &mut OutputWriter,
    mut tmp_writer: BufWriter<File>,
    append_file_buffer: &mut [u8],
    configuration: &Configuration,
) -> BufWriter<File> {
    if let Some(mut current_contig_line) = current_contig_line.take() {
        match &mut current_contig_line {
            Wtdbg2CtgLayLine::Contig { length, .. } => {
                *length = *current_offset + *current_last_edge_length;
                *current_offset = 0;
                *current_last_edge_length = 0;
                output_writer
                    .write_all(current_contig_line.to_string().as_bytes())
                    .unwrap();
                output_writer.write_all(&[b'\n']).unwrap();

                // Append the tmp file to the actual file, now that we know how long the decompressed contig is.
                let mut tmp_file = tmp_writer.into_inner().unwrap();
                tmp_file.seek(SeekFrom::Start(0)).unwrap();
                loop {
                    let length = tmp_file.read(append_file_buffer).unwrap();
                    if length > 0 {
                        output_writer
                            .write_all(&append_file_buffer[..length])
                            .unwrap();
                    } else {
                        break;
                    }
                }
                tmp_file.set_len(0).unwrap();
                tmp_file.seek(SeekFrom::Start(0)).unwrap();
                tmp_writer = BufWriter::with_capacity(configuration.io_buffer_size, tmp_file);
            }
            _ => unreachable!("Contig line is not a contig line {current_contig_line:?}"),
        }
    }

    tmp_writer
}

fn reverse_complement<
    IntoIter: DoubleEndedIterator<Item = u8>,
    DnaIterator: IntoIterator<Item = u8, IntoIter = IntoIter>,
>(
    dna: DnaIterator,
) -> Vec<u8> {
    dna.into_iter()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => panic!("Unknown dna character: {other}"),
        })
        .rev()
        .collect()
}

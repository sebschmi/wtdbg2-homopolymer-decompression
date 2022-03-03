use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::slice;
use bio::io::fasta;

struct FileSlice {
    offset: u64,
    len: usize,
}

pub struct FastaSequenceIndex {
    file: File,
    index: HashMap<Vec<u8>, FileSlice>,
}

impl FastaSequenceIndex {
    pub fn build<P1: AsRef<Path>, P2: AsRef<Path>>(input_file: P1, tmp_file: P2) -> Self {
        let reader =
            fasta::Reader::with_capacity(64 * 1024 * 1024, File::open(input_file).unwrap());
        let mut writer = BufWriter::with_capacity(64 * 1024 * 1024, File::create(tmp_file).unwrap());

        let mut index = HashMap::new();

        let mut offset = 0;
        for record in reader.records() {
            let record = record.unwrap();

            let len = record.seq().len();
            writer.write_all(record.seq()).unwrap();
            // Write delimiter character to catch errors.
            writer.write_all(&[b'\n']).unwrap();

            let previous = index.insert(record.id().as_bytes().to_vec(), FileSlice {offset, len});
            offset += u64::try_from(len).unwrap() + 1;
            assert!(previous.is_none());
        }

        Self {
            file: writer.into_inner().unwrap(),
            index,
        }
    }



    pub fn get_sequence(&mut self, id: &Vec<u8>, output: &mut Vec<u8>) {
        let file_slice = self.index.get(id).unwrap();
        output.clear();
        output.reserve(file_slice.len);

        let buffer = output.as_mut_ptr();
        let capacity = output.capacity();

        self.file.read_exact_at(unsafe {slice::from_raw_parts_mut(buffer, file_slice.len)}, file_slice.offset).unwrap();
        *output = unsafe {Vec::from_raw_parts(buffer, file_slice.len, capacity)};
    }
}
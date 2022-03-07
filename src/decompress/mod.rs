pub fn decompress(offset: usize, limit: usize, sequence: &[u8]) -> (usize, usize) {
    // Find offset.
    // Use a block to ensure the next block does not accidentally reuse any variable.
    let shifted_offset = {
        let mut shifted_offset = 0;
        let mut current_offset = 0;
        for character_window in sequence.windows(2) {
            if current_offset == offset {
                break;
            }

            // Safety: windows of size 2.
            if unsafe { character_window.get_unchecked(0) != character_window.get_unchecked(1) } {
                //if character_window.get(0).unwrap() != character_window.get(1).unwrap() {
                current_offset += 1;
            }
            shifted_offset += 1;
        }
        // The windowed iteration cannot recognise an offset after the end of the sequence.
        if current_offset != offset {
            shifted_offset += 1;
        }
        shifted_offset
    };

    // Find limit.
    let mut shifted_limit = shifted_offset;
    let mut current_limit = offset;
    for character_window in sequence.windows(2).skip(shifted_offset) {
        if current_limit == limit {
            break;
        }

        // Safety: windows of size 2.
        if unsafe { character_window.get_unchecked(0) != character_window.get_unchecked(1) } {
            //if character_window.get(0).unwrap() != character_window.get(1).unwrap() {
            current_limit += 1;
        }
        shifted_limit += 1;
    }
    // The windowed iteration cannot recognise a limit after the end of the sequence.
    if current_limit != limit {
        shifted_limit += 1;
    }
    (shifted_offset, shifted_limit)
}

#[cfg(test)]
mod tests {
    use crate::decompress;

    #[test]
    fn test_decompress() {
        let sequence = vec![0, 0, 1, 1, 2, 3, 3, 3, 4, 5];
        let tests = [
            (2, 4, 4, 8),
            (0, 0, 0, 0),
            (0, 1, 0, 2),
            (0, 2, 0, 4),
            (1, 2, 2, 4),
            (2, 2, 4, 4),
            (1, 3, 2, 5),
            (2, 3, 4, 5),
            (3, 4, 5, 8),
            (3, 5, 5, 9),
            (4, 5, 8, 9),
            (5, 5, 9, 9),
            (0, 6, 0, 10),
            (4, 6, 8, 10),
            (5, 6, 9, 10),
            (6, 6, 10, 10),
        ];
        for (offset, limit, shifted_offset, shifted_limit) in tests {
            let (decompressed_offset, decompressed_limit) = decompress(offset, limit, &sequence);
            assert_eq!((decompressed_offset, decompressed_limit), (shifted_offset, shifted_limit), "({offset}, {limit}): expected ({shifted_offset}, {shifted_limit}) but got ({decompressed_offset}, {decompressed_limit})");
        }
    }
}

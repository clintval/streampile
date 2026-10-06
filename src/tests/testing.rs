use noodles::bam;
use noodles::core::Position;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value;

use super::name;
use crate::testing::{Frag, Pair, READ_GROUP_ID, SAMPLE_NAME, SamBuilder, Strand, default_contigs};

fn start(record: &RecordBuf) -> Option<usize> {
    record.alignment_start().map(usize::from)
}

fn mate_start(record: &RecordBuf) -> Option<usize> {
    record.mate_alignment_start().map(usize::from)
}

fn tag(record: &RecordBuf, tag: Tag) -> Option<&Value> {
    record.data().get(&tag)
}

#[test]
fn test_a_builder_has_fgbios_contigs_and_read_group_in_coordinate_order() {
    let builder = SamBuilder::new();
    let header = builder.header();
    let contigs: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();
    assert_eq!(contigs, default_contigs());
    assert_eq!(
        (contigs[0].as_str(), contigs[24].as_str(), contigs.len()),
        ("chr1", "chrM", 25)
    );
    let text = {
        let mut writer = noodles::sam::io::Writer::new(Vec::new());
        writer.write_header(header).unwrap();
        String::from_utf8(writer.into_inner()).unwrap()
    };
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"), "{text}");
    assert!(
        text.contains(&format!("@RG\tID:{READ_GROUP_ID}\tSM:{SAMPLE_NAME}\n")),
        "{text}"
    );
}

/// fgbio's `SamBuilder.addPair`, with mate fields as htsjdk's `SamPairUtil.setMateInfo` sets
/// them.
#[test]
fn test_add_pair_sets_mate_information_as_htsjdk_does() {
    let mut builder = SamBuilder::new().read_length(40);
    let pair = builder.add_pair(Pair::at(100, 140));
    let (r1, r2) = (&pair[0], &pair[1]);
    assert_eq!((r1.template_length(), r2.template_length()), (80, -80));
    assert_eq!((mate_start(r1), mate_start(r2)), (start(r2), start(r1)));
    assert_eq!(
        (
            r1.mate_reference_sequence_id(),
            r2.mate_reference_sequence_id()
        ),
        (Some(0), Some(0))
    );
    let paired = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED;
    assert_eq!(
        r1.flags(),
        paired | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED
    );
    assert_eq!(
        r2.flags(),
        paired | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED
    );
    for record in [r1, r2] {
        assert_eq!(tag(record, Tag::MATE_CIGAR), Some(&Value::from("40M")));
        assert_eq!(
            tag(record, Tag::MATE_MAPPING_QUALITY),
            Some(&Value::from(60_u8))
        );
        assert_eq!(
            tag(record, Tag::READ_GROUP),
            Some(&Value::from(READ_GROUP_ID))
        );
        assert_eq!(record.mapping_quality(), MappingQuality::new(60));
        assert_eq!(
            (record.sequence().len(), record.quality_scores().len()),
            (40, 40)
        );
        assert!(
            record
                .quality_scores()
                .as_ref()
                .iter()
                .all(|&quality| quality == 30)
        );
    }
    assert_eq!(builder.records().len(), 2);
    assert_ne!(r1.sequence(), r2.sequence());
}

#[test]
fn test_add_pair_takes_the_reads_it_is_given() {
    let mut builder = SamBuilder::new();
    let pair = builder.add_pair(Pair {
        name: Some("q".into()),
        bases1: Some("ACGTA".into()),
        bases2: Some("TTGCA".into()),
        quals1: Some(vec![10, 20, 30, 40, 50]),
        cigar1: Some("2S3M".into()),
        cigar2: Some("3M2I".into()),
        contig2: Some(1),
        mapq1: 5,
        strand1: Strand::Minus,
        strand2: Strand::Minus,
        attrs: vec![(Tag::ALIGNMENT_HIT_COUNT, Value::from(2_u8))],
        ..Pair::at(10, 20)
    });
    let (r1, r2) = (&pair[0], &pair[1]);
    assert_eq!(r1.name().map(|name| name.to_vec()), Some(b"q".to_vec()));
    assert_eq!(r1.quality_scores().as_ref(), [10, 20, 30, 40, 50]);
    assert_eq!(r2.quality_scores().as_ref(), [30; 5]);
    assert_eq!(
        (r1.reference_sequence_id(), r2.reference_sequence_id()),
        (Some(0), Some(1))
    );
    assert_eq!(r1.template_length(), 0);
    assert!(!r1.flags().is_properly_segmented() && r1.flags().is_mate_reverse_complemented());
    assert_eq!(tag(r1, Tag::MATE_CIGAR), Some(&Value::from("3M2I")));
    assert_eq!(tag(r2, Tag::MATE_CIGAR), Some(&Value::from("2S3M")));
    assert_eq!(tag(r2, Tag::MATE_MAPPING_QUALITY), Some(&Value::from(5_u8)));
    assert_eq!(tag(r2, Tag::ALIGNMENT_HIT_COUNT), Some(&Value::from(2_u8)));
}

/// htsjdk's `SamPairUtil.setMateInfo` when one read is unmapped: it takes its mate's position,
/// and only it has the mate's mapping quality and CIGAR.
#[test]
fn test_an_unmapped_read_of_a_pair_takes_its_mates_position() {
    let mut builder = SamBuilder::new().read_length(10);
    let pair = builder.add_pair(Pair {
        unmapped2: true,
        ..Pair::at(100, 200)
    });
    let (mapped, unmapped) = (&pair[0], &pair[1]);
    assert!(unmapped.flags().is_unmapped() && mapped.flags().is_mate_unmapped());
    assert_eq!(
        (start(unmapped), mate_start(mapped), mate_start(unmapped)),
        (Some(100), Some(100), Some(100))
    );
    assert_eq!(unmapped.reference_sequence_id(), Some(0));
    assert_eq!(unmapped.cigar().as_ref(), []);
    assert_eq!(
        (
            tag(mapped, Tag::MATE_CIGAR),
            tag(mapped, Tag::MATE_MAPPING_QUALITY)
        ),
        (None, None)
    );
    assert_eq!(tag(unmapped, Tag::MATE_CIGAR), Some(&Value::from("10M")));
    assert_eq!(
        tag(unmapped, Tag::MATE_MAPPING_QUALITY),
        Some(&Value::from(60_u8))
    );
    assert_eq!(
        (mapped.template_length(), unmapped.template_length()),
        (0, 0)
    );
    assert!(!mapped.flags().is_properly_segmented());

    let both = builder.add_pair(Pair::default());
    for record in &both {
        assert!(record.flags().is_unmapped() && record.flags().is_mate_unmapped());
        assert_eq!((start(record), mate_start(record)), (None, None));
        assert_eq!(tag(record, Tag::MATE_CIGAR), None);
    }
}

#[test]
fn test_add_frag_builds_an_unpaired_read_with_sequential_names() {
    let mut builder = SamBuilder::new().read_length(8).base_quality(20);
    let first = builder.add_frag(Frag::at(5)).remove(0);
    let second = builder.add_frag(Frag {
        strand: Strand::Minus,
        cigar: Some("4M2D4M".into()),
        ..Frag::at(9)
    });
    assert_eq!(
        first.name().map(|name| name.to_vec()),
        Some(b"0000".to_vec())
    );
    assert_eq!(
        second[0].name().map(|name| name.to_vec()),
        Some(b"0001".to_vec())
    );
    assert_eq!(first.flags(), Flags::empty());
    assert_eq!(second[0].flags(), Flags::REVERSE_COMPLEMENTED);
    assert_eq!(first.quality_scores().as_ref(), [20; 8]);
    let placed = builder.add_frag(Frag {
        unmapped: true,
        ..Frag::at(7)
    });
    assert!(placed[0].flags().is_unmapped());
    assert_eq!(
        (start(&placed[0]), placed[0].mapping_quality()),
        (Some(7), MappingQuality::new(0))
    );
}

#[test]
#[should_panic(expected = "the bases of q do not agree with its CIGAR on length")]
fn test_a_read_whose_bases_disagree_with_its_cigar_is_refused() {
    SamBuilder::new().add_frag(Frag {
        name: Some("q".into()),
        bases: Some("ACGT".into()),
        cigar: Some("5M".into()),
        ..Frag::at(1)
    });
}

#[test]
fn test_records_come_out_in_coordinate_order_with_unplaced_reads_last() {
    let mut builder = SamBuilder::new().read_length(10);
    builder.add_frag(Frag {
        name: Some("unplaced".into()),
        ..Frag::default()
    });
    builder.add_frag(Frag {
        name: Some("chr2".into()),
        contig: 1,
        ..Frag::at(5)
    });
    builder.add_pair(Pair {
        name: Some("pair".into()),
        ..Pair::at(50, 20)
    });
    builder.add_frag(Frag {
        name: Some("late".into()),
        ..Frag::at(20)
    });
    let names: Vec<String> = builder.to_bam_records().iter().map(name).collect();
    assert_eq!(names, ["pair", "late", "pair", "chr2", "unplaced"]);
    let added: Vec<String> = builder
        .records()
        .iter()
        .map(|record| String::from_utf8_lossy(record.name().unwrap()).into_owned())
        .collect();
    assert_eq!(added, ["unplaced", "chr2", "pair", "pair", "late"]);
}

#[test]
fn test_built_records_pile_up_in_memory_with_the_real_engine() {
    let mut builder = SamBuilder::new().read_length(10);
    builder.add_pair(Pair::filled(100, 105, 'A', 10));
    builder.add_frag(Frag::at(103));
    let mut pileups = builder.to_pileup_builder();
    let pileup = pileups.pileup("chr1", 105).unwrap();
    assert_eq!(pileup.unfiltered_depth(), 3);
    let flags: Vec<u16> = pileup.iter().map(|entry| entry.flags().bits()).collect();
    assert_eq!(flags, [99, 0, 147]);
    let ends: Vec<Option<usize>> = pileup
        .iter()
        .map(|entry| entry.template_end_distance().unwrap())
        .collect();
    assert_eq!(ends, [Some(8), None, Some(6)]);
    pileups.close().unwrap();
}

#[test]
fn test_write_bam_writes_the_header_and_records_in_coordinate_order() {
    let directory = std::env::temp_dir().join(format!("streampile-testing-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("reads.bam");
    let mut builder = SamBuilder::new().read_length(10);
    builder.add_frag(Frag::at(30));
    builder.add_pair(Pair::at(10, 20));
    builder.write_bam(&path).unwrap();
    let mut reader = bam::io::reader::Builder.build_from_path(&path).unwrap();
    let header = reader.read_header().unwrap();
    assert_eq!(header.reference_sequences().len(), 25);
    let starts: Vec<Option<Position>> = reader
        .records()
        .map(|record| record.unwrap().alignment_start().transpose().unwrap())
        .collect();
    assert_eq!(starts, [10, 20, 30].map(Position::new));
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn test_the_record_helpers_change_one_field() {
    let mut builder = SamBuilder::new().read_length(4);
    let record = builder.add_pair(Pair::at(10, 20)).remove(0);
    let moved = SamBuilder::with_mate_reference_sequence_id(record.clone(), 3);
    assert_eq!(moved.mate_reference_sequence_id(), Some(3));
    let bare = SamBuilder::without_mate_cigar(record.clone());
    assert_eq!(tag(&bare, Tag::MATE_CIGAR), None);
    let duplicate = SamBuilder::with_flags(record.clone(), 0x400);
    assert!(duplicate.flags().is_duplicate() && duplicate.flags().is_first_segment());
    let bases = SamBuilder::with_bases(record, "GGGG");
    assert_eq!(bases.sequence().as_ref(), b"GGGG");
}

#[test]
fn test_random_bases_are_as_many_as_the_cigar_reads() {
    let mut builder = SamBuilder::new().read_length(10);
    let pair = builder.add_pair(Pair {
        cigar1: Some("5S10M2D3M".into()),
        cigar2: Some("4M1I4M".into()),
        ..Pair::at(10, 20)
    });
    let lengths: Vec<usize> = pair.iter().map(|record| record.sequence().len()).collect();
    assert_eq!(lengths, [18, 9]);
    assert_eq!(builder.add_frag(Frag::at(30))[0].sequence().len(), 10);
}

# Test targets

Every public function, class and method declared in `src/` and every
command-line option chromap registers (`src/chromap_driver.cc`), each
assigned to one of four areas so four writers can work in parallel on
disjoint files. A row is ticked `[x]` when tests check its behaviour against
simulation truth or a hand-worked value, `[~]` when they cover part of it,
`[ ]` when nothing does. A test that only checks that chromap exits 0 is not
coverage. "Kept out" marks behaviour the base gets wrong: no suite test
asserts it, and the removed test is kept with a note in the search's
`known-bugs/` directory.

Private methods are not targets themselves; they are listed under the public
entry point that reaches them ("via ..."), so the row's tests must exercise
them.

## Areas and file ownership

| Area | Scope | Owns (create only these) |
|---|---|---|
| **A** | index build and load, minimizers, k-mer and seed encoding, sequence and FASTA/FASTQ parsing, read input files | `tests/test_a_*.py`, `tests/cpp/a_*.cc`, `tests/helpers_a.py` |
| **B** | mapping core: candidates, draft mappings, banded alignment and ksw, single- and paired-end best mappings, split alignment, MAPQ, multimappers, adapter trimming, presets, the candidate cache, thread equivalence | `tests/test_b_*.py`, `tests/cpp/b_*.cc`, `tests/helpers_b.py` |
| **C** | output and post-processing: BED, TagAlign, BEDPE, SAM and PAF records and writers, sorting, deduplication (bulk and cell level), Tn5 shift, `--low-mem` temp files and merging, multi-mapping allocation, summary CSV | `tests/test_c_*.py`, `tests/cpp/c_*.cc`, `tests/helpers_c.py` |
| **D** | single-cell and Hi-C: barcodes, whitelist correction, `--barcode-translate`, `--read-format`, feature/barcode matrix classes, pairs format and flipping, `--chr-order`; plus the command-line front end (help, version, every error path) | `tests/test_d_*.py`, `tests/cpp/d_*.cc`, `tests/helpers_d.py` |

Shared and read-only for writers: `chromap_helpers.py`, `chromap_build.py`,
`conftest.py`, `pytest.ini`, `cpp/chromap_test.h`, `test_harness.py`,
`test_pilot_*.py`, `cpp/pilot.cc`, `README.md`. A harness bug or a missing
general helper is reported back, not patched in place; helpers for one area
go in that area's `helpers_<area>.py`. Each writer ticks only its own rows
here.

Run budget per area (measured at 949043c, see README *Runtime*): at most 20
chromap mapping processes (about 0.3 s single-end, 0.5 s paired-end each)
and at most 2 `--low-mem` runs with BED or TagAlign output, presets `atac`
and `chip` included (about 3 s and 12 GB each, marked `lowmem`). Index builds
of the synthetic genome and C++ driver calls cost milliseconds. Share runs:
module-scoped fixtures, `shared=True`, and read sets that carry many designed
cases at once.

Test names below are `<file>::<test>` with `h` test_harness.py, `pc`
test_pilot_cli.py, `pu` test_pilot_unit.py.

## Excluded

| File | Why |
|---|---|
| `cxxopts.hpp` | vendored option parser (MIT); chromap's use of it is covered through the options below |
| `khash.h` | vendored hash table (attractive chaos); chromap's additions `kh_load`/`kh_save` (`kh_write`) are covered through `Index::Save`/`Index::Load` (A) and `KHASH_*` instantiations through their users |
| `kseq.h` | vendored FASTA/FASTQ reader; its behaviour is covered through `SequenceBatch` (A) |
| `ksw.cc`, `ksw.h`: `ksw_qinit`, `ksw_u8`, `ksw_i16`, `ksw_align`, `ksw_align2`, `ksw_extend`, `ksw_extend2`, `ksw_global`, `ksw_global2`, `main` (under `#ifdef _KSW_MAIN`) | unmodified klib code that chromap never calls. The chromap-specific `ksw_semi_global3` and its wrappers `ksw_semi_global`, `ksw_semi_global2` (char sequences, chromap's base table) are in scope (B): `ksw_semi_global3` computes every SAM CIGAR |
| `utils.h`: `GetRealTime`, `GetCPUTime` | wall/CPU clocks used only for log timestamps |

`mmcache.hpp` is in scope (B): it is chromap's own code on the mapping path
(the candidate cache queried for every read), and the CLI can only reach its
filled state after a 500,000-read batch, so its `Query`/`Update` contract is
tested directly through a C++ driver. The cache never changes the first
batch's results, which is all a CLI test of the suite maps.

## Functions and classes

### Area A

| Done | Area | File | Name(s) | What to check | Covered by |
|---|---|---|---|---|---|
| [ ] | A | utils.h | `CharToUint8`, `Uint8ToChar`, `char_to_uint8_table_`, `uint8_to_char_table_` | ACGT and acgt map to 0-3, every other byte to 4; 0-3 map to ACGT, 4-7 to N | |
| [~] | A | utils.h, sequence_batch.h | `GenerateSeedFromSequence`, `SequenceBatch::GenerateSeedFromSequenceAt` | 2-bit layout, first base highest; N encodes as A; positions past the end pad with A; start offset; lengths 1 and 32 | pu::test_pilot_cpp[seed_encoding_hand_worked] |
| [ ] | A | utils.h | `Hash64` | equals `chromap_helpers.hash64` for random keys and every mask 2k bits, k = 1..28; bijective on small masks | |
| [~] | A | utils.h, hit_utils.h | `GenerateMinimizer`, `HitToSequenceIndex`, `HitToSequencePosition`, `HitToStrand`, `AreTwoHitsOnTheSameStrand` | hit layout `((index << 32 | position) << 1) | strand`; round trips at index and position extremes | pu::test_pilot_cpp[hit_encoding_round_trip] |
| [~] | A | strand.h, minimizer.h | `Strand`, `Minimizer` (both constructors, `GetHash`, `GetHit`, `GetSequenceIndex`, `GetSequencePosition`, `GetSequenceStrand`, `operator<`) | getters decode the hit; `operator<` orders by hash, then hit | pu::test_pilot_cpp[hit_encoding_round_trip] |
| [~] | A | minimizer_generator.h/.cc | `MinimizerGenerator::MinimizerGenerator`, `GenerateMinimizers` | the k-mers that are smallest in some window of w k-mers, hit encoding, position order, sequences shorter than one window (`chromap_helpers.reference_minimizers`); N bases reset the window (and the last window minimizer before an N run, which minimap2-style winnowing can drop); lower case; even k (k-mers equal to their reverse complement are skipped); sequences shorter than k; appending to a non-empty vector | pu::test_minimizer_generator_matches_definition (N-free, odd k) |
| [ ] | A | sequence_batch.h/.cc | `SequenceBatch` (both constructors), `InitializeLoading`, `FinalizeLoading`, `LoadAllSequences`, `LoadBatch`, `LoadOneSequenceAndSaveAt`, `GetNumSequences`, `GetMaxBatchSize`, `GetNumBases`, `GetSequenceBatch`, `GetSequenceAt`, `GetSequenceLengthAt`, `GetSequenceNameAt`, `GetSequenceNameLengthAt`, `GetSequenceQualAt`, `GetSequenceIdAt`, `SwapSequenceBatch` | FASTA (multi-line, lower case, N) and FASTQ, plain and gzipped; names cut at the first whitespace; empty records skipped; ids count every record ever loaded; `LoadBatch` stops at the batch size and returns 0 at the end; `LoadOneSequenceAndSaveAt` refuses to overwrite an earlier slot; missing file and truncated file exit non-zero with chromap's messages (`ExitWithMessage`) | |
| [ ] | A | sequence_batch.h | `GetNegativeSequenceBatch`, `GetNegativeSequenceAt`, `PrepareNegativeSequenceAt`, `TrimSequenceAt`, `CorrectBaseAt`, `ReorderSequences`, `GetSequenceNsAt`, `IsNInSequenceAt` | reverse complement (N stays N, lower case complemented to upper case); trimming shortens sequence, quality and the reverse complement consistently and is a no-op for a longer length; `ReorderSequences` places sequence i at rank i; N positions in both orders | |
| [ ] | A | index_parameters.h | `IndexParameters` | defaults k 17, w 7 | |
| [ ] | A | index_utils.h | `GenerateHashInLookupTable`, `GenerateEntryValueInLookupTable`, `GenerateOffsetInOccurrenceTable`, `GenerateNumOccurrenceInOccurrenceTable`, `SequenceIndexAndPositionToCandidatePosition`, `GenerateCandidatePositionFromOccurrenceTableEntry`, `IsSingletonLookupKey`, `RepetitiveSeedStats`, `CandidatePositionWithListIndex`, `HeapMergeCandidatePositionLists`, `KHashFunctionForIndex`, `KHashEqForIndex` | bit layouts by hand; heap merge of sorted lists equals the sorted concatenation (empty lists, duplicates, one list) | |
| [ ] | A | index.h/.cc | `Index::Index` (both), `Construct`, `Save`, `Load`, `Destroy`, `Statistics`, `CheckIndex`, `GetKmerSize`, `GetWindowSize`, `GetLookupTableSize` | the index file (int32 k, int32 w, uint32 table size, khash dump: n_buckets, size, n_occupied, upper_bound, flags, keys, values; uint32 occurrence count, uint64 occurrences) parsed in Python holds every reference minimizer: singletons with key bit 0 set and the hit as value, others with `offset << 32 | count` and their hits sorted in the occurrence table; Save, Load, Save is byte-identical; `GetLookupTableSize` is the number of distinct minimizer hashes; one test may pin the md5 of the default genome's index ("pins the current output"; the search keeps the index format fixed) | |
| [ ] | A | chromap.h/.cc | `Chromap::Chromap(const IndexParameters &)`, `Chromap::ConstructIndex` | index building end to end through `-i` (options below); FASTA with N runs, lower case, a contig shorter than k, many contigs | |
| [ ] | A | chromap_driver.cc | `GetMatchedFilePaths` (both overloads, file-local) | `-1`/`-2`/`-b` take comma-separated lists and glob patterns, expanded in order and printed to stderr; a pattern that matches nothing exits non-zero (`glob() failed`) | |
| [ ] | A | chromap.h/.cc | via `MapSingleEndReads`/`MapPairedEndReads`: `LoadSingleEndReadsWithBarcodes`, `LoadPairedEndReadsWithBarcodes` | several input files mapped as one run (read ids continue across files); unequal read counts in `-1`/`-2` (and `-b`) exit with "Numbers of reads and barcodes don't match!"; parallel parsing (paired-end at `-t` >= 3) gives the same output as `-t 1`; reads shorter than `--min-read-length` (either mate) are dropped | |

### Area B

| Done | Area | File | Name(s) | What to check | Covered by |
|---|---|---|---|---|---|
| [ ] | B | candidate.h | `Candidate` (`GetReferenceSequenceIndex`, `GetReferenceSequencePosition`, `GetCount`, `operator<`) | position layout; order by count descending, then position | |
| [ ] | B | candidate_position_generating_config.h | `CandidatePositionGeneratingConfig` | `IsFrequentSeed`/`IsRepetitiveSeed` thresholds are inclusive; getters | |
| [ ] | B | index.h/.cc | `Index::GenerateCandidatePositions`, `Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand` (via `GenerateCandidatePositionFromHits`, `UpdateRepetitiveSeedStats`) | candidate start = reference position minus read position (same strand) or reference + read position - k + 1 (opposite strand), sorted; seeds at or above the frequency threshold skipped and counted as repetitive with the repetitive length rule; heap-merge and sort paths agree; mate-window search returns only positions within the search range of the best mate candidates and returns a negative count when the mate has too many candidates | |
| [ ] | B | candidate_processor.h/.cc | `CandidateProcessor::CandidateProcessor`, `GenerateCandidates`, `SupplementCandidates`, `ReduceCandidatesForPairedEndRead` (via `GenerateCandidatesOnOneStrand`, `MergeCandidates`, `GenerateCandidatesFromRepetitiveReadWithMateInfoOnOneStrand`, `ReduceCandidatesForPairedEndReadOnOneDirection`) | `-s` minimum seeds; second round with `-f` second value when every seed is frequent; supplementing a repetitive mate from its partner; pair reduction keeps candidates within the insert distance | |
| [ ] | B | mapping_metadata.h, paired_end_mapping_metadata.h | `MappingMetadata` (every public method, `PrintCandidates`), `PairedEndMappingMetadata` (every public method) | counters, sorting by position, buffer swap, `GetBestMappings` picks the vector for each strand combination; protected members are reachable from a test subclass (`using`), no access hacks | |
| [ ] | B | draft_mapping.h, draft_mapping_generator.h/.cc | `DraftMapping`, `DraftMappingGenerator::DraftMappingGenerator`, `GenerateDraftMappings` (via `IsValidCandidate`, `GenerateNonSplitDraftMappingSupportedByAllMinimizers`, `GenerateDraftMappingsOnOneStrandUsingSIMD`, `GenerateDraftMappingsOnOneStrand`) | edit distance of each draft mapping equals a Python banded edit distance; the three code paths chosen by `-e` (below 8: 8 SSE lanes; 8 to 15: 4 lanes; 16 and above: scalar) give the same mappings for reads all of them can map; candidates near contig ends | |
| [ ] | B | alignment.h/.cc | `GetLongestMatchLength`, `AdjustGapBeginning`, `GenerateNMAndMDTag`, `BandedAlignPatternToText`, `BandedAlignPatternToTextWithDropOff`, `BandedAlignPatternToTextWithDropOffFrom3End`, `BandedAlign4PatternsToText`, `BandedAlign8PatternsToText`, `BandedTraceback`, `BandedTracebackToEnd` | against a Python banded edit-distance reference on small hand-made and random cases (substitutions, insertions, deletions, errors above the threshold); the 4- and 8-lane versions agree with the scalar one lane by lane; NM/MD by hand | |
| [ ] | B | ksw.h/.cc | `ksw_semi_global3`, `ksw_semi_global`, `ksw_semi_global2` | score, CIGAR and start/end on small cases against a Python semi-global affine-gap reference with chromap's matrix (`-A`, `-B`, `-O`, `-E`) | |
| [ ] | B | mapping_in_memory.h | `MappingInMemory`, `PairedEndMappingInMemory` | fragment start, length and per-strand alignment lengths for both read-1 strands | |
| [ ] | B | mapping_parameters.h | `MappingParameters`, `MappingOutputFormat`, `GetNumVPULanes` | every default listed in the options table; lanes 8, 4, 4, 0 at `-e` 7, 8, 15, 16 | |
| [~] | B | mapping_generator.h/.cc | `MappingGenerator::MappingGenerator`, `GenerateBestMappingsForSingleEndRead`, `GenerateBestMappingsForPairedEndRead` (via `ProcessBestMappingsForSingleEndRead`, `GenerateBestMappingsForPairedEndReadOnOneDirection`, `ProcessBestMappingsForPairedEndReadOnOneDirection`, `GetRefStartEndPositionForReadFromMapping`, `GetMAPQForSingleEndRead`, `GetMAPQForPairedEndRead`, `EmplaceBack*MappingRecord`) | true positions of single-end and paired-end reads with substitutions and indels; MAPQ worked from the formulas for unique, second-best and repetitive cases (`dupmm`, `rep8`, `dup2`); `-n` reports that many best mappings (reservoir sampling seeded with 11); `--drop-repetitive-reads`; insert limit `-l`; split alignment (`--split-alignment`) of chimeric pairs | pc::test_single_end_reads_map_to_their_origin, pc::test_paired_end_fragments |
| [ ] | B | mmcache.hpp | `mm_cache` (`mm_cache`, `SetKmerLength`, `Query`, `Update`, `DirectUpdateWeight`, `GetMemoryBytes`, `GetUpdateThreshold`, `PrintStats`), `_mm_cache_entry`, `_mm_history` (utils.h) | a slot stores an entry once it has had at least 10 updates and the entry's fingerprint makes up at least a fifth of them, and stops changing after 100; `Query` of the same minimizers (or of the reverse-complemented read) returns the slot index and the stored candidates shifted to the query's positions, and -1 otherwise; `GetUpdateThreshold` arithmetic for single- and paired-end blocks | |
| [~] | B | chromap.h/.cc | `Chromap::Chromap(const MappingParameters &)`, `MapSingleEndReads<T>`, `MapPairedEndReads<T>` (via `TrimAdapterForPairedEndRead`, `OutputMappingStatistics`) | mapping loop end to end; adapter trimming of pairs with inserts shorter than the read (seeds of half `--min-read-length`, at most one mismatch); stderr statistics ("Number of reads", "Number of mapped reads", ...) match the truth; outputs identical at `-t 1` and `-t 4` for every output type and mode area B touches | pc::test_thread_count_does_not_change_output |

### Area C

| Done | Area | File | Name(s) | What to check | Covered by |
|---|---|---|---|---|---|
| [ ] | C | mapping.h | `Mapping` (interface) | implemented by every record class below | |
| [~] | C | bed_mapping.h | `MappingWithBarcode`, `MappingWithoutBarcode`, `PairedEndMappingWithBarcode`, `PairedEndMappingWithoutBarcode`: constructors, `operator<`, `operator==`, `IsSamePosition`, `GetBarcode`, `Tn5Shift`, `IsPositiveStrand`, `GetStartPosition`, `GetEndPosition` | sort keys and duplicate keys by hand (single-end duplicates share a start whatever the strand; paired-end ones share start and length; barcoded ones also the barcode); Tn5 shift +4/-5 | pu::test_pilot_cpp[tn5_shift_single_end_record], pu::test_pilot_cpp[tn5_shift_paired_end_record] |
| [~] | C | sam_mapping.h | `SAMMapping` (constructor, `operator<`, `operator==`, `IsSamePosition`, `GetBarcode`, `Tn5Shift`, `IsPositiveStrand`, `GenerateCigarString`, `GenerateIntTagString`, `GetAlignmentLength`, `GetSequenceLength`, `GetStartPosition`, `GetEndPosition`, `GetByteSize`, `WriteToFile`, `LoadFromFile`), `bam_cigar_table`, `bam_cigar_*` macros, `BAM_F*` flags | CIGAR strings, query and reference lengths by hand; the constructor reverses QUAL for reverse-strand records and truncates SEQ to the CIGAR's query length; Write/Load round trip | pc::test_sam_fields_from_truth |
| [ ] | C | paf_mapping.h | `PAFMapping`, `PairedPAFMapping` (same methods as above, Write/Load round trip) | paired-end PAF coordinates kept out (see known bug paired-paf-coordinates) | |
| [ ] | C | temp_mapping.h | `TempMappingFileHandle<T>` (`GetCurrentMapping`, `HasMappings`, `InitializeTempMappingLoading`, `FinalizeTempMappingLoading`, `LoadTempMappingBlock` and its PAF/SAM/pairs specializations, `Next`) | blocks smaller than one reference's mappings, references with no mappings, the last reference | |
| [~] | C | mapping_writer.h/.cc | `MappingWriter<T>` for every record type except pairs: constructor, `OutputHeader`, `OutputMappings`, `OutputTempMappings`, `ProcessAndOutputMappingsInLowMemory`, `OutputSummaryMetadata`, `UpdateSummaryMetadata`, `UpdateSpeicalCategorySummaryMetadata`, `AdjustSummaryPairedEndOverCount` (via `AppendMapping` specializations, `OutputMappingsInVector`, `OutputTempMapping`, `FindBestMappingIndexFromDuplicates`) | every line format: single-end BED/TagAlign (7 columns), barcoded BED (fragments, 5 columns), barcoded TagAlign (6 columns), paired-end BED (fragment), paired-end TagAlign (two lines, count on the second for bulk), SAM (RNEXT `=`/name/`*`, PNEXT, TLEN sign, NM, MD, CB), single-end PAF; MAPQ threshold applied after deduplication; k-way merge of several temp files gives the same output as one; temp files are removed; `-o /dev/stdout` | pc::test_single_end_reads_map_to_their_origin, pc::test_paired_end_fragments, pc::test_sam_fields_from_truth |
| [~] | C | mapping_processor.h | `MappingProcessor<T>` (`SortOutputMappings`, `ParallelSortOutputMappings`, `RemovePCRDuplicate`, `AllocateMultiMappings`, `ApplyTn5ShiftOnMappings`, `MoveMappingsInBuffersToMappingContainer`, `OutputMappingStatistics`; via `BuildAugmentedTree`, `GetNumOverlappedMappings`), `ReadIdLess`, `StackCell` (utils.h) | sort order; duplicate groups and their counts (capped at 255); which duplicate is written where both modes agree (identical duplicates); multi-mapping allocation weights by overlapping unique mappings within the allocation distance, seeded generator | pc::test_tn5_shift_on_fragments |
| [ ] | C | summary_metadata.h | `SummaryMetadata` (`inverse_logit`, `OutputCounts`, `Output`, `UpdateCount`, `UpdateNonWhitelistCount`, `SetBarcodeLength`, `AdjustPairedEndOverCount`), `SummaryMetadataField`, `_barcodeSummaryMetadata` | columns total, duplicate, unmapped, lowmapq, cachehit, fric, estfrip (inverse logit of the five `--frip-est-params` coefficients), numcacheslots; bulk row barcode is empty; per-barcode rows plus `non-whitelist` with a whitelist; paired-end SAM counts halved; row order is not specified (compare as sets) | |
| [ ] | C | chromap.h | `K_MinHash` (`add`, `compute_cardinality`) | k smallest distinct values kept; cardinality `k * range / max - 1`, 0 below k values | |

### Area D

| Done | Area | File | Name(s) | What to check | Covered by |
|---|---|---|---|---|---|
| [ ] | D | sequence_effective_range.h | `SequenceEffectiveRange` (`InitializeParsing`, `FinalizeParsing`, `ParseFormatStringAndAppendEffectiveRange`, `Replace`) | `r1:start:end[:strand]`, `bc:...`, end -1, several ranges per read, `-` strand reverse-complements sequence and reverses quality; malformed strings rejected | |
| [ ] | D | barcode_translator.h | `BarcodeTranslator` (`SetTranslateTable`, `Translate`) | `to<TAB or comma>from` lines, plain and gzipped; barcodes longer than one entry translated chunk by chunk and joined with `-`; a barcode missing from the table exits non-zero; without a table the barcode sequence itself | |
| [ ] | D | utils.h | `BarcodeWithQual` (`operator>`) | ordering by score, then positions and bases | |
| [ ] | D | chromap.h/.cc | via `MapSingleEndReads`/`MapPairedEndReads`: `ParseReadFormat`, `SampleInputBarcodesAndExamineLength`, `LoadBarcodeWhitelist`, `ComputeBarcodeAbundance`, `CorrectBarcodeAt`, `OutputBarcodeStatistics` | whitelist exact matches, one-mismatch correction (and two with `--bc-error-threshold 2`), N in the barcode, correction probability from base quality and whitelist abundance against `--bc-probability-threshold`, reads with uncorrectable barcodes dropped unless `--output-mappings-not-in-whitelist`; "Number of barcodes in whitelist" and "Number of corrected barcodes" on stderr; length checks and the 5% whitelist check | |
| [ ] | D | chromap.h/.cc | `PairedEndReadWithBarcodeIsDuplicate`, `UpdateBarcodeAbundance` | private and never called at 949043c; no target unless a public path reaches them | |
| [ ] | D | chromap.h/.cc | via `MapPairedEndReads`/`MapSingleEndReads`: `GenerateCustomRidRanks`, `RerankCandidatesRid` | `--chr-order` (paired-end) orders output by the custom list, contigs not listed follow in FASTA order, an unknown name exits non-zero; single-end `--chr-order` kept out (see known bug single-end-chr-order) | |
| [ ] | D | pairs_mapping.h | `PairsMapping` (constructor, `operator<`, `operator==`, `IsSamePosition`, `GetBarcode`, `Tn5Shift`, `GetPosition`, `GetStrand`, `IsPositiveStrand`, `GetStartPosition`, `GetEndPosition`, `GetByteSize`, `WriteToFile`, `LoadFromFile`) | 1-based 5' positions; Write/Load round trip including the packed MAPQ/unique/dups field | |
| [ ] | D | mapping_writer.h/.cc, mapping_generator.cc | `MappingWriter<PairsMapping>` (`OutputHeader`, `AppendMapping`, `OutputTempMapping`), `MappingGenerator<PairsMapping>::EmplaceBackPairedEndMappingRecord` | header lines (`## pairs format v1.0.0`, `#shape: upper triangle`, `#chromsize` in pairs order, `#columns`); each pair flipped so chrom1 has the lower rank (or the lower position on one chromosome); positions are each read's 5' end (1-based); `--pairs-natural-chr-order` | |
| [ ] | D | feature_barcode_matrix.h/.cc, feature_barcode_matrix_writer.h | `FeatureBarcodeMatrix` (constructor, `OutputFeatureMatrix`), `FeatureBarcodeMatrixWriter` (`InitializeMatrixOutput`, `OutputPeaks` x2, `AppendBarcodeOutput`, `WriteMatrixOutputHead`, `AppendMatrixOutput`, `FinalizeMatrixOutput`, `SetBarcodeTranslateTable`, `SetBarcodeLength`), `Peak` | not reachable from the command line (the peak options are never registered and `-p` is commented out); C++ only: peaks are runs of coverage at or above the cutoff, bins of `bin_size`, matrix lines and barcodes with `-1` | |
| [ ] | D | chromap_driver.h/.cc | `ChromapDriver::ParseArgsAndRun`, `main`; file-local `AddIndexingOptions`, `AddMappingOptions`, `AddInputOptions`, `AddOutputOptions`, `AddDevelopmentOptions`, `AddPeakOptions` | option registration, precedence (preset first, then explicit options; output format flags SAM > pairs > PAF > TagAlign > BED regardless of order), help and version, every error path (options table); an option that is not registered is rejected with a non-zero exit naming it (cxxopts throws and chromap does not catch it, so the process aborts: assert the exit is non-zero and the message, not the signal) | |
| [ ] | D | utils.h | `ExitWithMessage` | message on stderr, non-zero exit (255) | |

## Options

Defaults are from `src/mapping_parameters.h`, `src/index_parameters.h` and
the driver. "Covered by" lists tests that set the option to a non-default
value, or check its default effect by name.

| Done | Area | Option | Default | What to check | Covered by |
|---|---|---|---|---|---|
| [ ] | A | `-i`, `--build-index` | | builds an index; needs `-r` and `-o` | |
| [ ] | A | `--min-frag-length` | 30 | <= 60: k 17, w 7; <= 80: k 19, w 10; else k 23, w 11; `-k`/`-w` given as well win | |
| [ ] | A | `-k`, `--kmer` | 17 | index built with k; mapping reads k from the index | |
| [ ] | A | `-w`, `--window` | 7 | as `-k` | |
| [ ] | A | `-r`, `--ref` | | reference FASTA for indexing and mapping (D: missing) | |
| [~] | A | `-1`, `--read1` | | lists, globs, FASTA/FASTQ, gzip | pc::* |
| [~] | A | `-2`, `--read2` | | as `-1`; count mismatch | pc::test_paired_end_fragments |
| [ ] | A | `-x`, `--index` | | index loaded; k and w from the file | |
| [ ] | A | `--min-read-length` | 30 | shorter reads (either mate) dropped; also the adapter-trimming minimum overlap (B) | |
| [ ] | B | `--preset` | | `atac`: insert 2000, trim adapters, dedup, Tn5 shift, BED, low-mem; `chip`: insert 2000, dedup, BED, low-mem; `hic`: `-e 4`, `-q 1`, split alignment, low-mem, pairs; options given after it override; unknown value (D) | |
| [ ] | B | `--split-alignment` | off | chimeric pairs reported at both parts | |
| [ ] | B | `-e`, `--error-threshold` | 8 | reads with more errors unmapped; the SIMD lane paths | |
| [ ] | B | `-s`, `--min-num-seeds` | 2 | | |
| [ ] | B | `-f`, `--max-seed-frequencies` | 500,1000 | | |
| [ ] | B | `-l`, `--max-insert-size` | 1000 | pairs farther apart are not reported as a pair | |
| [~] | B | `-q`, `--MAPQ-threshold` | 30 | mappings below are not written; MAPQ 0 multimappers written at `-q 0` | pc::test_thread_count_does_not_change_output |
| [ ] | B | `--trim-adapters` | off | Nextera read-through trimmed so short inserts map | |
| [~] | B | `-t`, `--num-threads` | 1 | outputs identical at 1 and 4 threads (paired-end: above 5000 pairs to get two mapping tasks) | pc::test_thread_count_does_not_change_output |
| [ ] | B | `-A`, `--match-score` | 1 | SAM CIGAR and ksw scores | |
| [ ] | B | `-B`, `--mismatch-penalty` | 4 | | |
| [ ] | B | `-O`, `--gap-open-penalties` | 6,6 | | |
| [ ] | B | `-E`, `--gap-extension-penalties` | 1,1 | | |
| [ ] | B | `-n`, `--max-num-best-mappings` | 1 | reports up to n locations of a multimapper; capped by `--drop-repetitive-reads` with a warning | |
| [ ] | B | `--drop-repetitive-reads` | 500000 | | |
| [ ] | B | `--cache-size` | 4000003 | accepted in [2000000, 15000000], outputs unchanged (D: out of range) | |
| [ ] | B | `--cache-update-param` | 0.01 | accepted in [0, 1] (D: out of range) | |
| [ ] | B | `--debug-cache` | off | debug lines on stdout; output file unchanged | |
| [~] | C | `-o`, `--output` | | output path; `/dev/stdout` (temp files then go to the working directory) | pc::* |
| [~] | C | `--BED` | default | formats in the MappingWriter row | pc::test_single_end_reads_map_to_their_origin |
| [ ] | C | `--TagAlign` | | | |
| [~] | C | `--SAM` | | | pc::test_sam_fields_from_truth |
| [ ] | C | `--PAF` | | single-end; paired-end kept out | |
| [ ] | C | `--remove-pcr-duplicates` | off | | |
| [ ] | C | `--remove-pcr-duplicates-at-bulk-level` | on for single-cell | duplicates across barcodes collapse to the most supported barcode | |
| [ ] | C | `--remove-pcr-duplicates-at-cell-level` | | duplicates collapse within a barcode only | |
| [~] | C | `--Tn5-shift` | off | | pc::test_tn5_shift_on_fragments |
| [ ] | C | `--low-mem` | off | same groups and counts as in memory; temp files merged and removed | |
| [ ] | C | `--summary` | | | |
| [ ] | C | `--frip-est-params` | -1.0996;4.2391;3.0164e-05;-2.1087e-04;-5.5825e-05 | estfrip column; parsed only for paired-end runs, where anything but 5 numbers exits non-zero | |
| [ ] | C | `--turn-off-num-uniq-cache-slots` | off | drops the numcacheslots column | |
| [ ] | C | `--k-for-minhash` | 250 | accepted in [1, 1999] (D: out of range) | |
| [ ] | C | `--allocate-multi-mappings` | off | | |
| [ ] | C | `--multi-mapping-allocation-distance` | 0 | | |
| [ ] | C | `--multi-mapping-allocation-seed` | 11 | | |
| [ ] | D | `-b`, `--barcode` | | single-cell mode; fragments output | |
| [ ] | D | `--barcode-whitelist` | | error without `-b`; warning on stderr when `-b` comes without it | |
| [ ] | D | `--read-format` | | reads and barcodes cut from the given ranges | |
| [ ] | D | `--bc-error-threshold` | 1 | | |
| [ ] | D | `--bc-probability-threshold` | 0.9 | | |
| [ ] | D | `--skip-barcode-check` | off | skips the 5% whitelist check | |
| [ ] | D | `--output-mappings-not-in-whitelist` | off | | |
| [ ] | D | `--barcode-translate` | | | |
| [ ] | D | `--chr-order` | | paired-end only (single-end kept out) | |
| [ ] | D | `--pairs` | | paired-end only: single-end exits "No support for single-end HiC yet!" | |
| [ ] | D | `--pairs-natural-chr-order` | | | |
| [ ] | D | `-h`, `--help` | | usage on stderr, exit 0; also printed without `-i` or `-1` | |
| [ ] | D | `-v`, `--version` | | `0.3.3-r521` on stderr, exit 0 | |
| [ ] | D | errors | | missing `-r`, `-o`, `-x` with `-i` or `-1`; unknown `--preset`; `--cache-size`, `--cache-update-param`, `--k-for-minhash` out of range; `--frip-est-params` malformed (paired-end); `--barcode-whitelist` without `-b`; barcode lengths unequal in the input or whitelist, whitelist length differing from the reads, barcodes over 32 bp; fewer than 5% of barcodes in the whitelist; unregistered options (`--cell-by-bin`, `--bin-size`, `--depth-cutoff`, `--peak-min-length`, `--peak-merge-max-length`, `-p`) and malformed values | |

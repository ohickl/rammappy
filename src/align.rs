use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use rammap::align::index::Index as RustIndex;
use rammap::align::map::AlignFlags;
use rammap::align::occurrence_sidecar::{
    OccurrenceRecord, OccurrenceSidecarMetadata, OccurrenceSidecarWriter,
};
use rammap::align::partition::{map_partitioned_fasta_to_paf, PartitionedMapConfig};
use rammap::api::{
    strainxpress_sr_ava_config, Aligner as RustAligner, CigarOp as RustCigarOp,
    Mapping as RustMapping, Preset as RustPreset, Strand as RustStrand,
};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::path::PathBuf;
use std::sync::Arc;

fn digest32(value: &Bound<'_, PyBytes>, name: &str) -> PyResult<[u8; 32]> {
    value.as_bytes().try_into().map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(format!("{name} must contain exactly 32 bytes"))
    })
}

fn thread_pool(threads: usize) -> PyResult<Arc<ThreadPool>> {
    if threads == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "threads must be positive",
        ));
    }
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map(Arc::new)
        .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
}

const PARTITION_CAPABILITY_DESCRIPTOR: &str =
    "sx-native-partition-v1|xc-independent-occurrence-v1|occurrence-fasta-v1|project-bridge-admitted-v1";

/// The mapping presets available in `rammappy`.
///
/// These presets configure the aligner for different sequencing technologies
/// and use cases, tuning heuristics and scoring.
///
/// Examples:
///     >>> from rammappy import Index, Aligner, Preset
///     >>> index = Index.build([(b"target1", b"ATGC...")])
///     >>> aligner = Aligner(index, preset=Preset.MapOnt)
#[gen_stub_pyclass_enum]
#[pyclass(module = "rammappy._rammappy", eq, eq_int, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum Preset {
    MapOnt,
    MapHifi,
    Sr,
    StrainxpressSrAva,
    Splice,
    Asm5,
    Asm10,
    Asm20,
    MapPb,
}

impl From<Preset> for RustPreset {
    fn from(preset: Preset) -> Self {
        match preset {
            Preset::MapOnt => RustPreset::MapOnt,
            Preset::MapHifi => RustPreset::MapHifi,
            Preset::Sr => RustPreset::Sr,
            Preset::StrainxpressSrAva => RustPreset::StrainxpressSrAva,
            Preset::Splice => RustPreset::Splice,
            Preset::Asm5 => RustPreset::Asm5,
            Preset::Asm10 => RustPreset::Asm10,
            Preset::Asm20 => RustPreset::Asm20,
            Preset::MapPb => RustPreset::MapPb,
        }
    }
}

/// Strand orientation of an alignment.
///
/// Represents whether the query mapped to the forward or reverse strand of the target.
///
/// Attributes:
///     Forward: The forward strand.
///     Reverse: The reverse complement strand.
#[gen_stub_pyclass_enum]
#[pyclass(module = "rammappy._rammappy", eq, eq_int, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl From<RustStrand> for Strand {
    fn from(s: RustStrand) -> Self {
        match s {
            RustStrand::Forward => Strand::Forward,
            RustStrand::Reverse => Strand::Reverse,
        }
    }
}

impl From<Strand> for RustStrand {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Forward => RustStrand::Forward,
            Strand::Reverse => RustStrand::Reverse,
        }
    }
}

/// BAM CIGAR operation encodings as a Python Enum.
///
/// The values correspond to the official BAM specification.
#[gen_stub_pyclass_enum]
#[pyclass(eq, eq_int, module = "rammappy._rammappy")]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CigarOp {
    M = 0,
    I = 1,
    D = 2,
    N = 3,
    S = 4,
    H = 5,
    P = 6,
    EQ = 7,
    X = 8,
    B = 9,
}

impl From<u8> for CigarOp {
    fn from(op: u8) -> Self {
        match op {
            0 => CigarOp::M,
            1 => CigarOp::I,
            2 => CigarOp::D,
            3 => CigarOp::N,
            4 => CigarOp::S,
            5 => CigarOp::H,
            6 => CigarOp::P,
            7 => CigarOp::EQ,
            8 => CigarOp::X,
            9 => CigarOp::B,
            _ => CigarOp::M,
        }
    }
}

/// Structured CIGAR operation element (length and operation type).
///
/// Attributes:
///     len (int): Operation length.
///     op (CigarOp): Operation type enum.
#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", get_all, from_py_object)]
#[derive(Clone)]
pub struct CigarElement {
    pub len: u32,
    pub op: CigarOp,
}

impl From<RustCigarOp> for CigarElement {
    fn from(op: RustCigarOp) -> Self {
        CigarElement {
            len: op.len,
            op: op.op.into(),
        }
    }
}

/// Python representation of an alignment `Mapping`.
///
/// Mappings are lazily evaluated: the actual Rust-level objects are preserved
/// until accessed via the Python properties.
///
/// Attributes:
///     target_name (bytes): The name of the target sequence.
///     target_id (int): Target sequence numeric index.
///     target_len (int): Target sequence length.
///     target_start (int): Start coordinate on the target.
///     target_end (int): End coordinate on the target.
///     query_start (int): Start coordinate on the query.
///     query_end (int): End coordinate on the query.
///     strand (Strand): Orientation of the alignment.
///     score (int): Alignment score.
///     mapq (int): Mapping quality (0-255).
///     is_primary (bool): True if this is the primary alignment.
///     is_supplementary (bool): True if this is a supplementary alignment.
///     is_spliced (bool): True if this alignment contains splice junctions.
///     trans_strand (Strand | None): Transcript strand for splice alignments, if known.
///     matches (int): Number of matching bases.
///     block_len (int): Alignment block length.
///     edit_distance (int): Edit distance (NM tag).
///     divergence (float): Sequence divergence (0.0 = identical).
///     cigar (bytes | None): The CIGAR string, if requested during alignment.
///     cigar_ops (list[CigarElement] | None): Structured CIGAR operations, if requested.
///     cs (bytes | None): CS tag string, if requested.
///     md (bytes | None): MD tag string, if requested.

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy")]
pub struct Mapping {
    // Hold the underlying Rust mapping object directly.
    inner: RustMapping,
}

#[gen_stub_pymethods]
#[pymethods]
impl Mapping {
    /// Return the target name as Python `bytes`.
    #[getter]
    fn target_name<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, self.inner.target_name.as_bytes())
    }

    #[getter]
    fn target_start(&self) -> usize {
        self.inner.target_start
    }

    #[getter]
    fn target_end(&self) -> usize {
        self.inner.target_end
    }

    #[getter]
    fn target_len(&self) -> usize {
        self.inner.target_len
    }

    #[getter]
    fn query_start(&self) -> usize {
        self.inner.query_start
    }

    #[getter]
    fn query_end(&self) -> usize {
        self.inner.query_end
    }

    #[getter]
    fn strand(&self) -> Strand {
        self.inner.strand.into()
    }

    #[getter]
    fn score(&self) -> i32 {
        self.inner.score
    }

    #[getter]
    fn mapq(&self) -> i32 {
        self.inner.mapq
    }

    #[getter]
    fn is_primary(&self) -> bool {
        self.inner.is_primary
    }

    #[getter]
    fn is_supplementary(&self) -> bool {
        self.inner.is_supplementary
    }

    #[getter]
    fn is_spliced(&self) -> bool {
        self.inner.is_spliced
    }

    #[getter]
    fn trans_strand(&self) -> Option<Strand> {
        self.inner.trans_strand.map(|s| s.into())
    }

    #[getter]
    fn matches(&self) -> usize {
        self.inner.matches
    }

    #[getter]
    fn block_len(&self) -> usize {
        self.inner.block_len
    }

    #[getter]
    fn edit_distance(&self) -> u32 {
        self.inner.edit_distance
    }

    #[getter]
    fn target_id(&self) -> usize {
        self.inner.target_id
    }

    #[getter]
    fn divergence(&self) -> f64 {
        self.inner.divergence
    }

    /// Returns the optional CIGAR string as a lazy byte array.
    #[getter]
    fn cigar<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.inner
            .cigar
            .as_ref()
            .map(|s| PyBytes::new(py, s.as_bytes()))
    }

    /// Returns the structured CIGAR operations.
    #[getter]
    fn cigar_ops(&self) -> Option<Vec<CigarElement>> {
        self.inner
            .cigar_ops
            .as_ref()
            .map(|ops| ops.iter().map(|&op| op.into()).collect())
    }

    /// Returns the optional cs string as a lazy byte array.
    #[getter]
    fn cs<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.inner
            .cs
            .as_ref()
            .map(|s| PyBytes::new(py, s.as_bytes()))
    }

    /// Returns the optional MD string as a lazy byte array.
    #[getter]
    fn md<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyBytes>> {
        self.inner
            .md
            .as_ref()
            .map(|s| PyBytes::new(py, s.as_bytes()))
    }

    fn __str__(&self) -> String {
        let strand = match self.inner.strand {
            RustStrand::Forward => "+",
            RustStrand::Reverse => "-",
        };
        let primary = if self.inner.is_primary { "*" } else { "" };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
            self.inner.query_start,
            self.inner.query_end,
            strand,
            self.inner.target_name,
            self.inner.target_len,
            self.inner.target_start,
            self.inner.target_end,
            self.inner.score,
            primary
        )
    }

    fn __repr__(&self) -> String {
        let strand = match self.inner.strand {
            RustStrand::Forward => "+",
            RustStrand::Reverse => "-",
        };
        format!(
            "<Mapping: target='{}' [{}:{}] query=[{}:{}] strand='{}' mapq={} score={}>",
            self.inner.target_name,
            self.inner.target_start,
            self.inner.target_end,
            self.inner.query_start,
            self.inner.query_end,
            strand,
            self.inner.mapq,
            self.inner.score
        )
    }
}

/// A lazy iterator that provides `Mapping` objects.
///
/// Instead of allocating a list, we hold an iterator of Rust mappings
/// and materialize Python wrapper objects only when requested via `next()`.

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy")]
pub struct MappingIterator {
    iter: std::vec::IntoIter<RustMapping>,
}

#[gen_stub_pymethods]
#[pymethods]
impl MappingIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Mapping> {
        slf.iter.next().map(|m| Mapping { inner: m })
    }
}

/// The `Index` object represents a genomic sequence index.
///
/// It holds the internal Rust Index for alignment. You can construct an index
/// from a collection of sequences, or load it from a previously saved file.
///
/// Examples:
///     >>> from rammappy import Index
///     >>> index = Index.build([(b"target1", b"ATGC...")])
///     >>> index.save("my_index.mmi")
///     >>> loaded_index = Index.load("my_index.mmi")

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy")]
#[derive(Clone)]
pub struct Index {
    inner: RustIndex,
}

unsafe impl Send for Index {}
unsafe impl Sync for Index {}

#[gen_stub_pymethods]
#[pymethods]
impl Index {
    /// Build an index from target sequences.
    ///
    /// Args:
    ///     seqs (list[tuple[bytes, bytes]]): A list of tuples containing `(name, sequence)` as bytes.
    ///     w (int): Window size. Defaults to 10.
    ///     k (int): K-mer size. Defaults to 15.
    ///     is_hpc (bool): Homopolymer compressed. Defaults to False.
    ///     max_occ (int): Maximum occurrences. Defaults to 50000.
    ///
    /// Returns:
    ///     Index: The built index.
    #[staticmethod]
    #[pyo3(signature = (seqs, w=10, k=15, is_hpc=false, max_occ=50000, threads=1))]
    fn build(
        py: Python<'_>,
        seqs: Vec<(Bound<'_, PyBytes>, Bound<'_, PyBytes>)>,
        w: usize,
        k: usize,
        is_hpc: bool,
        max_occ: usize,
        threads: usize,
    ) -> PyResult<Self> {
        let rust_seqs = seqs
            .into_iter()
            .map(|(name, seq)| {
                (
                    String::from_utf8_lossy(name.as_bytes()).to_string(),
                    seq.as_bytes().to_vec(),
                )
            })
            .collect();
        let pool = thread_pool(threads)?;
        let inner =
            py.detach(move || pool.install(|| RustIndex::build(rust_seqs, w, k, is_hpc, max_occ)));
        Ok(Index { inner })
    }

    /// Build an index by streaming a FASTA/FASTQ file without materializing
    /// all target records in Python. If `occurrence_counts` is supplied, the
    /// native builder writes a versioned pre-cap occurrence sidecar while the
    /// index is finalized.
    #[staticmethod]
    #[pyo3(signature = (path, w=10, k=15, is_hpc=false, max_occ=50000, threads=1, occurrence_counts=None))]
    fn build_fasta(
        py: Python<'_>,
        path: PathBuf,
        w: usize,
        k: usize,
        is_hpc: bool,
        max_occ: usize,
        threads: usize,
        occurrence_counts: Option<PathBuf>,
    ) -> PyResult<Self> {
        let path = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path string"))?
            .to_string();
        let pool = thread_pool(threads)?;
        let inner = py.detach(move || {
            pool.install(|| {
                let index = if let Some(sidecar_path) = occurrence_counts.as_ref() {
                    let file = std::fs::File::create(sidecar_path)
                        .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?;
                    let mut sidecar = OccurrenceSidecarWriter::new(
                        file,
                        OccurrenceSidecarMetadata {
                            bucket_bits: 10u32.min((2 * k) as u32),
                            shard_id: 0,
                            shard_count: 1,
                            parameter_digest: [0; 32],
                            target_digest: [0; 32],
                        },
                    )
                    .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?;
                    let mut sidecar_error = None;
                    let index = RustIndex::build_fasta_with_occurrence_counts(
                        &path,
                        w,
                        k,
                        is_hpc,
                        max_occ,
                        |bucket, hash, count| {
                            if sidecar_error.is_none() {
                                if let Err(error) = sidecar.write_record(OccurrenceRecord {
                                    bucket,
                                    hash,
                                    count,
                                }) {
                                    sidecar_error = Some(error);
                                }
                            }
                        },
                    )
                    .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?;
                    if let Some(error) = sidecar_error {
                        return Err(pyo3::exceptions::PyIOError::new_err(error.to_string()));
                    }
                    sidecar
                        .finish()
                        .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?
                        .sync_all()
                        .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?;
                    index
                } else {
                    RustIndex::build_fasta(&path, w, k, is_hpc, max_occ)
                        .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))?
                };
                Ok(Index { inner: index })
            })
        });
        inner
    }

    /// Load an index from file.
    ///
    /// Args:
    ///     path (os.PathLike): The file path to load the index from.
    ///
    /// Returns:
    ///     Index: The loaded index.
    #[staticmethod]
    fn load(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path string"))?
            .to_string();
        py.detach(move || RustIndex::load(&path_str))
            .map(|idx| Index { inner: idx })
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Save the index to a file.
    ///
    /// Args:
    ///     path (os.PathLike): The file path to save the index to.
    fn save(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path string"))?
            .to_string();
        py.detach(move || self.inner.save(&path_str))
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Strip sequences from the index to save memory.
    ///
    /// This removes the actual sequence bytes from memory, which is useful when
    /// you only need to perform mapping and do not need base-level alignment (CIGAR).
    ///
    /// Examples:
    ///     >>> index.strip_sequences()
    fn strip_sequences(&mut self) {
        self.inner.strip_sequences();
    }

    #[getter]
    fn kmer_size(&self) -> usize {
        self.inner.kmer_size
    }

    #[getter]
    fn window_size(&self) -> usize {
        self.inner.window_size
    }

    #[getter]
    fn homopolymer_compressed(&self) -> bool {
        self.inner.homopolymer_compressed
    }

    /// Returns the sequence names in the index.
    ///
    /// Returns:
    ///     list[str]: The list of sequence names.

    #[getter]
    fn seq_names(&self) -> Vec<String> {
        self.inner.seqs.iter().map(|s| s.name.clone()).collect()
    }

    /// Returns the sequence for a given target name.
    ///
    /// Args:
    ///     name (str): The name of the target sequence.
    ///     start (int, optional): The 0-based start coordinate. Defaults to 0.
    ///     end (int, optional): The 0-based end coordinate. Defaults to the end of the sequence.
    ///
    /// Returns:
    ///     str: The requested sequence.
    ///
    /// Raises:
    ///     ValueError: If the sequence name is not found in the index.
    #[pyo3(signature = (name, start=None, end=None))]
    fn seq(&self, name: &str, start: Option<usize>, end: Option<usize>) -> PyResult<String> {
        let rid = self
            .inner
            .seqs
            .iter()
            .position(|s| s.name == name)
            .ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Sequence not found: {}", name))
            })?;

        let seq_len = self.inner.seqs[rid].len;
        let s = start.unwrap_or(0).min(seq_len);
        let e = end.unwrap_or(seq_len).min(seq_len);

        if s >= e {
            return Ok(String::new());
        }

        let nt4 = self.inner.get_region_nt4(rid, s, e);
        let ascii: String = nt4
            .into_iter()
            .map(|b| RustIndex::NT4_TO_ASCII[b as usize] as char)
            .collect();

        Ok(ascii)
    }
}

/// The `Aligner` orchestrates the alignment process.
///
/// It encapsulates the alignment configuration and the index to map query sequences against reference targets.
///
/// Examples:
///     >>> from rammappy import Index, Aligner, Preset
///     >>> index = Index.build([(b"target1", b"ATGC...")])
///     >>> aligner = Aligner(index, preset=Preset.MapOnt)
///     >>> for mapping in aligner.map(b"query1", b"ATGC..."):
///     ...     print(mapping.score)

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy")]
pub struct Aligner {
    inner: RustAligner,
    pool: Arc<ThreadPool>,
}

unsafe impl Send for Aligner {}
unsafe impl Sync for Aligner {}

#[gen_stub_pymethods]
#[pymethods]
impl Aligner {
    /// Describe the native partition contract and its source-bound project
    /// adapter admission level.
    #[staticmethod]
    fn partitioned_capability_descriptor() -> &'static str {
        PARTITION_CAPABILITY_DESCRIPTOR
    }

    /// Build an aligner with the StrainXpress short-read all-vs-all settings.
    ///
    /// This mirrors the parent project's Minimap2 contract: short-read mode,
    /// all chains, no exact diagonal, and the explicit scoring/filter values
    /// used by its oracle command. The index must already use k=21 and w=11.
    #[staticmethod]
    #[pyo3(signature = (index, threads=1))]
    fn from_strainxpress_sr_ava(index: &Index, threads: usize) -> PyResult<Self> {
        if threads == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "threads must be positive",
            ));
        }
        if index.inner.kmer_size != 21
            || index.inner.window_size != 11
            || index.inner.homopolymer_compressed
        {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "StrainXpress short-read all-vs-all requires k=21, w=11, and hpc=false",
            ));
        }
        let aligner = Self::new(
            index,
            Some(Preset::StrainxpressSrAva),
            true,
            false,
            false,
            threads,
        )?;
        let expected = AlignFlags::SHORT_READ
            | AlignFlags::NO_DIAG
            | AlignFlags::ALL_CHAINS
            | AlignFlags::NO_DUAL
            | AlignFlags::NO_LJOIN
            | AlignFlags::OUT_CIGAR;
        if aligner.inner.options().flags != expected {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "rammap core returned an unexpected StrainXpress flag set",
            ));
        }
        Ok(aligner)
    }

    /// Create a new aligner instance using an already built index.
    ///
    /// Args:
    ///     index (Index): The built index object.
    ///     preset (Preset): The preset configuration (e.g. `Preset.MapOnt`). Defaults to `Preset.MapOnt`.
    ///     do_cigar (bool): Whether to compute CIGAR strings. Defaults to `True`.
    ///     do_cs (bool): Whether to compute `cs` tags. Defaults to `True`.
    ///     do_md (bool): Whether to compute `md` tags. Defaults to `True`.
    ///
    /// Returns:
    ///     Aligner: The initialized aligner object.
    #[new]
    #[pyo3(signature = (index, preset=Some(Preset::MapOnt), do_cigar=true, do_cs=true, do_md=true, threads=1))]
    fn new(
        index: &Index,
        preset: Option<Preset>,
        do_cigar: bool,
        do_cs: bool,
        do_md: bool,
        threads: usize,
    ) -> PyResult<Self> {
        let preset_enum: RustPreset = preset.unwrap_or(Preset::MapOnt).into();
        let pool = thread_pool(threads)?;

        let mut buf = Vec::new();
        index
            .inner
            .save_part(&mut buf)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        let mut cursor = std::io::Cursor::new(buf);
        let mut inner = RustAligner::from_index_reader(&mut cursor, preset_enum)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

        {
            let cfg = inner.output_config_mut();
            cfg.do_cigar = do_cigar;
            cfg.do_cs = do_cs;
            cfg.do_md = do_md;
        }

        Ok(Aligner { inner, pool })
    }

    /// Create an aligner from a FASTA file.
    ///
    /// Args:
    ///     path (os.PathLike): The file path to the FASTA file.
    ///     preset (Preset): The preset configuration (e.g. `Preset.MapOnt`). Defaults to `Preset.MapOnt`.
    ///
    /// Returns:
    ///     Aligner: The initialized aligner object.
    #[staticmethod]
    #[pyo3(signature = (path, preset=Some(Preset::MapOnt), threads=1))]
    fn from_fasta(
        py: Python<'_>,
        path: PathBuf,
        preset: Option<Preset>,
        threads: usize,
    ) -> PyResult<Self> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path string"))?
            .to_string();
        let preset_enum: RustPreset = preset.unwrap_or(Preset::MapOnt).into();
        let pool = thread_pool(threads)?;
        let build_pool = pool.clone();
        py.detach(move || build_pool.install(|| RustAligner::from_fasta(&path_str, preset_enum)))
            .map(|inner| Aligner { inner, pool })
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Create an aligner from an index file.
    ///
    /// Args:
    ///     path (os.PathLike): The file path to the saved index file.
    ///     preset (Preset): The preset configuration (e.g. `Preset.MapOnt`). Defaults to `Preset.MapOnt`.
    ///
    /// Returns:
    ///     Aligner: The initialized aligner object.
    #[staticmethod]
    #[pyo3(signature = (path, preset=Some(Preset::MapOnt), threads=1))]
    fn from_index(
        py: Python<'_>,
        path: PathBuf,
        preset: Option<Preset>,
        threads: usize,
    ) -> PyResult<Self> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path string"))?
            .to_string();
        let preset_enum: RustPreset = preset.unwrap_or(Preset::MapOnt).into();
        let pool = thread_pool(threads)?;
        let build_pool = pool.clone();
        py.detach(move || build_pool.install(|| RustAligner::from_index(&path_str, preset_enum)))
            .map(|inner| Aligner { inner, pool })
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Map a query FASTA/FASTQ against target FASTA shards through the native
    /// raw-candidate spool and one global finalization pass.
    ///
    /// The three 32-byte digests are caller-authenticated identities for the
    /// immutable parameters, target catalog, and query stream. The method
    /// returns `(shard_count, query_count, mid_occ, output_bytes)`.
    #[pyo3(signature = (target_paths, query_path, output_path, spool_dir, parameter_digest, target_digest, query_digest, index_max_occ=50000, mid_occ_frac=0.0002, resume=false))]
    fn map_partitioned_fasta_to_paf(
        &self,
        py: Python<'_>,
        target_paths: Vec<PathBuf>,
        query_path: PathBuf,
        output_path: PathBuf,
        spool_dir: PathBuf,
        parameter_digest: &Bound<'_, PyBytes>,
        target_digest: &Bound<'_, PyBytes>,
        query_digest: &Bound<'_, PyBytes>,
        index_max_occ: usize,
        mid_occ_frac: f32,
        resume: bool,
    ) -> PyResult<(u32, u64, usize, u64)> {
        let parameter_digest = digest32(parameter_digest, "parameter_digest")?;
        let target_digest = digest32(target_digest, "target_digest")?;
        let query_digest = digest32(query_digest, "query_digest")?;
        let config = PartitionedMapConfig {
            target_paths,
            query_path,
            output_path,
            spool_dir,
            k: self.inner.index().kmer_size,
            w: self.inner.index().window_size,
            is_hpc: self.inner.index().homopolymer_compressed,
            index_max_occ,
            mid_occ_frac,
            options: self.inner.options().clone(),
            output: self.inner.output_config().clone(),
            parameter_digest,
            target_digest,
            query_digest,
            resume,
        };
        py.detach(move || {
            map_partitioned_fasta_to_paf(&config)
                .map(|receipt| {
                    (
                        receipt.shard_count,
                        receipt.query_count,
                        receipt.mid_occ,
                        receipt.output_bytes,
                    )
                })
                .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))
        })
    }

    /// Map independent occurrence queries against target FASTA shards through
    /// the native resumable transaction without constructing an `Aligner`.
    ///
    /// This is deliberately a fixed StrainXpress short-read all-vs-all
    /// operation. The caller supplies authenticated identities for the exact
    /// parameters, target projection, and query stream. The returned tuple is
    /// `(shard_count, query_count, mid_occ, output_bytes)`; the durable native
    /// manifest remains the source of the complete transaction receipt.
    #[staticmethod]
    #[pyo3(signature = (target_paths, query_path, output_path, spool_dir, parameter_digest, target_digest, query_digest, threads=1, index_max_occ=50000, mid_occ_frac=0.0002, resume=false))]
    fn map_partitioned_fasta_to_paf_resumable(
        py: Python<'_>,
        target_paths: Vec<PathBuf>,
        query_path: PathBuf,
        output_path: PathBuf,
        spool_dir: PathBuf,
        parameter_digest: &Bound<'_, PyBytes>,
        target_digest: &Bound<'_, PyBytes>,
        query_digest: &Bound<'_, PyBytes>,
        threads: usize,
        index_max_occ: usize,
        mid_occ_frac: f32,
        resume: bool,
    ) -> PyResult<(u32, u64, usize, u64)> {
        let parameter_digest = digest32(parameter_digest, "parameter_digest")?;
        let target_digest = digest32(target_digest, "target_digest")?;
        let query_digest = digest32(query_digest, "query_digest")?;
        let pool = thread_pool(threads)?;
        // The partition transaction recalibrates mid_occ from authenticated
        // occurrence sidecars.  The preset constructor nevertheless requires
        // a positive initial value, so use its locked short-read default as a
        // bootstrap value before the transaction overwrites it.
        let (options, output) = strainxpress_sr_ava_config(21, 11, false, 1000)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        let config = PartitionedMapConfig {
            target_paths,
            query_path,
            output_path,
            spool_dir,
            k: 21,
            w: 11,
            is_hpc: false,
            index_max_occ,
            mid_occ_frac,
            options,
            output,
            parameter_digest,
            target_digest,
            query_digest,
            resume,
        };
        py.detach(move || {
            pool.install(|| map_partitioned_fasta_to_paf(&config))
                .map(|receipt| {
                    (
                        receipt.shard_count,
                        receipt.query_count,
                        receipt.mid_occ,
                        receipt.output_bytes,
                    )
                })
                .map_err(|error| pyo3::exceptions::PyIOError::new_err(error.to_string()))
        })
    }

    /// Maps a single query sequence sequentially to the targets.
    ///
    /// Args:
    ///     query_name (bytes): The name of the query sequence.
    ///     query_seq (bytes): The query sequence.
    ///
    /// Returns:
    ///     MappingIterator: An iterator over the generated mappings.
    fn map(
        &self,
        py: Python<'_>,
        query_name: &Bound<'_, PyBytes>,
        query_seq: &Bound<'_, PyBytes>,
    ) -> MappingIterator {
        let name_bytes = query_name.as_bytes();
        let seq_bytes = query_seq.as_bytes();

        let map_result = py.detach(move || {
            let query_name_str = String::from_utf8_lossy(name_bytes);
            self.inner.map_seq(query_name_str.as_ref(), seq_bytes)
        });

        MappingIterator {
            iter: map_result.mappings.into_iter(),
        }
    }

    /// Performs highly parallelized batch alignments mapping over multiple queries.
    ///
    /// Bypasses the GIL to utilize multiple threads for parallelism (via Rayon).
    ///
    /// Args:
    ///     queries (list[tuple[bytes, bytes]]): A list of tuples containing `(name, sequence)` as bytes.
    ///
    /// Returns:
    ///     list[MappingIterator]: A list of iterators, one for each query sequence.
    #[pyo3(signature = (queries))]
    fn map_batch(
        &self,
        py: Python<'_>,
        queries: Vec<(Bound<'_, PyBytes>, Bound<'_, PyBytes>)>,
    ) -> PyResult<Vec<MappingIterator>> {
        // Keep owning Python references alive until every worker has finished.
        // The detached workers only receive raw pointers because Python objects
        // cannot be accessed without the GIL.
        struct RawQuery {
            name_ptr: *const u8,
            name_len: usize,
            seq_ptr: *const u8,
            seq_len: usize,
        }

        // Safety: RawQuery is manually marked as Send and Sync so it can cross thread boundaries.
        // We only use the pointers while the GIL is temporarily released, meaning the Python
        unsafe impl Send for RawQuery {}
        unsafe impl Sync for RawQuery {}

        let mut owners = Vec::with_capacity(queries.len());
        let mut raw_queries = Vec::with_capacity(queries.len());
        for (name, seq) in queries {
            let name = name.unbind();
            let seq = seq.unbind();
            let name_bytes = name.bind(py);
            let seq_bytes = seq.bind(py);
            raw_queries.push(RawQuery {
                name_ptr: name_bytes.as_bytes().as_ptr(),
                name_len: name_bytes.as_bytes().len(),
                seq_ptr: seq_bytes.as_bytes().as_ptr(),
                seq_len: seq_bytes.as_bytes().len(),
            });
            owners.push((name, seq));
        }

        // Release the GIL via `detach` to allow other Python threads to execute concurrently.
        let iterators = py.detach(move || {
            let _keep_alive = &owners;
            self.pool.install(|| {
                let _ = _keep_alive;
                raw_queries
                    .par_iter()
                    .map(|raw_q| {
                        // Safety: Reconstructing the slice is safe because we know the pointer is valid
                        let name_bytes =
                            unsafe { std::slice::from_raw_parts(raw_q.name_ptr, raw_q.name_len) };
                        let seq_bytes =
                            unsafe { std::slice::from_raw_parts(raw_q.seq_ptr, raw_q.seq_len) };

                        let query_name_str = String::from_utf8_lossy(name_bytes);
                        let map_result = self.inner.map_seq(&query_name_str, seq_bytes);

                        MappingIterator {
                            iter: map_result.mappings.into_iter(),
                        }
                    })
                    .collect()
            })
        });
        Ok(iterators)
    }

    /// Map queries into a compact little-endian numeric buffer.
    ///
    /// The returned pair is `(records, offsets)`. `records` is a fixed-width
    /// buffer with 93-byte records; `offsets[i]..offsets[i + 1]` identifies the
    /// byte range for query `i`. The record fields are query ID, target ID,
    /// target length, query start/end, target start/end, matches, block
    /// length, edit distance, score, map quality, and strand. Python callers
    /// can wrap `records` with `numpy.frombuffer` without per-mapping Python
    /// objects or a copy.
    #[pyo3(signature = (queries))]
    fn map_batch_packed<'py>(
        &self,
        py: Python<'py>,
        queries: Vec<(Bound<'_, PyBytes>, Bound<'_, PyBytes>)>,
    ) -> PyResult<(Bound<'py, PyBytes>, Vec<u64>)> {
        struct RawQuery {
            name_ptr: *const u8,
            name_len: usize,
            seq_ptr: *const u8,
            seq_len: usize,
        }
        unsafe impl Send for RawQuery {}
        unsafe impl Sync for RawQuery {}

        let mut owners = Vec::with_capacity(queries.len());
        let mut raw_queries = Vec::with_capacity(queries.len());
        for (name, seq) in queries {
            let name = name.unbind();
            let seq = seq.unbind();
            let name_bytes = name.bind(py);
            let seq_bytes = seq.bind(py);
            raw_queries.push(RawQuery {
                name_ptr: name_bytes.as_bytes().as_ptr(),
                name_len: name_bytes.as_bytes().len(),
                seq_ptr: seq_bytes.as_bytes().as_ptr(),
                seq_len: seq_bytes.as_bytes().len(),
            });
            owners.push((name, seq));
        }

        let packed = py.detach(move || {
            let _keep_alive = &owners;
            self.pool.install(|| {
                let _ = _keep_alive;
                let per_query: Vec<Vec<RustMapping>> = raw_queries
                    .par_iter()
                    .map(|raw_q| {
                        let name_bytes =
                            unsafe { std::slice::from_raw_parts(raw_q.name_ptr, raw_q.name_len) };
                        let seq_bytes =
                            unsafe { std::slice::from_raw_parts(raw_q.seq_ptr, raw_q.seq_len) };
                        let query_name_str = String::from_utf8_lossy(name_bytes);
                        self.inner.map_seq(&query_name_str, seq_bytes).mappings
                    })
                    .collect();

                let mut records = Vec::new();
                let mut offsets = Vec::with_capacity(per_query.len() + 1);
                offsets.push(0);
                for (query_id, mappings) in per_query.into_iter().enumerate() {
                    for mapping in mappings {
                        records.extend_from_slice(&(query_id as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.target_id as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.target_len as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.query_start as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.query_end as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.target_start as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.target_end as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.matches as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.block_len as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.edit_distance as u64).to_le_bytes());
                        records.extend_from_slice(&(mapping.score as i64).to_le_bytes());
                        records.extend_from_slice(&(mapping.mapq as i32).to_le_bytes());
                        records.push(match mapping.strand {
                            RustStrand::Forward => 0,
                            RustStrand::Reverse => 1,
                        });
                    }
                    offsets.push(records.len() as u64);
                }
                (records, offsets)
            })
        });

        Ok((PyBytes::new(py, &packed.0), packed.1))
    }

    /// Load splice junctions from a BED file.
    ///
    /// Args:
    ///     path (os.PathLike | str): Path to the BED file.
    #[pyo3(signature = (path))]
    fn load_junctions_bed(&mut self, path: PathBuf) -> PyResult<()> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path for BED file"))?;
        self.inner.load_junctions_bed(path_str).map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("Failed to load junctions: {}", e))
        })
    }

    /// Load splice junctions from a SPSC file.
    ///
    /// Args:
    ///     path (os.PathLike | str): Path to the SPSC file.
    ///     scale (float | None): Optional scaling factor.
    #[pyo3(signature = (path, scale=None))]
    fn load_junctions_spsc(&mut self, path: PathBuf, scale: Option<f32>) -> PyResult<()> {
        let path_str = path
            .to_str()
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("Invalid path for SPSC file"))?;
        self.inner
            .load_junctions_spsc(path_str, scale)
            .map_err(|e| {
                pyo3::exceptions::PyIOError::new_err(format!("Failed to load junctions: {}", e))
            })
    }

    /// Get the current mapping options.
    ///
    /// Returns:
    ///     MapOptions: A copy of the current mapping options.
    #[getter]
    fn options(&self, py: Python<'_>) -> PyResult<crate::options::PyMapOptions> {
        crate::options::PyMapOptions::from_map(py, self.inner.options().clone())
    }

    /// Set the mapping options.
    ///
    /// Args:
    ///     opts (MapOptions): The new mapping options to apply.
    #[setter]
    fn set_options(&mut self, opts: &Bound<'_, crate::options::PyMapOptions>) -> PyResult<()> {
        let options = opts.borrow().into_map(opts.py());
        self.inner.output_config_mut().do_cigar = options.flags.contains(AlignFlags::OUT_CIGAR);
        *self.inner.options_mut() = options;
        Ok(())
    }

    #[getter]
    fn seq_names(&self) -> Vec<String> {
        self.inner
            .index()
            .seqs
            .iter()
            .map(|s| s.name.clone())
            .collect()
    }

    /// Returns the sequence for a given target name.
    ///
    /// Args:
    ///     name (str): The name of the target sequence.
    ///     start (int, optional): The 0-based start coordinate. Defaults to 0.
    ///     end (int, optional): The 0-based end coordinate. Defaults to the end of the sequence.
    ///
    /// Returns:
    ///     str: The requested sequence.
    ///
    /// Raises:
    ///     ValueError: If the sequence name is not found in the index.
    #[pyo3(signature = (name, start=None, end=None))]
    fn seq(&self, name: &str, start: Option<usize>, end: Option<usize>) -> PyResult<String> {
        let index = self.inner.index();
        let rid = index
            .seqs
            .iter()
            .position(|s| s.name == name)
            .ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!("Sequence not found: {}", name))
            })?;

        let seq_len = index.seqs[rid].len;
        let s = start.unwrap_or(0).min(seq_len);
        let e = end.unwrap_or(seq_len).min(seq_len);

        if s >= e {
            return Ok(String::new());
        }

        let nt4 = index.get_region_nt4(rid, s, e);
        let ascii: String = nt4
            .into_iter()
            .map(|b| RustIndex::NT4_TO_ASCII[b as usize] as char)
            .collect();

        Ok(ascii)
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    if let Err(e) = m.add_class::<Preset>() {
        println!("Error adding Preset: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<Strand>() {
        println!("Error adding Strand: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<CigarOp>() {
        println!("Error adding CigarOp: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<CigarElement>() {
        println!("Error adding CigarElement: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<Index>() {
        println!("Error adding Index: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<Aligner>() {
        println!("Error adding Aligner: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<Mapping>() {
        println!("Error adding Mapping: {:?}", e);
        return Err(e);
    }
    if let Err(e) = m.add_class::<MappingIterator>() {
        println!("Error adding MappingIterator: {:?}", e);
        return Err(e);
    }

    Ok(())
}

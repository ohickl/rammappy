use pyo3::prelude::*;
use pyo3_stub_gen::derive::gen_stub_pyclass;
use rammap::align::map::{
    AlignFlags, AlignmentParams, ChainingParams, FilteringParams, MapOptions, PairedEndParams,
    ScoringParams, SeedingParams,
};

/// Stable Python constants for the behavioral bits carried by rammap.
#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "AlignFlags")]
pub struct PyAlignFlags;

#[pymethods]
impl PyAlignFlags {
    #[classattr]
    const NO_DIAG: u64 = 0x001;
    #[classattr]
    const NO_DUAL: u64 = 0x002;
    #[classattr]
    const NO_QUAL: u64 = 0x010;
    #[classattr]
    const OUT_CIGAR: u64 = 0x020;
    #[classattr]
    const SPLICE: u64 = 0x080;
    #[classattr]
    const SPLICE_FOR: u64 = 0x100;
    #[classattr]
    const SPLICE_REV: u64 = 0x200;
    #[classattr]
    const NO_LJOIN: u64 = 0x400;
    #[classattr]
    const INDEPEND_SEG: u64 = 0x800;
    #[classattr]
    const SHORT_READ: u64 = 0x1000;
    #[classattr]
    const FRAG_MODE: u64 = 0x2000;
    #[classattr]
    const NO_PRINT_2ND: u64 = 0x4000;
    #[classattr]
    const EQX: u64 = 0x8000;
    #[classattr]
    const LONG_CIGAR: u64 = 0x10000;
    #[classattr]
    const SOFTCLIP: u64 = 0x20000;
    #[classattr]
    const SPLICE_FLANK: u64 = 0x40000;
    #[classattr]
    const FOR_ONLY: u64 = 0x100000;
    #[classattr]
    const REV_ONLY: u64 = 0x200000;
    #[classattr]
    const HEAP_SORT: u64 = 0x400000;
    #[classattr]
    const ALL_CHAINS: u64 = 0x800000;
    #[classattr]
    const NO_END_FLT: u64 = 0x1000000;
    #[classattr]
    const HARD_MASK_LEVEL: u64 = 0x2000000;
    #[classattr]
    const SAM_HIT_ONLY: u64 = 0x4000000;
    #[classattr]
    const PAF_NO_HIT: u64 = 0x8000000;
    #[classattr]
    const NO_HASH_NAME: u64 = 0x40000000;
    #[classattr]
    const RMQ_CHAIN: u64 = 0x80000000;
    #[classattr]
    const QSTRAND: u64 = 0x100000000;
    #[classattr]
    const NO_INV: u64 = 0x200000000;
    #[classattr]
    const SPLICE_OLD: u64 = 0x400000000;
    #[classattr]
    const SECONDARY_SEQ: u64 = 0x800000000;
    #[classattr]
    const WEAK_PAIRING: u64 = 0x4000000000;
    #[classattr]
    const SR_RNA: u64 = 0x8000000000;
    #[classattr]
    const OUT_JUNC: u64 = 0x10000000000;
    #[classattr]
    const PAR_CHAIN: u64 = 0x20000000000;
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "SeedingParams")]
#[derive(Clone)]
pub struct PySeedingParams {
    #[pyo3(get, set)]
    pub mid_occ: usize,
    #[pyo3(get, set)]
    pub max_occ: usize,
    #[pyo3(get, set)]
    pub max_max_occ: usize,
    #[pyo3(get, set)]
    pub occ_dist: i32,
    #[pyo3(get, set)]
    pub q_occ_frac: f32,
    #[pyo3(get, set)]
    pub min_mid_occ: i32,
    #[pyo3(get, set)]
    pub max_mid_occ: i32,
}

impl From<SeedingParams> for PySeedingParams {
    fn from(c: SeedingParams) -> Self {
        Self {
            mid_occ: c.mid_occ,
            max_occ: c.max_occ,
            max_max_occ: c.max_max_occ,
            occ_dist: c.occ_dist,
            q_occ_frac: c.q_occ_frac,
            min_mid_occ: c.min_mid_occ,
            max_mid_occ: c.max_mid_occ,
        }
    }
}

impl From<PySeedingParams> for SeedingParams {
    fn from(p: PySeedingParams) -> Self {
        Self {
            mid_occ: p.mid_occ,
            max_occ: p.max_occ,
            max_max_occ: p.max_max_occ,
            occ_dist: p.occ_dist,
            q_occ_frac: p.q_occ_frac,
            min_mid_occ: p.min_mid_occ,
            max_mid_occ: p.max_mid_occ,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "ChainingParams")]
#[derive(Clone)]
pub struct PyChainingParams {
    #[pyo3(get, set)]
    pub min_cnt: i32,
    #[pyo3(get, set)]
    pub min_chain_score: i32,
    #[pyo3(get, set)]
    pub max_gap: i32,
    #[pyo3(get, set)]
    pub max_gap_ref: i32,
    #[pyo3(get, set)]
    pub max_dist_x: i32,
    #[pyo3(get, set)]
    pub max_dist_y: i32,
    #[pyo3(get, set)]
    pub bandwidth: i32,
    #[pyo3(get, set)]
    pub bandwidth_long: i32,
    #[pyo3(get, set)]
    pub max_chain_skip: i32,
    #[pyo3(get, set)]
    pub max_chain_iter: i32,
    #[pyo3(get, set)]
    pub chn_pen_gap: f32,
    #[pyo3(get, set)]
    pub chn_pen_skip: f32,
    #[pyo3(get, set)]
    pub chain_gap_scale: f32,
    #[pyo3(get, set)]
    pub rmq_rescue_size: i32,
    #[pyo3(get, set)]
    pub rmq_rescue_ratio: f32,
    #[pyo3(get, set)]
    pub rmq_inner_dist: i32,
    #[pyo3(get, set)]
    pub rmq_size_cap: i32,
}

impl From<ChainingParams> for PyChainingParams {
    fn from(c: ChainingParams) -> Self {
        Self {
            min_cnt: c.min_cnt,
            min_chain_score: c.min_chain_score,
            max_gap: c.max_gap,
            max_gap_ref: c.max_gap_ref,
            max_dist_x: c.max_dist_x,
            max_dist_y: c.max_dist_y,
            bandwidth: c.bandwidth,
            bandwidth_long: c.bandwidth_long,
            max_chain_skip: c.max_chain_skip,
            max_chain_iter: c.max_chain_iter,
            chn_pen_gap: c.chn_pen_gap,
            chn_pen_skip: c.chn_pen_skip,
            chain_gap_scale: c.chain_gap_scale,
            rmq_rescue_size: c.rmq_rescue_size,
            rmq_rescue_ratio: c.rmq_rescue_ratio,
            rmq_inner_dist: c.rmq_inner_dist,
            rmq_size_cap: c.rmq_size_cap,
        }
    }
}

impl From<PyChainingParams> for ChainingParams {
    fn from(p: PyChainingParams) -> Self {
        Self {
            min_cnt: p.min_cnt,
            min_chain_score: p.min_chain_score,
            max_gap: p.max_gap,
            max_gap_ref: p.max_gap_ref,
            max_dist_x: p.max_dist_x,
            max_dist_y: p.max_dist_y,
            bandwidth: p.bandwidth,
            bandwidth_long: p.bandwidth_long,
            max_chain_skip: p.max_chain_skip,
            max_chain_iter: p.max_chain_iter,
            chn_pen_gap: p.chn_pen_gap,
            chn_pen_skip: p.chn_pen_skip,
            chain_gap_scale: p.chain_gap_scale,
            rmq_rescue_size: p.rmq_rescue_size,
            rmq_rescue_ratio: p.rmq_rescue_ratio,
            rmq_inner_dist: p.rmq_inner_dist,
            rmq_size_cap: p.rmq_size_cap,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "ScoringParams")]
#[derive(Clone)]
pub struct PyScoringParams {
    #[pyo3(get, set)]
    pub match_score: i32,
    #[pyo3(get, set)]
    pub mismatch_penalty: i32,
    #[pyo3(get, set)]
    pub gap_open: i32,
    #[pyo3(get, set)]
    pub gap_extend: i32,
    #[pyo3(get, set)]
    pub gap_open2: i32,
    #[pyo3(get, set)]
    pub gap_extend2: i32,
    #[pyo3(get, set)]
    pub transition: i32,
    #[pyo3(get, set)]
    pub ambig_penalty: i32,
    #[pyo3(get, set)]
    pub noncanon_penalty: i32,
    #[pyo3(get, set)]
    pub junc_bonus: i32,
    #[pyo3(get, set)]
    pub junc_pen: i32,
}

impl From<ScoringParams> for PyScoringParams {
    fn from(c: ScoringParams) -> Self {
        Self {
            match_score: c.match_score,
            mismatch_penalty: c.mismatch_penalty,
            gap_open: c.gap_open,
            gap_extend: c.gap_extend,
            gap_open2: c.gap_open2,
            gap_extend2: c.gap_extend2,
            transition: c.transition,
            ambig_penalty: c.ambig_penalty,
            noncanon_penalty: c.noncanon_penalty,
            junc_bonus: c.junc_bonus,
            junc_pen: c.junc_pen,
        }
    }
}

impl From<PyScoringParams> for ScoringParams {
    fn from(p: PyScoringParams) -> Self {
        Self {
            match_score: p.match_score,
            mismatch_penalty: p.mismatch_penalty,
            gap_open: p.gap_open,
            gap_extend: p.gap_extend,
            gap_open2: p.gap_open2,
            gap_extend2: p.gap_extend2,
            transition: p.transition,
            ambig_penalty: p.ambig_penalty,
            noncanon_penalty: p.noncanon_penalty,
            junc_bonus: p.junc_bonus,
            junc_pen: p.junc_pen,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "AlignmentParams")]
#[derive(Clone)]
pub struct PyAlignmentParams {
    #[pyo3(get, set)]
    pub zdrop: i32,
    #[pyo3(get, set)]
    pub zdrop_inv: i32,
    #[pyo3(get, set)]
    pub end_bonus: i32,
    #[pyo3(get, set)]
    pub max_sw_mat: i64,
    #[pyo3(get, set)]
    pub min_dp_max: i32,
    #[pyo3(get, set)]
    pub min_dp_len: i32,
    #[pyo3(get, set)]
    pub anchor_ext_len: i32,
    #[pyo3(get, set)]
    pub anchor_ext_shift: i32,
    #[pyo3(get, set)]
    pub max_clip_ratio: f32,
}

impl From<AlignmentParams> for PyAlignmentParams {
    fn from(c: AlignmentParams) -> Self {
        Self {
            zdrop: c.zdrop,
            zdrop_inv: c.zdrop_inv,
            end_bonus: c.end_bonus,
            max_sw_mat: c.max_sw_mat,
            min_dp_max: c.min_dp_max,
            min_dp_len: c.min_dp_len,
            anchor_ext_len: c.anchor_ext_len,
            anchor_ext_shift: c.anchor_ext_shift,
            max_clip_ratio: c.max_clip_ratio,
        }
    }
}

impl From<PyAlignmentParams> for AlignmentParams {
    fn from(p: PyAlignmentParams) -> Self {
        Self {
            zdrop: p.zdrop,
            zdrop_inv: p.zdrop_inv,
            end_bonus: p.end_bonus,
            max_sw_mat: p.max_sw_mat,
            min_dp_max: p.min_dp_max,
            min_dp_len: p.min_dp_len,
            anchor_ext_len: p.anchor_ext_len,
            anchor_ext_shift: p.anchor_ext_shift,
            max_clip_ratio: p.max_clip_ratio,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "FilteringParams")]
#[derive(Clone)]
pub struct PyFilteringParams {
    #[pyo3(get, set)]
    pub best_n: i32,
    #[pyo3(get, set)]
    pub pri_ratio: f32,
    #[pyo3(get, set)]
    pub mask_level: f32,
    #[pyo3(get, set)]
    pub mask_len: i32,
    #[pyo3(get, set)]
    pub is_splice: bool,
    #[pyo3(get, set)]
    pub alt_drop: f32,
    #[pyo3(get, set)]
    pub seed: i32,
    #[pyo3(get, set)]
    pub chain_skip_scale: f32,
    #[pyo3(get, set)]
    pub max_qlen: i32,
    #[pyo3(get, set)]
    pub jump_min_match: i32,
}

impl From<FilteringParams> for PyFilteringParams {
    fn from(c: FilteringParams) -> Self {
        Self {
            best_n: c.best_n,
            pri_ratio: c.pri_ratio,
            mask_level: c.mask_level,
            mask_len: c.mask_len,
            is_splice: c.is_splice,
            alt_drop: c.alt_drop,
            seed: c.seed,
            chain_skip_scale: c.chain_skip_scale,
            max_qlen: c.max_qlen,
            jump_min_match: c.jump_min_match,
        }
    }
}

impl From<PyFilteringParams> for FilteringParams {
    fn from(p: PyFilteringParams) -> Self {
        Self {
            best_n: p.best_n,
            pri_ratio: p.pri_ratio,
            mask_level: p.mask_level,
            mask_len: p.mask_len,
            is_splice: p.is_splice,
            alt_drop: p.alt_drop,
            seed: p.seed,
            chain_skip_scale: p.chain_skip_scale,
            max_qlen: p.max_qlen,
            jump_min_match: p.jump_min_match,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "PairedEndParams")]
#[derive(Clone)]
pub struct PyPairedEndParams {
    #[pyo3(get, set)]
    pub max_frag_len: i32,
    #[pyo3(get, set)]
    pub pe_ori: i32,
    #[pyo3(get, set)]
    pub pe_bonus: i32,
}

impl From<PairedEndParams> for PyPairedEndParams {
    fn from(c: PairedEndParams) -> Self {
        Self {
            max_frag_len: c.max_frag_len,
            pe_ori: c.pe_ori,
            pe_bonus: c.pe_bonus,
        }
    }
}

impl From<PyPairedEndParams> for PairedEndParams {
    fn from(p: PyPairedEndParams) -> Self {
        Self {
            max_frag_len: p.max_frag_len,
            pe_ori: p.pe_ori,
            pe_bonus: p.pe_bonus,
        }
    }
}

#[gen_stub_pyclass]
#[pyclass(module = "rammappy._rammappy", name = "MapOptions")]
pub struct PyMapOptions {
    #[pyo3(get, set)]
    pub seeding: Py<PySeedingParams>,
    #[pyo3(get, set)]
    pub chaining: Py<PyChainingParams>,
    #[pyo3(get, set)]
    pub scoring: Py<PyScoringParams>,
    #[pyo3(get, set)]
    pub alignment: Py<PyAlignmentParams>,
    #[pyo3(get, set)]
    pub filtering: Py<PyFilteringParams>,
    #[pyo3(get, set)]
    pub pairing: Py<PyPairedEndParams>,
    /// Raw bitset preserving every rammap behavioral flag.
    #[pyo3(get, set)]
    pub flags: u64,
    #[pyo3(get, set)]
    pub mini_batch_size: i64,
}

impl PyMapOptions {
    pub fn from_map(py: Python<'_>, c: MapOptions) -> PyResult<Self> {
        Ok(Self {
            seeding: Py::new(py, PySeedingParams::from(c.seeding))?,
            chaining: Py::new(py, PyChainingParams::from(c.chaining))?,
            scoring: Py::new(py, PyScoringParams::from(c.scoring))?,
            alignment: Py::new(py, PyAlignmentParams::from(c.alignment))?,
            filtering: Py::new(py, PyFilteringParams::from(c.filtering))?,
            pairing: Py::new(py, PyPairedEndParams::from(c.pairing))?,
            flags: c.flags.bits(),
            mini_batch_size: c.mini_batch_size,
        })
    }

    pub fn into_map(&self, py: Python<'_>) -> MapOptions {
        MapOptions {
            seeding: self.seeding.bind(py).borrow().clone().into(),
            chaining: self.chaining.bind(py).borrow().clone().into(),
            scoring: self.scoring.bind(py).borrow().clone().into(),
            alignment: self.alignment.bind(py).borrow().clone().into(),
            filtering: self.filtering.bind(py).borrow().clone().into(),
            pairing: self.pairing.bind(py).borrow().clone().into(),
            flags: AlignFlags::from_bits_retain(self.flags),
            mini_batch_size: self.mini_batch_size,
        }
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySeedingParams>()?;
    m.add_class::<PyChainingParams>()?;
    m.add_class::<PyScoringParams>()?;
    m.add_class::<PyAlignmentParams>()?;
    m.add_class::<PyFilteringParams>()?;
    m.add_class::<PyPairedEndParams>()?;
    m.add_class::<PyAlignFlags>()?;
    m.add_class::<PyMapOptions>()?;
    Ok(())
}

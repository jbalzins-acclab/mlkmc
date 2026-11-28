# Changelog

All notable changes on the `perf-opt` branch of `mlkmc` are documented in this file.

## [Unreleased]

### Fixed
- **Prevent `minim_cache` cross-contamination and isolated-monomer poisoning (`kmc.f90`):**
  - **Root cause:** `create_fingerprint` hashed only the neighbor shell counts around `final_id` (excluding `site_id`) without encoding `jump_shell`, and `get_rates()` / `get_rates_swiftmd_batch()` called `check_minim_cache` *before* `only_substrate_lae`. Whenever any longer or multi-barrier jump (`1nn`, `over-top`, `2nn`, `zigzag`, `3nn`, or detachment) onto an empty terrace `fcc` or `hcp` site failed final-state minimization (`too large displacement`) or NEB validation (`last_neb_forbid` / `reverse_barrier <= 0`), it poisoned the bare-terrace monomer fingerprint `(1, 0, 0, 0, 3, ...)` in `minim_cache`. Every isolated monomer's `fcc <-> hcp` terrace hop then hit `MINIM_RESULT_FAIL` before reaching `only_substrate_lae` and was marked `EVENT_FORBIDDEN (-3)` (`hop_rate = 0`), leaving `470/470` tabGAP monomers and `430/431` EAM monomers frozen at 500 K and collapsing tabGAP `total_rate` down to `dep_rate` (`100,440 Hz`).
  - **Fix:**
    - Evaluate `only_substrate_lae` *before* `check_minim_cache` in `get_rates()`, `get_rates_swiftmd_batch()`, and `execute_jump_event()`; bypass `minim_cache` reads and writes (`cache_minim=(only_substrate_jump == SUBSTRATE_NOT)`) for bare-substrate `fcc <-> hcp` jumps.
    - Encode `jump_shell` into `fingerprint(atom_type, 0) = 1 + 10 * jump_shell` in `create_fingerprint()` (`1` when omitted for exchange events) so failures on longer jump shells can never poison shorter `fcc <-> hcp` hops or exchanges.

### Changed
- **Carved-out-only NEB path dumping by default (`swiftmd_mod.f90`):**
  - Default `MLKMC_SWIFTMD_SAVE_NEB` to `"neb_paths"` when unset (`0`/`off`/`false`/`none` disables) and default `MLKMC_SWIFTMD_SAVE_NEB_FULL` to `.false.` so `<prefix>_active.xyz` and `highlight_<prefix>_active.xyz` are saved without writing multi-GB `<prefix>_full.bin` full-slab trajectories unless opted in (`MLKMC_SWIFTMD_SAVE_NEB_FULL=1`).
  - Relabel substrate atoms (`frozen_substrate` + `substrate`) to `Pd` (`Z = 46`) in `highlight_<prefix>_active.xyz` (with deposited adatoms in true species `Ag`, jumping atom `Au`, exchange middle partner `Cu`) for clear contrast in OVITO.

### Added
- **SwiftMD persistent C-ABI engine & NEB backend (`swiftmd_mod.f90`, `mlkmc.f90`, `kmc.f90`, `MAKE/makefile.defs`):**
  - Support `eam/fs` and `tabgap` (`hybrid/overlay tabgap eam/fs`) potentials on CPU (`gpu=0`) and Metal GPU (`gpu=1`) through `libSwiftMDC`, with automatic backend selection and `backend` / `MLKMC_BACKEND` (`swiftmd` | `lammps`) overrides.
  - Configurable NEB controls via input file and environment variables (`MLKMC_SWIFTMD_ACTIVE_REGION`, `MLKMC_SWIFTMD_GPU`, `MLKMC_SWIFTMD_OPTIMIZER`, `MLKMC_SWIFTMD_FTOL`, `MLKMC_SWIFTMD_CLIMB`, `MLKMC_SWIFTMD_ENDPOINT_MODE`, `MLKMC_SWIFTMD_PATH_INIT`, `MLKMC_SWIFTMD_DYNEB`, `MLKMC_SWIFTMD_AUTO_CPU`, `MLKMC_SWIFTMD_NEB_BATCH`, `MLKMC_SWIFTMD_DEBUG`).
  - Experimental batched NEB scheduling path (`get_rates_swiftmd_batch`, `prepare_barrier_or_queue`, `flush_neb_batch`) with deduplication/aliasing.
  - NEB quality gates: `MLKMC_NEB_REJECT_FMAX` (`<= 0` disabled) keeps high-residual unconverged bands out of ML training sets (`train_ok`), and `MLKMC_NEB_FORBID_FMAX` (default `0.5 eV/A`) forbids and caches broken/exploded unconverged NEBs alongside `reverse_barrier <= 0` events.
- **JGAP single-center SOAP backend (`jgap_soap_mod.f90`, `kmc.f90`, `mlkmc.f90`):**
  - Persistent `jgap_c` SOAP evaluator operating directly on contiguous LAE buffers (`get_jump_lae_buf`) with `MLKMC_JGAP_SOAP_VERIFY` verification mode against QUIP SOAP.
- **Hybrid potential parsing (`mlkmc.f90`):**
  - Widen `pot_name` / `pot_style` buffers and support multi-file `hybrid` / `hybrid/overlay` `pair_coeff` generation for LAMMPS fallback and SwiftMD handoff.

## [7a79405] - 2026-11-25

### Fixed
- Align `size_up_int_2D`, `size_up_real_2D`, `size_up_double_2D`, `size_up_real_3D`, `size_up_int_4D`, and `size_up_logical` in `utils.f90` with 1D exponential (`2x`) growth so co-resized arrays (`predictor%x` vs `%y`, etc.) cannot desync and corrupt the heap.
- Clear leftover `site_blocked` entries at `get_rates()` entry (`kmc.f90`) and guard `first_dynamic_site` in `delete_site()` (`atoms.f90`).

## [0ab4537] - 2026-11-25

### Changed
- Raise default Fortran optimization flags from `-O2` to `-O3 -march=native` without `-ffast-math` across `MAKE/makefile.defs` and target makefiles (preserving IEEE-754 reproducibility).

## [356dcd7] - 2026-11-25

### Changed
- Optimize KMC rate evaluation, lockedness bookkeeping, dynamic-site tracking, and system setup while preserving bit-identical output (`~29.4x` total speedup, `~64.1x` non-LAMMPS over 300k steps):
  - Park non-NEB MPI workers on a single `mpi_bcast(worker_cmd)` loop in `run_kmc()`.
  - Replace per-step `O(n_sites)` `nn_tmp`/`n_near_tmp` allocations and strided `t_site` copies with `O(1)` local shell math, persistent scratch arrays, sparse `blocked_sites` clearing, lazy config construction, and sorted `deposited_atoms` summation.
  - Replace `O(n_sites)` `lock_checked` sweeps with a sparse `lock_dirty_sites` list (`atoms_types.f90`, `atoms.f90`, `kmc.f90`).
  - Restrict per-step `inactive_counter` scans to dynamic sites (`>= first_dynamic_site`).
  - Pre-allocate `system%sites` in `create_system()` and use `2x` growth + `move_alloc` in `size_up_sites`.

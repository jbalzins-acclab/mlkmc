! MLKMC - Fortran iso_c_binding interface to SwiftMD C ABI (swiftmd.h)
!
! Provides persistent engine management, zero-reallocation coordinate/species
! caching across KMC steps, single-structure minimization, and NEB barrier calls.

module swiftmd_module
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env, only : int8, int32, int64, real32, real64
  use constants_module
  use mlkmc_atoms_types_module
  use libAtoms_module, only : ElementName
  implicit none

  private
  public :: t_swiftmd_ctx
  public :: swiftmd_init_ctx
  public :: swiftmd_destroy_ctx
  public :: swiftmd_sync_atoms
  public :: swiftmd_minimize_cfg
  public :: swiftmd_neb_cfg
  public :: swiftmd_neb_batch_cfg
  public :: swiftmd_is_supported_potential

  interface
    function c_swiftmd_create(potential_path, eam_path, species, is_fixed, n_atoms, cell_3x3, pbc_3, use_gpu) &
        bind(C, name="swiftmd_create") result(handle)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: potential_path
      type(c_ptr), value :: eam_path
      type(c_ptr), intent(in) :: species(*)
      integer(c_int), intent(in) :: is_fixed(*)
      integer(c_int), value :: n_atoms
      real(c_double), intent(in) :: cell_3x3(9)
      integer(c_int), intent(in) :: pbc_3(3)
      integer(c_int), value :: use_gpu
      type(c_ptr) :: handle
    end function c_swiftmd_create

    subroutine c_swiftmd_destroy(handle) bind(C, name="swiftmd_destroy")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine c_swiftmd_destroy

    function c_swiftmd_update_atoms(handle, species, is_fixed, n_atoms) &
        bind(C, name="swiftmd_update_atoms") result(rc)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      type(c_ptr), intent(in) :: species(*)
      integer(c_int), intent(in) :: is_fixed(*)
      integer(c_int), value :: n_atoms
      integer(c_int) :: rc
    end function c_swiftmd_update_atoms

    function c_swiftmd_set_neb_endpoint_mode(handle, endpoint_mode) &
        bind(C, name="swiftmd_set_neb_endpoint_mode") result(rc)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int), value :: endpoint_mode
      integer(c_int) :: rc
    end function c_swiftmd_set_neb_endpoint_mode

    function c_swiftmd_set_active_region(handle, radius_angstrom) &
        bind(C, name="swiftmd_set_active_region") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), value :: radius_angstrom
      integer(c_int) :: rc
    end function c_swiftmd_set_active_region

    function c_swiftmd_set_neb_optimizer(handle, optimizer) &
        bind(C, name="swiftmd_set_neb_optimizer") result(rc)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int), value :: optimizer
      integer(c_int) :: rc
    end function c_swiftmd_set_neb_optimizer

    function c_swiftmd_set_neb_options(handle, path_init, weighted_springs, dyneb) &
        bind(C, name="swiftmd_set_neb_options") result(rc)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int), value :: path_init, weighted_springs, dyneb
      integer(c_int) :: rc
    end function c_swiftmd_set_neb_options

    function c_swiftmd_set_auto_cpu_refinement(handle, enabled) &
        bind(C, name="swiftmd_set_auto_cpu_refinement") result(rc)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int), value :: enabled
      integer(c_int) :: rc
    end function c_swiftmd_set_auto_cpu_refinement

    function c_swiftmd_minimize(handle, positions_inout, etol, ftol, dmax, max_steps, &
        energy_out, fmax_out, steps_out, converged_out) bind(C, name="swiftmd_minimize") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(inout) :: positions_inout(*)
      real(c_double), value :: etol, ftol, dmax
      integer(c_int), value :: max_steps
      type(c_ptr), value :: energy_out, fmax_out, steps_out, converged_out
      integer(c_int) :: rc
    end function c_swiftmd_minimize

    function c_swiftmd_neb(handle, init_positions, final_positions, n_images, spring_k, climb, &
        etol, ftol, dmax, max_steps, pe_out, relaxed_init_out, relaxed_final_out, &
        fmax_out, steps_out, converged_out) bind(C, name="swiftmd_neb") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: init_positions(*)
      real(c_double), intent(in) :: final_positions(*)
      integer(c_int), value :: n_images
      real(c_double), value :: spring_k
      integer(c_int), value :: climb
      real(c_double), value :: etol, ftol, dmax
      integer(c_int), value :: max_steps
      type(c_ptr), value :: pe_out, relaxed_init_out, relaxed_final_out, fmax_out, steps_out, converged_out
      integer(c_int) :: rc
    end function c_swiftmd_neb

    function c_swiftmd_neb_images(handle, init_positions, final_positions, n_images, spring_k, climb, &
        etol, ftol, dmax, max_steps, pe_out, relaxed_init_out, relaxed_final_out, &
        fmax_out, steps_out, converged_out, images_out) bind(C, name="swiftmd_neb_images") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: init_positions(*)
      real(c_double), intent(in) :: final_positions(*)
      integer(c_int), value :: n_images
      real(c_double), value :: spring_k
      integer(c_int), value :: climb
      real(c_double), value :: etol, ftol, dmax
      integer(c_int), value :: max_steps
      type(c_ptr), value :: pe_out, relaxed_init_out, relaxed_final_out, fmax_out, steps_out, converged_out
      type(c_ptr), value :: images_out
      integer(c_int) :: rc
    end function c_swiftmd_neb_images

    function c_swiftmd_neb_batch(handle, n_bands, init_positions_batch, final_positions_batch, &
        n_images, spring_k, climb, etol, ftol, dmax, max_steps, &
        barriers_fwd_out, barriers_rev_out, pe_out, fmax_out, steps_out, converged_out) &
        bind(C, name="swiftmd_neb_batch") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      integer(c_int), value :: n_bands
      real(c_double), intent(in) :: init_positions_batch(*)
      real(c_double), intent(in) :: final_positions_batch(*)
      integer(c_int), value :: n_images
      real(c_double), value :: spring_k
      integer(c_int), value :: climb
      real(c_double), value :: etol, ftol, dmax
      integer(c_int), value :: max_steps
      real(c_double), intent(out) :: barriers_fwd_out(*)
      real(c_double), intent(out) :: barriers_rev_out(*)
      type(c_ptr), value :: pe_out
      real(c_double), intent(out) :: fmax_out(*)
      integer(c_int), intent(out) :: steps_out(*)
      integer(c_int), intent(out) :: converged_out(*)
      integer(c_int) :: rc
    end function c_swiftmd_neb_batch

    function c_swiftmd_last_error() bind(C, name="swiftmd_last_error") result(err_ptr)
      import :: c_ptr
      type(c_ptr) :: err_ptr
    end function c_swiftmd_last_error
  end interface

  type t_swiftmd_ctx
    type(c_ptr) :: handle = c_null_ptr
    logical :: enabled = .false.
    integer :: n_substrate = 0
    integer :: n_atoms = 0
    integer :: capacity = 0
    integer :: use_gpu = 0
    ! 0 = free endpoints (default): end images relax under their true forces during the NEB, exactly like
    !     LAMMPS 'fix neb' without the 'end' keyword. mlkmc hands over UNRELAXED lattice-site endpoints,
    !     so the endpoints must relax; 1 = fixed endpoints would pin the band to unrelaxed states.
    integer :: endpoint_mode = 0
    integer :: climb = 1
    double precision :: active_region = 15.0d0
    double precision :: ftol_override = -1.0d0
    double precision :: min_etol = 0.0d0
    double precision :: neb_etol = 0.0d0
    real(c_double) :: cell_3x3(9) = 0.0_c_double
    integer(c_int) :: pbc_3(3) = [1_c_int, 1_c_int, 0_c_int]
    character(kind=c_char, len=1), allocatable :: species_chars(:,:)
    type(c_ptr), allocatable :: species_ptrs(:)
    integer(c_int), allocatable :: is_fixed(:)
    real(c_double), allocatable :: pos_init(:,:)
    real(c_double), allocatable :: pos_final(:,:)
    real(c_double), allocatable :: pe_buf(:)
    logical :: substrate_pos_loaded = .false.
    logical :: sub_mask_dirty = .false.
    logical :: debug = .false.
    logical :: use_neb_batch = .false.
    integer(int64) :: neb_batch_calls = 0_int64
    integer(int64) :: neb_batch_bands = 0_int64
    real(c_double), allocatable :: batch_init(:,:,:)
    real(c_double), allocatable :: batch_final(:,:,:)
    ! NEB path saving (MLKMC_SWIFTMD_SAVE_NEB=<prefix>, default "neb_paths"):
    !   <prefix>_active.xyz             active-region cutout, full slab box, true species only
    !   highlight_<prefix>_active.xyz   active cutout centered at (0,0,0) with tight box for OVITO;
    !                                   substrate relabelled Pd (Z=46), adatoms Ag, jump Au, exchange Cu
    !   <prefix>_full.bin               (opt-in via MLKMC_SWIFTMD_SAVE_NEB_FULL=1 only; huge)
    logical :: save_neb = .false.
    logical :: save_neb_full = .false.
    integer :: neb_unit_full = -1
    integer :: neb_unit_active = -1
    integer :: neb_unit_highlight = -1
    integer :: neb_saved = 0
    real(c_double), allocatable :: img_buf(:,:)
  end type t_swiftmd_ctx

  integer, parameter :: NEB_MARK_Z_SUBSTRATE = 46 ! Pd
  integer, parameter :: NEB_MARK_Z_JUMP = 79      ! Au
  integer, parameter :: NEB_MARK_Z_MIDDLE = 29    ! Cu

contains

  subroutine set_c_string(str, buf)
    character(len=*), intent(in) :: str
    character(kind=c_char, len=1), intent(out) :: buf(:)
    integer :: i, n

    buf = c_null_char
    n = min(len_trim(str), size(buf) - 1)
    do i = 1, n
      buf(i) = str(i:i)
    end do
    buf(n + 1) = c_null_char
  end subroutine set_c_string

  subroutine print_swiftmd_error(prefix)
    character(len=*), intent(in) :: prefix
    type(c_ptr) :: p
    character(kind=c_char), pointer :: chars(:)
    integer :: n, i
    character(len=512) :: msg

    p = c_swiftmd_last_error()
    if (.not. c_associated(p)) then
      write(0, "(a)") "*** " // trim(prefix) // ": (no error message)"
      return
    end if
    call c_f_pointer(p, chars, [511])
    msg = ""
    n = 0
    do i = 1, 511
      if (chars(i) == c_null_char) exit
      n = i
      msg(i:i) = chars(i)
    end do
    write(0, "(a)") "*** " // trim(prefix) // ": " // trim(msg(1:n))
  end subroutine print_swiftmd_error

  logical function swiftmd_is_supported_potential(pot_style, pot_files, n_pot_files, &
                                                  substrate_type, potential_path, eam_path)
    character(len=*), intent(in) :: pot_style
    character(len=*), intent(in) :: pot_files(:)
    integer, intent(in) :: n_pot_files, substrate_type
    character(len=*), intent(out) :: potential_path, eam_path

    integer :: i, l
    character(len=256) :: tok, style_l

    swiftmd_is_supported_potential = .false.
    potential_path = ""
    eam_path = ""

    if (substrate_type == SUBSTRATE_TYPE_WEAK .or. n_pot_files < 1) return
    style_l = adjustl(pot_style)

    if (index(style_l, "eam/fs") > 0 .and. index(style_l, "tabgap") == 0 .and. &
        index(style_l, "hybrid") == 0) then
      if (n_pot_files == 1) then
        potential_path = trim(pot_files(1))
        swiftmd_is_supported_potential = .true.
      end if
      return
    end if

    if (index(style_l, "tabgap") > 0) then
      do i = 1, n_pot_files
        tok = trim(adjustl(pot_files(i)))
        l = len_trim(tok)
        if ((l >= 7 .and. tok(l-6:l) == ".eam.fs") .or. (l >= 4 .and. tok(l-3:l) == ".eam")) then
          eam_path = tok
        else
          potential_path = tok
        end if
      end do
      if (len_trim(potential_path) > 0) then
        swiftmd_is_supported_potential = .true.
      end if
    end if
  end function swiftmd_is_supported_potential

  subroutine ensure_capacity(ctx, needed)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    integer, intent(in) :: needed

    integer :: new_cap, i
    character(kind=c_char, len=1), allocatable, target :: new_sp(:,:)
    integer(c_int), allocatable :: new_fixed(:)
    real(c_double), allocatable :: new_pos_init(:,:), new_pos_final(:,:)

    if (needed <= ctx%capacity .and. allocated(ctx%species_chars)) return

    new_cap = max(needed + 256, max(1024, ctx%capacity * 2))
    allocate(new_sp(8, new_cap))
    new_sp = c_null_char
    allocate(new_fixed(new_cap))
    new_fixed = 0_c_int
    allocate(new_pos_init(3, new_cap))
    new_pos_init = 0.0_c_double
    allocate(new_pos_final(3, new_cap))
    new_pos_final = 0.0_c_double

    if (allocated(ctx%species_chars) .and. ctx%n_atoms > 0) then
      new_sp(:, 1:ctx%n_atoms) = ctx%species_chars(:, 1:ctx%n_atoms)
      new_fixed(1:ctx%n_atoms) = ctx%is_fixed(1:ctx%n_atoms)
      new_pos_init(:, 1:ctx%n_atoms) = ctx%pos_init(:, 1:ctx%n_atoms)
      new_pos_final(:, 1:ctx%n_atoms) = ctx%pos_final(:, 1:ctx%n_atoms)
    end if

    call move_alloc(new_sp, ctx%species_chars)
    call move_alloc(new_fixed, ctx%is_fixed)
    call move_alloc(new_pos_init, ctx%pos_init)
    call move_alloc(new_pos_final, ctx%pos_final)

    if (allocated(ctx%species_ptrs)) deallocate(ctx%species_ptrs)
    allocate(ctx%species_ptrs(new_cap))
    do i = 1, new_cap
      ctx%species_ptrs(i) = c_loc(ctx%species_chars(1, i))
    end do
    ctx%capacity = new_cap
  end subroutine ensure_capacity

  subroutine fill_substrate_metadata(ctx, system)
    type(t_swiftmd_ctx), intent(inout) :: ctx
    type(t_system), intent(in) :: system

    integer :: i, f_start, f_end, s_start, s_end

    f_start = system%frozen_substrate%first_id
    f_end = system%frozen_substrate%last_id
    s_start = system%substrate%first_id
    s_end = system%substrate%last_id

    do i = f_start, f_end
      call set_c_string(ElementName(system%frozen_substrate%Z(i)), ctx%species_chars(:, i))
      ctx%is_fixed(i) = 1_c_int
      ctx%pos_init(:, i) = real(system%frozen_substrate%pos(:, i), c_double)
    end do
    do i = s_start, s_end
      call set_c_string(ElementName(system%substrate%Z(i)), ctx%species_chars(:, i))
      ctx%is_fixed(i) = 0_c_int
      ctx%pos_init(:, i) = real(system%substrate%pos(:, i), c_double)
    end do
    ctx%pos_final(:, 1:system%n_substrate) = ctx%pos_init(:, 1:system%n_substrate)
    ctx%substrate_pos_loaded = .true.
  end subroutine fill_substrate_metadata

  subroutine swiftmd_init_ctx(ctx, potential_path, eam_path, system, use_gpu, active_region, error)
    type(t_swiftmd_ctx), intent(inout) :: ctx
    character(len=*), intent(in) :: potential_path, eam_path
    type(t_system), intent(in) :: system
    integer, intent(in) :: use_gpu
    double precision, intent(in) :: active_region
    integer, intent(out) :: error

    integer :: i, site_id, atom_id, n_tot, dbg_len, dbg_stat
    integer(c_int) :: rc
    character(kind=c_char, len=1), target :: pot_c(512), eam_c(512)
    character(len=16) :: dbg_env
    type(c_ptr) :: eam_ptr

    error = 0
    call get_environment_variable("MLKMC_SWIFTMD_DEBUG", dbg_env, length=dbg_len, status=dbg_stat)
    ctx%debug = (dbg_stat == 0 .and. dbg_len > 0)
    ctx%use_gpu = use_gpu
    ctx%active_region = active_region
    ctx%n_substrate = system%n_substrate
    n_tot = system%n_substrate + system%n_deposited

    call ensure_capacity(ctx, n_tot)
    ctx%n_atoms = n_tot

    call fill_substrate_metadata(ctx, system)
    do i = 1, system%n_deposited
      site_id = system%deposited_atoms(i)
      atom_id = system%sites(site_id)%id
      call set_c_string(ElementName(system%sites(site_id)%Z), ctx%species_chars(:, atom_id))
      ctx%is_fixed(atom_id) = 0_c_int
      ctx%pos_init(:, atom_id) = real(system%sites(site_id)%cartesian_coords(:), c_double)
    end do

    ctx%cell_3x3 = reshape(real(system%lattice(1:3, 1:3), c_double), [9])
    do i = 1, 3
      if (system%periodic(i)) then
        ctx%pbc_3(i) = 1_c_int
      else
        ctx%pbc_3(i) = 0_c_int
      end if
    end do

    call set_c_string(potential_path, pot_c)
    if (len_trim(eam_path) > 0) then
      call set_c_string(eam_path, eam_c)
      eam_ptr = c_loc(eam_c(1))
    else
      eam_ptr = c_null_ptr
    end if

    ctx%handle = c_swiftmd_create(c_loc(pot_c(1)), eam_ptr, ctx%species_ptrs, ctx%is_fixed, &
                                  int(n_tot, c_int), ctx%cell_3x3, ctx%pbc_3, int(use_gpu, c_int))
    if (.not. c_associated(ctx%handle)) then
      call print_swiftmd_error("swiftmd_create failed")
      error = 1
      return
    end if

    block
      character(len=32) :: ep_val
      integer :: ep_len, ep_stat
      call get_environment_variable("MLKMC_SWIFTMD_ENDPOINT_MODE", ep_val, length=ep_len, status=ep_stat)
      if (ep_stat == 0 .and. ep_len > 0) read(ep_val, *) ctx%endpoint_mode
    end block
    write(6, "(a,i0,a)") "# SwiftMD NEB endpoint mode: ", ctx%endpoint_mode, " (0 = free as in LAMMPS fix neb, 1 = fixed)"
    rc = c_swiftmd_set_neb_endpoint_mode(ctx%handle, int(ctx%endpoint_mode, c_int))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_set_neb_endpoint_mode failed")
      error = 1
      return
    end if

    rc = c_swiftmd_set_active_region(ctx%handle, real(ctx%active_region, c_double))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_set_active_region failed")
      error = 1
      return
    end if

    block
      character(len=32) :: env_val
      integer :: env_len, env_stat, ios, val, pinit, wspr, dyn
      call get_environment_variable("MLKMC_SWIFTMD_AUTO_CPU", env_val, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_val, *, iostat=ios) val
        if (ios == 0) then
          rc = c_swiftmd_set_auto_cpu_refinement(ctx%handle, int(val, c_int))
          if (rc /= 0_c_int) then
            call print_swiftmd_error("swiftmd_set_auto_cpu_refinement failed")
            error = 1
            return
          end if
        end if
      end if
      call get_environment_variable("MLKMC_SWIFTMD_OPTIMIZER", env_val, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_val, *, iostat=ios) val
        if (ios == 0) then
          rc = c_swiftmd_set_neb_optimizer(ctx%handle, int(val, c_int))
          if (rc /= 0_c_int) then
            call print_swiftmd_error("swiftmd_set_neb_optimizer failed")
            error = 1
            return
          end if
        end if
      end if
      pinit = 1; wspr = 0; dyn = 1
      ios = 0
      call get_environment_variable("MLKMC_SWIFTMD_PATH_INIT", env_val, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_val, *, iostat=ios) val
        if (ios == 0) pinit = val
      end if
      call get_environment_variable("MLKMC_SWIFTMD_DYNEB", env_val, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        read(env_val, *, iostat=ios) val
        if (ios == 0) dyn = val
      end if
      call get_environment_variable("MLKMC_SWIFTMD_NEB_BATCH", env_val, length=env_len, status=env_stat)
      if (env_stat == 0 .and. env_len > 0) then
        ctx%use_neb_batch = any(env_val(1:1) == ["1", "T", "t", "Y", "y"])
      end if
      if (pinit /= 1 .or. dyn /= 1) then
        rc = c_swiftmd_set_neb_options(ctx%handle, int(pinit, c_int), int(wspr, c_int), int(dyn, c_int))
        if (rc /= 0_c_int) then
          call print_swiftmd_error("swiftmd_set_neb_options failed")
          error = 1
          return
        end if
      end if
    end block

    if (.not. allocated(ctx%pe_buf)) allocate(ctx%pe_buf(64))

    call open_neb_save_files(ctx, error)
    if (error /= 0) return

    ctx%enabled = .true.
  end subroutine swiftmd_init_ctx

  subroutine open_neb_save_files(ctx, error)
    type(t_swiftmd_ctx), intent(inout) :: ctx
    integer, intent(out) :: error

    character(len=512) :: env_val, full_val, prefix, fname, hfname
    integer :: env_len, env_stat, full_len, full_stat, ios, slash_pos
    logical :: existed

    error = 0
    call get_environment_variable("MLKMC_SWIFTMD_SAVE_NEB", env_val, length=env_len, status=env_stat)
    if (env_stat /= 0 .or. env_len == 0) then
      ! Default ON: user always needs NEB paths saved unless explicitly disabled.
      prefix = "neb_paths"
    else
      prefix = trim(adjustl(env_val))
      select case (trim(prefix))
      case ("0", "off", "OFF", "false", "FALSE", "F", "f", "none", "NONE")
        return
      case ("1")
        prefix = "neb_paths"
      end select
    end if

    ctx%save_neb_full = .false.
    call get_environment_variable("MLKMC_SWIFTMD_SAVE_NEB_FULL", full_val, length=full_len, status=full_stat)
    if (full_stat == 0 .and. full_len > 0) then
      select case (trim(adjustl(full_val)))
      case ("1", "on", "ON", "true", "TRUE", "T", "t")
        ctx%save_neb_full = .true.
      end select
    end if

    ! Append to existing files so restarts from a checkpoint keep earlier paths.
    if (ctx%save_neb_full) then
      fname = trim(prefix)//"_full.bin"
      inquire(file=trim(fname), exist=existed)
      open(newunit=ctx%neb_unit_full, file=trim(fname), access="stream", form="unformatted", &
           status="unknown", position="append", action="write", iostat=ios)
      if (ios /= 0) then
        write(0, "(a)") "*** Error: couldn't open "//trim(fname)//" for NEB path saving"
        error = 1
        return
      end if
      if (.not. existed) then
        ! File header: magic, format version, cell (3x3, Angstrom), pbc flags
        write(ctx%neb_unit_full) "MLKMCNEB", int(1, int32), real(ctx%cell_3x3, real64), int(ctx%pbc_3, int32)
      end if
    end if

    fname = trim(prefix)//"_active.xyz"
    open(newunit=ctx%neb_unit_active, file=trim(fname), form="formatted", status="unknown", &
         position="append", action="write", iostat=ios)
    if (ios /= 0) then
      write(0, "(a)") "*** Error: couldn't open "//trim(fname)//" for NEB path saving"
      error = 1
      return
    end if

    slash_pos = index(trim(prefix), "/", back=.true.)
    if (slash_pos > 0) then
      hfname = prefix(1:slash_pos)//"highlight_"//trim(prefix(slash_pos+1:))//"_active.xyz"
    else
      hfname = "highlight_"//trim(prefix)//"_active.xyz"
    end if
    open(newunit=ctx%neb_unit_highlight, file=trim(hfname), form="formatted", status="unknown", &
         position="append", action="write", iostat=ios)
    if (ios /= 0) then
      write(0, "(a)") "*** Error: couldn't open "//trim(hfname)//" for NEB path saving"
      error = 1
      return
    end if

    ctx%save_neb = .true.
    if (ctx%save_neb_full) then
      write(6, "(a)") "# Saving NEB paths to "//trim(prefix)//"_full.bin, "//trim(fname)//", and "//trim(hfname)
    else
      write(6, "(a)") "# Saving NEB paths to "//trim(fname)//" and "//trim(hfname)//" (carved-out active region only)"
    end if
  end subroutine open_neb_save_files

  ! Builds per-atom true Z for the current configuration.
  subroutine neb_species(ctx, system, init_cfg, z)
    type(t_swiftmd_ctx), intent(in) :: ctx
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: init_cfg
    integer, intent(out) :: z(:)

    integer :: i

    z = 0
    do i = system%frozen_substrate%first_id, system%frozen_substrate%last_id
      z(i) = system%frozen_substrate%Z(i)
    end do
    do i = system%substrate%first_id, system%substrate%last_id
      z(i) = system%substrate%Z(i)
    end do
    do i = init_cfg%first_id, min(init_cfg%last_id, ctx%n_atoms)
      z(i) = init_cfg%Z(i)
    end do
  end subroutine neb_species

  ! Full-slab record layout (native little-endian, no padding):
  !   int32 record_id, int64 kmc_step, float64 kmc_time,
  !   int32 n_atoms, int32 n_images, int32 jumping_atom, int32 middle_atom (1-based, 0 = none),
  !   int32 site_id, int32 final_id, int32 converged, int32 neb_steps,
  !   float64 fmax, float64 barrier, float64 reverse_barrier,
  !   float64 energies(n_images), int8 Z(n_atoms), float32 pos(3, n_atoms, n_images)
  subroutine save_neb_record(ctx, system, init_cfg, jumping_atom_i, middle_atom_i, n_images, &
                             barrier, reverse_barrier, fmax, steps, converged, &
                             kmc_step, kmc_time, site_id, final_id)
    type(t_swiftmd_ctx), intent(inout) :: ctx
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: init_cfg
    integer, intent(in) :: jumping_atom_i, middle_atom_i, n_images, steps, converged
    integer, intent(in) :: kmc_step, site_id, final_id
    double precision, intent(in) :: barrier, reverse_barrier, fmax, kmc_time

    integer :: n, r, i, k, n_sel, idx
    integer, allocatable :: z(:), z_hl(:), sel_ids(:)
    logical, allocatable :: sel(:)
    double precision, allocatable :: hpos(:,:,:)
    double precision :: box(3), d(3), dc(3), r2, c_first(3), c_last(3), c_mid(3)
    double precision :: bmin(3), bmax(3), bmid(3), hbox(3), horig(3)
    double precision, parameter :: HBOX_PAD = 2.0d0
    character(len=1024) :: meta, header, hheader
    character(len=3) :: sym

    n = ctx%n_atoms
    ctx%neb_saved = ctx%neb_saved + 1
    allocate(z(n), z_hl(n))
    call neb_species(ctx, system, init_cfg, z)
    z_hl = z
    do i = system%frozen_substrate%first_id, system%frozen_substrate%last_id
      if (i >= 1 .and. i <= n) z_hl(i) = NEB_MARK_Z_SUBSTRATE
    end do
    do i = system%substrate%first_id, system%substrate%last_id
      if (i >= 1 .and. i <= n) z_hl(i) = NEB_MARK_Z_SUBSTRATE
    end do
    if (jumping_atom_i > 0) z_hl(jumping_atom_i) = NEB_MARK_Z_JUMP
    if (middle_atom_i > 0) z_hl(middle_atom_i) = NEB_MARK_Z_MIDDLE

    ! --- full slab, single precision (only when explicitly enabled via MLKMC_SWIFTMD_SAVE_NEB_FULL=1) ---
    if (ctx%save_neb_full .and. ctx%neb_unit_full >= 0) then
      write(ctx%neb_unit_full) int(ctx%neb_saved, int32), int(kmc_step, int64), kmc_time, &
        int(n, int32), int(n_images, int32), int(jumping_atom_i, int32), int(middle_atom_i, int32), &
        int(site_id, int32), int(final_id, int32), int(converged, int32), int(steps, int32), &
        fmax, barrier, reverse_barrier, real(ctx%pe_buf(1:n_images), real64), int(z, int8), &
        real(ctx%img_buf(:, 1:n*n_images), real32)
      flush(ctx%neb_unit_full)
    end if

    ! --- active region: atoms within R of the hopping atom's start or end position ---
    box = [ctx%cell_3x3(1), ctx%cell_3x3(5), ctx%cell_3x3(9)]
    r2 = ctx%active_region**2
    if (r2 <= 0.0d0) r2 = 15.0d0**2
    c_first = ctx%img_buf(:, jumping_atom_i)
    c_last = ctx%img_buf(:, (n_images-1)*n + jumping_atom_i)
    dc = c_last - c_first
    do k = 1, 3
      if (ctx%pbc_3(k) /= 0 .and. box(k) > 0.0d0) dc(k) = dc(k) - box(k)*anint(dc(k)/box(k))
    end do
    c_mid = c_first + 0.5d0 * dc

    allocate(sel(n))
    sel = .false.
    do i = 1, n
      d = ctx%img_buf(:, i) - c_first
      do k = 1, 3
        if (ctx%pbc_3(k) /= 0 .and. box(k) > 0.0d0) d(k) = d(k) - box(k)*anint(d(k)/box(k))
      end do
      if (sum(d*d) <= r2) then
        sel(i) = .true.
        cycle
      end if
      d = ctx%img_buf(:, (n_images-1)*n + i) - c_last
      do k = 1, 3
        if (ctx%pbc_3(k) /= 0 .and. box(k) > 0.0d0) d(k) = d(k) - box(k)*anint(d(k)/box(k))
      end do
      if (sum(d*d) <= r2) sel(i) = .true.
    end do
    n_sel = count(sel)
    allocate(sel_ids(n_sel), hpos(3, n_sel, n_images))
    idx = 0
    do i = 1, n
      if (.not. sel(i)) cycle
      idx = idx + 1
      sel_ids(idx) = i
    end do

    ! Unwrap active region across minimum-image PBC around c_mid and center bounding box at (0,0,0).
    bmin = huge(1.0d0)
    bmax = -huge(1.0d0)
    do r = 1, n_images
      do idx = 1, n_sel
        i = sel_ids(idx)
        d = ctx%img_buf(:, (r-1)*n + i) - c_mid
        do k = 1, 3
          if (ctx%pbc_3(k) /= 0 .and. box(k) > 0.0d0) d(k) = d(k) - box(k)*anint(d(k)/box(k))
        end do
        hpos(:, idx, r) = d
        bmin = min(bmin, d)
        bmax = max(bmax, d)
      end do
    end do
    if (n_sel > 0) then
      bmid = 0.5d0 * (bmin + bmax)
      do r = 1, n_images
        do idx = 1, n_sel
          hpos(:, idx, r) = hpos(:, idx, r) - bmid
        end do
      end do
      hbox = (bmax - bmin) + 2.0d0 * HBOX_PAD
    else
      hbox = 2.0d0 * HBOX_PAD
    end if
    horig = -0.5d0 * hbox

    do r = 1, n_images
      meta = " neb_id="//trim(istr(ctx%neb_saved))//" image="//trim(istr(r))//" n_images="//trim(istr(n_images))// &
        " energy="//trim(fstr(ctx%pe_buf(r), "(f24.8)"))// &
        " rel_energy="//trim(fstr(ctx%pe_buf(r) - ctx%pe_buf(1), "(f24.8)"))// &
        " barrier="//trim(fstr(barrier, "(f24.8)"))// &
        " reverse_barrier="//trim(fstr(reverse_barrier, "(f24.8)"))// &
        " kmc_step="//trim(istr(kmc_step))//" kmc_time="//trim(fstr(kmc_time, "(es24.10)"))// &
        " jumping_atom="//trim(istr(jumping_atom_i))//" middle_atom="//trim(istr(middle_atom_i))// &
        " site_id="//trim(istr(site_id))//" final_id="//trim(istr(final_id))// &
        " converged="//trim(istr(converged))//" neb_steps="//trim(istr(steps))

      header = 'Lattice="'//trim(fstr(ctx%cell_3x3(1), "(f20.6)"))
      do k = 2, 9
        header = trim(header)//" "//trim(fstr(ctx%cell_3x3(k), "(f20.6)"))
      end do
      header = trim(header)//'" Properties=species:S:1:pos:R:3:id:I:1 pbc="'// &
        trim(merge("T", "F", ctx%pbc_3(1) /= 0))//" "//trim(merge("T", "F", ctx%pbc_3(2) /= 0))//" "// &
        trim(merge("T", "F", ctx%pbc_3(3) /= 0))//'"'//trim(meta)

      hheader = 'Lattice="'//trim(fstr(hbox(1), "(f12.4)"))//" 0.0 0.0 0.0 "// &
        trim(fstr(hbox(2), "(f12.4)"))//" 0.0 0.0 0.0 "//trim(fstr(hbox(3), "(f12.4)"))//'" Origin="'// &
        trim(fstr(horig(1), "(f12.4)"))//" "//trim(fstr(horig(2), "(f12.4)"))//" "// &
        trim(fstr(horig(3), "(f12.4)"))//'" Properties=species:S:1:pos:R:3:id:I:1 pbc="F F F"'//trim(meta)

      write(ctx%neb_unit_active, "(i0)") n_sel
      write(ctx%neb_unit_active, "(a)") trim(header)
      write(ctx%neb_unit_highlight, "(i0)") n_sel
      write(ctx%neb_unit_highlight, "(a)") trim(hheader)
      do idx = 1, n_sel
        i = sel_ids(idx)
        sym = adjustl(ElementName(z(i)))
        write(ctx%neb_unit_active, "(a,3(1x,f11.5),1x,i0)") trim(sym), ctx%img_buf(:, (r-1)*n + i), i
        sym = adjustl(ElementName(z_hl(i)))
        write(ctx%neb_unit_highlight, "(a,3(1x,f11.5),1x,i0)") trim(sym), hpos(:, idx, r), i
      end do
    end do
    flush(ctx%neb_unit_active)
    flush(ctx%neb_unit_highlight)
  end subroutine save_neb_record

  function fstr(x, fmt) result(s)
    double precision, intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=32) :: s
    write(s, fmt) x
    s = adjustl(s)
  end function fstr

  function istr(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    write(s, "(i0)") i
  end function istr

  subroutine swiftmd_destroy_ctx(ctx)
    type(t_swiftmd_ctx), intent(inout) :: ctx

    if (c_associated(ctx%handle)) then
      call c_swiftmd_destroy(ctx%handle)
      ctx%handle = c_null_ptr
    end if
    ctx%enabled = .false.
    if (ctx%neb_unit_full >= 0) then
      close(ctx%neb_unit_full)
      ctx%neb_unit_full = -1
    end if
    if (ctx%neb_unit_active >= 0) then
      close(ctx%neb_unit_active)
      ctx%neb_unit_active = -1
    end if
    if (ctx%neb_unit_highlight >= 0) then
      close(ctx%neb_unit_highlight)
      ctx%neb_unit_highlight = -1
    end if
    ctx%save_neb = .false.
    ctx%save_neb_full = .false.
    if (allocated(ctx%img_buf)) deallocate(ctx%img_buf)
    if (allocated(ctx%species_ptrs)) deallocate(ctx%species_ptrs)
    if (allocated(ctx%species_chars)) deallocate(ctx%species_chars)
    if (allocated(ctx%is_fixed)) deallocate(ctx%is_fixed)
    if (allocated(ctx%pos_init)) deallocate(ctx%pos_init)
    if (allocated(ctx%pos_final)) deallocate(ctx%pos_final)
    if (allocated(ctx%pe_buf)) deallocate(ctx%pe_buf)
    if (allocated(ctx%batch_init)) deallocate(ctx%batch_init)
    if (allocated(ctx%batch_final)) deallocate(ctx%batch_final)
    ctx%n_atoms = 0
    ctx%capacity = 0
    ctx%substrate_pos_loaded = .false.
    ctx%sub_mask_dirty = .false.
  end subroutine swiftmd_destroy_ctx

  subroutine restore_substrate_mask(ctx, system, error)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    type(t_system), intent(in) :: system
    integer, intent(out) :: error

    integer(c_int) :: rc

    error = 0
    if (.not. ctx%sub_mask_dirty) return
    ctx%is_fixed(system%substrate%first_id:system%substrate%last_id) = 0_c_int
    rc = c_swiftmd_update_atoms(ctx%handle, ctx%species_ptrs, ctx%is_fixed, int(ctx%n_atoms, c_int))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_update_atoms (restore) failed")
      error = 1
      return
    end if
    ctx%sub_mask_dirty = .false.
  end subroutine restore_substrate_mask

  subroutine apply_minimize_active_mask(ctx, system, center_pos, error)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    type(t_system), intent(in) :: system
    real(kind=rk), intent(in) :: center_pos(3)
    integer, intent(out) :: error

    integer :: i, s_start, s_end
    real(c_double) :: r2, dx, dy, dz, lx, ly, lz, cx, cy, cz
    integer(c_int) :: rc

    error = 0
    if (ctx%active_region <= 0.0d0) return

    s_start = system%substrate%first_id
    s_end = system%substrate%last_id
    if (s_end < s_start) return

    r2 = real(ctx%active_region * ctx%active_region, c_double)
    cx = real(center_pos(1), c_double)
    cy = real(center_pos(2), c_double)
    cz = real(center_pos(3), c_double)
    lx = ctx%cell_3x3(1)
    ly = ctx%cell_3x3(5)
    lz = ctx%cell_3x3(9)

    do i = s_start, s_end
      dx = ctx%pos_init(1, i) - cx
      dy = ctx%pos_init(2, i) - cy
      dz = ctx%pos_init(3, i) - cz
      if (ctx%pbc_3(1) /= 0_c_int .and. lx > 0.0_c_double) dx = dx - lx * nint(dx / lx)
      if (ctx%pbc_3(2) /= 0_c_int .and. ly > 0.0_c_double) dy = dy - ly * nint(dy / ly)
      if (ctx%pbc_3(3) /= 0_c_int .and. lz > 0.0_c_double) dz = dz - lz * nint(dz / lz)
      if (dx*dx + dy*dy + dz*dz > r2) then
        ctx%is_fixed(i) = 1_c_int
      else
        ctx%is_fixed(i) = 0_c_int
      end if
    end do

    rc = c_swiftmd_update_atoms(ctx%handle, ctx%species_ptrs, ctx%is_fixed, int(ctx%n_atoms, c_int))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_update_atoms (active mask) failed")
      error = 1
      return
    end if
    ctx%sub_mask_dirty = .true.
  end subroutine apply_minimize_active_mask

  subroutine swiftmd_sync_atoms(ctx, system, atoms_cfg, error)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: atoms_cfg
    integer, intent(out) :: error

    integer :: i, n_tot, first_id, last_id
    integer(c_int) :: rc

    error = 0
    n_tot = system%n_substrate + system%n_deposited
    if (n_tot == ctx%n_atoms) return

    call ensure_capacity(ctx, n_tot)
    if (.not. ctx%substrate_pos_loaded) then
      call fill_substrate_metadata(ctx, system)
    end if

    first_id = atoms_cfg%first_id
    last_id = atoms_cfg%last_id
    do i = first_id, last_id
      call set_c_string(ElementName(atoms_cfg%Z(i)), ctx%species_chars(:, i))
      ctx%is_fixed(i) = 0_c_int
    end do
    ctx%n_atoms = n_tot

    rc = c_swiftmd_update_atoms(ctx%handle, ctx%species_ptrs, ctx%is_fixed, int(n_tot, c_int))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_update_atoms failed")
      error = 1
    end if
  end subroutine swiftmd_sync_atoms

  subroutine swiftmd_minimize_cfg(ctx, system, final_config, center_pos, etol, ftol, dmax, max_steps, &
                                  relaxed_pos, error)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: final_config
    real(kind=rk), intent(in) :: center_pos(3)
    double precision, intent(in) :: etol, ftol, dmax
    integer, intent(in) :: max_steps
    real(kind=rk), intent(out) :: relaxed_pos(:,:)
    integer, intent(out) :: error

    integer :: first_id, last_id, n_tot
    integer(c_int) :: rc
    real(c_double), target :: e_out, fmax_out
    integer(c_int), target :: steps_out, conv_out

    call swiftmd_sync_atoms(ctx, system, final_config, error)
    if (error /= 0) return

    call apply_minimize_active_mask(ctx, system, center_pos, error)
    if (error /= 0) return

    n_tot = ctx%n_atoms
    if (ctx%n_substrate > 0) then
      ctx%pos_final(:, 1:ctx%n_substrate) = ctx%pos_init(:, 1:ctx%n_substrate)
    end if
    first_id = final_config%first_id
    last_id = final_config%last_id
    if (last_id >= first_id) then
      ctx%pos_final(:, first_id:last_id) = real(final_config%pos(:, first_id:last_id), c_double)
    end if

    rc = c_swiftmd_minimize(ctx%handle, ctx%pos_final, &
                            real(etol, c_double), real(ftol, c_double), real(dmax, c_double), &
                            int(max_steps, c_int), c_loc(e_out), c_loc(fmax_out), &
                            c_loc(steps_out), c_loc(conv_out))
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_minimize failed")
      error = 1
      return
    end if

    if (last_id >= first_id) then
      relaxed_pos(:, 1:(last_id - first_id + 1)) = real(ctx%pos_final(:, first_id:last_id), rk)
    end if
  end subroutine swiftmd_minimize_cfg

  subroutine swiftmd_neb_cfg(ctx, system, init_cfg, jumping_atom_i, diff_init_target, &
                             middle_atom_i, diff_init_middle, diff_middle_target, &
                             n_images, spring_k, climb, etol, ftol, dmax, max_steps, &
                             barrier, reverse_barrier, error, kmc_step, kmc_time, site_id, final_id, &
                             fmax_res, converged_res, jump_disp_res)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    integer, optional, intent(in) :: kmc_step, site_id, final_id
    double precision, optional, intent(out) :: fmax_res
    double precision, optional, intent(out) :: jump_disp_res  ! jumping atom: first vs last image (A)
    integer, optional, intent(out) :: converged_res
    double precision, optional, intent(in) :: kmc_time
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: init_cfg
    integer, intent(in) :: jumping_atom_i, middle_atom_i
    real(kind=rk), intent(in) :: diff_init_target(3), diff_init_middle(3), diff_middle_target(3)
    integer, intent(in) :: n_images, climb, max_steps
    double precision, intent(in) :: spring_k, etol, ftol, dmax
    double precision, intent(out) :: barrier, reverse_barrier
    integer, intent(out) :: error

    integer :: first_id, last_id
    integer(c_int) :: rc
    real(c_double), target :: fmax_out
    integer(c_int), target :: steps_out, conv_out
    integer(int64) :: clk0, clk1, clk_rate
    double precision :: jd_dbg

    call system_clock(clk0, clk_rate)
    call swiftmd_sync_atoms(ctx, system, init_cfg, error)
    if (error /= 0) return

    call restore_substrate_mask(ctx, system, error)
    if (error /= 0) return

    if (ctx%n_substrate > 0) then
      ctx%pos_final(:, 1:ctx%n_substrate) = ctx%pos_init(:, 1:ctx%n_substrate)
    end if
    first_id = init_cfg%first_id
    last_id = init_cfg%last_id
    if (last_id >= first_id) then
      ctx%pos_init(:, first_id:last_id) = real(init_cfg%pos(:, first_id:last_id), c_double)
      ctx%pos_final(:, first_id:last_id) = ctx%pos_init(:, first_id:last_id)
    end if

    if (middle_atom_i == 0) then
      ctx%pos_final(:, jumping_atom_i) = ctx%pos_init(:, jumping_atom_i) + real(diff_init_target, c_double)
    else
      ctx%pos_final(:, jumping_atom_i) = ctx%pos_init(:, jumping_atom_i) + real(diff_init_middle, c_double)
      ctx%pos_final(:, middle_atom_i) = ctx%pos_init(:, middle_atom_i) + real(diff_middle_target, c_double)
    end if

    if (.not. allocated(ctx%pe_buf)) then
      allocate(ctx%pe_buf(max(64, n_images)))
    else if (size(ctx%pe_buf) < n_images) then
      deallocate(ctx%pe_buf)
      allocate(ctx%pe_buf(max(n_images * 2, 64)))
    end if

    if (ctx%save_neb) then
      if (allocated(ctx%img_buf)) then
        if (size(ctx%img_buf, 2) < ctx%n_atoms*n_images) deallocate(ctx%img_buf)
      end if
      if (.not. allocated(ctx%img_buf)) allocate(ctx%img_buf(3, (ctx%n_atoms + 256)*n_images))
      rc = c_swiftmd_neb_images(ctx%handle, ctx%pos_init, ctx%pos_final, &
                         int(n_images, c_int), real(spring_k, c_double), int(climb, c_int), &
                         real(etol, c_double), real(ftol, c_double), real(dmax, c_double), &
                         int(max_steps, c_int), c_loc(ctx%pe_buf(1)), c_null_ptr, c_null_ptr, &
                         c_loc(fmax_out), c_loc(steps_out), c_loc(conv_out), c_loc(ctx%img_buf(1, 1)))
    else
      ! Only endpoints needed for the "relaxed back to start" displacement check (avoids copying all R images).
      if (allocated(ctx%img_buf)) then
        if (size(ctx%img_buf, 2) < ctx%n_atoms*2) deallocate(ctx%img_buf)
      end if
      if (.not. allocated(ctx%img_buf)) allocate(ctx%img_buf(3, (ctx%n_atoms + 256)*2))
      rc = c_swiftmd_neb(ctx%handle, ctx%pos_init, ctx%pos_final, &
                         int(n_images, c_int), real(spring_k, c_double), int(climb, c_int), &
                         real(etol, c_double), real(ftol, c_double), real(dmax, c_double), &
                         int(max_steps, c_int), c_loc(ctx%pe_buf(1)), &
                         c_loc(ctx%img_buf(1, 1)), c_loc(ctx%img_buf(1, ctx%n_atoms + 1)), &
                         c_loc(fmax_out), c_loc(steps_out), c_loc(conv_out))
    end if
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_neb failed")
      error = 1
      return
    end if

    barrier = maxval(ctx%pe_buf(1:n_images)) - ctx%pe_buf(1)
    reverse_barrier = maxval(ctx%pe_buf(1:n_images)) - ctx%pe_buf(n_images)
    if (present(fmax_res)) fmax_res = real(fmax_out, kind(1.0d0))
    if (present(converged_res)) converged_res = int(conv_out)
    if (present(jump_disp_res)) then
      block
        double precision :: d(3), box(3)
        integer :: k, final_off
        box = [ctx%cell_3x3(1), ctx%cell_3x3(5), ctx%cell_3x3(9)]
        if (ctx%save_neb) then
          final_off = (n_images - 1)*ctx%n_atoms
        else
          final_off = ctx%n_atoms
        end if
        d = ctx%img_buf(:, final_off + jumping_atom_i) - ctx%img_buf(:, jumping_atom_i)
        do k = 1, 3
          if (ctx%pbc_3(k) /= 0 .and. box(k) > 0.0d0) d(k) = d(k) - box(k)*anint(d(k)/box(k))
        end do
        jump_disp_res = sqrt(sum(d*d))
      end block
    end if

    if (ctx%save_neb) then
      block
        integer :: k_step, k_site, k_final
        double precision :: k_time
        k_step = 0; k_site = 0; k_final = 0; k_time = 0.0d0
        if (present(kmc_step)) k_step = kmc_step
        if (present(site_id)) k_site = site_id
        if (present(final_id)) k_final = final_id
        if (present(kmc_time)) k_time = kmc_time
        call save_neb_record(ctx, system, init_cfg, jumping_atom_i, middle_atom_i, n_images, &
                             barrier, reverse_barrier, real(fmax_out, kind(1.0d0)), int(steps_out), &
                             int(conv_out), k_step, k_time, k_site, k_final)
      end block
    end if
    if (ctx%debug) then
      call system_clock(clk1)
      jd_dbg = -1.0d0
      if (present(jump_disp_res)) jd_dbg = jump_disp_res
      write(6, "(a,i0,a,es10.3,a,i0,a,f8.3,a,f9.5,a,f9.5,a,f7.3,a)") "# [SwiftMD NEB] steps=", steps_out, &
        " fmax=", fmax_out, " conv=", conv_out, " wall=", real(clk1 - clk0, kind(1.0d0))/real(clk_rate, kind(1.0d0)), &
        " s barrier=", barrier, " rev=", reverse_barrier, " jump=", jd_dbg, " A"
    end if
  end subroutine swiftmd_neb_cfg

  subroutine swiftmd_neb_batch_cfg(ctx, system, init_cfg, n_bands, jumping_atoms, middle_atoms, &
                                   disp1, disp2, n_images, spring_k, climb, etol, ftol, dmax, &
                                   max_steps, barriers_fwd, barriers_rev, fmax_arr, steps_arr, &
                                   conv_arr, error)
    type(t_swiftmd_ctx), target, intent(inout) :: ctx
    type(t_system), intent(in) :: system
    type(t_atoms), intent(in) :: init_cfg
    integer, intent(in) :: n_bands, n_images, climb, max_steps
    integer, intent(in) :: jumping_atoms(n_bands), middle_atoms(n_bands)
    real(kind=rk), intent(in) :: disp1(3, n_bands), disp2(3, n_bands)
    double precision, intent(in) :: spring_k, etol, ftol, dmax
    double precision, intent(out) :: barriers_fwd(n_bands), barriers_rev(n_bands), fmax_arr(n_bands)
    integer, intent(out) :: steps_arr(n_bands), conv_arr(n_bands)
    integer, intent(out) :: error

    integer :: first_id, last_id, b, cur_band_cap
    integer(c_int) :: rc
    real(c_double), allocatable :: c_bfwd(:), c_brev(:), c_fmax(:)
    integer(c_int), allocatable :: c_steps(:), c_conv(:)

    error = 0
    if (n_bands <= 0) return

    call swiftmd_sync_atoms(ctx, system, init_cfg, error)
    if (error /= 0) return

    call restore_substrate_mask(ctx, system, error)
    if (error /= 0) return

    first_id = init_cfg%first_id
    last_id = init_cfg%last_id
    if (last_id >= first_id) then
      ctx%pos_init(:, first_id:last_id) = real(init_cfg%pos(:, first_id:last_id), c_double)
    end if

    cur_band_cap = 0
    if (allocated(ctx%batch_init)) then
      if (size(ctx%batch_init, 2) < ctx%n_atoms .or. size(ctx%batch_init, 3) < n_bands) then
        deallocate(ctx%batch_init, ctx%batch_final)
      else
        cur_band_cap = size(ctx%batch_init, 3)
      end if
    end if
    if (cur_band_cap == 0) then
      cur_band_cap = max(32, n_bands)
      allocate(ctx%batch_init(3, ctx%n_atoms, cur_band_cap))
      allocate(ctx%batch_final(3, ctx%n_atoms, cur_band_cap))
    end if

    do b = 1, n_bands
      ctx%batch_init(:, 1:ctx%n_atoms, b) = ctx%pos_init(:, 1:ctx%n_atoms)
      ctx%batch_final(:, 1:ctx%n_atoms, b) = ctx%pos_init(:, 1:ctx%n_atoms)
      if (middle_atoms(b) == 0) then
        ctx%batch_final(:, jumping_atoms(b), b) = ctx%pos_init(:, jumping_atoms(b)) + real(disp1(:, b), c_double)
      else
        ctx%batch_final(:, jumping_atoms(b), b) = ctx%pos_init(:, jumping_atoms(b)) + real(disp1(:, b), c_double)
        ctx%batch_final(:, middle_atoms(b), b) = ctx%pos_init(:, middle_atoms(b)) + real(disp2(:, b), c_double)
      end if
    end do

    allocate(c_bfwd(n_bands), c_brev(n_bands), c_fmax(n_bands), c_steps(n_bands), c_conv(n_bands))
    rc = c_swiftmd_neb_batch(ctx%handle, int(n_bands, c_int), ctx%batch_init, ctx%batch_final, &
                             int(n_images, c_int), real(spring_k, c_double), int(climb, c_int), &
                             real(etol, c_double), real(ftol, c_double), real(dmax, c_double), &
                             int(max_steps, c_int), c_bfwd, c_brev, c_null_ptr, &
                             c_fmax, c_steps, c_conv)
    if (rc /= 0_c_int) then
      call print_swiftmd_error("swiftmd_neb_batch failed")
      error = 1
      return
    end if

    barriers_fwd(1:n_bands) = real(c_bfwd(1:n_bands), kind(1.0d0))
    barriers_rev(1:n_bands) = real(c_brev(1:n_bands), kind(1.0d0))
    fmax_arr(1:n_bands) = real(c_fmax(1:n_bands), kind(1.0d0))
    steps_arr(1:n_bands) = int(c_steps(1:n_bands))
    conv_arr(1:n_bands) = int(c_conv(1:n_bands))
    ctx%neb_batch_calls = ctx%neb_batch_calls + 1_int64
    ctx%neb_batch_bands = ctx%neb_batch_bands + int(n_bands, int64)

    if (ctx%debug) then
      write(6, "(a,i0,a,es10.3)") "# [SwiftMD NEB batch] n_bands=", n_bands, &
        " max_fmax=", maxval(fmax_arr(1:n_bands))
    end if
  end subroutine swiftmd_neb_batch_cfg

end module swiftmd_module

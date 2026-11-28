! MLKMC - Fortran iso_c_binding interface to JGAP persistent single-center SOAP C ABI (jgap_c.h)
!
! Replaces QUIP Atoms + calc_connect_atom + soap_calc with zero-allocation single-center
! evaluation directly from contiguous LAE coordinate/Z buffers.

module jgap_soap_module
  use, intrinsic :: iso_c_binding
  implicit none

  private
  public :: t_jgap_soap_ctx
  public :: jgap_soap_init_ctx
  public :: jgap_soap_eval_buf
  public :: jgap_soap_destroy_ctx

  interface
    function c_jgap_soap_create_from_quip_str(quip_soap_str) &
        bind(C, name="jgap_soap_create_from_quip_str") result(handle)
      import :: c_ptr
      type(c_ptr), value :: quip_soap_str
      type(c_ptr) :: handle
    end function c_jgap_soap_create_from_quip_str

    function c_jgap_soap_dim(handle) bind(C, name="jgap_soap_dim") result(dim)
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: dim
    end function c_jgap_soap_dim

    function c_jgap_soap_eval_lae(handle, cell_3x3, pbc_3, pos_3xn, z_n, n_atoms, &
        center_idx_0based, desc_out) bind(C, name="jgap_soap_eval_lae") result(rc)
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      real(c_double), intent(in) :: cell_3x3(9)
      integer(c_int), intent(in) :: pbc_3(3)
      real(c_double), intent(in) :: pos_3xn(*)
      integer(c_int), intent(in) :: z_n(*)
      integer(c_int), value :: n_atoms
      integer(c_int), value :: center_idx_0based
      real(c_double), intent(out) :: desc_out(*)
      integer(c_int) :: rc
    end function c_jgap_soap_eval_lae

    subroutine c_jgap_soap_destroy(handle) bind(C, name="jgap_soap_destroy")
      import :: c_ptr
      type(c_ptr), value :: handle
    end subroutine c_jgap_soap_destroy

    function c_jgap_last_error() bind(C, name="jgap_last_error") result(err_ptr)
      import :: c_ptr
      type(c_ptr) :: err_ptr
    end function c_jgap_last_error
  end interface

  type t_jgap_soap_ctx
    type(c_ptr) :: handle = c_null_ptr
    logical :: enabled = .false.
    logical :: verify = .false.
    integer :: dim = 0
    real(c_double) :: cell_3x3(9) = 0.0_c_double
    integer(c_int) :: pbc_3(3) = [1_c_int, 1_c_int, 1_c_int]
    integer(8) :: eval_count = 0
    double precision :: max_verify_diff = 0.0d0
  end type t_jgap_soap_ctx

contains

  subroutine print_jgap_error(prefix)
    character(len=*), intent(in) :: prefix
    type(c_ptr) :: p
    character(kind=c_char), pointer :: chars(:)
    integer :: i, n
    character(len=512) :: msg

    p = c_jgap_last_error()
    if (.not. c_associated(p)) then
      write(0, "(a)") "*** "//trim(prefix)//": (no error message)"
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
    write(0, "(a)") "*** "//trim(prefix)//": "//trim(msg(1:n))
  end subroutine print_jgap_error

  subroutine jgap_soap_init_ctx(ctx, desc_str, lattice, error)
    type(t_jgap_soap_ctx), intent(inout) :: ctx
    character(len=*), intent(in) :: desc_str
    double precision, intent(in) :: lattice(3, 3)
    integer, intent(out) :: error

    character(len=32) :: env_val
    integer :: env_len, env_stat, i, n
    character(kind=c_char, len=1), allocatable, target :: c_str(:)

    error = 0
    ctx%enabled = .false.
    ctx%verify = .false.
    ctx%dim = 0
    ctx%eval_count = 0
    ctx%max_verify_diff = 0.0d0

    call get_environment_variable("MLKMC_JGAP_SOAP", env_val, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      env_val = adjustl(env_val)
      if (trim(env_val) == "0" .or. trim(env_val) == "F" .or. trim(env_val) == "f" .or. &
          trim(env_val) == "false" .or. trim(env_val) == "FALSE" .or. trim(env_val) == "off") then
        return
      end if
    end if

    call get_environment_variable("MLKMC_JGAP_SOAP_VERIFY", env_val, length=env_len, status=env_stat)
    if (env_stat == 0 .and. env_len > 0) then
      env_val = adjustl(env_val)
      if (trim(env_val) == "1" .or. trim(env_val) == "T" .or. trim(env_val) == "t" .or. &
          trim(env_val) == "true" .or. trim(env_val) == "TRUE") then
        ctx%verify = .true.
      end if
    end if

    n = len_trim(desc_str)
    allocate(c_str(n + 1))
    do i = 1, n
      c_str(i) = desc_str(i:i)
    end do
    c_str(n + 1) = c_null_char

    ctx%handle = c_jgap_soap_create_from_quip_str(c_loc(c_str(1)))
    deallocate(c_str)
    if (.not. c_associated(ctx%handle)) then
      call print_jgap_error("jgap_soap_create_from_quip_str failed")
      error = 1
      return
    end if

    ctx%dim = int(c_jgap_soap_dim(ctx%handle))
    if (ctx%dim <= 0) then
      call print_jgap_error("jgap_soap_dim returned <= 0")
      call c_jgap_soap_destroy(ctx%handle)
      ctx%handle = c_null_ptr
      error = 1
      return
    end if

    ctx%cell_3x3 = reshape(real(lattice(1:3, 1:3), c_double), [9])
    ! Match QUIP initialise(lae, 0, lattice), which uses 3D minimum-image PBC inside the cell.
    ctx%pbc_3 = [1_c_int, 1_c_int, 1_c_int]
    ctx%enabled = .true.
  end subroutine jgap_soap_init_ctx

  subroutine jgap_soap_eval_buf(ctx, pos_3xn, z_n, n_atoms, center_idx_0based, desc_out, error)
    type(t_jgap_soap_ctx), intent(inout) :: ctx
    real(c_double), intent(in) :: pos_3xn(3, *)
    integer(c_int), intent(in) :: z_n(*)
    integer, intent(in) :: n_atoms, center_idx_0based
    double precision, intent(out) :: desc_out(*)
    integer, intent(out) :: error

    integer(c_int) :: rc

    error = 0
    rc = c_jgap_soap_eval_lae(ctx%handle, ctx%cell_3x3, ctx%pbc_3, &
                              pos_3xn, z_n, int(n_atoms, c_int), &
                              int(center_idx_0based, c_int), desc_out)
    if (rc /= 0_c_int) then
      call print_jgap_error("jgap_soap_eval_lae failed")
      error = int(rc)
      return
    end if
    ctx%eval_count = ctx%eval_count + 1
  end subroutine jgap_soap_eval_buf

  subroutine jgap_soap_destroy_ctx(ctx)
    type(t_jgap_soap_ctx), intent(inout) :: ctx

    if (c_associated(ctx%handle)) then
      call c_jgap_soap_destroy(ctx%handle)
      ctx%handle = c_null_ptr
    end if
    ctx%enabled = .false.
    ctx%verify = .false.
    ctx%dim = 0
  end subroutine jgap_soap_destroy_ctx

end module jgap_soap_module

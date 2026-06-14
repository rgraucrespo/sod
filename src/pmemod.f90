!*******************************************************************************
!    Copyright (c) 2026 Ricardo Grau-Crespo, Said Hamad, Salvador R.G. Balestra
!
!    This file is part of the SOD package.
!
!    SOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SOD.  If not, see <http://www.gnu.org/licenses/>.
!
!*******************************************************************************

module pmemod
  use iso_fortran_env, only: real64, int64, error_unit
  use structwriters
  use insod_reader,    only: insod_t, read_insod
  use ensemble_io,     only: read_ensemble, write_ensemble, read_energies_file
  implicit none
  private

  integer, parameter :: max_motif_order = 4
  integer, parameter :: n_calib_configs = 9
  real(real64), parameter :: assign_tol = 1.0e-8_real64

  integer, save :: nop = 0
  integer, save :: npos = 0
  integer, save :: target_species = 0
  integer, save :: target_atini = 0
  integer, save :: target_atfin = 0
  integer, save :: max_low_order = 0
  integer, save :: max_high_order = 0
  logical, save :: high_base_loaded = .false.
  logical, save :: model_initialized = .false.
  real(real64), save :: v0_low = 0.0_real64
  real(real64), save :: v0_high = 0.0_real64
  real(real64), save :: eps_low(0:4)  = [0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]
  real(real64), save :: eps_high(0:4) = [0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]
  real(real64), save :: alpha_hybrid = 2.0_real64   ! power-law exponent for PMEh weighting
  real(real64), save :: eta_hybrid   = 0.0_real64   ! log-scale asymmetry; >0 favours low-end
  integer, save :: pme_choice = 2            ! 0=PME0, 1=PME1, 2=PMEh (default)
  integer, save :: pme_order_cap = max_motif_order  ! user-specified cap on expansion order
  integer, save :: pme_target_level = 0            ! set during initialize; guards reference loading
  logical, save :: pme_calibration_applied = .false.
  integer, save :: calib_config_list(9) = 0  ! config indices from INPME calib_config_list
  integer, save :: n_calib_list = 0          ! how many entries are populated
  integer,      save :: n_calib_used   = 0           ! configs actually used in last fit
  real(real64), save :: rms_calib_low  = 0.0_real64  ! RMS of last low-side fit
  real(real64), save :: rms_calib_high = 0.0_real64  ! RMS of last high-side fit

  integer, allocatable, save :: eqmatrix(:, :)
  real(real64), allocatable, save :: V1_low(:), V2_low(:, :)
  integer, allocatable, save :: V3_i_low(:), V3_j_low(:), V3_k_low(:)
  real(real64), allocatable, save :: V3_val_low(:)
  integer(int64), allocatable, save :: V3_key_low(:)
  integer, allocatable, save :: V4_i_low(:), V4_j_low(:), V4_k_low(:), V4_l_low(:)
  real(real64), allocatable, save :: V4_val_low(:)
  integer(int64), allocatable, save :: V4_key_low(:)
  real(real64), allocatable, save :: V1_high(:), V2_high(:, :)
  integer, allocatable, save :: V3_i_high(:), V3_j_high(:), V3_k_high(:)
  real(real64), allocatable, save :: V3_val_high(:)
  integer(int64), allocatable, save :: V3_key_high(:)
  integer, allocatable, save :: V4_i_high(:), V4_j_high(:), V4_k_high(:), V4_l_high(:)
  real(real64), allocatable, save :: V4_val_high(:)
  integer(int64), allocatable, save :: V4_key_high(:)
  logical, allocatable, save :: work_is_ge(:)
  integer, allocatable, save :: work_holes(:)

  ! Compact V residual cache — one value per inequivalent reference config per level.
  ! Populated during V fitting (load_*_side_model) or loaded from INPME V block.
  integer, save :: nv1_low  = 0, nv2_low  = 0, nv3_low  = 0, nv4_low  = 0
  integer, save :: nv1_high = 0, nv2_high = 0, nv3_high = 0, nv4_high = 0
  real(real64), allocatable, save :: v1res_low(:),  v2res_low(:)
  real(real64), allocatable, save :: v3res_low(:),  v4res_low(:)
  real(real64), allocatable, save :: v1res_high(:), v2res_high(:)
  real(real64), allocatable, save :: v3res_high(:), v4res_high(:)
  logical, save :: v_loaded_from_cache = .false.
  character(len=256), save :: pme_model_filename = 'pme.model'

  ! Set .false. when the model is initialised geometry-only (uniform sampling
  ! without reference ENERGIES). The V expansion is not built, so energies must
  ! not be evaluated.
  logical, save :: pme_energies_avail = .true.

  ! INSOD data cache — populated once by build_target_layout_from_inputs
  logical, save :: insod_cache_valid = .false.
  integer, save :: insod_filer_cache = -1
  integer, save :: insod_nsp_cache = 0
  character(len=4), allocatable, save :: insod_symbols_cache(:)
  character(len=5), save :: insod_newsym_new_cache = '', insod_newsym_old_cache = ''
  logical, save :: insod_has_special_cache = .false.

  public :: pme_initialize_model, pme_finalize_model, pme_get_target_level_from_insod
  public :: pme_write_level_outputs, pme_print_model_summary, pme_evaluate_configuration
  public :: pme_get_npos, pme_get_max_low_order, pme_get_max_high_order
  public :: pme_has_high_side
  public :: pme_energies_available
  public :: pme_preload_recalibration
  public :: pme_get_target_atini, pme_get_v0
  public :: pme_get_epsilon
  public :: pme_get_alpha_eta
  public :: pme_get_choice
  public :: canonicalize_config
  public :: apply_epsilon_energy
  public :: apply_epsilon_energy_decomp
  public :: mc_insertion_sort
  public :: compute_config_degeneracy
  public :: format_level_directory, write_pme_model_template, write_pme_model_calibrated
  public :: pme_set_model_filename
  public :: pme_variant_dir_from_model
  public :: pme_get_newsymbol

contains

  subroutine pme_initialize_model(target_level, energy_optional)
    !  Build the PME effective Hamiltonian for the given target level.
    !  If energy_optional is present and .true. (uniform sampling) and the
    !  low-side reference ENERGIES are absent, initialise geometry only:
    !  the V expansion is skipped and pme_energies_avail is set .false., so
    !  the caller can sample configurations without evaluating energies.
    integer, intent(in) :: target_level
    logical, intent(in), optional :: energy_optional
    logical :: v_ok, recon_ok, opt

    opt = .false.
    if (present(energy_optional)) opt = energy_optional

    call pme_finalize_model()
    pme_target_level = target_level
    call read_eqmatrix_file('EQMATRIX')
    call build_target_layout_from_inputs(target_species)
    if (npos <= 0) then
      write(error_unit,'(A)') 'Error: invalid EQMATRIX header while initialising PME Hamiltonian.'
      stop 1
    end if
    if (target_level < 0 .or. target_level > npos) then
      write(error_unit,'(A,I0,A,I0,A)') 'Error: requested level ', target_level, ' is outside 0..', npos, '.'
      stop 1
    end if

    ! Geometry-only mode: no Hamiltonian to build, no energies to evaluate.
    if (opt .and. .not. low_side_refs_present()) then
      pme_energies_avail = .false.
      model_initialized  = .true.
      return
    end if
    pme_energies_avail = .true.

    call read_pme_model_v_block(trim(pme_model_filename), v_ok)

    if (v_ok .and. nv1_low > 0) then
      call reconstruct_v_low_from_cache(recon_ok)
      if (.not. recon_ok) then
        call load_low_side_model()
      else
        v_loaded_from_cache = .true.
      end if
    else
      call load_low_side_model()
    end if

    if (v_ok .and. nv1_high > 0) then
      call reconstruct_v_high_from_cache(recon_ok)
      if (.not. recon_ok) call load_high_side_model()
    else
      call load_high_side_model()
    end if

    if (v_loaded_from_cache) then
      write(*,'(A,A,A)') ' > V parameters loaded from ', trim(pme_model_filename), &
        ' energy-terms cache; reference ENERGIES not re-read.'
    end if

    call note_target_energies_if_present(target_level)

    model_initialized = .true.
  end subroutine pme_initialize_model

  subroutine pme_finalize_model()
    if (allocated(eqmatrix)) deallocate(eqmatrix)
    if (allocated(V1_low)) deallocate(V1_low)
    if (allocated(V2_low)) deallocate(V2_low)
    if (allocated(V3_i_low)) deallocate(V3_i_low, V3_j_low, V3_k_low, V3_val_low)
    if (allocated(V3_key_low)) deallocate(V3_key_low)
    if (allocated(V4_i_low)) deallocate(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low)
    if (allocated(V4_key_low)) deallocate(V4_key_low)
    if (allocated(V1_high)) deallocate(V1_high)
    if (allocated(V2_high)) deallocate(V2_high)
    if (allocated(V3_i_high)) deallocate(V3_i_high, V3_j_high, V3_k_high, V3_val_high)
    if (allocated(V3_key_high)) deallocate(V3_key_high)
    if (allocated(V4_i_high)) deallocate(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high)
    if (allocated(V4_key_high)) deallocate(V4_key_high)
    if (allocated(work_is_ge)) deallocate(work_is_ge)
    if (allocated(work_holes)) deallocate(work_holes)
    if (allocated(v1res_low))  deallocate(v1res_low)
    if (allocated(v2res_low))  deallocate(v2res_low)
    if (allocated(v3res_low))  deallocate(v3res_low)
    if (allocated(v4res_low))  deallocate(v4res_low)
    if (allocated(v1res_high)) deallocate(v1res_high)
    if (allocated(v2res_high)) deallocate(v2res_high)
    if (allocated(v3res_high)) deallocate(v3res_high)
    if (allocated(v4res_high)) deallocate(v4res_high)
    nv1_low  = 0;  nv2_low  = 0;  nv3_low  = 0;  nv4_low  = 0
    nv1_high = 0;  nv2_high = 0;  nv3_high = 0;  nv4_high = 0
    v_loaded_from_cache = .false.
    pme_energies_avail = .true.
    nop = 0
    npos = 0
    target_species = 0
    target_atini = 0
    target_atfin = 0
    max_low_order = 0
    max_high_order = 0
    high_base_loaded = .false.
    model_initialized = .false.
    v0_low = 0.0_real64
    v0_high = 0.0_real64
    eps_low(0)   = 0.0_real64; eps_low(1:4)  = 1.0_real64
    eps_high(0)  = 0.0_real64; eps_high(1:4) = 1.0_real64
    alpha_hybrid = 2.0_real64
    eta_hybrid   = 0.0_real64
    pme_choice = 2
    pme_order_cap = max_motif_order
    pme_target_level = 0
    n_calib_used   = 0
    rms_calib_low  = 0.0_real64
    rms_calib_high = 0.0_real64
  end subroutine pme_finalize_model

  subroutine pme_set_model_filename(filename)
    character(len=*), intent(in) :: filename
    pme_model_filename = trim(filename)
  end subroutine pme_set_model_filename

  integer function pme_get_npos() result(value)
    value = npos
  end function pme_get_npos

  integer function pme_get_max_low_order() result(value)
    value = max_low_order
  end function pme_get_max_low_order

  integer function pme_get_max_high_order() result(value)
    value = max_high_order
  end function pme_get_max_high_order

  logical function pme_is_initialized() result(value)
    value = model_initialized
  end function pme_is_initialized

  logical function pme_has_high_side() result(value)
    value = high_base_loaded
  end function pme_has_high_side

  logical function pme_energies_available() result(value)
    !  .false. when the model was initialised geometry-only (uniform sampling
    !  without reference ENERGIES); callers must not evaluate configuration
    !  energies in that case.
    value = pme_energies_avail
  end function pme_energies_available

  logical function low_side_refs_present() result(present_all)
    !  The minimal references required to build the low-side expansion are
    !  n00/ENERGIES, n01/ENERGIES and n02/ENERGIES (see load_low_side_model).
    character(len=32) :: dir
    logical :: ex
    integer :: lvl
    present_all = .true.
    do lvl = 0, 2
      call format_level_directory(lvl, dir)
      inquire(file=trim(dir)//'/ENERGIES', exist=ex)
      if (.not. ex) then
        present_all = .false.
        return
      end if
    end do
  end function low_side_refs_present

  integer function pme_get_choice() result(value)
    value = pme_choice
  end function pme_get_choice

  integer function pme_get_target_species() result(value)
    value = target_species
  end function pme_get_target_species

  integer function pme_get_nop() result(value)
    value = nop
  end function pme_get_nop

  integer function pme_get_target_atini() result(value)
    value = target_atini
  end function pme_get_target_atini

  subroutine pme_get_v0(v0_low_out, v0_high_out)
    real(real64), intent(out) :: v0_low_out, v0_high_out
    v0_low_out  = v0_low
    v0_high_out = v0_high
  end subroutine pme_get_v0

  subroutine pme_get_newsymbol(orig_sym_out, sym_new_out, sym_rem_out)
    character(len=*), intent(out) :: orig_sym_out, sym_new_out, sym_rem_out
    integer :: sp
    sp = target_species
    orig_sym_out = ''
    if (allocated(insod_symbols_cache) .and. sp >= 1 .and. sp <= size(insod_symbols_cache)) &
      orig_sym_out = trim(insod_symbols_cache(sp))
    sym_new_out = trim(insod_newsym_new_cache)
    sym_rem_out = trim(insod_newsym_old_cache)
  end subroutine pme_get_newsymbol

  subroutine pme_get_epsilon(eps_low_out, eps_high_out)
    !  Return the per-order renormalisation parameters epsilon (PME §3.2).
    real(real64), intent(out) :: eps_low_out(0:4), eps_high_out(0:4)
    eps_low_out  = eps_low(0:4)
    eps_high_out = eps_high(0:4)
  end subroutine pme_get_epsilon

  subroutine pme_get_alpha_eta(alpha_out, eta_out)
    real(real64), intent(out) :: alpha_out, eta_out
    alpha_out = alpha_hybrid
    eta_out   = eta_hybrid
  end subroutine pme_get_alpha_eta

  subroutine pme_print_model_summary()
    integer :: k
    character(len=32) :: tag

    write(*,*)
    write(*,'(A)') '  > Effective Hamiltonian summary:'
    write(*,'(A,I6)')   '    Target sites   (npos) : ', npos
    write(*,'(A,I6)')   '    Symmetry ops   (nop)  : ', nop
    write(*,*)
    if (.not. pme_energies_avail) then
      write(*,'(A)') '    Geometry-only: no PME Hamiltonian built (uniform sampling'
      write(*,'(A)') '    without reference ENERGIES). Energies will not be evaluated.'
      write(*,*)
      return
    end if
    write(*,'(A)')      '    LOW-SIDE expansion  (x → 0):'
    write(*,'(A,F22.10,A)') '      V₀_low               : ', v0_low, ' eV'
    write(*,'(A,I2)')   '      Maximum motif order: ', max_low_order
    if (v_loaded_from_cache) then
      write(*,'(A,A,A)') '      V parameters: loaded from ', trim(pme_model_filename), ' cache'
    else
      write(*,'(A)', advance='no') '      Training levels used: '
      call format_level_directory(0, tag)
      write(*,'(A)', advance='no') trim(tag)
      do k = 1, max_low_order
        call format_level_directory(k, tag)
        write(*,'(A,A)', advance='no') ', ', trim(tag)
      end do
      write(*,*)
    end if
    write(*,*)
    if (high_base_loaded) then
      write(*,'(A)')    '    HIGH-SIDE expansion (x → 1):'
      write(*,'(A,F22.10,A)') '      V₀_high              : ', v0_high, ' eV'
      write(*,'(A,I2)') '      Maximum motif order: ', max_high_order
      if (v_loaded_from_cache) then
        write(*,'(A,A,A)') '      V parameters: loaded from ', trim(pme_model_filename), ' cache'
      else
        write(*,'(A)', advance='no') '      Training levels used: '
        call format_level_directory(npos, tag)
        write(*,'(A)', advance='no') trim(tag)
        do k = 1, max_high_order
          call format_level_directory(npos - k, tag)
          write(*,'(A,A)', advance='no') ', ', trim(tag)
        end do
        write(*,*)
      end if
    else
      write(*,'(A)') '    HIGH-SIDE expansion: not available (no n{npos}/ENERGIES found).'
    end if
    write(*,*)
  end subroutine pme_print_model_summary

  subroutine pme_get_blend_weights(level, weight_low, weight_high)
    integer, intent(in) :: level
    real(real64), intent(out) :: weight_low, weight_high

    call require_model_initialized('pme_get_blend_weights')
    call validate_level(level, 'pme_get_blend_weights')
    if (.not. high_base_loaded) then
      weight_low = 1.0_real64
      weight_high = 0.0_real64
      return
    end if
    call level_blend_weights(level, npos, weight_low, weight_high)
  end subroutine pme_get_blend_weights

  subroutine pme_load_level_ensemble(level, mm, conf, omega, ok)
    integer, intent(in) :: level
    integer, intent(out) :: mm
    integer, allocatable, intent(out) :: conf(:, :), omega(:)
    logical, intent(out) :: ok
    character(len=32) :: level_dir
    integer :: level_from_file, total_sites

    call require_model_initialized('pme_load_level_ensemble')
    call validate_level(level, 'pme_load_level_ensemble')
    call format_level_directory(level, level_dir)
    call read_ensemble_level(trim(level_dir)//'/ENSEMBLE', level_from_file, total_sites, mm, conf, omega, ok)
    if (.not. ok) return
    if (level_from_file /= level .or. total_sites /= npos) then
      ok = .false.
      if (allocated(conf)) deallocate(conf)
      if (allocated(omega)) deallocate(omega)
    end if
  end subroutine pme_load_level_ensemble

  subroutine pme_load_level_dataset(level, mm, conf, omega, energies, ok)
    integer, intent(in) :: level
    integer, intent(out) :: mm
    integer, allocatable, intent(out) :: conf(:, :), omega(:)
    real(real64), allocatable, intent(out) :: energies(:)
    logical, intent(out) :: ok

    call require_model_initialized('pme_load_level_dataset')
    call validate_level(level, 'pme_load_level_dataset')
    call read_level_data(level, mm, conf, omega, energies, ok)
  end subroutine pme_load_level_dataset

  subroutine pme_build_ordered_configuration(level, conf, ok)
    integer, intent(in) :: level
    integer, allocatable, intent(out) :: conf(:)
    logical, intent(out) :: ok
    integer :: i

    call require_model_initialized('pme_build_ordered_configuration')
    ok = .false.
    if (level < 0 .or. level > npos) return
    allocate(conf(level))
    do i = 1, level
      conf(i) = i
    end do
    ok = .true.
  end subroutine pme_build_ordered_configuration

  logical function pme_configuration_is_valid(conf_row) result(ok)
    integer, intent(in) :: conf_row(:)
    integer :: i, level

    ok = .false.
    if (.not. model_initialized) return
    level = size(conf_row)
    if (level < 0 .or. level > npos) return
    do i = 1, level
      if (conf_row(i) < 1 .or. conf_row(i) > npos) return
      if (i > 1) then
        if (conf_row(i) <= conf_row(i - 1)) return
      end if
    end do
    ok = .true.
  end function pme_configuration_is_valid

  subroutine pme_sort_configuration(conf_row)
    integer, intent(inout) :: conf_row(:)
    integer :: i, j, tmp

    do i = 1, size(conf_row) - 1
      do j = i + 1, size(conf_row)
        if (conf_row(j) < conf_row(i)) then
          tmp = conf_row(i)
          conf_row(i) = conf_row(j)
          conf_row(j) = tmp
        end if
      end do
    end do
  end subroutine pme_sort_configuration

  subroutine pme_evaluate_subset(conf_row, energy, energy_low, energy_high, low_terms, high_terms)
    integer, intent(in) :: conf_row(:)
    real(real64), intent(out) :: energy
    real(real64), intent(out), optional :: energy_low, energy_high
    real(real64), intent(out), optional :: low_terms(4), high_terms(4)
    real(real64) :: local_energy_low, local_energy_high
    real(real64) :: local_low_terms(4), local_high_terms(4)

    call require_model_initialized('pme_evaluate_subset')
    if (.not. pme_configuration_is_valid(conf_row)) then
      write(error_unit,'(A)') 'Error: invalid configuration passed to pme_evaluate_subset.'
      stop 1
    end if

    call pme_evaluate_configuration(conf_row, size(conf_row), energy, local_energy_low, local_energy_high, &
      local_low_terms, local_high_terms)
    if (present(energy_low)) energy_low = local_energy_low
    if (present(energy_high)) energy_high = local_energy_high
    if (present(low_terms)) low_terms = local_low_terms
    if (present(high_terms)) high_terms = local_high_terms
  end subroutine pme_evaluate_subset

  subroutine pme_get_target_level_from_insod(target_level)
    integer, intent(out) :: target_level
    type(insod_t) :: d
    call read_insod('INSOD', d)
    target_level = d%nsubs_max
    if (d%nsubs_min /= d%nsubs_max) &
      write(*,*) ' > INSOD range detected; using upper endpoint as PME target level:', target_level
  end subroutine pme_get_target_level_from_insod


  subroutine read_target_species_from_insod(species_index)
    integer, intent(out) :: species_index
    type(insod_t) :: d
    call read_insod('INSOD', d)
    species_index = d%sptarget(1)
    if (species_index <= 0) stop 'Error: invalid sptarget in INSOD.'
  end subroutine read_target_species_from_insod

  subroutine build_target_layout_from_inputs(species_index)
    integer, intent(out) :: species_index
    integer :: unit_sgo, ios, nsp, nat0, at0, op1, nop1
    integer :: na, nb, nc, sp, cumnatsp, at1r, at1, at1i, nat1r, nat1
    integer :: nat_super
    integer, allocatable :: natsp0(:), natsp1(:), natsp(:), spat0(:), spat1(:), spat1r(:)
    real(real64), allocatable :: coords0(:, :), coords1(:, :), coords1r(:, :)
    real(real64), allocatable :: mgroup1(:, :, :), vgroup1(:, :)
    real(real64) :: prod
    logical :: found
    real(real64), parameter :: tol0 = 1.0e-3_real64
    type(insod_t) :: d

    call read_insod('INSOD', d)
    nsp            = d%nsp
    nat0           = d%nat0
    na = d%na; nb = d%nb; nc = d%nc
    species_index  = d%sptarget(1)

    ! Populate module-level cache for write_pme_model_template
    if (allocated(insod_symbols_cache)) deallocate(insod_symbols_cache)
    allocate(insod_symbols_cache(nsp))
    insod_symbols_cache(1:nsp) = d%symbol
    insod_nsp_cache            = nsp
    insod_newsym_new_cache     = d%newsymbol(1,1)
    insod_newsym_old_cache     = d%newsymbol(1,2)
    insod_has_special_cache    = (insod_newsym_new_cache(1:1) == '@' .or. &
                                  insod_newsym_new_cache(1:1) == '%' .or. &
                                  insod_newsym_old_cache(1:1) == '@' .or. &
                                  insod_newsym_old_cache(1:1) == '%')
    insod_filer_cache          = d%filer
    insod_cache_valid          = .true.

    allocate(natsp0(nsp), natsp1(nsp), natsp(nsp))
    natsp0 = d%natsp0
    allocate(coords0(nat0, 3), spat0(nat0))
    coords0(1:nat0, 1:3) = d%coords0

    open(newunit=unit_sgo, file='SGO', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A)') 'Error: could not open SGO while building PME target layout.'
      stop 1
    end if

    cumnatsp = 0
    do sp = 1, nsp
      do at0 = cumnatsp + 1, cumnatsp + natsp0(sp)
        spat0(at0) = sp
      end do
      cumnatsp = cumnatsp + natsp0(sp)
    end do

    read(unit_sgo,*,iostat=ios)
    if (ios /= 0) stop 'Error: malformed SGO header.'
    read(unit_sgo,*,iostat=ios) op1
    if (ios /= 0) stop 'Error: malformed SGO operator header.'
    nop1 = 0
    do while (op1 > 0)
      nop1 = max(nop1, op1)
      do at0 = 1, 3
        read(unit_sgo,*,iostat=ios)
        if (ios /= 0) stop 'Error: malformed SGO operator.'
      end do
      read(unit_sgo,*,iostat=ios) op1
      if (ios /= 0) stop 'Error: malformed SGO operator terminator.'
    end do
    rewind(unit_sgo)
    read(unit_sgo,*)
    allocate(mgroup1(nop1, 3, 3), vgroup1(nop1, 3))
    mgroup1 = 0.0_real64
    vgroup1 = 0.0_real64
    read(unit_sgo,*) op1
    do while (op1 > 0)
      do at0 = 1, 3
        read(unit_sgo,*,iostat=ios) mgroup1(op1, at0, 1:3), vgroup1(op1, at0)
        if (ios /= 0) stop 'Error: malformed SGO operator body.'
      end do
      read(unit_sgo,*,iostat=ios) op1
      if (ios /= 0) stop 'Error: malformed SGO operator terminator.'
    end do
    close(unit_sgo)

    allocate(coords1r(nat0 * nop1, 3), spat1r(nat0 * nop1))
    at1r = 0
    do at0 = 1, nat0
      do op1 = 1, nop1
        at1r = at1r + 1
        coords1r(at1r, 1:3) = matmul(mgroup1(op1, 1:3, 1:3), coords0(at0, 1:3)) + vgroup1(op1, 1:3)
        coords1r(at1r, 1) = wrap_fractional(coords1r(at1r, 1))
        coords1r(at1r, 2) = wrap_fractional(coords1r(at1r, 2))
        coords1r(at1r, 3) = wrap_fractional(coords1r(at1r, 3))
        spat1r(at1r) = spat0(at0)
      end do
    end do
    nat1r = at1r

    allocate(coords1(nat1r, 3), spat1(nat1r))
    nat1 = 0
    do at1r = 1, nat1r
      found = .false.
      do at1i = 1, nat1
        prod = sum((coords1r(at1r, 1:3) - coords1(at1i, 1:3))**2)
        if (prod <= tol0) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        nat1 = nat1 + 1
        coords1(nat1, 1:3) = coords1r(at1r, 1:3)
        spat1(nat1) = spat1r(at1r)
      end if
    end do

    natsp1 = 0
    do at1 = 1, nat1
      natsp1(spat1(at1)) = natsp1(spat1(at1)) + 1
    end do
    nat_super = na * nb * nc
    natsp = nat_super * natsp1

    target_atini = 1
    do sp = 1, species_index - 1
      target_atini = target_atini + natsp(sp)
    end do
    target_atfin = target_atini + natsp(species_index) - 1
    if (natsp(species_index) <= 0) then
      write(error_unit,'(A)') 'Error: target species has zero multiplicity in the supercell.'
      stop 1
    end if
    if (npos > 0 .and. natsp(species_index) /= npos) then
      write(error_unit,'(A,I0,A,I0,A)') 'Error: target-site count from INSOD/SGO (', natsp(species_index), &
        ') does not match EQMATRIX target count (', npos, ').'
      stop 1
    end if

    deallocate(natsp0, natsp1, natsp, spat0, spat1, spat1r)
    deallocate(coords0, coords1, coords1r, mgroup1, vgroup1)
  end subroutine build_target_layout_from_inputs

  subroutine pme_write_level_outputs(target_level)
    integer, intent(in) :: target_level
    integer :: mm, m, choice
    integer, allocatable :: conf(:, :), omega(:)
    integer, allocatable :: conf_by_level(:, :)
    logical :: ok, have_model, copy_ok
    character(len=32) :: level_dir
    character(len=256) :: model_path, model_tmp_path, level_model_tmp_path
    real(real64) :: energy, energy_low, energy_high, low_terms(4), high_terms(4)
    ! Per-config contribution cache
    real(real64), allocatable :: all_low_terms(:,:), all_high_terms(:,:)
    real(real64), allocatable :: pred_low(:), pred_high(:)

    call require_model_initialized('pme_write_level_outputs')
    call validate_level(target_level, 'pme_write_level_outputs')
    call format_level_directory(target_level, level_dir)
    call pme_load_level_ensemble(target_level, mm, conf, omega, ok)
    if (.not. ok) then
      write(error_unit,'(A,1X,A)') 'Error: could not read target ENSEMBLE file', trim(level_dir)//'/ENSEMBLE'
      stop 1
    end if

    if (target_level > 0) then
      allocate(conf_by_level(target_level, mm))
      conf_by_level = transpose(conf(:, 1:target_level))
    else
      allocate(conf_by_level(1, mm))
      conf_by_level = 0
    end if
    if (allocated(conf)) deallocate(conf)

    ! --- Pass 1: evaluate all configs and cache per-order contributions ---
    allocate(all_low_terms(4, mm), all_high_terms(4, mm))
    allocate(pred_low(mm), pred_high(mm))
    all_low_terms  = 0.0_real64
    all_high_terms = 0.0_real64

    do m = 1, mm
      if (target_level > 0) then
        call pme_evaluate_configuration(conf_by_level(1:target_level, m), target_level, &
          energy, energy_low, energy_high, low_terms, high_terms)
      else
        energy_low  = v0_low
        energy_high = v0_high
        low_terms   = 0.0_real64
        high_terms  = 0.0_real64
      end if
      all_low_terms(:, m)  = low_terms
      all_high_terms(:, m) = high_terms
      pred_low(m)  = energy_low
      pred_high(m) = energy_high
    end do

    ! --- Determine µ scaling factors ---
    model_tmp_path = trim(pme_model_filename)//'.tmp'
    level_model_tmp_path = trim(level_dir)//'/'//trim(pme_model_filename)//'.tmp'
    call pme_preload_recalibration(target_level, ok)
    model_path = trim(pme_model_filename)
    inquire(file=trim(model_path), exist=have_model)

    if (.not. have_model .and. target_level > 0 .and. mm > 0) then
      call write_pme_model_template(model_tmp_path, all_low_terms, all_high_terms, mm, target_level)
      write(*,'(A,A)') ' > Model template written to ', trim(model_tmp_path)
      call copy_text_file(model_tmp_path, level_model_tmp_path, copy_ok)
      if (copy_ok) write(*,'(A,A)') ' > Model template copy written to ', trim(level_model_tmp_path)
      write(*,'(A)') '   Run: sod_gener.sh -choose <calib_config_list indices>, run DFT,'
      write(*,'(A,A,A)') '   add energies to nXX/ENERGIES, rename ', trim(model_tmp_path), &
        ' to '//trim(pme_model_filename)//', re-run pmesod.'
    end if

    if (have_model .and. pme_calibration_applied) then
      call write_pme_model_calibrated(model_tmp_path)
      write(*,'(A,A)') ' > Calibrated model written to ', trim(model_tmp_path)
      write(*,'(A,A,A)') '   Rename ', trim(model_tmp_path), ' to '//trim(pme_model_filename)//' to use on future runs.'
    end if

    if (have_model) then
      call write_pme_variant_outputs(pme_choice, target_level, level_dir, &
        all_low_terms, all_high_terms, mm)
    else
      if (high_base_loaded) then
        do choice = 0, 2
          call write_pme_variant_outputs(choice, target_level, level_dir, &
            all_low_terms, all_high_terms, mm)
        end do
      else
        call write_pme_variant_outputs(0, target_level, level_dir, &
          all_low_terms, all_high_terms, mm)
      end if
    end if

    if (allocated(omega))          deallocate(omega)
    if (allocated(all_low_terms))  deallocate(all_low_terms)
    if (allocated(all_high_terms)) deallocate(all_high_terms)
    if (allocated(pred_low))       deallocate(pred_low)
    if (allocated(pred_high))      deallocate(pred_high)
    if (allocated(conf_by_level))  deallocate(conf_by_level)
  end subroutine pme_write_level_outputs

  subroutine write_pme_variant_outputs(choice, target_level, level_dir, all_low_terms, all_high_terms, mm)
    integer, intent(in) :: choice, target_level, mm
    character(len=*), intent(in) :: level_dir
    real(real64), intent(in) :: all_low_terms(:,:), all_high_terms(:,:)

    integer :: m, unit_energy, ios_exec
    character(len=8) :: variant_dir
    character(len=256) :: output_dir, energy_path
    real(real64) :: energy, e_low_cal, e_high_cal, weight_low, weight_high

    if (choice == 1 .and. .not. high_base_loaded) return
    if (choice == 2 .and. .not. high_base_loaded) return

    call pme_variant_directory(choice, variant_dir)
    output_dir = trim(level_dir)//'/'//trim(variant_dir)
    call execute_command_line('mkdir -p '//trim(output_dir), exitstat=ios_exec)
    if (ios_exec /= 0) then
      write(error_unit,'(A,A)') 'Error: could not create PME output directory ', trim(output_dir)
      stop 1
    end if

    ! Copy parent nXX/ENSEMBLE so that sod_stat.sh can be run from this variant directory
    block
      character(len=256) :: src_ensemble, dst_ensemble
      logical :: have_parent_ensemble
      integer :: cp_stat
      src_ensemble = trim(level_dir)//'/ENSEMBLE'
      dst_ensemble = trim(output_dir)//'/ENSEMBLE'
      inquire(file=trim(src_ensemble), exist=have_parent_ensemble)
      if (have_parent_ensemble) then
        call execute_command_line('cp '//trim(src_ensemble)//' '//trim(dst_ensemble), exitstat=cp_stat)
        if (cp_stat /= 0) &
          write(*,'(A,A)') ' > Warning: could not copy ENSEMBLE to ', trim(output_dir)
      end if
    end block

    energy_path = trim(output_dir)//'/ENERGIES'
    open(newunit=unit_energy, file=trim(energy_path), status='replace', action='write')
    do m = 1, mm
      call compute_pme_variant_energy(choice, target_level, all_low_terms(:, m), all_high_terms(:, m), &
        energy, e_low_cal, e_high_cal)
      write(unit_energy,'(F20.10)') energy
    end do
    close(unit_energy)

    if (choice == 2 .and. high_base_loaded) then
      call level_blend_weights(target_level, npos, weight_low, weight_high)
      write(*,'(A,A,A,F7.5,A,F7.5,A,F6.3,A,F6.3,A,A)') ' > ', trim(variant_dir), &
        ' energies written to ', weight_low, '×E_low + ', weight_high, '×E_high (alpha=', &
        alpha_hybrid, ', eta=', eta_hybrid, ') → ', trim(energy_path)
    else
      write(*,'(A,A,A,A)') ' > ', trim(variant_dir), ' energies written to ', trim(energy_path)
    end if
  end subroutine write_pme_variant_outputs

  subroutine compute_pme_variant_energy(choice, level_val, low_terms, high_terms, energy, e_low_cal, e_high_cal)
    integer, intent(in) :: choice, level_val
    real(real64), intent(in) :: low_terms(4), high_terms(4)
    real(real64), intent(out) :: energy, e_low_cal, e_high_cal
    integer :: eff_low_ord, eff_high_ord
    real(real64) :: weight_low, weight_high

    eff_low_ord  = min(max_low_order,  pme_order_cap)
    eff_high_ord = min(max_high_order, pme_order_cap)

    if (eff_low_ord >= 1) then
      e_low_cal = v0_low + eps_low(0) + dot_product(eps_low(1:eff_low_ord), low_terms(1:eff_low_ord))
    else
      e_low_cal = v0_low + eps_low(0)
    end if

    if (high_base_loaded .and. eff_high_ord >= 1) then
      e_high_cal = v0_high + eps_high(0) + dot_product(eps_high(1:eff_high_ord), high_terms(1:eff_high_ord))
    else
      e_high_cal = e_low_cal
    end if

    select case (choice)
    case (0)
      energy = e_low_cal
    case (1)
      energy = e_high_cal
    case default
      if (high_base_loaded) then
        call level_blend_weights(level_val, npos, weight_low, weight_high)
        energy = weight_low * e_low_cal + weight_high * e_high_cal
      else
        energy = e_low_cal
      end if
    end select
  end subroutine compute_pme_variant_energy

  subroutine pme_variant_directory(choice, dirname)
    integer, intent(in) :: choice
    character(len=*), intent(out) :: dirname

    select case (choice)
    case (0)
      dirname = 'PME0'
    case (1)
      dirname = 'PME1'
    case default
      dirname = 'PMEh'
    end select
  end subroutine pme_variant_directory

  subroutine pme_variant_dir_from_model(filename, dirname, ok)
    !  Read just the PME choice (first data line) from a pme.model file and
    !  return the corresponding variant directory name ('PME0'/'PME1'/'PMEh').
    !  Lets standalone tools (e.g. mcstatsod) locate the variant subdirectory
    !  without building the full Hamiltonian.  ok = .false. if the file is
    !  missing or the choice line is unreadable (dirname defaults to 'PMEh').
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: dirname
    logical,          intent(out) :: ok
    integer :: unit_in, ios, choice
    character(len=512) :: line
    logical :: lk

    ok = .false.
    dirname = 'PMEh'
    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    call next_data_line(unit_in, line, lk)
    close(unit_in)
    if (.not. lk) return
    read(line, *, iostat=ios) choice
    if (ios /= 0 .or. choice < 0 .or. choice > 2) return
    call pme_variant_directory(choice, dirname)
    ok = .true.
  end subroutine pme_variant_dir_from_model

  subroutine pme_evaluate_configuration(conf_row, level, energy, energy_low, energy_high, low_terms, high_terms)
    integer, intent(in) :: conf_row(:)
    integer, intent(in) :: level
    real(real64), intent(out) :: energy, energy_low, energy_high
    real(real64), intent(out) :: low_terms(4), high_terms(4)
    integer :: p, q, r, s, hole_count, eff_low_ord, eff_high_ord
    real(real64) :: weight_low, weight_high

    eff_low_ord  = min(max_low_order,  pme_order_cap)
    eff_high_ord = min(max_high_order, pme_order_cap)

    low_terms = 0.0_real64
    high_terms = 0.0_real64
    energy_low = v0_low

    if (level > 0) then
      do p = 1, level
        low_terms(1) = low_terms(1) + V1_low(conf_row(p))
      end do
      if (level > 1 .and. eff_low_ord >= 2) then
        do p = 1, level - 1
          do q = p + 1, level
            low_terms(2) = low_terms(2) + V2_low(conf_row(p), conf_row(q))
          end do
        end do
      end if
      if (level > 2 .and. eff_low_ord >= 3 .and. allocated(V3_i_low)) then
        do p = 1, level - 2
          do q = p + 1, level - 1
            do r = q + 1, level
              low_terms(3) = low_terms(3) + get_sorted_triplet_value(conf_row(p), conf_row(q), conf_row(r), V3_key_low, V3_val_low)
            end do
          end do
        end do
      end if
      if (level > 3 .and. eff_low_ord >= 4 .and. allocated(V4_i_low)) then
        do p = 1, level - 3
          do q = p + 1, level - 2
            do r = q + 1, level - 1
              do s = r + 1, level
                low_terms(4) = low_terms(4) + get_sorted_quad_value(conf_row(p), conf_row(q), conf_row(r), conf_row(s), &
                  V4_key_low, V4_val_low)
              end do
            end do
          end do
        end do
      end if
    end if
    energy_low = energy_low + sum(low_terms(1:eff_low_ord))

    energy_high = energy_low
    if (high_base_loaded) then
      hole_count = npos - level
      call ensure_work_buffers()
      call build_hole_configuration(conf_row, level, work_holes, hole_count)
      high_terms = 0.0_real64
      energy_high = v0_high
      if (hole_count > 0) then
        do p = 1, hole_count
          high_terms(1) = high_terms(1) + V1_high(work_holes(p))
        end do
        if (hole_count > 1 .and. eff_high_ord >= 2) then
          do p = 1, hole_count - 1
            do q = p + 1, hole_count
              high_terms(2) = high_terms(2) + V2_high(work_holes(p), work_holes(q))
            end do
          end do
        end if
        if (hole_count > 2 .and. eff_high_ord >= 3 .and. allocated(V3_i_high)) then
          do p = 1, hole_count - 2
            do q = p + 1, hole_count - 1
              do r = q + 1, hole_count
                high_terms(3) = high_terms(3) + get_sorted_triplet_value(work_holes(p), work_holes(q), work_holes(r), &
                  V3_key_high, V3_val_high)
              end do
            end do
          end do
        end if
        if (hole_count > 3 .and. eff_high_ord >= 4 .and. allocated(V4_i_high)) then
          do p = 1, hole_count - 3
            do q = p + 1, hole_count - 2
              do r = q + 1, hole_count - 1
                do s = r + 1, hole_count
                  high_terms(4) = high_terms(4) + &
                    get_sorted_quad_value(work_holes(p), work_holes(q), &
                    work_holes(r), work_holes(s), V4_key_high, V4_val_high)
                end do
              end do
            end do
          end do
        end if
      end if
      energy_high = energy_high + sum(high_terms(1:eff_high_ord))
    end if

    select case (pme_choice)
    case (0)
      energy = energy_low
    case (1)
      if (high_base_loaded) then
        energy = energy_high
      else
        energy = energy_low
      end if
    case default  ! 2 = PMEh
      if (high_base_loaded) then
        call level_blend_weights(level, npos, weight_low, weight_high)
        energy = weight_low * energy_low + weight_high * energy_high
      else
        energy = energy_low
      end if
    end select
  end subroutine pme_evaluate_configuration

  subroutine ensure_work_buffers()
    if (.not. allocated(work_is_ge)) allocate(work_is_ge(npos))
    if (size(work_is_ge) /= npos) then
      deallocate(work_is_ge)
      allocate(work_is_ge(npos))
    end if
    if (.not. allocated(work_holes)) allocate(work_holes(max(1, npos)))
    if (size(work_holes) < max(1, npos)) then
      deallocate(work_holes)
      allocate(work_holes(max(1, npos)))
    end if
  end subroutine ensure_work_buffers

  subroutine read_eqmatrix_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_in, ios, op

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,'(A,1X,A)') 'Error: could not open EQMATRIX file', trim(filename)
      stop 1
    end if

    read(unit_in,*,iostat=ios) nop, npos
    if (ios /= 0 .or. nop <= 0 .or. npos <= 0) then
      close(unit_in)
      write(error_unit,'(A)') 'Error: invalid EQMATRIX header.'
      stop 1
    end if
    allocate(eqmatrix(nop, npos))
    do op = 1, nop
      read(unit_in,*,iostat=ios) eqmatrix(op, 1:npos)
      if (ios /= 0) then
        close(unit_in)
        write(error_unit,'(A)') 'Error: invalid EQMATRIX row.'
        stop 1
      end if
    end do
    close(unit_in)
  end subroutine read_eqmatrix_file

  subroutine load_low_side_model()
    integer :: level, low_max_level, mm, m, op, i1, i2, i3, i4
    integer, allocatable :: conf(:, :), omega(:)
    real(real64), allocatable :: energies(:)
    logical :: ok
    real(real64) :: residual

    call read_v0_from_level(0, v0_low, ok)
    if (.not. ok) then
      write(error_unit,'(A)') 'Error: n00/ENERGIES is required to build the low-side expansion.'
      stop 1
    end if

    call read_level_data(1, mm, conf, omega, energies, ok)
    if (.not. ok) then
      write(error_unit,'(A)') 'Error: n01/ENSEMBLE and n01/ENERGIES are required to build the low-side expansion.'
      stop 1
    end if
    allocate(V1_low(npos))
    V1_low = 0.0_real64
    allocate(v1res_low(mm))
    do m = 1, mm
      residual = energies(m) - v0_low
      v1res_low(m) = residual
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        call assign_scalar(V1_low(i1), residual)
      end do
    end do
    nv1_low = mm
    max_low_order = 1
    call cleanup_level_data(conf, omega, energies)

    call read_level_data(2, mm, conf, omega, energies, ok)
    if (.not. ok) then
      write(error_unit,'(A)') 'Error: n02/ENSEMBLE and n02/ENERGIES are required to build the low-side expansion.'
      stop 1
    end if
    allocate(V2_low(npos, npos))
    V2_low = 0.0_real64
    allocate(v2res_low(mm))
    do m = 1, mm
      i1 = eqmatrix(1, conf(m, 1))
      i2 = eqmatrix(1, conf(m, 2))
      call sort_pair(i1, i2)
      residual = energies(m) - v0_low - V1_low(i1) - V1_low(i2)
      v2res_low(m) = residual
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        i2 = eqmatrix(op, conf(m, 2))
        call sort_pair(i1, i2)
        call assign_pair_term(V2_low, i1, i2, residual)
      end do
    end do
    nv2_low = mm
    max_low_order = 2
    call cleanup_level_data(conf, omega, energies)

    low_max_level = max_motif_order
    if (pme_target_level > 0) low_max_level = min(max_motif_order, pme_target_level - 1)
    do level = 3, low_max_level
      call read_level_data(level, mm, conf, omega, energies, ok)
      if (.not. ok) exit
      if (level == 4 .and. allocated(V3_i_low) .and. .not. allocated(V3_key_low)) then
        call prepare_triplet_lookup(V3_i_low, V3_j_low, V3_k_low, V3_val_low, V3_key_low)
      end if
      if (level == 3) allocate(v3res_low(mm))
      if (level == 4) allocate(v4res_low(mm))
      do m = 1, mm
        if (level == 3) then
          i1 = eqmatrix(1, conf(m, 1))
          i2 = eqmatrix(1, conf(m, 2))
          i3 = eqmatrix(1, conf(m, 3))
          call sort_triplet(i1, i2, i3)
          residual = energies(m) - v0_low - V1_low(i1) - V1_low(i2) - V1_low(i3) &
                   - V2_low(i1, i2) - V2_low(i1, i3) - V2_low(i2, i3)
          v3res_low(m) = residual
          do op = 1, nop
            i1 = eqmatrix(op, conf(m, 1))
            i2 = eqmatrix(op, conf(m, 2))
            i3 = eqmatrix(op, conf(m, 3))
            call sort_triplet(i1, i2, i3)
            call assign_triplet_term(V3_i_low, V3_j_low, V3_k_low, V3_val_low, i1, i2, i3, residual)
          end do
        else
          i1 = eqmatrix(1, conf(m, 1))
          i2 = eqmatrix(1, conf(m, 2))
          i3 = eqmatrix(1, conf(m, 3))
          i4 = eqmatrix(1, conf(m, 4))
          call sort_quad(i1, i2, i3, i4)
          residual = energies(m) - v0_low
          residual = residual - V1_low(i1) - V1_low(i2) - V1_low(i3) - V1_low(i4)
          residual = residual - V2_low(i1, i2) - V2_low(i1, i3) - V2_low(i1, i4) &
                   - V2_low(i2, i3) - V2_low(i2, i4) - V2_low(i3, i4)
          residual = residual - get_sorted_triplet_value(i1, i2, i3, V3_key_low, V3_val_low)
          residual = residual - get_sorted_triplet_value(i1, i2, i4, V3_key_low, V3_val_low)
          residual = residual - get_sorted_triplet_value(i1, i3, i4, V3_key_low, V3_val_low)
          residual = residual - get_sorted_triplet_value(i2, i3, i4, V3_key_low, V3_val_low)
          v4res_low(m) = residual
          do op = 1, nop
            i1 = eqmatrix(op, conf(m, 1))
            i2 = eqmatrix(op, conf(m, 2))
            i3 = eqmatrix(op, conf(m, 3))
            i4 = eqmatrix(op, conf(m, 4))
            call sort_quad(i1, i2, i3, i4)
            call assign_quad_term(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, i1, i2, i3, i4, residual)
          end do
        end if
      end do
      if (level == 3) nv3_low = mm
      if (level == 4) nv4_low = mm
      max_low_order = level
      call cleanup_level_data(conf, omega, energies)
    end do

    call prepare_triplet_lookup(V3_i_low, V3_j_low, V3_k_low, V3_val_low, V3_key_low)
    call prepare_quad_lookup(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, V4_key_low)
  end subroutine load_low_side_model

  subroutine load_high_side_model()
    integer :: ge_count, hole_count, mm, m, op
    integer :: i1, i2, i3, i4
    integer, allocatable :: conf(:, :), omega(:), holes(:, :)
    real(real64), allocatable :: energies(:)
    logical :: ok
    real(real64) :: residual

    call read_v0_from_level(npos, v0_high, ok)
    if (.not. ok) then
      high_base_loaded = .false.
      return
    end if

    high_base_loaded = .true.
    allocate(V1_high(npos))
    allocate(V2_high(npos, npos))
    V1_high = 0.0_real64
    V2_high = 0.0_real64
    max_high_order = 0

    do hole_count = 1, min(max_motif_order, npos)
      ge_count = npos - hole_count
      if (pme_target_level > 0 .and. ge_count <= pme_target_level) exit
      call read_level_data(ge_count, mm, conf, omega, energies, ok)
      if (.not. ok) exit
      if (hole_count == 4 .and. allocated(V3_i_high) .and. .not. allocated(V3_key_high)) then
        call prepare_triplet_lookup(V3_i_high, V3_j_high, V3_k_high, V3_val_high, V3_key_high)
      end if
      allocate(holes(mm, max(1, hole_count)))
      call build_hole_configurations(conf, mm, ge_count, hole_count, holes)
      if (hole_count == 1) allocate(v1res_high(mm))
      if (hole_count == 2) allocate(v2res_high(mm))
      if (hole_count == 3) allocate(v3res_high(mm))
      if (hole_count == 4) allocate(v4res_high(mm))

      do m = 1, mm
        i1 = eqmatrix(1, holes(m, 1))
        if (hole_count == 1) then
          residual = energies(m) - v0_high
          v1res_high(m) = residual
        else if (hole_count == 2) then
          i2 = eqmatrix(1, holes(m, 2))
          call sort_pair(i1, i2)
          residual = energies(m) - v0_high - V1_high(i1) - V1_high(i2)
          v2res_high(m) = residual
        else if (hole_count == 3) then
          i2 = eqmatrix(1, holes(m, 2))
          i3 = eqmatrix(1, holes(m, 3))
          call sort_triplet(i1, i2, i3)
          residual = energies(m) - v0_high - V1_high(i1) - V1_high(i2) - V1_high(i3) &
                   - V2_high(i1, i2) - V2_high(i1, i3) - V2_high(i2, i3)
          v3res_high(m) = residual
        else
          i2 = eqmatrix(1, holes(m, 2))
          i3 = eqmatrix(1, holes(m, 3))
          i4 = eqmatrix(1, holes(m, 4))
          call sort_quad(i1, i2, i3, i4)
          residual = energies(m) - v0_high
          residual = residual - V1_high(i1) - V1_high(i2) - V1_high(i3) - V1_high(i4)
          residual = residual - V2_high(i1, i2) - V2_high(i1, i3) - V2_high(i1, i4) &
                   - V2_high(i2, i3) - V2_high(i2, i4) - V2_high(i3, i4)
          residual = residual - get_sorted_triplet_value(i1, i2, i3, V3_key_high, V3_val_high)
          residual = residual - get_sorted_triplet_value(i1, i2, i4, V3_key_high, V3_val_high)
          residual = residual - get_sorted_triplet_value(i1, i3, i4, V3_key_high, V3_val_high)
          residual = residual - get_sorted_triplet_value(i2, i3, i4, V3_key_high, V3_val_high)
          v4res_high(m) = residual
        end if
        do op = 1, nop
          i1 = eqmatrix(op, holes(m, 1))
          if (hole_count == 1) then
            call assign_scalar(V1_high(i1), residual)
          else if (hole_count == 2) then
            i2 = eqmatrix(op, holes(m, 2))
            call sort_pair(i1, i2)
            call assign_pair_term(V2_high, i1, i2, residual)
          else if (hole_count == 3) then
            i2 = eqmatrix(op, holes(m, 2))
            i3 = eqmatrix(op, holes(m, 3))
            call sort_triplet(i1, i2, i3)
            call assign_triplet_term(V3_i_high, V3_j_high, V3_k_high, V3_val_high, i1, i2, i3, residual)
          else
            i2 = eqmatrix(op, holes(m, 2))
            i3 = eqmatrix(op, holes(m, 3))
            i4 = eqmatrix(op, holes(m, 4))
            call sort_quad(i1, i2, i3, i4)
            call assign_quad_term(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, i1, i2, i3, i4, residual)
          end if
        end do
      end do

      if (hole_count == 1) nv1_high = mm
      if (hole_count == 2) nv2_high = mm
      if (hole_count == 3) nv3_high = mm
      if (hole_count == 4) nv4_high = mm
      max_high_order = hole_count
      deallocate(holes)
      call cleanup_level_data(conf, omega, energies)
    end do

    call prepare_triplet_lookup(V3_i_high, V3_j_high, V3_k_high, V3_val_high, V3_key_high)
    call prepare_quad_lookup(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, V4_key_high)
  end subroutine load_high_side_model

  subroutine read_v0_from_level(level, energy_value, ok)
    integer, intent(in) :: level
    real(real64), intent(out) :: energy_value
    logical, intent(out) :: ok
    character(len=32) :: level_dir
    real(real64) :: energies(1)
    integer :: n_missing

    call format_level_directory(level, level_dir)
    call read_energies_file(trim(level_dir)//'/ENERGIES', 1, energies, ok, n_missing)
    if (ok .and. n_missing == 0) then
      energy_value = energies(1)
    else
      ok = .false.
      energy_value = 0.0_real64
    end if
  end subroutine read_v0_from_level

  subroutine read_level_data(level, mm, conf, omega, energies, ok)
    integer, intent(in) :: level
    integer, intent(out) :: mm
    integer, allocatable, intent(out) :: conf(:, :), omega(:)
    real(real64), allocatable, intent(out) :: energies(:)
    logical, intent(out) :: ok
    character(len=32) :: level_dir
    integer :: level_from_file, total_sites

    call format_level_directory(level, level_dir)
    call read_ensemble_level(trim(level_dir)//'/ENSEMBLE', level_from_file, total_sites, mm, conf, omega, ok)
    if (.not. ok) return
    if (level_from_file /= level) then
      ok = .false.
      call cleanup_level_data(conf, omega, energies)
      return
    end if
    if (total_sites /= npos) then
      ok = .false.
      call cleanup_level_data(conf, omega, energies)
      return
    end if
    allocate(energies(mm))
    block
      integer :: n_missing
      call read_energies_file(trim(level_dir)//'/ENERGIES', mm, energies, ok, n_missing)
      if (ok .and. n_missing > 0) then
        write(error_unit,'(A,I0,A,A)') 'Error: missing energies for ', n_missing, &
          ' configuration(s) in ', trim(level_dir)//'/ENERGIES'
        ok = .false.
      end if
    end block
    if (.not. ok) call cleanup_level_data(conf, omega, energies)
  end subroutine read_level_data

  subroutine read_level_conf_only(level, mm, conf, omega, ok)
    integer, intent(in)  :: level
    integer, intent(out) :: mm
    integer, allocatable, intent(out) :: conf(:, :), omega(:)
    logical, intent(out) :: ok
    character(len=32) :: level_dir
    integer :: level_from_file, total_sites

    call format_level_directory(level, level_dir)
    call read_ensemble_level(trim(level_dir)//'/ENSEMBLE', level_from_file, total_sites, mm, conf, omega, ok)
    if (.not. ok) return
    if (level_from_file /= level .or. total_sites /= npos) then
      ok = .false.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      mm = 0
    end if
  end subroutine read_level_conf_only

  subroutine cleanup_level_data(conf, omega, energies)
    integer, allocatable, intent(inout) :: conf(:, :), omega(:)
    real(real64), allocatable, intent(inout) :: energies(:)

    if (allocated(conf)) deallocate(conf)
    if (allocated(omega)) deallocate(omega)
    if (allocated(energies)) deallocate(energies)
  end subroutine cleanup_level_data

  subroutine read_ensemble_level(filename, level, total_sites, mm, conf, omega, ok)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: level, total_sites, mm
    integer, allocatable, intent(out) :: conf(:, :), omega(:)
    logical, intent(out) :: ok

    integer :: ntarget_l
    integer, allocatable :: nsfl(:), npst(:), ic(:,:)
    real(real64) :: tsampling_l

    ok = .false.
    level = 0
    total_sites = 0
    mm = 0

    call read_ensemble(filename, ntarget_l, nsfl, npst, mm, ic, omega, tsampling_l, ok)
    if (.not. ok) return
    level = nsfl(1)
    total_sites = npst(1)
    allocate(conf(mm, max(1, level)))
    conf = 0
    if (level > 0) conf(1:mm, 1:level) = ic(1:mm, 1:level)
    deallocate(ic, nsfl, npst)
    if (level > 0) then
      call localize_configurations(conf, mm, level, ok)
      if (.not. ok) deallocate(conf, omega)
    end if
  end subroutine read_ensemble_level

  subroutine localize_configurations(conf, mm, level, ok)
    integer, intent(in) :: mm, level
    integer, intent(inout) :: conf(:,:)
    logical, intent(out) :: ok
    integer :: m, pos, raw_index

    ok = .true.
    if (target_atini <= 0 .or. target_atfin < target_atini) then
      ok = .false.
      return
    end if

    do m = 1, mm
      do pos = 1, level
        raw_index = conf(m, pos)
        if (raw_index < target_atini .or. raw_index > target_atfin) then
          ok = .false.
          return
        end if
        conf(m, pos) = raw_index - target_atini + 1
        if (conf(m, pos) < 1 .or. conf(m, pos) > npos) then
          ok = .false.
          return
        end if
      end do
    end do
  end subroutine localize_configurations


  subroutine build_hole_configurations(conf, mm, ge_count, hole_count, holes)
    integer, intent(in) :: mm, ge_count, hole_count
    integer, intent(in) :: conf(mm, max(1, ge_count))
    integer, intent(out) :: holes(mm, max(1, hole_count))
    integer :: m, i, idx

    call ensure_work_buffers()
    do m = 1, mm
      work_is_ge = .false.
      if (ge_count > 0) then
        do i = 1, ge_count
          work_is_ge(conf(m, i)) = .true.
        end do
      end if
      idx = 0
      do i = 1, npos
        if (.not. work_is_ge(i)) then
          idx = idx + 1
          if (idx <= hole_count) holes(m, idx) = i
        end if
      end do
    end do
  end subroutine build_hole_configurations

  subroutine build_hole_configuration(conf_row, level, holes, hole_count)
    integer, intent(in) :: conf_row(:)
    integer, intent(in) :: level, hole_count
    integer, intent(out) :: holes(:)
    integer :: i, idx

    call ensure_work_buffers()
    work_is_ge = .false.
    if (level > 0) then
      do i = 1, level
        work_is_ge(conf_row(i)) = .true.
      end do
    end if
    idx = 0
    do i = 1, npos
      if (.not. work_is_ge(i)) then
        idx = idx + 1
        if (idx <= hole_count) holes(idx) = i
      end if
    end do
  end subroutine build_hole_configuration

  subroutine assign_scalar(slot, value)
    real(real64), intent(inout) :: slot
    real(real64), intent(in) :: value
    if (abs(slot) <= assign_tol) then
      slot = value
    else if (abs(slot - value) > assign_tol) then
      slot = 0.5_real64 * (slot + value)
    end if
  end subroutine assign_scalar

  subroutine assign_pair_term(table, i, j, value)
    real(real64), intent(inout) :: table(:, :)
    integer, intent(in) :: i, j
    real(real64), intent(in) :: value
    real(real64) :: merged

    merged = table(i, j)
    if (abs(merged) <= assign_tol) then
      merged = value
    else if (abs(merged - value) > assign_tol) then
      merged = 0.5_real64 * (merged + value)
    end if
    table(i, j) = merged
    table(j, i) = merged
  end subroutine assign_pair_term

  subroutine assign_triplet_term(ai, aj, ak, aval, i, j, k, value)
    integer, allocatable, intent(inout) :: ai(:), aj(:), ak(:)
    real(real64), allocatable, intent(inout) :: aval(:)
    integer, intent(in) :: i, j, k
    real(real64), intent(in) :: value
    integer :: idx, nold
    integer, allocatable :: tmp_i(:), tmp_j(:), tmp_k(:)
    real(real64), allocatable :: tmp_v(:)

    idx = find_triplet_index(i, j, k, ai, aj, ak)
    if (idx > 0) then
      call assign_scalar(aval(idx), value)
      return
    end if

    if (.not. allocated(ai)) then
      allocate(ai(1), aj(1), ak(1), aval(1))
      ai(1) = i
      aj(1) = j
      ak(1) = k
      aval(1) = value
      return
    end if

    nold = size(ai)
    allocate(tmp_i(nold + 1), tmp_j(nold + 1), tmp_k(nold + 1), tmp_v(nold + 1))
    tmp_i(1:nold) = ai
    tmp_j(1:nold) = aj
    tmp_k(1:nold) = ak
    tmp_v(1:nold) = aval
    tmp_i(nold + 1) = i
    tmp_j(nold + 1) = j
    tmp_k(nold + 1) = k
    tmp_v(nold + 1) = value
    call move_alloc(tmp_i, ai)
    call move_alloc(tmp_j, aj)
    call move_alloc(tmp_k, ak)
    call move_alloc(tmp_v, aval)
  end subroutine assign_triplet_term

  subroutine assign_quad_term(ai, aj, ak, al, aval, i, j, k, l, value)
    integer, allocatable, intent(inout) :: ai(:), aj(:), ak(:), al(:)
    real(real64), allocatable, intent(inout) :: aval(:)
    integer, intent(in) :: i, j, k, l
    real(real64), intent(in) :: value
    integer :: idx, nold
    integer, allocatable :: tmp_i(:), tmp_j(:), tmp_k(:), tmp_l(:)
    real(real64), allocatable :: tmp_v(:)

    idx = find_quad_index(i, j, k, l, ai, aj, ak, al)
    if (idx > 0) then
      call assign_scalar(aval(idx), value)
      return
    end if

    if (.not. allocated(ai)) then
      allocate(ai(1), aj(1), ak(1), al(1), aval(1))
      ai(1) = i
      aj(1) = j
      ak(1) = k
      al(1) = l
      aval(1) = value
      return
    end if

    nold = size(ai)
    allocate(tmp_i(nold + 1), tmp_j(nold + 1), tmp_k(nold + 1), tmp_l(nold + 1), tmp_v(nold + 1))
    tmp_i(1:nold) = ai
    tmp_j(1:nold) = aj
    tmp_k(1:nold) = ak
    tmp_l(1:nold) = al
    tmp_v(1:nold) = aval
    tmp_i(nold + 1) = i
    tmp_j(nold + 1) = j
    tmp_k(nold + 1) = k
    tmp_l(nold + 1) = l
    tmp_v(nold + 1) = value
    call move_alloc(tmp_i, ai)
    call move_alloc(tmp_j, aj)
    call move_alloc(tmp_k, ak)
    call move_alloc(tmp_l, al)
    call move_alloc(tmp_v, aval)
  end subroutine assign_quad_term

  integer function find_triplet_index(i, j, k, ai, aj, ak) result(idx)
    integer, intent(in) :: i, j, k
    integer, allocatable, intent(in) :: ai(:), aj(:), ak(:)
    integer :: n

    idx = 0
    if (.not. allocated(ai)) return
    do n = 1, size(ai)
      if (ai(n) == i .and. aj(n) == j .and. ak(n) == k) then
        idx = n
        return
      end if
    end do
  end function find_triplet_index

  integer function find_quad_index(i, j, k, l, ai, aj, ak, al) result(idx)
    integer, intent(in) :: i, j, k, l
    integer, allocatable, intent(in) :: ai(:), aj(:), ak(:), al(:)
    integer :: n

    idx = 0
    if (.not. allocated(ai)) return
    do n = 1, size(ai)
      if (ai(n) == i .and. aj(n) == j .and. ak(n) == k .and. al(n) == l) then
        idx = n
        return
      end if
    end do
  end function find_quad_index

  subroutine prepare_triplet_lookup(ai, aj, ak, aval, akey)
    integer, allocatable, intent(inout) :: ai(:), aj(:), ak(:)
    real(real64), allocatable, intent(inout) :: aval(:)
    integer(int64), allocatable, intent(inout) :: akey(:)
    integer :: n

    if (allocated(akey)) deallocate(akey)
    if (.not. allocated(ai)) return
    n = size(ai)
    if (n <= 0) return
    allocate(akey(n))
    do n = 1, size(ai)
      akey(n) = pack_triplet_key(ai(n), aj(n), ak(n))
    end do
    call sort_triplet_lookup(ai, aj, ak, aval, akey, 1, size(ai))
  end subroutine prepare_triplet_lookup

  subroutine prepare_quad_lookup(ai, aj, ak, al, aval, akey)
    integer, allocatable, intent(inout) :: ai(:), aj(:), ak(:), al(:)
    real(real64), allocatable, intent(inout) :: aval(:)
    integer(int64), allocatable, intent(inout) :: akey(:)
    integer :: n

    if (allocated(akey)) deallocate(akey)
    if (.not. allocated(ai)) return
    n = size(ai)
    if (n <= 0) return
    allocate(akey(n))
    do n = 1, size(ai)
      akey(n) = pack_quad_key(ai(n), aj(n), ak(n), al(n))
    end do
    call sort_quad_lookup(ai, aj, ak, al, aval, akey, 1, size(ai))
  end subroutine prepare_quad_lookup

  recursive subroutine sort_triplet_lookup(ai, aj, ak, aval, akey, left, right)
    integer, intent(inout) :: ai(:), aj(:), ak(:)
    real(real64), intent(inout) :: aval(:)
    integer(int64), intent(inout) :: akey(:)
    integer, intent(in) :: left, right
    integer :: i, j, ti
    real(real64) :: tv
    integer(int64) :: pivot, tk

    if (left >= right) return
    pivot = akey((left + right) / 2)
    i = left
    j = right
    do
      do while (akey(i) < pivot)
        i = i + 1
      end do
      do while (akey(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tk = akey(i); akey(i) = akey(j); akey(j) = tk
        ti = ai(i); ai(i) = ai(j); ai(j) = ti
        ti = aj(i); aj(i) = aj(j); aj(j) = ti
        ti = ak(i); ak(i) = ak(j); ak(j) = ti
        tv = aval(i); aval(i) = aval(j); aval(j) = tv
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call sort_triplet_lookup(ai, aj, ak, aval, akey, left, j)
    if (i < right) call sort_triplet_lookup(ai, aj, ak, aval, akey, i, right)
  end subroutine sort_triplet_lookup

  recursive subroutine sort_quad_lookup(ai, aj, ak, al, aval, akey, left, right)
    integer, intent(inout) :: ai(:), aj(:), ak(:), al(:)
    real(real64), intent(inout) :: aval(:)
    integer(int64), intent(inout) :: akey(:)
    integer, intent(in) :: left, right
    integer :: i, j, ti
    real(real64) :: tv
    integer(int64) :: pivot, tk

    if (left >= right) return
    pivot = akey((left + right) / 2)
    i = left
    j = right
    do
      do while (akey(i) < pivot)
        i = i + 1
      end do
      do while (akey(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tk = akey(i); akey(i) = akey(j); akey(j) = tk
        ti = ai(i); ai(i) = ai(j); ai(j) = ti
        ti = aj(i); aj(i) = aj(j); aj(j) = ti
        ti = ak(i); ak(i) = ak(j); ak(j) = ti
        ti = al(i); al(i) = al(j); al(j) = ti
        tv = aval(i); aval(i) = aval(j); aval(j) = tv
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call sort_quad_lookup(ai, aj, ak, al, aval, akey, left, j)
    if (i < right) call sort_quad_lookup(ai, aj, ak, al, aval, akey, i, right)
  end subroutine sort_quad_lookup

  integer(int64) function pack_triplet_key(i, j, k) result(key)
    integer, intent(in) :: i, j, k
    integer(int64) :: base

    base = int(npos, int64)
    key = ((int(i - 1, int64) * base + int(j - 1, int64)) * base) + int(k, int64)
  end function pack_triplet_key

  integer(int64) function pack_quad_key(i, j, k, l) result(key)
    integer, intent(in) :: i, j, k, l
    integer(int64) :: base

    base = int(npos, int64)
    key = (((int(i - 1, int64) * base + int(j - 1, int64)) * base + int(k - 1, int64)) * base) + int(l, int64)
  end function pack_quad_key

  integer function binary_search_key(key, akey) result(idx)
    integer(int64), intent(in) :: key
    integer(int64), allocatable, intent(in) :: akey(:)
    integer :: left, right, mid

    idx = 0
    if (.not. allocated(akey)) return
    left = 1
    right = size(akey)
    do while (left <= right)
      mid = left + (right - left) / 2
      if (akey(mid) == key) then
        idx = mid
        return
      else if (akey(mid) < key) then
        left = mid + 1
      else
        right = mid - 1
      end if
    end do
  end function binary_search_key

  real(real64) function get_sorted_triplet_value(i, j, k, akey, aval) result(value)
    integer, intent(in) :: i, j, k
    integer(int64), allocatable, intent(in) :: akey(:)
    real(real64), allocatable, intent(in) :: aval(:)
    integer :: idx

    value = 0.0_real64
    if (.not. allocated(akey)) return
    idx = binary_search_key(pack_triplet_key(i, j, k), akey)
    if (idx > 0) value = aval(idx)
  end function get_sorted_triplet_value

  real(real64) function get_sorted_quad_value(i, j, k, l, akey, aval) result(value)
    integer, intent(in) :: i, j, k, l
    integer(int64), allocatable, intent(in) :: akey(:)
    real(real64), allocatable, intent(in) :: aval(:)
    integer :: idx

    value = 0.0_real64
    if (.not. allocated(akey)) return
    idx = binary_search_key(pack_quad_key(i, j, k, l), akey)
    if (idx > 0) value = aval(idx)
  end function get_sorted_quad_value

  subroutine format_level_directory(level, dirname)
    integer, intent(in) :: level
    character(len=*), intent(out) :: dirname
    character(len=32) :: tag
    integer :: width

    width = max(2, count_digits(max(1, npos)))
    call zero_pad_integer(level, width, tag)
    dirname = 'n'//trim(tag)
  end subroutine format_level_directory

  integer function count_digits(value) result(ndig)
    integer, intent(in) :: value
    integer :: tmp

    tmp = abs(value)
    ndig = 1
    do while (tmp >= 10)
      tmp = tmp / 10
      ndig = ndig + 1
    end do
  end function count_digits

  subroutine zero_pad_integer(value, width, text)
    integer, intent(in) :: value, width
    character(len=*), intent(out) :: text
    character(len=32) :: raw
    integer :: pad

    write(raw,'(I0)') value
    pad = max(0, width - len_trim(raw))
    text = repeat('0', pad)//trim(raw)
  end subroutine zero_pad_integer

  subroutine sort_pair(i, j)
    integer, intent(inout) :: i, j
    integer :: tmp
    if (i > j) then
      tmp = i
      i = j
      j = tmp
    end if
  end subroutine sort_pair

  subroutine sort_triplet(i, j, k)
    integer, intent(inout) :: i, j, k
    call sort_pair(i, j)
    call sort_pair(j, k)
    call sort_pair(i, j)
  end subroutine sort_triplet

  subroutine sort_quad(i, j, k, l)
    integer, intent(inout) :: i, j, k, l
    call sort_pair(i, j)
    call sort_pair(k, l)
    call sort_pair(i, k)
    call sort_pair(j, l)
    call sort_pair(j, k)
  end subroutine sort_quad

  subroutine next_data_line(unit_in, line, ok)
    integer, intent(in) :: unit_in
    character(len=*), intent(out) :: line
    logical, intent(out) :: ok
    integer :: ios

    ok = .false.
    line = ''
    do
      read(unit_in,'(A)',iostat=ios) line
      if (ios /= 0) return
      line = adjustl(trim(line))
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      ok = .true.
      return
    end do
  end subroutine next_data_line

  pure real(real64) function wrap_fractional(x) result(corx)
    real(real64), intent(in) :: x
    real(real64), parameter :: tol1 = 1.0e-4_real64

    corx = x - floor(x)
    if (corx > (1.0_real64 - tol1)) corx = 0.0_real64
  end function wrap_fractional

  subroutine require_model_initialized(context)
    character(len=*), intent(in) :: context

    if (.not. model_initialized) then
      write(error_unit,'(A,1X,A)') 'Error: PME Hamiltonian is not initialized in', trim(context)
      stop 1
    end if
  end subroutine require_model_initialized

  subroutine validate_level(level, context)
    integer, intent(in) :: level
    character(len=*), intent(in) :: context

    if (level < 0 .or. level > npos) then
      write(error_unit,'(A,1X,A,1X,A,I0,A,I0,A)') 'Error:', trim(context), 'received level', level, &
        ' outside 0..', npos, '.'
      stop 1
    end if
  end subroutine validate_level

  subroutine level_blend_weights(level, total_sites, weight_low, weight_high)
    !  Piecewise power-law hybrid weights controlled by alpha_hybrid and eta_hybrid.
    !
    !  Edge regions (x <= xref or x >= 1-xref) use the corresponding end-member
    !  with full weight, where xref = pme_order_cap / total_sites.
    !
    !  Central region: u = (x-xref)/(1-2*xref), then
    !    w_low  = rho*(1-u)^alpha / (rho*(1-u)^alpha + u^alpha),  rho = exp(eta)
    !    w_high = 1 - w_low
    !
    !  alpha > 1 guarantees continuous first derivatives at the boundaries.
    integer, intent(in) :: level, total_sites
    real(real64), intent(out) :: weight_low, weight_high
    real(real64) :: x, xref, span, u, rho, a0, aN

    if (total_sites <= 0) then
      weight_low = 1.0_real64; weight_high = 0.0_real64
      return
    end if

    x    = real(level, real64) / real(total_sites, real64)
    xref = real(pme_order_cap, real64) / real(total_sites, real64)
    span = 1.0_real64 - 2.0_real64 * xref

    if (span <= 0.0_real64) then
      if (x < 0.5_real64) then
        weight_low = 1.0_real64; weight_high = 0.0_real64
      else if (x > 0.5_real64) then
        weight_low = 0.0_real64; weight_high = 1.0_real64
      else
        weight_low = 0.5_real64; weight_high = 0.5_real64
      end if
    else if (x <= xref) then
      weight_low = 1.0_real64; weight_high = 0.0_real64
    else if (x >= 1.0_real64 - xref) then
      weight_low = 0.0_real64; weight_high = 1.0_real64
    else
      u   = (x - xref) / span
      rho = exp(eta_hybrid)
      a0  = rho * (1.0_real64 - u)**alpha_hybrid
      aN  = u**alpha_hybrid
      weight_low  = a0 / (a0 + aN)
      weight_high = aN / (a0 + aN)
    end if
  end subroutine level_blend_weights

  ! ---------------------------------------------------------------------------
  ! Recalibration helpers
  ! ---------------------------------------------------------------------------

  subroutine fit_epsilon(v0, contributions, unique_count, active_order, &
      ref_indices, ref_energies, ref_count, eps_values, rms_residual, success)
    ! Fits eps_0..eps_K for E(c) = V0 + eps_0 + eps_1*T_1(c) + ... + eps_K*T_K(c).
    ! eps_0 is an additive offset (default 0); eps_1..eps_K are multiplicative (default 1).
    ! Constant predictors (T_1 when all sites are equivalent) are excluded; eps_0 is always
    ! included since it is the intercept term.  If the remaining fit is still unstable
    ! (collinear predictors), the lowest-variance active predictor is progressively pruned.
    ! Pruned predictors keep their default eps value (0 for eps_0, 1 for eps_k, k>=1).
    real(real64), intent(in) :: v0
    real(real64), intent(in) :: contributions(:, :)
    integer, intent(in) :: unique_count, active_order, ref_count
    integer, intent(in) :: ref_indices(:)
    real(real64), intent(in) :: ref_energies(:)
    real(real64), intent(out) :: eps_values(:)
    real(real64), intent(out) :: rms_residual
    logical, intent(out) :: success
    integer :: n_vars, n_active, iref, iord, jord, config_idx, ia, ja, min_ia
    integer :: active_map(max_motif_order+1)
    logical :: is_active(max_motif_order+1)
    real(real64) :: normal(max_motif_order+1, max_motif_order+1)
    real(real64) :: rhs(max_motif_order+1), solution(max_motif_order+1)
    real(real64) :: t_i, t_j, predicted, residual_sum
    real(real64) :: t_mean, t_var, t_var_arr(max_motif_order+1)
    real(real64) :: e_corr(n_calib_configs)
    logical :: stable

    n_vars         = active_order + 1
    eps_values     = 1.0_real64
    eps_values(1)  = 0.0_real64    ! eps_0 additive default
    rms_residual = huge(1.0_real64)
    success      = .false.
    if (ref_count < 1) return
    do iref = 1, ref_count
      if (ref_indices(iref) < 1 .or. ref_indices(iref) > unique_count) return
    end do

    ! Identify predictors with non-negligible variance across calibration configs.
    ! V0 (iord=1) is always included as an intercept: it corrects the absolute energy
    ! offset even though it is constant at a single composition.
    ! Constant T_k (iord>=2, e.g. T1 for symmetric systems) are excluded; they are
    ! collinear with V0 and cannot be independently resolved.
    n_active = 0
    t_var_arr = 0.0_real64
    do iord = 1, n_vars
      t_mean = 0.0_real64
      do iref = 1, ref_count
        config_idx = ref_indices(iref)
        if (iord == 1) then
          t_mean = t_mean + 1.0_real64
        else
          t_mean = t_mean + contributions(iord-1, config_idx)
        end if
      end do
      t_mean = t_mean / real(ref_count, real64)
      t_var = 0.0_real64
      do iref = 1, ref_count
        config_idx = ref_indices(iref)
        if (iord == 1) then
          t_i = 1.0_real64 - t_mean
        else
          t_i = contributions(iord-1, config_idx) - t_mean
        end if
        t_var = t_var + t_i * t_i
      end do
      t_var_arr(iord) = t_var
      ! V0 is always active; T_k (k>=1) are active only if they vary across calib configs
      if (iord == 1 .or. &
          t_var > 1.0e-10_real64 * t_mean**2 * real(ref_count, real64) + 1.0e-30_real64) then
        n_active = n_active + 1
        active_map(n_active) = iord
      end if
    end do

    ! Progressive pruning: fit active set; if unstable, drop lowest-variance predictor.
    ! Inactive predictors keep eps=1; their contribution is subtracted from the target
    ! before building the normal equations (corrected-target formulation).
    do
      if (n_active < 1 .or. ref_count < n_active) return

      ! Mark which variable indices are active
      is_active = .false.
      do ia = 1, n_active
        is_active(active_map(ia)) = .true.
      end do

      ! Corrected target: (E_DFT - V0) minus inactive T_k contributions (eps_k=1)
      ! V0 is always subtracted (fixed, not a parameter); eps_0 is always active so
      ! its default (0) never contributes to the correction.
      do iref = 1, ref_count
        config_idx = ref_indices(iref)
        e_corr(iref) = ref_energies(iref) - v0
        do iord = 2, n_vars
          if (.not. is_active(iord)) &
            e_corr(iref) = e_corr(iref) - contributions(iord-1, config_idx)
        end do
      end do

      ! Build n_active x n_active normal equations against corrected target
      normal   = 0.0_real64
      rhs      = 0.0_real64
      solution = 0.0_real64
      do iref = 1, ref_count
        config_idx = ref_indices(iref)
        do ia = 1, n_active
          iord = active_map(ia)
          if (iord == 1) then
            t_i = 1.0_real64
          else
            t_i = contributions(iord-1, config_idx)
          end if
          rhs(ia) = rhs(ia) + t_i * e_corr(iref)
          do ja = 1, n_active
            jord = active_map(ja)
            if (jord == 1) then
              t_j = 1.0_real64
            else
              t_j = contributions(jord-1, config_idx)
            end if
            normal(ia, ja) = normal(ia, ja) + t_i * t_j
          end do
        end do
      end do

      solution(1:n_active) = rhs(1:n_active)
      call solve_small_linear_system(normal, solution, n_active, success)
      if (.not. success) return

      ! Stability check: |eps| > 100 signals near-collinear predictors
      stable = .true.
      do ia = 1, n_active
        if (abs(solution(ia)) > 1.0e2_real64) then
          stable = .false.
          exit
        end if
      end do
      if (stable) exit

      if (n_active <= 1) then
        success = .false.
        return
      end if

      ! Drop the active predictor with the smallest variance and retry
      min_ia = 1
      do ia = 2, n_active
        if (t_var_arr(active_map(ia)) < t_var_arr(active_map(min_ia))) min_ia = ia
      end do
      active_map(min_ia:n_active-1) = active_map(min_ia+1:n_active)
      n_active = n_active - 1
    end do

    do ia = 1, n_active
      eps_values(active_map(ia)) = solution(ia)
    end do

    residual_sum = 0.0_real64
    do iref = 1, ref_count
      config_idx = ref_indices(iref)
      predicted = v0 + eps_values(1)
      do iord = 1, active_order
        predicted = predicted + eps_values(iord+1) * contributions(iord, config_idx)
      end do
      residual_sum = residual_sum + (predicted - ref_energies(iref))**2
    end do
    rms_residual = sqrt(residual_sum / real(ref_count, real64))
  end subroutine fit_epsilon

  subroutine solve_small_linear_system(matrix, rhs, n, success)
    real(real64), intent(inout) :: matrix(:, :), rhs(:)
    integer, intent(in) :: n
    logical, intent(out) :: success
    real(real64), parameter :: eps = 1.0e-12_real64
    real(real64) :: factor, pivot_abs, temp_val
    real(real64) :: temp_row(n)
    integer :: irow, kcol, pivot_row

    success = .false.
    if (size(matrix, 1) < n .or. size(matrix, 2) < n .or. size(rhs) < n) return

    do kcol = 1, n
      pivot_row = kcol
      pivot_abs = abs(matrix(kcol, kcol))
      do irow = kcol + 1, n
        if (abs(matrix(irow, kcol)) > pivot_abs) then
          pivot_abs = abs(matrix(irow, kcol))
          pivot_row = irow
        end if
      end do
      if (pivot_abs <= eps) return

      if (pivot_row /= kcol) then
        temp_row = matrix(kcol, 1:n)
        matrix(kcol, 1:n) = matrix(pivot_row, 1:n)
        matrix(pivot_row, 1:n) = temp_row
        temp_val = rhs(kcol)
        rhs(kcol) = rhs(pivot_row)
        rhs(pivot_row) = temp_val
      end if

      do irow = kcol + 1, n
        factor = matrix(irow, kcol) / matrix(kcol, kcol)
        matrix(irow, kcol:n) = matrix(irow, kcol:n) - factor * matrix(kcol, kcol:n)
        rhs(irow) = rhs(irow) - factor * rhs(kcol)
      end do
    end do

    do irow = n, 1, -1
      if (abs(matrix(irow, irow)) <= eps) return
      if (irow < n) rhs(irow) = rhs(irow) - dot_product(matrix(irow, irow+1:n), rhs(irow+1:n))
      rhs(irow) = rhs(irow) / matrix(irow, irow)
    end do

    success = .true.
  end subroutine solve_small_linear_system

  subroutine fit_partial_epsilon(v0, contributions, n_calib_fit, full_order, &
      eps_fixed, calib_energies, eps_out, rms_residual, success)
    ! Fits the first n_calib_fit epsilons: eps_0..eps_{n_calib_fit-1}
    ! for the Hamiltonian: E = V0 + eps_0 + sum_{k=1}^{full_order} eps_k*T_k(c)
    ! eps_fixed and eps_out are size full_order+1 (index 1=eps_0 .. full_order+1=eps_{full_order}).
    ! Keeps eps_fixed(n_calib_fit+1..full_order+1) fixed; fits eps_0..eps_{n_calib_fit-1}.
    real(real64), intent(in)  :: v0
    real(real64), intent(in)  :: contributions(:,:)  ! (4, n_configs); first n_calib_fit used
    integer,      intent(in)  :: n_calib_fit, full_order
    real(real64), intent(in)  :: eps_fixed(:)         ! full_order+1 values (1-based: eps_0..eps_K)
    real(real64), intent(in)  :: calib_energies(:)
    real(real64), intent(out) :: eps_out(:)            ! full_order+1 values on output
    real(real64), intent(out) :: rms_residual
    logical,      intent(out) :: success

    integer      :: n_total, iref, iord, jord
    real(real64) :: normal(max_motif_order+1, max_motif_order+1)
    real(real64) :: rhs(max_motif_order+1), solution(max_motif_order+1)
    real(real64) :: fixed_sum, predicted, residual_sum, t_i, t_j

    n_total = full_order + 1
    eps_out      = eps_fixed
    rms_residual = huge(1.0_real64)
    success      = .false.

    if (n_calib_fit < 1 .or. n_total > max_motif_order+1) return
    if (n_calib_fit > n_total) return
    if (n_calib_fit > size(calib_energies) .or. n_calib_fit > size(contributions, 2)) return

    normal = 0.0_real64
    rhs    = 0.0_real64

    do iref = 1, n_calib_fit
      ! Fixed contribution: upper epsilons eps_{n_calib_fit}..eps_{full_order}
      ! in local 1-based: indices n_calib_fit+1..n_total
      fixed_sum = 0.0_real64
      do iord = n_calib_fit+1, n_total
        if (iord == 1) then
          t_i = 1.0_real64
        else
          t_i = contributions(iord-1, iref)
        end if
        fixed_sum = fixed_sum + eps_fixed(iord) * t_i
      end do
      ! Build normal equations for fitted epsilons (indices 1..n_calib_fit)
      do iord = 1, n_calib_fit
        if (iord == 1) then
          t_i = 1.0_real64
        else
          t_i = contributions(iord-1, iref)
        end if
        rhs(iord) = rhs(iord) + t_i * (calib_energies(iref) - v0 - fixed_sum)
        do jord = 1, n_calib_fit
          if (jord == 1) then
            t_j = 1.0_real64
          else
            t_j = contributions(jord-1, iref)
          end if
          normal(iord, jord) = normal(iord, jord) + t_i * t_j
        end do
      end do
    end do

    solution(1:n_calib_fit) = rhs(1:n_calib_fit)
    call solve_small_linear_system(normal, solution, n_calib_fit, success)
    if (.not. success) then
      eps_out = eps_fixed
      return
    end if

    eps_out = eps_fixed
    do iord = 1, n_calib_fit
      eps_out(iord) = solution(iord)
    end do

    residual_sum = 0.0_real64
    do iref = 1, n_calib_fit
      predicted = v0 + eps_out(1)
      do iord = 2, n_total
        predicted = predicted + eps_out(iord) * contributions(iord-1, iref)
      end do
      residual_sum = residual_sum + (predicted - calib_energies(iref))**2
    end do
    rms_residual = sqrt(residual_sum / real(n_calib_fit, real64))
  end subroutine fit_partial_epsilon

  subroutine read_calib_data(target_level, calib_list, n_list, n_use, &
      calib_low_contribs, calib_high_contribs, calib_energies, n_out, ok)
    !  Reads nXX/ENSEMBLE and nXX/ENERGIES, then selects the configurations at
    !  calib_list(1:n_use) for epsilon calibration.
    integer,      intent(in)  :: target_level
    integer,      intent(in)  :: calib_list(:)   ! config indices (1-based) from INPME
    integer,      intent(in)  :: n_list           ! entries available in calib_list
    integer,      intent(in)  :: n_use            ! how many to actually use (≤ n_list)
    real(real64), allocatable, intent(out) :: calib_low_contribs(:,:)   ! (4, n_out)
    real(real64), allocatable, intent(out) :: calib_high_contribs(:,:)  ! (4, n_out)
    real(real64), allocatable, intent(out) :: calib_energies(:)
    integer,      intent(out) :: n_out
    logical,      intent(out) :: ok

    character(len=256) :: energies_path
    character(len=32)  :: level_dir
    integer :: mm_all, i, idx, count, n_missing
    integer, allocatable :: conf_all(:,:), omega_all(:)
    real(real64), allocatable :: energies_all(:)
    logical, allocatable :: e_mask(:)
    real(real64) :: energy, energy_low, energy_high, low_terms(4), high_terms(4)
    logical :: read_ok

    ok    = .false.
    n_out = 0

    if (n_use <= 0 .or. n_list <= 0) return

    call pme_load_level_ensemble(target_level, mm_all, conf_all, omega_all, read_ok)
    if (.not. read_ok .or. mm_all <= 0) then
      write(*,'(A)') ' > Warning: could not read target ENSEMBLE for calibration.'
      return
    end if

    call format_level_directory(target_level, level_dir)
    energies_path = trim(level_dir)//'/ENERGIES'
    allocate(energies_all(mm_all), e_mask(mm_all))
    call read_energies_file(energies_path, mm_all, energies_all, read_ok, n_missing, e_mask)
    if (.not. read_ok) then
      write(*,'(A,A)') ' > Warning: could not read ', trim(energies_path)
      deallocate(conf_all, omega_all, energies_all, e_mask)
      return
    end if

    allocate(calib_low_contribs(4, n_use))
    allocate(calib_high_contribs(4, n_use))
    allocate(calib_energies(n_use))
    count = 0

    do i = 1, n_use
      idx = calib_list(i)
      if (idx < 1 .or. idx > mm_all) then
        write(*,'(A,I0,A,I0,A)') ' > Error: calib_config_list index ', idx, &
          ' out of range (ENSEMBLE has ', mm_all, ' configs).'
        stop 1
      end if
      if (.not. e_mask(idx)) then
        write(*,'(A,I0,A)') ' > Warning: energy missing for calib config ', idx, '; skipping.'
        cycle
      end if
      count = count + 1
      calib_energies(count) = energies_all(idx)
      if (target_level > 0) then
        call pme_evaluate_configuration(conf_all(idx, 1:target_level), target_level, &
          energy, energy_low, energy_high, low_terms, high_terms)
      else
        low_terms  = 0.0_real64
        high_terms = 0.0_real64
      end if
      calib_low_contribs(:, count)  = low_terms
      calib_high_contribs(:, count) = high_terms
    end do

    deallocate(conf_all, omega_all, energies_all, e_mask)
    n_out = count
    ok    = (count > 0)
  end subroutine read_calib_data

  subroutine bisect_calibration_selection(energies, mm, selected_indices, n_selected)
    !  Select up to n_calib_configs (9) configurations using recursive bisection of the
    !  sorted energy list.  Every prefix is representative: the first 3 give min/max/median,
    !  and each additional entry halves the largest remaining gap.
    !  Output order: A, B, C, D, E, F, G, H, I (bisection order, not ascending energy).
    real(real64), intent(in) :: energies(:)
    integer, intent(in) :: mm
    integer, intent(out) :: selected_indices(:)
    integer, intent(out) :: n_selected

    integer, allocatable :: sorted_idx(:)
    logical, allocatable :: seen(:)
    integer :: pos(9), n_pos, i, j, k
    real(real64) :: key_e

    n_selected = 0
    selected_indices = 0
    if (mm <= 0) return

    allocate(sorted_idx(mm), seen(mm))
    seen = .false.

    ! Insertion sort: sorted_idx(i) gives original config index, ascending by energy
    do i = 1, mm
      sorted_idx(i) = i
    end do
    do i = 2, mm
      k = sorted_idx(i)
      key_e = energies(k)
      j = i - 1
      do while (j >= 1 .and. energies(sorted_idx(j)) > key_e)
        sorted_idx(j + 1) = sorted_idx(j)
        j = j - 1
      end do
      sorted_idx(j + 1) = k
    end do

    ! Bisection positions (1-based into sorted_idx):
    ! A=1, B=mm, C=mid(1,mm), D=mid(1,C), E=mid(C,mm),
    ! F=mid(1,D), G=mid(C,E), H=mid(D,C), I=mid(E,mm)
    n_pos = min(n_calib_configs, mm)
    if (mm >= 1) pos(1) = 1
    if (mm >= 2) pos(2) = mm
    if (mm >= 3) pos(3) = (1 + mm) / 2
    if (mm >= 4) pos(4) = (1 + pos(3)) / 2
    if (mm >= 5) pos(5) = (pos(3) + mm) / 2
    if (mm >= 6) pos(6) = (1 + pos(4)) / 2
    if (mm >= 7) pos(7) = (pos(3) + pos(5)) / 2
    if (mm >= 8) pos(8) = (pos(4) + pos(3)) / 2
    if (mm >= 9) pos(9) = (pos(5) + mm) / 2

    do i = 1, n_pos
      j = pos(i)
      if (j < 1 .or. j > mm) cycle
      if (seen(j)) cycle
      seen(j) = .true.
      n_selected = n_selected + 1
      selected_indices(n_selected) = sorted_idx(j)
    end do

    deallocate(sorted_idx, seen)
  end subroutine bisect_calibration_selection

  subroutine write_pme_model_template(filename, all_low_terms, all_high_terms, mm, target_level)
    !  Writes SODPROJECT/INPME.tmp with the new 7-data-line format:
    !    PME choice, PMEorder, n_calib, calib_config_list, eps_low, eps_high, alpha eta.
    !  No calib/ subdirectory is created; calibration configs are generated via
    !  sod_gener.sh -choose and energies are read from nXX/ENERGIES.
    character(len=*), intent(in) :: filename
    real(real64), intent(in) :: all_low_terms(:, :), all_high_terms(:, :)
    integer, intent(in) :: mm, target_level

    integer :: unit_tmp, idx, n_shared, variant_choice, order_cap_out
    integer :: shared_idx(n_calib_configs)
    character(len=32) :: level_dir
    real(real64), allocatable :: selection_energy(:)
    real(real64) :: e_low_tmp, e_high_tmp

    call format_level_directory(target_level, level_dir)

    allocate(selection_energy(mm))
    do idx = 1, mm
      call compute_pme_variant_energy(pme_choice, target_level, all_low_terms(:, idx), all_high_terms(:, idx), &
        selection_energy(idx), e_low_tmp, e_high_tmp)
    end do
    call bisect_calibration_selection(selection_energy, mm, shared_idx, n_shared)
    deallocate(selection_energy)

    variant_choice = pme_choice
    if (.not. high_base_loaded) variant_choice = 0
    order_cap_out = min(pme_order_cap, max(max_low_order, max_high_order))
    order_cap_out = max(2, min(max_motif_order, order_cap_out))

    open(newunit=unit_tmp, file=trim(filename), status='replace', action='write')
    write(unit_tmp,'(A)') '# SOD PME calibration/control file.'
    write(unit_tmp,'(A,A,A)') '# Rename this file to ', trim(pme_model_filename), ' to enable calibration.'
    write(unit_tmp,'(A,A,A)') '# Target level directory: ', trim(level_dir), '/'
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PME Hamiltonian: 0 (PME0 low-side); 1 (PME1 high-side); 2 (PMEh hybrid)'
    write(unit_tmp,'(I0)') variant_choice
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PMEorder cap (2, 3, or 4)'
    write(unit_tmp,'(I0)') order_cap_out
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# n_calib: 0 no calibration (manual epsilons); 1..9 number of configurations for calibration'
    write(unit_tmp,'(I0)') n_shared
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# calib_config_list: configuration indices; only first n_calib used'
    write(unit_tmp,'(*(1X,I0))') shared_idx(1:n_shared)
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# epsilon_low (epsilon_0 ... epsilon_PMEorder)'
    write(unit_tmp,'(*(1X,F10.6))') 0.0_real64, eps_low(1:order_cap_out)
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# epsilon_high (epsilon_N ... epsilon_{N-PMEorder})'
    write(unit_tmp,'(*(1X,F10.6))') 0.0_real64, eps_high(1:order_cap_out)
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PMEh_alpha (>1.0 for smooth blending, default 2.0) and PMEh_eta (default 0.0; >0 favours low-end, <0 favours high-end)'
    write(unit_tmp,'(F10.6,1X,F10.6)') 2.0_real64, 0.0_real64
    write(unit_tmp,'(A)') ''
    call write_v_block_section(unit_tmp, variant_choice)
    close(unit_tmp)
  end subroutine write_pme_model_template

  ! Writes INPME.tmp with n_calib=0 and the fitted eps/alpha/eta from module state,
  ! plus the V block.  Called automatically after a successful calibration run.
  subroutine write_pme_model_calibrated(filename)
    character(len=*), intent(in) :: filename
    integer :: unit_tmp, eff_low_ord, eff_high_ord
    eff_low_ord  = min(max_low_order,  pme_order_cap)
    eff_high_ord = min(max_high_order, pme_order_cap)
    open(newunit=unit_tmp, file=trim(filename), status='replace', action='write')
    write(unit_tmp,'(A)') '# SOD PME calibration/control file (calibrated — n_calib=0).'
    write(unit_tmp,'(A,A,A)') '# Rename to ', trim(pme_model_filename), ' to use these parameters directly.'
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PME Hamiltonian: 0 (PME0 low-side); 1 (PME1 high-side); 2 (PMEh hybrid)'
    write(unit_tmp,'(I0)') pme_choice
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PMEorder cap (2, 3, or 4)'
    write(unit_tmp,'(I0)') pme_order_cap
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# n_calib: 0 no calibration (manual epsilons); 1..9 number of configurations for calibration'
    write(unit_tmp,'(I0)') 0
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# calib_config_list: configuration indices (kept for reference; not used when n_calib=0)'
    if (n_calib_list > 0) then
      write(unit_tmp,'(*(1X,I0))') calib_config_list(1:n_calib_list)
    else
      write(unit_tmp,'(A)') ' 0'
    end if
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# epsilon_low (epsilon_0 ... epsilon_PMEorder)'
    write(unit_tmp,'(*(1X,F20.14))') eps_low(0:eff_low_ord)
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# epsilon_high (epsilon_N ... epsilon_{N-PMEorder})'
    write(unit_tmp,'(*(1X,F20.14))') eps_high(0:eff_high_ord)
    write(unit_tmp,'(A)') ''
    write(unit_tmp,'(A)') '# PMEh_alpha (>1.0 for smooth blending, default 2.0) and PMEh_eta (default 0.0; >0 favours low-end, <0 favours high-end)'
    write(unit_tmp,'(F20.14,1X,F20.14)') alpha_hybrid, eta_hybrid
    write(unit_tmp,'(A)') ''
    call write_v_block_section(unit_tmp, pme_choice)
    close(unit_tmp)
  end subroutine write_pme_model_calibrated

  ! Writes the "# PME energy terms" V-block section to an already-open unit.
  ! Uses the compact V residual arrays currently in module state.
  ! 3 columns (level m v) for PME0/PME1; 4 columns (level m vl vh) for PMEh.
  subroutine write_v_block_section(unit_out, pme_ch)
    integer, intent(in) :: unit_out, pme_ch
    integer :: m

    if (nv1_low < 1 .and. nv1_high < 1) return

    write(unit_out,'(A)') '# PME energy terms'
    write(unit_out,'(A)') ''
    if (pme_ch == 2) then
      write(unit_out,'(I2,1X,I2,2(1X,F22.14))') 0, 1, v0_low, v0_high
    else if (pme_ch == 0) then
      write(unit_out,'(I2,1X,I2,1X,F22.14)') 0, 1, v0_low
    else
      write(unit_out,'(I2,1X,I2,1X,F22.14)') 0, 1, v0_high
    end if

    do m = 1, max(nv1_low, nv1_high)
      if (pme_ch == 2) then
        block
          real(real64) :: vl, vh
          if (m <= nv1_low)  then; vl = v1res_low(m);  else; vl = 0.0_real64; end if
          if (m <= nv1_high) then; vh = v1res_high(m); else; vh = 0.0_real64; end if
          write(unit_out,'(I2,1X,I2,2(1X,F22.14))') 1, m, vl, vh
        end block
      else if (pme_ch == 0 .and. m <= nv1_low) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 1, m, v1res_low(m)
      else if (pme_ch == 1 .and. m <= nv1_high) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 1, m, v1res_high(m)
      end if
    end do

    do m = 1, max(nv2_low, nv2_high)
      if (pme_ch == 2) then
        block
          real(real64) :: vl, vh
          if (m <= nv2_low)  then; vl = v2res_low(m);  else; vl = 0.0_real64; end if
          if (m <= nv2_high) then; vh = v2res_high(m); else; vh = 0.0_real64; end if
          write(unit_out,'(I2,1X,I2,2(1X,F22.14))') 2, m, vl, vh
        end block
      else if (pme_ch == 0 .and. m <= nv2_low) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 2, m, v2res_low(m)
      else if (pme_ch == 1 .and. m <= nv2_high) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 2, m, v2res_high(m)
      end if
    end do

    do m = 1, max(nv3_low, nv3_high)
      if (pme_ch == 2) then
        block
          real(real64) :: vl, vh
          if (m <= nv3_low)  then; vl = v3res_low(m);  else; vl = 0.0_real64; end if
          if (m <= nv3_high) then; vh = v3res_high(m); else; vh = 0.0_real64; end if
          write(unit_out,'(I2,1X,I2,2(1X,F22.14))') 3, m, vl, vh
        end block
      else if (pme_ch == 0 .and. m <= nv3_low) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 3, m, v3res_low(m)
      else if (pme_ch == 1 .and. m <= nv3_high) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 3, m, v3res_high(m)
      end if
    end do

    do m = 1, max(nv4_low, nv4_high)
      if (pme_ch == 2) then
        block
          real(real64) :: vl, vh
          if (m <= nv4_low)  then; vl = v4res_low(m);  else; vl = 0.0_real64; end if
          if (m <= nv4_high) then; vh = v4res_high(m); else; vh = 0.0_real64; end if
          write(unit_out,'(I2,1X,I2,2(1X,F22.14))') 4, m, vl, vh
        end block
      else if (pme_ch == 0 .and. m <= nv4_low) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 4, m, v4res_low(m)
      else if (pme_ch == 1 .and. m <= nv4_high) then
        write(unit_out,'(I2,1X,I2,1X,F22.14)') 4, m, v4res_high(m)
      end if
    end do
  end subroutine write_v_block_section

  subroutine read_pme_model(filename, pme_ch, order_cap, n_calib, &
      calib_list_out, n_calib_list_out, alpha_out, eta_out, eps_low_out, eps_high_out, ok)
    !  INPME format (7 data lines; comment lines starting with '#' and blank lines are skipped):
    !    data line 1: PME choice          (integer: 0=PME0, 1=PME1, 2=PMEh)
    !    data line 2: PMEOrder cap        (integer: 2, 3, or 4; other values are a fatal error)
    !    data line 3: n_calib             (integer: 0=no calib; 1..9=fit n configs; other values are fatal)
    !    data line 4: calib_config_list   (space-separated integers; only first n_calib used)
    !    data line 5: eps_low values      (order_cap+1 reals: epsilon_0 .. epsilon_K)
    !    data line 6: eps_high values     (order_cap+1 reals: epsilon_0 .. epsilon_K)
    !    data line 7: PMEh_alpha PMEh_eta  (two reals on one line; applied only for PMEh)
    character(len=*), intent(in) :: filename
    integer,      intent(out) :: pme_ch, order_cap, n_calib
    integer,      intent(out) :: calib_list_out(9)
    integer,      intent(out) :: n_calib_list_out
    real(real64), intent(out) :: alpha_out, eta_out
    real(real64), intent(out) :: eps_low_out(0:4), eps_high_out(0:4)
    logical,      intent(out) :: ok

    integer :: unit_in, ios, i, n_calib_list_tmp, tmp_list(9)
    character(len=512) :: line
    logical :: lk

    ok               = .false.
    pme_ch           = 2
    order_cap        = max_motif_order
    n_calib          = 0
    calib_list_out   = 0
    n_calib_list_out = 0
    alpha_out        = 2.0_real64
    eta_out          = 0.0_real64
    eps_low_out(0)   = 0.0_real64; eps_low_out(1:)  = 1.0_real64
    eps_high_out(0)  = 0.0_real64; eps_high_out(1:) = 1.0_real64

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    ! Line 1: PME choice
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) pme_ch
    if (ios /= 0) then; close(unit_in); return; end if
    if (pme_ch < 0 .or. pme_ch > 2) then; close(unit_in); return; end if

    ! Line 2: PMEOrder cap (must be 2, 3, or 4)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) order_cap
    if (ios /= 0) then; close(unit_in); return; end if
    if (order_cap < 2 .or. order_cap > 4) then
      write(*,'(A,A,A,I0,A)') ' > Error: ', trim(pme_model_filename), ' PMEorder = ', order_cap, '; must be 2, 3, or 4.'
      stop 1
    end if

    ! Line 3: n_calib (0..9 only)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) n_calib
    if (ios /= 0) then; close(unit_in); return; end if
    if (n_calib < 0 .or. n_calib > 9) then
      write(*,'(A,A,A,I0,A)') ' > Error: ', trim(pme_model_filename), ' n_calib = ', n_calib, '; must be 0..9.'
      stop 1
    end if

    ! Line 4: calib_config_list (space-separated integers; try to read up to 9)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    n_calib_list_tmp = 0
    do i = 1, 9
      read(line, *, iostat=ios) tmp_list(1:i)
      if (ios /= 0) exit
      n_calib_list_tmp = i
    end do
    if (n_calib_list_tmp == 0) then; close(unit_in); return; end if
    if (n_calib > 0 .and. n_calib_list_tmp < n_calib) then
      write(*,'(A,I0,A,I0,A)') ' > Warning: calib_config_list has ', n_calib_list_tmp, &
        ' entries but n_calib = ', n_calib, '; reducing n_calib to match.'
      n_calib = n_calib_list_tmp
    end if
    if (n_calib > 0 .and. n_calib_list_tmp > n_calib) then
      write(*,'(A,I0,A,I0,A)') ' > Note: calib_config_list has ', n_calib_list_tmp, &
        ' entries; only first ', n_calib, ' will be used.'
    end if
    calib_list_out(1:n_calib_list_tmp) = tmp_list(1:n_calib_list_tmp)
    n_calib_list_out = n_calib_list_tmp

    ! Line 5: eps_low (order_cap+1 reals: epsilon_0 .. epsilon_K)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) eps_low_out(0:order_cap)
    if (ios /= 0) then; close(unit_in); return; end if

    ! Line 6: eps_high (order_cap+1 reals: epsilon_0 .. epsilon_K)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) eps_high_out(0:order_cap)
    if (ios /= 0) then; close(unit_in); return; end if

    ! Line 7: PMEh_alpha and PMEh_eta (two reals on one line)
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) alpha_out, eta_out
    if (ios /= 0) then; close(unit_in); return; end if
    if (alpha_out <= 0.0_real64) then
      write(*,'(A,F10.4,A)') ' > Error: PMEh_alpha = ', alpha_out, ' must be positive.'
      stop 1
    end if
    if (alpha_out <= 1.0_real64) then
      write(*,'(A,F10.4,A)') ' > Warning: PMEh_alpha = ', alpha_out, &
        '; alpha > 1 is required for continuous first derivatives at the blending boundaries.'
    end if

    close(unit_in)
    ok = .true.
  end subroutine read_pme_model

  ! Reads the optional V-block section from INPME (lines after the 7 mandatory lines).
  ! Populates module-level compact V residual arrays and v0_low / v0_high.
  ! Format: one line per reference config — "level m v_low v_high" (PMEh) or
  !         "level m v_side" (PME0 or PME1, 3 columns).  pme_choice on line 1
  !         determines the column count.
  subroutine read_pme_model_v_block(filename, ok)
    character(len=*), intent(in) :: filename
    logical, intent(out) :: ok

    integer :: unit_in, ios, i, lev, idx, pme_ch_local
    character(len=512) :: line
    logical :: lk
    real(real64) :: vl, vh
    real(real64), allocatable :: tmp1(:), tmp2(:), tmp3(:), tmp4(:)
    real(real64), allocatable :: tmp1h(:), tmp2h(:), tmp3h(:), tmp4h(:)
    integer :: n1l, n2l, n3l, n4l, n1h, n2h, n3h, n4h
    integer, parameter :: max_m = 9999

    ok = .false.
    n1l = 0;  n2l = 0;  n3l = 0;  n4l = 0
    n1h = 0;  n2h = 0;  n3h = 0;  n4h = 0
    allocate(tmp1(max_m), tmp2(max_m), tmp3(max_m), tmp4(max_m))
    allocate(tmp1h(max_m), tmp2h(max_m), tmp3h(max_m), tmp4h(max_m))

    open(newunit=unit_in, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      deallocate(tmp1, tmp2, tmp3, tmp4, tmp1h, tmp2h, tmp3h, tmp4h)
      return
    end if

    ! Read line 1 to get pme_choice; skip lines 2–7
    call next_data_line(unit_in, line, lk)
    if (.not. lk) then; close(unit_in); return; end if
    read(line, *, iostat=ios) pme_ch_local
    if (ios /= 0 .or. pme_ch_local < 0 .or. pme_ch_local > 2) then
      close(unit_in)
      deallocate(tmp1, tmp2, tmp3, tmp4, tmp1h, tmp2h, tmp3h, tmp4h)
      return
    end if
    do i = 2, 7
      call next_data_line(unit_in, line, lk)
      if (.not. lk) then
        close(unit_in)
        deallocate(tmp1, tmp2, tmp3, tmp4, tmp1h, tmp2h, tmp3h, tmp4h)
        return
      end if
    end do

    ! Read V-block lines (data lines 8+)
    do
      call next_data_line(unit_in, line, lk)
      if (.not. lk) exit
      vl = 0.0_real64;  vh = 0.0_real64
      if (pme_ch_local == 2) then
        read(line, *, iostat=ios) lev, idx, vl, vh
        if (ios /= 0) cycle
      else
        read(line, *, iostat=ios) lev, idx, vl
        if (ios /= 0) cycle
        if (pme_ch_local == 1) then
          vh = vl;  vl = 0.0_real64
        end if
      end if
      if (lev < 0 .or. lev > 4 .or. idx < 1 .or. idx > max_m) cycle
      if (lev == 0 .and. idx == 1) then
        if (pme_ch_local == 0 .or. pme_ch_local == 2) v0_low  = vl
        if (pme_ch_local == 1 .or. pme_ch_local == 2) v0_high = vh
      else if (lev == 1) then
        if (pme_ch_local == 0 .or. pme_ch_local == 2) then; tmp1(idx)  = vl;  n1l = max(n1l, idx);  end if
        if (pme_ch_local == 1 .or. pme_ch_local == 2) then; tmp1h(idx) = vh;  n1h = max(n1h, idx);  end if
      else if (lev == 2) then
        if (pme_ch_local == 0 .or. pme_ch_local == 2) then; tmp2(idx)  = vl;  n2l = max(n2l, idx);  end if
        if (pme_ch_local == 1 .or. pme_ch_local == 2) then; tmp2h(idx) = vh;  n2h = max(n2h, idx);  end if
      else if (lev == 3) then
        if (pme_ch_local == 0 .or. pme_ch_local == 2) then; tmp3(idx)  = vl;  n3l = max(n3l, idx);  end if
        if (pme_ch_local == 1 .or. pme_ch_local == 2) then; tmp3h(idx) = vh;  n3h = max(n3h, idx);  end if
      else if (lev == 4) then
        if (pme_ch_local == 0 .or. pme_ch_local == 2) then; tmp4(idx)  = vl;  n4l = max(n4l, idx);  end if
        if (pme_ch_local == 1 .or. pme_ch_local == 2) then; tmp4h(idx) = vh;  n4h = max(n4h, idx);  end if
      end if
    end do
    close(unit_in)

    ! Transfer into module compact arrays if any V data was found
    if (n1l > 0 .or. n1h > 0) then
      ok = .true.
      if (n1l > 0) then; allocate(v1res_low(n1l));   v1res_low  = tmp1(1:n1l);   nv1_low  = n1l;  end if
      if (n2l > 0) then; allocate(v2res_low(n2l));   v2res_low  = tmp2(1:n2l);   nv2_low  = n2l;  end if
      if (n3l > 0) then; allocate(v3res_low(n3l));   v3res_low  = tmp3(1:n3l);   nv3_low  = n3l;  end if
      if (n4l > 0) then; allocate(v4res_low(n4l));   v4res_low  = tmp4(1:n4l);   nv4_low  = n4l;  end if
      if (n1h > 0) then; allocate(v1res_high(n1h));  v1res_high = tmp1h(1:n1h);  nv1_high = n1h;  end if
      if (n2h > 0) then; allocate(v2res_high(n2h));  v2res_high = tmp2h(1:n2h);  nv2_high = n2h;  end if
      if (n3h > 0) then; allocate(v3res_high(n3h));  v3res_high = tmp3h(1:n3h);  nv3_high = n3h;  end if
      if (n4h > 0) then; allocate(v4res_high(n4h));  v4res_high = tmp4h(1:n4h);  nv4_high = n4h;  end if
    end if
    deallocate(tmp1, tmp2, tmp3, tmp4, tmp1h, tmp2h, tmp3h, tmp4h)
  end subroutine read_pme_model_v_block

  ! Reconstructs V1_low..V4_low + max_low_order from the compact cache.
  ! Reads n0k/ENSEMBLE for site patterns; skips ENERGIES entirely.
  subroutine reconstruct_v_low_from_cache(ok)
    logical, intent(out) :: ok

    integer :: mm, m, op, i1, i2, i3, i4
    integer, allocatable :: conf(:, :), omega(:)
    logical :: ens_ok

    ok = .false.
    if (nv1_low < 1) return

    call read_level_conf_only(1, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv1_low) then
      write(*,'(A,I0,A,I0,A)') ' > Warning: n01/ENSEMBLE has ', mm, &
        ' configs but V cache has ', nv1_low, '; falling back to reference ENERGIES.'
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(V1_low(npos))
    V1_low = 0.0_real64
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        V1_low(i1) = v1res_low(m)
      end do
    end do
    max_low_order = 1
    deallocate(conf, omega)

    if (nv2_low < 1) then; ok = .true.; return; end if
    call read_level_conf_only(2, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv2_low) then
      write(*,'(A)') ' > Warning: n02/ENSEMBLE mismatch with V cache; falling back to reference ENERGIES.'
      deallocate(V1_low)
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(V2_low(npos, npos))
    V2_low = 0.0_real64
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        i2 = eqmatrix(op, conf(m, 2))
        call sort_pair(i1, i2)
        call assign_pair_term(V2_low, i1, i2, v2res_low(m))
      end do
    end do
    max_low_order = 2
    deallocate(conf, omega)

    if (nv3_low < 1) then; ok = .true.; return; end if
    call read_level_conf_only(3, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv3_low) then
      write(*,'(A)') ' > Warning: n03/ENSEMBLE mismatch with V cache; using order-2 low model.'
      ok = .true.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        i2 = eqmatrix(op, conf(m, 2))
        i3 = eqmatrix(op, conf(m, 3))
        call sort_triplet(i1, i2, i3)
        call assign_triplet_term(V3_i_low, V3_j_low, V3_k_low, V3_val_low, i1, i2, i3, v3res_low(m))
      end do
    end do
    max_low_order = 3
    deallocate(conf, omega)

    if (nv4_low < 1) then
      call prepare_triplet_lookup(V3_i_low, V3_j_low, V3_k_low, V3_val_low, V3_key_low)
      call prepare_quad_lookup(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, V4_key_low)
      ok = .true.
      return
    end if
    call prepare_triplet_lookup(V3_i_low, V3_j_low, V3_k_low, V3_val_low, V3_key_low)
    call read_level_conf_only(4, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv4_low) then
      write(*,'(A)') ' > Warning: n04/ENSEMBLE mismatch with V cache; using order-3 low model.'
      call prepare_quad_lookup(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, V4_key_low)
      ok = .true.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, conf(m, 1))
        i2 = eqmatrix(op, conf(m, 2))
        i3 = eqmatrix(op, conf(m, 3))
        i4 = eqmatrix(op, conf(m, 4))
        call sort_quad(i1, i2, i3, i4)
        call assign_quad_term(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, i1, i2, i3, i4, v4res_low(m))
      end do
    end do
    max_low_order = 4
    deallocate(conf, omega)
    call prepare_triplet_lookup(V3_i_low, V3_j_low, V3_k_low, V3_val_low, V3_key_low)
    call prepare_quad_lookup(V4_i_low, V4_j_low, V4_k_low, V4_l_low, V4_val_low, V4_key_low)
    ok = .true.
  end subroutine reconstruct_v_low_from_cache

  ! Reconstructs V1_high..V4_high + max_high_order from the compact cache.
  subroutine reconstruct_v_high_from_cache(ok)
    logical, intent(out) :: ok

    integer :: mm, m, op, i1, i2, i3, i4
    integer, allocatable :: conf(:, :), omega(:), holes(:, :)
    logical :: ens_ok

    ok = .false.
    if (nv1_high < 1) return

    high_base_loaded = .true.
    allocate(V1_high(npos), V2_high(npos, npos))
    V1_high = 0.0_real64;  V2_high = 0.0_real64
    max_high_order = 0

    call read_level_conf_only(npos - 1, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv1_high) then
      write(*,'(A)') ' > Warning: high-side n-1/ENSEMBLE mismatch with V cache; falling back to reference ENERGIES.'
      if (allocated(V1_high)) deallocate(V1_high)
      if (allocated(V2_high)) deallocate(V2_high)
      high_base_loaded = .false.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(holes(mm, 1))
    call build_hole_configurations(conf, mm, npos - 1, 1, holes)
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, holes(m, 1))
        V1_high(i1) = v1res_high(m)
      end do
    end do
    max_high_order = 1
    deallocate(conf, omega, holes)

    if (nv2_high < 1) then; ok = .true.; return; end if
    call read_level_conf_only(npos - 2, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv2_high) then
      write(*,'(A)') ' > Warning: high-side n-2/ENSEMBLE mismatch; using order-1 high model.'
      ok = .true.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(holes(mm, 2))
    call build_hole_configurations(conf, mm, npos - 2, 2, holes)
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, holes(m, 1))
        i2 = eqmatrix(op, holes(m, 2))
        call sort_pair(i1, i2)
        call assign_pair_term(V2_high, i1, i2, v2res_high(m))
      end do
    end do
    max_high_order = 2
    deallocate(conf, omega, holes)

    if (nv3_high < 1) then; ok = .true.; return; end if
    call read_level_conf_only(npos - 3, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv3_high) then
      write(*,'(A)') ' > Warning: high-side n-3/ENSEMBLE mismatch; using order-2 high model.'
      ok = .true.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(holes(mm, 3))
    call build_hole_configurations(conf, mm, npos - 3, 3, holes)
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, holes(m, 1))
        i2 = eqmatrix(op, holes(m, 2))
        i3 = eqmatrix(op, holes(m, 3))
        call sort_triplet(i1, i2, i3)
        call assign_triplet_term(V3_i_high, V3_j_high, V3_k_high, V3_val_high, i1, i2, i3, v3res_high(m))
      end do
    end do
    max_high_order = 3
    deallocate(conf, omega, holes)

    if (nv4_high < 1) then
      call prepare_triplet_lookup(V3_i_high, V3_j_high, V3_k_high, V3_val_high, V3_key_high)
      call prepare_quad_lookup(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, V4_key_high)
      ok = .true.
      return
    end if
    call prepare_triplet_lookup(V3_i_high, V3_j_high, V3_k_high, V3_val_high, V3_key_high)
    call read_level_conf_only(npos - 4, mm, conf, omega, ens_ok)
    if (.not. ens_ok .or. mm /= nv4_high) then
      write(*,'(A)') ' > Warning: high-side n-4/ENSEMBLE mismatch; using order-3 high model.'
      call prepare_quad_lookup(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, V4_key_high)
      ok = .true.
      if (allocated(conf))  deallocate(conf)
      if (allocated(omega)) deallocate(omega)
      return
    end if
    allocate(holes(mm, 4))
    call build_hole_configurations(conf, mm, npos - 4, 4, holes)
    do m = 1, mm
      do op = 1, nop
        i1 = eqmatrix(op, holes(m, 1))
        i2 = eqmatrix(op, holes(m, 2))
        i3 = eqmatrix(op, holes(m, 3))
        i4 = eqmatrix(op, holes(m, 4))
        call sort_quad(i1, i2, i3, i4)
        call assign_quad_term(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, i1, i2, i3, i4, v4res_high(m))
      end do
    end do
    max_high_order = 4
    deallocate(conf, omega, holes)
    call prepare_triplet_lookup(V3_i_high, V3_j_high, V3_k_high, V3_val_high, V3_key_high)
    call prepare_quad_lookup(V4_i_high, V4_j_high, V4_k_high, V4_l_high, V4_val_high, V4_key_high)
    ok = .true.
  end subroutine reconstruct_v_high_from_cache

  subroutine copy_text_file(source_path, dest_path, ok)
    character(len=*), intent(in) :: source_path, dest_path
    logical, intent(out) :: ok

    integer :: unit_in, unit_out, ios
    character(len=1024) :: line

    ok = .false.
    open(newunit=unit_in, file=trim(source_path), status='old', action='read', iostat=ios)
    if (ios /= 0) return

    open(newunit=unit_out, file=trim(dest_path), status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      close(unit_in)
      return
    end if

    do
      read(unit_in, '(A)', iostat=ios) line
      if (ios /= 0) exit
      write(unit_out, '(A)') trim(line)
    end do

    close(unit_in)
    close(unit_out)
    ok = .true.
  end subroutine copy_text_file

  ! ---------------------------------------------------------------------------
  ! FILER / structure-file helpers for INPME.tmp embedding
  ! ---------------------------------------------------------------------------

  subroutine read_filer_from_insod(filer_value, ok)
    !  Scans INSOD for the '# FILER' comment marker and reads the integer on
    !  the next non-comment, non-blank data line.
    integer, intent(out) :: filer_value
    logical, intent(out) :: ok
    integer :: unit_in, ios
    character(len=512) :: line
    logical :: found_filer_comment

    ok = .false.
    filer_value = -1
    found_filer_comment = .false.

    if (insod_cache_valid) then
      filer_value = insod_filer_cache
      ok = .true.
      return
    end if

    open(newunit=unit_in, file='INSOD', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    do
      read(unit_in, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (found_filer_comment) then
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read(line, *, iostat=ios) filer_value
        if (ios == 0) ok = .true.
        exit
      end if
      ! Detect the FILER comment line (starts with '#' and contains 'FILER')
      if (line(1:1) == '#' .and. index(line, 'FILER') > 0) then
        found_filer_comment = .true.
      end if
    end do

    close(unit_in)
  end subroutine read_filer_from_insod

  subroutine structure_file_for_filer(filer, fname)
    !  Returns the structure-file name for a given FILER code.
    !  Returns an empty string for FILER types without a single primary file.
    integer, intent(in) :: filer
    character(len=*), intent(out) :: fname

    select case (filer)
    case (1)
      fname = 'input.gin'
    case (11)
      fname = 'POSCAR'
    case (12)
      fname = 'castep.cell'
    case (13)
      fname = 'pw.in'
    case default
      fname = ''
    end select
  end subroutine structure_file_for_filer

  ! ---------------------------------------------------------------------------
  ! INSOD species-info reader
  ! ---------------------------------------------------------------------------

  subroutine read_insod_species_info(nsp_out, symbols_out, newsym_new, newsym_old, has_special, ok)
    !  Reads nsp, species symbols, and newsymbol line from INSOD by sequential parsing.
    !  newsym_new = first token of newsymbol line (the new/substituting species).
    !  newsym_old = second token (the remaining species on the target site).
    !  has_special = .true. if either newsymbol starts with '@' (molecule) or '%' (vacancy).
    integer, intent(out) :: nsp_out
    character(len=4), allocatable, intent(out) :: symbols_out(:)
    character(len=5), intent(out) :: newsym_new, newsym_old
    logical, intent(out) :: has_special, ok

    ok          = .false.
    nsp_out     = 0
    newsym_new  = ''
    newsym_old  = ''
    has_special = .false.

    if (insod_cache_valid) then
      nsp_out    = insod_nsp_cache
      newsym_new = insod_newsym_new_cache
      newsym_old = insod_newsym_old_cache
      has_special = insod_has_special_cache
      if (allocated(insod_symbols_cache)) then
        allocate(symbols_out(insod_nsp_cache))
        symbols_out = insod_symbols_cache
      end if
      ok = .true.
      return
    end if

    block
      type(insod_t) :: d
      call read_insod('INSOD', d)
      nsp_out = d%nsp
      allocate(symbols_out(nsp_out))
      symbols_out(1:nsp_out) = d%symbol
      newsym_new = d%newsymbol(1,1)
      newsym_old = d%newsymbol(1,2)
    end block
    has_special = (newsym_new(1:1) == '@' .or. newsym_new(1:1) == '%' .or. &
                   newsym_old(1:1) == '@' .or. newsym_old(1:1) == '%')
    ok = .true.
  end subroutine read_insod_species_info

  ! ---------------------------------------------------------------------------
  ! Supercell CIF reader
  ! ---------------------------------------------------------------------------

  subroutine read_supercell_cif_data(nsp, symbols, nat_out, coords_out, spat_out, &
      a, b, c, al, be, ga, ok)
    !  Reads supercell.cif: cell parameters + fractional coords + species assignment.
    !  On success allocates coords_out(nat_out,3) and spat_out(nat_out).
    integer, intent(in) :: nsp
    character(len=4), intent(in) :: symbols(:)
    integer, intent(out) :: nat_out
    real(real64), allocatable, intent(out) :: coords_out(:, :)
    integer, allocatable, intent(out) :: spat_out(:)
    real(real64), intent(out) :: a, b, c, al, be, ga
    logical, intent(out) :: ok

    integer, parameter :: nat_max = 30000
    real(real64), allocatable :: coords_buf(:, :)
    integer, allocatable :: spat_buf(:)
    integer :: unit_in, ios, nat, sp
    character(len=200) :: cifline
    character(len=20) :: atmlabel
    character(len=4) :: atmsymbol
    logical :: in_atom_loop

    ok  = .false.
    nat = 0
    a   = 1.0_real64;  b  = 1.0_real64;  c  = 1.0_real64
    al  = 90.0_real64; be = 90.0_real64; ga = 90.0_real64
    in_atom_loop = .false.

    allocate(coords_buf(nat_max, 3), spat_buf(nat_max))
    open(newunit=unit_in, file='supercell.cif', status='old', action='read', iostat=ios)
    if (ios /= 0) then; deallocate(coords_buf, spat_buf); return; end if

    do
      read(unit_in, '(A)', iostat=ios) cifline
      if (ios /= 0) exit
      cifline = adjustl(cifline)
      if      (cifline(1:14) == '_cell_length_a') then
        read(cifline(15:), *, iostat=ios) a
      else if (cifline(1:14) == '_cell_length_b') then
        read(cifline(15:), *, iostat=ios) b
      else if (cifline(1:14) == '_cell_length_c') then
        read(cifline(15:), *, iostat=ios) c
      else if (cifline(1:17) == '_cell_angle_alpha') then
        read(cifline(18:), *, iostat=ios) al
      else if (cifline(1:16) == '_cell_angle_beta') then
        read(cifline(17:), *, iostat=ios) be
      else if (cifline(1:17) == '_cell_angle_gamma') then
        read(cifline(18:), *, iostat=ios) ga
      else if (cifline(1:18) == '_atom_site_fract_z') then
        in_atom_loop = .true.
      else if (in_atom_loop .and. len_trim(cifline) > 0 .and. &
               cifline(1:1) /= '_' .and. cifline(1:1) /= '#' .and. &
               cifline(1:4) /= 'loop' .and. cifline(1:4) /= 'data' .and. &
               cifline(1:1) /= "'") then
        nat = nat + 1
        if (nat > nat_max) then
          write(error_unit,'(A,I0,A)') ' Warning: supercell.cif has more than ', nat_max, &
            ' atoms; truncating structure read.'
          nat = nat - 1
          exit
        end if
        read(cifline, *, iostat=ios) atmlabel, atmsymbol, &
          coords_buf(nat,1), coords_buf(nat,2), coords_buf(nat,3)
        if (ios /= 0) then
          nat = nat - 1
          cycle
        end if
        spat_buf(nat) = 1
        do sp = 1, nsp
          if (trim(atmsymbol) == trim(symbols(sp))) then
            spat_buf(nat) = sp
            exit
          end if
        end do
      end if
    end do
    close(unit_in)

    if (nat == 0) then; deallocate(coords_buf, spat_buf); return; end if
    nat_out = nat
    allocate(coords_out(nat_out, 3), spat_out(nat_out))
    coords_out = coords_buf(1:nat_out, :)
    spat_out   = spat_buf(1:nat_out)
    deallocate(coords_buf, spat_buf)
    ok = .true.
  end subroutine read_supercell_cif_data

  ! ---------------------------------------------------------------------------
  ! Lattice-vector utility
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  ! Structure-file dispatcher (writers are in structwriters module)
  ! ---------------------------------------------------------------------------

  subroutine write_config_structure(confdir, conf_local, level_val, nat, coords, spat, &
      symbols, a, b, c, al, be, ga, sptarget_val, newsym_new, newsym_old, &
      filer_val, config_idx, ok)
    !  Builds sym_arr(nat) from substitution info and delegates to structwriters.
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: conf_local(:), level_val, nat, sptarget_val, filer_val, config_idx
    integer, intent(in) :: spat(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    character(len=4), intent(in) :: symbols(:)
    character(len=5), intent(in) :: newsym_new, newsym_old
    logical, intent(out) :: ok

    character(len=5), allocatable :: sym_arr(:)
    integer :: at, sp, local_idx

    ok = .true.
    if (filer_val < 0 .or. filer_val == 2) return  ! FILER=-1 or LAMMPS: no structure file

    ! Build resolved symbol array
    allocate(sym_arr(nat))
    do at = 1, nat
      sp = spat(at)
      if (sp /= sptarget_val) then
        sym_arr(at) = trim(symbols(sp))
      else
        local_idx = at - target_atini + 1
        if (any(conf_local(1:level_val) == local_idx)) then
          sym_arr(at) = trim(newsym_new)
        else
          sym_arr(at) = trim(newsym_old)
        end if
      end if
    end do

    select case (filer_val)
    case (0)
      call sw_write_cif_p1(confdir, nat, sym_arr, coords, a, b, c, al, be, ga, ok)
    case (1)
      call sw_write_gulp(confdir, nat, sym_arr, coords, a, b, c, al, be, ga, config_idx, ok)
    case (11)
      call sw_write_poscar(confdir, nat, sym_arr, coords, a, b, c, al, be, ga, ok)
    case (12)
      call sw_write_castep(confdir, nat, sym_arr, coords, a, b, c, al, be, ga, config_idx, ok)
    case (13)
      call sw_write_qe(confdir, nat, sym_arr, coords, a, b, c, al, be, ga, config_idx, ok)
    end select

    deallocate(sym_arr)
  end subroutine write_config_structure

  ! Old structure-file writers (write_config_cif_p1, write_config_gulp,
  ! write_config_poscar, write_config_castep, write_config_qe, lat_vectors)
  ! have been moved to the structwriters module (src/structwriters.f90).

  ! ---------------------------------------------------------------------------
  !  Canonicalization: map a configuration to its lexicographically minimal
  !  form under all symmetry operations in eqmatrix(nop, npos).
  ! ---------------------------------------------------------------------------

  subroutine canonicalize_config(conf, level, canonical)
    integer, intent(in)  :: conf(:), level
    integer, intent(out) :: canonical(level)
    integer :: op, i, mapped(level)
    logical :: is_less

    canonical = conf(1:level)  ! start with identity
    do op = 1, nop
      ! Apply operation: map each site index
      do i = 1, level
        mapped(i) = eqmatrix(op, conf(i))
      end do
      call mc_insertion_sort(mapped, level)
      ! Compare lexicographically: is mapped < canonical?
      is_less = .false.
      do i = 1, level
        if (mapped(i) < canonical(i)) then
          is_less = .true.; exit
        else if (mapped(i) > canonical(i)) then
          exit
        end if
      end do
      if (is_less) canonical = mapped
    end do
  end subroutine canonicalize_config

  ! ---------------------------------------------------------------------------
  !  Compute symmetry degeneracy Ω = nop / |Stab(C)| for a canonical config.
  !  |Stab(C)| is the number of symmetry operations that leave C invariant.
  ! ---------------------------------------------------------------------------

  subroutine compute_config_degeneracy(canonical, level, omega)
    integer, intent(in)  :: canonical(:), level
    integer, intent(out) :: omega
    integer :: op, i, stab_count
    integer :: mapped(level)

    stab_count = 0
    do op = 1, nop
      do i = 1, level
        mapped(i) = eqmatrix(op, canonical(i))
      end do
      call mc_insertion_sort(mapped, level)
      if (all(mapped(1:level) == canonical(1:level))) stab_count = stab_count + 1
    end do
    omega = nop / stab_count
  end subroutine compute_config_degeneracy

  ! ---------------------------------------------------------------------------
  !  Apply epsilon renormalisation (§3.2) and alpha/eta power-law hybrid weighting
  !  to the per-order Hamiltonian contributions → calibrated energy.
  !  When epsilon = 1.0, alpha = 2.0 and eta = 0.0 (defaults), result is identical to
  !  raw pme_evaluate_configuration with symmetric power-law blending.
  ! ---------------------------------------------------------------------------

  subroutine apply_epsilon_energy(level_val, low_terms, high_terms, energy)
    integer, intent(in) :: level_val
    real(real64), intent(in) :: low_terms(4), high_terms(4)
    real(real64), intent(out) :: energy
    real(real64) :: e_low_cal, e_high_cal, weight_low, weight_high
    integer :: eff_low_ord, eff_high_ord

    eff_low_ord  = min(max_low_order,  pme_order_cap)
    eff_high_ord = min(max_high_order, pme_order_cap)

    if (eff_low_ord >= 1) then
      e_low_cal = v0_low + eps_low(0) + dot_product(eps_low(1:eff_low_ord), low_terms(1:eff_low_ord))
    else
      e_low_cal = v0_low + eps_low(0)
    end if

    select case (pme_choice)
    case (0)
      energy = e_low_cal
    case (1)
      if (high_base_loaded .and. eff_high_ord >= 1) then
        e_high_cal = v0_high + eps_high(0) + dot_product(eps_high(1:eff_high_ord), high_terms(1:eff_high_ord))
        energy = e_high_cal
      else
        energy = e_low_cal
      end if
    case default  ! 2 = PMEh
      if (high_base_loaded .and. eff_high_ord >= 1) then
        e_high_cal = v0_high + eps_high(0) + dot_product(eps_high(1:eff_high_ord), high_terms(1:eff_high_ord))
        call level_blend_weights(level_val, npos, weight_low, weight_high)
        energy = weight_low * e_low_cal + weight_high * e_high_cal
      else
        energy = e_low_cal
      end if
    end select
  end subroutine apply_epsilon_energy

  subroutine apply_epsilon_energy_decomp(level_val, low_terms, high_terms, e_blend, e_low, e_high)
    !  Apply epsilon renormalisation and return all three energies for diagnostics.
    integer, intent(in) :: level_val
    real(real64), intent(in) :: low_terms(4), high_terms(4)
    real(real64), intent(out) :: e_blend, e_low, e_high
    real(real64) :: weight_low, weight_high
    integer :: eff_low_ord, eff_high_ord

    eff_low_ord  = min(max_low_order,  pme_order_cap)
    eff_high_ord = min(max_high_order, pme_order_cap)

    if (eff_low_ord >= 1) then
      e_low = v0_low + eps_low(0) + dot_product(eps_low(1:eff_low_ord), low_terms(1:eff_low_ord))
    else
      e_low = v0_low + eps_low(0)
    end if

    if (high_base_loaded .and. eff_high_ord >= 1) then
      e_high = v0_high + eps_high(0) + dot_product(eps_high(1:eff_high_ord), high_terms(1:eff_high_ord))
    else
      e_high = e_low
    end if

    select case (pme_choice)
    case (0)
      e_blend = e_low
    case (1)
      e_blend = e_high
    case default  ! 2 = PMEh
      if (high_base_loaded) then
        call level_blend_weights(level_val, npos, weight_low, weight_high)
        e_blend = weight_low * e_low + weight_high * e_high
      else
        e_blend = e_low
      end if
    end select
  end subroutine apply_epsilon_energy_decomp

  subroutine check_target_vs_reference(target_level, pme_ch, eff_low_ord, eff_high_ord)
    integer, intent(in) :: target_level, pme_ch, eff_low_ord, eff_high_ord
    logical :: low_conflict, high_conflict
    low_conflict  = (pme_ch == 0 .or. pme_ch == 2) .and. (target_level <= eff_low_ord)
    high_conflict = (pme_ch == 1 .or. pme_ch == 2) .and. high_base_loaded .and. &
                    (target_level >= npos - eff_high_ord)
    if (low_conflict) then
      write(error_unit,'(A,I0,A,I0,A)') &
        'Error: target n', target_level, ' is a low-side reference level (n00..n', eff_low_ord, &
        '). PME cannot predict its own training data.'
      stop 1
    end if
    if (high_conflict) then
      write(error_unit,'(A,I0,A,I0,A,I0,A)') &
        'Error: target n', target_level, ' is a high-side reference level (n', npos - eff_high_ord, &
        '..n', npos, '). PME cannot predict its own training data.'
      stop 1
    end if
  end subroutine check_target_vs_reference

  subroutine note_target_energies_if_present(tgt)
    integer, intent(in) :: tgt
    character(len=32) :: level_dir
    character(len=300) :: energies_path
    logical :: ex
    if (tgt <= 0) return
    call format_level_directory(tgt, level_dir)
    energies_path = trim(level_dir)//'/ENERGIES'
    inquire(file=trim(energies_path), exist=ex)
    if (ex) write(*,'(A,A,A)') &
      ' Note: Reference ENERGIES found for the target level (', trim(energies_path), &
      '); PME predictions can be benchmarked against them.'
  end subroutine note_target_energies_if_present

  ! ---------------------------------------------------------------------------
  !  Pre-load recalibration from SODPROJECT/INPME and nXX/ENERGIES.
  !  Sets module vars: pme_choice, pme_order_cap, eps_low, eps_high, alpha_hybrid, eta_hybrid.
  !  Must be called BEFORE the MC walk or energy output so that Metropolis
  !  acceptance uses calibrated Hamiltonian energies.
  !  If INPME is absent, sensible defaults apply (PMEh if high side available,
  !  PME0 otherwise; eps=1, alpha=2, eta=0).
  ! ---------------------------------------------------------------------------

  subroutine pme_preload_recalibration(target_level, ok)
    integer, intent(in)  :: target_level
    logical, intent(out) :: ok

    character(len=32)  :: level_dir
    character(len=256) :: model_path, model_snapshot_path
    logical  :: have_model, fit_ok, fit_ok_loo, copy_ok
    integer  :: pme_ch, order_cap, n_calib_req, eff_low_ord, eff_high_ord, n_use, n_calib_read
    integer  :: mkdir_status
    integer  :: iref, jloo, loo_n
    real(real64) :: alpha_read, eta_read, rms_low, rms_high, rms_tmp
    real(real64) :: eps_low_read(0:4), eps_high_read(0:4), eps_loo(0:4)
    integer  :: calib_list_read(9), n_calib_list_read
    integer  :: loo_idx(n_calib_configs)
    real(real64), allocatable :: calib_low_contribs(:,:), calib_high_contribs(:,:), calib_energies(:)
    integer, allocatable :: seq_idx(:)
    real(real64) :: rms_before_low, rms_before_high, rms_before_pmeh
    real(real64) :: rms_loo_low, rms_loo_high, rms_after_pmeh
    real(real64) :: weight_low_calib, weight_high_calib
    real(real64) :: e_low_bef, e_high_bef, e_low_aft, e_high_aft, e_loo_low, e_loo_high, resid_sum
    logical :: low_fit_ok, high_fit_ok, use_loo_low, use_loo_high, use_loo_pmeh

    ok = .true.
    eps_low(0)   = 0.0_real64; eps_low(1:4)  = 1.0_real64
    eps_high(0)  = 0.0_real64; eps_high(1:4) = 1.0_real64
    alpha_hybrid = 2.0_real64
    eta_hybrid   = 0.0_real64
    pme_choice    = 2
    pme_order_cap = max_motif_order
    pme_calibration_applied = .false.

    call format_level_directory(target_level, level_dir)
    model_path = trim(pme_model_filename)
    model_snapshot_path = trim(level_dir)//'/'//trim(pme_model_filename)
    inquire(file=trim(model_path), exist=have_model)

    if (.not. have_model) then
      if (high_base_loaded) then
        pme_choice = 2
        write(*,'(A,A,A)') ' > No ', trim(pme_model_filename), ' found; using PMEh with eps=1, alpha=2, eta=0.'
      else
        pme_choice = 0
        write(*,'(A,A,A)') ' > No ', trim(pme_model_filename), ' found, no high-side data: using PME0 with eps=1.'
      end if
      eff_low_ord   = min(max_low_order,  pme_order_cap)
      eff_high_ord  = min(max_high_order, pme_order_cap)
      write(*,'(A,I0,A,I0,A)') ' > Effective order applied: low=', eff_low_ord, ', high=', eff_high_ord, '.'
      call check_target_vs_reference(target_level, pme_choice, eff_low_ord, eff_high_ord)
      return
    end if

    call execute_command_line('mkdir -p '//trim(level_dir), exitstat=mkdir_status)
    call copy_text_file(model_path, model_snapshot_path, copy_ok)
    if (copy_ok) then
      write(*,'(A,A)') ' > Model snapshot written to ', trim(model_snapshot_path)
    else
      write(*,'(A,A)') ' > Warning: could not write model snapshot to ', trim(model_snapshot_path)
    end if

    call read_pme_model(model_path, pme_ch, order_cap, n_calib_req, &
      calib_list_read, n_calib_list_read, alpha_read, eta_read, eps_low_read, eps_high_read, ok)
    if (.not. ok) then
      write(*,'(A,A,A)') ' > Warning: could not parse ', trim(pme_model_filename), '; using defaults (PMEh, eps=1, alpha=2, eta=0).'
      ok = .true.
      return
    end if

    calib_config_list = calib_list_read
    n_calib_list      = n_calib_list_read

    pme_choice    = pme_ch
    pme_order_cap = order_cap
    eff_low_ord   = min(max_low_order,  pme_order_cap)
    eff_high_ord  = min(max_high_order, pme_order_cap)

    if (eff_low_ord < pme_order_cap .and. (pme_ch == 0 .or. pme_ch == 2)) then
      write(*,'(A,I0,A,I0,A)') ' > Warning: PMEOrder cap ', pme_order_cap, &
        ' exceeds available low-side orders ', eff_low_ord, '; capping.'
    end if
    if (high_base_loaded .and. eff_high_ord < pme_order_cap .and. (pme_ch == 1 .or. pme_ch == 2)) then
      write(*,'(A,I0,A,I0,A)') ' > Warning: PMEOrder cap ', pme_order_cap, &
        ' exceeds available high-side orders ', eff_high_ord, '; capping.'
    end if
    write(*,'(A,I0,A,I0,A)') ' > Effective order applied: low=', eff_low_ord, ', high=', eff_high_ord, '.'
    call check_target_vs_reference(target_level, pme_choice, eff_low_ord, eff_high_ord)

    if (n_calib_req == 0) then
      if (pme_ch == 0 .or. pme_ch == 2) eps_low(0:eff_low_ord)  = eps_low_read(0:eff_low_ord)
      if (pme_ch == 1 .or. pme_ch == 2) then
        if (high_base_loaded) eps_high(0:eff_high_ord) = eps_high_read(0:eff_high_ord)
      end if
      if (pme_ch == 2) then
        alpha_hybrid = alpha_read
        eta_hybrid   = eta_read
      end if
      pme_calibration_applied = .true.
      write(*,'(A,A,A)', advance='no') ' > eps read from ', trim(pme_model_filename), ' (n_calib=0).'
      if (pme_ch == 0 .or. pme_ch == 2) then
        write(*,'(A)', advance='no') '  ε_low ='
        do iref = 0, eff_low_ord
          write(*,'(1X,F9.5)', advance='no') eps_low(iref)
        end do
      end if
      if ((pme_ch == 1 .or. pme_ch == 2) .and. high_base_loaded) then
        write(*,'(A)', advance='no') '  ε_high ='
        do iref = 0, eff_high_ord
          write(*,'(1X,F9.5)', advance='no') eps_high(iref)
        end do
      end if
      write(*,*)
      return
    end if

    ! n_calib_req >= 1: fit from nXX/ENSEMBLE + nXX/ENERGIES using calib_config_list
    call read_calib_data(target_level, calib_config_list, n_calib_list, n_calib_req, &
      calib_low_contribs, calib_high_contribs, calib_energies, n_calib_read, fit_ok)
    if (.not. fit_ok .or. n_calib_read == 0) then
      write(*,'(A)') ' > Warning: calib data not found or all energies missing; using eps=1.'
      return
    end if

    pme_calibration_applied = .true.   ! pme.model was read; calibration attempted
    n_use = n_calib_read
    n_calib_used = n_use

    ! Sequential indices 1..n_use into the calib contribs array
    allocate(seq_idx(n_use))
    do iref = 1, n_use
      seq_idx(iref) = iref
    end do

    ! Pre-calibration RMSE (eps=1 throughout)
    call level_blend_weights(target_level, npos, weight_low_calib, weight_high_calib)
    resid_sum = 0.0_real64
    do iref = 1, n_use
      e_low_bef = v0_low + sum(calib_low_contribs(1:eff_low_ord, iref))
      resid_sum = resid_sum + (e_low_bef - calib_energies(iref))**2
    end do
    rms_before_low = sqrt(resid_sum / real(n_use, real64))

    resid_sum = 0.0_real64
    do iref = 1, n_use
      e_high_bef = v0_high + sum(calib_high_contribs(1:eff_high_ord, iref))
      resid_sum  = resid_sum + (e_high_bef - calib_energies(iref))**2
    end do
    rms_before_high = sqrt(resid_sum / real(n_use, real64))

    resid_sum = 0.0_real64
    do iref = 1, n_use
      e_low_bef  = v0_low  + sum(calib_low_contribs(1:eff_low_ord, iref))
      e_high_bef = v0_high + sum(calib_high_contribs(1:eff_high_ord, iref))
      resid_sum  = resid_sum + &
        (weight_low_calib*e_low_bef + weight_high_calib*e_high_bef - calib_energies(iref))**2
    end do
    rms_before_pmeh = sqrt(resid_sum / real(n_use, real64))

    low_fit_ok  = .false.
    high_fit_ok = .false.

    ! --- Low-side ---
    use_loo_low = .false.
    if ((pme_ch == 0 .or. pme_ch == 2) .and. eff_low_ord >= 1) then
      if (n_use <= eff_low_ord) then
        ! Underdetermined: keep upper eps from INPME, fit lower n_use epsilons (eps_0..eps_{n_use-1})
        eps_low(0:eff_low_ord) = eps_low_read(0:eff_low_ord)
        call fit_partial_epsilon(v0_low, calib_low_contribs, n_use, eff_low_ord, &
          eps_low_read(0:eff_low_ord), calib_energies, eps_low(0:eff_low_ord), rms_low, fit_ok)
      else
        ! Overdetermined or square LS: fit all eps_0..eps_{eff_low_ord}
        call fit_epsilon(v0_low, calib_low_contribs, n_calib_read, eff_low_ord, &
          seq_idx(1:n_use), calib_energies, n_use, eps_low(0:eff_low_ord), rms_low, fit_ok)
      end if
      if (fit_ok) then
        low_fit_ok = .true.
        ! LOO-CV: each fold has n_use-1 configs; requires n_use-1 > eff_low_ord
        if (n_use >= eff_low_ord + 2) then
          resid_sum  = 0.0_real64
          use_loo_low = .true.
          do jloo = 1, n_use
            loo_n = 0
            do iref = 1, n_use
              if (iref /= jloo) then; loo_n = loo_n + 1; loo_idx(loo_n) = iref; end if
            end do
            call fit_epsilon(v0_low, calib_low_contribs, n_calib_read, eff_low_ord, &
              loo_idx(1:loo_n), calib_energies, loo_n, eps_loo, rms_tmp, fit_ok_loo)
            if (.not. fit_ok_loo) then; use_loo_low = .false.; exit; end if
            e_loo_low = v0_low + eps_loo(0) + &
              dot_product(eps_loo(1:eff_low_ord), calib_low_contribs(1:eff_low_ord, jloo))
            resid_sum = resid_sum + (e_loo_low - calib_energies(jloo))**2
          end do
          if (use_loo_low) rms_loo_low = sqrt(resid_sum / real(n_use, real64))
        end if
        rms_calib_low = merge(rms_loo_low, rms_low, use_loo_low)
        if (use_loo_low) then
          write(*,'(A,I0,A,I0,A,F9.5,A,F9.5,A)', advance='no') &
            ' > Low-side eps (order=', eff_low_ord, ', ', n_use, ' configs, RMSE ', &
            rms_before_low, ' -> LOO-CV ', rms_loo_low, ' eV):'
        else
          write(*,'(A,I0,A,I0,A,F9.5,A,F9.5,A)', advance='no') &
            ' > Low-side eps (order=', eff_low_ord, ', ', n_use, ' configs, RMSE ', &
            rms_before_low, ' -> ', rms_low, ' eV):'
        end if
        do iref = 0, eff_low_ord
          write(*,'(1X,F9.5)', advance='no') eps_low(iref)
        end do
        write(*,*)
      else
        eps_low(0) = 0.0_real64; eps_low(1:4) = 1.0_real64
        write(*,'(A)') ' > Warning: low-side eps calibration unsuccessful; keeping eps=default.'
      end if
    end if

    ! --- High-side ---
    use_loo_high = .false.
    if (high_base_loaded .and. (pme_ch == 1 .or. pme_ch == 2) .and. eff_high_ord >= 1) then
      if (n_use <= eff_high_ord) then
        eps_high(0:eff_high_ord) = eps_high_read(0:eff_high_ord)
        call fit_partial_epsilon(v0_high, calib_high_contribs, n_use, eff_high_ord, &
          eps_high_read(0:eff_high_ord), calib_energies, eps_high(0:eff_high_ord), rms_high, fit_ok)
      else
        call fit_epsilon(v0_high, calib_high_contribs, n_calib_read, eff_high_ord, &
          seq_idx(1:n_use), calib_energies, n_use, eps_high(0:eff_high_ord), rms_high, fit_ok)
      end if
      if (fit_ok) then
        high_fit_ok = .true.
        ! LOO-CV: each fold has n_use-1 configs; requires n_use-1 > eff_high_ord
        if (n_use >= eff_high_ord + 2) then
          resid_sum   = 0.0_real64
          use_loo_high = .true.
          do jloo = 1, n_use
            loo_n = 0
            do iref = 1, n_use
              if (iref /= jloo) then; loo_n = loo_n + 1; loo_idx(loo_n) = iref; end if
            end do
            call fit_epsilon(v0_high, calib_high_contribs, n_calib_read, eff_high_ord, &
              loo_idx(1:loo_n), calib_energies, loo_n, eps_loo, rms_tmp, fit_ok_loo)
            if (.not. fit_ok_loo) then; use_loo_high = .false.; exit; end if
            e_loo_high = v0_high + eps_loo(0) + &
              dot_product(eps_loo(1:eff_high_ord), calib_high_contribs(1:eff_high_ord, jloo))
            resid_sum = resid_sum + (e_loo_high - calib_energies(jloo))**2
          end do
          if (use_loo_high) rms_loo_high = sqrt(resid_sum / real(n_use, real64))
        end if
        rms_calib_high = merge(rms_loo_high, rms_high, use_loo_high)
        if (use_loo_high) then
          write(*,'(A,I0,A,I0,A,F9.5,A,F9.5,A)', advance='no') &
            ' > High-side eps (order=', eff_high_ord, ', ', n_use, ' configs, RMSE ', &
            rms_before_high, ' -> LOO-CV ', rms_loo_high, ' eV):'
        else
          write(*,'(A,I0,A,I0,A,F9.5,A,F9.5,A)', advance='no') &
            ' > High-side eps (order=', eff_high_ord, ', ', n_use, ' configs, RMSE ', &
            rms_before_high, ' -> ', rms_high, ' eV):'
        end if
        do iref = 0, eff_high_ord
          write(*,'(1X,F9.5)', advance='no') eps_high(iref)
        end do
        write(*,*)
      else
        eps_high(0) = 0.0_real64; eps_high(1:4) = 1.0_real64
        write(*,'(A)') ' > Warning: high-side eps calibration unsuccessful; keeping eps=default.'
      end if
    end if

    ! --- PMEh combined score (after both sides fitted) ---
    if (pme_ch == 2 .and. high_base_loaded .and. (low_fit_ok .or. high_fit_ok)) then
      ! LOO-CV for PMEh: re-fit both sides per fold; requires both to be LOO-eligible
      use_loo_pmeh = low_fit_ok .and. high_fit_ok .and. &
                     n_use >= max(eff_low_ord, eff_high_ord) + 2
      if (use_loo_pmeh) then
        resid_sum = 0.0_real64
        do jloo = 1, n_use
          loo_n = 0
          do iref = 1, n_use
            if (iref /= jloo) then; loo_n = loo_n + 1; loo_idx(loo_n) = iref; end if
          end do
          call fit_epsilon(v0_low, calib_low_contribs, n_calib_read, eff_low_ord, &
            loo_idx(1:loo_n), calib_energies, loo_n, eps_loo, rms_tmp, fit_ok_loo)
          if (.not. fit_ok_loo) then; use_loo_pmeh = .false.; exit; end if
          e_loo_low = v0_low + eps_loo(0) + &
            dot_product(eps_loo(1:eff_low_ord), calib_low_contribs(1:eff_low_ord, jloo))
          call fit_epsilon(v0_high, calib_high_contribs, n_calib_read, eff_high_ord, &
            loo_idx(1:loo_n), calib_energies, loo_n, eps_loo, rms_tmp, fit_ok_loo)
          if (.not. fit_ok_loo) then; use_loo_pmeh = .false.; exit; end if
          e_loo_high = v0_high + eps_loo(0) + &
            dot_product(eps_loo(1:eff_high_ord), calib_high_contribs(1:eff_high_ord, jloo))
          resid_sum = resid_sum + &
            (weight_low_calib*e_loo_low + weight_high_calib*e_loo_high - calib_energies(jloo))**2
        end do
      end if
      if (use_loo_pmeh) then
        rms_after_pmeh = sqrt(resid_sum / real(n_use, real64))
        write(*,'(A,I0,A,F9.5,A,F9.5,A)') ' > PMEh LOO-CV (', n_use, ' configs): ', &
          rms_before_pmeh, ' -> ', rms_after_pmeh, ' eV'
      else
        resid_sum = 0.0_real64
        do iref = 1, n_use
          e_low_aft  = v0_low  + eps_low(0)  + &
            dot_product(eps_low(1:eff_low_ord),  calib_low_contribs(1:eff_low_ord, iref))
          e_high_aft = v0_high + eps_high(0) + &
            dot_product(eps_high(1:eff_high_ord), calib_high_contribs(1:eff_high_ord, iref))
          resid_sum  = resid_sum + &
            (weight_low_calib*e_low_aft + weight_high_calib*e_high_aft - calib_energies(iref))**2
        end do
        rms_after_pmeh = sqrt(resid_sum / real(n_use, real64))
        write(*,'(A,I0,A,F9.5,A,F9.5,A)') ' > PMEh RMSE (', n_use, ' configs): ', &
          rms_before_pmeh, ' -> ', rms_after_pmeh, ' eV'
      end if
    end if

    ! --- Alpha/eta (PMEh only, manual) ---
    if (pme_ch == 2) then
      alpha_hybrid = alpha_read
      eta_hybrid   = eta_read
    end if

    if (allocated(calib_low_contribs))  deallocate(calib_low_contribs)
    if (allocated(calib_high_contribs)) deallocate(calib_high_contribs)
    if (allocated(calib_energies))      deallocate(calib_energies)
    if (allocated(seq_idx))             deallocate(seq_idx)
  end subroutine pme_preload_recalibration

  subroutine mc_insertion_sort(arr, n)
    !  In-place insertion sort for a small integer array.
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: n
    integer :: i, j, key

    do i = 2, n
      key = arr(i)
      j   = i - 1
      do while (j >= 1)
        if (arr(j) <= key) exit
        arr(j+1) = arr(j)
        j = j - 1
      end do
      arr(j+1) = key
    end do
  end subroutine mc_insertion_sort

end module pmemod

!*******************************************************************************
!    Structure-file writers shared by genersod and pmemod.
!
!    Each writer receives a fully resolved configuration: sym_arr(nat) contains
!    the final atomic symbols (CIF names for non-GULP, CIF names that the GULP
!    writer maps via sod_type_map internally).
!
!    Part of the SOD package — GNU GPL v3+.
!*******************************************************************************

module structwriters
  use iso_fortran_env, only: real64, error_unit
  implicit none
  private

  public :: sw_write_gulp, sw_write_poscar, sw_write_cif_p1
  public :: sw_write_castep, sw_write_qe
  public :: sw_lat_vectors

  integer, parameter :: sw_max_tmpl = 100000
  integer, parameter :: sw_max_sym  = 100
  integer, parameter :: sw_max_map  = 50

contains

  ! ---------------------------------------------------------------------------
  !  Lattice vectors from cell parameters
  ! ---------------------------------------------------------------------------

  subroutine sw_lat_vectors(a, b, c, al, be, ga, lv)
    !  Converts a,b,c (angstrom) and alpha,beta,gamma (degrees) to 3x3 Cartesian matrix.
    !  lv(i,j): i=Cartesian component, j=vector (1=a, 2=b, 3=c).
    real(real64), intent(in) :: a, b, c, al, be, ga
    real(real64), intent(out) :: lv(3, 3)
    real(real64), parameter :: pi = acos(-1.0_real64)
    real(real64), parameter :: d2r = pi / 180.0_real64
    real(real64) :: cosa, cosb, cosg, sing, trm1

    cosa = merge(0.0_real64, cos(al*d2r), al == 90.0_real64)
    cosb = merge(0.0_real64, cos(be*d2r), be == 90.0_real64)
    if (ga == 90.0_real64) then
      cosg = 0.0_real64;  sing = 1.0_real64
    else
      cosg = cos(ga*d2r); sing = sin(ga*d2r)
    end if

    lv(1,1) = a;       lv(2,1) = 0.0_real64; lv(3,1) = 0.0_real64
    lv(1,2) = b*cosg;  lv(2,2) = b*sing;      lv(3,2) = 0.0_real64
    lv(1,3) = c*cosb
    lv(2,3) = c*(cosa - cosg*cosb)/sing
    trm1    = lv(2,3)/c
    lv(3,3) = c*sqrt(1.0_real64 - cosb**2 - trm1**2)
  end subroutine sw_lat_vectors

  ! ---------------------------------------------------------------------------
  !  CIF P1 writer (FILER=0)
  ! ---------------------------------------------------------------------------

  subroutine sw_write_cif_p1(confdir, nat, sym_arr, coords, &
      a, b, c, al, be, ga, ok)
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: nat
    character(len=*), intent(in) :: sym_arr(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    logical, intent(out) :: ok

    integer :: unit_out, at, ios
    character(len=16) :: lbl

    ok = .false.
    open(newunit=unit_out, file=trim(confdir)//'/configuration.cif', &
      status='replace', action='write', iostat=ios)
    if (ios /= 0) return

    write(unit_out,'(A)') 'data_SOD_configuration'
    write(unit_out,'(A,F12.6)') '_cell_length_a   ', a
    write(unit_out,'(A,F12.6)') '_cell_length_b   ', b
    write(unit_out,'(A,F12.6)') '_cell_length_c   ', c
    write(unit_out,'(A,F12.6)') '_cell_angle_alpha ', al
    write(unit_out,'(A,F12.6)') '_cell_angle_beta  ', be
    write(unit_out,'(A,F12.6)') '_cell_angle_gamma ', ga
    write(unit_out,'(A)') "_symmetry_space_group_name_H-M   'P 1'"
    write(unit_out,'(A)') 'loop_'
    write(unit_out,'(A)') '_atom_site_label'
    write(unit_out,'(A)') '_atom_site_type_symbol'
    write(unit_out,'(A)') '_atom_site_fract_x'
    write(unit_out,'(A)') '_atom_site_fract_y'
    write(unit_out,'(A)') '_atom_site_fract_z'

    do at = 1, nat
      write(lbl,'(A,I0)') trim(sym_arr(at)), at
      write(unit_out,'(A,2X,A,2X,3(F10.7,2X))') trim(lbl), trim(sym_arr(at)), &
        coords(at,1), coords(at,2), coords(at,3)
    end do

    close(unit_out)
    ok = .true.
  end subroutine sw_write_cif_p1

  ! ---------------------------------------------------------------------------
  !  GULP writer (FILER=1)
  ! ---------------------------------------------------------------------------

  subroutine sw_write_gulp(confdir, nat, sym_arr, coords, &
      a, b, c, al, be, ga, config_idx, ok)
    !  Reads template_input.gin, applies sod_type_map, detects shells (from
    !  template + library), replaces markers, writes confdir/input.gin.
    !  Copies any referenced library file to confdir/.
    !  sym_arr contains CIF-level symbols; mapping to GULP names is done here.
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: nat, config_idx
    character(len=*), intent(in) :: sym_arr(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    logical, intent(out) :: ok

    integer, parameter :: len_marker = 22  ! len('@configuration_number@')
    character(len=200), allocatable :: tmpl(:)
    integer :: unit_in, unit_out, ios, ntmpl, struc_idx, l, at, is, idx_tok
    integer :: nmap, im, ios_cp
    character(len=200) :: lbuf, outbuf_gulp, lib_file
    character(len=10) :: tok1, tok2, map_from(sw_max_map), map_to(sw_max_map)
    character(len=16) :: numstr
    logical :: ishell(sw_max_sym)
    character(len=5) :: all_sym(sw_max_sym), sym_mapped
    integer :: nall
    logical :: found_s

    ok = .false.

    ! --- Read template ---
    open(newunit=unit_in, file='template_input.gin', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    allocate(tmpl(sw_max_tmpl))
    ntmpl = 0; struc_idx = 0
    do
      read(unit_in, '(A)', iostat=ios) lbuf
      if (ios /= 0) exit
      ntmpl = ntmpl + 1
      if (ntmpl > sw_max_tmpl) then; close(unit_in); deallocate(tmpl); return; end if
      tmpl(ntmpl) = lbuf
      if (trim(adjustl(lbuf)) == '@configuration_structure@' .and. struc_idx == 0) struc_idx = ntmpl
    end do
    close(unit_in)
    if (struc_idx == 0) then; deallocate(tmpl); return; end if

    ! --- Parse sod_type_map directives ---
    nmap = 0
    do l = 1, ntmpl
      lbuf = adjustl(tmpl(l))
      if (lbuf(1:15) == '# sod_type_map ') then
        if (nmap < sw_max_map) then
          nmap = nmap + 1
          read(lbuf(16:), *, iostat=ios) map_from(nmap), map_to(nmap)
          if (ios /= 0) nmap = nmap - 1
        end if
      end if
    end do

    ! --- Build unique symbol list with GULP names ---
    nall = 0
    do at = 1, nat
      sym_mapped = trim(sym_arr(at))
      do im = 1, nmap
        if (trim(sym_mapped) == trim(map_from(im))) then
          sym_mapped = trim(map_to(im)); exit
        end if
      end do
      found_s = .false.
      do is = 1, nall
        if (trim(sym_mapped) == trim(all_sym(is))) then; found_s = .true.; exit; end if
      end do
      if (.not. found_s .and. nall < sw_max_sym) then
        nall = nall + 1
        all_sym(nall) = sym_mapped
      end if
    end do

    ! --- Detect shell species: scan template ---
    ishell = .false.
    do l = 1, ntmpl
      lbuf = adjustl(tmpl(l))
      if (len_trim(lbuf) == 0 .or. lbuf(1:1) == '#') cycle
      tok1 = ''; tok2 = ''
      read(lbuf, *, iostat=ios) tok1, tok2
      if (ios == 0 .and. trim(tok2) == 'shel') then
        do is = 1, nall
          if (trim(tok1) == trim(all_sym(is))) then; ishell(is) = .true.; exit; end if
        end do
      end if
    end do

    ! --- Detect shell species: scan library file if referenced ---
    do l = 1, ntmpl
      lbuf = adjustl(tmpl(l))
      if (lbuf(1:1) == '#' .or. len_trim(lbuf) == 0) cycle
      tok1 = ''
      read(lbuf, *, iostat=ios) tok1
      if (trim(tok1) == 'library') then
        lib_file = ''
        read(lbuf(len_trim(tok1)+2:), *, iostat=ios) lib_file
        if (ios == 0 .and. len_trim(lib_file) > 0) then
          open(newunit=unit_in, file=trim(lib_file), status='old', action='read', iostat=ios)
          if (ios == 0) then
            do
              read(unit_in, '(A)', iostat=ios) lbuf
              if (ios /= 0) exit
              lbuf = adjustl(lbuf)
              if (len_trim(lbuf) == 0 .or. lbuf(1:1) == '#') cycle
              tok1 = ''; tok2 = ''
              read(lbuf, *, iostat=ios) tok1, tok2
              if (ios == 0 .and. trim(tok2) == 'shel') then
                do is = 1, nall
                  if (trim(tok1) == trim(all_sym(is))) then; ishell(is) = .true.; exit; end if
                end do
              end if
            end do
            close(unit_in)
          end if
        end if
        exit
      end if
    end do

    ! --- Write output input.gin ---
    write(numstr, '(I0)') config_idx

    open(newunit=unit_out, file=trim(confdir)//'/input.gin', &
      status='replace', action='write', iostat=ios)
    if (ios /= 0) then; deallocate(tmpl); return; end if

    do l = 1, ntmpl
      if (l == struc_idx) then
        write(unit_out,'(A)') 'cell'
        write(unit_out,'(6(F10.4,2X))') a, b, c, al, be, ga
        write(unit_out,'(A)') 'frac'
        do at = 1, nat
          sym_mapped = trim(sym_arr(at))
          do im = 1, nmap
            if (trim(sym_mapped) == trim(map_from(im))) then
              sym_mapped = trim(map_to(im)); exit
            end if
          end do
          write(unit_out,'(A,2X,A,2X,3(F11.7,2X))') trim(sym_mapped), 'core', &
            coords(at,1), coords(at,2), coords(at,3)
          do is = 1, nall
            if (trim(sym_mapped) == trim(all_sym(is)) .and. ishell(is)) then
              write(unit_out,'(A,2X,A,2X,3(F11.7,2X))') trim(sym_mapped), 'shel', &
                coords(at,1), coords(at,2), coords(at,3)
              exit
            end if
          end do
        end do
      else
        outbuf_gulp = tmpl(l)
        lbuf = adjustl(outbuf_gulp)
        if (lbuf(1:15) == '# sod_type_map ') cycle
        do while (index(outbuf_gulp, '@configuration_number@') /= 0)
          idx_tok = index(outbuf_gulp, '@configuration_number@')
          outbuf_gulp = outbuf_gulp(1:idx_tok-1)//trim(numstr)//outbuf_gulp(idx_tok+len_marker:)
        end do
        write(unit_out,'(A)') trim(outbuf_gulp)
      end if
    end do

    close(unit_out)

    ! --- Copy library files ---
    do l = 1, ntmpl
      lbuf = adjustl(tmpl(l))
      if (lbuf(1:1) == '#') cycle
      tok1 = ''
      read(lbuf, *, iostat=ios) tok1
      if (trim(tok1) == 'library') then
        lib_file = ''
        read(lbuf(len_trim(tok1)+2:), *, iostat=ios) lib_file
        if (ios == 0 .and. len_trim(lib_file) > 0) then
          call execute_command_line('cp '//trim(lib_file)//' '//trim(confdir)//'/', &
            exitstat=ios_cp)
        end if
      end if
    end do

    deallocate(tmpl)
    ok = .true.
  end subroutine sw_write_gulp

  ! ---------------------------------------------------------------------------
  !  VASP POSCAR writer (FILER=11)
  ! ---------------------------------------------------------------------------

  subroutine sw_write_poscar(confdir, nat, sym_arr, coords, &
      a, b, c, al, be, ga, ok)
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: nat
    character(len=*), intent(in) :: sym_arr(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    logical, intent(out) :: ok

    integer :: unit_out, at, ios, nsp_uniq, isp
    character(len=5) :: uniq_sym(sw_max_sym)
    integer :: uniq_cnt(sw_max_sym)
    real(real64) :: lv(3, 3)
    logical :: found
    character(len=256) :: outbuf

    ok = .false.

    ! Collect unique species in order of first appearance
    nsp_uniq = 0; uniq_cnt = 0
    do at = 1, nat
      found = .false.
      do isp = 1, nsp_uniq
        if (trim(sym_arr(at)) == trim(uniq_sym(isp))) then
          uniq_cnt(isp) = uniq_cnt(isp) + 1
          found = .true.; exit
        end if
      end do
      if (.not. found .and. nsp_uniq < sw_max_sym) then
        nsp_uniq = nsp_uniq + 1
        uniq_sym(nsp_uniq) = trim(sym_arr(at))
        uniq_cnt(nsp_uniq) = 1
      end if
    end do

    call sw_lat_vectors(a, b, c, al, be, ga, lv)

    open(newunit=unit_out, file=trim(confdir)//'/POSCAR', &
      status='replace', action='write', iostat=ios)
    if (ios /= 0) return

    write(unit_out,'(A)') 'Generated by SOD'
    write(unit_out,'(F10.6)') 1.0_real64
    write(unit_out,'(3(F12.6,2X))') lv(1,1), lv(2,1), lv(3,1)
    write(unit_out,'(3(F12.6,2X))') lv(1,2), lv(2,2), lv(3,2)
    write(unit_out,'(3(F12.6,2X))') lv(1,3), lv(2,3), lv(3,3)
    outbuf = trim(uniq_sym(1))
    do isp = 2, nsp_uniq
      outbuf = trim(outbuf)//' '//trim(uniq_sym(isp))
    end do
    write(unit_out,'(A)') trim(outbuf)
    do isp = 1, nsp_uniq
      write(unit_out,'(I0,A)',advance='no') uniq_cnt(isp), ' '
    end do
    write(unit_out,*)
    write(unit_out,'(A)') 'Direct'
    do isp = 1, nsp_uniq
      do at = 1, nat
        if (trim(sym_arr(at)) == trim(uniq_sym(isp))) then
          write(unit_out,'(3(F12.8,2X))') coords(at,1), coords(at,2), coords(at,3)
        end if
      end do
    end do

    close(unit_out)
    ok = .true.
  end subroutine sw_write_poscar

  ! ---------------------------------------------------------------------------
  !  CASTEP writer (FILER=12)
  ! ---------------------------------------------------------------------------

  subroutine sw_write_castep(confdir, nat, sym_arr, coords, &
      a, b, c, al, be, ga, config_idx, ok)
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: nat, config_idx
    character(len=*), intent(in) :: sym_arr(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    logical, intent(out) :: ok

    character(len=200), allocatable :: tmpl(:)
    integer, parameter :: len_marker = 22
    integer :: unit_in, unit_out, ios, ntmpl, struc_idx, l, at, idx_tok
    character(len=200) :: lbuf, outbuf
    character(len=16) :: numstr
    real(real64) :: lv(3, 3)

    ok = .false.
    open(newunit=unit_in, file='template_castep.cell', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    allocate(tmpl(sw_max_tmpl))
    ntmpl = 0; struc_idx = 0
    do
      read(unit_in, '(A)', iostat=ios) lbuf
      if (ios /= 0) exit
      ntmpl = ntmpl + 1
      if (ntmpl > sw_max_tmpl) then; close(unit_in); deallocate(tmpl); return; end if
      tmpl(ntmpl) = lbuf
      if (trim(adjustl(lbuf)) == '@configuration_structure@' .and. struc_idx == 0) struc_idx = ntmpl
    end do
    close(unit_in)
    if (struc_idx == 0) then; deallocate(tmpl); return; end if

    call sw_lat_vectors(a, b, c, al, be, ga, lv)
    write(numstr, '(I0)') config_idx

    open(newunit=unit_out, file=trim(confdir)//'/castep.cell', &
      status='replace', action='write', iostat=ios)
    if (ios /= 0) then; deallocate(tmpl); return; end if

    do l = 1, ntmpl
      if (l == struc_idx) then
        write(unit_out,'(A)') '%BLOCK lattice_cart'
        write(unit_out,'(3(F12.6,2X))') lv(1,1), lv(2,1), lv(3,1)
        write(unit_out,'(3(F12.6,2X))') lv(1,2), lv(2,2), lv(3,2)
        write(unit_out,'(3(F12.6,2X))') lv(1,3), lv(2,3), lv(3,3)
        write(unit_out,'(A)') '%ENDBLOCK lattice_cart'
        write(unit_out,'(A)') '%BLOCK positions_frac'
        do at = 1, nat
          write(unit_out,'(A,2X,3(F11.7,2X))') trim(sym_arr(at)), &
            coords(at,1), coords(at,2), coords(at,3)
        end do
        write(unit_out,'(A)') '%ENDBLOCK positions_frac'
      else
        outbuf = tmpl(l)
        do while (index(outbuf, '@configuration_number@') /= 0)
          idx_tok = index(outbuf, '@configuration_number@')
          outbuf = outbuf(1:idx_tok-1)//trim(numstr)//outbuf(idx_tok+len_marker:)
        end do
        write(unit_out,'(A)') trim(outbuf)
      end if
    end do

    close(unit_out)
    deallocate(tmpl)
    ok = .true.
  end subroutine sw_write_castep

  ! ---------------------------------------------------------------------------
  !  Quantum ESPRESSO writer (FILER=13)
  ! ---------------------------------------------------------------------------

  subroutine sw_write_qe(confdir, nat, sym_arr, coords, &
      a, b, c, al, be, ga, config_idx, ok)
    character(len=*), intent(in) :: confdir
    integer, intent(in) :: nat, config_idx
    character(len=*), intent(in) :: sym_arr(:)
    real(real64), intent(in) :: coords(:, :), a, b, c, al, be, ga
    logical, intent(out) :: ok

    character(len=200), allocatable :: tmpl(:)
    integer, parameter :: len_marker = 22
    integer :: unit_in, unit_out, ios, ntmpl, struc_idx, l, at, idx_tok
    character(len=200) :: lbuf, outbuf
    character(len=16) :: numstr
    real(real64) :: lv(3, 3)

    ok = .false.
    open(newunit=unit_in, file='template_pw.in', status='old', action='read', iostat=ios)
    if (ios /= 0) return

    allocate(tmpl(sw_max_tmpl))
    ntmpl = 0; struc_idx = 0
    do
      read(unit_in, '(A)', iostat=ios) lbuf
      if (ios /= 0) exit
      ntmpl = ntmpl + 1
      if (ntmpl > sw_max_tmpl) then; close(unit_in); deallocate(tmpl); return; end if
      tmpl(ntmpl) = lbuf
      if (trim(adjustl(lbuf)) == '@configuration_structure@' .and. struc_idx == 0) struc_idx = ntmpl
    end do
    close(unit_in)
    if (struc_idx == 0) then; deallocate(tmpl); return; end if

    call sw_lat_vectors(a, b, c, al, be, ga, lv)
    write(numstr, '(I0)') config_idx

    open(newunit=unit_out, file=trim(confdir)//'/pw.in', &
      status='replace', action='write', iostat=ios)
    if (ios /= 0) then; deallocate(tmpl); return; end if

    do l = 1, ntmpl
      if (l == struc_idx) then
        write(unit_out,'(A)') 'CELL_PARAMETERS {angstrom}'
        write(unit_out,'(3(F12.6,2X))') lv(1,1), lv(2,1), lv(3,1)
        write(unit_out,'(3(F12.6,2X))') lv(1,2), lv(2,2), lv(3,2)
        write(unit_out,'(3(F12.6,2X))') lv(1,3), lv(2,3), lv(3,3)
        write(unit_out,'(A)') 'ATOMIC_POSITIONS {crystal}'
        do at = 1, nat
          write(unit_out,'(A,2X,3(F11.7,2X))') trim(sym_arr(at)), &
            coords(at,1), coords(at,2), coords(at,3)
        end do
      else
        outbuf = tmpl(l)
        do while (index(outbuf, '@configuration_number@') /= 0)
          idx_tok = index(outbuf, '@configuration_number@')
          outbuf = outbuf(1:idx_tok-1)//trim(numstr)//outbuf(idx_tok+len_marker:)
        end do
        write(unit_out,'(A)') trim(outbuf)
      end if
    end do

    close(unit_out)
    deallocate(tmpl)
    ok = .true.
  end subroutine sw_write_qe

end module structwriters

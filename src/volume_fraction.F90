module volume_fraction

  use global
  use constants
  use mesh,             only : get_mesh_indices
  use random_lcg,       only : prn, ! set_particle_seed
  use geometry,         only : find_cell, check_cell_overlap
  use geometry_header,  only : Cell, BASE_UNIVERSE 
  use particle_header,  only : deallocate_coord, Particle
  use material_header,  only : Material
  use output,           only : write_message
  use error,            only : fatal_error
  use output_interface

  implicit none

contains

!=========================================================================================
! Volume fractions module does all of the work needed to create a vector of volume 
! fractions the same size as the user defined mesh for said fractions. It creates random
! locations in geometry and finds the material identifiers, whether or not they are 
! fissionable, and then does the math needed to find the fractions based on the number
! of particles run. The more particles the more accurate the volume fractions.
!=========================================================================================

  subroutine run_volfrac()

    real(8)         :: xyz(3)  ! random location in geometry
    integer         :: ijk(3)  ! indces in ufs mesh
    logical         :: in_mesh ! whether point specified is in mesh
    integer         :: volfrac_dim(4) ! dimensions of the volume fraction array
    type(Particle)  :: p
    integer         :: i ! counter variable
    logical         :: found_cell
    type(Cell), pointer :: c => null()
    type(Material), pointer :: m => null()
    
    !Allocate and initialize particle
     call p % intialize()
     p % coord % uvw = (/ 0.5, 0.5, 0.5 /)
     p % coord % universe = BASE_UNIVERSE

    ! Set particle seed
    !call set_particle_seed()

    !Loop over particles
    do i = 1, ufs_vol_res

      !Random location
      xyz(1) = prn()*(ufs_mesh % width(1)*ufs_mesh % dimension(1)) + ufs_mesh % lower_left(1)
      xyz(2) = prn()*(ufs_mesh % width(2)*ufs_mesh % dimension(2)) + ufs_mesh % lower_left(2)
      xyz(3) = prn()*(ufs_mesh % width(3)*ufs_mesh % dimension(3)) + ufs_mesh % lower_left(3)
      p % coord % xyz = xyz
      call deallocate_coord(p % coord0 % next)
      p % coord => p % coord0

      !Identify cell and material
      call find_cell(p, found_cell)
      if (check_overlaps) call check_cell_overlap(p)
      if (.not. found_cell) then
        print*, 'Overlap found, adding particle' 
        ufs_vol_res = ufs_vol_res + 1
        return
      else !Cell is found, then figure out material id of celli
        c => cells(p % coord % cell)
        m => materials( c % material)
      end if

      !If not fissionable, return to next particle, else sort into user defined UFS mesh
      if (.not. m % mat_fissionable) then
        return
      else
        ! Find location in UFS mesh
        call get_mesh_indices(ufs_mesh, xyz, ijk, in_mesh)
        if (.not. in_mesh) then
          message = "Location is not in UFS mesh"
          call fatal_error()
        end if
        ! Bank it into the volume fraction mesh
        volume_frac(1, ijk(1), ijk(2), ijk(3)) = volume_frac(1, ijk(1), ijk(2), ijk(3)) + 1
      end if
    
    end do
    
    !Find volume fractions (math)
    volume_frac=volume_frac/(sum(volume_frac))

    !Write out to xml file
    call write_volfrac_xml()
 
    !print out to user that it is complete
    message = 'Volume fraction file has been written'
    call write_message(1)

  end subroutine run_volfrac()

!========================================================================================
! write_volfrac_xml() creates the xml file containging the volume fractions
! that can be used as an input when the ufs method chosen is approximation
!========================================================================================
    
  subroutine write_volfrac_xml()

    character(len=11) :: node_tag_open
    character(len=12) :: node_tag_close
    integer           :: z ! counter for z loop
    integer           :: y ! counter for y loop
    integer           :: i ! counter for inner loops
    integer           :: xdim ! shorthand for dimensions
    integer           :: ydim 
    integer           :: zdim

    ! Define some variables for readability
    xdim = ufs_mesh % dimension(1)  ! number of data in x direction
    ydim = ufs_mesh % dimension(2)  ! " y "
    zdim = ufs_mesh % dimension(3)  ! "z "

    ! Open the file for writing
    open(UNIT=UNIT_VOLFRAC, FILE='volume_fractions.xml', ACTION='readwrite', &
         STATUS='replace', ACCESS='stream')

    ! Write header to file
    WRITE(UNIT_VOLFRAC, 100)
    100 FORMAT(' ','<volume_fracions>')
    
    ! Write node name
    node_tag_open='<fractions>'
    WRITE(UNIT_VOLFRAC, 110) node_tag_open
    110 FORMAT(' ', A13)
    
    ! Write data
    do z = 1, zdim
      
      do y = 1, ydim
        ! Write out entire horizontal line at once
        WRITE(UNIT_VOLFRAC, '(1X, f12.5, (8.5))') &
             volume_frac(1,1,y,z), (volume_frac(1,i,y,z), i=2, xdim)
      end do
      ! Add a line between Z to show XvsY slice
      WRITE(UNIT_VOLFRAC, '(A)') new_line(' ') ! 2005 fortran?
    end do 
    
    ! Write closing node tag
    node_tag_close='</fractions>'
    WRITE(UNIT_VOLFRAC, 190) node_tag_close
    190 FORMAT(' ', A14)

    ! Write closing tag
    WRITE(UNIT_VOLFRAC, 200) 
    200 FORMAT(' ', '</volume_fractions>')
    
    ! Close the file
    close(UNIT=UNIT_VOLFRAC)

  end subroutine write_volfrac_xml()

end module volume_fraction

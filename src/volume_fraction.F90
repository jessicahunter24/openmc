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

    ! Open the file for writing
    open(UNIT=UNIT_VOLFRAC, FILE='volume_fractions.xml', ACTION='readwrite',STATUS='replace', ACCESS='stream')
    ! Write header to file
    ! Write node name
    ! Write data
    ! Write closing tag
    
    ! Close the file
    close(UNIT=UNIT_VOLFRAC)

    !create the file
    ! call file_create( , "volume_fractions.xml")
    !volfrac_dim(1) = 1
    !volfrac_dim(2:4) = ufs_mesh % dimension
    !call write_double_3Darray( , datavarname, "data name", length(4))
    !call file_open( , "Volume_fractions.xml", "w)" !also r for reeeeeead ;)  

  end subroutine write_volfrac_xml()

end module volume_fraction

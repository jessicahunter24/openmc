module volume_fraction

  use global
  use mesh,             only : get_mesh_indices
  use random_lcg,       only : prn, ! set_particle_seed
  use geometry,         only : find_cell, check_cell_overlap
  use geometry_header,  only : Cell, BASE_UNIVERSE 
  use particle_header,  only : deallocate_coord, Particle
  use material_header,  only : Material
  use output,           only : write_message
  use error,            only : fatal_error

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

    !write out to volume fraction xml file 

    !Print to user that it has been completed.
     ! see header above

  end subroutine run_volfrac()

end module volume_fraction

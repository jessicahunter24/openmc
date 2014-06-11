module volume_fraction

  use global

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

    integer   :: particlenum !Number of particles to run simulation

    !Pull in the data from the user on the mesh size and how many particles to run
    ufs_mesh % lower_left(3)
    ufs_mesh % upper_right(3)
    ufs_mesh % width(3)
    ufs_mesh % dimension(3)
    ufs_vol_res
    
    !Allocate arrays

    !Loop over particles

      !Random location

      !Identify cell

      !Identify material

      !Check if fissionable

      !If not fissionable, return to next particle, else sort into user defined UFS mesh
    
    !Find volume fractions (math)
    
    !write out to volume fraction xml file 

    !Print to user that it has been completed.


  end subroutine run_volfrac()

end module volume_fraction

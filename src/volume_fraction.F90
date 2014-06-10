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


  end subroutine run_volfrac()

end module volume_fraction

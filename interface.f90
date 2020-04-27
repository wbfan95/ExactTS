! An interface connecting ExactTS and your PES
! 
! Initializing your potential energy surface
subroutine TSinit
call prepot
end subroutine

! Hock your PES here
! Unit: Hartree and Angstrom, please! 
subroutine TSenergy(Natoms, coord, energy)
implicit none
integer, intent(in) :: Natoms
real(8), intent(in) :: coord(3, Natoms)
!real(8) coordT(Natoms, 3)
real(8), intent(out) :: energy
real(8) CART, PENGYGS
INTEGER NATOM
PARAMETER (NATOM=25)
COMMON/USROCM/ PENGYGS
COMMON/USRICM/ CART(NATOM,3)

cart = transpose(coord) / 0.5291772d0
call pot
energy = PENGYGS

end subroutine

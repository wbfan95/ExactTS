! ********************************************************
! Program: ExactTS
! Author: Wenbin, FAN (fanwenbin@shu.edu.cn, langzihuigu@qq.com)
!
! Release: Dec. 7, 2019
! Update: Jul. 19, 2020
!
! Hartree and Angstrom are used here! 
! ********************************************************

program main
!use lapack95
!use fakeG
implicit none
integer Natoms
real(8), allocatable :: coord(:,:), coord1D(:) ! coordinates 3D and 1D
real(8), allocatable :: grad(:,:), grad1D(:)
real(8), allocatable :: hessian(:,:), hessianT(:,:) ! Hessian
integer, allocatable :: ipiv(:) ! for matrix inversion
real(8), allocatable :: coord1DNew(:)
character(4), allocatable :: eleName(:) ! element names
real(8), allocatable :: eleMass(:)
real(8) rand
real(8) energy
!real(8) error ! difference between current and the last Cartesian geometries. 
integer step, convIF
integer i, j, k, pinvInfo
real(8) rr, mc(2), R, x1, y1, x2, y2, gam, rmgh1, rmgh2 ! output for bonds


! PES initialization
call TSinit

! read coordinates from TS.xyz ! angstrom
open(999, file='TS.xyz', status='old', action='read') ! xyz -> 999
read(999, *) Natoms
allocate(coord(3, Natoms), eleName(Natoms))
read(999, *)
do i = 1, Natoms
read(999, *) eleName(i), coord(:, i)
!# write(*,'(3F11.6)') coord(:,i)
end do
close(999)

! output file ! out -> 688, ene -> 363
open(688, file='opt.xyz', status='unknown', action='write')
open(363, file='energy.txt', status='unknown', action='write')
open(373, file='freq_out.txt', action='write') ! fre -> 373
open(566, file='bond_out.txt', action='write') ! bon -> 566

! allocate
allocate(eleMass(Natoms))
allocate(coord1D(3 * Natoms))
coord1D = reshape(coord, (/Natoms * 3/)) ! reshape coordinates to 1D array
allocate(grad(3, Natoms), grad1D(3 * Natoms))
allocate(hessian(Natoms * 3, Natoms * 3), &
         hessianT(Natoms * 3, Natoms * 3), &
         ipiv(Natoms * 3))
allocate(coord1DNew(3 * Natoms))

! Get mass
call eleName2Mass(Natoms, eleName, eleMass)

!error = 1.d0
step = 1
convIF = 0
do while (convIF .lt. 3.5d0)
!do while (step .lt. 3)
write(*,'(A12,I5)') 'Steps: ', step
coord = reshape(coord1D, (/3, Natoms/))
call TSenergy(Natoms, coord, energy)
write(363,*) energy * 627.503
write(*,'(A12,F11.6,A)') 'Energy: ', energy * 627.503, ' kcal/mol'

! 1st derivation
call gradCalc(Natoms, coord, grad)
grad1D = reshape(grad, (/Natoms * 3/))

! 2nd derivation / Hessian matrix
call hessianCalc(Natoms, coord1D, hessian)
call vibAna(Natoms, eleMass, hessian, 373)
! Inversion of Hessian matrix
!hessianT = hessian
!call getrf(hessianT, ipiv)
!call getri(hessianT, ipiv)
!#if (maxval(hessianT - hessian) .lt. 1E-9) then
!#write(*,*) 'Inversion of Hessian matrix failed! '
!#!!stop
!#end if
call pinv(Natoms * 3, hessian, hessianT, pinvInfo)
if (pinvInfo .ne. 0) stop
!# write(*,*) matmul(hessianT, hessian) ! Check the identity

! Evolution: x' = x - dx * Hessian^T
coord1DNew = coord1D - matmul(hessianT, grad1D) *0.01d0!* (rand() * 0.2d0 + 0.1d0)
call convCrit(Natoms, coord1D, coord1DNew, grad1D, convIF)
!# write(*,*) shape(coord1D)
!# write(*,*) shape(hessianT)
! Export results
call cleanGeom(Natoms, eleMass, coord1DNew)
call writeXYZ(Natoms, eleName, reshape(coord1DNew, (/3, Natoms/)), 688)

! write R, r, gamma
rr = dsqrt((coord1D(4) - coord1D(1)) ** 2d0 + (coord1D(5) - coord1D(2)) ** 2d0)
mc(1) = (coord1D(4) * 24.305 + coord1D(1) * 1.008) / (24.305 + 1.008)
mc(2) = (coord1D(5) * 24.305 + coord1D(2) * 1.008) / (24.305 + 1.008)
R = dsqrt((coord1D(7) - mc(1)) ** 2d0 + (coord1D(8) - mc(2)) ** 2d0)
x1 = coord1D(1) - coord1D(4)
y1 = coord1D(2) - coord1D(5)
x2 = coord1D(7) - mc(1)
y2 = coord1D(8) - mc(2)
gam = dacos((x1 * x2 + y1 * y2)/dsqrt((x1*x1 + y1*y1) * (x2*x2 + y2*y2))) * 180d0 / dacos(-1d0)
rmgh1 = dsqrt((coord1D(4)-coord1D(1))**2d0 + (coord1D(5)-coord1D(2))**2d0 + (coord1D(6)-coord1D(3))**2d0)
rmgh2 = dsqrt((coord1D(4)-coord1D(7))**2d0 + (coord1D(5)-coord1D(8))**2d0 + (coord1D(6)-coord1D(9))**2d0)
write(566, '(5F15.6)') R/0.529177210903d0, rr/0.529177210903d0, gam, rmgh1/0.529177210903d0, rmgh2/0.529177210903d0

! loop
coord1D = coord1DNew
step = step + 1
write(*,*) ''
!error = 1d-7
end do ! while loop

call fakeG(Natoms, eleName, coord)

end program

subroutine convCrit(Natoms, coord1D, coord1DNew, grad1D, convIF)
implicit none
integer Natoms, i, j, k
real(8) coord1D(Natoms * 3), coord1DNew(Natoms * 3)
real(8) grad1D(Natoms * 3)
real(8) maxF, rmsF, maxR, rmsR
integer :: convIF

convIF = 0
maxF = 0d0
rmsF = 0d0
maxR = 0d0
rmsR = 0d0

! max force
maxF = maxval(abs(grad1D))
if (maxF .lt. 2d-6) convIF = convIF + 1

! RMS force
do i = 1, Natoms * 3
rmsF = rmsF + grad1D(i) ** 2d0
end do
rmsF = dsqrt(rmsF) / (Natoms * 3d0)
if (rmsF .lt. 1d-6) convIF = convIF + 1

! max displacement
maxR = maxval(abs(coord1D - coord1DNew))
if (maxR .lt. 6d-6) convIF = convIF + 1

! RMS displacement
do i = 1, Natoms * 3
rmsR = rmsR + (coord1D(i) - coord1DNew(i)) ** 2d0
end do
rmsR = dsqrt(rmsR) / (Natoms * 3d0)
if (rmsR .lt. 4d-6) convIF = convIF + 1

! ouput
write(*,'(2(A12,E11.2,A4))') 'Max force: ', maxF, '', 'RMS force: ', rmsF
write(*,'(2(A12,E11.2,A4))') 'Max disp.: ', maxR, '', 'RMS disp.: ', rmsR

end subroutine

! Vibrational analysis
subroutine vibAna(Natoms, eleMass, hessian, fileID)
use lapack95
implicit none
integer Natoms
real(8) hessian(3 * Natoms, 3 * Natoms), hessianMW(3 * Natoms, 3 * Natoms) ! mass-weighted
real(8) eleMass(Natoms), eigenvalue(3 * Natoms), freq(3 * Natoms)!, wavenumber(3 * Natoms)
real(8) tmp ! temp value for sorting freq
integer i, j, m, n
integer fileID

call convMWhessian(Natoms, eleMass, hessian, hessianMW)

! Eigenvalue
call syev(hessianMW, eigenvalue)
!# write(*,'(6F11.6)') eigenvalue

call calcFreq(Natoms, eigenvalue, freq)

write(fileID, '(300F11.2)') freq

! sort modes by absolute value
do i = 1, Natoms * 3
do j = i + 1, Natoms * 3
if (abs(freq(i))>abs(freq(j))) then
tmp=freq(i)
freq(i)=freq(j)
freq(j)=tmp
end if
end do
end do

! output freq
write(*,'(A)') 'Freq(cm-1): '
write(*,'(6F12.2)') freq(7:)

end subroutine

subroutine convMWhessian(Natoms, eleMass, hessian, hessianMW)
implicit none
integer :: i,j,k,l,m,n,Natoms
real(8) hessian(3 * Natoms, 3 * Natoms), hessianMW(3 * Natoms, 3 * Natoms)
real(8) eleMass(Natoms)

hessianMW = 0d0
! Convert to mass-weighted Hessian
do j = 0, Natoms - 1
do i = 0, Natoms - 1
! 9 elements for each atom ! 3 * 3 matrix, ie. xx, xy, xz, yx, yy, yz, zx, zy, zz
do n = 1, 3
do m = 1, 3
hessianMW(i*3+m, j*3+n) = hessian(i*3+m, j*3+n) / dsqrt(eleMass(i+1) * eleMass(j+1))
end do
end do
end do
end do

end subroutine

subroutine calcFreq(Natoms, eigenvalue, freq)
implicit none
integer :: Natoms, i, j, k
real(8) eigenvalue(3 * Natoms), freq(3 * Natoms)
real(8) amu2kg, a2m, au2J
real(8), parameter :: pi = 4.d0*atan(1.d0)

! To frequency ! Below 13 lines were referred Sobereva's `Hess2freq`. 
amu2kg=1.66053878D-27
a2m=1D-10
au2J=4.35974434D-18
eigenvalue=eigenvalue*au2J/a2m**2/amu2kg !convert force constant from a.u. to SI
do i = 1, Natoms * 3
if (eigenvalue(i)<0) then
freq(i) = -dsqrt(abs(eigenvalue(i))) / (2 * pi)
else
freq(i) = dsqrt(eigenvalue(i)) / (2 * pi)
end if
end do

freq = freq/2.99792458D10

end subroutine

! Shift centroid to origin
subroutine cleanGeom(Natoms, eleMass, coord)
integer Natoms
real(8) coord(3 * Natoms), coord3(3, Natoms), cm(3) ! center of mass
real(8) eleMass(Natoms), sumMass(3), totalMass

coord3 = reshape(coord, (/3, Natoms/))
!# coord3(1,:) = coord3(1,:) + 1d0
!# write(*,*) 'Coordintaes'
!# write(*,'(3F11.6)') coord3

! Find center of mass
sumMass = 0d0
do j = 1, Natoms
do i = 1, 3
sumMass(i) = sumMass(i) + coord3(i, j) * eleMass(j)
!# write(*,'(3F11.6)') coord3(i, j) * eleMass(j)
end do
end do
!# write(*,*) 'Sum of mass * r'
!# write(*,'(3F11.6)') sumMass


! total mass
totalMass = sum(eleMass)
!# write(*,*) 'Total mass: '
!# write(*,'(3F11.6)') totalMass

! center of mass
cm = sumMass / totalMass
!# write(*,*) 'Center of mass: '
!# write(*,'(3F11.6)') cm

! shift CM to origin
do j = 1, Natoms
coord3(:, j) = coord3(:, j) - cm
end do
!# write(*,*) 'Coordintaes shifted '
!# write(*,'(3F11.6)') coord3

coord = reshape(coord3, (/Natoms/))

end subroutine

! Get each mass from element name
subroutine eleName2Mass(Natoms, eleName, eleMass)
implicit none
integer Natoms
character(4) :: allName(22) = &
(/'H', 'He',   'Li',   'Be',   'B', &
 'C',  'N',    'O',    'F',    'Ne', &
 'Na', 'Mg',   'Al',   'Si',   'P', &
 'S',  'Cl',   'Ar',   'K',    'Ca', &
 'D', 'C13'/)
real(8) :: allMass(22) = &
(/1.008d0,      4.002602d0, 6.94d0,       9.0121831d0,    10.81d0, &
 12.000d0,      14.007d0,   15.999d0,     18.998403163d0, 20.1797d0, &
 22.98976928d0, 24.305d0,   26.9815385d0, 28.085d0,       30.973761998d0, &
 32.06d0,       35.45d0,    39.948d0,     39.0983d0,      40.078d0, &
 2.0141d0,      13.00335483521d0/)
character(4) eleName(Natoms)
real(8), intent(out) :: eleMass(Natoms)
integer i, j

do i = 1, Natoms
do j = 1, 22
!# write(*,*) i, j, eleName(i), allName(j)
if (trim(eleName(i)) == trim(allName(j))) then
eleMass(i) = allMass(j)
exit
!else
!write(*, *) '[ERROR] Element not supported: ', eleName(i)
end if
end do
end do
!# write(*,*) 'Element mass: '
!# write(*,*) eleMass
end subroutine

! Export coordinates to .xyz file
subroutine writeXYZ(Natoms, eleName, coord, fileID)
implicit none
integer Natoms, fileID, i
character(4) eleName(Natoms)
character(64) info
real(8) coord(3, Natoms)

write(fileID,'(I4)') Natoms
write(fileID,'(A)') 'TEST'
do i = 1, Natoms
write(fileID,'(A,3F18.10)') eleName(i), coord(:,i)
end do
end subroutine

! Calculate the Hessian matrix
subroutine hessianCalc(Natoms, coord1D, hessian)
implicit none
integer Natoms
real(8) coord1D(3 * Natoms)
real(8) hessian(Natoms * 3, Natoms * 3)
real(8) energy2nd

integer i, j

do j = 1, Natoms * 3
do i = 1, Natoms * 3
call der2nd(Natoms, coord1D, i, j, energy2nd)
hessian(i, j) = energy2nd
!# write(*,*) hessian(i,j)
end do
end do

!# ! check Hessian: $H_{i,j} = H_{j,i}$ ! check passed! 
!# do j = 1, Natoms * 3
!# do i = 1, Natoms * 3
!# write(*,'(3F11.6)') hessian(i, j), hessian(j, i), hessian(i, j) - hessian(j, i)
!# end do
!# end do

end subroutine

! Calculate the 2nd derivation
subroutine der2nd(Natoms, coord1D, i, j, energy)
implicit none
integer Natoms, i, j
real(8) coord1D(3 * Natoms), energy
real(8) :: delta = 1d-5, delta4squared = 4d-10
real(8) e1, e2, e3, e4 ! See ref: William, Fortran Reciept, Page 182
real(8) bkup1, bkup2 ! backup energy

! e1: x + delta, y + delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) + delta
coord1D(j) = coord1D(j) + delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e1)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e2: x + delta, y - delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) + delta
coord1D(j) = coord1D(j) - delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e2)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e3: x - delta, y + delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) - delta
coord1D(j) = coord1D(j) + delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e3)
coord1D(i) = bkup1
coord1D(j) = bkup2

! e4: x - delta, y - delta
bkup1 = coord1D(i)
bkup2 = coord1D(j)
coord1D(i) = coord1D(i) - delta
coord1D(j) = coord1D(j) - delta
call TSenergy(Natoms, reshape(coord1D, (/3, Natoms/)), e4)
coord1D(i) = bkup1
coord1D(j) = bkup2

energy = (e1 - e2 - e3 + e4) / delta4squared

end subroutine

! Calculate the gradient
subroutine gradCalc(Natoms, coord, grad)
implicit none
integer Natoms
real(8) coord(3, Natoms)
real(8) grad(3, Natoms)
real(8) eneFor, eneBak, coordBackup ! forward, backward, backup
real(8) :: delta = 1d-5 ! displacement
integer i, j

do j = 1, Natoms
do i = 1, 3 ! x,y,z
! finite displacement
coordBackup = coord(i,j)
coord(i,j) = coord(i,j) + delta
call TSenergy(Natoms, coord, eneFor)
coord(i,j) = coord(i,j) - 2d0 * delta
call TSenergy(Natoms, coord, eneBak)
coord(i,j) = coordBackup
! calculate gradient
grad(i,j) = (eneFor - eneBak) / (2d0 * delta)
end do
!# write(*,'(3F11.6)') grad(:,j)
end do
end subroutine

subroutine pinv(N, A, pinvA, info)
use lapack95, only: gesvd
implicit none
integer N
real(8), intent(in) :: A(N,N)
real(8), intent(out) :: pinvA(N,N)
real(8) U(N,N), S(N), VT(N,N)!, WW(3)
real(8) S1(N,N), V(N,N), VS1(N,N)
real(8) work(N*N), lwork
integer, intent(out) :: info
integer i

call dgesvd('A', 'A', N, N, A, N, S, U, N, VT, N, work, N*N, info)

!#write(*,*) 'A'
!#write(*,'(6F11.4)') A

! pseudo inversion for S
S1 = 0d0
do i = 1, N
if (S(i) .gt. 1d-10) S1(i,i) = 1/S(i)
end do

! A-1 = V * S-1 * U^T
VS1 = matmul(transpose(VT), S1)
pinvA = matmul(VS1, transpose(U))

end subroutine

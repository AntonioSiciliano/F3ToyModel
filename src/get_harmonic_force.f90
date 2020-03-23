

! Subroutine to calculate the force on each atom from the
! force constants. The force constants are given as input.        

subroutine get_harmonic_force_from_fc(fc,u,forces, natom)

  implicit none

  double precision, dimension(3,3,natom,natom), intent(in) :: fc
  double precision, dimension(natom,3), intent(in) :: u
  double precision, dimension(natom,3), intent(out) :: forces 

  integer :: natom
  integer :: i, j, k, l 

  !natom = size(u(:,1))

  do i = 1, natom
    do j = 1, 3
      forces(i,j) = 0.0d0
      do k = 1, natom
        do l = 1, 3     
          forces(i,j) = forces(i,j) - fc(j,l,i,k) * u(k,l)
        end do
      end do
    end do
  end do
        
end subroutine get_harmonic_force_from_fc

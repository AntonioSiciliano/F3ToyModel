
! This subroutine gives the contribution to the free energy of a 
! harmonic potential. What it is calculated is 1/2 u * phi * u 

subroutine v2_0_harmonic (u,phi,f, nat)

  implicit none

  integer, intent(in) :: nat
  double precision, dimension(nat,3), intent(in) :: u
  double precision, dimension(3,3,nat,nat), intent(in) :: phi
  double precision, intent(out) :: f

  integer :: i, j, alpha, beta
  double precision, dimension(:), allocatable :: u_vec, u_prod
  double precision, dimension(:,:), allocatable :: phitwo
  
  !nat = size(u(:,1))

  allocate(phitwo(3*nat,3*nat))
  allocate(u_vec(3*nat))
  allocate(u_prod(3*nat))

  do i = 1, nat
    do alpha = 1, 3
      u_vec(3*(i-1)+alpha) = u(i,alpha)
    end do
  end do

  do i = 1, nat
    do j = 1, nat
      do alpha = 1, 3
        do beta = 1, 3
          phitwo(3*(i-1)+alpha,3*(j-1)+beta) = phi(alpha,beta,i,j)
        end do   
      end do
    end do
  end do

  call dgemv('N',3*nat,3*nat,1.0d0,phitwo,3*nat,u_vec,1,0.0d0,u_prod,1)

  f = 0.5d0 * dot_product(u_vec,u_prod)

  deallocate(phitwo, u_vec, u_prod)

end subroutine v2_0_harmonic


! This subroutine assigns nearest neighbors if an atom in a 
! supercell in a fcc lattice

subroutine assign_NN (tau_sc, at_sc, ityp_sc, nn_at, nn_vect, natsc)

  implicit none

  integer, dimension(natsc), intent(in) :: ityp_sc
  double precision, dimension(3,natsc), intent(in) :: tau_sc
  double precision, dimension(3,3), intent(in) :: at_sc
  integer, dimension(natsc,6), intent(out) :: nn_at  
  double precision, dimension(natsc,6,3), intent(out) :: nn_vect
  integer :: natsc

  integer,parameter :: far =2
  double precision, dimension((2*far+1)**3,3) :: t_vectors
  double precision, dimension(3) :: vect
  integer :: i, k, j, l
  double precision :: nn_dist, dist

  !natsc = size(tau_sc(1,:))

  ! Distance between nearest neighbors

  nn_dist = 0.5d0

  ! Create supercell lattice vectors

  l = 1

  do i = -far,far  
    do j = -far,far    
      do k = -far,far   
        t_vectors(l,:) = dble(i) * at_sc(:,1) + dble(j) * at_sc(:,2) + dble(k) * at_sc(:,3)
        l = l + 1    
      end do
    end do
  end do

  nn_at = 0

  do i = 1, natsc
    k = 0  
    do j = 1, natsc
      do l = 1, (2*far+1)**3
        vect(:) = tau_sc(:,i)-tau_sc(:,j)-t_vectors(l,:)
        dist = dsqrt(dot_product(vect(:),vect(:))) 
        if ( abs(dist - nn_dist) .lt. 1.0d-5 ) then
          k = k + 1
          nn_at(i,k) = j
          nn_vect(i,k,:) = vect(:)
        end if   
      end do 
    end do
    if (k .ne. 6) then
      print *, ' Not all or too many NN assigned', k
      stop
    end if
  end do

end subroutine assign_NN

subroutine assign_PM(nat_sc, index_vect, index_mat, slv_idx)
  implicit none

  integer,intent(IN)   :: nat_sc
  double precision, intent(IN) :: index_vect(nat_sc,6,3)
  integer,intent(IN)   :: index_mat(nat_sc,6)
  integer,intent(OUT) :: slv_idx(nat_sc,6)

  double precision :: dist(3)
  integer :: i,j, jj, is
  double precision ::  special(3,6), proj, aux

  special(:,1) = (/ 1.d0, 0.d0, 0.d0 /)
  special(:,2) = (/-1.d0, 0.d0, 0.d0 /)
  special(:,3) = (/ 0.d0, 1.d0, 0.d0 /)
  special(:,4) = (/ 0.d0,-1.d0, 0.d0 /)
  special(:,5) = (/ 0.d0, 0.d0, 1.d0 /)
  special(:,6) = (/ 0.d0, 0.d0,-1.d0 /)


  do i = 1,nat_sc
  write(*,*) "neighbours of atom", i
  do is = 1,6 ! six special directions
    proj = 0.d0
    do jj = 1,6 ! first neighbours
      j = index_mat(i,jj)
!       dist(:) = (coords(j,:)+atom_transl(i,j,:)) - coords(i,:)
!       dist(:) = (coords(j,:)+(tau_sc(:,j)-tau_sc(:,i)) - coords(i,:)
      dist(:) = index_vect(i,jj,:)
      aux = SUM( dist*special(:,is) )
      IF( aux > proj ) THEN
        proj = aux
        slv_idx(i,is) = j
      ENDIF

    enddo
   write(*,'(2(a,i3),2(a,3f10.4))') "type", is, " idx", slv_idx(i, is) 
  enddo
  enddo

end subroutine assign_PM




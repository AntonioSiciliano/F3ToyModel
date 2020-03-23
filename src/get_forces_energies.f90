SUBROUTINE get_forces_energies(phi_sc_harmonic,u_vector,type_cal,ityp_sc,at_sc, tau_sc, &
                               b,c,p2,p3,p4,p5,p6,p4x,p3x,p4f,p4g, &
                               forces, v, n_random, nat_sc)
  
  implicit none

  integer :: n_random, nat_sc
  DOUBLE PRECISION, dimension(3,3,nat_sc,nat_sc), intent(in) :: phi_sc_harmonic 
  DOUBLE PRECISION, dimension(n_random,nat_sc,3), intent(in) :: u_vector
  CHARACTER (len=6), intent(in) :: type_cal  
  INTEGER,  dimension(nat_sc), intent(in) :: ityp_sc
  DOUBLE PRECISION, dimension(3,3), intent(in) :: at_sc
  DOUBLE PRECISION, dimension(3, nat_sc), intent(in) :: tau_sc
  DOUBLE PRECISION, intent(in) :: b,c,p2, p3, p3x, p4, p4x, p4f, p4g, p6, p5
  DOUBLE PRECISION, dimension(n_random,nat_sc,3), intent(out) :: forces  
  DOUBLE PRECISION, dimension(n_random), intent(out) :: v        

  
  integer, dimension(:,:), allocatable :: nn_at, slv_idx
  DOUBLE PRECISION, dimension(:,:,:), allocatable :: nn_vect
  DOUBLE PRECISION :: v_anharmonic
  DOUBLE PRECISION, dimension(3) :: dispvect

  
  ! Get integer numbers
  integer i1, i2, i3

  ! Get the nearest neighbors 
  print *, "BEFORE ALLOCATION"
  call flush()

  ALLOCATE(nn_at(nat_sc,6))
  ALLOCATE(nn_vect(nat_sc,6,3))

  print *, "ASSIGN NEAR NEIGHBOUR:"
  print *, "TAU_SC:"
  do i1 =1 , nat_sc 
    print *, tau_sc(:, i1)
  enddo
  print *, "AT_SC:"
  do i1 =1 , 3 
    print *, at_sc(:, i1)
  enddo

  print *, "ITYP_SC: ", ityp_sc(:)
  call flush()
  ! NN : TODO, pass to NN a parameter that is the nearest neighbour distance.
  !      This can be done in the initialization function and passed through it.
  call assign_NN (tau_sc,at_sc,ityp_sc,nn_at, nn_vect, nat_sc)
  
  ALLOCATE(slv_idx(nat_sc,6))

  print *, "ASSIGN PM"
  call flush()
  call assign_PM (nat_sc, nn_vect, nn_at, slv_idx)

  print *, "AFTER ALLOCATION"
  call flush()

  ! Get forces for each configuration

  forces=0.d0

  CONFIGURATIONS_LOOP_2 : &
  DO i1 = 1, n_random

    print *, "N CONF:", i1 
    call flush()
  
   !-------------------------------------
   ! CAREFUL:
   ! Not including effective charges here
   call get_harmonic_force_from_fc(phi_sc_harmonic,u_vector(i1,:,:),forces(i1,:,:), nat_sc)

    ! Add anharmonic bit to some of the atoms on top
    ANHARMONIC_FORCE_FIELDS : &
    IF(type_cal == 'pdhxx' ) THEN
      do i2 = 1, nat_sc
        !!!!!!!!!
        ! CAREFUL
        ! 
        ! We assume that the anharmonic atom is the atom 2
        if (ityp_sc(i2) .eq. 2) then
          dispvect(:) = u_vector(i1,i2,:)
          do i3 = 1, 6
            dispvect(:) = dispvect(:) - u_vector(i1,nn_at(i2,i3),:) / 6.d0
          end do
          forces(i1,i2,1) = forces(i1,i2,1) - 4.d0*b*dispvect(1)**3.d0 &
                            - 2.d0*c*dispvect(1)*(dispvect(2)**2.d0+dispvect(3)**2.d0)
          forces(i1,i2,2) = forces(i1,i2,2) - 4.d0*b*dispvect(2)**3.d0 &
                            - 2.d0*c*dispvect(2)*(dispvect(1)**2.d0+dispvect(3)**2.d0)
          forces(i1,i2,3) = forces(i1,i2,3) - 4.d0*b*dispvect(3)**3.d0 &
                            - 2.d0*c*dispvect(3)*(dispvect(1)**2.d0+dispvect(2)**2.d0)
  !      else if (ityp_sc(i2) .eq. 1) then
          do i3 = 1, 6
  !          dispvect(:) = u_vector(i1,nn_at(i2,i3),:) 
  !          do i4 = 1, 6
  !            dispvect(:) = dispvect(:) - u_vector(i1,nn_at(nn_at(i2,i3),i4),:) / 6.d0
  !          end do 
            forces(i1,nn_at(i2,i3),1) = forces(i1,nn_at(i2,i3),1) + (4.d0*b*dispvect(1)**3.d0 &
                              + 2.d0*c*dispvect(1)*(dispvect(2)**2.d0+dispvect(3)**2.d0)) / 6.d0
            forces(i1,nn_at(i2,i3),2) = forces(i1,nn_at(i2,i3),2) + (4.d0*b*dispvect(2)**3.d0 &
                              + 2.d0*c*dispvect(2)*(dispvect(1)**2.d0+dispvect(3)**2.d0)) / 6.d0
            forces(i1,nn_at(i2,i3),3) = forces(i1,nn_at(i2,i3),3) + (4.d0*b*dispvect(3)**3.d0 &
                              + 2.d0*c*dispvect(3)*(dispvect(1)**2.d0+dispvect(2)**2.d0)) / 6.d0
          end do
        end if
         !print*, forces(1,i2,1),forces(1,i2,2),forces(1,i2,3)
      end do
        
    ELSE IF (type_cal == 'pbtex') THEN ! ANHARMONIC_FORCE_FIELDS 
      !
      do i2 = 1, nat_sc

! forces,  x component !
  forces(i1,i2,1)= forces(i1,i2,1) &
        +p2*(-2.d0/sqrt(2.d0)**2) * ( &                ! Anharmonic contribution
        ( (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1)) & 
         +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1)))&
          ) &
        +p3*(-3.d0/sqrt(2.d0)**3)*&    ! 3rd order contribution
         ( (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**2 &
          -(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**2)&
         +p4*(-1.d0) * ( &                ! Anharmonic contribution
          ( (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**3 & 
           +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**3)&
           ) &
         +p5*(-5.d0/sqrt(2.d0)**5) * ( &                ! Anharmonic contribution
          ( (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**4 & 
           -(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**4)&
           ) &
         +p6*(-6.d0/sqrt(2.d0)**6) * ( &                ! Anharmonic contribution
          ( (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**5 & 
           +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**5)&
           )
   ! including crossed terms  !
  
    forces(i1,i2,1)= forces(i1,i2,1) &
    - 0.5d0*p4x * ( &
       (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))&
        *( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2 +(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 ) &  ! First term       
      +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1)) &
        *( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2 +(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 ) &  ! Second term
      +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1)) &  ! third term
      +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1)) &  ! fourth term
      +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**2 * (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1)) &  ! fifth term
      +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**2 * (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1)) &  ! sixth  term  
      ) &
     -(1.d0/dsqrt(2.d0)**3)*p3x*( &
            (-1)*( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2  +(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 ) &
           -(-1)*( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2  +(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 ) &
           +2*(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))*(u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1)) &
           -2*(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))*(u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1)) &
           +2*(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))*(u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1)) &
           -2*(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))*(u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1)) &
           ) &
     -0.25*p4f*( &
          2*( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1)) &
             *(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 )&
         +2*( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1)) &
             *(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 )&
         +2*( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1)) &
             *(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 )&
         +2*( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1)) &
             *(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 )&
          ) &
     -0.25*p4g*( &
         +4*(u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**3 &
         +4*(u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**3 &
         +4*(u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**3 &
         +4*(u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**3 &
          ) 

      ! end x component    !

! forces, y component !
    forces(i1,i2,2)= forces(i1,i2,2) &
        +p2*(-2.d0/sqrt(2.d0)**2) * ( &                ! Anharmonic contribution
        ( (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))& 
         +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2)))&
          ) &
         +p3*(-3.d0/sqrt(2.d0)**3)*&    ! 3rd order contribution
         ( (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**2 &
          -(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**2)&
         +p4*(-4.d0/sqrt(2.0)**4) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**3 &
          +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**3)&
           ) &
         +p5*(-5.d0/sqrt(2.d0)**5) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**4 &
          -(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**4)&
           ) &
         +p6*(-6.d0/sqrt(2.d0)**6) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**5 &
          +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**5)&
           )
   ! including crossed terms  !
    forces(i1,i2,2)= forces(i1,i2,2) &
     - 0.5d0*p4x * ( &
      (u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))          &
       *( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2   &
         +(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 ) &  ! First term       
      +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2)) &
       *( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2   &
         +(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 ) &  ! Second term
      +(u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2)) &  ! third term
      +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2)) &  ! fourth term
      +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**2 * (u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2)) &  ! fifth term
      +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**2 * (u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2)) &  ! sixth  term  
      ) &
     -(1.d0/dsqrt(2.d0)**3)*p3x*( &
            (-1)*( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2  +(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 ) &
           -(-1)*( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2  +(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 ) &
           +2*(u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))*(u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2)) &
           -2*(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))*(u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2)) &
           +2*(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))*(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2)) &
           -2*(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))*(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2)) &
           ) &
     -0.25*p4f*( &
          2*( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))    * (u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 )&
         +2*( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))    * (u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 )&
         +2*( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))    )&
         +2*( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))    )&
          )&
     -0.25*p4g*( &
         +4*(u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**3 &
         +4*(u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**3 &
         +4*(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**3 &
         +4*(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**3 &
          ) 

! end y component  !

! forces, z component !
    forces(i1,i2,3)= forces(i1,i2,3) &
        +p2*(-2.d0/sqrt(2.d0)**2)*&    ! 3rd order contribution
         ( (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))&
          +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3)))&
        +p3*(-3.d0/sqrt(2.d0)**3)*&    ! 3rd order contribution
         ( (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**2&
          -(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**2)&
         +p4*(-1.d0) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**3 &
          +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**3)&
           ) &
         +p5*(-5.d0/sqrt(2.d0)**5) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**4 &
          -(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**4)&
           ) &
         +p6*(-6.d0/sqrt(2.d0)**6) * ( &                ! Anharmonic contribution
         ( (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**5 &
          +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**5)&
           )
   ! including crossed terms  !
    forces(i1,i2,3)= forces(i1,i2,3) &
    - 0.5d0*p4x * ( &   
       (u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3)) &
        *((u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 +(u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2 ) &  ! First term       
      +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3)) &
        *((u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 +(u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2 ) &  ! Second term
      +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3)) &  ! third term
      +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3)) &  ! fourth term
      +(u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3)) &  ! fifth term
      +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3)) &  ! sixth  term  
      ) &
     -(1.d0/dsqrt(2.d0)**3)*p3x*( &
            (-1)*( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2  +(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 ) &
           -(-1)*( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2  +(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 ) & 
           +2*(u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))*(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3)) &
           -2*(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))*(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3)) &
           +2*(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))*(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3)) &
           -2*(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))*(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3)) &
           ) &
     -0.25*p4f*( &
          2*( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3)) )&
         +2*( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3)) )&
         +2*( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3)) )&
         +2*( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3)) )&
          )&
     -0.25*p4g*( &
         +4*(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**3 &
         +4*(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**3 &
         +4*(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**3 &
         +4*(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**3 &
          ) 

      
    ! end z component
    enddo
    
    ELSE IF (type_cal == "harmx") THEN 
      print*, "skipping 3rd and 4th order terms"
    ELSE
      print*, "what?"
      STOP 3
    ENDIF &
    ANHARMONIC_FORCE_FIELDS 
    
    END DO &
    CONFIGURATIONS_LOOP_2 
        
  print *, ''
  print *, ' ****************************************'
  print *, ' *                                      *'
  print *, ' *   GET THE FREE ENERGY FOR EACH       *'
  print *, ' *           CONFIGURATION              *'
  print *, ' *                                      *'
  print *, ' ****************************************'
  print *, ''

  !--------------------------------------------------------------!
  ! THIS NEEDS TO BE CORRECTED IF EFFECTIVE CHARGES ARE PRESENT  ! 
  ! AND ALSO FOR THE OTHER FORCE FIELDS                           !
  !                                                              !
  do i1 = 1, n_random
    ! Get harmonic contribution of the potential to the free energy
  v_anharmonic = 0.0d0
    call v2_0_harmonic (u_vector(i1,:,:),phi_sc_harmonic,v_anharmonic, nat_sc)
    do i2 = 1, nat_sc
      !!!!!!!!!
      ! CAREFUL
      ! 
      ! We assume that the anharmonic atom is the atom 2
      if (ityp_sc(i2) .eq. 2) then
        dispvect(:) = u_vector(i1,i2,:)
        do i3 = 1, 6
          dispvect(:) = dispvect(:) - u_vector(i1,nn_at(i2,i3),:) / 6.d0
        end do
        

        if (type_cal .eq. 'pdhxx' ) then
        v_anharmonic = v_anharmonic + &
                       b*(dispvect(1)**4 + dispvect(2)**4 + dispvect(3)**4) + &
                       c*(dispvect(1)**2 * dispvect(2)**2 + &
                          dispvect(1)**2 * dispvect(3)**2 + &
                          dispvect(2)**2 * dispvect(3)**2 )


       elseif (type_cal .eq. 'pbtex' ) then
       
       
       v_anharmonic = v_anharmonic &
     -(0.5d0)*p2*(   &
           (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**2 &   ! x slave modes 3order
          +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**2 &
          +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**2 &   ! y slave modes 3order
          +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**2 &
          +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**2 &   ! z slave modes 3order
          +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**2 &
          )&
     -(0.5d0/dsqrt(2.d0))*p3*(   &
           (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**3 &   ! x slave modes 3order
          -(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**3 &
          +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**3 &   ! y slave modes 3order
          -(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**3 &
          +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**3 &   ! z slave modes 3order
          -(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**3 &
          )&
     -(0.5d0/dsqrt(2.d0))*p3x*( &
          (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))      &
           *( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2  &
             +(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 )&
         -(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))      &    
           *( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2  &
             +(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))      &
           *( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 )&
         -(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))      &    
           *( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))      &
           *( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 )&
         -(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))      &    
           *( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 )&
          ) &
     -0.25d0*p4*( &
          (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**4 &   ! x slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**4 & 
         +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**4 &   ! y slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**4 &
         +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**4 &   ! z slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**4 &
         ) &
     -0.25d0*p4x*( &
          (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**2      &
           *( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2  &
             +(u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**2      &    
           *( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2  &
             +(u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**2      &
           *( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**2      &    
           *( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 )&
         +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**2      &
           *( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 )&
         +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**2      &    
           *( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2  &
             +(u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 )&
          )&
     -0.25*p4f*( &
          ( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**2 )&
         +( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**2 * (u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**2 )&
         +( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**2 )&
         +( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**2 )&
         +( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**2 )&
         +( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**2 * (u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**2 )&
          ) &
     -0.25*p4g*( &
          ( (u_vector(i1,slv_idx(i2,1),2)-u_vector(i1,i2,2))**4 + (u_vector(i1,slv_idx(i2,1),3)-u_vector(i1,i2,3))**4 )&
         +( (u_vector(i1,slv_idx(i2,2),2)-u_vector(i1,i2,2))**4 + (u_vector(i1,slv_idx(i2,2),3)-u_vector(i1,i2,3))**4 )&
         +( (u_vector(i1,slv_idx(i2,3),1)-u_vector(i1,i2,1))**4 + (u_vector(i1,slv_idx(i2,3),3)-u_vector(i1,i2,3))**4 )&
         +( (u_vector(i1,slv_idx(i2,4),1)-u_vector(i1,i2,1))**4 + (u_vector(i1,slv_idx(i2,4),3)-u_vector(i1,i2,3))**4 )&
         +( (u_vector(i1,slv_idx(i2,5),1)-u_vector(i1,i2,1))**4 + (u_vector(i1,slv_idx(i2,5),2)-u_vector(i1,i2,2))**4 )&
         +( (u_vector(i1,slv_idx(i2,6),1)-u_vector(i1,i2,1))**4 + (u_vector(i1,slv_idx(i2,6),2)-u_vector(i1,i2,2))**4 )&
          ) &
     -p5/(sqrt(2.d0)**5)*( &
          (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**5 &   ! x slave modes anharmonic
         -(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**5 & 
         +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**5 &   ! y slave modes anharmonic
         -(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**5 &
         +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**5 &   ! z slave modes anharmonic
         -(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**5 &
         ) &
     -p6/(sqrt(2.d0)**6)*( &
          (u_vector(i1,slv_idx(i2,1),1)-u_vector(i1,i2,1))**6 &   ! x slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,2),1)-u_vector(i1,i2,1))**6 & 
         +(u_vector(i1,slv_idx(i2,3),2)-u_vector(i1,i2,2))**6 &   ! y slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,4),2)-u_vector(i1,i2,2))**6 &
         +(u_vector(i1,slv_idx(i2,5),3)-u_vector(i1,i2,3))**6 &   ! z slave modes anharmonic
         +(u_vector(i1,slv_idx(i2,6),3)-u_vector(i1,i2,3))**6 &
         ) 
        end if
        
      end if
    end do
    v(i1) = v_anharmonic
  end do

END SUBROUTINE get_forces_energies



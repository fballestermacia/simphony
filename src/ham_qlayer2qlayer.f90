  subroutine ham_qlayer2qlayer(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), allocatable :: Hij(:, :, :)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     !     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: H00new(Ndim,Ndim)

     ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
     !     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Ndim,Ndim)

     allocate(Hij(Num_wann,Num_wann,-ijmax:ijmax))

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
       
        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*new_ia+ k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(1:Num_wann, 1:Num_wann, inew_ic )&
           =Hij(1:Num_wann, 1:Num_wann, inew_ic)&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

     H00new=0.0d0
     H01new=0.0d0

     ! nslab's principle layer 
     ! H00new
     do i=1,Np
     do j=1,Np
        if (abs(i-j).le.(ijmax)) then
          H00new(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Hij(:,:,j-i)
        endif
     enddo
     enddo

     ! H01new
     do i=1,Np
     do j=Np+1,Np*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(:,:,j-i)
        endif
     enddo
     enddo

   !   do i=1,Ndim
   !   do j=1,Ndim
   !      if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
   !     !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
   !     !stop
   !      endif

   !   enddo
   !   enddo

     deallocate(Hij)

  return
  end subroutine ham_qlayer2qlayer

  subroutine ham_qlayer2qlayer_bak(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), allocatable :: Hij(:, :, :)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     !     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: H00new(Ndim,Ndim)

     ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
     !     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Ndim,Ndim)

     allocate(Hij(-ijmax:ijmax,Num_wann,Num_wann))

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
       
        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*new_ia+ k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

     H00new=0.0d0
     H01new=0.0d0

     ! nslab's principle layer 
     ! H00new
     do i=1,Np
     do j=1,Np
        if (abs(i-j).le.(ijmax)) then
          H00new(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Hij(j-i,:,:)
        endif
     enddo
     enddo

     ! H01new
     do i=1,Np
     do j=Np+1,Np*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(j-i,:,:)
        endif
     enddo
     enddo

     do i=1,Ndim
     do j=1,Ndim
        if(abs(H00new(i,j)-conjg(H00new(j,i))).ge.1e-4)then
       !  write(stdout,*)'there are something wrong with ham_qlayer2qlayer'
       !stop
        endif

     enddo
     enddo

     deallocate(Hij)

  return
  end subroutine ham_qlayer2qlayer_bak



  subroutine ham_qlayer2qlayer_LOTO(k,H00new,H01new)
     ! This subroutine caculates Hamiltonian between
     ! slabs For surface state calculation
     ! Aug/01/2017 by QS Wu
     ! Copyright (c) 2017 QuanSheng Wu. All rights reserved.

      !        July/2/2024 by Francisco M. Ballester
      !          - LO-TO implementation.


     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic
     integer  :: ii,jj,pp,qq, counter, iia

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)
     
     ! wave vector transformed as if it where 3D system for lrange calculation
     real(Dp) :: k3d(3), keps(3) = (/eps12,eps12,eps12/)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers

     complex(dp) :: ratio

     complex(Dp), allocatable :: Hij(:, :, :)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers

     !     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: H00new(Ndim,Ndim)

     ! H01 Hamiltonian between next-nearest neighbour-quintuple-layers
     !     complex(Dp),allocatable,intent(out) :: H01new(:,:)
     complex(Dp),intent(out) :: H01new(Ndim,Ndim)

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     complex(Dp),allocatable :: nac_correction(:, :) 

     real(dp) :: temp1(2), temp2, zag(3), zbg(3), &
                  R1(3), R2(3),R3(3), R12_cross(3),angle_t
     real(dp) :: constant_t, qeq
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     real(dp), allocatable :: pos_cart_ic(:,:)

     !> check if we are exactly at gamma
     logical :: atGamma=.false.

     real(dp) :: nac_q
     real(dp), external :: norm, angle

     allocate(Hij(-ijmax:ijmax,Num_wann,Num_wann))
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))
     allocate(nac_correction(Num_wann, Num_wann))
     allocate(pos_cart_ic(3,  Origin_cell%Num_atoms))
     mat1 = 0d0
     mat2 = 0d0
     nac_correction= 0d0
     atGamma = .false.

     
     
     !>  add loto splitting term
     temp1(1:2)= (/0.0,0.0/)
     
     k3d(:) = k(1)*Cell_defined_by_surface%reciprocal_lattice(1,:) + k(2)*Cell_defined_by_surface%reciprocal_lattice(2,:) ! ESTO NO ME CONVENCE
     k3d = k3d*Cell_defined_by_surface%cell_parameters(1)/(twopi)
     
     if (abs((k3d(1)**2+k3d(2)**2+k3d(3)**2)).le.eps12)then  !> skip k=0
         atGamma=.true.

         qeq = (keps(1)*(Diele_Tensor(1,1)*keps(1)+Diele_Tensor(1,2)*keps(2)+Diele_Tensor(1,3)*keps(3))+    &
            keps(2)*(Diele_Tensor(2,1)*keps(1)+Diele_Tensor(2,2)*keps(2)+Diele_Tensor(2,3)*keps(3))+    &
            keps(3)*(Diele_Tensor(3,1)*keps(1)+Diele_Tensor(3,2)*keps(2)+Diele_Tensor(3,3)*keps(3)))

         constant_t= 2.0d0*4.0d0*Pi/Origin_cell%CellVolume
          do pp = 1,Origin_cell%Num_atoms
            do qq = 1,Origin_cell%Num_atoms
               do ii=1,3
                  zag(ii) = keps(1)*Born_Charge(pp,1,ii) +  keps(2)*Born_Charge(pp,2,ii) + keps(3)*Born_Charge(pp,3,ii)
                  
                  zbg(ii) = keps(1)*Born_Charge(qq,1,ii) +  keps(2)*Born_Charge(qq,2,ii) + keps(3)*Born_Charge(qq,3,ii)

               end do
               
               do ii=1,3
                  do jj=1,3
                     
                     nac_q= constant_t*zag(ii)*zbg(jj)/qeq


                     mat2(3*(pp-1)+ii,3*(qq-1)+jj) = nac_q


                  enddo  ! jj
               enddo  ! ii
            enddo ! qq
          enddo ! pp
      endif
     
   !   nac_correction= 0d0
   !   call long_range_phonon_interaction(0,0,0,k3d(:),.false.,1.0d0,mat1,  &
   !        Cell_defined_by_surface%Atom_position_cart/Cell_defined_by_surface%cell_parameters(1),  &
   !        Born_Charge(:,:,:), Cell_defined_by_surface%reciprocal_lattice*Cell_defined_by_surface%cell_parameters(1)/(twopi), &
   !        Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
     
   !   do ii=1,Num_wann
   !      do jj=1, Num_wann
   !          pp = Origin_cell%spinorbital_to_atom_index(ii)
   !          qq = Origin_cell%spinorbital_to_atom_index(jj)
   !          nac_correction(ii,jj) = (mat1(ii,jj) + mat2(ii,jj)*(108.97077184367376*eV2Hartree)**2)/SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
   !      end do
   !   end do
           
     counter = 0
     !> First check the max number of blocks in the hamiltonian matrix so as to add the long range interaction
     ! counter is the normlization factor, if the supercell is not excesively large, it should coincide with Nrpts, but just to be sure
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)
        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
        
         inew_ic= int(new_ic)
         if (abs(new_ic).le.ijmax)then
            counter = counter + 1
         end if
     end do

     R1=Cell_defined_by_surface%Rua
      R2=Cell_defined_by_surface%Rub
      R3=Cell_defined_by_surface%Ruc

      !> R12_cross=R1xR2
      call cross_product(R1, R2, R12_cross)

      !> angle of R12_cross and R3
      angle_t= angle(R12_cross, R3)
      angle_t= angle_t*pi/180d0

      ratio= Vacuum_thickness_in_Angstrom/cos(angle_t)/norm(R3)
     
 
     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
            pos_cart_ic=0d0
            do iia=1, Origin_cell%Num_atoms
               pos_cart_ic(:, iia)= Origin_cell%Atom_position_cart(:, iia)+ R3*(new_ic-1d0+ratio/2d0)
            enddo
            
            
            mat1 = 0d0
            nac_correction= 0d0
            call long_range_phonon_interaction(0,0,0,k3d(:),.false.,1.0d0,mat1,  &
                  pos_cart_ic/Origin_cell%cell_parameters(1),  &
                  Born_Charge(:,:,:), Cell_defined_by_surface%reciprocal_lattice*Cell_defined_by_surface%cell_parameters(1)/(twopi), &
                  Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
            
            do ii=1,Num_wann
               do jj=1, Num_wann
                     pp = Origin_cell%spinorbital_to_atom_index(ii)
                     qq = Origin_cell%spinorbital_to_atom_index(jj)
                     nac_correction(ii,jj) = (mat1(ii,jj) + mat2(ii,jj)*(108.97077184367376*eV2Hartree)**2)/SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
               end do
            end do

           kdotr=k(1)*new_ia+k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           = Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           + (HmnR(1:Num_wann,1:Num_wann,iR)+ nac_correction(1:Num_wann, 1:Num_wann)/counter)*ratio/ndegen(iR) 
        endif

     enddo

     H00new=0.0d0
     H01new=0.0d0

     ! nslab's principle layer 
     ! H00new
     do i=1,Np
     do j=1,Np
        if (abs(i-j).le.(ijmax)) then
          H00new(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Hij(j-i,:,:)
        endif
     enddo
     enddo

     ! H01new
     do i=1,Np
     do j=Np+1,Np*2
        if (j-i.le.ijmax) then
           H01new(Num_wann*(i-1)+1:Num_wann*i,&
               Num_wann*(j-1-Np)+1:Num_wann*(j-Np))=Hij(j-i,:,:)
        endif
     enddo
     enddo

      ! preserve previous k point for NAC calculation in Gamma so as to not have unnecessary discontinuities
     if ((k3d(1).ne.0.0d0) .or. (k3d(2).ne.0.0d0) .or. (k3d(3).ne.0.0d0)) then
         
        keps(:) = k3d
        do ii=1,3
           if (abs(keps(ii)).le.eps12/1000)then
               keps(ii) = 0.0d0
           end if
        end do  
     end if

     deallocate(mat1)
     deallocate(mat2)
     deallocate(nac_correction)
     deallocate(Hij)
     deallocate(pos_cart_ic)

  return
  end subroutine ham_qlayer2qlayer_LOTO


  subroutine ham_qlayer2qlayer_velocity(k,Vij_x,Vij_y)
     ! This subroutine caculates velocity matrix between
     ! slabs  
     ! 06/Aug/2018 by QS Wu
     ! Copyright (c) 2018 QuanSheng Wu. All rights reserved.

     use para
     implicit none

     ! loop index
     real(dp) :: ia, ib, ic
     integer :: iR,  inew_ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr, r0(3), r1(3)

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Vij_x(-ijmax:ijmax,Num_wann,Num_wann)
     complex(Dp), intent(out) :: Vij_y(-ijmax:ijmax,Num_wann,Num_wann)

     Vij_x=0.0d0; Vij_y= 0d0;
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
        r0= new_ia*Rua_newcell+ new_ib* Rub_newcell

        !> rotate the vector from the original coordinate to the new coordinates
        !> which defined like this: x is along R1', z is along R1'xR2', y is along z x y
        call rotate(r0, r1)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr= k(1)*new_ia+ k(2)*new_ib
           ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

           Vij_x(inew_ic, 1:Num_wann, 1:Num_wann ) &
           = Vij_x(inew_ic, 1:Num_wann, 1:Num_wann ) &
           + zi*r1(1)*HmnR(:,:,iR)*ratio/ndegen(iR)
           Vij_y(inew_ic, 1:Num_wann, 1:Num_wann ) &
           = Vij_y(inew_ic, 1:Num_wann, 1:Num_wann ) &
           + zi*r1(2)*HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
     enddo

     return
  end subroutine ham_qlayer2qlayer_velocity


  subroutine ham_qlayer2qlayer2(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! slabs  
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: iR, inew_ic
     real(dp) :: ia, ib, ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
           kdotr=k(1)*new_ia+k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif

     enddo

  return
  end subroutine ham_qlayer2qlayer2


!   subroutine long_range_phonon_interaction_real_space(R,tau,Ham,zeu,natoms,map2atoms)
!    ! This subroutine computes the rigid-ion (long-range) term 
!    ! of the dynamical matrix in real space, following equation (70)
!    ! of X. Gonze et al, PRB 55. 10355 (1997).
!    ! It is intendet to be used in a given supercell, using the irreducible vectors written
!    ! that are later used for interpolation
!    !
!    ! History
!    !     
!    !       Jan/22/2025 by Francesc Ballester
!    !

!    use para

!    implicit none

!    integer :: natoms, map2atoms(natoms)
!    complex(Dp) :: Ham(3*natoms,3*natoms)
!    real(Dp) &
!          R(3),                               & ! irVec
!          tau(3,natoms),                       & ! atom position
!          zeu(Origin_cell%Num_atoms,3,3)        ! effective charges
   


!    !  local variables

!    integer :: i,j,ii,jj,na,nb                      !  loop integers
!    real(Dp) :: einv(3,3), sqrtedet                 !  dielectric tensor inverse and determinant
!    real(Dp) :: d(3), Delta(3), bigD                !  variables used for computation
!    real(Dp) :: unitConversor, sumalphabeta         !  just for clarity
!    real(Dp), external ::  det3


!    unitConversor = (108.97077184367376*eV2Hartree)**2

!    sqrtedet = SQRT(det3(Diele_Tensor))
!    einv(:,:) = Diele_Tensor(:,:)
   

!    call inv_r(3,einv)

   
!    do na=1, natoms
!       do nb=1, natoms

!          d(:) = R(:) + tau(:,nb) - tau(:,na) ! eq 68

!          if (abs((d(1)**2+d(2)**2+d(3)**2)).ge.eps12)then ! skip the d=0 case
!             do i=1, 3
!                Delta(i) = einv(i,1)*d(1) + einv(i,2)*d(2) + einv(i,3)*d(3)
!             end do
!             bigD = SQRT(Delta(1)*d(1) + Delta(2)*d(2) + Delta(3)*d(3))

!             do i=1, 3
!                do j=1, 3
!                   sumalphabeta = 0.0d0
!                   do ii=1, 3
!                      do jj=1, 3

!                         sumalphabeta = sumalphabeta +                                     &
!                                        (zeu(map2atoms(na),i,ii)*zeu(map2atoms(nb),j,jj)/sqrtedet*   &
!                                        (einv(ii,jj)/bigD**3 - 3*Delta(ii)*Delta(jj)/bigD**5))


!                      end do
!                   end do
!                   Ham(3*(na-1)+i,3*(nb-1)+j) = Ham(3*(na-1)+i,3*(nb-1)+j) + sumalphabeta/SQRT(Atom_Mass(na)*Atom_Mass(nb))*unitConversor
!                end do
!             end do
!          end if

!       end do
!    end do

!   return
!   end subroutine long_range_phonon_interaction_real_space

!   subroutine asr_real_space(nirs,irvecs,tau,Ham,zeu,natoms,map2atoms)
!    ! Tmimimimimimimimi
!    ! omimimimimimimimi
!    ! omimimimimimimimi
!    ! Imimimimimimimimi
!    ! tmimimimimimimimi
!    !
!    ! History
!    !     
!    !       Jan/23/2025 by Francesc Ballester
!    !

!    use para

!    implicit none

!    integer :: nirs, natoms, map2atoms(natoms), irvecs(3,nirs)
!    complex(Dp) :: Ham(3*natoms,3*natoms)
!    real(Dp) &
!          tau(3,natoms),                       & ! atom position
!          zeu(Origin_cell%Num_atoms,3,3)        ! effective charges
   


!    !  local variables

!    integer :: i,j,ii,jj,na,nb,inew_ic                 !   
!    real(Dp) ::  bigR(3),ia,ib,ic,new_ia,new_ib,new_ic              !  
!    complex(Dp) :: sumforasr(3*natoms,3*natoms)

!    !real(Dp), external ::  det3
!    sumforasr = 0.0d0
!    do i=1, nirs
!       ia=irvec(1,i)
!       ib=irvec(2,i)
!       ic=irvec(3,i)


      
!       bigR = ia*Origin_cell%Rua + ib*Origin_cell%Rub + ic*Origin_cell%Ruc
!       bigR = bigR/Origin_cell%cell_parameters(1)
!       call long_range_phonon_interaction_real_space(bigR,tau,sumforasr,zeu,natoms,map2atoms)
!    end do

!    do na=1, natoms
!       do nb=1, natoms
!          if (na.ne.nb) then
!             do i=1, 3
!                do j=1, 3
!                   Ham(3*(na-1)+i,3*(na-1)+j) = Ham(3*(na-1)+i,3*(na-1)+j) + sumforasr(3*(na-1)+i,3*(nb-1)+j)
!                end do
!             end do
!          end if
!       end do
!    end do

   
!   return
!   end subroutine asr_real_space


   subroutine FT_long_range_to_R(r,nkft1,nkft2,nkft3, DofR, tau, zeu, rec_lattice, natoms, map2atoms)
      !
      !
      !
      !
      !
      !

      use para

      !implicit none

      integer :: natoms, map2atoms(natoms),nkft1,nkft2,nkft3
      complex(Dp) :: DofR(3*natoms,3*natoms)
      real(Dp) &
            r(3),                               & ! irVec
            tau(3,natoms),                       & ! atom position
            rec_lattice(3,3),                   & ! rec_lattice = Origin_cell%reciprocal_lattice*Origin_cell%cell_parameters(1)/(twopi)
            zeu(Origin_cell%Num_atoms,3,3)        ! effective charges
      


      !  local variables

      integer :: i,j,ii,jj,na,nb, n1, n2, n3, totalnknumber                   !  loop integers
      real(Dp) :: unitConversor, qdotr, qingrid(3)                !  
      complex(Dp) :: sumoverq(3*natoms,3*natoms), dummy(3*natoms,3*natoms), expqdotr

      unitConversor = (108.97077184367376*eV2Hartree)**2
      
      totalnknumber = 0
      sumoverq = 0.0d0
      do n1=-nkft1,nkft1
         do n2=-nkft2,nkft2
            do n3=-nkft3,nkft3
               qingrid(:) = rec_lattice(:,1)*(n1)/dble(nkft1) + &
                           rec_lattice(:,2)*(n2)/dble(nkft2)  + &
                           rec_lattice(:,3)*(n3)/dble(nkft3)

               if (abs((qingrid(1)**2+qingrid(2)**2+qingrid(3)**2)).ge.eps12)then
                  dummy = 0.0d0
                  qdotr = (n1)/dble(nkft1)*r(1) + (n2)/dble(nkft2)*r(2) + (n3)/dble(nkft3)*r(3)!qingrid(1)*r(1) + qingrid(2)*r(2) + qingrid(3)*r(3)
                  call long_range_phonon_interaction(0,0,0,qingrid,.false.,1.0d0,dummy,tau,zeu,rec_lattice,natoms, map2atoms)
                  expqdotr = cos(2d0*pi*qdotr)-zi*sin(2d0*pi*qdotr)!exp(-zi*qdotr*twopi)
                  sumoverq = sumoverq+dummy*expqdotr
                  totalnknumber = totalnknumber + 1

               end if

            end do
         end do
      end do
      
      DofR = sumoverq/dble(totalnknumber)!*2d0*(twopi)**3/Origin_cell%cellvolume
      
      do ii=1, 3*natoms
         do jj=1, 3*natoms
            i = Origin_cell%spinorbital_to_atom_index(ii)
            j = Origin_cell%spinorbital_to_atom_index(jj)
            
            DofR(ii,jj) = DofR(ii,jj)/SQRT(Atom_Mass(i)*Atom_Mass(j)) 
            !write(*,*) DofR(ii,jj)
            !write(*,*) i,j,SQRT(Atom_Mass(i)*Atom_Mass(j)), Atom_Mass(i), Atom_Mass(j)
         end do
      end do
      




      
      return
   end subroutine FT_long_range_to_R



  subroutine ham_qlayer2qlayer2_LOTO(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! slabs for phonon system with LOTO correction
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     ! Modified and corrected by Francesc Ballester
     ! September 2024

     use para

     implicit none

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic
     integer  :: ii,jj,pp,qq, iia

     ! 

     ! new index used to sign irvec     
     integer  :: inew_ic
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(dp),intent(in) :: k(2)
     real(Dp) :: k3d(3)
     complex(dp) :: ratio

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin
     complex(dp), intent(out) :: Hij(-ijmax:ijmax,Num_wann,Num_wann)

     !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
     complex(Dp),allocatable :: nac_correction(:, :) 

     real(dp) ::  constant_t, R(3)
     complex(dp), allocatable :: mat1(:, :)
     complex(dp), allocatable :: mat2(:, :)
     real(dp), allocatable :: pos_cart_ic(:, :)
     integer :: counter
 
     !> k times Born charge
     logical :: atGamma
     real(dp) :: qeq, temp1(3), temp2, zag(3), zbg(3),  keps(3) = (/eps12,eps12,0.0d0/), &
                  R1(3), R2(3),R3(3), R12_cross(3),angle_t
     complex(dp) :: nac_q
     real(dp), external :: norm, angle
     

     
     allocate(mat1(Num_wann, Num_wann))
     allocate(mat2(Num_wann, Num_wann))
     allocate(nac_correction(Num_wann, Num_wann))
     allocate(pos_cart_ic(3,  Origin_cell%Num_atoms))

     mat1 = 0d0
     mat2 = 0d0
     nac_correction= 0d0
     atGamma = .false.

     
     
     !>  add loto splitting term
     temp1(1:2)= (/0.0,0.0/)
     
     k3d(:) = k(1)*Cell_defined_by_surface%reciprocal_lattice(1,:) + k(2)*Cell_defined_by_surface%reciprocal_lattice(2,:) ! ESTO NO ME CONVENCE
     
     k3d = k3d!*Cell_defined_by_surface%cell_parameters(1)/(twopi)
     
     if (abs((k3d(1)**2+k3d(2)**2+k3d(3)**2)).le.eps12)then  !> skip k=0
         atGamma=.true.

         qeq = (keps(1)*(Diele_Tensor(1,1)*keps(1)+Diele_Tensor(1,2)*keps(2)+Diele_Tensor(1,3)*keps(3))+    &
            keps(2)*(Diele_Tensor(2,1)*keps(1)+Diele_Tensor(2,2)*keps(2)+Diele_Tensor(2,3)*keps(3))+    &
            keps(3)*(Diele_Tensor(3,1)*keps(1)+Diele_Tensor(3,2)*keps(2)+Diele_Tensor(3,3)*keps(3)))

         constant_t= 2.0d0*4.0d0*Pi/Origin_cell%CellVolume
          do pp = 1,Origin_cell%Num_atoms
            do qq = 1,Origin_cell%Num_atoms
               do ii=1,3
                  zag(ii) = keps(1)*Born_Charge(pp,1,ii) +  keps(2)*Born_Charge(pp,2,ii) + keps(3)*Born_Charge(pp,3,ii)
                  
                  zbg(ii) = keps(1)*Born_Charge(qq,1,ii) +  keps(2)*Born_Charge(qq,2,ii) + keps(3)*Born_Charge(qq,3,ii)

               end do
               
               do ii=1,3
                  do jj=1,3
                     
                     nac_q= constant_t*zag(ii)*zbg(jj)/qeq


                     mat2(3*(pp-1)+ii,3*(qq-1)+jj) = nac_q


                  enddo  ! jj
               enddo  ! ii
            enddo ! qq
          enddo ! pp
      endif

      
   !   mat1 = 0d0
   !   nac_correction= 0d0
   !   call long_range_phonon_interaction(0,0,0,k3d(:),.false.,1.0d0,mat1,  &
   !        Cell_defined_by_surface%Atom_position_cart/Cell_defined_by_surface%cell_parameters(1),  &
   !        Born_Charge(:,:,:), Cell_defined_by_surface%reciprocal_lattice*Cell_defined_by_surface%cell_parameters(1)/(twopi), &
   !        Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
     
   !   do ii=1,Num_wann
   !      do jj=1, Num_wann
   !          pp = Origin_cell%spinorbital_to_atom_index(ii)
   !          qq = Origin_cell%spinorbital_to_atom_index(jj)
   !          nac_correction(ii,jj) = (mat1(ii,jj) + mat2(ii,jj)*(108.97077184367376*eV2Hartree)**2)/SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
   !      end do
   !   end do
      
     Hij=0.0d0
   !    R1=Cell_defined_by_surface%Rua
   !    R2=Cell_defined_by_surface%Rub
   !    R3=Cell_defined_by_surface%Ruc

   !    !> R12_cross=R1xR2
   !    call cross_product(R1, R2, R12_cross)

   !    !> angle of R12_cross and R3
   !    angle_t= angle(R12_cross, R3)
   !    angle_t= angle_t*pi/180d0

   !    ratio= Vacuum_thickness_in_Angstrom/cos(angle_t)/norm(R3)
   !   counter = 0
     !> First check the max number of blocks in the hamiltonian matrix so as to add the long range interaction
     ! counter is the normalization factor, if the supercell is not excesively large, it should coincide with Nrpts, but just to be sure
   !   do iR=1,Nrpts
   !      ia=irvec(1,iR)
   !      ib=irvec(2,iR)
   !      ic=irvec(3,iR)
   !      !> new lattice
   !      call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
        
   !       inew_ic= int(new_ic)
   !       if (abs(new_ic).le.ijmax)then
   !          counter = counter + 1
   !       end if
   !   end do


     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)


        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
        
        inew_ic= int(new_ic)
        if (abs(new_ic).le.ijmax)then
            ! pos_cart_ic=0d0
            ! do iia=1, Origin_cell%Num_atoms
            !    pos_cart_ic(:, iia)= Origin_cell%Atom_position_cart(:, iia)+ R3*(new_ic-1d0+ratio/2d0)
            ! enddo

            mat1 = 0d0
            !R =(/new_ia, new_ib, new_ic/)
            ! R(1) = new_ia*Cell_defined_by_surface%lattice(1,1) + new_ib*Cell_defined_by_surface%lattice(1,2) + new_ic*Cell_defined_by_surface%lattice(1,3)
            ! R(2) = new_ia*Cell_defined_by_surface%lattice(2,1) + new_ib*Cell_defined_by_surface%lattice(2,2) + new_ic*Cell_defined_by_surface%lattice(2,3)
            ! R(3) = new_ia*Cell_defined_by_surface%lattice(3,1) + new_ib*Cell_defined_by_surface%lattice(3,2) + new_ic*Cell_defined_by_surface%lattice(3,3)


            ! R = new_ia*Origin_cell%Rua + new_ib*Origin_cell%Rub + new_ic*Origin_cell%Ruc
            ! R = R/Origin_cell%cell_parameters(1)

            ! call FT_long_range_to_R(R,3,3,3,mat1,     &
            !                         Cell_defined_by_surface%Atom_position_cart/Cell_defined_by_surface%cell_parameters(1),  &
            !                         Born_Charge(:,:,:), Cell_defined_by_surface%reciprocal_lattice*Cell_defined_by_surface%cell_parameters(1)/(twopi), &
            !                         Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
           


            ! call long_range_phonon_interaction_real_space(R, Origin_cell%Atom_position_cart(:, :)/Origin_cell%cell_parameters(1), mat1, Born_Charge(:,:,:), &
            !       Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
            
            ! if ((ia*ia+ib*ib+ic*ic).eq.0)then
            !    call asr_real_space(Nrpts,irvec,Origin_cell%Atom_position_cart(:, :)/Origin_cell%cell_parameters(1), mat1, Born_Charge(:,:,:), &
            !       Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
            ! end if


            ! mat1 = 0d0
            ! nac_correction= 0d0
            ! call long_range_phonon_interaction(0,0,0,k3d(:),.false.,1.0d0,mat1,  &
            !       pos_cart_ic/Origin_cell%cell_parameters(1),  &
            !       Born_Charge(:,:,:), Cell_defined_by_surface%reciprocal_lattice*Cell_defined_by_surface%cell_parameters(1)/(twopi), &
            !       Origin_cell%Num_atoms, Origin_cell%spinorbital_to_atom_index(::3))
            
            ! do ii=1,Num_wann
            !    do jj=1, Num_wann
            !          pp = Origin_cell%spinorbital_to_atom_index(ii)
            !          qq = Origin_cell%spinorbital_to_atom_index(jj)
            !          nac_correction(ii,jj) = (mat1(ii,jj) + mat2(ii,jj)*(108.97077184367376*eV2Hartree)**2)/SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
            !    end do
            ! end do

           kdotr=k(1)*new_ia+k(2)*new_ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)
            

            Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
            = Hij(inew_ic, 1:Num_wann, 1:Num_wann )&
            + (HmnR(:,:,iR))*ratio/ndegen(iR)  !+nac_correction(1:Num_wann, 1:Num_wann)/counter
         endif

     enddo

     ! preserve previous k point for NAC calculation in Gamma so as to not have unnecessary discontinuities
     if ((k3d(1).ne.0.0d0) .or. (k3d(2).ne.0.0d0) .or. (k3d(3).ne.0.0d0)) then
         
        keps(:) = k3d
        do ii=1,3
           if (abs(keps(ii)).le.eps12/1000)then
               keps(ii) = 0.0d0
           end if
        end do  
     end if
     

     
     deallocate(mat1)
     deallocate(mat2)
     deallocate(nac_correction)
     deallocate(pos_cart_ic)
  return
  end subroutine ham_qlayer2qlayer2_LOTO

  subroutine ham_qlayer2qlayerribbon(k,Hij)
     ! This subroutine caculates Hamiltonian between
     ! ribbon
     ! 4/23/2010 by QS Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.

     use para

     implicit none

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

     complex(dp) :: ratio

     complex(Dp), intent(out) :: Hij(-ijmax:ijmax, &
        -ijmax:ijmax,Num_wann,Num_wann)

     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)


        inew_ia= int(new_ia)
        inew_ib= int(new_ib)
        if (abs(new_ia).le.ijmax)then
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ic
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
        endif

     enddo

     return
  end subroutine ham_qlayer2qlayerribbon


  subroutine latticetransform(a, b, c, x, y, z)
     !> use Umatrix to get the new representation of a vector in new basis
     !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
     use para
     implicit none

     real(dp), intent(in)  :: a, b, c
     real(dp), intent(out) :: x, y, z

     real(dp), allocatable :: Uinv(:, :)

     allocate(Uinv(3, 3))
     Uinv= Umatrix

     call inv_r(3, Uinv)

     x= a*Uinv(1, 1)+ b*Uinv(2, 1)+ c*Uinv(3, 1)
     y= a*Uinv(1, 2)+ b*Uinv(2, 2)+ c*Uinv(3, 2)
     z= a*Uinv(1, 3)+ b*Uinv(2, 3)+ c*Uinv(3, 3)

     return
  end subroutine latticetransform
  
    subroutine transform_r(v0, x, y, z)
     !> use Umatrix to get the new representation of a vector in new basis
     !> R= a*R1+b*R2+c*R3= x*R1'+y*R2'+z*R3'
     use para
     implicit none

     real(dp), intent(in)  :: v0(3)
     real(dp), intent(out) :: x, y, z

     real(dp), allocatable :: Uinv(:, :)

     allocate(Uinv(3, 3))
     Uinv= Umatrix

     call inv_r(3, Uinv)

     x= v0(1)*Uinv(1, 1)+ v0(2)*Uinv(2, 1)+ v0(3)*Uinv(3, 1)
     y= v0(1)*Uinv(1, 2)+ v0(2)*Uinv(2, 2)+ v0(3)*Uinv(3, 2)
     z= v0(1)*Uinv(1, 3)+ v0(2)*Uinv(2, 3)+ v0(3)*Uinv(3, 3)

     return
  end subroutine 

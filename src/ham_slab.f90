  subroutine ham_slab(k,Hamk_slab)
     ! This subroutine is used to caculate Hamiltonian for 
     ! slab system. 
     ! 
     ! History  
     !        4/18/2010 by Quansheng Wu
     ! Copyright (c) 2010 QuanSheng Wu. All rights reserved.
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2, qq, pp, ii, jj

     ! wave vector in 2d
     real(Dp), intent(in) :: k(2)      

     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

     ! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :)

     ! for LO-TO correction
     real(Dp) :: k3d(3)
     real(dp) :: temp1(3), temp2, keps(2)
     real(dp) ::  constant_t
     complex(dp), allocatable :: mat1(:, :), mat2(:,:)
 
     !> k times Born charge
     real(dp), allocatable :: kBorn(:, :)
     real(dp) :: nac_q

     allocate( Hij(-ijmax:ijmax,Num_wann,Num_wann))
     allocate( mat1(Num_wann,Num_wann))
     allocate( mat2(Num_wann,Num_wann))
     allocate(kBorn(Origin_cell%Num_atoms, 3))
     !mat1 = 0.0d0

     call ham_qlayer2qlayer2(k,Hij)

     Hamk_slab=0.0d0 
     ! i1 column index
     do i1=1, nslab
        ! i2 row index
        do i2=1, nslab
          if (abs(i2-i1).le.ijmax)then
            Hamk_slab((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                      (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann )&
            = Hij(i1-i2,1:Num_wann,1:Num_wann) 
          endif 
        enddo ! i2
     enddo ! i1
     write(*,*) Num_wann!,Num_atoms_slab
     !> deal with phonon system
     if (LOTO_correction) then
         k3d(:) = k(1)*Umatrix(1,:) + k(2)*Umatrix(2,:) ! ESTO NO ME CONVENCE

         !call ham_qlayer2qlayer2(k,Hij)
         !call ham_qlayer2qlayer2_LOTO(k,Hij)
         call long_range_phonon_interaction(0,0,0,k3d(:),.false.,1.0d0,mat1)
         
         temp1(1:2)= (/0.0,0.0/)
         keps(1:2) = (/eps12,eps12/)
         temp2= 0.0
         if (abs((k(1)**2+k(2)**2)).le.eps12) then  !> skip k=0
            !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
            !k3d(1:3) = keps(1)*Umatrix(1,:) + keps(2)*Umatrix(2,:)
            
            do qq= 1, 3
               temp1(qq)= keps(1)*Diele_Tensor(qq, 1)+keps(2)*Diele_Tensor(qq, 2)!+k3d(3)*Diele_Tensor(qq,3)
            enddo
            temp2= keps(1)*temp1(1)+ keps(2)*temp1(2) !+ k3d(3)*temp1(3)
            constant_t= 2.0d0*4.0d0*Pi/(temp2*Origin_cell%CellVolume)
            
            do ii=1, Origin_cell%Num_atoms
               do pp=1, 3
                  kBorn(ii, pp)=  keps(1)*Born_Charge(ii,1,pp)+keps(2)*Born_Charge(ii,2,pp)!+k3d(3)*Born_Charge(ii,3,pp)
               enddo
            enddo
            
            do ii= 1,Origin_cell%Num_atoms
               do pp= 1, 3
                  do jj= 1, Origin_cell%Num_atoms
                     do qq= 1,3
                        nac_q = zzero
                        nac_q= kBorn(jj, qq)*kBorn(ii, pp)*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
      
                        mat1(3*(ii-1)+pp,3*(jj-1)+qq) =mat1(3*(ii-1)+pp,3*(jj-1)+qq) + nac_q*(108.97077184367376*eV2Hartree)**2
                        
                     enddo  ! qq
                  enddo  ! jj
               enddo ! pp
            enddo  ! ii

         end if
      endif

     
     
     
     ! check hermitcity

     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
         !write(stdout,*)'there are something wrong with Hamk_slab'
         !stop
        endif 
     enddo
     enddo

     deallocate( Hij)
     deallocate( mat1)
     deallocate( kBorn)
  return
  end subroutine ham_slab


subroutine ham_slab_sparseHR(nnzmax, k, acoo,jcoo,icoo)
   !> Calculate slab hamiltonian with the sparse hr format
   !> Dec 17 2018 EPFL
   !> QuanSheng Wu (wuquansheng@gmail.com)
   use para
   implicit none

   !> input: nnzmax is the maximum number of non-zeros entries 
   !> output: nnzmax is the number of non-zeros entries of acoo
   integer, intent(inout) :: nnzmax
   real(dp), intent(in) :: k(3)

   !> output hamiltonian stored as COO sparse matrix format
   complex(dp), intent(inout) :: acoo(nnzmax)
   integer, intent(inout) :: jcoo(nnzmax)
   integer, intent(inout) :: icoo(nnzmax)

   ! loop index
   integer :: i1, i2, ncoo, iR, ims

   ! index used to sign irvec
   real(dp) :: ia,ib,ic

   integer :: inew_ic

   !> new index used to sign irvec
   real(dp) :: new_ia,new_ib,new_ic

   !> wave vector k times lattice vector R
   real(dp) :: kdotr
   complex(dp) :: ratio, tmp

   acoo=zzero
   ncoo=0
   tmp=0d0

   ! i1 column index, sweep over slab along the third vectors in the SURFACE card
   do i1=1, Nslab
      ! i2 row index
      do i2=1, Nslab
         if (abs(i2-i1)> ijmax) cycle

         !> sum over R points to get H(k1, k2)
         do ims=1,splen
            ia= hirv(1, ims)
            ib= hirv(2, ims)
            ic= hirv(3, ims)

            !> new lattice
            call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

            !> Fourier transform confined on the surface plane
            inew_ic= int(new_ic)
            if (inew_ic /= (i2-i1)) cycle

            !> exp(i k.R)
            kdotr= k(1)*new_ia+ k(2)*new_ib
            ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)

            tmp=hacoo(ims)*ratio/ndegen(iR)
            if(abs(tmp) > 1e-6) ncoo=ncoo+1
            icoo(ncoo)= hicoo(ims)+ (i1-1)*Num_wann
            jcoo(ncoo)= hjcoo(ims)+ (i2-1)*Num_wann
            acoo(ncoo)= acoo(ncoo)+ tmp
         enddo ! iR

      enddo ! i2
   enddo ! i1

   if (ncoo>nnzmax) STOP ' ERROR: please increase nnzmax in the subroutine ham_slab_sparseHR'

   nnzmax= ncoo

   return
end subroutine ham_slab_sparseHR

  subroutine ham_slab_parallel_B(k,Hamk_slab)
     ! This subroutine is used to caculate Hamiltonian for 
     ! slab system . 
     !> for slab with in-plane magnetic field
     !> the magnetic phase are chosen like this
     !> phi_ij= alpha*[By*(xj-xi)*(zi+zj)-Bx*(yj-yi)*(zi+zj)] 
     !> x, z are in unit of Angstrom, Bx, By are in unit of Tesla
     !> History :
     !        9/21/2015 by Quansheng Wu @ETH Zurich
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2

     ! wave vector in 2d
     real(Dp), intent(inout) :: k(2)      

     ! loop index
     integer :: iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic
     integer :: ia1, ia2

     integer :: istart1, istart2
     integer :: iend1, iend2

     integer :: inew_ic

     !> nwann= Num_wann/2
     integer :: nwann
     
     integer, allocatable :: orbital_start(:)

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  
     real(dp) :: phase
     complex(dp) :: ratio
     complex(dp) :: fac

     real(dp) :: Rp1(3)
     real(dp) :: Rp2(3)
     real(dp) :: R1(3)
     real(dp) :: R2(3)
     real(dp) :: Ri(3)
     real(dp) :: Rj(3)
     real(dp) :: tau1(3)
     real(dp) :: tau2(3)


     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

     nwann= Num_wann/2
     allocate( orbital_start(Origin_cell%Num_atoms+ 1))
     orbital_start= 0
     orbital_start(1)= 1
     do i1=1, Origin_cell%Num_atoms
        orbital_start(i1+1)= orbital_start(i1)+ Origin_cell%nprojs(i1)
     enddo

     Hamk_slab=0.0d0 
     ! i1 column index
     do i1=1, Nslab
        ! i2 row index
        do i2=1, Nslab
           !> sum over R points to get H(k1, k2)
           do iR=1, Nrpts
              ia=irvec(1,iR)
              ib=irvec(2,iR)
              ic=irvec(3,iR)
      
              !> new lattice
              call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      
              inew_ic= int(new_ic)
              if (inew_ic /= (i2-i1)) cycle
      
              !> exp(i k.R)
              kdotr= k(1)*new_ia+ k(2)*new_ib
              ratio= cos(2d0*pi*kdotr)+ zi*sin(2d0*pi*kdotr)
      
              R1= (i1-1)*Ruc_new
              R2= new_ia*Rua_new+ new_ib*Rub_new+ (i2-1)*Ruc_new
      
              do ia1=1, Origin_cell%Num_atoms
              do ia2=1, Origin_cell%Num_atoms
                 R1= Origin_cell%Atom_position_cart(:, ia1)
                 R2= Origin_cell%Atom_position_cart(:, ia2)
                 call rotate(R1, tau1)
                 call rotate(R2, tau2)
      
                
                 Ri= Rp1+ tau1
                 Rj= Rp2+ tau2
      
                 phase= alpha*By*(Rj(3)+Ri(3))*(Rj(1)-Ri(1))  &
                      - alpha*Bx*(Rj(3)+Ri(3))*(Rj(2)-Ri(2))
                 fac= cos(phase)+ zi*sin(phase)
      
                !write(*, '(a, 4i5   )') 'iR, ia ib ic', ir, ia, ib, ic
                !write(*, '(a, 4f10.5)') '            ', new_ia, new_ib, new_ic
                !write(*, '(a, 3f10.5)') 'Ri', Ri
                !write(*, '(a, 3f10.5)') 'Rj', Rj
                !write(*, '(a, 3f10.5)') 'phase', phase
      
                 istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                 iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                 istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                 iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                 
                 Hamk_slab( istart1:iend1, istart2:iend2) &
                 = Hamk_slab( istart1:iend1, istart2:iend2) &
                 + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                 istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                 !> there is soc term in the hr file
                 if (soc>0) then
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2)
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1)
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
      
                    istart1= (i2-1)*Num_wann+ orbital_start(ia1) + Nwann
                    iend1= (i2-1)*Num_wann+ orbital_start(ia1+1)- 1 + Nwann 
                    istart2= (i1-1)*Num_wann+ orbital_start(ia2) + Nwann
                    iend2= (i1-1)*Num_wann+ orbital_start(ia2+1)- 1 + Nwann
                    
                    Hamk_slab( istart1:iend1, istart2:iend2) &
                    = Hamk_slab( istart1:iend1, istart2:iend2) &
                    + HmnR( istart1- (i2-1)*Num_wann:iend1- (i2-1)*Num_wann, &
                    istart2- (i1-1)*Num_wann:iend2- (i1-1)*Num_wann, iR)*ratio/ndegen(iR)* fac
                 endif ! soc
              enddo ! ia2
              enddo ! ia1
           enddo ! iR
        enddo ! i2
     enddo ! i1

 

     !> check hermitcity
     do i1=1,nslab*Num_wann
     do i2=1,nslab*Num_wann
        if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
          write(stdout,*)'there is something wrong with Hamk_slab'
          stop
        endif 
     enddo
     enddo

  return
  end subroutine ham_slab_parallel_B


  subroutine ham_slab(k,Hamk_slab)
     ! This subroutine is used to caculate Hamiltonian for 
     ! slab system. 
     ! 
     ! History  
     !        4/18/2010 by Quansheng Wu
     !       12/09/2024 modified by Francesc Ballester to include LO-TO splitting
  
     use para
     implicit none

     ! loop index  
     integer :: i1, i2, qq, pp, ii, jj, mm, nn

     ! wave vector in 2d
     real(Dp), intent(in) :: k(2)      

     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::Hamk_slab(Num_wann*nslab,Num_wann*nslab) 

     ! the factor 2 is induced by spin
     complex(Dp), allocatable :: Hij(:, :, :)

     ! for LO-TO correction
     real(Dp) :: k3d(3), R1(3), R2(3), R3(3), R12_cross(3),R3_slab(3)
     real(dp) :: temp1(3), temp2, zag(3), zbg(3), keps(3) = (/eps12,eps12,0.0d0/)
     real(dp) ::  constant_t, ratio, angle_t, volume_slab, R(3)
     integer ::  Num_atoms_slab, ia, ib, ic, iR
     complex(dp) :: mat1(Num_wann,Num_wann)
     real(dp), external :: norm, angle
     integer, allocatable :: maptoprimitive(:), orbitaltoatom(:)
     real(dp), allocatable :: pos_cart(:,:)
 
     !> k times Born charge
     real(dp) :: qeq
     complex(dp) :: nac_q

     !> check if we are exactly at gamma
      logical :: atGamma=.false.

     allocate( Hij(-ijmax:ijmax,Num_wann,Num_wann))
     
     !mat1 = 0.0d0
     if (LOTO_correction) then
      call ham_qlayer2qlayer2_LOTO(k,Hij)
     else
      call ham_qlayer2qlayer2(k,Hij)
     end if



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
     
      
      !> DEBUG
      ! if ((k(1).eq.0.5d0).and.(k(2).eq.0.5d0))then
      !    write(*,*) 'Writing Hamiltonian'
      !    outfileindex= outfileindex+ 1
      !    open(unit=outfileindex, file='HamSlabatK_small.dat')
      !    do i1=1, Num_wann*nslab
      !       do i2=1, Num_wann*nslab
      !          write(outfileindex, *) REALPART(Hamk_slab(i1,i2)), IMAGPART(Hamk_slab(i1,i2))
      !       end do
      !    end do
      !    write(outfileindex , *)''
      !    close(outfileindex)
      ! end if
     
     
     ! check hermitcity

   !   do i1=1,nslab*Num_wann
   !   do i2=1,nslab*Num_wann
   !      if(abs(Hamk_slab(i1,i2)-conjg(Hamk_slab(i2,i1))).ge.1e-6)then
   !       write(stdout,*)'there is something wrong with Hamk_slab'
   !       !stop
   !      endif 
   !   enddo
   !   enddo

     deallocate( Hij)
     
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




subroutine make_translational_invariant(Ham)
     ! This subroutine is used to apply the ASR for
     ! slab systems by projecting out the non-traslational-invariant part. 
     !
     ! We use the scheme based on Eq. (81) of X. Gonze et al, PRB 50. 13035 (1994) 
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, m1, m2,islab,jslab, lamd, i, j, na, nb, nasr

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(inout) ::Ham(Num_wann*nslab,Num_wann*nslab) 

     complex(Dp), allocatable :: igHig(:,:)


     if (.not. allocated(ntiH)) then
      allocate(ntiH(Num_wann*nslab,Num_wann*nslab))
      allocate(igHig(Num_wann*nslab,Num_wann*nslab))
      ntiH = 0.0d0
      igHig = 0.0d0

      call Ham_slab_realSpace((/0.0d0,0.0d0/), igHig, 1)
      ! call ham_slab((/0.0d0,0.0d0/), igHig)
      ! call ASR_R_G(igHig)

      if (LOTO_correction)then
         LOTO_correction = .false.
         call ham_slab((/0.0d0,0.0d0/), ntiH)
         LOTO_correction = .true.
      else
         call ham_slab((/0.0d0,0.0d0/), ntiH)
      endif

      ntiH = ntiH - igHig
      deallocate(igHig)
     endif

      ! na = 3*Origin_cell%Num_atoms

      ! do islab=1, Nslab
      !    do jslab=1, Nslab
      !       do kappa=1, Origin_cell%Num_atoms
      !          do kappap=1, Origin_cell%Num_atoms
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
      !             sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
      !          enddo
      !       enddo
      !    enddo
      ! enddo



      Ham = Ham - ntiH


      ! do islab=1, Nslab
      !    do jslab=1, Nslab
      !       do kappa=1, Origin_cell%Num_atoms
      !          do kappap=1, Origin_cell%Num_atoms
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
      !             sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
      !          enddo
      !       enddo
      !    enddo
      ! enddo


end subroutine make_translational_invariant


subroutine make_translational_invariant_ribbon(Ham)
     ! This subroutine is used to apply the ASR for
     ! slab systems by projecting out the non-traslational-invariant part. 
     !
     ! We use the scheme based on Eq. (81) of X. Gonze et al, PRB 50. 13035 (1994) 
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, m1, m2,islab,jslab, lamd, i, j, na, nb, nasr

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(inout) ::Ham(Num_wann*nslab1*nslab2,Num_wann*nslab1*nslab2) 

     complex(Dp), allocatable :: igHig(:,:)


     if (.not. allocated(ntiH)) then
      allocate(ntiH(Num_wann*nslab1*nslab2,Num_wann*nslab1*nslab2))
      allocate(igHig(Num_wann*nslab1*nslab2,Num_wann*nslab1*nslab2))
      ntiH = 0.0d0
      igHig = 0.0d0

      ! call Ham_slab_realSpace((/0.0d0,0.0d0/), igHig, 1)
      call ham_ribbon((/0.0d0,0.0d0/), igHig)
      call ASR_R_G_ribbon(igHig)

      call ham_ribbon((/0.0d0,0.0d0/), ntiH)

      ntiH = ntiH - igHig
      deallocate(igHig)
     endif

      ! na = 3*Origin_cell%Num_atoms

      ! do islab=1, Nslab
      !    do jslab=1, Nslab
      !       do kappa=1, Origin_cell%Num_atoms
      !          do kappap=1, Origin_cell%Num_atoms
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
      !             sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
      !          enddo
      !       enddo
      !    enddo
      ! enddo



      Ham = Ham - ntiH


      ! do islab=1, Nslab
      !    do jslab=1, Nslab
      !       do kappa=1, Origin_cell%Num_atoms
      !          do kappap=1, Origin_cell%Num_atoms
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
      !             Ham(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
      !             sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
      !          enddo
      !       enddo
      !    enddo
      ! enddo


end subroutine make_translational_invariant_ribbon


subroutine Ham_slab_realSpace(k,slabham, Gproj) 
     ! This subroutine generates a set of slab Hamiltonians
     ! in real space
     ! prior to Fourier interpolation onto reciprocal space
     ! 19/06/2026 F Ballester
     
     
     use para

     implicit none

     real(dp), intent(in):: k(2)
     complex(dp), intent(inout) :: slabham(nslab*Num_wann,nslab*Num_wann)
     integer, intent(in) :: Gproj

     ! loop index
     integer :: iR, inew_ic, i1, i2, iRsearch
     real(dp) :: ia, ib, ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic, kdotr, ndegen_slab(nrpts)

     ! eigenvalues for debugging
     real(dp) :: eigen(nslab*Num_wann)
     complex(dp) :: ratio

     ! intermediate array so the math checks out
     complex(dp) :: HijR(-ijmax:ijmax,Num_wann,Num_wann,nrpts),&
                     HmnR_SLAB(Num_wann*Nslab,Num_wann*Nslab,nrpts), &
                     irvec_SLAB(2,nrpts)


     real(Dp) :: tx(nslab*Num_wann), ty(nslab*Num_wann), tz(nslab*Num_wann)

      HijR = 0.0d0
      HmnR_SLAB = 0.0d0
      irvec_SLAB = 0.0d0
      ndegen_slab = 1.0d0
      eigen = 0.0d0
      slabham = 0.0d0

      do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)

        irvec_SLAB(1,iR) = new_ia
        irvec_SLAB(2,iR) = new_ib
        inew_ic= int(new_ic)
         ! print *, irvec(:, iR), new_ia, new_ib, new_ic

         HijR(inew_ic, 1:Num_wann, 1:Num_wann,iR)&
         =HijR(inew_ic, 1:Num_wann, 1:Num_wann,iR)&
         +HmnR(:,:,iR)!*ratio/ndegen(iR)
      enddo

      do iR=1,Nrpts
         ! i1 column index
         do i1=1, nslab
            ! i2 row index
            do i2=1, nslab
               if (abs(i2-i1).le.ijmax)then
                  HmnR_SLAB((i2-1)*Num_wann+1:(i2-1)*Num_wann+Num_wann,&
                           (i1-1)*Num_wann+1:(i1-1)*Num_wann+Num_wann, iR )&
                  = HijR(i1-i2,1:Num_wann,1:Num_wann,iR) 
               endif 
            enddo ! i2
         enddo ! i1
         if (Gproj.eq.1) call ASR_R_G(HmnR_SLAB(:,:,iR))
     enddo


      do iR=1,Nrpts
           kdotr=k(1)*irvec_SLAB(1,iR)+k(2)*irvec_SLAB(2,iR)
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           slabham(:,:) =slabham(:,:)&
           +HmnR_SLAB(:,:,iR)*ratio/ndegen_slab(iR)
      enddo


   !   call ASR_R_G(slabham)
   !   tx=0.0d0
   !   tx(1::3) = 1.0d0
   !   ty=0.0d0
   !   ty(2::3) = 1.0d0
   !   tz=0.0d0
   !   tz(3::3) = 1.0d0
   ! !   tx = MATMUL(slabham,tx)
   !   ty = MATMUL(slabham,ty)
   !   tz = MATMUL(slabham,tz)
   !   print *,'tx'
   !   print *, tx
   !   print *,'ty'
   !   print *, ty
   !   print *,'tz'
   !   print *, tz
  

      
   !    print *, 'start'
   !   do i1 =1, nslab*num_wann
   !    print *, slabham(i1,:) 
   !   enddo
   !   print *, 'end'

   !    call eigensystem_c('V', 'L', Num_wann*Nslab,  slabham , eigen)   

   !    print *, 'start'
   !    print *, eigen
   !   print *, 'end'

end subroutine Ham_slab_realSpace


subroutine ASR_R_G(Hs)
     ! This subroutine projects out the translational part
     ! of the slab Hamiltonian in real space 
     ! prior to Fourier interpolation onto reciprocal space
     ! 22/06/2026 F Ballester
     
     
     use para

     implicit none

     complex(dp), intent(inout) :: Hs(nslab*Num_wann,nslab*Num_wann)

     integer :: a,b, islab, jslab, kappa, na, kappap
     
     real(dp) :: G(nslab*Num_wann,nslab*Num_wann), tx(nslab*Num_wann), ty(nslab*Num_wann), tz(nslab*Num_wann), eye(nslab*Num_wann,nslab*Num_wann)

     G = 0.0d0
     eye = 0.0d0

     tx=0.0d0
     tx(1::3) = 1.0d0
     ty=0.0d0
     ty(2::3) = 1.0d0
     tz=0.0d0
     tz(3::3) = 1.0d0

     do a=1, nslab*Num_wann
      eye(a,a) = 1.0d0
      do b=1, nslab*Num_wann
         G(a,b) =3.0d0/dble(nslab*Num_wann)*&
         (tx(a)*tx(b)+ty(a)*ty(b)+tz(a)*tz(b))
      enddo
     enddo
   !    print *, 'start G'
   !   do a =1, nslab*num_wann
   !    print *, G(a,:) 
   !   enddo
   !   print *, 'end G'

     eye(:,:) = eye(:,:) - G(:,:)



      na = 3*Origin_cell%Num_atoms

      do islab=1, Nslab
         do jslab=1, Nslab
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Hs(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
                  Hs(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
      enddo

     Hs = MATMUL(eye,Hs)
     Hs = MATMUL(Hs,eye)




      do islab=1, Nslab
         do jslab=1, Nslab
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Hs(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
                  Hs(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
      enddo

   !   tx = MATMUL(Hs,tx)
   !   ty = MATMUL(Hs,ty)
   !   tz = MATMUL(Hs,tz)
   !   print *,'tx'
   !   print *, tx
   !   print *,'ty'
   !   print *, ty
   !   print *,'tz'
   !   print *, tz


end subroutine ASR_R_G


subroutine ASR_R_G_ribbon(Hs)
     ! This subroutine projects out the translational part
     ! of the slab Hamiltonian in real space 
     ! prior to Fourier interpolation onto reciprocal space
     ! 22/06/2026 F Ballester
     
     
     use para

     implicit none

     complex(dp), intent(inout) :: Hs(nslab1*nslab2*Num_wann,nslab1*nslab2*Num_wann)

     integer :: a,b, islab1, jslab1, islab2, jslab2, kappa, na, kappap
     
     real(dp) :: G(nslab1*nslab2*Num_wann,nslab1*nslab2*Num_wann), tx(nslab1*nslab2*Num_wann), ty(nslab1*nslab2*Num_wann), tz(nslab1*nslab2*Num_wann), eye(nslab1*nslab2*Num_wann,nslab1*nslab2*Num_wann)

     G = 0.0d0
     eye = 0.0d0

     tx=0.0d0
     tx(1::3) = 1.0d0
     ty=0.0d0
     ty(2::3) = 1.0d0
     tz=0.0d0
     tz(3::3) = 1.0d0

     do a=1, nslab1*nslab2*Num_wann
      eye(a,a) = 1.0d0
      do b=1, nslab1*nslab2*Num_wann
         G(a,b) = 3.0d0/dble(nslab1*nslab2*Num_wann)*&
         (tx(a)*tx(b)+ty(a)*ty(b)+tz(a)*tz(b))
      enddo
     enddo
   !    print *, 'start G'
   !   do a =1, nslab*num_wann
   !    print *, G(a,:) 
   !   enddo
   !   print *, 'end G'

     eye(:,:) = eye(:,:) - G(:,:)



      na = 3*Origin_cell%Num_atoms

      do islab1=1, Nslab1
         do jslab1=1, Nslab2

         do islab2=1, Nslab1
         do jslab2=1, Nslab2
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Hs(na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+1:na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+3,&
                  na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+1:na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+3) =&
                  Hs(na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+1:na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+3,&
                  na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+1:na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+3) *&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
         enddo
         enddo
      enddo

     Hs = MATMUL(eye,Hs)
     Hs = MATMUL(Hs,eye)





      do islab1=1, Nslab1
         do jslab1=1, Nslab2

         do islab2=1, Nslab1
         do jslab2=1, Nslab2
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Hs(na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+1:na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+3,&
                  na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+1:na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+3) =&
                  Hs(na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+1:na*(islab1-1)*nslab2+na*(jslab1)+3*(kappa-1)+3,&
                  na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+1:na*(islab2-1)*nslab2+na*(jslab2)+3*(kappa-1)+3) /&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
         enddo
         enddo
      enddo

   !   tx = MATMUL(Hs,tx)
   !   ty = MATMUL(Hs,ty)
   !   tz = MATMUL(Hs,tz)
   !   print *,'tx'
   !   print *, tx
   !   print *,'ty'
   !   print *, ty
   !   print *,'tz'
   !   print *, tz


end subroutine ASR_R_G_ribbon



subroutine apply_ASR_slab(Ham_to_ASR)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on Eq. (81) of X. Gonze et al, PRB 50. 13035 (1994) 
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, m1, m2,islab,jslab, lamd, i, j, na, nb, nasr

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(inout) ::Ham_to_ASR(Num_wann*nslab,Num_wann*nslab) 

     ! Hamiltonian of slab at q=0
     complex(Dp), allocatable :: Ham0(:,:), Ham_no(:,:)

     ! Sum of Hamiltonian of slab at q=0
     complex(Dp), allocatable :: sum(:,:)

     ! Displacement vector for debugging
     real(Dp), allocatable :: onesvec(:), resultsofvec(:)

     real(Dp) :: k0(2)

     double precision :: sumd, Qdd

     allocate(Ham0(Num_wann*nslab,Num_wann*nslab), Ham_no(Num_wann*nslab,Num_wann*nslab))
     allocate(onesvec(Num_wann*nslab), resultsofvec(Num_wann*nslab))
     

     !!!!!!!!!!!!!!!!!!!!!!!
     ! TODO: ADD MASSES
     !!!!!!!!!!!!!!!!!!!!!!
     
     k0 = 0.0d0
     Ham0 = 0.0d0
     Ham_no = 0.0d0
     ! Get H(q=0)
     call ham_slab(k0,Ham0)
     Ham_no = Ham0

   !   call apply_ASR_slab_iterative_Gamma(Ham0)

   !    Ham0 = Ham_no-Ham0 
      ! resultsofvec=0.0d0
      ! call eigensystem_c('V', 'L', Num_wann*Nslab, Ham0, resultsofvec)  
      ! call force_positive_definite_slab(Ham0, resultsofvec)

     onesvec = 0.0d0
     onesvec(3::3) = 1.0d0
   ! !   print *, onesvec
     resultsofvec = matmul(Ham_to_ASR, onesvec)

   !   print *, 'vector pre-ASR'
   !   print *, resultsofvec 
   !   print *, 'end of vector pre-ASR'

   na = 3*Origin_cell%Num_atoms

   do islab=1, Nslab
      do jslab=1, Nslab
         do kappa=1, Origin_cell%Num_atoms
            do kappap=1, Origin_cell%Num_atoms
               Ham0(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
               Ham0(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
               sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))

               Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
               Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
               sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
            enddo
         enddo
      enddo
   enddo


     do ialpha=1, 3
      do ibeta=1, 3
         do kappa=1, Origin_cell%Num_atoms*Nslab
               sumd = 0.0d0 
               Qdd = 0.0d0
               do kappapp=1, Origin_cell%Num_atoms*Nslab
                 sumd =&
                 sumd + Ham0(3*(kappa-1)+ialpha,3*(kappapp-1)+ibeta)! + Ham0(3*(kappa-1)+ibeta,3*(kappapp-1)+ialpha))/2.0d0
                 Qdd = Qdd !+ (Ham0(3*(kappa-1)+ialpha,3*(kappapp-1)+ibeta) - Ham0(3*(kappa-1)+ibeta,3*(kappapp-1)+ialpha))
               enddo
               Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) =Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) +  DCMPLX(-(sumd+Qdd), 0.d0) 
               ! Ham_no(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) =Ham_no(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) +  DCMPLX(-(sumd+Qdd), 0.d0) 
         enddo
      enddo
     end do


   do islab=1, Nslab
      do jslab=1, Nslab
         do kappa=1, Origin_cell%Num_atoms
            do kappap=1, Origin_cell%Num_atoms
               Ham0(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
               Ham0(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
               sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))

               Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
               Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
               sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
            enddo
         enddo
      enddo
   enddo

   !   Ham0 = Ham_no

   !   do ialpha=1, 3
   !    do ibeta=1, 3
   !       do kappa=1, Origin_cell%Num_atoms*Nslab
   !             sumd = 0.0d0 
   !             Qdd = 0.0d0
   !             do kappapp=1, Origin_cell%Num_atoms*Nslab
   !               sumd =&
   !               sumd + Ham0(3*(kappa-1)+ialpha,3*(kappapp-1)+ibeta)! + Ham0(3*(kappa-1)+ibeta,3*(kappapp-1)+ialpha))/2.0d0

   !               Qdd = Qdd !+ (Ham0(3*(kappa-1)+ialpha,3*(kappapp-1)+ibeta) - Ham0(3*(kappa-1)+ibeta,3*(kappapp-1)+ialpha))
   !             enddo
   !             Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) =Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) +  DCMPLX(-(sumd+Qdd), 0.d0) 
   !             Ham_no(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) =Ham_no(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) +  DCMPLX(-(sumd+Qdd), 0.d0) 
   !       enddo
   !    enddo
   !   end do

      ! do i=1, Num_wann*nslab
      !    do j=1, ialpha - 1
      !       Ham_to_ASR(i,j) = 0.5d0* (Ham_to_ASR(i,j)+CONJG(Ham_to_ASR(j,i)))
      !       Ham_to_ASR(j,i) = CONJG(Ham_to_ASR(i,j))
      !    end do
      ! end do

      ! do ialpha=1, 3
      ! do ibeta=ialpha, 3
      !    do kappa=1, Origin_cell%Num_atoms*Nslab
      !       do kappap=kappa, Origin_cell%Num_atoms*Nslab
      !          Ham_to_ASR(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = &
      !             conjg(Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta))
      !       enddo
      !    enddo
      ! enddo
      ! enddo
   
   !   print *, 'chamk in asr'
   !      do j=1, Nslab* Num_wann
   !        print *, Ham_to_ASR(j,:)
   !      end do
       
   !      print *, 'end of chamk in asr'
   !   onesvec = 0.0d0
   !   onesvec(3::3) = 1.0d0
   ! ! !   print *, onesvec
   !   resultsofvec = matmul(Ham_to_ASR, onesvec)

   !   print *, 'vector'
   !   print *, resultsofvec 
   !   print *, 'end of vector'

   !   print *, 'hamtoasr'
   !   print *, dble(Ham_to_ASR)
   !   print *, 'end of hamtoasr'
     
   
     deallocate(Ham0, Ham_no)
     deallocate(onesvec,resultsofvec)
   return
end subroutine apply_ASR_slab

subroutine apply_ASR_slab_iterative_Gamma(Ham_to_ASR)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on https://doi.org/10.1016/j.cpc.2011.04.019
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, niters, iiter, na, islab, jslab

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(inout) ::Ham_to_ASR(Num_wann*nslab,Num_wann*nslab) 

     complex(Dp) :: sum0
     

    na = 3*Origin_cell%Num_atoms

      do islab=1, Nslab
         do jslab=1, Nslab
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
                  Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) *&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
      enddo
     niters=100
     do iiter=1,niters

      !> naive ASR
      do ialpha=1,3
         do ibeta=1,3
            do kappa=1,Origin_cell%Num_atoms*nslab
               sum0=0.0d0
               do kappap =1, Origin_cell%Num_atoms*nslab
                  sum0 = sum0 + Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
               sum0 = sum0/(Origin_cell%Num_atoms*nslab)
               do kappap =1, Origin_cell%Num_atoms*nslab
                  Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
               end do
            enddo
         enddo
      enddo

      !> Symmetrize
      do kappa=1, Origin_cell%Num_atoms*nslab
         do kappap=kappa, Origin_cell%Num_atoms*nslab
            sum0=0.0d0
            do ialpha=1,3
               do ibeta=1,3
               sum0 = (Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) + conjg(Ham_to_ASR(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha)))*0.5d0
               Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = sum0
               Ham_to_ASR(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = sum0
               end do
            end do
         end do
      end do

      
     enddo
     !> Symmetric ASR
      do ialpha=1,3
         do ibeta=1,3
            do kappa=1,Origin_cell%Num_atoms*nslab
               sum0=0.0d0
               do kappap =1, Origin_cell%Num_atoms*nslab
                  sum0 = sum0 + Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
               sum0 = sum0/(Origin_cell%Num_atoms*nslab-kappa+1)
               do kappap =1, Origin_cell%Num_atoms*nslab
                  Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
                  Ham_to_ASR(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
            enddo
         enddo
      enddo

      do islab=1, Nslab
         do jslab=1, Nslab
            do kappa=1, Origin_cell%Num_atoms
               do kappap=1, Origin_cell%Num_atoms
                  Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) =&
                  Ham_to_ASR(na*(islab-1)+3*(kappa-1)+1:na*(islab-1)+3*(kappa-1)+3,na*(jslab-1)+3*(kappap-1)+1:na*(jslab-1)+3*(kappap-1)+3) /&
                  sqrt(Atom_Mass(kappa)*Atom_Mass(kappap))
               enddo
            enddo
         enddo
      enddo


   return
end subroutine apply_ASR_slab_iterative_Gamma


subroutine apply_ASR_RealSpace(irorigin)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on Eq. (82) of X. Gonze et al, PRB 50. 13035 (1994) 
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, irab

     integer, intent(in) :: irorigin ! index of R(ir)=(0, 0, 0)

     ! Hamiltonian of system to apply the ASR to
   !   complex(Dp),intent(out) ::Ham_to_ASR(Num_wann*nslab,Num_wann*nslab, Nrpts) 


     ! Sum of Hamiltonian
     complex(Dp), allocatable :: sum(:,:)

     allocate(sum(Num_wann,3))

     sum = 0.0d0
     do kappapp=1, Origin_cell%Num_atoms
      do kappa=1, Origin_cell%Num_atoms
         do ialpha=1, 3
            do ibeta=1, 3
               do irab=1, Nrpts
                  if (.not.(irab.eq.irorigin).or.(.not.(kappa.eq.kappapp))) then
                     sum(3*(kappa-1)+ialpha,ibeta) = sum(3*(kappa-1)+ialpha,ibeta) + HmnR(3*(kappa-1)+ialpha,3*(kappapp-1)+ibeta, irab) 
                  endif
               end do
            end do
         enddo
      enddo
     enddo

   
     do kappa=1, Origin_cell%Num_atoms
      do ialpha=1, 3
         do ibeta=1, 3
            HmnR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta,irorigin) = -1.0d0*sum(3*(kappa-1)+ialpha,ibeta)
         end do
      enddo
     enddo


     

     deallocate(sum)
   return
end subroutine apply_ASR_RealSpace



subroutine apply_ASR_slab_RealSpace(Ham_to_ASR)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on Eq. (82) of X. Gonze et al, PRB 50. 13035 (1994) 
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, irab


   !   Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(out) ::Ham_to_ASR(-ijmax:ijmax,Num_wann,Num_wann) 


     ! Sum of Hamiltonian
     complex(Dp) :: sum

   !   allocate(sum(Num_wann))

     do kappa=1, Origin_cell%Num_atoms
      do ialpha=1, 3
         do ibeta=1,3
            sum = 0.0d0
            do irab=-ijmax, ijmax
               do kappapp=1, Origin_cell%Num_atoms
               if ((.not.(kappa.eq.kappapp)).and.(.not.(irab.eq.0))) then
                  sum = sum + Ham_to_ASR(irab,3+kappa+ialpha,3*kappapp+ibeta) 
               endif
               enddo
            enddo
            Ham_to_ASR(0,3+kappa+ialpha,3*kappa+ibeta) = -sum
         enddo
      enddo
     enddo



   return
end subroutine apply_ASR_slab_RealSpace


subroutine apply_ASR_slab_iterative(Ham_to_ASR)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on https://doi.org/10.1016/j.cpc.2011.04.019
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, niters, iiter

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(out) ::Ham_to_ASR(Num_wann*nslab,Num_wann*nslab) 

     ! Hamiltonian of slab at q=0
     complex(Dp), allocatable :: Ham0(:,:)

     ! Sum of Hamiltonian of slab at q=0
     complex(Dp), allocatable :: sum(:,:)

     complex(Dp) :: sum0

     allocate(Ham0(Num_wann*nslab,Num_wann*nslab))
     allocate(sum(Num_wann*nslab,3))
   
     Ham0 = 0.0d0
     ! Get H(q=0)
     call ham_slab((/0.0d0, 0.0d0/),Ham0)


     niters=100
     do iiter=1,niters

      !> naive ASR
      do ialpha=1,3
         do ibeta=1,3
            do kappa=1,Origin_cell%Num_atoms*nslab
               sum0=0.0d0
               do kappap =1, Origin_cell%Num_atoms*nslab
                  sum0 = sum0 + Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
               sum0 = sum0/(Origin_cell%Num_atoms*nslab)
               do kappap =1, Origin_cell%Num_atoms*nslab
                  Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
               end do
            enddo
         enddo
      enddo

      !> Symmetrize
      do kappa=1, Origin_cell%Num_atoms*nslab
         do kappap=kappa, Origin_cell%Num_atoms*nslab
            sum0=0.0d0
            do ialpha=1,3
               do ibeta=1,3
               sum0 = (Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) + conjg(Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha)))*0.5d0
               Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = sum0
               Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = sum0
               end do
            end do
         end do
      end do

      
     enddo
     !> Symmetric ASR
      do ialpha=1,3
         do ibeta=1,3
            do kappa=1,Origin_cell%Num_atoms*nslab
               sum0=0.0d0
               do kappap =1, Origin_cell%Num_atoms*nslab
                  sum0 = sum0 + Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
               sum0 = sum0/(Origin_cell%Num_atoms*nslab-kappa+1)
               do kappap =1, Origin_cell%Num_atoms*nslab
                  Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
                  Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
            enddo
         enddo
      enddo




   
     !> ASR on the rest?????
     sum = 0.0d0
     do kappapp=1, Origin_cell%Num_atoms*nslab
      do ibeta=1, 3
         sum(:,ibeta) = sum(:,ibeta) + (Ham0(:,3*(kappapp-1)+ibeta))
      end do
     enddo
     !> skip kappap by imposing delta_kappa,kappap
     do kappa=1, Origin_cell%Num_atoms*nslab
         do ialpha=1, 3
            do ibeta=1, 3
               Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) = Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) - sum(3*(kappa-1)+ialpha,ibeta)
            end do
         enddo
     enddo

     deallocate(Ham0)
     deallocate(sum)
   return
end subroutine apply_ASR_slab_iterative


subroutine apply_ASR_slab_iterative_RealSpace(Ham_to_ASR)
     ! This subroutine is used to apply the ASR for
     ! slab systems. 
     !
     ! We use the scheme based on https://doi.org/10.1016/j.cpc.2011.04.019
     ! 
     ! History  
     !       24/04/2026  Francesc Ballester 

  
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, niters, iiter

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(out) ::Ham_to_ASR(-ijmax:ijmax,Num_wann,Num_wann) 

     ! Hamiltonian of slab at q=0
     complex(Dp), allocatable :: Ham0(:,:)

     ! Sum of Hamiltonian of slab at q=0
     complex(Dp), allocatable :: sum(:,:)

     complex(Dp) :: sum0

     allocate(Ham0(Num_wann,Num_wann))
     allocate(sum(Num_wann,3))
   
     ! Get H(q=0)
   !   call ham_slab((/0.0d0, 0.0d0/),Ham0)
     Ham0= 0.0d0
     Ham0(:,:) = Ham_to_ASR(0,:,:)

     niters=1
     do iiter=1,niters

      !> naive ASR
      do ialpha=1,3
         do ibeta=1,3
            do kappa=1,Origin_cell%Num_atoms
               sum0=0.0d0
               do kappap =1, Origin_cell%Num_atoms
                  sum0 = sum0 + Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
               end do
               sum0 = sum0/(Origin_cell%Num_atoms)
               do kappap =1, Origin_cell%Num_atoms
                  Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
               end do
            enddo
         enddo
      enddo

      !> Symmetrize
      do kappa=1, Origin_cell%Num_atoms
         do kappap=kappa, Origin_cell%Num_atoms
            sum0=0.0d0
            do ialpha=1,3
               do ibeta=1,3
               sum0 = (Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) + Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha))*0.5d0
               Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = sum0
               Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = sum0
               end do
            end do
         end do
      end do

      
     enddo
     !> Symmetric ASR
      ! do ialpha=1,3
      !    do ibeta=1,3
      !       do kappa=1,Origin_cell%Num_atoms
      !          sum0=0.0d0
      !          do kappap =1, Origin_cell%Num_atoms
      !             sum0 = sum0 + Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
      !          end do
      !          sum0 = sum0/(Origin_cell%Num_atoms-kappa+1)
      !          do kappap =1, Origin_cell%Num_atoms
      !             Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta) - sum0
      !             Ham0(3*(kappap-1)+ibeta,3*(kappa-1)+ialpha) = Ham0(3*(kappa-1)+ialpha,3*(kappap-1)+ibeta)
      !          end do
      !       enddo
      !    enddo
      ! enddo

      Ham_to_ASR(0,:,:) = Ham0(:,:)



   
   !   !> ASR on the rest?????
   !   sum = 0.0d0
   !   do kappapp=1, Origin_cell%Num_atoms
   !    do ibeta=1, 3
   !       sum(:,ibeta) = sum(:,ibeta) + (Ham0(:,3*(kappapp-1)+ibeta))
   !    end do
   !   enddo
   !   !> skip kappap by imposing delta_kappa,kappap
   !   do kappa=1, Origin_cell%Num_atoms
   !       do ialpha=1, 3
   !          do ibeta=1, 3
   !             Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) = Ham_to_ASR(3*(kappa-1)+ialpha,3*(kappa-1)+ibeta) - sum(3*(kappa-1)+ialpha,ibeta)
   !          end do
   !       enddo
   !   enddo
     deallocate(Ham0)
     deallocate(sum)
   return
end subroutine apply_ASR_slab_iterative_RealSpace


subroutine force_positive_definite_slab(Ham_to_force, omega2)
     use para
     implicit none

     ! loop index  
     integer :: ialpha, ibeta, kappa, kappap, kappapp, m1, m2,islab,jslab, lamd, i, j, na, nb, nasr, a, b, mu

     ! Hamiltonian of slab system to apply the ASR to
     complex(Dp),intent(inout) ::Ham_to_force(Num_wann*nslab,Num_wann*nslab)
     real(Dp), intent(in) :: omega2(Num_wann*nslab)

    
     complex(Dp), allocatable :: Hamdummy(:,:)

     allocate(Hamdummy(Num_wann*nslab,Num_wann*nslab))
     
     Hamdummy = 0.0d0

     do i=1, nslab
      do j=1, nslab
         do kappa=1, Origin_cell%Num_atoms
            do kappap=1, Origin_cell%Num_atoms
               do ialpha=1, 3
                  do ibeta=1, 3
                     a = 3*Origin_cell%Num_atoms*(i-1)+3*(kappa-1)+ialpha
                     b = 3*Origin_cell%Num_atoms*(j-1)+3*(kappap-1)+ibeta
                     do mu=1,Num_wann*nslab
                        Hamdummy(a,b) = Hamdummy(a,b) +&
                        abs(omega2(mu))*Ham_to_force(a,mu)*conjg(Ham_to_force(b,mu))
                     enddo
                  end do
               end do
            end do
         end do
      end do
     end do

     do a=1, Num_wann*nslab
      do b=1, Num_wann*nslab
         Ham_to_force(a,b) = Hamdummy(a,b)
      end do
     end do


   !   call eigensystem_c('V', 'L', Num_wann*Nslab, Ham_to_force, omega2) 

   deallocate(Hamdummy)

end subroutine force_positive_definite_slab
! subroutine remove_Trans(Ham_to_ASR)
!      ! This subroutine is used to apply the ASR for
!      ! slab systems by projecting out directly the translational vectors. 
!      !
!      ! 
!      ! History  
!      !       24/04/2026  Francesc Ballester 

  
!      use para
!      implicit none

!      ! loop index  
!      integer :: ialpha, ibeta, kappa, kappap, kappapp  

!      ! Hamiltonian of slab system to apply the ASR to
!      complex(Dp),intent(out) ::Ham_to_ASR(Num_wann,Num_wann) 

!      ! translational projection matrix
!      complex(Dp), allocatable :: transproj(:,:)

!      ! translational vector
!      complex(Dp), allocatable :: v1(:)

!      allocate(transproj(Num_wann,Num_wann))
!      allocate(v1(Num_wann))
   
!      transproj = 0.0d0

!      do ialpha=1, Num_wann
!       transproj(ialpha,ialpha) = 1.0d0
!      end do

!      do ialpha=1, 3
!       v1 = 0.0d0
!       do kappa=1, Origin_cell%Num_atoms
!          v1(3*(kappa-1)+ialpha) = 1.0d0/sqrt(dble(Origin_cell%Num_atoms))
!       enddo
!       do kappap = 1,Num_wann
!          do kappapp=1, Num_wann
!             transproj(kappap,kappapp) = transproj(kappap,kappapp) - (v1(kappap) * v1(kappapp))
!          enddo
!       enddo 
!      enddo

!      Ham_to_ASR = MATMUL(Ham_to_ASR, transproj)

!      Ham_to_ASR = MATMUL(transproj, Ham_to_ASR)

!      deallocate(transproj)
!      deallocate(v1)

! end subroutine remove_Trans
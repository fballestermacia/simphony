! This subroutine is used to caculate Hamiltonian for
! bulk system .

! History
!        May/29/2011 by Quansheng Wu

subroutine ham_bulk_atomicgauge(k,Hamk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system with the consideration of the atom's position
   !
   ! History
   !
   !        May/29/2011 by Quansheng Wu
   !  Atomic gauge Guan Yifei 2019
   !  Lattice gauge Hl
   !  Atomic gauge Ha= U* Hl U 
   !  where U = e^ik.wc(i) on diagonal

   use para
   implicit none

   integer :: i1,i2,iR

   ! wave vector in 3d
   real(Dp) :: k(3), kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)
   complex(dp), allocatable :: mat1(:, :)
   real(dp), external :: norm

   allocate(mat1(Num_wann, Num_wann))

   Hamk_bulk=0d0
  !do iR=1, Nrpts
  !   kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
  !   factor= exp(pi2zi*kdotr)

  !   Hamk_bulk(:, :)= Hamk_bulk(:, :) &
  !      + HmnR(:, :, iR)*factor/ndegen(iR)
  !enddo ! iR
 
  !mat1=0d0
  !do i1=1,Num_wann
  !   pos0=Origin_cell%wannier_centers_direct(:, i1)
  !   kdotr= k(1)*pos0(1)+ k(2)*pos0(2)+ k(3)*pos0(3)
  !   mat1(i1,i1)= exp(pi2zi*kdotr)
  !enddo
  !Hamk_bulk=matmul(conjg(mat1),matmul(Hamk_bulk,mat1))


   !> the first atom in home unit cell
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

            Hamk_bulk(i1, i2)= Hamk_bulk(i1, i2) &
               + HmnR(i1, i2, iR)*factor/ndegen(iR)
         enddo ! i1
      enddo ! i2
   enddo ! iR

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
            write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
            !stop
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_atomicgauge


subroutine valley_k_atomicgauge(k,valley_k)
   ! This subroutine performs the Fourier transform of avalley operator
   ! History
   !        Nov/5/2023 by Quansheng Wu

   use para
   implicit none

   integer :: i1,i2,iR

   ! wave vector in 3d
   real(Dp) :: k(3), kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   ! Hamiltonian of bulk system
   complex(Dp),intent(out) :: valley_k(Num_wann, Num_wann)
   real(dp), external :: norm
   

   valley_k= 0d0
   !> the first atom in home unit cell
   do iR=1, Nrpts_valley
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec_valley(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

            valley_k(i1, i2)= valley_k(i1, i2) &
               + valley_operator_R(i1, i2, iR)*factor
         enddo ! i1
      enddo ! i2
   enddo ! iR

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(valley_k(i1,i2)-conjg(valley_k(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', valley_k(i1, i2)
            write(stdout,*)'value at (i2, i1)', valley_k(i2, i1)
            !stop
         endif
      enddo
   enddo


   return
end subroutine valley_k_atomicgauge

subroutine d2Hdk2_atomicgauge(k, DHDk2_wann)
   !> second derivatve of H(k)
   use para, only : Nrpts, irvec, crvec, Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi2zi, zi, twopi, pi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> second derivate of H(k)
   complex(dp), intent(out) :: DHDk2_wann(Num_wann, Num_wann, 3, 3)

   integer :: iR, i1, i2, i, j

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   DHDk2_wann= 0d0
   !> the first atom in home unit cell
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)
            pos_direct= pos_direct+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)

            do j=1, 3
               do i=1, 3
                  DHDk2_wann(i1, i2, i, j)=DHDk2_wann(i1, i2, i, j) &
                     -pos_cart(i)*pos_cart(j)*HmnR(i1, i2, iR)*ratio
               enddo ! j 
            enddo ! i
         enddo ! i1
      enddo ! i2
   enddo ! iR

   return
end subroutine d2Hdk2_atomicgauge

subroutine dHdk_atomicgauge(k, velocity_Wannier)
   !> Velocity operator in Wannier basis using atomic gauge
   !> First derivate of H(k); dH(k)/dk
   use para, only : Nrpts, irvec, Origin_cell, HmnR, ndegen, &
       Num_wann, dp, Rcut, pi2zi,  &
       zi, soc, zzero, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator in Wannier basis using atomic gauge 
   complex(dp), intent(out) :: velocity_Wannier(Num_wann, Num_wann, 3)

   integer :: iR, i1, i2, i

   real(dp) :: pos1(3), pos2(3), pos_cart(3), pos_direct(3)
   real(dp) :: kdotr, dis
   complex(dp) :: ratio
   real(dp), external :: norm

   velocity_Wannier= zzero
   do iR=1, Nrpts
      do i2=1, Num_wann
         pos2= Origin_cell%wannier_centers_direct(:, i2)
         !> the second atom in unit cell R
         do i1=1, Num_wann
            !> the first atom in home unit cell
            pos1= Origin_cell%wannier_centers_direct(:, i1)
            pos_direct= irvec(:, iR)+ pos2- pos1

            call direct_cart_real(pos_direct, pos_cart, Origin_cell%lattice)

            dis= norm(pos_cart)
            if (dis> Rcut) cycle

            kdotr=k(1)*pos_direct(1) + k(2)*pos_direct(2) + k(3)*pos_direct(3)
            ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))/ndegen(iR)

            do i=1, 3
               velocity_Wannier(i1, i2, i)=velocity_Wannier(i1, i2, i)+ &
                  zi*pos_cart(i)*HmnR(i1, i2, iR)*ratio
            enddo ! i

         enddo ! i2
      enddo ! i1
   enddo ! iR

   return
end subroutine dHdk_atomicgauge

subroutine dHdk_atomicgauge_Ham(k, eigvec, Vmn_Ham)
   !> Velocity operator in Hamiltonian basis using atomic gauge
   !> see https://www.wanniertools.org/theory/tight-binding-model/
   use para, only : Num_wann, dp
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> eigenvectors of H, H*eigvec(:, n)= E(n)*eigvec(:, n)
   complex(dp), intent(in) :: eigvec(Num_wann, Num_wann)

   !> velocity operator in the diagonalized Hamiltonian basis using lattice gauge
   complex(dp), intent(out) :: Vmn_Ham(Num_wann, Num_wann, 3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), allocatable :: Vmn_wann(:, :, :)

   integer :: i

   allocate(Vmn_wann(Num_wann, Num_wann, 3))
   Vmn_Ham= 0d0
   call dHdk_atomicgauge(k, Vmn_wann)
   do i=1, 3
      call rotation_to_Ham_basis(eigvec, Vmn_wann(:, :, i), Vmn_Ham(:, :, i))
   enddo
   deallocate(Vmn_wann)

   return
end subroutine dHdk_atomicgauge_Ham

subroutine ham_bulk_latticegauge(k,Hamk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system without the consideration of the atom's position
   !
   ! History
   !
   !        May/29/2011 by Quansheng Wu

   use para, only : dp, pi2zi, HmnR, ndegen, nrpts, irvec, Num_wann, stdout, twopi, zi
   implicit none

   ! loop index
   integer :: i1,i2,iR

   real(dp) :: kdotr

   complex(dp) :: factor

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   Hamk_bulk=0d0
   do iR=1, Nrpts
      kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
      factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

      Hamk_bulk(:, :)= Hamk_bulk(:, :)+ HmnR(:, :, iR)*factor/ndegen(iR)
   enddo ! iR
   
   !call mat_mul(Num_wann, mirror_z, Hamk_bulk, mat1)
   !call mat_mul(Num_wann, mat1, mirror_z, mat2)
   !Hamk_bulk= (Hamk_bulk+ mat2)/2d0

   ! check hermitcity
  !do i1=1, Num_wann
  !   do i2=1, Num_wann
  !      if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
  !         write(stdout,*)'there is something wrong with Hamk_bulk'
  !         write(stdout,*)'i1, i2', i1, i2
  !         write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
  !         write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
  !         !stop
  !      endif
  !   enddo
  !enddo

   return
end subroutine ham_bulk_latticegauge

subroutine dHdk_latticegauge_wann(k, velocity_Wannier)
   use para, only : Nrpts, irvec, crvec, Origin_cell, &
      HmnR, ndegen, Num_wann, zi, pi2zi, dp, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), intent(out) :: velocity_Wannier(Num_wann, Num_wann, 3)

   integer :: iR, i

   real(dp) :: kdotr
   complex(dp) :: ratio

   do i=1, 3
      do iR= 1, Nrpts
         kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
         ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
         velocity_Wannier(:, :, i)= velocity_Wannier(:, :, i)+ &
            zi*crvec(i, iR)*HmnR(:,:,iR)*ratio/ndegen(iR)
      enddo ! iR
   enddo

   return
end subroutine dHdk_latticegauge_wann

subroutine dHdk_latticegauge_Ham(k, eigval, eigvec, Vmn_Ham)
   use para, only : Nrpts, irvec, crvec, Origin_cell, &
      HmnR, ndegen, Num_wann, zi, pi2zi, dp, zzero, twopi
   implicit none

   !> momentum in 3D BZ
   real(dp), intent(in) :: k(3)

   real(dp), intent(in) :: eigval(Num_wann)
   complex(dp), intent(in) :: eigvec(Num_wann, Num_wann)

   !> velocity operator in the diagonalized Hamiltonian basis using lattice gauge
   complex(dp), intent(out) :: Vmn_Ham(Num_wann, Num_wann, 3)

   !> velocity operator in Wannier basis using lattice gauge
   complex(dp), allocatable :: Vmn_wann(:, :, :), Vmn_Ham0(:, :, :)

   !> wnm=En-Em
   complex(dp), allocatable :: wnm(:, :), temp(:, :)

   !> it's a diagonal matrix for each axis
   complex(dp), allocatable :: itau(:, :, :)

   integer :: i, m, n, l

   allocate(Vmn_wann(Num_wann, Num_wann, 3))
   allocate(Vmn_Ham0(Num_wann, Num_wann, 3))
   Vmn_wann= 0d0; Vmn_Ham0= 0d0

   call dHdk_latticegauge_wann(k, Vmn_wann)
   do i=1, 3
      call rotation_to_Ham_basis(eigvec, Vmn_wann(:, :, i), Vmn_Ham0(:, :, i))
   enddo
 
   allocate(wnm(Num_wann, Num_wann))
   allocate(temp(Num_wann, Num_wann))
   allocate(itau(Num_wann, Num_wann, 3))
   wnm=0d0; temp=0d0; itau= zzero

   !\  -i\tau
   do i=1, 3
      do m=1, Num_wann
         itau(m, m, i)= -zi*Origin_cell%wannier_centers_cart(i, m)
      enddo
   enddo

   do m=1, Num_wann
      do n=1, Num_wann
         wnm(m, n)= eigval(n)-eigval(m)
      enddo
   enddo
  
   !> \sum_l (-i*\tau_{l})*conjg(vec(l, m))*vec(l, n)
   do i=1, 3
      call mat_mul(Num_wann, itau(:, :, i), eigvec, temp)
      call mat_mul(Num_wann, conjg(transpose(eigvec)), temp, Vmn_Ham(:, :, i))
   enddo
 
   do i=1, 3
      do n=1, Num_wann
         do m=1, Num_wann
            temp(m, n) = wnm(m, n)*Vmn_Ham(m, n, i)
         enddo
      enddo
      Vmn_Ham(:, :, i)= temp(:, :)
   enddo


   Vmn_Ham= Vmn_Ham+ Vmn_Ham0
  !Vmn_Ham = Vmn_Ham0

!  print *, Vmn_Ham(1, 1, 1)

  !Vmn_Ham= 0d0
  !do m=1, Num_wann
  !   do n=1, Num_wann
  !      do i=1, 3
  !         do l=1, Num_wann
  !         Vmn_Ham(m, n, i)= Vmn_Ham(m, n, i)- &
  !            (eigval(n)- eigval(m))*zi*Origin_cell%wannier_centers_cart(i, l)*conjg(eigvec(l, m))*eigvec(l, n)
  !         enddo ! summation over l
  !      enddo
  !   enddo
  !enddo
  !Vmn_Ham= Vmn_Ham+ Vmn_Ham0

!  print *, Vmn_Ham(1, 1, 1)
!  stop

   return
end subroutine dHdk_latticegauge_Ham

subroutine long_range_phonon_interaction(nfr1,nfr2,nfr3,k,loto_2d,sign,Hamk_bulk,tau,zeu) 
   ! This subroutine computes the rigid-ion (long-range) term for q
   ! X. Gonze et al, PRB 50. 13035 (1994) 
   ! the Ewald parameter alpha must be large enough to
   ! have negligible r-space contribution
   !
   ! Note that this subroutine is based on the rgd_blk subroutine found on
   ! SSCHA's CellConstructor, which is based on the subroutine with the same name found 
   ! on Quantum ESPRESSO's source code
   !
   ! History
   !        
   !        June/13/2024 by Francesc Ballester
   use para

   implicit none

   ! Diele_Tensor(3,3)   = epsil              Dielectric tensor
   ! Born_Charge(Origin_cell$Num_atoms,3,3) = zeu     Born effective charges  in CC (3,3,Origin_cell%Num_Atoms)
   ! Origin_cell%lattice          = at    Lattice vectors
   ! Origin_cell%reciprocal_lattice  = bg   Reciprocal lattice vectors
   ! Origin_cell%CellVolume  = Omega
   ! alat = Origin_cell%cell_parameters(1) ! In Bohrs
   ! Origin_cell%Atom_position_direct(3,Origin_cell%Num_atoms) = tau  Atom positions in direct coordinates
   integer :: nfr1, nfr2, nfr3   !  FFT grid             
   complex(Dp) :: Hamk_bulk(Num_wann, Num_wann)
   real(Dp) &
        k(3),           &! q-vector
        sign,          &! sign=+/-1.0 ==> add/subtract rigid-ion term
        tau(3,Origin_cell%Num_atoms),  & ! atom position
        zeu(Origin_cell%Num_atoms,3,3) ! Effective charges
   logical :: loto_2d ! 2D LOTO correction 

   
   ! local variables

   real(Dp):: geg, gp2, r                    !  <q+G| epsil | q+G>,  For 2d loto: gp2, r
   integer :: na,nb, i,j, m1, m2, m3
   integer :: nr1x, nr2x, nr3x
   real(Dp) :: alph, fac,g1,g2,g3, facgd, arg, gmax, e2
   !real(DP) :: tau(3,Origin_cell%Num_atoms) 
   real(Dp) :: zag(3),zbg(3),zcg(3), fnat(3), reff(2,2)
   complex(dp) :: facg, fmtx(3,3)
   real(Dp) :: q(3) !new k
   real(Dp) :: rec_lattice(3,3) 
   !
   ! alph is the Ewald parameter, geg is an estimate of G^2
   ! such that the G-space sum is convergent for that alph
   ! very rough estimate: geg/4/alph > gmax = 14
   ! (exp (-14) = 10^-6)
   !
   complex(Dp) :: dummydebug(3,Origin_cell%Num_atoms,3,Origin_cell%Num_atoms),dummydebug2(Num_wann,Num_wann), dummyFC(Num_wann,Num_wann,1)

   integer :: CartToOrb(3) !map cartesian axis to orbital ordering


   CartToOrb=(/1,2,3/)
   dummydebug = 0d0
   dummydebug2 = 0d0

   e2 = 2.0d0

   
   gmax= 14.0d0
   alph= 1.0d0
   geg = gmax*alph*4.0d0

   rec_lattice = Origin_cell%reciprocal_lattice*Origin_cell%cell_parameters(1)/(twopi) !This value is the same as QE
   
   q(1) = k(1)*rec_lattice(1,1) + k(2)*rec_lattice(1,2) + k(3)*rec_lattice(1,3)
   q(2) = k(1)*rec_lattice(2,1) + k(2)*rec_lattice(2,2) + k(3)*rec_lattice(2,3)
   q(3) = k(1)*rec_lattice(3,1) + k(2)*rec_lattice(3,2) + k(3)*rec_lattice(3,3) !The value of q is the same as QE   
   
   do i=1,3
      if (abs(q(i)).le.eps12/1000)then
         q(i) = 0.0d0
      end if
   end do   
   
   !write(*,*) rec_lattice
   
   if (nfr1 == 1) then
      nr1x=0
   else
      nr1x = int ( sqrt (geg) / &
                   (sqrt (rec_lattice (1, 1) **2 + rec_lattice (2, 1) **2 + rec_lattice (3, 1) **2) )) + 1
   endif
   if (nfr2 == 1) then
      nr2x=0
   else
      nr2x = int ( sqrt (geg) / &
                   ( sqrt (rec_lattice (1, 2) **2 + rec_lattice (2, 2) **2 + rec_lattice (3, 2) **2) )) + 1
   endif
   if (nfr3 == 1) then
      nr3x=0
   else
      nr3x = int ( sqrt (geg) / &
                   (sqrt (rec_lattice (1, 3) **2 + rec_lattice (2, 3) **2 + rec_lattice (3, 3) **2) )) + 1
   endif
   
   if (loto_2d) then 
      fac = sign*4*twopi/Origin_cell%CellVolume*0.5d0/rec_lattice(3,3)!*alat !CHECK IF THE ANG2BOHRS SHOULD BE HERE OR NOT
      reff=0.0d0
      do i=1,2
         do j=1,2
            reff(i,j)=Diele_Tensor(i,j)*0.5d0*twopi/rec_lattice(3,3) ! (eps)*c/2 in 2pi/a units
         enddo
      enddo
      do i=1,2
         reff(i,i)=reff(i,i)-0.5d0*twopi/rec_lattice(3,3) ! (-1)*c/2 in 2pi/a units
      enddo 
   else
     fac = sign/Origin_cell%CellVolume*e2*2.0d0*twopi
   endif

   do m1 = -nr1x,nr1x
      do m2 = -nr2x,nr2x
         do m3 = -nr3x,nr3x
            
            g1 = m1*rec_lattice(1,1) + m2*rec_lattice(1,2) + m3*rec_lattice(1,3)
            g2 = m1*rec_lattice(2,1) + m2*rec_lattice(2,2) + m3*rec_lattice(2,3)
            g3 = m1*rec_lattice(3,1) + m2*rec_lattice(3,2) + m3*rec_lattice(3,3)

            if (loto_2d) then 
               geg = g1**2 + g2**2 + g3**2
               r=0.0d0
               gp2=g1**2+g2**2
               if (gp2>1.0d-8) then
                  r=g1*reff(1,1)*g1+g1*reff(1,2)*g2+g2*reff(2,1)*g1+g2*reff(2,2)*g2
                  r=r/gp2
               endif
            else
                geg = (g1*(Diele_Tensor(1,1)*g1+Diele_Tensor(1,2)*g2+Diele_Tensor(1,3)*g3)+      &
                   g2*(Diele_Tensor(2,1)*g1+Diele_Tensor(2,2)*g2+Diele_Tensor(2,3)*g3)+      &
                   g3*(Diele_Tensor(3,1)*g1+Diele_Tensor(3,2)*g2+Diele_Tensor(3,3)*g3))
            endif
           
            if (geg > 0.0d0 .and. geg/alph/4.0d0 < gmax ) then
               
               if (loto_2d) then 
                 facgd = fac*exp(-geg/alph/4.0d0)/SQRT(geg)/(1.0+r*SQRT(geg)) 
               else
                 facgd = fac*exp(-geg/alph/4.0d0)/geg
               endif
               !
            
               !arg = twopi* (g1 * (dble(irvec(1,iR))) + g2 * (dble(irvec(2,iR))) + g3 * (dble(irvec(3,iR))))
               do na = 1,Origin_cell%Num_atoms
                  zag(:)=g1*zeu(na,1,:)+g2*zeu(na,2,:)+g3*zeu(na,3,:)
                  fnat(:) = 0.d0
                  do nb = 1,Origin_cell%Num_atoms
                     arg = 2.d0 * Pi* (g1 * (tau(1,na)-tau(1,nb))+             &
                                       g2 * (tau(2,na)-tau(2,nb))+             &
                                       g3 * (tau(3,na)-tau(3,nb)))
                     zcg(:) = g1*zeu(nb,1,:) + g2*zeu(nb,2,:) + g3*zeu(nb,3,:)

                     fnat(:) = fnat(:) + zcg(:)*cos(arg)
                  end do
                  fmtx = MATMUL(RESHAPE(zag, (/3,1/)), RESHAPE(fnat, (/1,3/)))
                  !fmtx = (fmtx + TRANSPOSE(fmtx)) / 2.d0
                  !if (abs((q(1)**2+q(2)**2+q(3)**2)).le.eps12)then
                  !   write(*,*) 'q',q, 'on-site', - facgd * fmtx(i,j)
                  !end if
                  do j=1,3
                     do i=1,3
                        dummydebug(i,na,j,na) = dummydebug(i,na,j,na) - facgd * fmtx(i,j)!zag(CartToOrb(i)) * fnat(CartToOrb(j)) !* (108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(na)*Atom_Mass(na))


                        !Hamk_bulk(3*(na-1)+j,3*(na-1)+i) = Hamk_bulk(3*(na-1)+j,3*(na-1)+i) - &
                        !facgd * zag(i) * fnat(j) * (108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(na)*Atom_Mass(na))
                     end do
                  end do
                  
               end do
            
               
            end if
            
            g1 = g1 + q(1)
            g2 = g2 + q(2)
            g3 = g3 + q(3)

            if (loto_2d) then
               geg = g1**2+g2**2+g3**2
               r=0.0d0
               gp2=g1**2+g2**2
               if (gp2>1.0d-8) then
                  r=g1*reff(1,1)*g1+g1*reff(1,2)*g2+g2*reff(2,1)*g1+g2*reff(2,2)*g2
                  r=r/gp2
               endif
            else
               geg = (g1*(Diele_Tensor(1,1)*g1+Diele_Tensor(1,2)*g2+Diele_Tensor(1,3)*g3)+      &
                   g2*(Diele_Tensor(2,1)*g1+Diele_Tensor(2,2)*g2+Diele_Tensor(2,3)*g3)+      &
                   g3*(Diele_Tensor(3,1)*g1+Diele_Tensor(3,2)*g2+Diele_Tensor(3,3)*g3))
            endif
            
            if (geg > 0.0_DP .and. geg/alph/4.0_DP < gmax ) then
               
               if (loto_2d) then 
                 facgd = fac*exp(-geg/alph/4.0d0)/SQRT(geg)/(1.0+r*SQRT(geg))
               else
                 facgd = fac*exp(-geg/alph/4.0d0)/geg
               endif
               !
            
               do nb = 1,Origin_cell%Num_atoms
                  
                  zbg(:)=g1*zeu(nb,1,:)+g2*zeu(nb,2,:)+g3*zeu(nb,3,:)
                  do na = 1,Origin_cell%Num_atoms
                     zag(:)=g1*zeu(na,1,:)+g2*zeu(na,2,:)+g3*zeu(na,3,:) 
                  
                     
                     arg = 2.d0*Pi* (g1 * (tau(1,na)-tau(1,nb))+             & 
                                       g2 * (tau(2,na)-tau(2,nb))+             &
                                       g3 * (tau(3,na)-tau(3,nb)))
                     
                     
                     
                     facg = facgd * exp(zi*arg)!cmplx(cos(arg), sin(arg), kind=dp)
                     !fmtx = MATMUL(RESHAPE(zag, (/3,1/)), RESHAPE(zbg, (/1,3/)))
                     !fmtx = (fmtx + TRANSPOSE(fmtx)) / 2.d0
                     !if (abs((q(1)**2+q(2)**2+q(3)**2)).le.eps12)then
                     !   write(*,*) 'q',q, 'full', - facgd * fmtx(i,j)
                     !end if
                     do j=1,3
                        do i=1,3
                           
                           dummydebug(i,na,j,nb) = dummydebug(i,na,j,nb) + &
                            facg * zag(i) * zbg(j)!*(108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(na)*Atom_Mass(nb))
                           !Hamk_bulk(3*(na-1)+j,3*(nb-1)+i) = Hamk_bulk(3*(na-1)+j,3*(nb-1)+i) + &
                            !facg * zag(j) * zbg(i)*(108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(na)*Atom_Mass(nb))
                        end do
                     end do
                  end do
               end do
            end if
         end do
      end do
   end do

   !dummydebug2 = RESHAPE(dummydebug,(/Num_wann,Num_wann/))
   !dummydebug2 = 0.0d0
   do nb=1,Origin_cell%Num_atoms
      do na=1,Origin_cell%Num_atoms
         do j=1,3
            do i=1,3
               dummydebug2(3*(na-1)+i,3*(nb-1)+j) = dummydebug(i,na,j,nb)*(108.97077184367376*eV2Hartree)**2!/sqrt(Atom_Mass(na)*Atom_Mass(nb))
            end do
         end do
      end do
   end do


   !do na=1,Num_wann
   !   do nb=1, Num_wann
   !      dummydebug2(na,nb) = dummydebug2(na,nb)*(108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(na))*Atom_Mass(Origin_cell%spinorbital_to_atom_index(nb)))
   !   end do
   !end do

   !dummyFC(:,:,1) = dummydebug2
   !> Impose ASR on NAC
   !call impose_ASR2(k,dummydebug2,Num_wann,1)


   
   !write(*,*) 'start', SUM(dummydebug2,1)
   !if (abs(k(1)**2+k(2)**2+k(3)**2) < eps12) then
   !   write(*,*) q
   !   write(*,*) dummydebug(3,4,1,5)
   !end if
   !dummydebug2 = dummydebug2 - SUM(dummydebug2,2)/Num_wann
   !write(*,*) 'start', SUM(dummydebug2,2)
   do na=1,Num_wann
      do nb=1, Num_wann
         Hamk_bulk(na,nb) =  dummydebug2(na,nb) + Hamk_bulk(na,nb)
         !*(108.97077184367376*eV2Hartree)**2/sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(na))*Atom_Mass(Origin_cell%spinorbital_to_atom_index(nb)))
      end do
   end do
   !write(*,*) 'dummydebug', real(dummydebug(2,1,2,1))
   !write(*,*) 'dummydebug2',real(dummydebug2(2,2))
   return
end subroutine long_range_phonon_interaction

subroutine impose_ASR2(k, ham,nwann, slabs)
   ! This subroutine imposes the Acoustic Sum Rule to a Hamiltonian,
   ! based on CellConstructor's python code
   use para
   
   implicit none

   integer :: nwann, slabs
   complex(Dp), intent(out) :: ham(nwann*slabs, nwann*slabs)
   real(Dp) :: k(3)


   real(dp), external :: norm

   !> local variables
   real(Dp) :: bg(3,3), &
   Q_proj(nwann*slabs,nwann*slabs), &! ASR projector
   v1(nwann*slabs),  &
   rec_lattice(3,3),    &
   lattice_ang(3,3),    &
   q(3),                &
   tol,                 &
   f_q2i(3),            &
   f_q2,                &
   mask(3),             &
   atdotq,              &
   normvec,             &
   sigma,               &    
   sqrtmass      
   
   complex(Dp) :: asrCorrection(nwann*slabs,nwann*slabs), dummyham(nwann*slabs,nwann*slabs)

   integer :: N_i(3)

   !> loop indices
   integer :: i,j,l
   
   
   

   tol = 1e-8

   q = 0.0d0
   rec_lattice = Origin_cell%reciprocal_lattice*Origin_cell%cell_parameters(1)/(twopi) !This value is the same as QE
   
   q(1) = k(1)*rec_lattice(1,1) + k(2)*rec_lattice(1,2) + k(3)*rec_lattice(1,3)
   q(2) = k(1)*rec_lattice(2,1) + k(2)*rec_lattice(2,2) + k(3)*rec_lattice(2,3)
   q(3) = k(1)*rec_lattice(3,1) + k(2)*rec_lattice(3,2) + k(3)*rec_lattice(3,3) !The value of q is the same as QE

   q = -q(:)*Ang2Bohr/Origin_cell%cell_parameters(1)
   lattice_ang =Origin_cell%lattice/Ang2Bohr

   
   !> Create the ASR projector
   Q_proj = 0.0d0
   do i =1,3
      v1 = 0.0d0
      normvec = 0.0d0
      do j=1,nwann*slabs,3
         v1(j+(i-1))=1.0d0
         normvec = normvec + 1.0d0
      end do
      v1 = v1/sqrt(normvec)
      
      Q_proj = Q_proj + spread(source = v1, dim = 2, ncopies = nwann*slabs) * spread(source = v1, dim = 1, ncopies = nwann*slabs)
   end do
   !write(*,*) Q_proj
   do i=1,3
      N_i(i) = 2*abs(maxval(irvec(i,:))-minval(irvec(i,:))) + 1
      if (MOD(N_i(i),2) .eq. 0) then 
         N_i(i) = N_i(i) + 1
      end if
   end do
   
   f_q2i = 1.0d0
   
   mask = (/(lattice_ang(1,1)*q(1)+lattice_ang(2,1)*q(2)+lattice_ang(3,1)*q(3)),      &
           (lattice_ang(1,2)*q(1)+lattice_ang(2,2)*q(2)+lattice_ang(3,2)*q(3)),      &
           (lattice_ang(1,3)*q(1)+lattice_ang(2,3)*q(2)+lattice_ang(3,3)*q(3))/)
   
   sigma = 2e-1
   !write(*,*) abs(sin(Pi*mask)), q
   !> We use the mask to avoid division by zero because the limit 0/0 is 1 in this case
   do i=1,3
      if (abs(sin(Pi*mask(i))) > tol) then
         atdotq = lattice_ang(1,i)*q(1)+lattice_ang(2,i)*q(2)+lattice_ang(3,i)*q(3)
         f_q2i(i) = sin(N_i(i)*Pi*atdotq)/(N_i(i)*sin(Pi*atdotq))
      end if
   end do
   f_q2 = f_q2i(1)*f_q2i(2)*f_q2i(3)
   !f_q2 = exp( - (q(1)**2+q(2)**2+q(3)**2) / (2 * sigma**2))
   
   asrCorrection = 0.0d0
   !> Go to proper units
   do i=1, nwann*slabs
      do j=1, nwann*slabs
         if (MOD(i,nwann).eq.0)then
            sqrtmass = sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(nwann)))
         else 
            sqrtmass = sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(MOD(i,nwann))))
         end if
         if (MOD(j,nwann).eq.0)then
            sqrtmass = sqrtmass *sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(nwann)))
         else 
            sqrtmass = sqrtmass*sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(MOD(j,nwann))))
         end if
         dummyham(i,j) = ham(i,j)/(108.97077184367376*eV2Hartree)**2*sqrtmass
      end do
   end do
   
   !> impose the acoustic sum rule
   do i=1, nwann*slabs
      do j=1, nwann*slabs
         do l=1, nwann*slabs
            asrCorrection(i,j) = asrCorrection(i,j) - dummyham(i,l)*Q_proj(j,l) * f_q2 !/sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(i))*Atom_Mass(Origin_cell%spinorbital_to_atom_index(j)))
         end do
      end do
   end do
   
   dummyham = dummyham + asrCorrection
   asrCorrection = 0.0d0
   do i=1, nwann*slabs
      do j=1, nwann*slabs
         do l=1, nwann*slabs
            asrCorrection(i,j) = asrCorrection(i,j) - dummyham(l,j)*Q_proj(i,l) * f_q2 !/sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(i))*Atom_Mass(Origin_cell%spinorbital_to_atom_index(j)))
         end do
      end do
   end do
   
   dummyham = dummyham + asrCorrection
   !> Go back to WT units
   do i=1, nwann*slabs
      do j=1, nwann*slabs
         if (MOD(i,nwann).eq.0)then
            sqrtmass = sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(nwann)))
         else 
            sqrtmass = sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(MOD(i,nwann))))
         end if
         if (MOD(j,nwann).eq.0)then
            sqrtmass = sqrtmass *sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(nwann)))
         else 
            sqrtmass = sqrtmass*sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(MOD(j,nwann))))
         end if
         ham(i,j) = dummyham(i,j)*(108.97077184367376*eV2Hartree)**2/sqrtmass
      end do
   end do
   !write(*,*) dummyham - Hamk_bulk, q
end subroutine impose_ASR2

subroutine impose_ASR_on_eff_charges (nat, tau, zeu)
   !-----------------------------------------------------------------------
   !
   ! Copyright (C) 2001-2012 Quantum ESPRESSO group
   ! This file is distributed under the terms of the
   ! GNU General Public License. See the file `License'
   ! in the root directory of the present distribution,
   ! or http://www.gnu.org/copyleft/gpl.txt .
   !
   ! Based on the QE code, modified by Francesc Ballester
   !
   !-----------------------------------------------------------------------
   implicit none
   integer, intent(in) :: nat
   double precision, intent(in) :: tau(3,nat)
   double precision, intent(inout) :: zeu(nat,3,3) !changed from (3,3,nat)
   !double complex, intent(inout) :: dyn(3,3,nat,nat)
   !
   integer :: i,j,n,m,p,k,l,q,r,na, nb, na1, i1, j1
   double precision, allocatable:: dynr_new(:,:,:,:,:), zeu_new(:,:,:)
   !double precision, allocatable :: u(:,:,:,:,:)
   ! These are the "vectors" associated with the sum rules
   !
   !integer u_less(6*3*nat)!,n_less,i_less
   ! indices of the vectors u that are not independent to the preceding ones,
   ! n_less = number of such vectors, i_less = temporary parameter
   !
   integer, allocatable :: ind_v(:,:,:)
   double precision, allocatable :: v(:,:)
   ! These are the "vectors" associated with symmetry conditions, coded by
   ! indicating the positions (i.e. the four indices) of the non-zero elements
   ! (there should be only 2 of them) and the value of that element.
   ! We do so in order to use limit the amount of memory used.
   !
   !double precision, allocatable :: w(:,:,:,:), x(:,:,:,:)
   double precision sum, scal, norm2
   ! temporary vectors and parameters
   !
   double precision, allocatable :: zeu_u(:,:,:,:)
   ! These are the "vectors" associated with the sum rules on effective charges
   !
   integer zeu_less(6*3),nzeu_less,izeu_less
   ! indices of the vectors zeu_u that are not independent to the preceding
   ! ones, nzeu_less = number of such vectors, izeu_less = temporary parameter
   !
   double precision, allocatable :: zeu_w(:,:,:), zeu_x(:,:,:)

   character(len=10) :: asr='crystal'
   ! temporary vectors
   !
   ! Initialization
   ! n is the number of sum rules to be considered (if asr.ne.'simple')
   ! 'axis' is the rotation axis in the case of a 1D system (i.e. the rotation
   !  axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
   !
   !if ( (asr.ne.'simple') .and. (asr.ne.'crystal') .and. (asr.ne.'one-dim') &
   !                       .and.(asr.ne.'zero-dim')) then
   !   call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
   !endif
   if(asr.eq.'crystal') n=3
   if(asr.eq.'one-dim') then
   !   write(6,'("asr rotation axis in 1D system= ",I4)') axis
      n=4
   endif
   if(asr.eq.'zero-dim') n=6
   !
   ! ASR on effective charges
   !
   if(asr.eq.'simple') then
      do i=1,3
         do j=1,3
            sum=0.0d0
            do na=1,nat
               sum = sum + zeu(na,i,j)
            end do
            do na=1,nat
               zeu(na,i,j) = zeu(na,i,j) - sum/nat
            end do
         end do
      end do
   else
      ! generating the vectors of the orthogonal of the subspace to project
      ! the effective charges matrix on
      !
      allocate ( zeu_new(3,3,nat) )
      allocate ( zeu_u(6*3,3,3,nat) )
      zeu_u(:,:,:,:)=0.0d0
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu_new(i,j,na)=zeu(na,i,j)
            enddo
         enddo
      enddo
      !
      p=0
      do i=1,3
         do j=1,3
            ! These are the 3*3 vectors associated with the
            ! translational acoustic sum rules
            p=p+1
            zeu_u(p,i,j,:)=1.0d0
            !
         enddo
      enddo
      !
      !if (n.eq.4) then
      !   do i=1,3
            ! These are the 3 vectors associated with the
            ! single rotational sum rule (1D system)
      !      p=p+1
      !      do na=1,nat
      !         zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
      !         zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
      !      enddo
            !
      !   enddo
      !endif
      
      if (n.eq.6) then
         do i=1,3
            do j=1,3
               ! These are the 3*3 vectors associated with the
               ! three rotational sum rules (0D system - typ. molecule)
               p=p+1
               do na=1,nat
                  zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                  zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
               enddo
               !
            enddo
         enddo
      endif
      !
      ! Gram-Schmidt orthonormalization of the set of vectors created.
      !
      allocate ( zeu_w(3,3,nat), zeu_x(3,3,nat) )
      nzeu_less=0
      do k=1,p
         zeu_w(:,:,:)=zeu_u(k,:,:,:)
         zeu_x(:,:,:)=zeu_u(k,:,:,:)
         do q=1,k-1
            r=1
            do izeu_less=1,nzeu_less
               if (zeu_less(izeu_less).eq.q) r=0
            enddo
            if (r.ne.0) then
               call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
               zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
            endif
         enddo
         call sp_zeu(zeu_w,zeu_w,nat,norm2)
         if (norm2.gt.1.0d-16) then
            zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
         else
            nzeu_less=nzeu_less+1
            zeu_less(nzeu_less)=k
         endif
      enddo
      !
      !
      ! Projection of the effective charge "vector" on the orthogonal of the
      ! subspace of the vectors verifying the sum rules
      !
      zeu_w(:,:,:)=0.0d0
      do k=1,p
         r=1
         do izeu_less=1,nzeu_less
            if (zeu_less(izeu_less).eq.k) r=0
         enddo
         if (r.ne.0) then
            zeu_x(:,:,:)=zeu_u(k,:,:,:)
            call sp_zeu(zeu_x,zeu_new,nat,scal)
            zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
         endif
      enddo
      !
      ! Final substraction of the former projection to the initial zeu, to get
      ! the new "projected" zeu
      !
      zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
      call sp_zeu(zeu_w,zeu_w,nat,norm2)
      !write(6,'(5x,"Acoustic Sum Rule: || Z*(ASR) - Z*(orig)|| = ",E15.6)') &
      !     SQRT(norm2)
      !
      ! Check projection
      !
      !write(6,'("Check projection of zeu")')
      !do k=1,p
      !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
      !  call sp_zeu(zeu_x,zeu_new,nat,scal)
      !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
      !enddo
      !
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu(na,i,j)=zeu_new(i,j,na)
            enddo
         enddo
      enddo
      deallocate (zeu_w, zeu_x)
      deallocate (zeu_u)
      deallocate (zeu_new)
   endif

end subroutine impose_ASR_on_eff_charges



subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
   !-----------------------------------------------------------------------
   !
   ! does the scalar product of two effective charges matrices zeu_u and zeu_v
   ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
   !
   implicit none
   integer i,j,na,nat
   double precision zeu_u(3,3,nat)
   double precision zeu_v(3,3,nat)
   double precision scal
   !
   !
   scal=0.0d0
   do i=1,3
     do j=1,3
       do na=1,nat
         scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
       enddo
     enddo
   enddo
   !
   return
   !
 end subroutine sp_zeu

subroutine ham_bulk_LOTO(k,Hamk_bulk)
   ! This subroutine caculates Hamiltonian for
   ! bulk system without the consideration of the atom's position
   ! with the LOTO correction for phonon system
   !
   ! History
   !
   !        July/15/2017 by TianTian Zhang
   !
   !        May/8/2024 by Francisco M. Ballester
   !          - Bug correction and propper LO-TO implementation.

   use para
   implicit none

   ! loop index
   integer :: i1,i2,ia,ib,ic,iR
   integer  :: ii,jj,mm,nn,pp,qq,xx,yy,zz

   real(dp) :: kdotr

   ! wave vector in 2d
   real(Dp) :: k(3), keps(3) = (/eps12,eps12,eps12/)

   ! coordinates of R vector
   real(Dp) :: R(3)
   complex(dp) :: ratio

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   !> see eqn. (3) in J. Phys.: Condens. Matter 22 (2010) 202201
   !complex(Dp),allocatable :: nac_correction(:, :)

   complex(dp), allocatable :: mat1(:, :)
   complex(dp), allocatable :: mat2(:, :)
   real(dp) :: temp1(3)=(/0.0,0.0,0.0/)
   real(dp) :: temp2=0.0
   real(dp) :: temp3(30),constant_t
   real(dp) ::A_ii(3)=(/0.0,0.0,0.0/)
   real(dp) ::A_jj(3)=(/0.0,0.0,0.0/), zag(3), zbg(3), qeq

   !> k times Born charge
   real(dp), allocatable :: kBorn(:, :)
   real(dp), allocatable :: tau(:,:), zeu(:,:,:)

   complex(dp) :: nac_q

   integer :: CartToOrb(3)

   !> check if we are exactly at gamma
   logical :: atGamma=.false.
   

   complex(Dp), allocatable :: dummydynmat(:,:,:)
   complex(Dp), allocatable :: dummydynmat2(:,:,:)

   allocate(kBorn(Origin_cell%Num_atoms, 3))
   allocate(mat1(Num_wann, Num_wann))
   allocate(mat2(Num_wann, Num_wann))
   allocate(dummydynmat(Num_wann,Num_wann,1))
   allocate(dummydynmat2(Num_wann,Num_wann,1))
   allocate(tau(3,Origin_cell%Num_atoms))
   allocate(zeu(Origin_cell%Num_atoms,3,3))
   
   !mat2 = 0d0
   !nac_correction= zzero

   CartToOrb=(/1,2,3/)

   Hamk_bulk=0d0
   dummydynmat2=0d0

   tau(:,:) = Origin_cell%Atom_position_cart/Origin_cell%cell_parameters(1)
   !write(*,*) Origin_cell%Atom_name
   ! do ia = 1,Origin_cell%Num_atoms
   !    do ib=1,3
   !       do ic=1,3
   !          tau(ib,ia) = tau(ib,ia) + (Origin_cell%Atom_position_direct(ic,ia)*Origin_cell%lattice(ib,ic))/Origin_cell%cell_parameters(1) !Now we have the same units as QE
   !       end do
   !    end do
   ! end do
   
   zeu(:,:,:) = Born_Charge(:,:,:)

   !write(*,*) 'start',Origin_cell%proj_name

   !write(*,*) 'start', SUM(zeu,1)

   !call impose_ASR_on_eff_charges(Origin_cell%Num_atoms,tau,zeu)
   !write(*,*) 'start',zeu - Born_Charge
   !>  add loto splitting term
   temp1(1:3)= (/0.0,0.0,0.0/)
   temp2= 0.0
   atGamma=.false.
   
   
   
   

   do iR=1, Nrpts
      R(1)=dble(irvec(1,iR))
      R(2)=dble(irvec(2,iR))
      R(3)=dble(irvec(3,iR))
      kdotr=k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
      
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      do mm = 1, Num_wann
         pp = Origin_cell%spinorbital_to_atom_index(mm)
         ii = Origin_cell%spinorbital_to_projector_index(mm)
         do nn = 1, Num_wann
            qq = Origin_cell%spinorbital_to_atom_index(nn)
            jj = Origin_cell%spinorbital_to_projector_index(nn)
            Hamk_bulk(mm, nn)= Hamk_bulk(mm, nn) + (HmnR(3*(pp-1)+ii,3*(qq-1)+jj,iR))*ratio/ndegen(iR)*SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
         end do
      end do
   enddo ! iR

   !call ham_bulk_latticegauge(k,Hamk_bulk)
   
   ! Add long-range interaction
   !Hamk_bulk = 0.0d0
   call long_range_phonon_interaction(0,0,0,k,.false.,1.0d0,Hamk_bulk,tau,zeu)
   
   if (abs((k(1)**2+k(2)**2+k(3)**2)).le.eps12)then  !> skip k=0
      atGamma=.true.

      qeq = (keps(1)*(Diele_Tensor(1,1)*keps(1)+Diele_Tensor(1,2)*keps(2)+Diele_Tensor(1,3)*keps(3))+    &
            keps(2)*(Diele_Tensor(2,1)*keps(1)+Diele_Tensor(2,2)*keps(2)+Diele_Tensor(2,3)*keps(3))+    &
            keps(3)*(Diele_Tensor(3,1)*keps(1)+Diele_Tensor(3,2)*keps(2)+Diele_Tensor(3,3)*keps(3)))
      
      constant_t= 2.0d0*4.0d0*Pi/Origin_cell%CellVolume

   endif
   
   
   
   
   
   mat2 = 0.0d0
   if (atGamma) then
      do pp = 1,Origin_cell%Num_atoms
         do qq = 1,Origin_cell%Num_atoms
            do ii=1,3
               zag(ii) = keps(1)*zeu(pp,1,ii) +  keps(2)*zeu(pp,2,ii) + keps(3)*zeu(pp,3,ii)
               
               zbg(ii) = keps(1)*zeu(qq,1,ii) +  keps(2)*zeu(qq,2,ii) + keps(3)*zeu(qq,3,ii)

            end do
            do ii=1,3
               do jj=1,3
                  
                  nac_q= constant_t*zag(ii)*zbg(jj)/qeq!/sqrt(Atom_Mass(pp)*Atom_Mass(qq))!kBorn(jj, CartToOrb(qq))*kBorn(ii, CartToOrb(pp))*constant_t/sqrt(Atom_Mass(ii)*Atom_Mass(jj))
                  !write(*,*) real(nac_q)
                  mat2(3*(pp-1)+ii,3*(qq-1)+jj) = nac_q!*(108.97077184367376*eV2Hartree)**2!/SQRT(Atom_Mass(pp)*Atom_Mass(qq))
               
               enddo  ! jj
            enddo  ! ii
         enddo ! qq
      enddo  ! pp
      
   end if
   

   !call impose_ASR2(k,mat2,Num_wann,1)


   do ii=1,Num_wann
      do jj=1, Num_wann
         Hamk_bulk(ii,jj) = Hamk_bulk(ii,jj) +mat2(ii,jj)*(108.97077184367376*eV2Hartree)**2!/sqrt(Atom_Mass(Origin_cell%spinorbital_to_atom_index(na))*Atom_Mass(Origin_cell%spinorbital_to_atom_index(nb)))
      end do
   end do

   do ii=1,Num_wann
      do jj=1, Num_wann
         pp = Origin_cell%spinorbital_to_atom_index(ii)
         qq = Origin_cell%spinorbital_to_atom_index(jj)
         Hamk_bulk(ii,jj) = Hamk_bulk(ii,jj)/SQRT(Atom_Mass(pp)*Atom_Mass(qq)) 
      end do
   end do


   
   
   

   ! Apply ASR
   !call impose_ASR2(k,Hamk_bulk,Num_wann, 1)

   ! preserve previous k point for NAC calculation in Gamma so as to not have unnecessary discontinuities
   if ((k(1).ne.0.0d0) .or. (k(2).ne.0.0d0) .or. (k(3).ne.0.0d0)) then
      !rec_lattice = Origin_cell%reciprocal_lattice*Origin_cell%cell_parameters(1)/(twopi)
      keps(1) = k(1)*Origin_cell%reciprocal_lattice(1,1) + k(2)*Origin_cell%reciprocal_lattice(1,2) + k(3)*Origin_cell%reciprocal_lattice(1,3)
      keps(2) = k(1)*Origin_cell%reciprocal_lattice(2,1) + k(2)*Origin_cell%reciprocal_lattice(2,2) + k(3)*Origin_cell%reciprocal_lattice(2,3)
      keps(3) = k(1)*Origin_cell%reciprocal_lattice(3,1) + k(2)*Origin_cell%reciprocal_lattice(3,2) + k(3)*Origin_cell%reciprocal_lattice(3,3)
      keps = keps*Origin_cell%cell_parameters(1)/(twopi)
      do ii=1,3
         if (abs(keps(ii)).le.eps12/1000)then
            keps(ii) = 0.0d0
         end if
      end do  
   end if
   
   

   deallocate(kBorn)
   deallocate(mat1)
   deallocate(mat2)
   deallocate(dummydynmat)
   deallocate(tau)
   return
end subroutine ham_bulk_LOTO

subroutine ham_bulk_kp(k,Hamk_bulk)
   ! > construct the kp model at K point of WC system
   ! > space group is 187. The little group is C3h
   ! Sep/10/2018 by Quansheng Wu

   use para, only : Num_wann, dp, stdout, zi
   implicit none

   integer :: i1,i2,ia,ib,ic,iR, nwann

   ! coordinates of R vector
   real(Dp) :: R(3), R1(3), R2(3), kdotr, kx, ky, kz
   real(dp) :: A1, A2, B1, B2, C1, C2, D1, D2
   real(dp) :: m1x, m2x, m3x, m4x, m1z, m2z, m3z, m4z
   real(dp) :: E1, E2, E3, E4

   complex(dp) :: factor, kminus, kplus

   real(dp), intent(in) :: k(3)

   ! Hamiltonian of bulk system
   complex(Dp),intent(out) ::Hamk_bulk(Num_wann, Num_wann)

   if (Num_wann/=4) then
      print *, "Error : in this kp model, num_wann should be 4"
      print *, 'Num_wann', Num_wann
      stop
   endif

   kx=k(1)
   ky=k(2)
   kz=k(3)
   E1= 1.25d0
   E2= 0.85d0
   E3= 0.25d0
   E4=-0.05d0

   A1= 0.10d0
   A2= 0.30d0
   B1= 0.05d0
   B2= 0.1d0

   C1= -1.211d0
   C2= 1.5d0
   D1=-0.6d0
   D2= 0.6d0

   m1x= -1.8d0
   m2x= -1.6d0
   m3x=  1.2d0
   m4x=  1.4d0
   m1z= 2d0
   m2z= 3d0
   m3z=-1d0
   m4z=-1d0

   kminus= kx- zi*ky
   kplus= kx+ zi*ky

   Hamk_bulk= 0d0
   !> kp
   Hamk_bulk(1, 1)= E1+ m1x*(kx*kx+ky*ky)+ m1z*kz*kz
   Hamk_bulk(2, 2)= E2+ m2x*(kx*kx+ky*ky)+ m2z*kz*kz
   Hamk_bulk(3, 3)= E3+ m3x*(kx*kx+ky*ky)+ m3z*kz*kz
   Hamk_bulk(4, 4)= E4+ m4x*(kx*kx+ky*ky)+ m4z*kz*kz

   Hamk_bulk(1, 2)=-zi*D1*kplus*kz
   Hamk_bulk(2, 1)= zi*D1*kminus*kz
   Hamk_bulk(3, 4)= zi*D2*kminus*kz
   Hamk_bulk(4, 3)=-zi*D2*kplus*kz

   Hamk_bulk(1, 4)= zi*C1*kz
   Hamk_bulk(2, 3)= zi*C2*kz
   Hamk_bulk(3, 2)=-zi*C2*kz
   Hamk_bulk(4, 1)=-zi*C1*kz

   Hamk_bulk(1, 3)=  A1*kplus+ B1*kminus*kminus !+ D*kplus*kplus*kplus
   Hamk_bulk(2, 4)=  A2*kminus+ B2*kplus*kplus !+ D*kminus*kminus*kminus
   Hamk_bulk(3, 1)= conjg(Hamk_bulk(1, 3))
   Hamk_bulk(4, 2)= conjg(Hamk_bulk(2, 4))

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
            write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
            !stop
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_kp

subroutine ham_bulk_coo_densehr(k,nnzmax, nnz, acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   ! History
   !        Dec/10/2018 by Quansheng Wu
   use para
   implicit none

   real(dp), intent(in) :: k(3)
   integer, intent(in) :: nnzmax
   integer, intent(out) :: nnz
   integer,intent(inout):: icoo(nnzmax),jcoo(nnzmax)
   complex(dp),intent(inout) :: acoo(nnzmax)

   ! loop index
   integer :: i1,i2,ia,ib,ic,iR, iz
   integer :: nwann

   real(dp) :: kdotr

   ! wave vector in 3D BZ
   ! coordinates of R vector
   real(Dp) :: R(3), R1(3), R2(3)

   complex(dp) :: factor

   ! Hamiltonian of bulk system
   complex(Dp), allocatable :: Hamk_bulk(:, :)

   allocate( Hamk_bulk(Num_wann, Num_wann))
   Hamk_bulk= 0d0

   do iR=1, Nrpts
      ia=irvec(1,iR)
      ib=irvec(2,iR)
      ic=irvec(3,iR)

      R(1)=dble(ia)
      R(2)=dble(ib)
      R(3)=dble(ic)
      kdotr=k(1)*R (1) + k(2)*R (2) + k(3)*R (3)
      factor= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))

      Hamk_bulk(:, :)=&
         Hamk_bulk(:, :) &
         + HmnR(:, :, iR)*factor/ndegen(iR)
   enddo ! iR

   !> transform Hamk_bulk into sparse COO format

   nnz= 0
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)).ge.1e-6)then
            nnz= nnz+ 1
            if (nnz>nnzmax) stop 'nnz is larger than nnzmax in ham_bulk_coo_densehr'
            icoo(nnz) = i1
            jcoo(nnz) = i2
            acoo(nnz) = Hamk_bulk(i1, i2)
         endif
      enddo
   enddo

   ! check hermitcity
   do i1=1, Num_wann
      do i2=1, Num_wann
         if(abs(Hamk_bulk(i1,i2)-conjg(Hamk_bulk(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with Hamk_bulk'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', Hamk_bulk(i1, i2)
            write(stdout,*)'value at (i2, i1)', Hamk_bulk(i2, i1)
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_coo_densehr


subroutine ham_bulk_coo_sparsehr_latticegauge(k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout):: icoo(splen),jcoo(splen)!,splen
   complex(dp),intent(inout) :: acoo(splen)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1,splen
      icoo(i)=hicoo(i)
      jcoo(i)=hjcoo(i)
      posij= hirv(1:3, i)
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*hacoo(i)
   end do

   return
end subroutine ham_bulk_coo_sparsehr_latticegauge

subroutine ham_bulk_coo_sparsehr(k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout):: icoo(splen),jcoo(splen)!,splen
   complex(dp),intent(inout) :: acoo(splen)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1, splen
      icoo(i)=hicoo(i)
      jcoo(i)=hjcoo(i)
      posij= hirv(1:3, i)+ Origin_cell%wannier_centers_direct(:, jcoo(i))- Origin_cell%wannier_centers_direct(:, icoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*hacoo(i)
   end do

   return
end subroutine ham_bulk_coo_sparsehr


subroutine overlap_bulk_coo_sparse(k, acoo, icoo, jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer,intent(inout) :: icoo(splen_overlap_input),jcoo(splen_overlap_input)
   complex(dp),intent(inout) :: acoo(splen_overlap_input)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1, splen_overlap_input
      icoo(i)=sicoo(i)
      jcoo(i)=sjcoo(i)
      posij= sirv(1:3, i)+ Origin_cell%wannier_centers_direct(:, sjcoo(i))- Origin_cell%wannier_centers_direct(:, sicoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*sacoo(i)
   end do

   return
end subroutine overlap_bulk_coo_sparse

subroutine ham_bulk_tpaw(mdim, k, Hamk_moire)
   ! This subroutine caculates Hamiltonian for
   ! a twist system
   !
   !  Lattice gauge Hl
   !  Atomic gauge Ha= U* Hl U 
   !  where U = e^ik.wc(i) on diagonal
   ! contact wuquansheng@gmail.com
   ! 2024 Jan 19

   use para
   implicit none

   !> leading dimension of Ham_moire
   integer, intent(in) :: mdim
   real(dp), intent(in) :: k(3)
   complex(dp), intent(inout) :: hamk_moire(mdim, mdim)

   ! wave vector in 3d
   real(Dp) :: kdotr, pos0(3), dis

   complex(dp) :: factor

   real(dp) :: pos(3), pos1(3), pos2(3), pos_cart(3), pos_direct(3)

   integer :: Num_gvectors

   !> dimension 3*Num_gvectors 
   real(dp), allocatable :: gvectors(:, :)

   integer :: i1, i2

   !> for a test, we assume a twist homo-bilayer system
   ! mdim= Num_gvectors* Folded_cell%NumberOfspinorbitals*2

   hamk_moire=0d0

   !> first, set G vectors

   !> construct T matrix
   !> the dimension of the T matrix is num_wann
   

   ! check hermitcity
   do i1=1, mdim
      do i2=1, mdim
         if(abs(hamk_moire(i1,i2)-conjg(hamk_moire(i2,i1))).ge.1e-6)then
            write(stdout,*)'there is something wrong with hamk_moire'
            write(stdout,*)'i1, i2', i1, i2
            write(stdout,*)'value at (i1, i2)', hamk_moire(i1, i2)
            write(stdout,*)'value at (i2, i1)', hamk_moire(i2, i1)
            !stop
         endif
      enddo
   enddo

   return
end subroutine ham_bulk_tpaw


subroutine valley_k_coo_sparsehr(nnz, k,acoo,icoo,jcoo)
   !> This subroutine use sparse hr format
   !> Here we use atomic gauge which means the atomic position is taken into account
   !> in the Fourier transformation
   use para
   implicit none

   real(dp) :: k(3), posij(3)
   real(dp) :: kdotr
   integer, intent(in) :: nnz
   integer,intent(inout) :: icoo(nnz),jcoo(nnz)
   complex(dp),intent(inout) :: acoo(nnz)
   complex(dp) ::  ratio
   integer :: i,j,ir

   do i=1,nnz
      ir= valley_operator_irv(i)
      icoo(i)=valley_operator_icoo(i)
      jcoo(i)=valley_operator_jcoo(i)
      posij=irvec_valley(:, ir)+ Origin_cell%wannier_centers_direct(:, jcoo(i))- Origin_cell%wannier_centers_direct(:, icoo(i))
      kdotr=posij(1)*k(1)+posij(2)*k(2)+posij(3)*k(3)
      ratio= (cos(twopi*kdotr)+zi*sin(twopi*kdotr))
      acoo(i)=ratio*valley_operator_acoo(i)
   enddo

   return
end subroutine valley_k_coo_sparsehr


subroutine rotation_to_Ham_basis(UU, mat_wann, mat_ham)
   !> this subroutine rotate the matrix from Wannier basis to Hamiltonian basis
   !> UU are the eigenvectors from the diagonalization of Hamiltonian
   !> mat_ham=UU_dag*mat_wann*UU
   use para, only : dp, Num_wann
   implicit none
   complex(dp), intent(in) :: UU(Num_wann, Num_wann)
   complex(dp), intent(in) :: mat_wann(Num_wann, Num_wann)
   complex(dp), intent(out) :: mat_ham(Num_wann, Num_wann)
   complex(dp), allocatable :: UU_dag(:, :), mat_temp(:, :)

   allocate(UU_dag(Num_wann, Num_wann), mat_temp(Num_wann, Num_wann))
   UU_dag= conjg(transpose(UU))

   call mat_mul(Num_wann, mat_wann, UU, mat_temp) 
   call mat_mul(Num_wann, UU_dag, mat_temp, mat_ham) 

   return
end subroutine rotation_to_Ham_basis




! ============================== INITIALIZATIONS ============================================

subroutine initialize_perm(R2, n_blocks,SClat,PBC)
   use para
   implicit none
   integer, intent(in)           :: R2(3,n_blocks),n_blocks
   integer, intent(in)           :: SCLat(3)
   logical, intent(in)           :: PBC
   !
   integer :: x(3), y(3),i_block,j_block
   logical :: found, Geq
  
   write(*,*) " "
   write(*,*) "Initialize indices for permutation...", perm_initialized

      
   allocate(P(n_blocks))

   do i_block = 1, n_blocks
               
       x=R2(:,i_block)
       found=.false.
       
       do j_block = 1, n_blocks
       !    
           y=R2(:,j_block)
           if (PBC) then
               Geq= ( ALL (mod(y+x,SCLat)==0) )
            else
               Geq= ( ALL (y+x ==0) )
           end if
           if ( Geq) then
               P(i_block)=j_block
               found=.true.
               cycle
           end if
       !    
       end do
       
       if ( .not. found ) then
       
           write(*,*) " ERROR: new vector found during the permutation symmetry initialization "
           write(*,*) "        the execution stops here. "
           write(*,*) "        If the 2ndFC is centered, try to repeat the centering with higher Far "            
           write(*,*) "        If the 2ndFC is not centered, try to repeat the ASR imposition with PBC=True "            
           stop
           
       end if
   
   end do

   perm_initialized =.true.

   write(*,*) "Initialize indices for permutation...done"
   write(*,*) " "
   
end subroutine initialize_perm

subroutine clear_all()
   use para
   implicit none

   ! This subroutine deallocates the permutations
   if (perm_initialized) then
       perm_initialized = .false.
       deallocate(P)
   end if
end subroutine clear_all

!==========================================================================================

!function Geq(v1,v2,Lat,PBC)
!   implicit none
!   logical    :: Geq
!   integer, dimension(:), intent(in)           :: v1,v2
!   integer, dimension(:), intent(in) :: Lat
!   logical, intent(in) :: PBC
   !    
!   if (PBC) then
!    Geq= ( ALL (mod(v1-v2,Lat)==0) )
!   else
!    Geq= ( ALL (v1-v2 ==0) )
!   end if
   !
   !   
!end function Geq

!=============================== PERM SYM =========================================

subroutine impose_perm_sym(FC,R2,SClat,PBC,verbose,FCvar,FC_sym,nat,n_blocks)
   use para
   implicit none
   !integer, parameter :: DP = selected_real_kind(14,200)
   complex(kind=DP), intent(in)     :: FC(3*nat,3*nat,n_blocks)
   integer, intent(in)           :: nat, n_blocks
   integer, intent(in)           :: R2(3,n_blocks)
   integer, intent(in)           :: SClat(3)
   logical, intent(in)           :: PBC, verbose
   !
   complex(kind=DP), intent(out)    :: FC_sym(3*nat,3*nat,n_blocks)
   real(kind=DP), intent(out) :: FCvar
   !
   integer                       :: jn1, jn2, i_block,nat3
   !

   
   if ( .not. perm_initialized ) call initialize_perm(R2,n_blocks,SClat,PBC)
         
   nat3=3*nat  
     
   do i_block = 1,n_blocks

   do jn1 = 1, nat3
   do jn2 = 1, nat3
   
   FC_sym(jn1,jn2, i_block)=    FC(jn1,jn2, i_block)    &
                              + FC(jn2,jn1, P(i_block))                                 

   end do
   end do
                                   
   end do

   FC_sym=FC_sym/2.0_dp


   FCvar=SUM(ABS(FC-FC_sym))/ SUM(ABS(FC)) 

   if (verbose) then
   write(*, * ) ""
   write(*,  "(' FC variation due to permut. symm.= 'e20.6)") FCvar
   write(*, * ) ""
   end if
   !
end subroutine impose_perm_sym

!================================ ASR ON 2nd INDEX ============================================

subroutine impose_ASR_2nd(FC,pow,SClat,PBC,verbose, &
                         FCvar,sum2nd,FC_asr,nat,n_blocks)
   use para
   implicit none
   !integer, parameter :: DP = selected_real_kind(14,200)
   complex(kind=DP), intent(in) :: FC(3*nat,3*nat, n_blocks)
   real(kind=DP), intent(in)    :: pow
   integer, intent(in) :: nat,n_blocks
   integer, intent(in) :: SCLat(3)
   logical, intent(in) :: PBC,verbose
   !
   real(kind=DP), intent(out) :: FCvar,sum2nd
   complex(kind=DP), intent(out) :: FC_asr(3*nat,3*nat, n_blocks)
   !
   integer :: nat3,jn1,j2,n2,i_block
   complex(kind=DP) :: d1,num,den,ratio,invpow

   !    
   nat3=3*nat
   
   FC_asr=0.0_dp
   
   d1=0.0_dp
   sum2nd=0.0_dp
   !
   do jn1=1,nat3
   do j2 =1,3
           !
           num=0.0_dp
           den=0.0_dp
           !
           if ( pow > 0.01_dp ) then
               !
               !
               do i_block=1,n_blocks
               do n2=1,nat
                   num=num+ FC(jn1 , j2+3*(n2-1), i_block)
                   den=den+ ABS(FC(jn1 , j2+3*(n2-1), i_block))**pow     
               end do
               end do
               !
               if ( real(den) > 0.0_dp ) then
                   ratio=num/den
                   !
                   do i_block=1,n_blocks
                   do n2=1,nat
                       FC_asr(jn1, j2+3*(n2-1), i_block)= &
                           FC(jn1, j2+3*(n2-1), i_block)- &
                           ratio*ABS(FC(jn1, j2+3*(n2-1), i_block))**pow
                   end do
                   end do
                   !  
               else ! no need to modify the FC
                   !
                   do i_block=1,n_blocks
                   do n2=1,nat
                       FC_asr(jn1, j2+3*(n2-1), i_block)= &
                           FC(jn1, j2+3*(n2-1), i_block)
                   end do
                   end do                
                   !           
               end if
               !
               !
           else
               !
               !
               do i_block=1,n_blocks
               do n2=1,nat
                   num=num+ FC(jn1 , j2+3*(n2-1), i_block)
                   den=den+1      
               end do
               end do            
               !
               ratio=num/den
               !
               do i_block=1,n_blocks
               do n2=1,nat
                   FC_asr(jn1, j2+3*(n2-1), i_block)= &
                       FC(jn1, j2+3*(n2-1), i_block)- ratio
               end do
               end do
               !  
           end if
           !            
           !
           d1 = d1 + den
           sum2nd = sum2nd + ABS(num) ! should be zero if the ASR were fulfilled
           !
           !
   end do
   end do
   !    
   FCvar=SUM(ABS(FC-FC_ASR))/ SUM(ABS(FC)) 
   sum2nd=sum2nd/SUM(ABS(FC))
   !
   if (verbose) then
   if ( pow > 0.01_dp ) then
   
       invpow=1.0_dp/pow
       
       write(*, * ) ""   
       write(*, "(' ASR imposition on 2nd index with pow= 'f5.3)") pow
       write(*, "(' Previous values: sum(|sum_2nd phi|)/sum(|phi|)=       'e20.6)" ) sum2nd
       write(*, "('                  sum(|phi|**pow)**(1/pow)/sum(|phi|)= 'e20.6)" ) d1**invpow/SUM(ABS(FC))
       write(*, "(' FC relative variation= 'e20.6)" ) FCvar    
       write(*, * ) "" 
   
   else

       write(*, * ) ""   
       write(*, "(' ASR imposition on 2nd index with pow= 0' )") 
       write(*, "(' Previous value: sum(|sum_2nd phi|)/sum(|phi|)= 'e20.6 )" ) sum2nd
       write(*, "(' FC relative variation= 'e20.6)" ) FCvar    
       write(*, * ) ""     
   
   
   end if
   end if
   !
   !
end subroutine impose_ASR_2nd


! ==================================  MAIN ==========================================


subroutine impose_ASR(FC,R2,pow,SClat,PBC,threshold,maxite,FC_out,verbose,nat,n_blocks)
   use para
   implicit none
   !integer, parameter :: DP = selected_real_kind(14,200)
   complex(kind=DP), intent(in) :: FC(3*nat,3*nat,n_blocks)
   real(kind=DP), intent(in) :: pow,threshold
   integer, intent(in)       :: maxite, nat, n_blocks
   integer, intent(in) :: R2(3,n_blocks)
   integer, intent(in) :: SCLat(3)
   logical, intent(in) :: PBC,verbose
   !
   complex(kind=DP), intent(out) :: FC_out(3*nat,3*nat,n_blocks)
   !
   complex(kind=DP)   :: FC_tmp(3*nat,3*nat,n_blocks)
   integer         :: ite, contr, iter, ios
   real(kind=DP)   :: FCvar, sum2nd
   logical :: converged

   if (verbose) then
   write(*,*) " "
   write(*, "(' Iterative ASR imposition with pow= 'f5.3)") pow
   write(*,*) " "
   end if

   call clear_all()
   call initialize_perm(R2, n_blocks,SClat,PBC)

   converged = .false.
   FC_tmp=FC

   ite=1
   if ( maxite == 0 ) then
   contr=-1
   else
   contr=1
   end if
   iter=ite*contr

   do while (iter < maxite)

      if (verbose) write(*,"(' Iter #' I5 '  ====')") ite
      if (verbose) write(*,*) ""
         call impose_ASR_2nd(FC_tmp,pow,SClat,PBC,.false.,FCvar,sum2nd,FC_out ,nat,n_blocks)
      if (verbose) write(*,"('         Sum on 2nd='  e20.6  '    Imp. ASR on 2nd:  => delta FC=' e20.6)") sum2nd,FCvar
         call impose_perm_sym(FC_out, R2, SClat,PBC,.false.,FCvar,FC_tmp,nat,n_blocks)
      if (verbose) write(*,"('                    '  20X    '    Imp. permut sym:  => delta FC=' e20.6)") FCvar
      if (verbose) write(*,*) ""

   !  check converg
   if ( sum2nd < threshold  .and. FCvar < threshold ) then
         write(*,*) " "
         write(*,"( ' * Convergence reached within threshold:' e20.6 )") threshold
         write(*,*) " "
         write(*,"( ' * Total FC relative variation:' e20.6 )") SUM(ABS(FC-FC_out))/ SUM(ABS(FC)) 
         converged = .True.
         EXIT
   end if
   
   !  check if there is STOP file
   OPEN(unit=100, file="STOP", status='OLD', iostat=ios)
   IF(ios==0) THEN
         CLOSE(100,status="DELETE")
         write(*,*) " File STOP found, the ASR execution terminates here"
         EXIT 
   ENDIF
   
   ite=ite+1
   iter=ite*contr
   end do


   if (.not. converged ) then
   write(*,*) " "
   write(*,"( ' Max number of iteration reached ('I6')' )") maxite
   write(*,"( ' Convergence not reached within threshold:' e20.6 )") threshold
   write(*,*) " "
   end if   


end subroutine impose_ASR



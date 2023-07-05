 MODULE MPI
#ifdef USE_MPI
include 'mpif.h'
#endif
integer myrank, nproc
integer ierr
integer*8 mpi_comm
integer np_loc
integer np_gb
integer partition(2)
END MODULE

    PROGRAM LagrRadBC
    
    USE IFPORT
	USE MPI

    implicit none

    ! Variables
    integer nstates, nparticles,nBparticles
    real, allocatable :: alpha(:,:),alphaWL(:,:), M(:),S(:),P(:),M0(:),S0(:),P0(:),P0_init(:)
    real, allocatable :: X(:),Y(:),Z(:) !particle coordinates
    real, allocatable :: U(:),V(:),W(:) !Particle states velocities
    real, allocatable :: Sconc(:),Kd_Var(:),bfrac(:)
    integer, allocatable :: Head(:),LL(:),Head_free(:),LL_free(:) 
    real K0, Kh, Kd0,Kdiff_Bott,Kpart
    real a12,Kd, Abds,Ars,rho,epsilon
    real DZ, Zstar
    real CB0
    real WlFlag
    
    real t0, tmax,dt,timeC
    integer nt,it
    integer ip,lab,lab0
    integer, allocatable :: Label(:)
    real r
    real*4 r4
    real, allocatable :: prob(:)
    integer ist,jst
    real, allocatable :: MP(:),SP(:)
    integer, allocatable :: NPART(:)
    integer InitState
    real dt_out,time_out
    integer i_out
    character*5 str
    character*40 filename
    real zmin,zmax,dx
    real mpart
    integer ix,nz
    real, allocatable :: Conc(:,:), Conc_GB(:,:)
    real p12
    real xx,mm1,mm2,ss,cc,cc1,cc2
    real, allocatable :: dep_rate(:), dep_time(:)
    integer npart_cur,i
    real lambda,KiSi
    real KDIFF, DKDIFF,DK0
    integer n_aver
    real p00,Ps0,pp0,pp1,pp2,pp3
    integer i1,iproc
    integer NtotConc
    integer, allocatable :: new_seed(:)
    integer iprev, ipnext,i0_ist, ilab
    integer Nbottom, n_Bottom_gb,n_new,n_new_water

    
    
#ifdef USE_MPI
     call mpi_init(ierr)
     ! Get number of processors
     call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
     call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
     mpi_comm = MPI_COMM_WORLD
     
     r = 1000.*myrank
     CALL DRANSET (r)
#endif 
    
    
    open(1,file = 'LagrRadBC.inp')
    read(1,*);read(1,*);read(1,*) 
    read(1,*) t0, tmax, dt
    read(1,*) nparticles,nBparticles
    read(1,*) nstates 
    read(1,*) InitState
    
    
   np_GB = nparticles 
   N_bottom_gb = nBparticles !50*NP_GB
   NBottom = N_bottom_gb
    
#ifdef USE_MPI
       

np_loc = np_GB/(nproc)
Nbottom = N_bottom_GB/nproc
i1 = 1
do iproc = 0, nproc-1
if (iproc == myrank)  partition(1) = i1

if (iproc/=nproc -1)then
    if (iproc == myrank)  partition(2) = i1+ np_loc -1
    i1 = i1+np_loc
else
    if (iproc == myrank)  partition(2) = NParticles
endif    
enddo
np_loc = partition(2) - partition(1) + 1

nparticles = NP_loc
#endif	
    
    
    
    allocate (Z(nparticles),X(nparticles),Y(nparticles))
    allocate (alpha(-nstates:nstates,-nstates:nstates),alphaWL(-nstates:nstates,-nstates:nstates))
    allocate (M(0:nstates),S(0:nstates),P(-nstates:nstates),M0(0:nstates),S0(0:nstates),P0(-nstates:nstates),P0_init(-nstates:nstates))
    allocate (U(-nstates:nstates),V(-nstates:nstates),W(-nstates:nstates))
    M = 0.;S = 0.;P = 0.;P0 = 0;P0_init = 0.;M0 = 0.;S0 = 0.
    U = 0.;V = 0.;W = 0.
    allocate(label(nparticles))
    label = 0
    allocate(prob(-nstates:nstates))
    prob = 0.
    allocate (MP(-nstates:nstates),SP(-nstates:nstates),npart(-nstates:nstates))
    MP = 0.
    SP = 0.
    npart = 0
    Allocate (Bfrac(nstates))
    Bfrac = 0.
    
    allocate (Head(-nstates:nstates), LL(nparticles))
    allocate (Head_free(-nstates:-1), LL_free(nparticles))
    
    
    n_aver = 0
    
    
    Kd0 = 0.0
    Kdiff_Bott = 0.
    
   
    
    
    allocate(Kd_var(nstates),SConc(nstates))
    read(1,*) (U(ist), ist = 0, nstates)
    read(1,*) (V(ist), ist = 0, nstates)
    read(1,*) (W(ist), ist = 0, nstates)
    read(1,*) (Kd_Var(ist), ist = 1, nstates)
    read(1,*) (SConc(ist), ist = 1, nstates)
    read(1,*) K0, Kh
    read(1,*) a12
    read(1,*) abds
    read(1,*) ars
    read(1,*) (Bfrac(ist), ist = 1, nstates)
    read(1,*) rho
    read(1,*) epsilon
    read(1,*) DZ
    read(1,*) Zstar
    read(1,*) lambda
    
    KiSi = 0.
    alpha = 0.
    do ist = 1, nstates
    KiSi = KiSi + a12*kd_var(ist)*Sconc(ist)
    alpha(ist,0) = a12*kd_var(ist)*Sconc(ist)
    alpha (0,ist) = a12
    enddo
        
    alpha (0,0) = -1.*KiSi
    
    do ist = 1, nstates
    alpha(ist,ist) = -1.*a12
    enddo
    
                
    read(1,*) DT_out
    read(1,*) CB0 
    read(1,*) zmin,zmax,dx
    close(1)
!    


Mpart = CB0*rho*(1-epsilon)*Zstar/N_bottom_GB

if (nproc ==1) then
open(11,file = 'balancePart.dat')
write(11,*) 'time, Ctot, Cw, CS, CSB'
endif

!R_ab matrix for BC water layer    
Call Alpha_calc_WaterLayer(nstates, DZ,Zstar, AlphaWL, Kd_var, Bfrac, W,abds, Ars,rho, epsilon)



!Initial particle distribution  
do i = 1, nparticles
    r = DRAND(0)
    if (Initstate >=0) then
        Z(i) = zmin + r*(zmax-zmin)
    else
        Z(i) = zmax + r*Zstar
    endif
enddo
    
    
    Time_out =t0 + Dt_out
    i_out = 1
    
    nz = int((zmax-zmin)/dx)
    allocate(Conc(nz,0:nstates),Conc_GB(nz,0:nstates))
    Conc = 0.
    
   
    !Initial state
    P0(InitState) = 1.
    P0_init(InitState) = 1.
    
    nt  = (tmax - t0)/dt
    timeC = t0

    Label = InitState
    Npart = 0;
    Npart(InitState) = Nbottom
    
    
!Initial states sorting
Head = 0; Head_free = 0;
LL = 0; LL_free = 0;
do  i = 1, nparticles
    lab = label(i)
    LL(i) = Head(lab)
    Head(lab) = i
    if (lab < 0) then
        LL_free(i) = Head_Free(lab)
        Head_free(lab) = i
    endif
enddo 
     
   npart_cur = nparticles 
     
    do it = 1, nt+1
     timeC = (it -1)*dt
        
water:     do ilab = 0, nstates
          ip = Head(ilab)
          iprev = 0
         do while (ip/=0)      
         
         ipnext = LL(ip) 
         
            lab = ilab
            P0(lab)=1.
            
            Call Bcond (Z(ip),zmin,zmax,DZ,Zstar,lab,WlFlag)    
            
            
            
            
                P0 = 0.
                P0(lab) = 1.
                lab0 = lab
                call Alpha_calc(nstates, Z(ip), M, Alpha, Kd_var, Sconc, a12)
                call SolveKolm (nstates,Alpha+WlFlag*AlphaWL,P0,dt,P)
                
                DK0 = 0. 
                 
                r = DRAND(0)
                prob(-nstates) = P(-nstates)
                do ist = -nstates +1, nstates
                    prob(ist) = prob(ist-1) + P(ist)
                enddo
                do ist =-nstates, nstates
                    if (r < prob(ist)) then
                        if (lab /= ist) then
                        
                           if (iprev == 0) then
                               Head(lab) = ipnext
                                                          
                           else
                             LL(iprev) = ipnext
                           endif
                           i0_ist = Head(ist)
                           Head(ist) = ip
                           LL(ip) = i0_ist
                           
                           Npart(lab) = npart(lab) -1
                           Npart(ist) = npart(ist) + 1
                           
                           if (ist == -1) Nbottom = nBottom + 1
                           
                           
                           
                        else
                        iprev = ip   
                        endif    
                    lab = ist
                    label(ip)=lab
                    exit
                    endif
                enddo
                
                If (lab0 >=0 .and. lab <0) then !went from water to bed
                    call random (r4);r = r4
                    !r = DRAND(0)
                    Z(ip) = zmax + r*Zstar
                elseif (lab0 < 0 .and. lab>=0) then !went from bed to water
                    call random (r4);r = r4
                    !r = DRAND(0)
                    Z(ip) = zmax - r*DZ
                endif
                
                
                
                call random(r4);r = r4
                r = 2.*r - 1.
                Kpart = K0
                If (lab <0) then
                    Dk0 = 0.
                    Kpart = Kdiff_Bott
                endif
                !Vartical motion
                Z(ip) = Z(ip) + (W(lab)+DK0)*dt + r*sqrt(6*Kpart*dt) 
                
                !Horizontal motion
                call random(r4);r = r4
                r = 2.*r - 1.
                x(ip) = x(ip) + U(lab)*dt + r*sqrt(6*Kh*dt)
                call random(r4);r = r4
                r = 2.*r - 1.
                y(ip) = y(ip) + V(lab)*dt + r*sqrt(6*Kh*dt)
                 
                
           
            
        
        ip = ipnext    
            
        enddo
        enddo water !lab
        
        !Bottom particles
        P0 = 0.
        lab = -1 ! for one-fraction sediments
        P0(lab) = 1.
        lab0 = lab
        WlFlag = 1.
        call Alpha_calc(nstates, zmax, M, Alpha, Kd_var, Sconc, a12)
        call SolveKolm (nstates,Alpha+WlFlag*AlphaWL,P0,dt,P)
        
        N_new_water = int(Npart(-1)*P(0))
        call random(r4) 
        if (Npart(-1)*P(0) > r4) N_new_water = N_new_water + 1
        
        n_new = 0
        
        ip = Head(lab) !Head_free(lab)
        iprev = 0
        do while (ip /= 0 .and. n_new < N_new_water)
        ipnext = LL(ip) !LL_free(ip) 
           if (iprev == 0) then
               !Head_free(lab) = ipnext
               Head(lab) = ipnext
           else
               !LL_free(iprev) = ipnext
               LL(iprev) = ipnext
           endif
           
           Npart(lab) = Npart(lab) -1
           Npart(0) = Npart(0) + 1    
           i0_ist = Head(0)
           Head(0) = ip
           LL(ip) = i0_ist
           n_new = N_new + 1
           label(ip) = 0
           call random (r4);r = r4
           Z(ip) = zmax - r*DZ
           ip = ipnext
           if (ip==0) then
            continue
           endif
        enddo
           
  
        
        
        timeC = timeC + dt
        

        
#ifdef USE_MPI
    if (myrank ==0) then
#endif    
     Print'(a7,f10.0,a4,f8.2,a5,a3,I8)', 'time = ' , real(timeC) , ' s, ' , real(timeC/86400.) , ' days ', 'n= ', npart_cur
     !write(*,'(a7,f10.0,a4,f8.2,a3,I8)') 'time = ', real(timeC), ' s, ', real(timeC/86400.) , ' days ', 'n= ', npart_cur
                
#ifdef USE_MPI
    endif
#endif


If (nproc == 1) then
    write(11,'(f20.5,4I)') timeC/86400., sum(Npart), Npart(0), Npart(1), Npart(-1)
endif
        
        if (abs(timeC - time_out) < 10.*86400.) then
            do ip = 1, nparticles
                ix = int((Z(ip)-zmin)/dx + 1.)
                lab = Label(ip)
                if (lab>=0) then
                    if (ix >0 .and. ix <= nz) Conc(ix,lab) = Conc(ix,lab) + 1
                endif
            enddo
            n_aver = n_aver + 1
        endif
            
            
        
        if (timeC >= time_out + 5.*86400. + 100.) then
#ifdef USE_MPI
         if (myrank ==0) then
#endif                       
            write(str,'(I5)') i_out
            i_out = i_out + 1
            filename = 'out\prof_'//trim(adjustl(str))//'.dat'
            open(2,file = filename)
            write(2,*) 'x, C0, C1, C2', ' time = ', time_out
           
#ifdef USE_MPI
         endif
#endif   
         
            time_out = time_out  + DT_out
            Conc = CONC /real(n_aver)
            n_aver = 0
            
#ifdef USE_MPI
NtotConc = nx*(nstates+1)
call MPI_allreduce(Conc,Conc_gb,NtotConc,MPI_DOUBLE,MPI_SUM,mpi_comm,ierr)
#endif
            
            
   
            Conc = 1.*Conc*mpart/dx ! / real(nparticles)
#ifdef USE_MPI            
            Conc_gb = Conc_gb*mpart/dx
#endif            
            !call Conc_gauss_sol(time,-alpha(1,1),u(2),K0,ss,mm1,mm2)
            
#ifdef USE_MPI
         if (myrank ==0) then
         Conc = Conc_GB
#endif            
            do ix = 2, nz
            xx=zmin+(ix-0.5)*dx
                
            write(2,'(4f20.8)') xx, (Conc(ix,ist), ist = 0, nstates)
            enddo
            close(2)
            
#ifdef USE_MPI
         endif
#endif            
            
           
            
            
        endif
        
    enddo
    
    
    
End program LagrRadBC
    
    
    
 

subroutine SolveKolm (nstates,A,P0,dt,P)
implicit none

integer nstates
real A(-nstates:nstates,-nstates:nstates), P0(-nstates:nstates),P(-nstates:nstates)
real dt

!llocal 
real AP(-nstates:nstates)
integer i,j

do i = -nstates, nstates
Ap(i) = 0.
    do j = -nstates, nstates
        AP(i) = AP(i) + A(i,j)*P0(j)
    enddo
enddo

P = P0 + dt*AP 
End subroutine
    
Subroutine Alpha_calc(nstates, Z,MM, Alpha, Kd_var, Sconc, a12)
                   !(nstates, Z+MM, Alpha, Kd_var, Sconc, a12)
implicit none

!input
integer nstates
real z
real alpha(-nstates:nstates,-nstates:nstates)
real Kd_var(nstates),Sconc(nstates)
real MM(0:nstates)
real a12
!local 
integer ist
real conc_var
real KiSi
real C_var

    KiSi = 0.
    alpha = 0.
    do ist = 1, nstates
     if (ist ==1) then
        C_var = Sconc(ist) !*conc_var(z+MM(ist),ist)
     else
        C_var = Sconc(ist)
     endif
    KiSi = KiSi + a12*kd_var(ist)*C_var
    alpha(ist,0) = a12*kd_var(ist)*C_var
    alpha (0,ist) = a12
    enddo
    
    alpha (0,0) = -1.*KiSi
    
    do ist = 1, nstates
    alpha(ist,ist) = -1.*a12
    enddo

    end subroutine
    
    
    Subroutine Alpha_calc_Decay(nstates, Z,MM, Alpha, Kd_var, Sconc, lambda)
                   !(nstates, Z+MM, Alpha, Kd_var, Sconc, a12)
implicit none

!input
integer nstates
real z
real alpha(-nstates:nstates,-nstates:nstates)
real Kd_var(nstates),Sconc(nstates)
real MM(0:nstates)
real lambda
!local 
integer ist
real conc_var
real KiSi
real C_var

    
    do ist =-nstates, nstates
        
      alpha(ist,ist) = lambda
    
    enddo
    
    

end subroutine

Subroutine Alpha_calc_WaterLayer(nstates, DZ,Zstar, Alpha, Kd_var, Bfrac, W,abds, Ars,rho, epsilon)
                   !(nstates, Z+MM, Alpha, Kd_var, Sconc, a12)
implicit none

!input
integer nstates
real DZ, Zstar
real alpha(-nstates:nstates,-nstates:nstates)
real Kd_var(nstates),Bfrac(nstates),W(-nstates:nstates)
real Abds, Ars,rho, epsilon

!local 
integer ist,jst
real KiSi


    KiSi = 0.
    alpha = 0.
    do ist = 1, nstates
       KiSi = KiSi + kd_var(ist)*Bfrac(ist)*rho
    enddo
    
    !Diffusion
    alpha (0,0) = -1.*Abds*(Zstar/DZ)*(1-epsilon)*KiSi
    do ist = 1,nstates
        alpha(0,-ist) = Abds
        alpha(-ist,0) = Abds*(Zstar/DZ)*(1-epsilon)*Bfrac(ist)*rho*Kd_Var(ist)
        alpha(-ist,-ist) = -Abds
    enddo
    
    !Deposition
    do ist = 1, nstates
       alpha(ist,ist) = -W(ist)/DZ
       alpha(-ist,ist) = W(ist)/DZ
    enddo
    
    !Bottom fraction transfer
    do ist = -nstates, -1
        alpha(ist,ist) = alpha(ist,ist)-Ars
    enddo
    
    do ist = 1,nstates !, -1
        do jst = 1, nstates !, -1
            alpha(-ist,-jst) = alpha(-ist,-jst) + Ars*rho*Bfrac(ist)*Kd_Var(ist)/KiSi
         enddo
     Enddo
          
    
      
    
    
end subroutine

Subroutine Bcond (x,xmin,xmax,DZ,Zstar,lab,WlFlag)

implicit none


real x, xmin, xmax
real DZ, Zstar
integer lab
real WlFlag

!local
real xx

if (lab >=0) then !in water
If (x > xmax) x = 2*xmax - x
If (x < xmin) x = 2*xmin - x

else !in bottom

    if (x < xmax) x = 2*xmax - x
    if (x > xmax + Zstar) x = 2*(xmax + zstar) - x
endif

xx = x - xmax
WlFlag = 0.
IF ( (xx + DZ)*(xx - Zstar) < 0.) then
 WlFlag = 1.
endif


END SUBroutine
    
   

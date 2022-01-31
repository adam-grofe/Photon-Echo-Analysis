program twoDir
use omp_lib
implicit none
! Author: Adam Grofe
! Created: 4-11-2016
! Purpose: Calculate 2D photon echo signal from a frequency trajectory
! Reference: J.H.Choi, K.W.Kwak, M.Cho. Computational Infrared and Two-Dimensional
!            infrared photon echo spectroscopy of both wild-type and double mutant
!            myoglobin-CO proteins. J.Phys.Chem.B, 2013, 117, 15462-15478.
!
! NOTE: In the reference, equations 9 have a sign problem.  Functions 3-6 need
!       to be plus instead of minus.  Meaning that the minus sign was removed for 
!       functions 3-6.

integer :: i,j,k,l,m,logfile,lspfile,minf,maxf,t1,t3,f1,f3,t2,cnt,t1_max,&
           t3_max,nthreads,npts,dphfile,rspfile,skip,nfile,nfreq,nconfig_max,nconfig_min
integer, dimension(:), allocatable :: nconfig

complex(8) :: im,phase,exp1,exp2
complex(8), dimension(:), allocatable :: phi
complex(8), dimension(:,:), allocatable :: psi_a, psi_b,phi_odd,phi_even,signal

real(8) :: dt,Tw,bw,min_freq,max_freq,df,c,pi,gamma_1,gamma_2,real_part,imag_part,&
           integ_1,integ_2,relax,norm,weight,max_sig,t_old, t_new
real(8), dimension(:,:), allocatable :: avg_freq,smoothed,fluc
real(8), dimension(:,:,:), allocatable :: freq
real(8), dimension(2) :: ts_dip

character(len=50) :: prmfile,logfile_name,lspfile_name,rspfile_name,dphfile_name
character(len=50),dimension(:),allocatable :: datfile

! READ THE PARAMETER FILE
call cpu_time(t_old)
call get_command_argument(number=1,value=prmfile)
open(30,file=prmfile,form='formatted',action='read')
read(30,*) nfreq
read(30,*) npts
read(30,*) dt
read(30,*) gamma_1 !First excited state lifetime
read(30,*) gamma_2 !Second excited state lifetime
read(30,*) Tw
read(30,*) bw
read(30,*) ts_dip(1)
read(30,*) ts_dip(2)
read(30,*) min_freq
read(30,*) max_freq
read(30,*) skip
read(30,*) nfile
allocate( datfile(nfile) )
allocate( nconfig(nfile) )
do i=1,nfile 
   read(30,*) datfile(i), nconfig(i)
enddo
close(30)
nconfig_max = maxval( nconfig,nfile)
nconfig_min = minval( nconfig,nfile)
if( npts .gt. nconfig_min ) then
    write(*,*) "ERROR: Number of Time Points is less than Number of Configurations"
    call exit(1)
endif

gamma_1 = 1.00D0 / gamma_1
gamma_2 = 1.00D0 / gamma_2

! MEMORY ALLOCATION
allocate( freq(nconfig_max,nfreq,nfile) )
allocate( fluc(nconfig_max,nfile) )
allocate( psi_a(npts,npts) )
allocate( psi_b(npts,npts) )
allocate( phi(6) )
allocate( phi_odd(npts,npts) )
allocate( phi_even(npts,npts) )
allocate( signal(npts,npts) )
allocate( smoothed(npts,npts) )

allocate( avg_freq(nfreq,nfile) )

! READ THE DATA FILE
do k=1,nfile
   open(30,file=datfile(k),form='formatted',action='read')
   do i=1, nconfig(k)
      read(30,*) (freq(i,j,k), j=1,nfreq)
   enddo
enddo

! WRITE PARAMETERS TO LOGFILE
logfile=54
logfile_name = adjustl(trim(prmfile) // '.log')
lspfile=55
lspfile_name = adjustl(trim(prmfile) // '.lsp')
rspfile=56
rspfile_name = adjustl(trim(prmfile) // '.rsp')
dphfile=55
dphfile_name = adjustl(trim(prmfile) // '.dph')
open(logfile, file=logfile_name, form='formatted',action='write')
   write(logfile,*) "Data files and the Number of Configurations per File: ", nconfig
do i=1,nfile
   write(logfile,*) datfile(i),nconfig(i)
enddo
write(logfile,*) "Number of frequencies per line:", nfreq
write(logfile,*) "Number of Points to evaluate:",npts
write(logfile,*) "Timestep(fs):", dt
write(logfile,*) "First Excited State Relaxation Rate:", gamma_1
write(logfile,*) "Second Excited State Relaxation Rate:", gamma_2
write(logfile,*) "Waiting Time:", Tw
write(logfile,*) "Band width for smoothing:",bw
write(logfile,*) "Transition Dipoles:", ts_dip(1), ts_dip(2)
write(logfile,*) "Minimum frequency in Fourier Transform:", min_freq
write(logfile,*) "Maximum Frequency in Fourier Transform:", max_freq
! START OPENMP
!$OMP parallel
!$ nthreads = omp_get_num_threads()
!$OMP master
!$ write(logfile, '(A17,I4,A8)') 'OpenMP is using', nthreads, 'threads'
!$OMP end master
!$OMP end parallel

! ALCULATE AVERAGE TRANSITION FREQ
avg_freq = 0.00D0
!$OMP parallel do &
!$OMP& default(shared) private(i,j,k) &
!$OMP& reduction(+:avg_freq)
do k=1,nfile
do i=1,nconfig(k)
   do j=1,nfreq
      avg_freq(j,k) = avg_freq(j,k) + freq(i,j,k)
   enddo
   enddo
enddo
!$OMP end parallel do
do i=1,nfreq
do j=1,nfile
   avg_freq(i,j) = avg_freq(i,j) / nconfig(j)
enddo
enddo

write(logfile,*)  "Average Frequency:"
do j=1,nfile
   write(logfile,*) (avg_freq(i,j), i=1,nfreq)
enddo
  
call cpu_time(t_new)
write(logfile,*)
write(logfile,'(A,F15.4,F15.4)') "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new

! CONVERT FREQ TO ANGULAR FREQ AND COMPUTE FOURIER CONSTANTS
c = 1.00D15 / 2.997924858D10
pi = acos(-1.00D0)
df = ( 1.00D0 / (npts*dt))*c
minf = floor( min_freq / df )
maxf = ceiling( max_freq / df )
im = cmplx(0.00D0,1.00D0)

if( minf .ge. npts ) then
    write(*,*) 'ERROR: Minimum Frequency is not in the allowed range'
    write(*,*) 'ERROR: Min Freq:', min_freq
    write(*,*) 'ERROR: Allowed Range:', df, npts*df
    write(logfile,*) 'ERROR: Minimum Frequency is not in the allowed range'
    write(logfile,*) 'ERROR: Min Freq:', min_freq
    write(logfile,*) 'ERROR: Allowed Range:', df, npts*df
    call exit(1)
endif
if( maxf .ge. npts ) maxf = npts
if( minf .eq. 0 ) minf = 1

write(logfile,*) "Minimum Frequency Integer:", minf, "Real:", minf*df
write(logfile,*) "Maximum Frequency Integer:", maxf, "Real:", maxf*df
call cpu_time(t_new)
write(logfile,*)
write(logfile,'(A,F15.4,F15.4)') "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new

freq = freq / c * 2.00D0 * pi ! FREQ = radians/femtosecon
avg_freq = avg_freq / c * 2.00D0 * pi
do j=1,nfile
do i=1,nconfig(j)
   fluc(i,j) = freq(i,1,j) - avg_freq(1,j)
enddo
enddo

! CALCULATE THE DEPHASING INDUCED LINEBROADENING TERMS
t2 = floor( TW / dt )  !NOTE: TW is the real value in fs and t2 is the integer equivalent
t1_max = 0
t3_max = 0
psi_a = cmplx(0.00D0,0.00D0)
psi_b = cmplx(0.00D0,0.00D0)
if( npts*dt+TW+npts*dt .gt. nconfig_min*dt) then
    write(*,*) "ERROR: Number of evaluation points exceeds allowed maximum"
    write(*,*) "ERROR: t1+t2+t3:",npts*dt+TW+npts*dt
    write(*,*) "ERROR: Time max:",nconfig_min*dt 
    call exit(1)
endif
fluc = fluc*dt

!$OMP parallel do default(shared) &
!$OMP& private(t1,t3,i,j,k,cnt,integ_1,integ_2)
do t3=1, npts
do t1=1, npts
   do k=1,nfile
       cnt = 0
       do i=1, nconfig(k)-t1-t2-t3,skip   !Average over Configurations
           integ_1 = 0.00D0
           do j=0,t1-1            !Perform First Integral
              integ_1 = integ_1 + fluc(i+j,k)!*dt
           enddo
           integ_2 = 0.00D0
           do j=t1+t2,t1+t2+t3    !Perform Second Integral
               integ_2 = integ_2 + fluc(i+j,k)!*dt
           enddo
           exp1 = exp( im*integ_1)
           exp2 = exp(-im*integ_2)
           psi_a(t1,t3) = psi_a(t1,t3) + exp1*exp2
           exp1 = conjg(exp1)
           psi_b(t1,t3) = psi_b(t1,t3) + exp1*exp2

           cnt = cnt + 1
       enddo

       if(cnt .gt. 1) psi_a(t1,t3) = cmplx( real(psi_a(t1,t3))/dble(cnt), imag(psi_a(t1,t3))/dble(cnt) )
       if(cnt .gt. 1) psi_b(t1,t3) = cmplx( real(psi_b(t1,t3))/dble(cnt), imag(psi_b(t1,t3))/dble(cnt) )
    enddo
    psi_a(t1,t3) = cmplx( real(psi_a(t1,t3))/dble(nfile), imag(psi_a(t1,t3))/dble(nfile) )
    psi_b(t1,t3) = cmplx( real(psi_b(t1,t3))/dble(nfile), imag(psi_b(t1,t3))/dble(nfile) )
enddo
enddo
!$OMP end parallel do
call cpu_time(t_new)
write(logfile,*) "Dephasing function complete"
write(logfile,*)
write(logfile,'(A,F15.4,F15.4)') "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new

! WRITE OUT THE PURE DEPHASING FUNCTION
!open(dphfile,file=dphfile_name, form='formatted', action='write')
!write(dphfile,'(2(A12),4(A15))') "T1", "T3", "Real Psi A", "Imag Psi A", "Real Psi B","Imag Psi B"
!do t1=1, npts 
!   do t2=1, npts
!      write(dphfile,'(2(F12.5),4(E15.5))') t1*dt, t2*dt,real(psi_a(t1,t2)),imag(psi_a(t1,t2)),&
!             real(psi_b(t1,t2)),imag(psi_b(t1,t2))
!   enddo
!   write(dphfile,*)
!enddo
!close(dphfile)

! CALCULATE RESPONSE FUNCTION
phi = cmplx(0.00D0,0.00D0)
phi_odd = cmplx(0.00D0,0.00D0)
phi_even = cmplx(0.00D0,0.00D0)
do i=1,nfreq
   do j=2,nfile
      avg_freq(i,1) = avg_freq(i,1) + avg_freq(i,j)
   enddo
   avg_freq(i,1) = avg_freq(i,1) / dble(nfile)
enddo
!$OMP parallel do default(shared) &
!$OMP& private(t1,t3,phase,relax,phi)
do t1=1, npts
   do t3=1, npts
      phase = exp( -im*cmplx(avg_freq(2,1)*dble(t3)*dt,0.00D0) + im*cmplx(avg_freq(1,1)*dble(t1)*dt,0.00D0) )
      relax = exp( -(gamma_1+gamma_2)/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(1) = cmplx(-ts_dip(1)**2*ts_dip(2)**2,0.00D0) * phase * psi_a(t1,t3) * cmplx(relax,0.00D0)

      phase = exp( -im*cmplx(avg_freq(2,1)*t3*dt,0.00D0) - im*cmplx(avg_freq(1,1)*t1*dt,0.00D0) )
      relax = exp( -(gamma_1+gamma_2)/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(2) = cmplx(-ts_dip(1)**2*ts_dip(2)**2,0.00D0) * phase * psi_b(t1,t3) * cmplx(relax,0.00D0)

      phase = exp( -im*cmplx(avg_freq(1,1)*t3*dt,0.00D0) + im*cmplx(avg_freq(1,1)*t1*dt,0.00D0) )
      relax = exp( -gamma_1/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(3) = cmplx( ts_dip(1)**4,0.00D0) * phase * psi_a(t1,t3) * cmplx(relax,0.00D0)

      phase = exp( -im*cmplx(avg_freq(1,1)*t3*dt,0.00D0) - im*cmplx(avg_freq(1,1)*t1*dt,0.00D0) )
      relax = exp( -gamma_1/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(4) = cmplx( ts_dip(1)**4,0.00D0) * phase * psi_b(t1,t3) * cmplx(relax,0.00D0)

      phase = exp( -im*cmplx(avg_freq(1,1)*t3*dt,0.00D0) + im*cmplx(avg_freq(1,1)*t1*dt,0.00D0) )
      relax = exp( -gamma_1/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(5) = cmplx( ts_dip(1)**4,0.00D0) * phase * psi_a(t1,t3) * cmplx(relax,0.00D0)

      phase = exp( -im*cmplx(avg_freq(1,1)*t3*dt,0.00D0) - im*cmplx(avg_freq(1,1)*t1*dt,0.00D0) )
      relax = exp( -gamma_1/2.00D0*t3*dt - gamma_1*Tw - gamma_1*t1*dt/2.00D0 )
      phi(6) = cmplx( ts_dip(1)**4,0.00D0) * phase * psi_b(t1,t3) * cmplx(relax,0.00D0)

      phi_odd(t1,t3)  = phi(1) + phi(3)  + phi(5)
      if( phi_odd(t1,t3) /= phi_odd(t1,t3) .or. isnan(real(phi_odd(t1,t3))) &
           .or. isnan(imag(phi_odd(t1,t3))) ) then
          write(*,*) "Found NAN in Odd Total Linebroadening Function:",t1,t3
          write(logfile,*) "Found NAN in Odd Total Linebroadening Function:",t1,t3
          call exit(1)
      endif

      phi_even(t1,t3) = phi(2) +  phi(4)  + phi(6) 
      if( phi_even(t1,t3) /= phi_even(t1,t3) .or. isnan(real(phi_even(t1,t3))) &
          .or. isnan(imag(phi_even(t1,t3)))  ) then
          write(*,*) "Found NAN in Even Total Linebroadening Function:",t1,t3
          write(logfile,*) "Found NAN in Even Total Linebroadening Function:",t1,t3
          call exit(1)
      endif
   enddo
enddo
!$OMP end parallel do
call cpu_time(t_new)
write(logfile,*) "Response function Complete"
write(logfile,*)
write(logfile,*) "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new

! WRITE OUT THE RESPONSE FUNCTION
!open( rspfile,file=rspfile_name, form='formatted', action='write')
!write(rspfile,'(2(A12),4(A15))') "T1", "T3", "Real Rsp Odd", "Imag Rsp Odd",&
!                        "Real Rsp Even", "Imag Rsp Even"
!do t1=1, npts 
!   do t2=1, npts
!      write(rspfile,'(2(F12.5),4(E15.5))') t1*dt, t2*dt,real(phi_odd(t1,t2)),imag(phi_odd(t1,t2)),&
!             real(phi_even(t1,t2)),imag(phi_even(t1,t2))
!   enddo
!   write(rspfile,*)
!enddo
!close(rspfile)

! PERFORM FOUIER TRANSFORM
signal = 0.00D0
df = df * 2.00D0 * pi / c
!$OMP parallel do default(shared) &
!$OMP private(f1,f3,t1,t3)
do f1=minf, maxf
do f3=minf, maxf
   do t1=1, npts
   do t3=1, npts
      signal(f1,f3) = signal(f1,f3) +  &
              exp( im*cmplx(dble(f3*df*t3*dt-f1*df*t1*dt),0.00D0) )*phi_odd(t1,t3) &
            + exp( im*cmplx(dble(f3*df*t3*dt+f1*df*t1*dt),0.00D0) )*phi_even(t1,t3) 

!      if( signal(f1,f3) /= signal(f1,f3) .or. isnan( real(signal(f1,f3))) ) then
!          write(*,*) "Found NAN in Signal Function:",f1,f3,t1,t3
!          write(logfile,*) "Found NAN in Signal Function:",f1,f3,t1,t3
!          call exit(1)
!      endif
   enddo
   enddo
enddo
enddo
!$OMP end parallel do
write(logfile,*) "Fourier Transform Complete"
call cpu_time(t_new)
write(logfile,*)
write(logfile,*) "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new
df = df * c / 2.00D0 / pi

! SMOOTH THE SIGNAL 
smoothed = 0.00D0
!$OMP parallel do default(shared) &
!$OMP& private(f1,f3,i,j,norm,weight)
do f1=minf, maxf
do f3=minf, maxf
   norm = 0.00D0
   do i=minf, maxf
   do j=minf, maxf
      weight = exp( -(df*(f1-i))**2-(df*(f3-j))**2/2.00D0/bw**2)
      norm = norm + weight
      smoothed(f1,f3) = smoothed(f1,f3) + real(signal(i,j))*weight
   enddo
   enddo
   smoothed(f1,f3) = smoothed(f1,f3) / norm
enddo
enddo
!$OMP end parallel do
write(logfile,*) "Smoothing Complete"
call cpu_time(t_new)
write(logfile,*)
write(logfile,*) "Time(min,sec):", (t_new-t_old)/60.0D0,t_new-t_old
write(logfile,*)
t_old = t_new

! NORMALIZE THE SPECTRA TO THE MAXIMUM
max_sig = 0.00D0
do f1=minf, maxf
do f3=minf, maxf
    if( abs(smoothed(f1,f3)) .gt. max_sig ) max_sig = abs(smoothed(f1,f3))
enddo
enddo
smoothed = smoothed / max_sig
signal = signal / max_sig
write(logfile,*) "Normalization Complete"
write(logfile,*) "Max signal:", max_sig

! WRITE OUT THE RESULT
open(lspfile, file=lspfile_name, form='formatted',action='write')
write(lspfile,'(5(A20))') 'Freq 1(cm-1)', 'Freq 1(cm-1)', 'Smoothed', 'Real Signal', 'Imag. Signal'
do f1=minf, maxf
   do f3=minf, maxf
      write(lspfile,'(2(F20.10),3(E20.10))') df*f1,df*f3,smoothed(f1,f3),real(signal(f1,f3)), &
                   imag(signal(f1,f3))
   enddo
   write(lspfile,*)
enddo
close(lspfile)

write(logfile,*) "Output written to:", lspfile_name
write(logfile,*) "NORMAL TERMINATION"



end program twoDir

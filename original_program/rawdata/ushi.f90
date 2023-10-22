MODULE nrtype
!    Symbolic names for kind types of 4-, 2-,and 1-byte integers:
      INTEGER, PARAMETER :: I4B=SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B=SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: I1B=SELECTED_INT_KIND(2)
!    SymbolicnarTlesfor kind types of single- and double-precision reals:
      INTEGER, PARAMETER :: SP=KIND(1.0)
      INTEGER, PARAMETER :: DP=KIND(1.0D0)
!    Symbolic names for kind types of single- and double-precision complex:
      INTEGER, PARAMETER :: SPC=KIND((1.0,1.0))
      INTEGER, PARAMETER :: DPC=KIND((1.0D0,1.0D0))
!    Symbolic name for kind type of default logical:
      INTEGER, PARAMETER ::LGT=KIND(.true.)
!    Frequently used mathematical constants (with precision to spare):
      REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
      REAL(SP), PARAMETER :: PI02=1.57079632679489661923132169163975144209858_sp
      REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
      REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
      REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
      REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
      REAL(DP), PARAMETER :: PI02_D=1.57079632679489661923132169163975144209858_dp
      REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
      REAL(DP), PARAMETER :: Grv=9.81_dp
!    Derived data types for sparse matrices, single and double precision:
      TYPE sprs2_sp
          INTEGER(I4B) :: n,len
          REAL(SP), DIMENSION(:), POINTER :: val
          INTEGER(I4B), DIMENSION(:), POINTER :: irow
          INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_sp
      TYPE sprs2_dP
          INTEGER(I4B) :: n,len   
          REAL(DP), DIMENSION(:), POINTER :: val
          INTEGER(I4B), DIMENSION(:), POINTER :: irow
          INTEGER(I4B), DIMENSION(:), POINTER :: jcol
      END TYPE sprs2_dp
END MODULE nrtype

MODULE variables
      USE nrtype
      INTEGER(I4B), PARAMETER               :: ndim=500, ndim2=6000
      INTEGER(I4B)                                       :: nmin, nmax
      INTEGER(I4B)                                       :: nelm1=30, nvtx1=31, counter1
      REAL(DP), DIMENSION(3,0:3,ndim)     :: X1, X10, X11, X00, X1dt, X1dt0, X1dt2, X1dt20, X1dt1, X1dt21, Flow1dt
      REAL(DP), DIMENSION(ndim,ndim)     :: Xs, Xsdt, Xsdt0, FXs
      REAL(DP), DIMENSION(24,24,ndim)    :: Mass1, Kl1, Kb1, Kall1
      REAL(DP), DIMENSION(3,0:1,4)            ::  rungeKob, rungeKab
      REAL(DP), DIMENSION(ndim2,ndim2)    :: M_st1, K_st1
      REAL(DP), DIMENSION(ndim2)          :: Fout1, Fext1
      REAL(DP), DIMENSION(ndim)            :: top_x, top_y
      REAL(DP), DIMENSION(ndim)            :: u1, u1x, u1y, u1z, ux1, uy1, dx1, dy1, dz1
      REAL(DP), DIMENSION(ndim)            :: Fd1, Fdx1, Fdy1, Fdz1, Fdx10, Fdy10, Fdz10
      REAL(DP), DIMENSION(ndim2)          :: Qe1, lk1, Me1, Qkl1, Qkb1
      REAL(DP), DIMENSION(ndim)     ::  Iyy1,  l_elm1, l_elm10, rho1, m_elm1, theta1, phi1, r1, Ap1
      REAL(DP), DIMENSION(ndim)     ::  d_out1, d_in1, A1
      REAL(DP), DIMENSION(ndim)     ::  vel, deg, dep_v, dep_d
      REAL(DP), DIMENSION(ndim)     ::  angleXs, angleXs0, cable1dt, cable1dt0, cable1dt2
      REAL(DP), DIMENSION(ndim)     ::  XabOA, XabAO, XobOA, XobAO
      REAL(DP)                                      :: E1=1.20d9, length1=100.d0, Cd=1.d0, tolerance1=1.d-7
      REAL(DP)                                      :: dt=1.d-4, t, emax1, alp=0.5d0, rhow=1030.d0, div
      REAL(DP)                                      :: Sx=0.2675d0, Sy=0.2675d0, Sz=0.19634954d0
      REAL(DP)                                      :: Fsdx_t, Fsdy_t, Fsdz_t
      REAL(DP)                                      :: FextXs_t, FextYs_t, FextZs_t
      REAL(DP)                                      :: xi, S1, S2, S3, S4
      REAL(DP)                                      :: down, rmd1
!======       Sinker      ============================================================================
      REAL(DP)                            :: mRs=30.0d0, volRs=0.0195122d0, T1x, T1y, T1z, Cds=-1.d0
      REAL(DP)                            :: CdxS, CdyS, CdzS      
      REAL(DP)                            :: PPAP
!============================================================================================================================================
!======  parametric_spline ============================================
      INTEGER(I4B), PARAMETER             :: ndim3=10000
      INTEGER(I4B)                                      :: nelm_spln, nvtx_spln, num
      REAL(DP), DIMENSION(0:ndim3)        ::  t_spln, tn_spln
      REAL(DP), DIMENSION(0:ndim3)        ::  x, x_n, Sx_spln, spln_a_x, spln_b_x, spln_c_x, spln_d_x, kn_x
      REAL(DP), DIMENSION(0:ndim3)        ::  y, y_n, Sy_spln, spln_a_y, spln_b_y, spln_c_y, spln_d_y, kn_y
!===============================================================
END MODULE variables

PROGRAM ROV
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, l, n, f
      REAL(DP)                            :: FZ1,FZ2, Error
                
      call Filesetting
      
      call InitialSettingsCable1
      call InitialSettingsSinker

     n=0
     do
       read(10,*, end=10)
       n=n+1
     end do
10  rewind(10)

     do i=1,n
        read(10,*) top_x(i), top_y(i)
     end do

     do i=1,n
        write(*,*) i, top_x(i), top_y(i)
     end do

!===========degree================= 
     write(*,*) 'degree'

     num=0
     do
       read(1,*, end=20)
       num=num+1
     end do
20  rewind(1)

     do i=1,num
       read(1,*) x(i-1),y(i-1)
       write(*,*) x(i-1),y(i-1)
     end do  
        
     call Parametric_SplineInterpolate

     do i=1,nvtx_spln
       if(0.d0<=x_n(i) .and. x_n(i)<270.d0) then
         deg(i)=(90.d0-x_n(i))/180.d0*PI
       else if (270.d0<=x_n(i) .and.  x_n(i)<=350.d0) then
         deg(i)=(90.d0-(x_n(i)-360.d0))/180.d0*PI
       else if (350.d0<x_n(i)) then
         deg(i)=(90.d0-(350.d0-360.d0))/180.d0*PI
       end if
       dep_d(i)=y_n(i)
     end do
     
     do i=1,nvtx_spln
       write(*,*) deg(i)/PI*180.d0, dep_d(i)
     end do 
!===========degree==================
!===========velocity==================       
    write(*,*) 'velocity'
    
     num=0
     do
       read(2,*, end=30)
       num=num+1
     end do
30  rewind(2)

     do i=1,num
       read(2,*) x(i-1),y(i-1)
       write(*,*) x(i-1),y(i-1)
     end do  

     call Parametric_SplineInterpolate

      do i=1,nvtx_spln
        vel(i)=x_n(i)
        dep_v(i)=y_n(i)
      end do
      pause
      
!==========velocity===================     

     call OutputCable1
     call OutputSinker     
     call Outputcurrent

     do i=1,2000001
        t=dt*real(i)
        counter1=0
        emax1=1.d0
        div=real(i)
        write(*,*) "t=",t  
          call Tidalcurrent
          call Extend_cable
          call rungeKuttaSinker
          call UpdateSinker
          call ANCFSinker
          
        PPAP=dsqrt((X1(1,0,1)-X1(1,0,nvtx1))**2+(X1(2,0,1)-X1(2,0,nvtx1))**2)
        write(*,*) "D=", PPAP
         
!=======OUTPUT=============================================
          if(mod(i,10000)==1) then  
            write(51,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',(nvtx1),' F=POINT'  
            write(61,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',10,' F=POINT'         
            call OutputCable1
            call OutputSinker 
            call Outputcurrent
          end if
!=========================================================
       end do
     
!=======OUTPUT=============================================
          write(51,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',(nvtx1),' F=POINT'  
          write(61,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',10,' F=POINT'        
                        
          call OutputCable1
          call OutputSinker 
     
          call Outputcurrent
!=========================================================

END PROGRAM ROV

SUBROUTINE Filesetting
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, l, n, f

      open(1,file='Tidalcurrent_degree.csv')
      open(2,file='Tidalcurrent_velocity.csv')
      open(10,file='Velocity_fukaemaru.csv')
      open(51,file='2K3D-Cable1.dat')
	  open(56,file='cable_simulation.dat')
      open(61,file='2K3D-Sinker.dat')
      open(71,file='Fext1.csv')
      open(73,file='T1.csv')
      open(81,file='Angle1.csv')
      open(91,file='Thrast1.csv')
      open(101,file='Me1.csv')
      open(111,file='Qkb1.csv')
      open(121,file='Qk1l.csv')
      open(131,file='Cable1dt2.csv')
      open(141,file='Cable1dt.csv')
      open(151,file='Sinkerdt.csv')
      open(161,file='Fsinker.csv')
      open(171,file='current1.csv')
      open(20,file='Fdx1.csv')
      open(21,file='Extend_cable1.csv')
      open(28,file='Fout1.csv')
      open(31,file='tansindo1.csv')
      open(77,file='output_spline_x.csv')    
      open(88,file='output_spline_y.csv')
      open(99,file='output_spline.csv')     

      write(71,*)'"t"', ',', '"Fext1(1)"', ',', '"Fext1(2)"', ',', '"Fext1(3)"', ',', '"Fext1(7)"', ',', '"Fext1(8)"', ',', '"Fext1(9)"' 
      write(73,*) 't', ',', 'T1_x', ',', 'T1_y', ',', 'T1_z'
      write(81,*) '"t"', ',', '"angleXs(1)"', ',', '"angleXs(2)"', ',', '"angleXs(3)"'
      write(101,*)'"t"', ',', '"Me1X"', ',', '"Me1Y"', ',', '"Me1Z"'
      write(111,*)'"t"', ',', '"Qkb1X"', ',', '"Qkb1Y"', ',', '"Qkb1Z"'
      write(121,*)'"t"', ',', '"Qkl1X"', ',', '"Qkl1Y"', ',', '"Qkl1Z"'
      write(131,*) 't', ',', 'X1dt2x', ',',  'X1dt2y', ',', 'X1dt2z'
      write(141,*) 't', ',', 'X1dtx', ',', 'X1dty', ',', 'X1dtz'
      write(151,*) 't', ',', 'Xsdtx', ',', 'Xsdty', ',', 'Xsdtz'
      write(161,*) 't', ',', 'FXs', ',', 'FYs', ',', 'FZs'
      write(171,*) 't', ',', 'Flow1dt_X', ',', 'Flow1dt_Y', ',', 'Flow1dt_Z'    
      write(31,*)    't', ',', 'X', ',', 'Y', ',', 'Z', ',', 'X', ',', 'Y', ',', 'Z'      
      write(21,*)    't', ',', 'down', ',', 'length1', ',', 'l_elm(i)'   
    
      write(51,*)'TITLE = "Cable1"'
      write(51,*)'VARIABLES ="X" "Y" "Z"'      
      write(51,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',(nvtx1),' F=POINT'
      
      write(61,*)'TITLE = "Sinker"'
      write(61,*)'VARIABLES ="X" "Y" "Z"'      
      write(61,'(a,f10.5,a,i4,a)') 'ZONE T="Time=',t,'[sec]" I=',10,' F=POINT'

END SUBROUTINE Filesetting


SUBROUTINE Parametric_SplineInterpolate
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, n
      REAL(DP)                            :: del, delmax, alp_spln=0.5d0, tol_spln=1.d-10

       n=num
      
      do i=0, n-1
        Sx_spln(i)=0.d0
        Sy_spln(i)=0.d0
      end do
      
      write(*,*) 'kn_x'
      do i=1,n-2
         kn_x(i)=6.d0*(x(i+1)-2*x(i)+x(i-1))
         write(*,*) i, kn_x(i)
      end do
      
      write(*,*) 'kn_y' 
      do i=1,n-2
         kn_y(i)=6.d0*(y(i+1)-2*y(i)+y(i-1))
         write(*,*) i, kn_y(i)
      end do

!---x---
  
      k=0
      delmax=1.d0
      
      do while(delmax>tol_spln)
        k=k+1
        delmax=0.d0
        do i=1, n-2
           del=(kn_x(i)-Sx_spln(i-1)-Sx_spln(i+1))/4.d0-Sx_spln(i)
           Sx_spln(i)=Sx_spln(i)+alp_spln*del
          if(abs(del) > delmax) then
            delmax=abs(del)
          end if
        end do
      end do      
      
      write(*,*) 'k=', k
      
      do i=0, n-1
         spln_d_x(i)=x(i)
         spln_b_x(i)=0.5d0*Sx_spln(i)
         spln_a_x(i)=(Sx_spln(i+1)-Sx_spln(i))/6.d0
         spln_c_x(i)=(x(i+1)-x(i))-(2.d0*Sx_spln(i)+Sx_spln(i+1))/6.d0
         
         write(*,'(i2,7(a1,f16.8))') i, ',', Sx_spln(i), ',', spln_a_x(i), ',', spln_b_x(i), ',', spln_c_x(i), ',', spln_d_x(i)  
         write(77,'(i2,7(a1,f16.8))') i, ',', Sx_spln(i), ',', spln_a_x(i), ',', spln_b_x(i), ',', spln_c_x(i), ',', spln_d_x(i)  
       end do
!--------
      
!---y---
      k=0
      delmax=1.d0
      do while(delmax>tol_spln)
        k=k+1
        delmax=0.d0
        do i=1, n-2
           del=(kn_y(i)-Sy_spln(i-1)-Sy_spln(i+1))/4.d0-Sy_spln(i)
           Sy_spln(i)=Sy_spln(i)+alp_spln*del
          if(abs(del) > delmax) then
            delmax=abs(del)
          end if
        end do
      end do    
     
      write(*,*) 'k=', k

      do i=0, n-1
        spln_d_y(i)=y(i)
        spln_b_y(i)=0.5d0*Sy_spln(i)
        spln_a_y(i)=(Sy_spln(i+1)-Sy_spln(i))/6.d0
        spln_c_y(i)=(y(i+1)-y(i))-(2.d0*Sy_spln(i)+Sy_spln(i+1))/6.d0
         
        write(*,'(i5,8(a1,f16.8))') i, ',',Sy_spln(i), ',', spln_a_y(i), ',', spln_b_y(i), ',', spln_c_y(i), ',', spln_d_y(i)   
        write(88,'(i5,7(a1,f16.8))') i, ',',Sy_spln(i), ',', spln_a_y(i), ',', spln_b_y(i), ',', spln_c_y(i), ',', spln_d_y(i)  
      end do
!-------
      nelm_spln=375
      
!~       write(*,*) 'nelm='
!~       read(*,*) nelm_spln
      
      do i=0,n-2
        tn_spln(i)=i
        do j=1,nelm_spln
          k=i*nelm_spln+j
          t_spln(j)=real(j-1)*1.d0/nelm_spln+tn_spln(i)
          x_n(k)=spln_a_x(i)*(t_spln(j)-tn_spln(i))**3+spln_b_x(i)*(t_spln(j)-tn_spln(i))**2+spln_c_x(i)*(t_spln(j)-tn_spln(i))+spln_d_x(i)
          y_n(k)=spln_a_y(i)*(t_spln(j)-tn_spln(i))**3+spln_b_y(i)*(t_spln(j)-tn_spln(i))**2+spln_c_y(i)*(t_spln(j)-tn_spln(i))+spln_d_y(i)
          write(*,'(i5,4(a1,f16.8))') k, ',', tn_spln(i), ',', t_spln(k), ',', x_n(k), ',', y_n(k)
          write(99,'(f16.8,2(a1,f16.8))') x_n(k), ',',  y_n(k)
        end do
      end do
      
      i=n-1
      k=i*nelm_spln+1
      nvtx_spln=k
           
      tn_spln(i)=i
      t_spln(j)=tn_spln(i)
      x_n(k)=spln_a_x(i)*(t_spln(j)-tn_spln(i))**3+spln_b_x(i)*(t_spln(j)-tn_spln(i))**2+spln_c_x(i)*(t_spln(j)-tn_spln(i))+spln_d_x(i)
      y_n(k)=spln_a_y(i)*(t_spln(j)-tn_spln(i))**3+spln_b_y(i)*(t_spln(j)-tn_spln(i))**2+spln_c_y(i)*(t_spln(j)-tn_spln(i))+spln_d_y(i)
      write(*,'(i5,4(a1,f16.8))') k, ',', tn_spln(k), ',', t_spln(k), ',', x_n(k), ',', y_n(k)
      write(99,'(f16.8,2(a1,f16.8))') x_n(k), ',', y_n(k)
      
END SUBROUTINE Parametric_SplineInterpolate


SUBROUTINE Tidalcurrent
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, ii, jj, n
     
     nmin=1
     nmax=nvtx_spln
     
     do i=1,nvtx1
        if (i==1) then
!~                 theta1(i)=-60.d0/180.d0*PI
                theta1(i)=deg(1)*PI/180.d0
                Flow1dt(1,0,i)=x_n(i)*cos(theta1(i))
                Flow1dt(2,0,i)=x_n(i)*sin(theta1(i))
                Flow1dt(3,0,i)=0.d0
        else
          do n=nmin,nmax-1
             if (y_n(n)>=X1(3,0,i) .and. X1(3,0,i)>=y_n(n+1)) then
!~                 theta1(i)=-60.d0/180.d0*PI
                theta1(i)=(((y_n(n+1)-y_n(n))/(y_n(n+1)-y_n(n)))*(X1(3,0,i)-y_n(n))+y_n(n))*PI/180.d0
                Flow1dt(1,0,i)=(((x_n(n+1)-x_n(n))/(y_n(n+1)-y_n(n)))*(X1(3,0,i)-y_n(n))+x_n(n))*cos(theta1(i))
                Flow1dt(2,0,i)=(((x_n(n+1)-x_n(n))/(y_n(n+1)-y_n(n)))*(X1(3,0,i)-y_n(n))+x_n(n))*sin(theta1(i))
                Flow1dt(3,0,i)=0.d0
             end if
          end do
        end if
     end do 
     
END SUBROUTINE Tidalcurrent


SUBROUTINE Outputcurrent
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, ii, jj, n
      
     write(171,*) t
     
     do i=1,nvtx1
       write(171,'(I2,a1,2(f16.8,a1),f16.8)') i, ',', Flow1dt(1,0,i), ',', Flow1dt(2,0,i), ',', Flow1dt(3,0,i)
     end do
    
END SUBROUTINE  Outputcurrent


SUBROUTINE InitialSettingsSinker
      USE nrtype
      USE variables
      IMPLICIT NONE

      Xs(1,0)=X1(1,0,nvtx1)
      Xs(2,0)=X1(2,0,nvtx1)
      Xs(3,0)=X1(3,0,nvtx1)-0.2675d0

      Xs(1,1)=X1(1,0,nvtx1)+0.250d0
      Xs(2,1)=X1(2,0,nvtx1)+0.250d0
      Xs(3,1)=X1(3,0,nvtx1)
      
      Xs(1,2)=X1(1,0,nvtx1)-0.25d0
      Xs(2,2)=X1(2,0,nvtx1)+0.25d0
      Xs(3,2)=X1(3,0,nvtx1)
      
      Xs(1,3)=X1(1,0,nvtx1)-0.250d0
      Xs(2,3)=X1(2,0,nvtx1)+0.250d0
      Xs(3,3)=X1(3,0,nvtx1)-0.535d0
      
      Xs(1,4)=X1(1,0,nvtx1)+0.250d0
      Xs(2,4)=X1(2,0,nvtx1)+0.250d0
      Xs(3,4)=X1(3,0,nvtx1)-0.535d0
      
      Xs(1,5)=X1(1,0,nvtx1)+0.250d0
      Xs(2,5)=X1(2,0,nvtx1)-0.250d0
      Xs(3,5)=X1(3,0,nvtx1)
      
      Xs(1,6)=X1(1,0,nvtx1)-0.250d0
      Xs(2,6)=X1(2,0,nvtx1)-0.250d0
      Xs(3,6)=X1(3,0,nvtx1)
      
      Xs(1,7)=X1(1,0,nvtx1)-0.250d0
      Xs(2,7)=X1(2,0,nvtx1)-0.250d0
      Xs(3,7)=X1(3,0,nvtx1)-0.535d0
      
      Xs(1,8)=X1(1,0,nvtx1)+0.250d0
      Xs(2,8)=X1(2,0,nvtx1)-0.250d0
      Xs(3,8)=X1(3,0,nvtx1)-0.535d0
      
      Xs(1,9)=X1(1,0,nvtx1)
      Xs(2,9)=X1(2,0,nvtx1)
      Xs(3,9)=X1(3,0,nvtx1)
      
      angleXs(1)=0.d0
      angleXs(2)=0.d0
      angleXs(3)=0.d0
      angleXs0(1)=0.d0
      angleXs0(2)=0.d0
      angleXs0(3)=0.d0

      Xsdt(1,0)=0.d0
      Xsdt(2,0)=0.d0
      Xsdt(3,0)=0.d0
      Xsdt(1,1)=0.d0
      Xsdt(2,1)=0.d0
      Xsdt(3,1)=0.d0
      Xsdt0(1,0)=0.d0
      Xsdt0(2,0)=0.d0
      Xsdt0(3,0)=0.d0
      Xsdt0(1,1)=0.d0
      Xsdt0(2,1)=0.d0
      Xsdt0(3,1)=0.d0
      
END SUBROUTINE InitialSettingsSinker


SUBROUTINE forthSinker
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i
      
      XabAO(1)=FextXs_t
      XabAO(2)=FextYs_t
      XabAO(3)=FextZs_t
      call Coordinate_Converter_AtoO_Sinker

      T1x=-XobAO(1)
      T1y=-XobAO(2)
      T1z=-XobAO(3)
      
      call Forceofcurrent1

      CdxS=0.5d0*rhow*Sx*Cds*(Xsdt(1,0)-Flow1dt(1,0,nvtx1))*abs(Xsdt(1,0)-Flow1dt(1,0,nvtx1))
      CdyS=0.5d0*rhow*Sy*Cds*(Xsdt(2,0)-Flow1dt(2,0,nvtx1))*abs(Xsdt(2,0)-Flow1dt(2,0,nvtx1))
      CdzS=0.d0!0.5d0*rhow*Sz*Cd*(X1dt(3,0,nvtx1)-Flow1dt(3,0,nvtx1))*abs(X1dt(3,0,nvtx1)-Flow1dt(3,0,nvtx1))
     
      FXs(1,0)=-(rhow*volRs-mRs)*Grv*dsin(angleXs(2))+T1x+Fsdx_t+CdxS
      FXs(2,0)= (rhow*volRs-mRs)*Grv*dsin(angleXs(1))*dcos(angleXs(2))+T1y+Fsdy_t+CdyS
      FXs(3,0)= (rhow*volRs-mRs)*Grv*dcos(angleXs(1))*dcos(angleXs(2))+T1z+Fsdz_t+CdzS
      FXs(1,1)=0.d0
      FXs(2,1)=0.d0
      FXs(3,1)=0.d0
     
END SUBROUTINE forthSinker


SUBROUTINE Forceofcurrent1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, ii, jj, n
	     
         i=nvtx1
         
         u1(i)=(X1dt1(1,0,i)-Flow1dt(1,0,i))**2+(X1dt1(2,0,i)-Flow1dt(2,0,i))**2
         
         if ((X1dt1(1,0,i)-Flow1dt(1,0,i))==0.d0) then
           phi1(i)=PI/2.d0
         else
           phi1(i)=atan(abs((X1dt1(2,0,i)-Flow1dt(2,0,i))/(X1dt1(1,0,i)-Flow1dt(1,0,i))))
         end if
         
         if ((X1dt1(1,0,i)-Flow1dt(1,0,i))<0.d0) then
            ux1(i)=-1
         else 
            ux1(i)=1
         end if
         
         if ((X1dt1(2,0,i)-Flow1dt(2,0,i))<0.d0) then
            uy1(i)=-1
         else 
            uy1(i)=1
         end if
            
          dx1(i)=X1(1,0,i+1)-X1(1,0,i)
          dy1(i)=X1(2,0,i+1)-X1(2,0,i)
          dz1(i)=X1(3,0,i+1)-X1(3,0,i)
          lk1(i)=dsqrt((dx1(i)**2+dy1(i)**2+dz1(i)**2))
          
          u1x(i)=-(X1dt1(1,0,i)-Flow1dt(1,0,i))
          u1y(i)=-(X1dt1(2,0,i)-Flow1dt(2,0,i))
          
!--------r1(i)=dsqrt(1-cos(theta)^2)=sin(theta)----------------
          r1(i)=dsqrt(1-(((u1x(i)*dx1(i))+(u1y(i)*dy1(i)))/(dsqrt(u1x(i)**2+u1y(i)**2)*lk1(i)))**2)
!---------------------------------------------------------------------------          

          Ap1(i)=d_out1(i)*lk1(i)*r1(i)
          Fd1(i)=-0.5d0*rhow*Cd*Ap1(i)*u1(i)

          Fdx1(i)=Fd1(i)*cos(phi1(i))*ux1(i)
          Fdy1(i)=Fd1(i)*sin(phi1(i))*uy1(i)
          Fdz1(i)=0.d0

          Fsdx_t=Fdx1(i)/2.d0
          Fsdy_t=Fdy1(i)/2.d0
          Fsdz_t=Fdz1(i)/2.d0

END SUBROUTINE Forceofcurrent1


SUBROUTINE rungeKuttaSinker
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k
      
        do k=1,4
          do j=0,0
            do i=1,3
           rungeKob(i,j,k)=0.d0
           rungeKab(i,j,k)=0.d0
            end do
          end do
        end do

        do j=0,0
          do i=1,3     
            Xsdt(i,j)=Xsdt0(i,j)
          end do
        end do
        do i=1,3
          angleXs(i)=angleXs0(i)
          XobOA(i)=Xsdt(i,0)
        end do
        call forthSinker
        call Coordinate_Converter_OtoA_Sinker
        rungeKob(1,0,1)=dt*(mRs*Xsdt(3,1)*Xsdt(2,0)-mRs*Xsdt(2,1)*Xsdt(3,0)+FXs(1,0))/mRs
        rungeKob(2,0,1)=dt*(mRs*Xsdt(1,1)*Xsdt(3,0)-mRs*Xsdt(3,1)*Xsdt(1,0)+FXs(2,0))/mRs
        rungeKob(3,0,1)=dt*(mRs*Xsdt(2,1)*Xsdt(1,0)-mRs*Xsdt(1,1)*Xsdt(2,0)+FXs(3,0))/mRs

        rungeKab(1,0,1)=dt*(XabOA(1))
        rungeKab(2,0,1)=dt*(XabOA(2))
        rungeKab(3,0,1)=dt*(XabOA(3))
    

        do j=0,0
          do i=1,3     
            Xsdt(i,j)=Xsdt0(i,j) + rungeKob(i,j,1)/2.d0
          end do
        end do
        do i=1,3
          angleXs(i)=angleXs0(i) + rungeKab(i,1,1)/2.d0
          XobOA(i)=Xsdt(i,0)
        end do
        call forthSinker
        call Coordinate_Converter_OtoA_Sinker
        rungeKob(1,0,2)=dt*(mRs*Xsdt(3,1)*Xsdt(2,0)-mRs*Xsdt(2,1)*Xsdt(3,0)+FXs(1,0))/mRs
        rungeKob(2,0,2)=dt*(mRs*Xsdt(1,1)*Xsdt(3,0)-mRs*Xsdt(3,1)*Xsdt(1,0)+FXs(2,0))/mRs
        rungeKob(3,0,2)=dt*(mRs*Xsdt(2,1)*Xsdt(1,0)-mRs*Xsdt(1,1)*Xsdt(2,0)+FXs(3,0))/mRs

        rungeKab(1,0,2)=dt*(XabOA(1))
        rungeKab(2,0,2)=dt*(XabOA(2))
        rungeKab(3,0,2)=dt*(XabOA(3))


        do j=0,0
          do i=1,3
            Xsdt(i,j)=Xsdt0(i,j) + rungeKob(i,j,2)/2.d0
          end do
        end do
        do i=1,3
          angleXs(i)=angleXs0(i) + rungeKab(i,1,2)/2.d0
          XobOA(i)=Xsdt(i,0)
        end do
        call forthSinker
        call Coordinate_Converter_OtoA_Sinker
        rungeKob(1,0,3)=dt*(mRs*Xsdt(3,1)*Xsdt(2,0)-mRs*Xsdt(2,1)*Xsdt(3,0)+FXs(1,0))/mRs
        rungeKob(2,0,3)=dt*(mRs*Xsdt(1,1)*Xsdt(3,0)-mRs*Xsdt(3,1)*Xsdt(1,0)+FXs(2,0))/mRs
        rungeKob(3,0,3)=dt*(mRs*Xsdt(2,1)*Xsdt(1,0)-mRs*Xsdt(1,1)*Xsdt(2,0)+FXs(3,0))/mRs

        rungeKab(1,0,3)=dt*(XabOA(1))
        rungeKab(2,0,3)=dt*(XabOA(2))
        rungeKab(3,0,3)=dt*(XabOA(3))

        
        do j=0,0
          do i=1,3     
            Xsdt(i,j)=Xsdt0(i,j) + rungeKob(i,j,3)
          end do
        end do
        do i=1,3
          angleXs(i)=angleXs0(i) + rungeKab(i,1,3)
          XobOA(i)=Xsdt(i,0)
        end do
        call forthSinker
        call Coordinate_Converter_OtoA_Sinker
        rungeKob(1,0,4)=dt*(mRs*Xsdt(3,1)*Xsdt(2,0)-mRs*Xsdt(2,1)*Xsdt(3,0)+FXs(1,0))/mRs
        rungeKob(2,0,4)=dt*(mRs*Xsdt(1,1)*Xsdt(3,0)-mRs*Xsdt(3,1)*Xsdt(1,0)+FXs(2,0))/mRs
        rungeKob(3,0,4)=dt*(mRs*Xsdt(2,1)*Xsdt(1,0)-mRs*Xsdt(1,1)*Xsdt(2,0)+FXs(3,0))/mRs

        rungeKab(1,0,4)=dt*(XabOA(1))
        rungeKab(2,0,4)=dt*(XabOA(2))
        rungeKab(3,0,4)=dt*(XabOA(3))

 
        do j=0,0
          do i=1,3     
            Xsdt(i,j)=Xsdt0(i,j) + (rungeKob(i,j,1)+2.d0*rungeKob(i,j,2)+2.d0*rungeKob(i,j,3)+rungeKob(i,j,4))/6.d0
          end do
        end do
        do i=1,3
          angleXs(i)=angleXs0(i) + (rungeKab(i,1,1)+2.d0*rungeKab(i,1,2)+2.d0*rungeKab(i,1,3)+rungeKab(i,1,4))/6.d0
          Xs(i,0)=Xs(i,0) + (rungeKab(i,0,1)+2.d0*rungeKab(i,0,2)+2.d0*rungeKab(i,0,3)+rungeKab(i,0,4))/6.d0
        end do

END SUBROUTINE rungeKuttaSinker


SUBROUTINE UpdateSinker
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j

      XobOA(1)=0.250d0
      XobOA(2)=0.250d0
      XobOA(3)=0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,1)=Xs(1,0)+XabOA(1)
      Xs(2,1)=Xs(2,0)+XabOA(2)
      Xs(3,1)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=-0.250d0
      XobOA(2)=0.250d0
      XobOA(3)=0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,2)=Xs(1,0)+XabOA(1)
      Xs(2,2)=Xs(2,0)+XabOA(2)
      Xs(3,2)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=-0.250d0
      XobOA(2)=0.250d0
      XobOA(3)=-0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,3)=Xs(1,0)+XabOA(1)
      Xs(2,3)=Xs(2,0)+XabOA(2)
      Xs(3,3)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=0.250d0
      XobOA(2)=0.250d0
      XobOA(3)=-0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,4)=Xs(1,0)+XabOA(1)
      Xs(2,4)=Xs(2,0)+XabOA(2)
      Xs(3,4)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=0.250d0
      XobOA(2)=-0.250d0
      XobOA(3)=0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,5)=Xs(1,0)+XabOA(1)
      Xs(2,5)=Xs(2,0)+XabOA(2)
      Xs(3,5)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=-0.250d0
      XobOA(2)=-0.250d0
      XobOA(3)=0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,6)=Xs(1,0)+XabOA(1)
      Xs(2,6)=Xs(2,0)+XabOA(2)
      Xs(3,6)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=-0.250d0
      XobOA(2)=-0.250d0
      XobOA(3)=-0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,7)=Xs(1,0)+XabOA(1)
      Xs(2,7)=Xs(2,0)+XabOA(2)
      Xs(3,7)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=0.250d0
      XobOA(2)=-0.250d0
      XobOA(3)=-0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      Xs(1,8)=Xs(1,0)+XabOA(1)
      Xs(2,8)=Xs(2,0)+XabOA(2)
      Xs(3,8)=Xs(3,0)+XabOA(3)
      
      XobOA(1)=0.d0
      XobOA(2)=0.d0
      XobOA(3)=0.2675d0
      call Coordinate_Converter_OtoA_Sinker
      cable1dt(1)=((Xs(1,0)+XabOA(1))-Xs(1,9))/dt
      cable1dt(2)=((Xs(2,0)+XabOA(2))-Xs(2,9))/dt
      cable1dt(3)=((Xs(3,0)+XabOA(3))-Xs(3,9))/dt
      Xs(1,9)=Xs(1,0)+XabOA(1)
      Xs(2,9)=Xs(2,0)+XabOA(2)
      Xs(3,9)=Xs(3,0)+XabOA(3)

      do j=0,1
        do i=1,3
          Xsdt0(i,j)=Xsdt(i,j)
        end do
      end do
      
      do i=1,3
       angleXs0(i)=angleXs(i)
      end do
        
END SUBROUTINE UpdateSinker


SUBROUTINE OutputSinker
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i
      
     do i=0,9
         write(61,*) Xs(1,i), Xs(2,i), Xs(3,i)
      end do      
      
      write(91,*) 't=', t 
      write(91,*) 'FXs(1,0)=', FXs(1,0)
      write(91,*) 'FXs(2,0)=', FXs(2,0)
      write(91,*) 'FXs(3,0)=', FXs(3,0)
      write(91,*) 'T1x=', T1x
      write(91,*) 'T1y=', T1y
      write(91,*) 'T1z=', T1z
      write(161,'(3(f16.8,a1),f16.8)') t, ',', FXs(1,0), ',', FXs(2,0), ',', FXs(3,0)
      write(81,'(3(f16.8,a1),f16.8)') t, ',', angleXs(1), ',', angleXs(2), ',', angleXs(3)
      write(151,'(3(f16.8,a1),f16.8)') t, ',', Xsdt(1,0), ',', Xsdt(2,0), ',', Xsdt(3,0)
     
END SUBROUTINE OutputSinker

SUBROUTINE ANCFSinker
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, l, n, f
      REAL(DP)                            :: FZ1,FZ2, Error
 
      do while(emax1>tolerance1)

        counter1=counter1+1
        
        do k=1,nelm1
          call MassMat1(k) 
          call KlMat1(k)
          call KbMat1(k)
          call KallMat1(k)
        end do
        
        call MassMat_TotStructure1
        call KMat_TotStructure1
        call KMat_fin1
        call QeMat1
        call Rev_QeMat1
        call GaussSeidel1
      end do
      
        call Update_velacc1
        call Cal_F1

       write(*,*)'counter1=', counter1  

END SUBROUTINE ANCFSinker

SUBROUTINE InitialSettingsCable1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, n
    
      do i=1,nelm1
        rho1(i)=1030.d0
        d_in1(i)=0.d0
        d_out1(i)=20.00d-3
        A1(i)=pi_d/4.d0*(d_out1(i)**2-d_in1(i)**2)
        Iyy1(i)=pi_d/64.d0*(d_out1(i)**4-d_in1(i)**4)
      end do
      
        rmd1=((mRs-rhow*volRs)*Grv*length1)/(A1(1)*E1)
        length1=length1+rmd1
      
      do i=1,nelm1 
        l_elm1(i)=length1/real(nelm1)
        m_elm1(i)=rho1(i)*l_elm1(i)*A1(i)
      end do
          
    do i=1,nvtx1
        X1(1,0,i)=0.d0
        X1(2,0,i)=0.d0
        X1(3,0,i)=0.d0-real(i-1)*l_elm1(i-1)
        X1(1,1,i)=0.d0
        X1(2,1,i)=0.d0
        X1(3,1,i)=-1.d0
        X10(1,0,i)=0.d0
        X10(2,0,i)=0.d0
        X10(3,0,i)=0.d0-real(i-1)*l_elm1(i-1)
        X10(1,1,i)=0.d0
        X10(2,1,i)=0.d0
        X10(3,1,i)=-1.d0
    
        X1dt(1,0,i)=0.d0
        X1dt(2,0,i)=0.d0
        X1dt(3,0,i)=0.d0
        X1dt(1,1,i)=0.d0
        X1dt(2,1,i)=0.d0
        X1dt(3,1,i)=0.d0
        
        X1dt0(1,0,i)=0.d0
        X1dt0(2,0,i)=0.d0
        X1dt0(3,0,i)=0.d0
        X1dt0(1,1,i)=0.d0
        X1dt0(2,1,i)=0.d0
        X1dt0(3,1,i)=0.d0
        
        X1dt2(1,0,i)=0.d0
        X1dt2(2,0,i)=0.d0
        X1dt2(3,0,i)=0.d0
        X1dt2(1,1,i)=0.d0
        X1dt2(2,1,i)=0.d0
        X1dt2(3,1,i)=0.d0
        
        X1dt20(1,0,i)=0.d0
        X1dt20(2,0,i)=0.d0
        X1dt20(3,0,i)=0.d0
        X1dt20(1,1,i)=0.d0
        X1dt20(2,1,i)=0.d0
        X1dt20(3,1,i)=0.d0
      end do
     
     FextXs_t=0.d0
     FextYs_t=0.d0
     FextZs_t=(rhow*volRs-mRs)*Grv
      
END SUBROUTINE InitialSettingsCable1

SUBROUTINE Extend_cable
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l
      
!       down=-0.5d0
       down=0.0d0
      
       length1=length1-down*dt     
         
      do i=1,nelm1
        l_elm10(i)=l_elm1(i)
        l_elm1(i)=length1/real(nelm1)
        m_elm1(i)=rho1(i)*l_elm1(i)*A1(i)
      end do 
         
      X11(1,0,1)=X1(1,0,1)
      X11(2,0,1)=X1(2,0,1)
      X11(3,0,1)=X1(3,0,1)  
      X00(1,0,1)=X10(1,0,1)
      X00(2,0,1)=X10(2,0,1)
      X00(3,0,1)=X10(3,0,1)  
         
      do i=2,nvtx1-1
        xi=(l_elm1(i)-l_elm10(i))/l_elm10(i)
        S1=1.d0-3.d0*xi**2+2.d0*xi**3
        S2=l_elm1(i)*(xi-2.d0*xi**2+xi**3)
        S3=3.d0*xi**2-2.d0*xi**3
        S4=l_elm1(i)*(-xi**2+xi**3)
        
        X11(1,0,i)=S1*X1(1,0,i)+S2*X1(1,1,i)+S3*X1(1,0,i+1)+S4*X1(1,1,i+1)
        X11(2,0,i)=S1*X1(2,0,i)+S2*X1(2,1,i)+S3*X1(2,0,i+1)+S4*X1(2,1,i+1)
        X11(3,0,i)=S1*X1(3,0,i)+S2*X1(3,1,i)+S3*X1(3,0,i+1)+S4*X1(3,1,i+1) 
        
        X00(1,0,i)=S1*X10(1,0,i)+S2*X10(1,1,i)+S3*X10(1,0,i+1)+S4*X10(1,1,i+1)
        X00(2,0,i)=S1*X10(2,0,i)+S2*X10(2,1,i)+S3*X10(2,0,i+1)+S4*X10(2,1,i+1)
        X00(3,0,i)=S1*X10(3,0,i)+S2*X10(3,1,i)+S3*X10(3,0,i+1)+S4*X10(3,1,i+1)         
      end do 
      
      X11(1,0,nvtx1)=Xs(1,9)
      X11(2,0,nvtx1)=Xs(2,9)
      X11(3,0,nvtx1)=Xs(3,9)   
      X00(1,0,nvtx1)=X10(1,0,nvtx1)
      X00(2,0,nvtx1)=X10(2,0,nvtx1)
      X00(3,0,nvtx1)=X10(3,0,nvtx1)

      do i=1,nvtx1
        do k=1,3
         X1(k,0,i)=X11(k,0,i)
         X10(k,0,i)=X00(k,0,i)         
        end do
      end do
      
END SUBROUTINE Extend_cable

SUBROUTINE MassMat1(k)
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l
      
        do i=1,12
          do j=1,12
            Mass1(i,j,k)=0.d0
          end do
        end do
      
        do i=1,3
          Mass1(i   ,i   ,k)=    (13.d0/ 35.d0)                             *m_elm1(k)
          Mass1(i   ,i+ 3,k)=  (11.d0/210.d0)        *l_elm1(k)       *m_elm1(k)
          Mass1(i   ,i+ 6,k)=  ( 9.d0/70.d0)                               *m_elm1(k)         
          Mass1(i   ,i+ 9,k)= -( 13.d0/ 420.d0)      *l_elm1(k)       *m_elm1(k)  
          
          Mass1(i+ 3,i   ,k)=    (11.d0/210.d0)       *l_elm1(k)       *m_elm1(k)
          Mass1(i+ 3,i+ 3,k)=  ( 1.d0/105.d0)        *l_elm1(k)**2   *m_elm1(k)
          Mass1(i+ 3,i+ 6,k)=  ( 13.d0/ 420.d0)     *l_elm1(k)        *m_elm1(k)
          Mass1(i+ 3,i+ 9,k)= -( 1.d0/ 140.d0)       *l_elm1(k)**2   *m_elm1(k)

          Mass1(i+ 6,i   ,k)=    ( 9.d0/ 70.d0)                              *m_elm1(k)
          Mass1(i+ 6,i+ 3,k)=  ( 13.d0/ 420.d0)      *l_elm1(k)       *m_elm1(k)
          Mass1(i+ 6,i+ 6,k)=  ( 13.d0/  35.d0)                           *m_elm1(k)
          Mass1(i+ 6,i+ 9,k)= -( 11.d0/  210.d0)     *l_elm1(k)       *m_elm1(k)

          Mass1(i+ 9,i   ,k)=   -( 13.d0/ 420.d0)      *l_elm1(k)        *m_elm1(k)
          Mass1(i+ 9,i+ 3,k)= -( 1.d0/ 140.d0)        *l_elm1(k)**2    *m_elm1(k)
          Mass1(i+ 9,i+ 6,k)= -( 11.d0/  210.d0)     *l_elm1(k)        *m_elm1(k)
          Mass1(i+ 9,i+ 9,k)=  ( 1.d0/  105.d0)       *l_elm1(k)**2    *m_elm1(k)
        end do
      
END SUBROUTINE MassMat1

SUBROUTINE MassMat_TotStructure1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      
      
      do i=1,6*nvtx1
        do j=1,6*nvtx1
          M_st1(i,j)=0.d0
        end do
      end do
      
      do i=1,6
        do j=1,12
          M_st1(i,j)=Mass1(i,j,1)
        end do
      end do
      
      do k=2,nelm1
        do i=1,6
          ii=6*(k-1)+i
          do l=1,3
            do j=1,6
              jj=6*(k-2)+6*(l-1)+j
              if(l==1) then
                M_st1(ii,jj)=Mass1(i+6,j,k-1)
              else if(l==2) then
                M_st1(ii,jj)=Mass1(i+6,j+6,k-1)+Mass1(i,j,k)
              else
                M_st1(ii,jj)=Mass1(i,j+6,k)
              end if
            end do
          end do
        end do
      end do
      
      do i=1,6
        ii=6*nelm1+i
        do j=1,12
          jj=6*(nelm1-1)+j
          M_st1(ii,jj)=Mass1(i+6,j,nelm1)
        end do
      end do

END SUBROUTINE MassMat_TotStructure1


SUBROUTINE KlMat1(k)
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l
      REAL(DP)                            :: epsl_d, d
      
        do i=1,12
          do j=1,12
            Kl1(i,j,k)=0.d0
          end do
        end do
      
        do i=1,3
          Kl1(i   ,i    ,k)= 1.d0
          Kl1(i   ,i+3,k)= 0.d0
          Kl1(i   ,i+6,k)=-1.d0
          Kl1(i   ,i+9,k)= 0.d0
         
          Kl1(i+ 3,i    ,k)= 0.d0
          Kl1(i+ 3,i+3,k)= 0.d0
          Kl1(i+ 3,i+6,k)= 0.d0
          Kl1(i+ 3,i+9,k)= 0.d0
         
          Kl1(i+6,i    ,k)=-1.d0
          Kl1(i+6,i+3,k)= 0.d0
          Kl1(i+6,i+6,k)= 1.d0
          Kl1(i+6,i+9,k)= 0.d0
        
          Kl1(i+9,i    ,k)= 0.d0
          Kl1(i+9,i+3,k)= 0.d0
          Kl1(i+9,i+6,k)= 0.d0
          Kl1(i+9,i+9,k)= 0.d0
        end do
     
        d=dsqrt((X1(1,0,k+1)-X1(1,0,k))**2+(X1(2,0,k+1)-X1(2,0,k))**2+(X1(3,0,k+1)-X1(3,0,k))**2)
        epsl_d=(d-l_elm1(k))/d
        
        do i=1,12
          do j=1,12
            Kl1(i,j,k)=E1*A1(k)/l_elm1(k)*epsl_d*Kl1(i,j,k)
          end do
        end do
      
END SUBROUTINE KlMat1

SUBROUTINE KbMat1(k)
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l
      
        do i=1,12
          do j=1,12
            Kb1(i,j,k)=0.d0
          end do
        end do
      
        do i=1,3
          Kb1(i   ,i   ,k)=      12.d0
          Kb1(i   ,i+ 3,k)=      6.d0*l_elm1(k)
          Kb1(i   ,i+ 6,k)= - 12.d0
          Kb1(i   ,i+ 9,k)=      6.d0*l_elm1(k) 
          
          Kb1(i+ 3,i   ,k)=      6.d0*l_elm1(k)
          Kb1(i+ 3,i+ 3,k)=   4.d0*l_elm1(k)**2
          Kb1(i+ 3,i+ 6,k)= - 6.d0*l_elm1(k)
          Kb1(i+ 3,i+ 9,k)=    2.d0*l_elm1(k)**2

          Kb1(i+ 6,i   ,k)=    -12.d0
          Kb1(i+ 6,i+ 3,k)= -  6.d0*l_elm1(k)
          Kb1(i+ 6,i+ 6,k)=  12.d0
          Kb1(i+ 6,i+ 9,k)= - 6.d0*l_elm1(k)

          Kb1(i+ 9,i   ,k)=       6.d0*l_elm1(k)
          Kb1(i+ 9,i+ 3,k)=    2.d0*l_elm1(k)**2
          Kb1(i+ 9,i+ 6,k)=  - 6.d0*l_elm1(k)
          Kb1(i+ 9,i+ 9,k)=    4.d0*l_elm1(k)**2 
        end do            
      
        do i=1,12
          do j=1,12
            Kb1(i,j,k)=E1*Iyy1(k)*Kb1(i,j,k)/(l_elm1(k)**3)
          end do
        end do
      
END SUBROUTINE KbMat1


SUBROUTINE KallMat1(k)
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l
      
        do i=1,12
          do j=1,12
            Kall1(i,j,k)=Kl1(i,j,k)+Kb1(i,j,k)
          end do
        end do

END SUBROUTINE KallMat1


SUBROUTINE KMat_TotStructure1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      
      do i=1,6*nvtx1
        do j=1,6*nvtx1
          K_st1(i,j)=0.d0
        end do
      end do
      
      do i=1,6
        do j=1,12
          K_st1(i,j)=Kall1(i,j,1)
        end do
      end do
      
      do k=2,nelm1
        do i=1,6
          ii=6*(k-1)+i
          do l=1,3
            do j=1,6
              jj=6*(k-2)+6*(l-1)+j
              if(l==1) then
                K_st1(ii,jj)=Kall1(i+6,j,k-1)
              else if(l==2) then
                K_st1(ii,jj)=Kall1(i+6,j+6,k-1)+Kall1(i,j,k)
              else
                K_st1(ii,jj)=Kall1(i,j+6,k)
              end if
            end do
          end do
        end do
      end do
      
      do i=1,6
        ii=6*nelm1+i
        do j=1,12
          jj=6*(nelm1-1)+j
          K_st1(ii,jj)=Kall1(i+6,j,nelm1)
        end do
      end do

END SUBROUTINE KMat_TotStructure1


SUBROUTINE KMat_fin1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      
      do i=1,6*nvtx1
        do j=1,6*nvtx1
          K_st1(i,j)=K_st1(i,j)+M_st1(i,j)*4.d0/(dt**2)
        end do
      end do

END SUBROUTINE KMat_fin1


SUBROUTINE QeMat1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      REAL(DP)                            :: Trq
      	
      do j=1,nvtx1
        do i=1,6
          ii=6*(j-1)+i
          if(j==1) then
            if(i==3) then
              Qe1(ii)=-0.5d0*m_elm1(j)*Grv+0.5d0*rhow*A1(j)*l_elm1(j)*Grv
            else
              Qe1(ii)= 0.0d0
            end if
          else if(j==nvtx1) then
            if(i==3) then
              Qe1(ii)=-0.5d0*m_elm1(j-1)*Grv+0.5d0*rhow*A1(j-1)*l_elm1(j-1)*Grv
            else
              Qe1(ii)= 0.0d0
            end if
          else
            if(i==3) then
              Qe1(ii)=-0.5d0*(m_elm1(j-1)+m_elm1(j))*Grv+0.5d0*rhow*A1(j-1)*l_elm1(j-1)*Grv+0.5d0*rhow*A1(j)*l_elm1(j)*Grv
            else
              Qe1(ii)= 0.0d0
            end if
          end if
        end do
      end do

      call Cal_velaccforRes1
      call Tidalcurrent
      call WaterRes1

END SUBROUTINE QeMat1


SUBROUTINE Cal_velaccforRes1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      REAL(DP)                            :: sum
      
     if(t==dt) then 
       do i=1,1
        do j=1,3
          do k=0,1
             X1dt1(j,k,i)=(X1(j,k,i)-X10(j,k,i))/dt
             X1dt21(j,k,i)=0.d0
           end do
         end do
       end do
     else 
      do i=1,1
        do j=1,2
          do k=0,1
            X1dt21(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt1(j,k,i)=X1dt0(j,k,i)+(X1dt21(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
        end do
     end do  
      do i=1,1
        do j=3,3
          do k=0,1
              X1dt1(j,k,i)=(X1(j,k,i)-X10(j,k,i))/dt
              X1dt21(j,k,i)=(X1dt1(j,k,i)-X1dt0(j,k,i))/dt
          end do
        end do
       end do
      end if  
      
      do i=2,nvtx1
        do j=1,3
          do k=0,1
            X1dt21(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt1(j,k,i)=X1dt0(j,k,i)+(X1dt21(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
        end do
      end do   

      X1dt1(1,0,nvtx1)=cable1dt(1)
      X1dt1(2,0,nvtx1)=cable1dt(2)
      X1dt1(3,0,nvtx1)=cable1dt(3)
      X1dt21(1,0,nvtx1)=(cable1dt(1)-X1dt0(1,0,nvtx1))/dt
      X1dt21(2,0,nvtx1)=(cable1dt(2)-X1dt0(2,0,nvtx1))/dt
      X1dt21(3,0,nvtx1)=(cable1dt(3)-X1dt0(3,0,nvtx1))/dt

      do i=nvtx1,nvtx1
        do j=1,3
          do k=1,1
            X1dt21(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt1(j,k,i)=X1dt0(j,k,i)+(X1dt21(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
        end do
      end do   

END SUBROUTINE Cal_velaccforRes1


SUBROUTINE WaterRes1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k, ii, jj, n

      do i=1,nelm1
         u1(i)=(X1dt1(1,0,i)-Flow1dt(1,0,i))**2+(X1dt1(2,0,i)-Flow1dt(2,0,i))**2
         
         if ((X1dt1(1,0,i)-Flow1dt(1,0,i))==0.d0) then
           phi1(i)=PI/2.d0
         else
           phi1(i)=atan(abs((X1dt1(2,0,i)-Flow1dt(2,0,i))/(X1dt1(1,0,i)-Flow1dt(1,0,i))))
         end if
         
         if ((X1dt1(1,0,i)-Flow1dt(1,0,i))<0.d0) then
            ux1(i)=-1
         else 
            ux1(i)=1
         end if
         
         if ((X1dt1(2,0,i)-Flow1dt(2,0,i))<0.d0) then
            uy1(i)=-1
         else 
            uy1(i)=1
         end if
            
          dx1(i)=X1(1,0,i+1)-X1(1,0,i)
          dy1(i)=X1(2,0,i+1)-X1(2,0,i)
          dz1(i)=X1(3,0,i+1)-X1(3,0,i)
          lk1(i)=dsqrt((dx1(i)**2+dy1(i)**2+dz1(i)**2))
          
          u1x(i)=-(X1dt1(1,0,i)-Flow1dt(1,0,i))
          u1y(i)=-(X1dt1(2,0,i)-Flow1dt(2,0,i))
          
!--------r1(i)=dsqrt(1-cos(theta)^2)=sin(theta)----------------
          r1(i)=dsqrt(1-(((u1x(i)*dx1(i))+(u1y(i)*dy1(i)))/(dsqrt(u1x(i)**2+u1y(i)**2)*lk1(i)))**2)
!---------------------------------------------------------------------------          
          Ap1(i)=d_out1(i)*lk1(i)*r1(i)
          Fd1(i)=-0.5d0*rhow*Cd*Ap1(i)*u1(i)

          Fdx10(i)=Fd1(i)*cos(phi1(i))*ux1(i)
          Fdy10(i)=Fd1(i)*sin(phi1(i))*uy1(i)
          Fdz10(i)=0.d0
      end do

      do i=1,nvtx1
         if (i==1) then
            Fdx1(i)=Fdx10(i)/2.d0
            Fdy1(i)=Fdy10(i)/2.d0
            Fdz1(i)=Fdz10(i)/2.d0         
         else if (i==nvtx1) then
            Fdx1(i)=Fdx10(i-1)/2.d0
            Fdy1(i)=Fdy10(i-1)/2.d0
            Fdz1(i)=Fdz10(i-1)/2.d0
         else
            Fdx1(i)=(Fdx10(i-1)+Fdx10(i))/2.d0
            Fdy1(i)=(Fdy10(i-1)+Fdy10(i))/2.d0
            Fdz1(i)=(Fdz10(i-1)+Fdz10(i))/2.d0    
         end if
      end do
            
      do i=1,nvtx1      
          Qe1(6*(i-1)+1)=Qe1(6*(i-1)+1)+Fdx1(i)
          Qe1(6*(i-1)+2)=Qe1(6*(i-1)+2)+Fdy1(i)
          Qe1(6*(i-1)+3)=Qe1(6*(i-1)+3)+Fdz1(i)     
       end do

END SUBROUTINE WaterRes1


SUBROUTINE Rev_QeMat1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      REAL(DP)                            :: sum
      
      do ii=1,6*nvtx1
        sum=0.d0
        do i=1,nvtx1
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k            
              sum=sum+M_st1(ii,jj)*X10(k,j,i)*4.d0/(dt**2)+M_st1(ii,jj)*(X1dt0(k,j,i)*4.d0/dt+X1dt20(k,j,i))
            end do
          end do
        end do
        Qe1(ii)=Qe1(ii)+sum
      end do
            
END SUBROUTINE Rev_QeMat1


SUBROUTINE GaussSeidel1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, i1, i2, j1, j2, k1, k2, l1, l2, ii, jj, kk
      REAL(DP)                            :: sum, theta, delta

      emax1=0.d0              
      
      i=int(t/30.d0)+1
     
      i1=1
      X1(1,0,i1)=X10(1,0,i1)+dt*top_x(i)
      X1(2,0,i1)=X10(2,0,i1)+dt*top_y(i)
      X1(3,0,i1)=0.d0+0.5d0*sin (2.d0*PI*(t/5.d0))
 
      do i1=1,1
        do j1=1,1
          do k1=1,3
            ii=6*(i1-1)+3*j1+k1
            sum=0.d0
            do i2=i1,i1+1
              do j2=0,1
                do k2=1,3
                  jj=6*(i2-1)+3*j2+k2
                  sum=sum+K_st1(ii,jj)*X1(k2,j2,i2)
                end do
              end do
            end do
            delta=(Qe1(ii)-sum)/K_st1(ii,ii)
            if(abs(delta)>emax1) then
              emax1=abs(delta)           
            end if
            X1(k1,j1,i1)=alp*delta+X1(k1,j1,i1)
          end do
        end do
      end do     
     
      do i1=2,nvtx1-1
        do j1=0,1
          do k1=1,3
            ii=6*(i1-1)+3*j1+k1
            sum=0.d0
            do i2=i1-1,i1+1
              do j2=0,1
                do k2=1,3
                  jj=6*(i2-1)+3*j2+k2
                  sum=sum+K_st1(ii,jj)*X1(k2,j2,i2)
                end do
              end do
            end do
            delta=(Qe1(ii)-sum)/K_st1(ii,ii)
            if(abs(delta)>emax1) then
              emax1=abs(delta)
            end if
            X1(k1,j1,i1)=alp*delta+X1(k1,j1,i1)
          end do
        end do
      end do
      
      i1=nvtx1
!========katan no heni no settei======       
      X1(1,0,i1)=Xs(1,9)
      X1(2,0,i1)=Xs(2,9)
      X1(3,0,i1)=Xs(3,9)
!============================

        do j1=1,1
          do k1=1,3
            ii=6*(i1-1)+3*j1+k1
            sum=0.d0
            do i2=i1-1,i1
              do j2=0,1
                do k2=1,3
                  jj=6*(i2-1)+3*j2+k2
                  sum=sum+K_st1(ii,jj)*X1(k2,j2,i2)
                end do
              end do
            end do
            delta=(Qe1(ii)-sum)/K_st1(ii,ii)
            if(abs(delta)>emax1) then
              emax1=abs(delta)
            end if
            X1(k1,j1,i1)=alp*delta+X1(k1,j1,i1)
          end do
        end do
      
END SUBROUTINE GaussSeidel1


SUBROUTINE Update_velacc1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj
      REAL(DP)                            :: sum
      
      
     if(t==dt) then 
       do i=1,1
        do j=1,3
          do k=0,1
             X1dt(j,k,i)=(X1(j,k,i)-X10(j,k,i))/dt
             X1dt2(j,k,i)=0.d0
           end do
         end do
       end do
     else 
      do i=1,1
        do j=1,2
          do k=0,1
            X1dt2(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt(j,k,i)=X1dt0(j,k,i)+(X1dt2(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
        end do
      end do  
      do i=1,1
        do j=3,3
          do k=0,1
              X1dt(j,k,i)=(X1(j,k,i)-X10(j,k,i))/dt
              X1dt2(j,k,i)=(X1dt(j,k,i)-X1dt0(j,k,i))/dt
          end do
        end do
      end do
     end if      
             
     do i=2,nvtx1-1
       do j=1,3
          do k=0,1
            X1dt2(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt(j,k,i)=X1dt0(j,k,i)+(X1dt2(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
       end do
     end do
       
      X1dt(1,0,nvtx1)=cable1dt(1)
      X1dt(2,0,nvtx1)=cable1dt(2)
      X1dt(3,0,nvtx1)=cable1dt(3)
      X1dt2(1,0,nvtx1)=(cable1dt(1)-X1dt0(1,0,nvtx1))/dt
      X1dt2(2,0,nvtx1)=(cable1dt(2)-X1dt0(2,0,nvtx1))/dt
      X1dt2(3,0,nvtx1)=(cable1dt(3)-X1dt0(3,0,nvtx1))/dt
      
      do i=nvtx1,nvtx1
        do j=1,3
          do k=1,1
            X1dt2(j,k,i)=(4.d0*X1(j,k,i)-4.d0*X10(j,k,i)-4.d0*X1dt0(j,k,i)*dt)/(dt**2)-X1dt20(j,k,i)
            X1dt(j,k,i)=X1dt0(j,k,i)+(X1dt2(j,k,i)+X1dt20(j,k,i))*dt/2.d0
          end do
        end do
      end do

      do i=1,nvtx1
        do j=1,3
          do k=0,1
            X10(j,k,i)=X1(j,k,i)
            X1dt0(j,k,i)=X1dt(j,k,i)
            X1dt20(j,k,i)=X1dt2(j,k,i)     
          end do
        end do
      end do

END SUBROUTINE Update_velacc1
 
 
SUBROUTINE Cal_F1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj, kk
      REAL(DP)                            :: sum
      
      do ii=7,9
       sum=0.d0
        do i=1,2
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k
              sum=sum+Mass1(ii,jj,nelm1)*X1dt2(k,j,nvtx1-2+i)
            end do
          end do
        end do
        Me1(ii)=sum
      end do
      
      do ii=7,9
       sum=0.d0
        do i=1,2
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k
              sum=sum+Kl1(ii,jj,nelm1)*X1(k,j,nvtx1-2+i)
            end do
          end do
        end do
        Qkl1(ii)=sum
      end do
      
      do ii=7,9
       sum=0.d0
        do i=1,2
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k
              sum=sum+Kb1(ii,jj,nelm1)*X1(k,j,nvtx1-2+i)
            end do
          end do
        end do
        Qkb1(ii)=sum
      end do
      
      do ii=1,3
       sum=0.d0
        do i=1,2
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k
              sum=sum+Mass1(ii,jj,1)*X1dt2(k,j,i)+Kall1(ii,jj,1)*X1(k,j,i)
            end do
          end do
        end do
        Fext1(ii)=sum
      end do
      
      do ii=7,9
       sum=0.d0
        do i=1,2
          do j=0,1
            do k=1,3
              jj=6*(i-1)+3*j+k
              sum=sum+Mass1(ii,jj,nelm1)*X1dt2(k,j,nvtx1-2+i)+Kall1(ii,jj,nelm1)*X1(k,j,nvtx1-2+i)
            end do
          end do
        end do
        Fext1(ii)=sum
      end do
    
!=================Qe==========================
       call  KMat_TotStructure1
       
       do jj=1,6*nvtx1          
          Fout1(jj)=0.d0 
        end do          
      
       do ii=1,6*nvtx1
       sum=0.d0
        do i=1,nvtx1
          do k=0,1
            do j=1,3
               jj=6*(i-1)+3*k+j
               sum=sum+M_st1(ii,jj)*X1dt2(j,k,i)+K_st1(ii,jj)*X1(j,k,i)
             end do
           end do
         end do
           Fout1(ii)=sum
       end do      
!==============================================      
!      write(*,*) Fext1(7),Fout1(6*nelm1+1)
!      write(*,*) Fext1(8),Fout1(6*nelm1+2)      
!      write(*,*) Fext1(9),Fout1(6*nelm1+3)
!      pause      
      
      FextXs_t=Fext1(7)
      FextYs_t=Fext1(8)
      FextZs_t=Fext1(9)
     
END SUBROUTINE Cal_F1


SUBROUTINE OutputCable1
      USE nrtype
      USE variables
      IMPLICIT NONE
      INTEGER(I4B)                        :: i, j, k,l, ii, jj, ipnt
      REAL(DP)                            :: dum1, dum2, absomg, omg, D
      
     do i=1,nvtx1
       if(mod(i,1)==0) then
        dum1=sqrt(X1(1,1,i)**2+X1(2,1,i)**2+X1(3,1,i)**2)
        if(i==1) then
          dum2=sqrt((X1(1,0,i+1)-X1(1,0,i))**2+(X1(2,0,i+1)-X1(2,0,i))**2+(X1(3,0,i+1)-X1(3,0,i))**2)/l_elm1(i)
        else
          dum2=sqrt((X1(1,0,i)-X1(1,0,i-1))**2+(X1(2,0,i)-X1(2,0,i-1))**2+(X1(3,0,i)-X1(3,0,i-1))**2)/l_elm1(i-1)
        end if
        write(51,*) X1(1,0,i), X1(2,0,i), X1(3,0,i)
        write(56,*) X1(1,0,i), X1(2,0,i), X1(3,0,i)
       else
       end if
      end do

     write(71,'(6(f16.8,a1),f16.8)') t, ',', Fext1(1),',',  Fext1(2), ',', Fext1(3), ',', Fext1(7), ',', Fext1(8), ',', Fext1(9)
     write(73,'(3(f16.8,a1),f16.8)') t, ',', FextXs_t,',', FextYs_t, ',', FextZs_t
     write(101,'(3(f16.8,a1),f16.8)') t, ',', Me1(7), ',', Me1(8), ',', Me1(9)
     write(111,'(3(f16.8,a1),f16.8)') t, ',', Qkb1(7), ',', Qkb1(8), ',', Qkb1(9)
     write(121,'(3(f16.8,a1),f16.8)') t, ',', Qkl1(7), ',', Qkl1(8), ',', Qkl1(9)
     write(131,'(3(f16.8,a1),f16.8)') t, ',', X1dt2(1,0,nvtx1), ',', X1dt2(2,0,nvtx1), ',', X1dt2(3,0,nvtx1)
     write(141,'(3(f16.8,a1),f16.8)') t, ',', X1dt(1,0,nvtx1), ',', X1dt(2,0,nvtx1), ',', X1dt(3,0,nvtx1)

     D=dsqrt((X1(1,0,1)-X1(1,0,nvtx1))**2+(X1(2,0,1)-X1(2,0,nvtx1))**2)

!-------------------------tansindo-------------------------------------------------------------------------
             write (31,'(7(f16.8,a1),f16.8)')  t, ',',   X1(1,0,1), ',',   X1(2,0,1), ',',   X1(3,0,1), ',',   X1(1,0,nvtx1), ',',   X1(2,0,nvtx1), ',',   X1(3,0,nvtx1), ',', D
!-------------------------------------------------------------------------------------------------------------

!-------------------------extend_cable-------------------------------------------------------------------------
             write (21,'(3(f16.8,a1),f16.8)')  t, ',',   down, ',',   length1, ',',   l_elm1(1)
!-------------------------------------------------------------------------------------------------------------

!-----------------Fout(1,0,i)--------------------------------------------------------------------------------
      write (28,*) 't=', ',', t
      do i=1,nvtx1
          write(28,'(3(f16.8,a1),f16.8)') i, ',', Fout1(6*(i-1)+1), ',',Fout1(6*(i-1)+2), ',', Fout1(6*(i-1)+3)
      end do
!--------------------------------------------------------------------------------------------------------------

!------------Fdx-----------------------------------------------------------       
       write (20,*) 't=', ',', t
       do i=1,nvtx1
              write(20,'(i2,a,2(f16.8,a1),f16.8)')  i, ',', X1dt(1,0,i), ',', Flow1dt(1,0,i), ',', Fdx1(i)
        end do    
!------------------------------------------------------------------------------   

END SUBROUTINE OutputCable1

SUBROUTINE Coordinate_Converter_AtoO_Sinker
      USE nrtype
      USE variables
      IMPLICIT NONE
    
      XobAO(1)= dcos(angleXs(2))*dcos(angleXs(3))*XabAO(1)+dcos(angleXs(2))*dsin(angleXs(3))*XabAO(2)-dsin(angleXs(2))*XabAO(3)
      XobAO(2)=(dsin(angleXs(1))*dsin(angleXs(2))*dcos(angleXs(3))-dcos(angleXs(1))*dsin(angleXs(3)))*XabAO(1)+(dsin(angleXs(1))*dsin(angleXs(2))*dsin(angleXs(3))+dcos(angleXs(1))*dcos(angleXs(3)))*XabAO(2)+dsin(angleXs(1))*dcos(angleXs(2))*XabAO(3)
      XobAO(3)=(dcos(angleXs(1))*dsin(angleXs(2))*dcos(angleXs(3))+dsin(angleXs(1))*dsin(angleXs(3)))*XabAO(1)+(dcos(angleXs(1))*dsin(angleXs(2))*dsin(angleXs(3))-dsin(angleXs(1))*dcos(angleXs(3)))*XabAO(2)+dcos(angleXs(1))*dcos(angleXs(2))*XabAO(3)
      
END SUBROUTINE Coordinate_Converter_AtoO_Sinker

SUBROUTINE Coordinate_Converter_OtoA_Sinker
      USE nrtype
      USE variables
      IMPLICIT NONE

      XabOA(1)=dcos(angleXs(2))*dcos(angleXs(3))*XobOA(1)+(dsin(angleXs(1))*dsin(angleXs(2))*dcos(angleXs(3))-dcos(angleXs(1))*dsin(angleXs(3)))*XobOA(2)+(dcos(angleXs(1))*dsin(angleXs(2))*dcos(angleXs(3))+dsin(angleXs(1))*dsin(angleXs(3)))*XobOA(3)
      XabOA(2)=dcos(angleXs(2))*dsin(angleXs(3))*XobOA(1)+(dsin(angleXs(1))*dsin(angleXs(2))*dsin(angleXs(3))+dcos(angleXs(1))*dcos(angleXs(3)))*XobOA(2)+(dcos(angleXs(1))*dsin(angleXs(2))*dsin(angleXs(3))-dsin(angleXs(1))*dcos(angleXs(3)))*XobOA(3)
      XabOA(3)=-dsin(angleXs(2))*XobOA(1)+dsin(angleXs(1))*dcos(angleXs(2))*XobOA(2)+dcos(angleXs(1))*dcos(angleXs(2))*XobOA(3)
      
END SUBROUTINE Coordinate_Converter_OtoA_Sinker

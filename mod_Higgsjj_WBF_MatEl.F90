module modHiggsjj_wbf
  implicit none

  integer, parameter  :: dp = selected_real_kind(15)
  real(dp), parameter :: pi =&
       & 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: sqrt2 = &
       &1.4142135623730950488016887242096980785696718753769_dp 
  real(dp), parameter :: one = 1.0_dp
  real(dp), parameter :: two = 2.0_dp
  real(dp), parameter :: three = 3.0_dp
  real(dp), parameter :: four = 4.0_dp
  real(dp), parameter :: zero = 0.0_dp
  complex(dp), parameter :: czero = (zero,zero)
  complex(dp), parameter :: ci = (zero,one)

  real(dp), parameter :: CF = four/three
  real(dp), parameter :: CA = three
  real(dp), parameter :: xn = three

  real(dp), parameter :: avegg = one/256.0_dp
  real(dp), parameter :: aveqg = one/96.0_dp
  real(dp), parameter :: aveqq = one/36.0_dp

  real(dp), parameter :: MW = 80.398_dp
  real(dp), parameter :: mwsq = mw**2
  integer, parameter :: nf = 5
  
!  real(dp), parameter :: couplfac = one ! -- for the final result

  !-- below: to check against MCFM
  real(dp), public, parameter :: Gf = 1.16639d-5
  real(dp), parameter :: gwsq = four * sqrt2 * MW**2 * Gf
  real(dp), parameter :: vev = one/sqrt(Gf*sqrt2)
  real(dp), parameter :: couplfac = gwsq**2/two * vev! -- to check against mcfm
  !-- end MCFM block

  private

  public :: EvalAmp_vbfH

contains

  !--- vvcoupl(1) -> ??
  !--- vvcoupl(2) -> ??
  !--- vvcoupl(3) -> ??
  !--- vvcoupl(4) -> ??

  !----- p1 and p2 used to get hadronic s
  !----- unphysical kinematics to match Markus notation
  !----- 0 -> P(-p1)+P(-p2) + j(p3) + j(p4) + H(p5)
  subroutine EvalAmp_vbfH(pin,vvcoupl,me2)
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: vvcoupl(1:4)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: xa, xb, p(4,5)
    real(dp) :: shad,etot,pztot,sqrts
    real(dp) :: ud_Hud_WW, uub_Hddb_WW
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    integer :: j,k
    integer, parameter :: pn(-5:5)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)


    !-- get hadronic energy
    shad = two * scr(pin(:,1),pin(:,2))
    sqrts = sqrt(shad)

    !-- get boost:
    !-- Etot = (xa+xb)sqrts/two
    !-- pztot = (xa-xb)sqrts/two
    etot = pin(1,3)+pin(1,4)+pin(1,5)
    pztot = pin(4,3)+pin(4,4)+pin(4,5)
    xa = (etot+pztot)/sqrts
    xb = (etot-pztot)/sqrts

    !-- all outgoing kinematics
    p(:,1) = -sqrts/two * (/xa,zero,zero,xa/)
    p(:,2) = -sqrts/two * (/xb,zero,zero,-xb/)
    p(:,3) = pin(:,3)
    p(:,4) = pin(:,4)
    p(:,5) = pin(:,5)

    call spinoru(4,p,za,zb,sprod)

    call me2_tree_qqqqH_WW(vvcoupl,1,3,2,4,za,zb,sprod,ud_hud_ww)
    call me2_tree_qqqqH_WW(vvcoupl,1,3,4,2,za,zb,sprod,uub_hddb_ww)

    ud_hud_ww = ud_hud_ww*aveqq
    uub_hddb_ww = uub_hddb_ww*aveqq

    do j=-(nf-1),nf-1
       do k=-(nf-1),nf-1
          if     ((j .gt. 0) .and. (k .lt. 0)) then
             if (pn(j) .eq. -pn(k)) me2(j,k)=uub_hddb_ww             
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
             if (pn(j) .eq. -pn(k)) me2(j,k)=uub_hddb_ww
          elseif ((j .gt. 0) .and. (k .gt. 0)) then
             if (pn(j)+pn(k) .eq. +3) me2(j,k)=ud_hud_ww
          elseif ((j .lt. 0) .and. (k .lt. 0)) then
             if (pn(j)+pn(k) .eq. -3) me2(j,k)=ud_hud_ww
          endif
       enddo
    enddo
        
  end subroutine EvalAmp_vbfH


  !-----------------------------------------------------------------

  ! Notation: 0 -> qbar(p1) q(p2) qbar(p3) q(p4) H
  ! Factored out: gwsq**2/two * vev
  subroutine me2_tree_qqqqH_WW(vvcoupl,j1,j2,j3,j4,za,zb,sprod,me2)
    complex(dp), intent(in) :: vvcoupl(4)
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    integer, intent(in) :: j1, j2, j3, j4
    real(dp), intent(out) :: me2
    real(dp) :: prefac, rme2
    complex(dp) :: amp(-1:1,-1:1,-1:1,-1:1)
    real(dp) :: q1sq, q2sq, prop1, prop2
    real(dp), parameter :: colf = xn**2

    me2 = zero
    rme2 = zero

    q1sq = sprod(j1,j2)
    q2sq = sprod(j3,j4)

    prop1 = one/(q1sq-mwsq)
    prop2 = one/(q2sq-mwsq)

    prefac = couplfac * prop1 * prop2
    prefac = prefac**2

    amp = A0Hqqqq_WW(j1,j2,j3,j4,za,zb,sprod)

!------ The W-boson contribution / only one helicity plays a role
    rme2 = rme2 + amp(-1,+1,-1,+1)*conjg(amp(-1,+1,-1,+1))

    me2 = rme2 * prefac * colf

    return

  end subroutine me2_tree_qqqqH_WW


  !----------------------------------------------------------------------------------------

  ! helicity amplitudes for 0-> q(j1) qbar(j2) q(j3) qbar(j4) (H->a(p4) a(p5))

  function A0Hqqqq_WW(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0Hqqqq_WW(-1:1,-1:1,-1:1,-1:1)
    integer :: i,j
    
    A0Hqqqq_WW = czero
    
    A0Hqqqq_WW(+1,-1,+1,-1) = zb(j1,j3)*za(j4,j2)
    A0Hqqqq_WW(+1,-1,-1,+1) = zb(j1,j4)*za(j3,j2)
    A0Hqqqq_WW(-1,+1,+1,-1) = za(j1,j4)*zb(j3,j2)
    A0Hqqqq_WW(-1,+1,-1,+1) = za(j1,j3)*zb(j4,j2)
    
    return
    
  end function A0Hqqqq_WW


  !--- auxiliary functions
    !--- spinor products
  subroutine spinoru(j,p,zza,zzb,ss1) 
    implicit none 
    integer, intent(in) :: j
    complex(dp), intent(out)  :: zza(j,j), zzb(j,j)
    real(dp), intent(out) :: ss1(j,j)
    real(dp), intent(in) :: p(4,j)
    integer :: i1,i2

    zza = czero
    zzb = czero
    ss1 = zero

    do i1=1,j
       do i2=i1+1,j
          
          zza(i1,i2) = aa(p(:,i1),p(:,i2))
          zzb(i1,i2) = bb(p(:,i1),p(:,i2)) 
          zza(i2,i1) = -zza(i1,i2)
          zzb(i2,i1) = -zzb(i1,i2)

          ss1(i1,i2) = zza(i1,i2)*zzb(i2,i1)
          ss1(i2,i1) = ss1(i1,i2)

       enddo
    enddo

    return 
  end subroutine spinoru


  function aa(p2,p1)
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: aa

    aa = sum(bspa(p2)*spa(p1))

  end function aa


  ! [1 2] 

  function bb(p2,p1)
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: bb
    
    bb = sum(bspb(p2)*spb(p1))
    
  end function bb


  ! | p > spinor

  function spa(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: spa(4)
      
    spa = u0(p,1)

  end function spa


  ! | p ] spinor

  function spb(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: spb(4)
      
    spb = u0(p,-1)

  end function spb


  ! < p | spinor

  function bspa(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspa(4)
  
    bspa = ubar0(p,-1)

  end function bspa


  ! [ p | spinor

  function bspb(p) 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspb(4)
    
    bspb = ubar0(p,1)

  end function bspb

  ! -- ubar spinor, massless
  function ubar0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    integer :: sgn
    logical :: flag_nan

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)
    

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in ubar -> p0=0'
       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = py/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = asin(rrr)
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif
    

    call get_NaN(theta,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ubar-th', p0,px,py,pz
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ubar-phi', p0,px,py,pz
       stop
    endif


    if (i.eq.1) then
       ubar0(1)=czero
       ubar0(2)=czero
       ubar0(3)=sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       ubar0(4)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
    elseif(i.eq.-1) then
       ubar0(1)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       ubar0(2)=-sqrt2*sqrt(cp0)*cmplx(abs(cos(theta/two)),0.0_dp,kind=dp)
       ubar0(3)=czero
       ubar0(4)=czero
    else
       stop 'ubar0: i out of range'
    endif

  end function ubar0


  ! -- u0  spinor, massless
  function u0(p,i)
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: u0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    integer :: sgn
    logical :: flag_nan
    
    flag_nan = .false.
    
    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)
    
    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
       pause
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = three*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif
    

    call get_NaN(phi,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif
    
    if (i.eq.-1) then
       u0(1)=czero
       u0(2)=czero
       u0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       u0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.1) then
       u0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       u0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       u0(3)=czero
       u0(4)=czero
    else
       stop 'u0: i out of range'
    endif

  end function u0

  subroutine get_NaN(value,flag_nan)
    implicit none
    real(dp), intent(in)  :: value
    logical, intent(out) :: flag_nan
    
    flag_nan = .false.
    
    if (.not.value.le.0.0_dp .and. .not.value.gt.0.0_dp ) then
       flag_nan =.true.
    endif
    
  end subroutine get_NaN

  function scr(p1,p2)   !scalar product of real vectors
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr

 !--- scalar products, 2*pi.pj for massless particles
  subroutine sprodu(j,p,sprod)
    implicit none
    integer, intent(in) :: j
    real(dp), intent(in) :: p(4,j)
    real(dp), intent(out) :: sprod(j,j)
    integer :: i1, i2

    sprod = zero

    do i1 = 1, j
       do i2 = i1+1, j
          sprod(i1,i2) = two * scr(p(:,i1),p(:,i2))
          sprod(i2,i1) = sprod(i1,i2)
       enddo
    enddo

    return

  end subroutine sprodu
  
 



end module modHiggsjj_wbf

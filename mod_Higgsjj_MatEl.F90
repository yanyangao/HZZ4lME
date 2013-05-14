module modHiggsjj
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

  real(dp), parameter :: nf = 5.0_dp

  real(dp), parameter :: avegg = one/256.0_dp
  real(dp), parameter :: aveqg = one/96.0_dp
  real(dp), parameter :: aveqq = one/36.0_dp

  private

  public :: EvalAmp_gg_jjH

contains

  !--- all amplitudes have gs^2 (as/(six*pi*vev)) factored out
  !--- ggcoupl(1) -> scalar
  !--- ggcoupl(2) -> pseudoscalar

  !----- p1 and p2 used to get hadronic s
  !----- P(p1)+P(p2) -> j(p3) + j(p4) + H(p5)
  subroutine EvalAmp_gg_jjH(pin,ggcoupl,me2)
    real(dp), intent(in) :: pin(4,5)
    complex(dp), intent(in) :: ggcoupl(1:3)
    real(dp), intent(out) :: me2(-5:5,-5:5)
    real(dp) :: xa, xb, p(4,5)
    real(dp) :: shad,etot,pztot,sqrts
    real(dp) :: gg_hgg, gg_hqa, qg_hqg, gq_hqg
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)

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

    call me2_ggggh(ggcoupl,1,2,3,4,za,zb,sprod,gg_hgg)
    call me2_qbqggh(ggcoupl,3,4,1,2,za,zb,sprod,gg_hqa)
    call me2_qbqggh(ggcoupl,1,3,2,4,za,zb,sprod,qg_hqg)
    call me2_qbqggh(ggcoupl,2,3,1,4,za,zb,sprod,gq_hqg)

    me2(0,0) = (gg_hgg/two + nf * gg_hqa)*avegg
    me2(1,0) = (qg_hqg)*aveqg
    me2(0,1) = (gq_hqg)*aveqg
        
  end subroutine EvalAmp_gg_jjH

  !--- 0 -> g(p1) g(p2) g(p3) g(p4) [H(p4)]
  subroutine me2_ggggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,me2gg)
    complex(dp), intent(in) :: ggcoupl(1:3)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2gg
    complex(dp) :: amp1234(-1:1,-1:1,-1:1), amp1324(-1:1,-1:1,-1:1)
    complex(dp) :: aPhi1234(1:2,-1:1,-1:1,-1:1), aPhi1324(1:2,-1:1,-1:1,-1:1)
    real(dp), parameter :: col1 = 8.0_dp * CA**3 * CF
    integer :: i2,i3,i4

    me2gg = zero

    aPhi1234 = A0phigggg_pxxx(j1,j2,j3,j4,za,zb,sprod)
    aPhi1324 = A0phigggg_pxxx(j1,j3,j2,j4,za,zb,sprod)

    amp1234(:,:,:) = ggcoupl(1) * (aPhi1234(1,:,:,:) + aPhi1234(2,:,:,:)) + & !-- scalar
         ggcoupl(2) * (-ci) * (aPhi1234(1,:,:,:) - aPhi1234(2,:,:,:)) !-- pseudoscalar

    amp1324(:,:,:) = ggcoupl(1) * (aPhi1324(1,:,:,:) + aPhi1324(2,:,:,:)) + & !-- scalar
         ggcoupl(2) * (-ci) * (aPhi1324(1,:,:,:) - aPhi1324(2,:,:,:)) !-- pseudoscalar

    do i2 = -1, 1, 2
       do i3 = -1, 1, 2
          do i4 = -1,1,2
             me2gg = me2gg + abs(amp1234(i2,i3,i4))**2
             me2gg = me2gg + abs(amp1324(i2,i3,i4))**2
             me2gg = me2gg + real(amp1234(i2,i3,i4)*conjg(amp1324(i3,i2,i4)),kind=dp)
          enddo
       enddo
    enddo

    !-- color factors and all
    me2gg = me2gg * col1

    !-- extra hels
    me2gg = me2gg * two

    return

  end subroutine me2_ggggh


  subroutine me2_qbqggh(ggcoupl,j1,j2,j3,j4,za,zb,sprod,me2q)
    complex(dp), intent(in) :: ggcoupl(1:3)
    integer, intent(in) :: j1,j2,j3,j4
    complex(dp), intent(in) :: za(4,4), zb(4,4)
    real(dp), intent(in) :: sprod(4,4)
    real(dp), intent(out) :: me2q
    complex(dp) :: amp1234(-1:1,-1:1), amp1243(-1:1,-1:1)
    complex(dp) :: aPhi1234(1:2,-1:1,-1:1), aPhi1243(1:2,-1:1,-1:1)
    real(dp), parameter :: colf1 = four * xn * CF**2
    real(dp), parameter :: colf2 = four * xn * CF * (CF - CA/2.0_dp)
    integer :: i3,i4

    me2q = zero

    aPhi1234 = A0phiqbqgg_mpxx(j1,j2,j3,j4,za,zb,sprod)
    aPhi1243 = A0phiqbqgg_mpxx(j1,j2,j4,j3,za,zb,sprod)

    amp1234(:,:) = ggcoupl(1) * (aPhi1234(1,:,:) + aPhi1234(2,:,:)) + & !-- scalar
         ggcoupl(2) * (-ci) * (aPhi1234(1,:,:) - aPhi1234(2,:,:)) !-- pseudoscalar

    amp1243(:,:) = ggcoupl(1) * (aPhi1243(1,:,:) + aPhi1243(2,:,:)) + & !-- scalar
         ggcoupl(2) * (-ci) * (aPhi1243(1,:,:) - aPhi1243(2,:,:)) !-- pseudoscalar

    do i3 = -1,1,2
       do i4 = -1,1,2
          me2q = me2q + colf1 * abs(amp1234(i3,i4))**2
          me2q = me2q + colf1 * abs(amp1243(i4,i3))**2
          me2q = me2q + colf2 * two * real(amp1234(i3,i4)*conjg(amp1243(i4,i3)),kind=dp)
       enddo
    enddo

    !-- extra hels
    me2q = me2q * two

    return

  end subroutine me2_qbqggh


!---------------------------------------------------------------------------

!--- phi-amplitudes: 0->qb(p1)q(p2)g(p3)g(p4)
!--- iphi = 1 --> phi
!--- iphi = 2 --> phid
!--- assume qb^-,q^+
  function A0phiqbqgg_mpxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phiqbqgg_mpxx(1:2,-1:1,-1:1)
    integer :: j1, j2, j3, j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: zab_1_ph_4, zab_3_ph_1, zab_3_ph_2, zab_3_ph_4
    complex(dp) :: zab_1_ph_2, zab_2_ph_4
    real(dp) :: s123, s412, qsq

    A0phiqbqgg_mpxx = czero

    s123 = sprod(j1,j2) + sprod(j1,j3) + sprod(j2,j3)
    s412 = sprod(j4,j1) + sprod(j4,j2) + sprod(j1,j2)

    zab_1_ph_4 = -za(j1,j2)*zb(j2,j4) - za(j1,j3)*zb(j3,j4)
    zab_3_ph_1 = -za(j3,j2)*zb(j2,j1) - za(j3,j4)*zb(j4,j1)
    zab_3_ph_2 = -za(j3,j1)*zb(j1,j2) - za(j3,j4)*zb(j4,j2)
    zab_3_ph_4 = -za(j3,j1)*zb(j1,j4) - za(j3,j2)*zb(j2,j4)
    zab_1_ph_2 = -za(j1,j3)*zb(j3,j2) - za(j1,j4)*zb(j4,j2)
    zab_2_ph_4 = -za(j2,j1)*zb(j1,j4) - za(j2,j3)*zb(j3,j4)

    qsq = s123 + sprod(j1,j4) + sprod(j2,j4) + sprod(j3,j4)

    A0phiqbqgg_mpxx(+1,+1,-1) = - za(j1,j4)**2 * za(j2,j4)/za(j1,j2)/za(j2,j3)/za(j3,j4) 
    A0phiqbqgg_mpxx(+2,+1,-1) = - zb(j2,j3)**2 * zb(j1,j3)/zb(j1,j2)/zb(j3,j4)/zb(j4,j1)

    A0phiqbqgg_mpxx(+1,-1,+1) = za(j1,j3)**3/za(j1,j2)/za(j3,j4)/za(j4,j1) 
    A0phiqbqgg_mpxx(+2,-1,+1) = zb(j2,j4)**3/zb(j1,j2)/zb(j2,j3)/zb(j3,j4)

    A0phiqbqgg_mpxx(+1,-1,-1) = qsq**2 * za(j1,j3)**3/s123/za(j1,j2)/zab_1_ph_4/zab_3_ph_4 &
         + zab_3_ph_1 * zab_3_ph_2**2/s412/zab_3_ph_4/zb(j2,j1)/zb(j4,j1) &
         - zab_1_ph_2**2/zab_1_ph_4/zb(j3,j2)/zb(j4,j3)

    A0phiqbqgg_mpxx(+2,+1,+1) = -zab_1_ph_2**2/za(j1,j4)/za(j3,j4)/zab_3_ph_2 + &
         zab_1_ph_4**2 * zab_2_ph_4/s123/za(j1,j2)/za(j2,j3)/zab_3_ph_4 &
         + qsq**2 * zb(j4,j2)**3/s412/zab_3_ph_2/zab_3_ph_4/zb(j2,j1)

    return

  end function A0phiqbqgg_mpxx


!--- phi-amplitudes: 0->g(p1)g(p2)g(p3)g(p4)
!--- iphi = 1 --> phi
!--- iphi = 2 --> phid
!--- assume h1=+1
  function A0phigggg_pxxx(j1,j2,j3,j4,za,zb,sprod)
    complex(dp) :: A0phigggg_pxxx(1:2,-1:1,-1:1,-1:1)
    integer :: j1, j2, j3, j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)

    A0phigggg_pxxx = czero

    A0phigggg_pxxx(+2,+1,+1,+1) = A0phiggggmmmm(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+2,+1,+1,-1) = A0phiggggpmmm(j4,j1,j2,j3,zb,za,sprod)

    A0phigggg_pxxx(+2,+1,-1,+1) = A0phiggggpmmm(j3,j4,j1,j2,zb,za,sprod)

    A0phigggg_pxxx(+1,+1,-1,-1) = A0phiggggmmpp(j3,j4,j1,j2,za,zb,sprod)
    A0phigggg_pxxx(+2,+1,-1,-1) = A0phiggggmmpp(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+2,-1,+1,+1) = A0phiggggpmmm(j2,j3,j4,j1,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,+1,-1) = A0phiggggmpmp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_pxxx(+2,-1,+1,-1) = A0phiggggmpmp(j1,j2,j3,j4,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,-1,+1) = A0phiggggmmpp(j2,j3,j4,j1,za,zb,sprod)
    A0phigggg_pxxx(+2,-1,-1,+1) = A0phiggggmmpp(j4,j1,j2,j3,zb,za,sprod)

    A0phigggg_pxxx(+1,-1,-1,-1) = A0phiggggpmmm(j1,j2,j3,j4,za,zb,sprod)

    return

  end function A0phigggg_pxxx


  function A0phiggggpmmm(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggpmmm
    
    real(dp) :: s3
    complex(dp) :: zab2
        
    s3(j1,j2,j3)=sprod(j1,j2)+sprod(j2,j3)+sprod(j3,j1)
    zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

    A0phiggggpmmm = &
         +(zab2(j3,j2,j4,j1)*za(j2,j4))**2/(s3(j1,j2,j4)*sprod(j1,j2)*sprod(j1,j4)) &
         +(zab2(j4,j2,j3,j1)*za(j2,j3))**2/(s3(j1,j2,j3)*sprod(j1,j2)*sprod(j2,j3)) &
         +(zab2(j2,j3,j4,j1)*za(j3,j4))**2/(s3(j1,j3,j4)*sprod(j1,j4)*sprod(j3,j4)) &
         -za(j2,j4)/(za(j1,j2)*zb(j2,j3)*zb(j3,j4)*za(j4,j1)) &
         *(-sprod(j2,j3)*zab2(j2,j3,j4,j1)/zb(j4,j1) &
         -sprod(j3,j4)*zab2(j4,j2,j3,j1)/zb(j1,j2) &
         -s3(j2,j3,j4)*za(j2,j4))
    
    return

  end function A0phiggggpmmm


  function A0phiggggmpmp(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggmpmp
    
    A0phiggggmpmp = za(j1,j3)**4 &
         /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
    
    return

  end function A0phiggggmpmp


  function A0phiggggmmpp(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4)
    complex(dp) :: A0phiggggmmpp
    
    A0phiggggmmpp = za(j1,j2)**4 &
         /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))
    
    return

  end function A0phiggggmmpp
             

  function A0phiggggmmmm(j1,j2,j3,j4,za,zb,sprod)
    integer :: j1,j2,j3,j4
    complex(dp) :: za(4,4), zb(4,4)
    real(dp) :: sprod(4,4),qsq
    complex(dp) :: A0phiggggmmmm
    
    qsq = sprod(j1,j2)+sprod(j1,j3)+sprod(j1,j4)+sprod(j2,j3)+sprod(j2,j4)+sprod(j3,j4)
    A0phiggggmmmm=qsq**2/(zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j4,j1))
    
    return
        
  end function A0phiggggmmmm



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





end module modHiggsjj

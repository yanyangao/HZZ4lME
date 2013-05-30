include 'variables.F90'

if( DecayMode1.le.3 ) then
   M_V = M_Z
   Ga_V= Ga_Z
elseif( (DecayMode1.ge.4) .and. (DecayMode1.le.6) ) then
   M_V = M_W
   Ga_V= Ga_W    
elseif( DecayMode1.eq.7 ) then
   M_V = 0d0
   Ga_V= 0d0    
endif



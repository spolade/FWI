module fwiindices !!!************ Beginning of Module ***********!!! 
! Module for calculation of Canadian Fire Weather Index (FWI)

! The FWI System outputs three moisture codes: Fine Fuel Moisture Code (FFMC); Duff Moisture Code (DMC); Drought Code (DC)
! & three fire behavior outputs: Initial Spread Index (ISI); Build-Up Index (BUI); Fire Weather Index (FWI)

 
! Moisture codes:  

! Subroutine FFMCcalc-  Input: Temp, relative humidity, wind speed, rain, and a previous FFMC value 
!                       Output Fine Fuel Moisture Code (FFMC)
!                       call FFMCcalc(17,42,25,0,85) = 87.692980092774448

!Subroutine DMCcalc-    Input: Temp, relative humidity, rain, previous DMC value, latitude, and month 
!                       OutPut: Duff Moisture Code (DMC)
!                       call DMCcalc(17,42,0,6,45.98,4) = 8.5450511359999997

!Subroutine DCcalc-     Input: Temp, rain, the previous DC value, latititude, month 
!                       Output: Drought Code (DC)
!                       call DCcalc(17,0,15,45.98,4) = 19.013999999999999


! Fire behavior outputs:

!Subroutine ISICalc-    Input: wind speed and current FFMC value 
!                       Output: Initial Spread Index (ISI)
!                       call ISICalc(25,87.692980092774448) = 10.853661073655068

!Subroutine BUICalc-    Input: the current DMC and DC values 
!                       Output: Build-Up Index (BUI)
!                       call BUICalc(8.5450511359999997,19.013999999999999) = 8.4904265358371838

!Subroutine FWICalc-    Input: the current ISI and BUI values 
!                       Output: Fire Weather Index (FWI)
!                       call FWICalc(4,17,42,25,0,85,6,15,45.98) = 10.096371392382368



implicit none
real, save :: Mo,Rf,Ed,Ew,M,Kl,Kw,K,Pr,Ra,Re,B,Mr,V,Rd,Dr,Qo
real, save :: Fw,Ff,P,Qr,D0,Kd,Ko,Fd  
private :: Mo,Rf,Ed,Ew,M,Kl,Kw,K,Pr,Ra,Re,B,Mr,V,Rd,Dr,Qo 
private :: Fw,Ff,P,Qr,D0,Kd,Ko,Fd 

Contains


! Subroutine FFMC calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates Fine Fuel Moisture Code
!Input:
!    T  --> TEMP is the 12:00 LST temperature in degrees celsius
!    H  --> RH is the 12:00 LST relative humidity in %
!    W  --> WIND is the 12:00 LST wind speed in km/h
!    Ro --> RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
!    Fo --> FFMCPrev is the previous day's FFMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FFMCcalc(T,H,W,Ro,Fo, ffmc) 
implicit none
	real, intent(in) :: T,H,W,Ro,Fo 
	real, intent(out) :: ffmc

        Mo = (147.2*(101.0 - Fo))/(59.5 + Fo)                                        !*Eq. 1*!
        
        if (Ro > 0.5) then
                rf = Ro - 0.5                                                            !*Eq. 2*!
                
                if(Mo <= 150.0) then
                        Mr = Mo + 42.5 * Rf * exp(-100.0/(251.0-Mo)) * (1.0 - exp(-6.93/Rf)) !*Eq. 3a*!
                else
                        Mr = (Mo + 42.5 * Rf * exp(-100.0/(251.0-Mo)) * (1.0 - exp(-6.93/Rf))) + & 
                        (.0015*(Mo - 150.0)**2)*sqrt(Rf)                                     !*Eq.3b*!
                endif
 
                if (Mr > 250.0) Mr = 250.0
                Mo = Mr 
        endif
        Ed =.942 * (H**.679) + (11.0*exp((H-100.0)/10.0))+0.18*(21.1-T) * & 
        (1.0 - 1.0/exp(.1150 * H))                                                   !*Eq. 4*!

        if(Mo > Ed) then
                Ko = .424*(1.0-(H/100.0)**1.7)+(.0694*sqrt(W))*(1.0- (H/100.0)**8)      !*Eq. 6a*!
                Kd = Ko * (.581*exp(.0365*T))                                           !*Eq. 6b*!
                M = Ed + (mo-ed)/10.0 ** Kd                                             !*Eq. 8*!
        else
                Ew = .618*(H**.753) + (10.0*exp((H-100.0)/10.0)) + .18*(21.1-T) * & 
                (1.0 - 1.0/exp(.115 * H))                                                !*Eq. 5*!
                if(Mo < Ew) then
                        Kl = .424*(1.0-((100.0 - H)/100.0)**1.7) + (.0694*sqrt(W)) * & 
                        (1.0 - ((100.0 - H)/100.0)**8)
                        Kw = Kl * (.581 * exp(.0365 * T))
                        M = Ew - (Ew - Mo)/10.0**Kw 
                else
                        M= Mo                                                               !*Eq. 7a*!
                endif                                                                   !*Eq. 7b*! 
        endif                                                                       !*Eq. 9*! 


        ffmc = (59.5 * (250.0 - m)) / (147.2 + m)                                   !*Eq. 10*!
        if (ffmc > 101.0) ffmc = 101.0
        if (ffmc <= 0.0) ffmc = 0.0
        return
end subroutine


! Subroutine DMC calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates Duff Moisture Code
!Input:
!    T   --> TEMP is the 12:00 LST temperature in degrees celsius
!    H   --> RH is the 12:00 LST relative humidity in %
!    Ro  --> RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
!    Po  --> prevvious day's DMC
!    I   --> month of Year (1..12)
!    Lat --> latitude in decimal degrees 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DMCcalc(T, H, Ro, Po, I, Lat, dmc)    
implicit none
        integer, intent(in)   :: I
	real,    intent(in)   :: T, H, Ro, Lat
	real,    intent(out)  :: dmc
	real :: Po, El
 
 
    !call DayLength(Lat, I, El)


    call daylengthDynamic(Lat, I, El)

    
	if(Ro <= 1.5) then
		Pr = Po
        
	else
		Re = 0.92*Ro - 1.27                                                     !*Eq. 11*!
		Mo = 20.0 + 280.0/exp(0.023*Po)                                         !*Eq. 12*!
		if(Po <= 33.0) then
			B = 100.0 /(0.5 + 0.3*Po)                                          !*Eq. 13a*!
		else
			B= 6.2 * log(Po) - 17.2
			if(Po-65.0 <= 0.0) B = 14.0 - 1.3*log(Po)                           !*Eq.13c*!
		endif
		                                                                        !*Eq.13b*!
		Mr = Mo + (1000.0*Re) / (48.77+B*Re)                                    !*Eq. 14*!
		Pr = 43.43 * (5.6348 - log(Mr-20.0))                                    !*Eq. 15*!
	    
	endif
	
	if(T >= -1.1) then
		K = 1.894*(T+1.1) * (100.0-H) * (El*0.0001)                             !*Eq. 15*!
	else
		K = 0.0    
		                                                                        !*Eq. 17*!
	endif
	
	if(Pr < 0.0) Pr = 0.0 
	dmc = Pr + K
	if(dmc <= 0.0) dmc = 0.0
	
	
	return
end subroutine

! Subroutine DC calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates Drought Code
!Input:
!    T   --> TEMP is the 12:00 LST temperature in degrees celsius
!    Ro  --> RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
!    Do  --> prevvious day's DC
!    I   --> month of Year (1..12)
!    Lat --> latitude in decimal degrees 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DCcalc (T, Ro, D0, I, Lat, dc)
implicit none
	real, intent(in)    :: T, Ro, Lat 
	real, intent(inout) :: D0
	real, intent(out)   :: dc
	real   :: Fl
	integer,intent(in)  :: I


    call DryingFactor(Lat, I,Fl)
	if (Ro > 2.8) then
		Rd = 0.83*Ro - 1.27                                                     !*Eq. 18*!
		Qo = 800.0 * exp(-D0/400.0)                                             !*Eq. 19*!
		Qr = Qo + 3.937*Rd                                                      !*Eq. 20*!   
		Dr = 400.0*log(800.0/Qr)                                                !*Eq. 21*!  
		if(Dr > 0.0) then
			D0 = Dr 
		else
			D0 = 0.0 
		endif
	endif
		
	if(T < -2.8) then 
		V = Fl 
	else
		V = (0.36*(T+2.8) + Fl) 
	endif
	if (V <= 0.0) V = 0.0                                                       !*Eq. 22*!
	dc = D0 + 0.5*V
	return
end subroutine



! Subroutine for ISI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates Initial Spread Index
!Input:
!    F  --> FFMC current day's FFMC
!    W  --> WIND is the 12:00 LST wind speed in kph
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ISICalc(F,W, isi) 
implicit none
	real,intent(in) :: F, W 
	real,intent(out):: isi
	 
	Mo = (147.2*(101.0-F)) / (59.5+F)                                           !*Eq. 1*!
	Fw = exp(0.05039*W)                                                         !*Eq. 24*!
	Ff = 91.9*exp(-0.1386*Mo) * (1.0+(Mo**5.31)/4.93e7)                         !*Eq. 25*!
	isi = 0.208 * Fw * Ff                                                       !*Eq. 26*!
	return 
end subroutine



! Subroutine for BUI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates Buildup Index
!Input:
!    P  --> current day's Duff Moisture Code
!    D  --> current day's Drought Code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BUICalc(P,D,bui) 
implicit none
	real, intent(in) :: P,D 
	real, intent(out):: bui
	if (P <= 0.4*D) then
		bui = (0.8*P*D) / (P+0.4*D)                                            !*Eq. 27a*!
	else
		bui = P-(1.0-0.8*D/(P+0.4*D))*(0.92+(0.0114*P)**1.7)                   !*Eq. 27b*!
	endif
	if (bui<0.0) bui = 0.0
	return
end subroutine


! Subroutine for FWI 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates FWI
!Input:
!    R  --> current day's ISI
!    U  --> current day's BUI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FWICalc(R,U,fwi)
implicit none
	real, intent(in) :: R,U 
	real, intent(out):: fwi
	!real :: Fd,B
	if (U <= 80.0) then
		Fd = 0.626*U**0.809 + 2.0                                              !*Eq. 28a*!
	else
		Fd= 1000.0/(25.0+108.64*exp(-0.023*U))                                 !*Eq. 28b*!
	endif
	
	B = 0.1 * R * Fd                                                            !*Eq. 29*!
	if(B > 1.0) then
		fwi = exp(2.72 * (0.434*log(B)) **0.647)                               !*Eq. 30a*!
	else
		fwi = B                                                                !*Eq. 30b*!
	endif
	return
end subroutine


!   Subroutine for the Drying Factor
!----------------------------------------------------------------------------------------------------
!   Day-length adjustment in DC
!   'NSH': the values for the Southern hemisphere are those for the northern shifted by 6 months
!   'NSHeq': the same idea is applied, but near the equator, one same value is applied for all months

!Input:
!    Latitude  --> latitude
!    Month  --> month
subroutine DryingFactor(Latitude, Month,Fl)
implicit none
	real, intent(in)    :: Latitude 
    integer, intent(in) :: Month 
	real, intent(out)   :: Fl
	real, save :: LfN(12),LfS(12), Lfeq(12)
    data LfN /-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6/
    data LfS /6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8/
    data Lfeq /1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4/

    if (Latitude > 20.0) Fl = LfN(Month)
    !if ((Latitude > -20.0) .and. (Latitude <= 20.0)) Fl = Lfeq

    IF (Latitude > -20.0 .AND. Latitude <= 20.0) Fl = Lfeq(Month)
    !END IF

    if (Latitude <= -20.0) Fl = LfS(Month)

    !if (Latitude > 0.0) then
    !    Fl = LfN(Month)
    !else    
    !    Fl = LfS(Month)
    !endif   
	return 
end subroutine


subroutine DayLength(Latitude, MONTH, El)
!    '''Approximates the length of the day given month and latitude'''
implicit none
	real, intent(in)    :: Latitude 
    integer, intent(in) :: Month 
	real, intent(out)   :: El
	real, save :: DayLength46N(12),DayLength20N(12),DayLength20S(12),DayLength40S(12)
	
	data DayLength46N/6.5,  7.5,  9.0, 12.8, 13.9, 13.9, 12.4, 10.9,  9.4,  8.0,  7.0,  6.0/ 
	data DayLength20N/7.9,  8.4,  8.9,  9.5,  9.9, 10.2, 10.1,  9.7,  9.1,  8.6,  8.1,  7.8/

	data DayLength20S/10.1,  9.6,  9.1,  8.5,  8.1,  7.8,  7.9,  8.3,  8.9,  9.4,  9.9, 10.2/ 
	data DayLength40S/11.5, 10.5,  9.2,  7.9,  6.8,  6.2,  6.5,  7.4,  8.7, 10.0, 11.2, 11.8/

    if(Latitude<= 90 .AND. Latitude > 33)           EL = DayLength46N(MONTH)
    if(Latitude <= 33 .AND. Latitude > 0.0)         El = DayLength20N(MONTH)
    if(Latitude <= 0.0 .AND. Latitude > -30.0)      El = DayLength20S(MONTH)
    if(Latitude <= -30.0 .AND. Latitude >= -90.0)   El = DayLength40S(MONTH)

    return 
end subroutine    
    

subroutine daylengthDynamic(lat_in, vdate, length)
implicit none
        integer, intent(in) :: vdate
	real, intent(in)    :: lat_in 
	real, intent(out)   :: length
        real   :: val_for_arccos, sun_dec, lat, sunset_hour_angle
        real(16), parameter :: PI_16 = 4 * atan (1.0_16)
        integer :: dayNo

        call dayOfYear(vdate, dayNo)

        lat = lat_in *PI_16/ 180  ! degree -> radian
        sun_dec = 0.409 * sin(2 * PI_16 / 365 * dayNo - 1.39)  ! equation 24
        ! preparing equation 25, with special cases handled
        val_for_arccos = -tan(lat) * tan(sun_dec)
        !val_for_arccos[np.where(val_for_arccos < -1)] = -1
        !val_for_arccos[np.where(val_for_arccos > 1)] = 1
        
        IF (abs(lat_in) < 90) THEN
         if (val_for_arccos < -1) val_for_arccos = -1
         if (val_for_arccos  > 1) val_for_arccos = 1

        ELSE
         if (val_for_arccos < -1) val_for_arccos = 1
         if (val_for_arccos  > 1) val_for_arccos = -1
        END IF

      
        
        sunset_hour_angle = acos(val_for_arccos)  ! equation 25
        length = 24 / PI_16 * sunset_hour_angle  !equation 34
        return 
end subroutine 


subroutine dayOfYear(vdate, dayNo)
implicit none
        integer, intent(in)  :: vdate
        logical :: year_leap
        integer, intent(out)  :: dayNo
        integer :: day, month, year,i , daysPreMonth(12)
        integer, parameter :: daysInMonth(12) = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]


        day = MOD(vdate/1,100)
        month = MOD(vdate/100,100)
        year = MOD(vdate/10000,100000000)

      if ( modulo( year, 4 ) == 0) THEN
         if ( modulo( year, 100 ) ==  0 .and. modulo( year,400 ) /= 0 ) then
            year_leap = .false.
         else
           year_leap  = .true.
         end if
     else
        year_leap  = .false.
     end if


      daysPreMonth(1) = 0
      do i = 2, 12
         daysPreMonth(i) = daysPreMonth(i-1) + daysInMonth(i-1)
      end do

        dayNo = daysPreMonth( month) + day
       if ( year_leap .and. month > 2 ) dayNo = dayNo + 1



        return
end subroutine


subroutine test_FireSeason_start(tasData, startTH, StartDate)
implicit none
         real :: tasData, startTH
       ! integer, dimension(:, :), allocatable :: StartDate, EndDate
       ! integer :: numLons, numLats, numTimes
       ! integer :: i
        logical :: StartDate

         !real, intent(in) :: tasData, startTH
         !logical, intent(out):: StartDate

       where (tasData > startTH)
             tasData = 1.0
       elsewhere
             tasData = 0.0
       endwhere


        if ( (sum(tasData)/3) ==  1 ) then
            StartDate = .true. 
         else
           StartDate  = .false.
         end if


        return
end subroutine


!Lawson equations:
!
!All of these equations take the current DMC and DC values and return moisture content as a % value
!
!LawsonEq1 - DMC National Standard and Coastal B.C. CWH (2.5-4 cm)^2
!          - LawsonEq1(8.5450511359999997)  = 250.7553985454235
!
!LawsonEq2 - Southern Interior B.C.3 (2-4 cm)^2
!          - LawsonEq2(8.5450511359999997)  = 194.93023948344205
!
!LawsonEq3 - Southern Yukon - Pine/White Spruce
!                             Feather moss, Sphagnum and Undifferentiated duff (2-4 cm)^2
!          - LawsonEq3(8.5450511359999997)  = 442.82109267231488
!
!LawsonEq4 - Southern Yukon - Pine/White Spruce
!                             Reindeer lichen (2-4 cm)^2
!          - LawsonEq4(8.5450511359999997)  = 746.02210402093272
!
!LawsonEq5 - Southern Yukon - White Spruce
!                             White spruce/feather moss (2-4 cm)^2
!          - LawsonEq5(8.5450511359999997)  = 853.2397847094652


!    '''National Standard and Best-fit Non-linear Regression Equations
!Linking DMC to Forest Floor Moisture Content in
!Coastal B.C., Southern Interior B.C. and Southern Yukon

!DMC National Standard and Coastal B.C. CWH (2.5-4 cm)^2

!LawsonEq1(8.5450511359999997)  = 250.7553985454235'''

subroutine LawsonEq1(DMC, Moisture)
implicit none
	real, intent(in)    :: DMC 
	real, intent(out)   :: Moisture
    Moisture = exp((DMC - 244.7) / (-43.4))+20.0
    

    	
	return 
end subroutine



!def LawsonEq2(DMC):
!    '''National Standard and Best-fit Non-linear Regression Equations
!Linking DMC to Forest Floor Moisture Content in
!Coastal B.C., Southern Interior B.C. and Southern Yukon
!Southern Interior B.C.3 (2-4 cm)^2
!LawsonEq2(8.5450511359999997)  = 194.93023948344205'''
!    return math.exp((DMC-223.9)/-41.7)+20.0

subroutine LawsonEq2(DMC, Moisture)
implicit none
	real, intent(in)    :: DMC 
	real, intent(out)   :: Moisture
    Moisture = exp((DMC-223.9)/ (-41.7))+20.0
	return 
end subroutine


!def LawsonEq3(DMC):
!    '''National Standard and Best-fit Non-linear Regression Equations
!Linking DMC to Forest Floor Moisture Content in
!Coastal B.C., Southern Interior B.C. and Southern Yukon
!
!Southern Yukon - Pine/White Spruce
!Feather moss, Sphagnum and Undifferentiated duff (2-4 cm)^2
!
!LawsonEq3(8.5450511359999997)  = 442.82109267231488'''
!    return math.exp((DMC-157.3)/-24.6)+20

subroutine LawsonEq3(DMC, Moisture)
implicit none
	real, intent(in)    :: DMC 
	real, intent(out)   :: Moisture
    Moisture = exp((DMC-157.3)/ (-24.6))+20
	return 
end subroutine

!def LawsonEq4(DMC):
!    '''National Standard and Best-fit Non-linear Regression Equations
!Linking DMC to Forest Floor Moisture Content in
!Coastal B.C., Southern Interior B.C. and Southern Yukon
!
!Southern Yukon - Pine/White Spruce
!Reindeer lichen (2-4 cm)^2
!
!LawsonEq4(8.5450511359999997)  = 746.02210402093272'''
!    return math.exp((DMC-106.7)/-14.9)+20.0
    
subroutine LawsonEq4(DMC, Moisture)
implicit none
	real, intent(in)    :: DMC 
	real, intent(out)   :: Moisture
    Moisture = exp((DMC-106.7)/ (-14.9))+20.0
	return 
end subroutine


!def LawsonEq5(DMC):
!    '''National Standard and Best-fit Non-linear Regression Equations
!Linking DMC to Forest Floor Moisture Content in
!Coastal B.C., Southern Interior B.C. and Southern Yukon
!
!Southern Yukon - White Spruce
!White spruce/feather moss (2-4 cm)^2
!
!LawsonEq5(8.5450511359999997)  = 853.2397847094652'''
!
!    return math.exp((DMC-149.6)/-20.9)
   
subroutine LawsonEq5(DMC, Moisture)
implicit none
	real, intent(in)    :: DMC 
	real, intent(out)   :: Moisture
    Moisture = exp((DMC-149.6)/ (-20.9))
	return 
end subroutine    



subroutine daylengthN(dayOfYear, lat, length)
implicit none
	integer, intent(in)    :: dayOfYear
	real, intent(in)    :: lat 
	real, intent(out)   :: length
    real   :: declinationOfEarth, latInRad, pi , deg2Rad, rad2deg, hourAngle
    
    
    pi = acos(-1.0)
    deg2Rad = pi/180.
    rad2deg = 180./pi
    
  !  print *, pi, deg2Rad, rad2deg
!    """Computes the length of the day (the time between sunrise and
!    sunset) given the day of the year and latitude of the location.
!    Function uses the Brock model for the computations.
!    For more information see, for example,
!    Forsythe et al., "A model comparison for daylength as a
!    function of latitude and day of year", Ecological Modelling,
!    1995.
!    Parameters
!    ----------
!    dayOfYear : int
!        The day of the year. 1 corresponds to 1st of January
!        and 365 to 31st December (on a non-leap year).
!    lat : float
!        Latitude of the location in degrees. Positive values
!        for north and negative for south.
!    Returns
!    -------
!   d : float
!        Daylength in hours.
!    """
    latInRad = deg2rad * lat
 !   print *, latInRad
    declinationOfEarth = 23.45*sin(deg2rad *(360.0*(283.0+dayOfYear)/365.0))
  !  print *, declinationOfEarth
    
	!if (-tan(latInRad) * tan(deg2rad * (declinationOfEarth)) <= -1.0) then
	!	hourAngle = 24
	!else 
	!	if (tan(latInRad) * tan(deg2rad * (declinationOfEarth)) >= 1.0) then
	!		hourAngle = 0.0
	!	else
			hourAngle = rad2deg *(acos(-tan(latInRad) * tan(deg2rad * (declinationOfEarth))))
	!	endif
	!endif
	
	if (tan(latInRad) * tan(deg2rad * (declinationOfEarth)) <= -1.0) then
		hourAngle = 24 
		!print *, hourAngle	
	endif
	
    if (tan(latInRad) * tan(deg2rad * (declinationOfEarth)) >= 1.0) then
    	hourAngle = 0.0
    	!print *, hourAngle	
 	endif 
	
    !print *, hourAngle		
    length  = 2.0*hourAngle/15.0
        
        
  	return 
end subroutine       
end module fwiindices 

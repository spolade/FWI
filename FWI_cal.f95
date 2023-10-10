program FWI_cal_ERA_grid
use netcdf
use fwiindices
implicit none
INCLUDE "cdi.inc"

integer :: ncid_tas, ncid_pr, ncid_d2m, ncid_fwi, ncid_UWind, ncid_VWind, ncid_ct
integer :: tasVarId,prVarId,u10VarId, v10VarId, d2mVarId,lonDimID,latDimId
integer :: latVarId, lonVarId, timeVarId, varid_array, varid_ffmc, varid_dmc, varid_dc, varid_lat, varid_lon
integer :: varid_time,varid_EndDate, varid_StartDate 
integer :: dimid_lon, dimid_lat, dimid_time
integer :: NDIMS, month1, day1
integer :: rec,vdate, vtime
integer :: asize, streamID, vlistID, taxisID
integer :: numLons, numLats, numTimes, loni,latj
integer :: len_in, len_out, len_ct, status
integer, dimension(nf90_max_var_dims) :: dimIDs
integer, dimension(3) :: arrdims

integer, dimension(:), allocatable :: start, count, time, time_in

real, dimension(:), allocatable :: lon, lat


real, dimension(:, :), allocatable :: tas, wspd, pr, d2m, dew, sat, rhum, FWI_grid, isi_grid, bui_grid, ffmc_grid, dmc_grid, dc_grid 

real, dimension(:, :), allocatable :: ffmc0, dmc0, dc0, u10, v10  !!! initial or yesterday’s FFMC, DMC, and DC values

real     :: ffmc, dmc, dc, Lat_in1 !!! FFMC, DMC, and DC to be calculated
real     :: isi, bui, fwi, dsr     !!! ISI, BUI, FWI, and DSR to be calculated


character(len=300) :: in_path, out_path, ct_path, temp_name, pr_name, uwind_name, vwind_name, d2m_name, out_name, ct_name
character(len=300) :: file_ct, file_temp, file_pr, file_UWind, file_VWind, file_d2m, file_out, val_in, val_out, val_ct


! Read DIR path and data file names from Python Wrappers
read (*,*) in_path, out_path, temp_name, pr_name, uwind_name, vwind_name, d2m_name, out_name, ct_path, ct_name
 
call get_environment_variable (in_path, val_in, len_in, status, .true.)
call get_environment_variable (out_path, val_out, len_out, status, .true.)
call get_environment_variable (ct_path, val_ct, len_ct, status, .true.)
 
 
 
file_ct  = trim(val_ct(1:len_ct))//trim(ct_name)
    
file_temp  = trim(val_in(1:len_in))//trim(temp_name)
file_pr    = trim(val_in(1:len_in))//trim(pr_name)
file_UWind = trim(val_in(1:len_in))//trim(uwind_name)
file_VWind = trim(val_in(1:len_in))//trim(vwind_name)
file_d2m   = trim(val_in(1:len_in))//trim(d2m_name)
 
file_out   = trim(val_out(1:len_out))//trim(out_name) 

write (*,*) 'temp file   = ', file_temp
write (*,*) 'output file = ', file_out

status = nf90_open(file_temp,  nf90_NoWrite, ncid_tas)
print *, status
IF (status == 0) THEN


status = nf90_open(file_pr,    nf90_NoWrite, ncid_pr)
status = nf90_open(file_UWind, nf90_NoWrite, ncid_UWind)
status = nf90_open(file_VWind, nf90_NoWrite, ncid_VWind)
status = nf90_open(file_d2m,   nf90_NoWrite, ncid_d2m)

! Open file with constant 
status = nf90_open(file_ct,   nf90_NoWrite, ncid_ct)


! Create ouput file
status = nf90_create(file_out, NF90_NETCDF4, ncid_fwi)
call check(status, 'open')



! using CDI for reading NetCDF time to standard time
streamID = streamOpenRead(file_temp)




vlistID  = streamInqVlist(streamID)
taxisID  = vlistInqTaxis(vlistID)


! Define the variable IDs


status = nf90_inq_varid(ncid_tas, "lat", latVarId)
status = nf90_inq_varid(ncid_tas, "lon", lonVarId)
status = nf90_inq_varid(ncid_tas, "time", timeVarId)

status = nf90_inq_varid(ncid_tas, "2t", tasVarId)
status = nf90_inq_varid(ncid_pr, "tp", prVarId)
status = nf90_inq_varid(ncid_UWind, "10u", u10VarId)
status = nf90_inq_varid(ncid_VWind, "10v", v10VarId)
status = nf90_inq_varid(ncid_d2m, "2d", d2mVarId)



! Get the variable dimentions 

   status = nf90_Inquire_Variable(ncid_tas, tasVarId, dimids = dimIDs, ndims = NDIMS)

   status = nf90_Inquire_Dimension(ncid_tas, dimIDs(1), len = numLons)
   status = nf90_Inquire_Dimension(ncid_tas, dimIDs(2), len = numLats)
   status = nf90_Inquire_Dimension(ncid_tas, dimIDs(3), len = numTimes)
   
   
   print *, numTimes, numLats, numLons 
  
  
! Read the constant from File
status = nf90_inq_varid(ncid_ct, "ffmc", varid_ffmc)
status = nf90_inq_varid(ncid_ct, "dmc", varid_dmc)
status = nf90_inq_varid(ncid_ct, "dc", varid_dc)

allocate(ffmc0(numLons, numLats))
allocate(dmc0(numLons, numLats))
allocate(dc0(numLons, numLats))

 status = nf90_get_var(ncid_ct, varid_ffmc, ffmc0)
 status = nf90_get_var(ncid_ct, varid_dmc, dmc0)
 status = nf90_get_var(ncid_ct, varid_dc, dc0)
status = nf90_close(ncid_ct)


 ! Read one time step at a time; allocate memory  


   allocate(start(NDIMS))
   allocate(count(NDIMS))
   allocate(tas(numLons, numLats))
   allocate(u10(numLons, numLats))
   allocate(v10(numLons, numLats))   
   allocate(wspd(numLons, numLats)) 
   allocate(d2m(numLons, numLats))
   allocate(pr(numLons, numLats))
   allocate(rhum(numLons, numLats))
   allocate(dew(numLons, numLats))
   allocate(sat(numLons, numLats))
   allocate(lat(numLats))
   allocate(lon(numLons))
   allocate(time_in(numTimes))


   ! allocate memory for the output varibale Lat and Lon, time dimention will be written in file    
    allocate(FWI_grid(numLons, numLats))
    !FWI_grid = 0.0
    allocate(isi_grid(numLons, numLats))
    allocate(bui_grid(numLons, numLats))
    allocate(ffmc_grid(numLons, numLats))
    allocate(dmc_grid(numLons, numLats))
    allocate(dc_grid(numLons, numLats))


    
    
  ! Read coordinate varibles from input file to write in output 
status = nf90_get_var(ncid_tas, lonVarId, lon)
status = nf90_get_var(ncid_tas, latVarId, lat)
status = nf90_get_var(ncid_tas, timeVarId, time_in)
  
 
 
! File constant 
! define the dimensions name and output variable 
! Create constant file
status = nf90_create(file_ct, NF90_NETCDF4, ncid_ct)
call check(status, 'open')

status = nf90_def_dim(ncid_ct, 'longitude', numLons, dimid_lon)
status = nf90_def_dim(ncid_ct, 'latitude', numLats, dimid_lat)

 status = nf90_def_var(ncid_ct, 'longitude', NF90_FLOAT, [dimid_lon], varid_lon)
 status = nf90_def_var(ncid_ct, 'latitude', NF90_FLOAT, [dimid_lat], varid_lat)

 
 status = nf90_def_var(ncid_ct, 'ffmc', NF90_FLOAT, [dimid_lon, dimid_lat], varid_ffmc)
 status = nf90_put_att(ncid_ct, varid_ffmc, "long_name", "ffmc0")

 status = nf90_def_var(ncid_ct, 'dmc', NF90_FLOAT, [dimid_lon, dimid_lat], varid_dmc)
 status = nf90_put_att(ncid_ct, varid_dmc, "long_name", "dmc0")

 status = nf90_def_var(ncid_ct, 'dc', NF90_FLOAT, [dimid_lon, dimid_lat], varid_dc)
 status = nf90_put_att(ncid_ct, varid_dc, "long_name", "dc0")

! Write the coordinate vaible in output file 
status = nf90_put_var(ncid_ct, varid_lon, lon)
status = nf90_put_var(ncid_ct, varid_lat, lat)

 
 
       
! define the dimensions name and output variable 
status = nf90_def_dim(ncid_fwi, 'longitude', numLons, dimid_lon)
status = nf90_def_dim(ncid_fwi, 'latitude', numLats, dimid_lat)
status = nf90_def_dim(ncid_fwi, 'time', nf90_unlimited, dimid_time)


 status = nf90_def_var(ncid_fwi, 'longitude', NF90_FLOAT, [dimid_lon], varid_lon)
 status = nf90_def_var(ncid_fwi, 'latitude', NF90_FLOAT, [dimid_lat], varid_lat)
 status = nf90_def_var(ncid_fwi, 'time', NF90_FLOAT, [dimid_time], varid_time)
 
 status = nf90_def_var(ncid_fwi, 'fwi', NF90_FLOAT, [dimid_lon, dimid_lat, dimid_time], varid_array)
 status = nf90_put_att(ncid_fwi, varid_array, "long_name", "fire weather index")



! Write the coordinate vaible in output file 
status = nf90_put_var(ncid_fwi, varid_lon, lon)
status = nf90_put_var(ncid_fwi, varid_lat, lat)
status = nf90_put_var(ncid_fwi, varid_time, time_in)

! Write attributes: copy from input file
status = nf90_put_att(ncid_fwi, NF90_GLOBAL, 'Created', 'by :Suraj Polade for ClimDT project')
!status = nf90_put_att(ncid_fwi, varid_lat, 'units', 'degree_north')
status = nf90_put_att(ncid_fwi, varid_array, '_FillValue', -2e8)

status = nf90_copy_att(ncid_tas, timeVarId, 'standard_name',ncid_fwi, varid_time)
status = nf90_copy_att(ncid_tas, timeVarId, 'units',ncid_fwi, varid_time)
status = nf90_copy_att(ncid_tas, timeVarId, 'calendar',ncid_fwi, varid_time)
status = nf90_copy_att(ncid_tas, timeVarId, 'axis',ncid_fwi, varid_time)

status = nf90_copy_att(ncid_tas, latVarId, 'standard_name',ncid_fwi, varid_lat)
status = nf90_copy_att(ncid_tas, latVarId, 'long_name',ncid_fwi, varid_lat)
status = nf90_copy_att(ncid_tas, latVarId, 'units',ncid_fwi, varid_lat)
status = nf90_copy_att(ncid_tas, latVarId, 'axis',ncid_fwi, varid_lat)

status = nf90_copy_att(ncid_tas, lonVarId, 'standard_name',ncid_fwi, varid_lon)
status = nf90_copy_att(ncid_tas, lonVarId, 'long_name',ncid_fwi, varid_lon)
status = nf90_copy_att(ncid_tas, lonVarId, 'units',ncid_fwi, varid_lon)
status = nf90_copy_att(ncid_tas, lonVarId, 'axis',ncid_fwi, varid_lon)


! done defining
 status = nf90_enddef(ncid_fwi)

   start = (/ 1, 1, 1 /)
   count = (/numLons, numLats, 1 /)
   

  ! one record at a time, read data from input file and calculate 
  do rec = 1, numTimes
     start(3) = rec
     !print *, start, count

     ! Read the time using CDI  
      status = streamInqTimestep(streamID, (rec-1))

     ! Get the verification date
      vdate = taxisInqVdate(taxisID)
      !print *, vdate
     
     ! Extract the month values from date  
      month1 = MOD(vdate/100,100)
      !print *, month1
    ! Extract the day values from date  
      day1 = MOD(vdate/1,100)
      !print *, day1
      
 IF (month1 < 3) THEN

    !!! First z initialize FFMC, DMC, and DC
    ffmc0    = 0.00001 
    dmc0     = 0.00001   
    dc0      = 0.00001                 
    FWI_grid = 0.0
   !!! First z initialize FFMC, DMC, and DC
   !ffmc0    = 85.0                  
   !dmc0     = 6.0                  
   !dc0      = 15.0     
 end if


  if (month1 >= 3) then

     status = nf90_get_var(ncid_tas, tasVarId, tas, start, count)
     status = nf90_get_var(ncid_pr, prVarId, pr, start, count)
     status = nf90_get_var(ncid_UWind, u10VarId, u10, start, count)
     status = nf90_get_var(ncid_VWind, v10VarId, v10, start, count)
     status = nf90_get_var(ncid_d2m, d2mVarId, d2m, start, count)

      

    ! Calculate wind speed from U and V component

      wspd = sqrt((u10 * u10) + (v10 * v10))
       
      !clear  u10, v10
       
    ! Calculate relative humidity from deqpoint temperature
     dew=6.11*exp(53.49-6808/d2m-5.09*log(d2m))
     sat=6.11*exp(53.49-6808/tas-5.09*log(tas))
     
     rhum = (dew/sat) *100 


    ! Correct the units 

    ! Temperature K == Deg C
     tas = tas - 273.15
  
    ! Wind speed m/s == km/h
      wspd = wspd * 3.6   
       
     !$OMP PARALLEL DO    
      do loni = 0, numLons
        do latj = 0, numLats
  
       if (rhum(loni,latj)>= 100.0) rhum(loni,latj)=100.0                    !!! make sure RHUM <= 100.0 
       call FFMCcalc(tas(loni,latj),rhum(loni,latj),wspd(loni,latj),pr(loni,latj),ffmc0(loni,latj),ffmc)   !!! calculate FFMC
       call DMCcalc(tas(loni,latj),rhum(loni,latj),pr(loni,latj),dmc0(loni,latj),vdate,lat(latj),dmc)       !!! calcualte DMC 
       call DCcalc (tas(loni,latj),pr(loni,latj),dc0(loni,latj),month1,lat(latj),dc)              !!! calcualte DC
       call ISIcalc (ffmc,wspd(loni,latj),isi)                    !!! calcualte ISI 
       call BUIcalc (dmc,dc,bui)                       !!! calculate BLL
       call FWIcalc (isi,bui,fwi)                      !!! calculate FWI
 
       FWI_grid(loni,latj)  = fwi  
       isi_grid(loni,latj)  = isi
       bui_grid(loni,latj)  = bui
       ffmc_grid(loni,latj) = ffmc 
       dmc_grid(loni,latj)  = dmc
       dc_grid(loni,latj)   = dc

      !print *, month1, tas(loni,latj), rhum(loni,latj), wspd(loni,latj), pr(loni,latj), FWI_grid(loni,latj)      
      end do ! numLons loop
     end do  ! numLats loop
     !$OMP END PARALLEL DO

       if (isnan(ffmc)) then
       else
        ffmc0 = ffmc_grid
       end if

       if (isnan(dmc)) then
       else
        dmc0 = dmc_grid
       end if


       if (isnan(dc)) then
       else
        dc0 = dc_grid
       end if

 

  end if ! if month > 3     
      status = nf90_put_var(ncid_fwi, varid_array, FWI_grid, start, count)      
      status = nf90_put_var(ncid_ct, varid_ffmc, ffmc_grid, start, count)
      status = nf90_put_var(ncid_ct, varid_dmc, dmc_grid, start, count)
      status = nf90_put_var(ncid_ct, varid_dc, dc_grid, start, count)

end do ! numTime loop



  status = nf90_close(ncid_tas)
  status = nf90_close(ncid_pr)
  status = nf90_close(ncid_UWind)
  status = nf90_close(ncid_VWind)
  status = nf90_close(ncid_d2m)
  status = nf90_close(ncid_ct)
  status = nf90_close(ncid_fwi)
call check(status, 'close')


   print *, numLons, numLats, numTimes
   print *, file_ct
 end if
  contains

subroutine check(status, operation)
    use netcdf
    implicit none
    integer, intent(in) :: status
    character(len=*), intent(in) :: operation
    if (status == NF90_NOERR) return
    print *, "Error encountered during ", operation
    print *, nf90_strerror(status)
    STOP 1
end subroutine check
  
end program FWI_cal_ERA_grid


!                    ****************
                     SUBROUTINE METEO
!                    ****************
!
     &(PATMOS,WINDX,WINDY,FUAIR,FVAIR,X,Y,AT,LT,NPOIN,VENT,ATMOS,
     & HN,TRA01,GRAV,ROEAU,NORD,PRIVE,FO1,FILES,LISTIN)
!
!***********************************************************************
! TELEMAC2D   V6P3                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES ATMOSPHERIC PRESSURE AND WIND VELOCITY FIELDS
!+               (IN GENERAL FROM INPUT DATA FILES).
!
!warning  CAN BE ADAPTED BY USER
!
!history  J-M HERVOUET (LNHE)
!+        02/01/2004
!+        V5P4
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        30/01/2013
!+        V6P3
!+   Now 2 options with an example for reading a file. 
!
!history  EHSAN SARHADI ZADEH (AnteaGroup, BELGIUM)
!+        17/01/2014
!+        V6P3
!+
!+   Variable atmospheric pressure in time & space
!+   Variable wind velocity fields in time & space
!+   Interpolation in time & space using 'FASP & INTERPOLATION' subroutine
!+   Variable input files in netCDF 
!+   Dataset type: NetCDF-3/CDM
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AT,LT          |-->| TIME, ITERATION NUMBER
!| ATMOS          |-->| YES IF PRESSURE TAKEN INTO ACCOUNT
!| FUAIR          |-->| VELOCITY OF WIND ALONG X, IF CONSTANT
!| FVAIR          |-->| VELOCITY OF WIND ALONG Y, IF CONSTANT
!| GRAV           |-->| GRAVITY ACCELERATION
!| HN             |-->| DEPTH
!| NORD           |-->| DIRECTION OF NORTH, COUNTER-CLOCK-WISE
!|                |   | STARTING FROM VERTICAL AXIS
!| NPOIN          |-->| NUMBER OF POINTS IN THE MESH
!| PATMOS         |<--| ATMOSPHERIC PRESSURE
!| PRIVE          |-->| USER WORKING ARRAYS (BIEF_OBJ BLOCK)
!| ROEAU          |-->| WATER DENSITY
!| TRA01          |-->| WORKING ARRAY
!| VENT           |-->| YES IF WIND TAKEN INTO ACCOUNT
!| WINDX          |<--| FIRST COMPONENT OF WIND VELOCITY
!| WINDY          |<--| SECOND COMPONENT OF WIND VELOCITY
!| X              |-->| ABSCISSAE OF POINTS
!| Y              |-->| ORDINATES OF POINTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D, ONLY : MARDAT,MARTIM
!      USE DECLARATIONS_TELEMAC2D, ONLY : T2DFO1,T2D_FILES,LISTIN
      USE netcdf
!      USE omp_lib
!
      IMPLICIT NONE
!      INCLUDE 'netcdf.inc'
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: LT,NPOIN,FO1
      LOGICAL, INTENT(IN)             :: ATMOS,VENT,LISTIN
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),HN(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: WINDX(NPOIN),WINDY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: PATMOS(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: TRA01(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: FUAIR,FVAIR,AT,GRAV,ROEAU,NORD
!      DOUBLE PRECISION, INTENT(IN)    :: DDC
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: PRIVE
      TYPE(BIEF_FILE), INTENT(IN)     :: FILES(*)
!      INTEGER          , INTENT(IN)    :: MARDAT(3),MARTIM(3)
      INTEGER YEAR,MONTH,DAY,HOUR,MINUTE,SEC,TS,STARTCDF
      DOUBLE PRECISION SHIFTTIME
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER MY_OPTION,UL,INPUT_FILE
      DOUBLE PRECISION P0,Z(1),FUAIR1,FUAIR2,FVAIR1,FVAIR2,COEF
      DOUBLE PRECISION UAIR,VAIR,AT1,AT2
      DOUBLE PRECISION, DIMENSION(NPOIN)   :: TRA00,TRA02,TRA03,TRA04
      DOUBLE PRECISION, DIMENSION(NPOIN)   :: PATMOS1,PATMOS2
!      DOUBLE PRECISION, DIMENSION(NPOIN)   :: PATMOS3
! 
!         BILINEAR INTERPOLATION 
!
!      DOUBLE PRECISION, DIMENSION(NPOIN)   :: L1,L2,L3,L4
!      INTEGER, DIMENSION(NPOIN)            :: NP1,NP2,NP3,NP4
!
!         BARYCENTRIC INTERPOLATION    
! 
      DOUBLE PRECISION, DIMENSION(NPOIN)   :: L1,L2
      INTEGER, DIMENSION(NPOIN)            :: NP1,NP2,NP3
      INTEGER, DIMENSION(NPOIN)            :: NP4
      
!
!-----------------------------------------------------------------------
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     DEFINITION OF ARRAYS FOR NetCDF INPUT FILE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER TIMEin,XCELLin,YCELLin,NPin
      DOUBLE PRECISION,ALLOCATABLE      :: Pin1(:),Pin2(:)
      DOUBLE PRECISION,ALLOCATABLE      :: Uin1(:),Uin2(:)
      DOUBLE PRECISION,ALLOCATABLE      :: Vin1(:),Vin2(:)
      DOUBLE PRECISION,ALLOCATABLE      :: LONGin(:),LATin(:)
      DOUBLE PRECISION,ALLOCATABLE      :: lon0(:),lat0(:),Xin(:),Yin(:)
      DOUBLE PRECISION                     pi,Radius
!      INTEGER, INTENT(IN) :: NPTFR,NBOR(NPTFR),KP1BOR(NPTFR)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     NAME OF THE NetCDF FILE 

      CHARACTER(len=*),PARAMETER::FILE_NAME="T2DBI2"
!      CHARACTER(len=*),PARAMETER::FILE_NAME="../HIRLAM.nc"
      INTEGER ncid,status     

      INTEGER NDIMS, NRECS , NLATS, NLONS

      CHARACTER (len = *), PARAMETER :: LAT_NAME = "y"
      CHARACTER (len = *), PARAMETER :: LON_NAME = "x"
      CHARACTER (len = *), PARAMETER :: REC_NAME = "time"
      INTEGER lon_dimid, lat_dimid, rec_dimid

!     TELL THE NetCDF LIBRARY WHERE TO READ DATA
      
      INTEGER,ALLOCATABLE      :: start(:), count(:), TIME_in(:)

!     COORDINATE VARIABLES

      DOUBLE PRECISION,ALLOCATABLE      ::     lats(:), lons(:)
      INTEGER lon_varid, lat_varid, time_varid 
      INTEGER ndims_in, nvars_in, ngatts_in, unlimdimid_in

!     VARIABLES NAME

      CHARACTER (len = *), PARAMETER :: WIND_V_NAME="wind_v"
      CHARACTER (len = *), PARAMETER :: WIND_U_NAME="wind_u"
      CHARACTER (len = *), PARAMETER :: PRES_NAME="p"
      CHARACTER (len = *), PARAMETER :: WIND_ABS_NAME="wind_uv_abs"

!     VARIABLES ID
            
      INTEGER v_varid, u_varid, p_varid, abs_varid
      INTEGER,ALLOCATABLE      :: dimids(:)

!     VARIABLE "UNITS" ATTRIBUTE

      CHARACTER (len = *), PARAMETER :: UNITS = "units"
      CHARACTER (len = *), PARAMETER :: WIND_V_UNITS = "m/s"
      CHARACTER (len = *), PARAMETER :: WIND_U_UNITS = "m/s"
      CHARACTER (len = *), PARAMETER :: WIND_ABS_UNITS = "m/s"
      CHARACTER (len = *), PARAMETER :: PRES_UNITS = "Pa"
      CHARACTER (len = *), PARAMETER :: LAT_UNITS = "degrees_north"
      CHARACTER (len = *), PARAMETER :: LON_UNITS = "degrees_east"
      CHARACTER (len = *), PARAMETER :: TIME_UNITS = "minutes 
     & since 1970-01-01 00:00:00 +00:00"

!     ARRAYS TO HOLD ONE TIMESTEP OF DATA

      DOUBLE PRECISION,ALLOCATABLE      ::   WIND_V_in(:,:)
      DOUBLE PRECISION,ALLOCATABLE      ::   WIND_U_in(:,:)
      DOUBLE PRECISION,ALLOCATABLE      ::   PRES_in(:,:)
!      DOUBLE PRECISION,ALLOCATABLE      ::   WIND_ABS_in(:,:)


!     LOOP INDICES
      INTEGER lat, lon, rec, I, J, K

!     CHECK THE UNITS ATTRIBUTES
      CHARACTER*80 WIND_V_UNITS_in, WIND_U_UNITS_in
      CHARACTER*80 WIND_ABS_UNITS_in, PRES_UNITS_in
      CHARACTER*80 LAT_UNITS_in, LON_UNITS_in, TIME_UNITS_in
      INTEGER att_len
!     ....
!     ....
!     ....

!-----------------------------------------------------------------------
!
!     DATA THAT YOU DECLARE AND READ HERE ONCE IN A FILE MAY HAVE TO BE
!     KEPT BECAUSE THIS SUBROUTINE IS CALLED AT EVERY TIME STEP.
!     WITHOUT THE SAVE COMMAND, ALL LOCAL DATA ARE FORGOTTEN IN THE NEXT
!     CALL.
!
      SAVE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!     CHOOSE YOUR OPTION !!!
!
!     1: CONSTANTS GIVEN BY THE KEYWORDS:
!        AIR PRESSURE (GIVEN HERE AS P0, NO KEYWORD)
!        WIND VELOCITY ALONG X (HERE FUAIR)
!        WIND VELOCITY ALONG Y (HERE FVAIR)
!        THEY WILL BE SET ONCE FOR ALL BEFORE THE FIRST ITERATION (LT=0)
!
!     2: WIND COMPONENTS & AIR PRESSURE VARIABLE IN SPACE & TIME 
!        IN THE FILE FO1_WIND DECLARED AS: 
!        FORMATTED DATA FILE 1             : HIRLAM.txt
!        OR
!        BINARY DATA FILE 2                : HIRLAM.nc
!        INTERPOLATION THE DATA FROM INPUT GRID TO MESH 
!
      MY_OPTION = 2
!
!-----------------------------------------------------------------------
!
!     CHOOSE YOUR OPTION !!!
!
!     1: INTERPOLATION THE DATA FROM ASCII INPUT FILE
!     2: INTERPOLATION THE DATA FROM BINARY INPUT FILE (NetCDF)
!
      INPUT_FILE = 2
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      IF(INPUT_FILE.EQ.2) THEN
!       
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     OPENING NetCDF INPUT FILE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!-----------------------------------------------------------------------
!
!     BEWARE, ONLY ONE TIME CALLING NetCDF UP AT FIRST TIMESTEP IS NEEDED
!
      IF(LT.EQ.0) THEN
!
!-----------------------------------------------------------------------
!
!     OPEN THE NetCDF FILE
!
!-----------------------------------------------------------------------
!
      status = NF90_OPEN(FILE_NAME, nf90_nowrite, ncid)

!     INQUIRY FUNCTIONS (NF90_INQ) FOR UNKNOWN NetCDF FILE
!     NF90_INQ TELLS HOW MANY VARIABLES AND DIMENSIONS ARE IN THE FILE
!       -VARIABLES
!       -DIMENSIONS
!       -GLOBAL ATTRIBUTES
!       -DIMENSION ID OF THE UNLIMITED DIMENSION
      status = NF90_INQUIRE(ncid, ndims_in, nvars_in, 
     & ngatts_in, unlimdimid_in)  

      NDIMS=ndims_in

      WRITE(LU,*) 'NDIMS' , ndims_in
!      IF (ndims_in==3) STOP
!
!     CHECK THE CONTENTS OF THE NetCDF FILE  
!     HOW MANY VALUES OF "LAT, LONG, & TIME" ARE THERE 
!     CHANGE THE SIZE OF ARRAYS 
!     

      status = nf90_inq_dimid(ncid, "time", rec_dimid)
      status = nf90_inquire_dimension(ncid, rec_dimid, len = NRECS)
      WRITE(LU,*) 'NRECS' , NRECS
!      IF (NRECS==249) STOP

      status = nf90_inq_dimid(ncid, "y", lat_dimid)
      status = nf90_inquire_dimension(ncid, lat_dimid, len = NLATS)
      WRITE(LU,*) 'NLATS' , NLATS
!      IF (NLATS==173) STOP

      status = nf90_inq_dimid(ncid, "x", lon_dimid)
      status = nf90_inquire_dimension(ncid, lon_dimid, len = NLONS)
      WRITE(LU,*) 'NLONS' , NLONS
!      IF (NLONS==201) STOP

!
!-----------------------------------------------------------------------
!
!     PRE-DIFINITION OF ARRAYS SIZE, IF THE SIZE IS NOT DEFINED IN NetCDF FILE 

!     READING 3D DATA, 173 X 201 LAT-LON GRID, WITH 249 TIMESTEPS

!      NDIMS = 3
!      NRECS = 249
!      NLATS = 173
!      NLONS = 201

      NP4=100.5

      ALLOCATE (start(NDIMS), count(NDIMS), dimids(NDIMS))
      ALLOCATE (TIME_in(NRECS))
      ALLOCATE (WIND_V_in(NLONS, NLATS))
      ALLOCATE (WIND_U_in(NLONS, NLATS))
      ALLOCATE (PRES_in(NLONS, NLATS))
!      ALLOCATE (WIND_ABS_in(NLONS, NLATS))
      ALLOCATE (lats(NLATS), lons(NLONS))
   
!
!-----------------------------------------------------------------------
!
!     ID OF LATITUDE AND LONGITUDE COORDINATE 
!
!-----------------------------------------------------------------------
!
      status = NF90_INQ_VARID(ncid, LAT_NAME, lat_varid)
      status = NF90_INQ_VARID(ncid, LON_NAME, lon_varid)
      status = NF90_INQ_VARID(ncid, REC_NAME, time_varid)

!
!-----------------------------------------------------------------------
!
!      CHECK THE "UNITS" ATTRIBUTE OF NetCDF VARIABLES 
!
!-----------------------------------------------------------------------
!
      status = NF90_GET_ATT(ncid, time_varid, UNITS, TIME_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,time_varid,UNITS,len=att_len)
!      IF (TIME_UNITS_in(1:att_len) /= TIME_UNITS) STOP 'TIME UNIT ERROR DETECTED !'
      IF (TIME_UNITS_in(1:att_len) /= TIME_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE TEMPS ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'TIME UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',TIME_UNITS
      WRITE(LU,*) ' netCDF Unit: ',TIME_UNITS_in(1:att_len)

      ENDIF

      status = NF90_GET_ATT(ncid, lat_varid, UNITS, LAT_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,lat_varid,UNITS,len=att_len)
!      IF (LAT_UNITS_in(1:att_len) /= LAT_UNITS) STOP 'LATITUDE UNIT ERROR DETECTED !'
      IF (LAT_UNITS_in(1:att_len) /= LAT_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE LATITUDE ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'LATITUDE UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',LAT_UNITS
      WRITE(LU,*) ' netCDF Unit: ',LAT_UNITS_in(1:att_len)

      ENDIF

      status = NF90_GET_ATT(ncid, lon_varid, UNITS, LON_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,lon_varid,UNITS,len=att_len)
!      IF (LON_UNITS_in(1:att_len) /= LON_UNITS) STOP 'LONGITUDE UNIT ERROR DETECTED !'      
      IF (LON_UNITS_in(1:att_len) /= LON_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE LONGITUDE ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'LONGITUDE UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',LON_UNITS
      WRITE(LU,*) ' netCDF Unit: ',LON_UNITS_in(1:att_len)

      ENDIF    
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     READ THE LATITUDE AND LONGITUDE COORDINATE
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      status = NF90_GET_VAR(ncid, lat_varid, lats)
      status = NF90_GET_VAR(ncid, lon_varid, lons)
      status = NF90_GET_VAR(ncid,time_varid,TIME_in)

          NPin=NLONS*NLATS

          ALLOCATE (Pin1(NPin),Pin2(NPin))
          ALLOCATE (Uin1(NPin),Uin2(NPin))
          ALLOCATE (Vin1(NPin),Vin2(NPin))
          ALLOCATE (LONGin(NPin),LATin(NPin))
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CONVERT THE NetCDF COORDINATE STYLE (long,lat) TO TELEMAC STYLE (NPOIN)
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!          
!          OPEN (50, FILE="temp_coordinate.txt")

!          DO J=1,NLATS
!          WRITE(50,*)    lons
!          END DO 

!          DO J=1,NLATS
!          lons=lats(J)
!          WRITE(50,*)    lons
!          END DO 

!          CLOSE (50)
!          OPEN (150, FILE="temp_coordinate.txt")

!          READ (150,*) LONGin, LATin
!          CLOSE (150)

      DO J=1,NLATS

          LONGin (1+NLONS*(J-1) : NLONS*(J)) = lons (1 : NLONS)
          
      END DO
          
      DO J=1,NLATS

          LATin (1+NLONS*(J-1) : NLONS*(J)) = lats (J)
          
      END DO

!      WRITE(LU,*) 
!      WRITE(LU,*) 'READ netCDF WIND & PRESSURE FILE'
!      WRITE(LU,*) '.nc ID', ncid
!      WRITE(LU,*) LON_NAME, LONGin
!      WRITE(LU,*) LAT_NAME, LATin

!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!         min and max of MERCATOR coordinates
!         φmax = ±85.05113° 
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
        IF(MINVAL(LATin).LT.-85.05113.OR.MAXVAL(LATin).GT.85.05113) THEN
!          IF(LNG.EQ.1) WRITE(LU,*) 'POINT OUTSIDE OF THE DOMAIN'
!          IF(LNG.EQ.1) WRITE(LU,*) 'SEGMENT INTERSECT DOMAIN BOUNDARY'
!          IF(LNG.EQ.2) WRITE(LU,*) 'POINT OUTSIDE OF THE DOMAIN'
!          IF(LNG.EQ.2) WRITE(LU,*) 'SEGMENT INTERSECT DOMAIN BOUNDARY'  
          IF(LNG.EQ.1) WRITE(LU,*) 'LATITUDE OUT OF RANGE'
          IF(LNG.EQ.1) WRITE(LU,*) 'φmax = ±85.05113°'
          IF(LNG.EQ.2) WRITE(LU,*) 'LATITUDE OUT OF RANGE'
          IF(LNG.EQ.2) WRITE(LU,*) 'φmax = ±85.05113°'
        END IF           
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!         function [X, Y] = lonlat2mercator_telemac(lon, lat)
!         % author : Olivier Gourgue
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
           ALLOCATE (lon0(NPin),lat0(NPin),Xin(NPin),Yin(NPin)) 

!         % lon and lat are in degrees (lonlat degrees WGS 84 - EPSG code 4326)
    
!         % "mercator_telemac" is the coordinate system required by Telemac if we
!         % want to take into account the sphericity of the domain
!         % It is not an official coordinate system (no EPSG code) but the
!         % transformation formulae from lonlat degrees WGS 84 are provided
    
!         % earth radius
           Radius = 6371000.D0
!         % number π = 3.1415926535897932384626433832795 
!           pi = acos(-1.)
           pi = 4.*atan(1.)
!         % arbitrary reference meridian (0° is Greenwich meridian)
           lon0 = 0.D0*pi/180.D0
           LONGin=LONGin*pi/180.D0
!         % arbitrary reference parallel
!         % '°' *pi/180 = 'Radians'
           lat0 = 50.D0*pi/360.D0+pi/4.D0
           LATin=LATin*pi/360.D0+pi/4.D0

           Xin=Radius*(LONGin-lon0)
           Yin=Radius*(log(tan(LATin))-log(tan(lat0)))
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CHECK THE SIMULATION DATE USING MARDAT AND MARTIM ARRAYS
!

       WRITE(LU,*) 'YEAR',MARDAT(1),'MONTH',MARDAT(2),'DAY',MARDAT(3)

       WRITE(LU,*) 'HOUR',MARTIM(1),'MINUTE',MARTIM(2),'SEC',MARTIM(3)

!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CONVERT SIMULATION DATE TO MJD (MODIFIED JULIAN DAYS SINCE 1858)     
! 
       CALL DATE_MJD( MARDAT(2),MARDAT(3),MARDAT(1),TS ) 
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CONVERT NetCDF DATE TO MJD (MODIFIED JULIAN DAYS SINCE 1858)       
! 
       CALL DATE_MJD( 1,1,1970,STARTCDF )
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CHECK THE NetCDF DATE 
!
       SHIFTTIME = 24*60*60*(TS - STARTCDF
!     &       + NetCDF(hr)/24.D0+NetCDF(min)/1440.D0+NetCDF(sec)/86400.D0
     &       + MARTIM(1)/24.D0+MARTIM(2)/1440.D0+MARTIM(3)/86400.D0)

            WRITE(LU,*) 'SHIFTTIME',SHIFTTIME


        IF (SHIFTTIME.LT.((TIME_in(1))*60)) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'NOMBRE DE FICHIERS NetCDF', NRECS 
      IF(LNG.EQ.1) WRITE(LU,*) 'DATE HORS netCDF LIMITE << netCDF'
      IF(LNG.EQ.2) WRITE(LU,*) 'NUMBER OF netCDF FILES', NRECS
      IF(LNG.EQ.2) WRITE(LU,*) 'DATE OUT OF netCDF RANGE << netCDF'

        ELSE IF (SHIFTTIME.GT.((TIME_in(NRECS))*60)) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'NOMBRE DE FICHIERS NetCDF', NRECS
      IF(LNG.EQ.1) WRITE(LU,*) 'DATE HORS netCDF LIMITE netCDF <<'
      IF(LNG.EQ.2) WRITE(LU,*) 'NUMBER OF netCDF FILES', NRECS
      IF(LNG.EQ.2) WRITE(LU,*) 'DATE OUT OF netCDF RANGE netCDF <<'

        END IF
!
!-----------------------------------------------------------------------
!              
!              'CALLING SUBROUTINE INTERPOLATION'
!              'FINDING WEIGHTING FACTORS (L1,L2,L3,L4) FOR'
!              'INTERPOLATING THE WIND & PRESSUER IN SPACE FROM FILE TO MESH'
!              'SINCE COMPUTATIONAL MESH & INPUT FILE DOMAIN ARE CONSTANT'
!              'BEWARE, ONLY ONE CALCULATION AT FIRST TIMESTEP IS NEEDED'
!
!              Wx=(L1*Wx1(NP1))+(L2*Wx2(NP2))+(L3*Wx3(NP3))+(L4*Wx4(NP4))
!              Wy=(L1*Wy1(NP1))+(L2*Wy2(NP2))+(L3*Wy3(NP3))+(L4*Wy4(NP4))
!              P =(L1*P(NP1))+(L2*P(NP2))+(L3*P(NP3))+(L4*P(NP4))
!
!-----------------------------------------------------------------------
! 
! 
!         BILINEAR INTERPOLATION 
!
!          CALL INTERPOLATION (X,Y,NPOIN,Xin,Yin,NPin,
!     &        L1,L2,L3,L4,NP1,NP2,NP3,NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
          CALL INTERPOLATION (X,Y,NPOIN,Xin,Yin,NPin,L1,L2,NP1,NP2,NP3)
 
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     ID OF NetCDF VARIABLES
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      status = NF90_INQ_VARID(ncid, WIND_V_NAME, v_varid)
      status = NF90_INQ_VARID(ncid, WIND_U_NAME, u_varid)
      status = NF90_INQ_VARID(ncid, PRES_NAME, p_varid)
!      status = NF90_INQ_VARID(ncid, WIND_ABS_NAME, abs_varid)

!     1ST RECORD OF NLATS*NLONS VALUES
!     STARTING AT BEGINNING OF THE RECORD (1, 1, REC) 

      count = (/ NLONS, NLATS, 1 /)
      start = (/ 1, 1, 1 /)
!      count(1) = NLONS
!      count(2) = NLATS
!      count(3) = 1
!      start(1) = 1
!      start(2) = 1

!
!-----------------------------------------------------------------------
!
!      CHECK THE "UNITS" ATTRIBUTE OF NetCDF VARIABLES 
!
!-----------------------------------------------------------------------
!      

      status = NF90_GET_ATT(ncid, u_varid, UNITS, WIND_U_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,u_varid,UNITS,len=att_len)
!      IF (WIND_U_UNITS_in(1:att_len) /= WIND_U_UNITS) STOP 'U UNIT ERROR DETECTED !' 
      IF (WIND_U_UNITS_in(1:att_len) /= WIND_U_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE VENT_U ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND_U UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',WIND_U_UNITS
      WRITE(LU,*) ' netCDF Unit: ',WIND_U_UNITS_in(1:att_len)

      ENDIF 

      status = NF90_GET_ATT(ncid, v_varid, UNITS, WIND_V_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,v_varid,UNITS,len=att_len)
!      IF (WIND_V_UNITS_in(1:att_len) /= WIND_V_UNITS) STOP 'V UNIT ERROR DETECTED !' 
      IF (WIND_V_UNITS_in(1:att_len) /= WIND_V_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE VENT_V ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND_V UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',WIND_V_UNITS
      WRITE(LU,*) ' netCDF Unit: ',WIND_V_UNITS_in(1:att_len)

      ENDIF

!      status = NF90_GET_ATT(ncid, abs_varid, UNITS, WIND_ABS_UNITS_in)
!      status = NF90_INQUIRE_ATTRIBUTE(ncid,abs_varid,UNITS,len=att_len)
!      IF (WIND_ABS_UNITS_in(1:att_len) /= WIND_ABS_UNITS) STOP 'ABS UNIT ERROR DETECTED !' 
!      IF (WIND_ABS_UNITS_in(1:att_len) /= WIND_ABS_UNITS) THEN

!      WRITE(LU,*) ' '
!      WRITE(LU,*) 'METEO'
!      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE VENT_ABS ERREUR DÉTECTÉE !'
!      IF(LNG.EQ.2) WRITE(LU,*) 'WIND_ABS UNIT ERROR DETECTED !'
!      WRITE(LU,*) 'TELEMAC Unit: ',WIND_ABS_UNITS
!      WRITE(LU,*) ' netCDF Unit: ',WIND_ABS_UNITS_in(1:att_len)

!      ENDIF

      status = NF90_GET_ATT(ncid, p_varid, UNITS, PRES_UNITS_in)
      status = NF90_INQUIRE_ATTRIBUTE(ncid,p_varid,UNITS,len=att_len)
!      IF (PRES_UNITS_in(1:att_len) /= PRES_UNITS) STOP 'PRESSURE UNIT ERROR DETECTED !' 
      IF (PRES_UNITS_in(1:att_len) /= PRES_UNITS) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'UNITÉ DE PRESSION ERREUR DÉTECTÉE !'
      IF(LNG.EQ.2) WRITE(LU,*) 'PRESSURE UNIT ERROR DETECTED !'
      WRITE(LU,*) 'TELEMAC Unit: ',PRES_UNITS
      WRITE(LU,*) ' netCDF Unit: ',PRES_UNITS_in(1:att_len)

      ENDIF 

!
!-----------------------------------------------------------------------
!         READING THE FIRST TWO SETS OF DATA FROM NetCDF FILE
!-----------------------------------------------------------------------
!
      rec = 1
      start(3) = rec

      status = NF90_GET_VAR(ncid,u_varid,WIND_U_in,start,count)
      status = NF90_GET_VAR(ncid,v_varid,WIND_V_in,start,count)
      status = NF90_GET_VAR(ncid,p_varid,PRES_in,start,count)
!      status = NF90_GET_VAR(ncid,abs_varid,WIND_ABS_in,start,count)
      
      AT1 = (TIME_in(rec))*60
!
!-----------------------------------------------------------------------
!     LOADING & CHECKING THE VARIABLES   
!-----------------------------------------------------------------------
!   

      Uin1=(/WIND_U_in(:,:)/)
      Vin1=(/WIND_V_in(:,:)/)
      Pin1=(/PRES_in(:,:)/)

            WRITE(LU,*) REC_NAME, rec
            WRITE(LU,*) WIND_U_NAME, Uin1(1:3)
            WRITE(LU,*) WIND_V_NAME, Vin1(1:3)
            WRITE(LU,*) PRES_NAME, Pin1(1:3)
            WRITE(LU,*) REC_NAME, TIME_in(1:3)
      
              IF(MAXVAL(Uin1)>200.D0) GO TO 500
              IF(MAXVAL(Vin1)>200.D0) GO TO 500
              IF(MINVAL(Uin1)<-200.D0) GO TO 500
              IF(MINVAL(Vin1)<-200.D0) GO TO 500           
              IF(MAXVAL(Pin1)>110000.D0) GO TO 500              
              IF(MAXVAL(Pin1)<90000.D0) GO TO 500
!
!-----------------------------------------------------------------------
!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
!  
!              CALL FASP(X,Y,TRA00,NPOIN,Xin,Yin,Uin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)                
!              CALL FASP(X,Y,TRA02,NPOIN,Xin,Yin,Vin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS1,NPOIN,Xin,Yin,Pin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA00=L1*Uin1(NP1)+L2*Uin1(NP2)+L3*Uin1(NP3)+L4*Uin1(NP4)
!      TRA02=L1*Vin1(NP1)+L2*Vin1(NP2)+L3*Vin1(NP3)+L4*Vin1(NP4)
!      PATMOS1=L1*Pin1(NP1)+L2*Pin1(NP2)+L3*Pin1(NP3)+L4*Pin1(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA00=L1*Uin1(NP1)+L2*Uin1(NP2)+(1-L1-L2)*Uin1(NP3)
      TRA02=L1*Vin1(NP1)+L2*Vin1(NP2)+(1-L1-L2)*Vin1(NP3)
      PATMOS1=L1*Pin1(NP1)+L2*Pin1(NP2)+(1-L1-L2)*Pin1(NP3)
         
!
!-----------------------------------------------------------------------
!
      rec = rec + 1
      start(3) = rec

      status = NF90_GET_VAR(ncid,u_varid,WIND_U_in,start,count)
      status = NF90_GET_VAR(ncid,v_varid,WIND_V_in,start,count)
      status = NF90_GET_VAR(ncid,p_varid,PRES_in,start,count)
!      status = NF90_GET_VAR(ncid,abs_varid,WIND_ABS_in,start,count)

      AT2 = (TIME_in(rec))*60
!
!-----------------------------------------------------------------------
!     LOADING & CHECKING THE VARIABLES       
!-----------------------------------------------------------------------
! 
      Uin2=(/WIND_U_in(:,:)/)
      Vin2=(/WIND_V_in(:,:)/)
      Pin2=(/PRES_in(:,:)/)


            WRITE(LU,*) REC_NAME, rec
            WRITE(LU,*) WIND_U_NAME, Uin2(1:3)
            WRITE(LU,*) WIND_V_NAME, Vin2(1:3)
            WRITE(LU,*) PRES_NAME, Pin2(1:3)
            WRITE(LU,*) REC_NAME, TIME_in(1:3)

              IF(MAXVAL(Uin2)>200.D0) GO TO 300
              IF(MAXVAL(Vin2)>200.D0) GO TO 300
              IF(MINVAL(Uin2)<-200.D0) GO TO 300
              IF(MINVAL(Vin2)<-200.D0) GO TO 300
              IF(MAXVAL(Pin2)>110000.D0) GO TO 300
              IF(MAXVAL(Pin2)<90000.D0) GO TO 300
!
!-----------------------------------------------------------------------
!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!              CALL FASP(X,Y,TRA03,NPOIN,Xin,Yin,Uin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,TRA04,NPOIN,Xin,Yin,Vin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS2,NPOIN,Xin,Yin,Pin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+L3*Uin2(NP3)+L4*Uin2(NP4)
!      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+L3*Vin2(NP3)+L4*Vin2(NP4)
!      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+L3*Pin2(NP3)+L4*Pin2(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+(1-L1-L2)*Uin2(NP3)
      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+(1-L1-L2)*Vin2(NP3)
      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+(1-L1-L2)*Pin2(NP3)
!
!-----------------------------------------------------------------------
!              
!            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,' UAIR=',TRA00,
!     &                                              ' VAIR=',TRA02
!            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,' UAIR=',TRA00,
!     &                                               ' VAIR=',TRA02   
!
!-----------------------------------------------------------------------
!           
      ENDIF
!     NEXT TIMESTEPS (LT>0) 
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!       INTERPOLATE THE WIND IN TIME
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+   
!
!     Nth RECORD OF NLATS*NLONS VALUES
!     NEXT RECORD UNTIL "NRECS"
!
600      CONTINUE
!
!-----------------------------------------------------------------------
!
         IF((SHIFTTIME+AT).GT.AT1.AND.(SHIFTTIME+AT).GT.AT2
     &      .AND.(SHIFTTIME+AT).LT.((TIME_in(NRECS))*60)) THEN

            AT1=AT2
            TRA00=TRA03
            TRA02=TRA04
            PATMOS1=PATMOS2
            rec = rec + 1
            start(3) = rec
!
!-----------------------------------------------------------------------
!
!          READING THE REMAINING SETS OF DATA FROM NetCDF FILE
!
!-----------------------------------------------------------------------
!
            status = NF90_GET_VAR(ncid,u_varid,WIND_U_in,start,count)
            status = NF90_GET_VAR(ncid,v_varid,WIND_V_in,start,count)
            status = NF90_GET_VAR(ncid,p_varid,PRES_in,start,count)
!            status = NF90_GET_VAR(ncid,abs_varid,WIND_ABS_in,start,count)

            AT2 = (TIME_in(rec))*60
!
!-----------------------------------------------------------------------
!     LOADING & CHECKING THE VARIABLES       
!-----------------------------------------------------------------------
! 
            Uin2=(/WIND_U_in(:,:)/)
            Vin2=(/WIND_V_in(:,:)/)
            Pin2=(/PRES_in(:,:)/)

            WRITE(LU,*) 'AT',(SHIFTTIME+AT),AT1,AT2
            WRITE(LU,*) 'TIME_in(rec)',TIME_in(rec)
            WRITE(LU,*) REC_NAME, rec
            WRITE(LU,*) WIND_U_NAME, Uin2(1:3)
            WRITE(LU,*) WIND_V_NAME, Vin2(1:3)
            WRITE(LU,*) PRES_NAME, Pin2(1:3)
            WRITE(LU,*) REC_NAME, TIME_in(1:3)

              IF(MAXVAL(Uin2)>200.D0) GO TO 300
              IF(MAXVAL(Vin2)>200.D0) GO TO 300
              IF(MINVAL(Uin2)<-200.D0) GO TO 300
              IF(MINVAL(Vin2)<-200.D0) GO TO 300
              IF(MAXVAL(Pin2)>110000.D0) GO TO 300
              IF(MAXVAL(Pin2)<90000.D0) GO TO 300
!
!-----------------------------------------------------------------------
!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!              CALL FASP(X,Y,TRA03,NPOIN,Xin,Yin,Uin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,TRA04,NPOIN,Xin,Yin,Vin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS2,NPOIN,Xin,Yin,Pin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+L3*Uin2(NP3)+L4*Uin2(NP4)
!      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+L3*Vin2(NP3)+L4*Vin2(NP4)
!      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+L3*Pin2(NP3)+L4*Pin2(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+(1-L1-L2)*Uin2(NP3)
      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+(1-L1-L2)*Vin2(NP3)
      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+(1-L1-L2)*Pin2(NP3)
!
!-----------------------------------------------------------------------
!
            GO TO 600
         ELSE
            IF(AT2-AT1.GE.0.0166) COEF=((SHIFTTIME+AT)-AT1)/(AT2-AT1)
            IF(AT2-AT1.LT.0.0166) COEF=0.D0
         ENDIF

!         WINDX=TRA00+COEF*(TRA03-TRA00)
!         WINDY=TRA02+COEF*(TRA04-TRA02)
!         PATMOS=PATMOS1+COEF*(PATMOS2-PATMOS1)

        IF ((SHIFTTIME+AT).LT.((TIME_in(1))*60)) COEF=0.D0
        IF ((SHIFTTIME+AT).GT.((TIME_in(NRECS))*60)) COEF=1.D0

        TRA01=TRA00+COEF*(TRA03-TRA00)
        CALL OV('X=Y     ',WINDX,TRA01,Z,UAIR,NPOIN)
        TRA01=TRA02+COEF*(TRA04-TRA02)
        CALL OV('X=Y     ',WINDY,TRA01,Z,VAIR,NPOIN) 
        TRA01=PATMOS1+COEF*(PATMOS2-PATMOS1)
        CALL OV('X=Y     ',PATMOS,TRA01,Z,VAIR,NPOIN)   

!          IF(LISTIN) THEN
!            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,
!     &                                              ' UAIR=',WINDX,
!     &                                              ' VAIR=',WINDY,
!     &                                              ' ATM=',PATMOS
!            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,
!     &                                              ' UAIR=',WINDX,
!     &                                              ' VAIR=',WINDY,
!     &                                              ' ATM=',PATMOS 
!          ENDIF 

!
!-----------------------------------------------------------------------
!
!     CLOSE THE NetCDF FILE
!      status = NF90_CLOSE(ncid) 
!
!-----------------------------------------------------------------------
!
!      WRITE(LU,*) REC_NAME, AT2
!      WRITE(LU,*) WIND_U_NAME, Uin2
!      WRITE(LU,*) WIND_V_NAME, Vin2
!      WRITE(LU,*) PRES_NAME, Pin2
!      WRITE(LU,*) REC_NAME, TIME_in

!      IF(LT.EQ.6) STOP
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!      CHOOSE INPUT FILE (NetCDF)
!      OPTION 2

       ELSE

!      CHOOSE INPUT FILE (ASCII)
!      OPTION 1
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     BEWARE, ONLY ONE CALCULATION AT FIRST TIMESTEP IS NEEDED
!
      IF(LT.EQ.0) THEN
!
!        UL=T2D_FILES(T2DFO1)%LU
        UL=FILES(FO1)%LU
!
!-----------------------------------------------------------------------
!
!       ATMOSPHERIC PRESSURE
!
        IF(ATMOS) THEN
          P0 = 100000.D0
          CALL OV( 'X=C     ' , PATMOS , Y , Z , P0 , NPOIN )
        ENDIF
!
!-----------------------------------------------------------------------
!
!       WIND : IN THIS CASE THE WIND IS CONSTANT,
!              VALUE GIVEN IN STEERING FILE.
!
!       MAY REQUIRE A ROTATION,
!       DEPENDING ON THE SYSTEM IN WHICH THE WIND VELOCITY WAS SUPPLIED
!
        IF(VENT) THEN
          CALL OV( 'X=C     ' , WINDX , Y , Z , FUAIR , NPOIN )
          CALL OV( 'X=C     ' , WINDY , Y , Z , FVAIR , NPOIN )
        ENDIF
!
        IF(MY_OPTION.EQ.2) THEN

!         READING THE INITIAL INFORMATION        

!         JUMPING ONE LINE OF COMMENTS
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200) TIMEin,XCELLin,YCELLin
          
          NPin=XCELLin*YCELLin

          ALLOCATE (Pin1(NPin),Pin2(NPin))
          ALLOCATE (Uin1(NPin),Uin2(NPin))
          ALLOCATE (Vin1(NPin),Vin2(NPin))
          ALLOCATE (LONGin(NPin),LATin(NPin))

!         READING THE COORDINATES

!         JUMPING ONE LINE OF COMMENTS
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200) LONGin    

!         JUMPING ONE LINE OF COMMENTS
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200) LATin 
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!         min and max of MERCATOR coordinates
!         φmax = ±85.05113° 
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
        IF(MINVAL(LATin).LT.-85.05113.OR.MAXVAL(LATin).GT.85.05113) THEN
!          IF(LNG.EQ.1) WRITE(LU,*) 'POINT OUTSIDE OF THE DOMAIN'
!          IF(LNG.EQ.1) WRITE(LU,*) 'SEGMENT INTERSECT DOMAIN BOUNDARY'
!          IF(LNG.EQ.2) WRITE(LU,*) 'POINT OUTSIDE OF THE DOMAIN'
!          IF(LNG.EQ.2) WRITE(LU,*) 'SEGMENT INTERSECT DOMAIN BOUNDARY'  
          IF(LNG.EQ.1) WRITE(LU,*) 'LATITUDE OUT OF RANGE'
          IF(LNG.EQ.1) WRITE(LU,*) 'φmax = ±85.05113°'
          IF(LNG.EQ.2) WRITE(LU,*) 'LATITUDE OUT OF RANGE'
          IF(LNG.EQ.2) WRITE(LU,*) 'φmax = ±85.05113°'
        END IF           
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!         function [X, Y] = lonlat2mercator_telemac(lon, lat)
!         % author : Olivier Gourgue
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
           ALLOCATE (lon0(NPin),lat0(NPin),Xin(NPin),Yin(NPin)) 

!         % lon and lat are in degrees (lonlat degrees WGS 84 - EPSG code 4326)
    
!         % "mercator_telemac" is the coordinate system required by Telemac if we
!         % want to take into account the sphericity of the domain
!         % It is not an official coordinate system (no EPSG code) but the
!         % transformation formulae from lonlat degrees WGS 84 are provided
    
!         % earth radius
           Radius = 6371000.D0
!         % number π = 3.1415926535897932384626433832795 
           pi = 4.*atan(1.)
!         % arbitrary reference meridian (0° is Greenwich meridian)
           lon0 = 0.D0*pi/180.D0
           LONGin=LONGin*pi/180.D0
!         % arbitrary reference parallel
!         % '°' *pi/180 = 'Radians'
           lat0 = 50.D0*pi/360.D0+pi/4.D0
           LATin=LATin*pi/360.D0+pi/4.D0

           Xin=Radius*(LONGin-lon0)
           Yin=Radius*(log(tan(LATin))-log(tan(lat0)))
!        
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!         JUMPING FOUR LINES OF COMMENTS
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200)
          READ(UL,*,ERR=100,END=200)

!         READING THE FIRST TWO SETS OF DATA FROM FILE
!
!-----------------------------------------------------------------------
!
          READ(UL,*,ERR=100,END=200) AT1,Uin1,Vin1,Pin1
!          READ(UL,*,ERR=100,END=200) AT1,TRA00,TRA02,PATMOS
              IF(MAXVAL(Uin1)>200.D0) GO TO 500
              IF(MAXVAL(Vin1)>200.D0) GO TO 500
              IF(MINVAL(Uin1)<-200.D0) GO TO 500
              IF(MINVAL(Vin1)<-200.D0) GO TO 500           
              IF(MAXVAL(Pin1)>110000.D0) GO TO 500              
              IF(MAXVAL(Pin1)<90000.D0) GO TO 500

!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'

!              CALL FASP(X,Y,TRA00,NPOIN,Xin,Yin,Uin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)                
!              CALL FASP(X,Y,TRA02,NPOIN,Xin,Yin,Vin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS1,NPOIN,Xin,Yin,Pin1,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA00=L1*Uin1(NP1)+L2*Uin1(NP2)+L3*Uin1(NP3)+L4*Uin1(NP4)
!      TRA02=L1*Vin1(NP1)+L2*Vin1(NP2)+L3*Vin1(NP3)+L4*Vin1(NP4)
!      PATMOS1=L1*Pin1(NP1)+L2*Pin1(NP2)+L3*Pin1(NP3)+L4*Pin1(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA00=L1*Uin1(NP1)+L2*Uin1(NP2)+(1-L1-L2)*Uin1(NP3)
      TRA02=L1*Vin1(NP1)+L2*Vin1(NP2)+(1-L1-L2)*Vin1(NP3)
      PATMOS1=L1*Pin1(NP1)+L2*Pin1(NP2)+(1-L1-L2)*Pin1(NP3)
!
!-----------------------------------------------------------------------
! 
          READ(UL,*,ERR=100,END=200) AT2,Uin2,Vin2,Pin2
!          READ(UL,*,ERR=100,END=200) AT2,TRA03,TRA04,PATMOS2
              IF(MAXVAL(Uin2)>200.D0) GO TO 300
              IF(MAXVAL(Vin2)>200.D0) GO TO 300
              IF(MINVAL(Uin2)<-200.D0) GO TO 300
              IF(MINVAL(Vin2)<-200.D0) GO TO 300
              IF(MAXVAL(Pin2)>110000.D0) GO TO 300
              IF(MAXVAL(Pin2)<90000.D0) GO TO 300

!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'

!              CALL FASP(X,Y,TRA03,NPOIN,Xin,Yin,Uin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,TRA04,NPOIN,Xin,Yin,Vin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS2,NPOIN,Xin,Yin,Pin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+L3*Uin2(NP3)+L4*Uin2(NP4)
!      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+L3*Vin2(NP3)+L4*Vin2(NP4)
!      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+L3*Pin2(NP3)+L4*Pin2(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+(1-L1-L2)*Uin2(NP3)
      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+(1-L1-L2)*Vin2(NP3)
      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+(1-L1-L2)*Pin2(NP3)
!
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!              
            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,' UAIR=',TRA00,
     &                                              ' VAIR=',TRA02
            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,' UAIR=',TRA00,
     &                                               ' VAIR=',TRA02   
!
!-----------------------------------------------------------------------
!       
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(MY_OPTION.EQ.2.AND.VENT) THEN

!-----------------------------------------------------------------------
!       'INTERPOLATE THE WIND IN TIME'
!-----------------------------------------------------------------------

400      CONTINUE
         IF((SHIFTTIME+AT).GT.AT1.AND.(SHIFTTIME+AT).GT.AT2) THEN
            AT1=AT2
            TRA00=TRA03
            TRA02=TRA04
            PATMOS1=PATMOS2
!
!          READING THE REMAINING SETS OF DATA FROM FILE
!
!-----------------------------------------------------------------------
!
            READ(UL,*,ERR=100,END=200) AT2,Uin2,Vin2,Pin2
!            READ(UL,*,ERR=100,END=200) AT2,TRA03,TRA04,PATMOS2
                IF(MAXVAL(Uin2)>200.D0) GO TO 300
                IF(MAXVAL(Vin2)>200.D0) GO TO 300
                IF(MINVAL(Uin2)<-200.D0) GO TO 300
                IF(MINVAL(Vin2)<-200.D0) GO TO 300
                IF(MAXVAL(Pin2)>110000.D0) GO TO 300
                IF(MAXVAL(Pin2)<90000.D0) GO TO 300

!              'CALLING SUBROUTINE FASP'
!              'INTERPOLATE THE WINDS IN SPACE FROM FILE TO MESH'

!              CALL FASP(X,Y,TRA03,NPOIN,Xin,Yin,Uin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,TRA04,NPOIN,Xin,Yin,Vin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!              CALL FASP(X,Y,PATMOS2,NPOIN,Xin,Yin,Pin2,NPin,NBOR,KP1BOR,
!     &                    NPTFR,0.D0)
!
!-----------------------------------------------------------------------
!              'USING WEIGHTING FACTORS (L1,L2,L3,L4)'
!              'INTERPOLATE THE WIND & PRESSURE IN SPACE FROM FILE TO MESH'
!-----------------------------------------------------------------------
! 
!         BILINEAR INTERPOLATION 
!
!      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+L3*Uin2(NP3)+L4*Uin2(NP4)
!      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+L3*Vin2(NP3)+L4*Vin2(NP4)
!      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+L3*Pin2(NP3)+L4*Pin2(NP4)
!
!         BARYCENTRIC INTERPOLATION    
! 
      TRA03=L1*Uin2(NP1)+L2*Uin2(NP2)+(1-L1-L2)*Uin2(NP3)
      TRA04=L1*Vin2(NP1)+L2*Vin2(NP2)+(1-L1-L2)*Vin2(NP3)
      PATMOS2=L1*Pin2(NP1)+L2*Pin2(NP2)+(1-L1-L2)*Pin2(NP3)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                
            GO TO 400
         ELSE
            IF(AT2-AT1.GE.0.0166) COEF=((SHIFTTIME+AT)-AT1)/(AT2-AT1)
            IF(AT2-AT1.LT.0.0166) COEF=0.D0
         ENDIF

!         WINDX=TRA00+COEF*(TRA03-TRA00)
!         WINDY=TRA02+COEF*(TRA04-TRA02)
!         PATMOS=PATMOS1+COEF*(PATMOS2-PATMOS1)

        TRA01=TRA00+COEF*(TRA03-TRA00)
        CALL OV('X=Y     ',WINDX,TRA01,Z,UAIR,NPOIN)
        TRA01=TRA02+COEF*(TRA04-TRA02)
        CALL OV('X=Y     ',WINDY,TRA01,Z,VAIR,NPOIN) 
        TRA01=PATMOS1+COEF*(PATMOS2-PATMOS1)
        CALL OV('X=Y     ',PATMOS,TRA01,Z,VAIR,NPOIN)   

!          IF(LISTIN) THEN
!            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,
!     &                                              ' UAIR=',WINDX,
!     &                                              ' VAIR=',WINDY,
!     &                                              ' ATM=',PATMOS
!            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,
!     &                                              ' UAIR=',WINDX,
!     &                                              ' VAIR=',WINDY,
!     &                                              ' ATM=',PATMOS 
!          ENDIF 

!-----------------------------------------------------------------------      
!
!       JUMPING TWO LINES OF COMMENTS
!
! 10      CONTINUE
!        IF(AT.GE.AT1.AND.AT.LT.AT2) THEN
!          IF(AT2-AT1.GT.1.D-6) THEN
!            COEF=(AT-AT1)/(AT2-AT1)
!          ELSE
!            COEF=0.D0
!          ENDIF
!          TRA05=TRA00+COEF*(TRA03-TRA00)
!          TRA06=TRA02+COEF*(TRA04-TRA02)
!          IF(LISTIN) THEN
!            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,' UAIR=',TRA00,
!     &                                              ' VAIR=',TRA02
!            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,' UAIR=',TRA00,
!     &                                               ' VAIR=',TRA02 
!          ENDIF
!        ELSE
!          AT1=AT2
!          TRA00=TRA03
!          TRA02=TRA04
!          READ(UL,*,ERR=100,END=200) AT2,TRA03,TRA04
!          GO TO 10
!        ENDIF
!        TRA01=TRA05
!        CALL OV('X=Y     ',WINDX,TRA01,Z,UAIR,NPOIN)
!        TRA01=TRA06
!        CALL OV('X=Y     ',WINDY,TRA01,Z,VAIR,NPOIN)    
!            IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,' UAIR=',WINDX,
!     &                                              ' VAIR=',WINDY
!            IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,' UAIR=',WINDX,
!     &                                               ' VAIR=',WINDY         
!
!-----------------------------------------------------------------------

      ENDIF
!
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
!       
!      CHOOSE INPUT FILE (ASCII)
      END IF
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! 
      RETURN
!
!-----------------------------------------------------------------------
! 
100   CONTINUE
      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'ERREUR DANS LE FICHIER DE VENT'
      IF(LNG.EQ.2) WRITE(LU,*) 'ERROR IN THE WIND & PRESSURE FILE'
      CALL PLANTE(1)
      STOP  
200   CONTINUE
      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'FIN PREMATUREE DU FICHIER DE VENT'
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND & PRESSURE FILE TOO SHORT'
      CALL PLANTE(1)
      STOP           
300   CONTINUE
      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT
      IF(LNG.EQ.1) WRITE(LU,*)  'U', MAXVAL(TRA03),MINVAL(TRA03)
      IF(LNG.EQ.1) WRITE(LU,*)  'V', MAXVAL(TRA04),MAXVAL(TRA04)
      IF(LNG.EQ.1) WRITE(LU,*)  'P', MAXVAL(PATMOS2),MAXVAL(PATMOS2)
      IF(LNG.EQ.1) WRITE(LU,*) 'VENT DU FICHIER HORS LIMITE -200< <200'
      IF(LNG.EQ.1) WRITE(LU,*) 'PRESSION HORS LIMITE 900< <110000'
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT
      IF(LNG.EQ.2) WRITE(LU,*)  'U', MAXVAL(TRA03),MINVAL(TRA03)
      IF(LNG.EQ.2) WRITE(LU,*)  'V', MAXVAL(TRA04),MAXVAL(TRA04)
      IF(LNG.EQ.2) WRITE(LU,*)  'P', MAXVAL(PATMOS2),MAXVAL(PATMOS2)
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND OUT OF RANGE -200< <200'
      IF(LNG.EQ.2) WRITE(LU,*) 'PRESSURE OUT OF RANGE 900< <110000'
      CALL PLANTE(1)
      STOP
500   CONTINUE
      WRITE(LU,*) ' '
      WRITE(LU,*) 'METEO'
      IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT
      IF(LNG.EQ.1) WRITE(LU,*)  'U', MAXVAL(TRA00),MINVAL(TRA00)
      IF(LNG.EQ.1) WRITE(LU,*)  'V', MAXVAL(TRA02),MAXVAL(TRA02)
      IF(LNG.EQ.1) WRITE(LU,*)  'P', MAXVAL(PATMOS1),MAXVAL(PATMOS1)
      IF(LNG.EQ.1) WRITE(LU,*) 'VENT DU FICHIER HORS LIMITE -200< <200'
      IF(LNG.EQ.1) WRITE(LU,*) 'PRESSION HORS LIMITE 900< <110000'
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT
      IF(LNG.EQ.2) WRITE(LU,*)  'U', MAXVAL(TRA00),MINVAL(TRA00)
      IF(LNG.EQ.2) WRITE(LU,*)  'V', MAXVAL(TRA02),MAXVAL(TRA02)
      IF(LNG.EQ.2) WRITE(LU,*)  'P', MAXVAL(PATMOS1),MAXVAL(PATMOS1)
      IF(LNG.EQ.2) WRITE(LU,*) 'WIND OUT OF RANGE -200< <200'
      IF(LNG.EQ.2) WRITE(LU,*) 'PRESSURE OUT OF RANGE 900< <110000'
      CALL PLANTE(1)
      STOP  
!
!-----------------------------------------------------------------------
!
      RETURN
      END
!                    ***************
                     SUBROUTINE FASP
!                    ***************
!
     &(X,Y,ZF,NPOIN,XRELV,YRELV,ZRELV,NP,NBOR,KP1BOR,NPTFR,DM)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    INTERPOLATES THE BOTTOM ELEVATIONS FROM A SET OF
!+                POINTS ON THE MESH NODES.
!
!history  J-M HERVOUET (LNHE)
!+        20/03/08
!+        V5P9
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  EHSAN SARHADI ZADEH (AnteaGroup, BELGIUM)
!+        23/04/2014
!+        V6P3
!+   Modification for netCDF wind interpolating 
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| DM             |-->| MINIMUM DISTANCE TO BOUNDARY TO ACCEPT A POINT
!| KP1BOR         |-->| GIVES THE NEXT BOUNDARY POINT IN A CONTOUR
!| NBOR           |-->| GLOBAL NUMBER OF BOUNDARY POINTS
!| NP             |-->| NUMBER OF BATHYMETRY POINTS
!| NPOIN          |-->| NUMBER OF POINTS IN THE MESH
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| X,Y            |-->| MESH COORDINATES
!| XRELV          |-->| ABCISSAE OF BATHYMETRY POINTS
!| YRELV          |-->| ORDINATES OF BATHYMETRY POINTS
!| ZF             |<--| INTERPOLATED BATHYMETRY 
!| ZRELV          |-->| ELEVATIONS OF BATHYMETRY POINTS
!| OK             |<--| IF YES, NO CROSSING OF BOUNDARY
!|                |   | IF NO, POINT OUTSIDE THE DOMAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF, EX_FASP => FASP
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPOIN,NP,NPTFR
      INTEGER, INTENT(IN) :: NBOR(NPTFR),KP1BOR(NPTFR)
      DOUBLE PRECISION, INTENT(IN)  :: X(NPOIN),Y(NPOIN),DM
      DOUBLE PRECISION, INTENT(IN)  :: XRELV(NP),YRELV(NP),ZRELV(NP)
      DOUBLE PRECISION, INTENT(OUT) :: ZF(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N,INUM,I
!
      DOUBLE PRECISION DIST1,DIST2,DIST3,DIST4
      DOUBLE PRECISION ZCADR1,ZCADR2,ZCADR3,ZCADR4
      DOUBLE PRECISION DIFX,DIFY,DIST,X1,Y1,X2,Y2,X3,Y3,X4,Y4
      DOUBLE PRECISION ZNUM,ZDEN
!
      LOGICAL OK1,OK2,OK3,OK4
!
!-----------------------------------------------------------------------
!
!  LOOP ON THE MESH NODES:
!
      DO 100 I = 1 , NPOIN
!
!     INTERPOLATES THE BOTTOM FROM 4 QUADRANTS
!
! ---->  INITIALISES:
!
!     min and max of MERCATOR coordinates
!     φmax = ±85.05113° 
!
      DIST1=625.D12
      DIST2=625.D12
      DIST3=625.D12
      DIST4=625.D12
!
      OK1 = .FALSE.
      OK2 = .FALSE.
      OK3 = .FALSE.
      OK4 = .FALSE.
!
      ZCADR1=0.D0
      ZCADR2=0.D0
      ZCADR3=0.D0
      ZCADR4=0.D0
!
! --------->  LOOP ON THE SET OF POINTS (THERE ARE NP):
      DO 30 N=1,NP
           DIFX = XRELV(N)-X(I)
           DIFY = YRELV(N)-Y(I)
           DIST = DIFX*DIFX + DIFY*DIFY
!
             IF ( DIST.LT.1.D-6 ) DIST=1.D-6
! ->QUADRANT 1 :
               IF( DIFX.LE.0.D0.AND.DIFY.LE.0.D0) THEN
                 IF(DIST.LE.DIST1)THEN
                      X1=XRELV(N)
                      Y1=YRELV(N)
                      DIST1=DIST
                      ZCADR1=ZRELV(N)
                      OK1 = .TRUE.
                 ENDIF
! ->QUADRANT 2 :
              ELSE IF( DIFX.GE.0.D0.AND.DIFY.LE.0.D0) THEN
                 IF(DIST.LE.DIST2)THEN
                      X2=XRELV(N)
                      Y2=YRELV(N)
                      DIST2=DIST
                      ZCADR2=ZRELV(N)
                      OK2 = .TRUE.
                 ENDIF
! ->QUADRANT 3 :
              ELSE IF( DIFX.GE.0.D0.AND.DIFY.GE.0.D0) THEN
                 IF(DIST.LE.DIST3)THEN
                      X3=XRELV(N)
                      Y3=YRELV(N)
                      DIST3=DIST
                      ZCADR3=ZRELV(N)
                      OK3 = .TRUE.
                 ENDIF
! ->QUADRANT 4 :
              ELSE IF( DIFX.LE.0.D0.AND.DIFY.GE.0.D0) THEN
                 IF(DIST.LE.DIST4)THEN
                      X4=XRELV(N)
                      Y4=YRELV(N)
                      DIST4=DIST
                      ZCADR4=ZRELV(N)
                      OK4 = .TRUE.
                 ENDIF
              ENDIF
 30        CONTINUE
!
! --------->  END OF LOOP ON THE SET OF POINTS
!
!-----------------------------------------------------------------------
!
!+            CHECKS HERE THAT THIS POINT IS NOT OUTSIDE OF THE
!+                DOMAIN, I.E. CHECKS THAT THE SEGMENT LINKING (X,Y)
!+                AND (XR,YR) DOES NOT INTERSECT WITH THE DOMAIN
!+                BOUNDARY.
!
!-----------------------------------------------------------------------
!
!
!      IF(OK1) CALL CROSFR(X(I),Y(I),X1,Y1,X,Y,NPOIN,NBOR,KP1BOR,
!     &                    NPTFR,DM,OK1)
!      IF(OK2) CALL CROSFR(X(I),Y(I),X2,Y2,X,Y,NPOIN,NBOR,KP1BOR,
!     &                    NPTFR,DM,OK2)
!      IF(OK3) CALL CROSFR(X(I),Y(I),X3,Y3,X,Y,NPOIN,NBOR,KP1BOR,
!     &                    NPTFR,DM,OK3)
!      IF(OK4) CALL CROSFR(X(I),Y(I),X4,Y4,X,Y,NPOIN,NBOR,KP1BOR,
!     &                    NPTFR,DM,OK4)
!
!
!-----------------------------------------------------------------------
!
         ZNUM = 0.D0
         ZDEN = 0.D0
         INUM = 0
         IF(OK1) THEN
          ZNUM = ZNUM + ZCADR1/DIST1
          ZDEN = ZDEN + 1.D0/DIST1
          INUM = INUM + 1
         ENDIF
         IF(OK2) THEN
          ZNUM = ZNUM + ZCADR2/DIST2
          ZDEN = ZDEN + 1.D0/DIST2
          INUM = INUM + 1
         ENDIF
         IF(OK3) THEN
          ZNUM = ZNUM + ZCADR3/DIST3
          ZDEN = ZDEN + 1.D0/DIST3
          INUM = INUM + 1
         ENDIF
         IF(OK4) THEN
          ZNUM = ZNUM + ZCADR4/DIST4
          ZDEN = ZDEN + 1.D0/DIST4
          INUM = INUM + 1
         ENDIF
!
         IF(INUM.NE.0) THEN
!         ZF : WIND AT THE POINT
          ZF(I)=ZNUM/ZDEN
         ELSE        
!         ZF : NO WIND AT THE POINT
          ZF(I) = -1.D6
         ENDIF
!
100   CONTINUE
!
!-----------------------------------------------------------------------
!
      RETURN
      END

!                    *****************
                     SUBROUTINE PROSOU
!                    *****************
!
     &(FU,FV,SMH,    UN,VN,HN,GRAV,NORD,
     & FAIR,WINDX,WINDY,VENT,HWIND,CORIOL,FCOR,
     & SPHERI,YASMH,COSLAT,SINLAT,AT,LT,DT,
     & NREJET,NREJEU,DSCE,ISCE,T1,MESH,MSK,MASKEL,
     & MAREE,MARDAT,MARTIM,PHI0,OPTSOU,
     & COUROU,NPTH,VARCL,NVARCL,VARCLA,UNSV2D,
     & FXWAVE,FYWAVE,RAIN,RAIN_MMPD,PLUIE,
     & T2D_FILES,T2DBI1,BANDEC,OPTBAN,
     & NSIPH,ENTSIP,SORSIP,DSIP,USIP,VSIP,
     & NBUSE,ENTBUS,SORBUS,DBUS,UBUS,VBUS,
     & TYPSEUIL,NWEIRS,
     & NPSING,NDGA1,NDGB1,NBOR)
!
!***********************************************************************
! TELEMAC2D   V6P3                                   21/08/2010
!***********************************************************************
!
!brief    PREPARES THE SOURCE TERMS IN THE CONTINUITY EQUATION
!+                AND IN THE DYNAMIC EQUATIONS. ARE TAKEN INTO ACCOUNT :
!+
!+              - WIND
!+
!+              - CORIOLIS FORCE
!+
!+              - TIDAL FORCE
!+
!+              - SOURCES AND SINKS
!+
!+              - WEIRS (IF TYPSEUIL=2)
!code
!+    RESPECTIVE TERMS ARE:
!+    ==========================
!+
!+     * WIND
!+       ---------
!+                                 1                         2      2
!+                FU           =  --- * F    * U    * SQRT( U    + V   )
!+                  VENT           H     AIR    AIR          AIR    AIR
!+
!+                                 1                         2      2
!+                FV           =  --- * F    * V    * SQRT( U    + V   )
!+                  VENT           H     AIR    AIR          AIR    AIR
!+
!+           WHERE :
!+                  UAIR   :  WIND VELOCITY ALONG X
!+                  VAIR   :  WIND VELOCITY ALONG Y
!+                  FAIR   :  AIR FRICTION COEFFICIENT
!+
!+     * CORIOLIS FORCE
!+       ---------------------
!+
!+                FU           =  + FCOR * V
!+                  CORIOLIS
!+
!+                FV           =  - FCOR * U
!+                  CORIOLIS
!+
!+           WHERE :
!+                  U       :  FLOW VELOCITY ALONG X
!+                  V       :  FLOW VELOCITY ALONG Y
!+                  FCOR    :  CORIOLIS PARAMETER
!
!note     BOTTOM FRICTION IS TAKEN INTO ACCOUNT IN THE PROPAGATION
!+         THROUGH CALL TO FROTXY, IT IS SEMI-IMPLICIT.
!note  IF SOURCES OR SINKS TERMS ARE ADDED TO THE CONTINUITY EQUATION,
!+         IT IS IDENTIFIED WITH VARIABLE YASMH (SET TO TRUE).
!note  SOURCE TERMS FU AND FV ARE FIRST COMPUTED IN P1.
!+         THEY ARE THEN EXTENDED TO QUASI-BUBBLE IF REQUIRED.
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  J-M HERVOUET (LNHE)
!+        20/02/2012
!+        V6P2
!+   Rain-evaporation added (after initial code provided by O. Boutron, 
!+   Tour du Valat and O. Bertrand, Artelia-group).
!
!history  C.COULET (ARTELIA)
!+        23/05/2012
!+        V6P2
!+   Modification for culvert management
!+   Addition of Tubes management
!
!history  R. KOPMANN (EDF R&D, LNHE)
!+        16/04/2013
!+        V6P3
!+   Adding the file format in calls to FIND_IN_SEL.
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        21/05/2013
!+        V6P3
!+   Possibility of negative depths taken into account when limiting
!+   evaporation on dry zones.
!
!history  C.COULET / A.REBAI / E.DAVID (ARTELIA)
!+        07/06/2013
!+        V6P3
!+   Modification for new treatment of weirs
!
!history  EHSAN SARHADI ZADEH (AnteaGroup, BELGIUM)
!+        23/04/2014
!+        V6P3
!+   Modification for reduced drag coefficient for high wind speeds
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| α (ALFA)       |-->| CONSTANT CHARNOCK PARAMETER
!| β (BETA)       |-->| VARIABLE CHARNOCK PARAMETER
!| κ (KAPPA)      |-->| von Karman CONSTANT
!| ω (OMEGA)      |-->| SPRAY IMPACT
!| ρair           |-->| DEFAULT AIR DENSITY ~ 1.29 kg/m3 
!| ρwater         |-->| DEFAULT WATER DENSITY ~ 1030.0 kg/m3
!| Acr            |-->| CRITICAL TERMINAL FALL VELOCITY OF DROPLET WITH RADIUS OF ~ 80 μm
!| AT             |-->| TIME
!| BANDEC         |-->| IF YES, TIDAL FLATS OR DRY ZONES
!| Cdrag          |-->| VARIABLE SURFACE DRAG COEFFICIENT FOR WIND
!| CORIOL         |-->| IF YES, CORIOLIS FORCE
!| COSLAT         |-->| COSINUS OF LATITUDE (SPHERICAL COORDINATES)
!| COUROU         |-->| IF YES, WAVE DRIVEN CURRENTS TAKEN INTO ACCOUNT
!| DENSCOEF       |-->| ρair / ρwater
!| DSIP           |-->| DISCHARGE OF TUBES
!| DSCE           |-->| DISCHARGE OF POINT SOURCES
!| DSIP           |-->| DISCHARGE OF CULVERTS
!| DT             |-->| TIME STEP IN SECONDS
!| ENTBUS         |-->| INDICES OF ENTRY OF TUBES IN GLOBAL NUMBERING
!| ENTSIP         |-->| INDICES OF ENTRY OF CULVERTS IN GLOBAL NUMBERING
!| FAIR           |-->| FRICTION COEFFICIENT FOR WIND
!| FCOR           |-->| CORIOLIS PARAMETER
!| FU             |<->| SOURCE TERMS ON VELOCITY U
!| FV             |<->| SOURCE TERMS ON VELOCITY V
!| FXWAVE         |<->| FORCING OF WAVES ALONG X
!| FYWAVE         |<->| FORCING OF WAVES ALONG Y
!| GRAV (g)       |-->| GRAVITY
!| HN             |-->| DEPTH AT TIME T(N)
!| HWIND          |-->| MINIMUM DEPTH FOR TAKING WIND INTO ACCOUNT
!| ISCE           |-->| NEAREST POINTS TO SOURCES
!| LT             |-->| TIME STEP NUMBER
!| MARDAT         |-->| DATE (YEAR, MONTH,DAY)
!| MAREE          |-->| IF YES, TAKES THE TIDAL FORCE INTO ACCOUNT
!| MARTIM         |-->| TIME (HOUR, MINUTE,SECOND)
!| MASKEL         |-->| MASKING OF ELEMENTS
!|                |   | =1. : NORMAL   =0. : MASKED ELEMENT
!| MESH           |-->| MESH STRUCTURE
!| MSK            |-->| IF YES, THERE IS MASKED ELEMENTS.
!| NBOR           |-->| INDICES OF BOUNDARY POINTS
!| NBUSE          |-->| NUMBER OF TUBES
!| NORD           |-->| DIRECTION OF NORTH WITH RESPECT TO Y AXIS
!|                |   | (TRIGONOMETRIC SENSE) IN DEGREES.
!| NPSING         |-->| NUMBER OF POINTS OF WEIRS
!| NPTH           |-->| RECORD NUMBER IN THE WAVE CURRENTS FILE
!| NREJET         |-->| NUMBER OF POINT SOURCES
!| NREJEU         |-->| NUMBER OF POINT SOURCES WITH GIVEN VELOCITY
!|                |   | IF NREJEU=0 VELOCITY OF SOURCES IS TAKEN EQUAL
!|                |   | TO VELOCITY.
!| NSIPH          |-->| NUMBER OF CULVERTS
!| NDGA1,NDGB1    |-->| INDICES OF POINTS OF WEIRS
!| NVARCL         |-->| NUMBER OF CLANDESTINE VARIABLES
!| NWEIRS         |-->| NUMBER OF WEIRS
!| OPTBAN         |-->| OPTION FOR THE TREATMENT OF TIDAL FLATS
!| OPTSOU         |-->| OPTION FOR THE TREATMENT OF SOURCES
!| PHI0           |-->| LATITUDE OF ORIGIN POINT
!| PLUIE          |-->| BIEF_OBJ STRUCTURE WITH RAIN OR EVAPORATION.
!| RAIN           |-->| IF YES, RAIN OR EVAPORATION TAKEN INTO ACCOUNT
!| RAIN_MMPD      |-->| RAIN OR EVAPORATION IN MM PER DAY
!| SINLAT         |-->| SINUS OF LATITUDE (SPHERICAL COORDINATES)
!| SMH            |-->| SOURCE TERM IN CONTINUITY EQUATION
!| SORBUS         |-->| INDICES OF TUBES EXITS IN GLOBAL NUMBERING
!| SORSIP         |-->| INDICES OF CULVERTS EXISTS IN GLOBAL NUMBERING
!| SPHERI         |-->| IF TRUE : SPHERICAL COORDINATES
!| T1             |<->| WORK BIEF_OBJ STRUCTURE
!| T2D_FILES      |-->| BIEF_FILE STRUCTURE WITH ALL TELEMAC-2D FILES
!| T2D_BI1        |-->| RANK OF BINARY FILE 1
!| TYPSEUIL       |-->| TYPE OS WEIRS (ONLY TYPSEUIL=2 IS MANAGE HERE)
!| UBUS           |-->| VELOCITY U AT TUBE EXTREMITY
!| UNSV2D         |-->| INVERSE OF INTEGRALS OF TEST FUNCTIONS
!| USIP           |-->| VELOCITY U AT CULVERT EXTREMITY
!| VBUS           |-->| VELOCITY V AT TUBE EXTREMITY
!| VARCL          |<->| BLOCK OF CLANDESTINE VARIABLES
!| VARCLA         |-->| NAMES OF CLANDESTINE VARIABLES
!| VENT           |-->| IF YES, WIND IS TAKEN INTO ACCOUNT
!| VSIP           |-->| VELOCITY V AT CULVERT EXTREMITY
!| WINDX          |-->| FIRST COMPONENT OF WIND VELOCITY
!| WINDY          |-->| SECOND COMPONENT OF WIND VELOCITY
!| YASMH          |<->| IF TRUE SMH IS TAKEN INTO ACCOUNT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
!     FOR SEEING OTHER VARIABLES IN DECLARATIONS_TELEMAC2D:
      USE DECLARATIONS_TELEMAC2D, ONLY : QWA,QWB
      USE INTERFACE_TELEMAC2D, EX_PROSOU => PROSOU
      USE M_COUPLING_ESTEL3D
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     WORKING ARRAYS
!
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: T1
!
!-----------------------------------------------------------------------
!
!     VECTORS
!
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: FU,FV,SMH,FXWAVE,FYWAVE
      TYPE(BIEF_OBJ)   , INTENT(IN)    :: MASKEL,UN,VN,HN,UNSV2D
      TYPE(BIEF_OBJ)   , INTENT(IN)    :: WINDX,WINDY,COSLAT,SINLAT
!
!-----------------------------------------------------------------------
!
!     MESH STRUCTURE
!
      TYPE(BIEF_MESH)  , INTENT(INOUT) :: MESH
!
!-----------------------------------------------------------------------
!
      INTEGER          , INTENT(IN)    :: NVARCL,LT,NREJET,NREJEU,OPTSOU
      INTEGER          , INTENT(IN)    :: NPTH,T2DBI1,OPTBAN
      INTEGER          , INTENT(IN)    :: MARDAT(3),MARTIM(3)
      INTEGER          , INTENT(IN)    :: ISCE(NREJET)
      DOUBLE PRECISION , INTENT(IN)    :: HWIND,AT,FAIR,FCOR
      DOUBLE PRECISION , INTENT(IN)    :: DSCE(NREJET)
      DOUBLE PRECISION , INTENT(IN)    :: GRAV,NORD,PHI0,RAIN_MMPD,DT
      CHARACTER(LEN=32), INTENT(IN)    :: VARCLA(NVARCL)
      LOGICAL          , INTENT(IN)    :: VENT,MAREE,CORIOL,SPHERI,MSK
      LOGICAL          , INTENT(IN)    :: COUROU,RAIN,BANDEC
      LOGICAL          , INTENT(INOUT) :: YASMH
      TYPE(BIEF_OBJ)   , INTENT(INOUT) :: VARCL,PLUIE
      TYPE(BIEF_FILE)  , INTENT(IN)    :: T2D_FILES(*)
!
      INTEGER          , INTENT(IN)    :: NSIPH
      INTEGER          , INTENT(IN)    :: ENTSIP(NSIPH),SORSIP(NSIPH)
      DOUBLE PRECISION , INTENT(IN)    :: DSIP(NSIPH)
      DOUBLE PRECISION , INTENT(IN)    :: USIP(NSIPH,2),VSIP(NSIPH,2)
!
      INTEGER          , INTENT(IN)    :: NBUSE
      INTEGER          , INTENT(IN)    :: ENTBUS(NBUSE),SORBUS(NBUSE)
      DOUBLE PRECISION , INTENT(IN)    :: DBUS(NBUSE)
      DOUBLE PRECISION , INTENT(IN)    :: UBUS(NBUSE,2),VBUS(NBUSE,2)
!
      INTEGER          , INTENT(IN)    :: TYPSEUIL,NWEIRS
      TYPE(BIEF_OBJ)   , INTENT(IN)    :: NBOR,NPSING,NDGA1,NDGB1
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N,I,IELMU,IELMH,IELM1,NPOIN,IR,ERR,NP
!
      DOUBLE PRECISION PI,WROT,WD,ATH,RAIN_MPS,SURDT
      DOUBLE PRECISION DENSCOEF,Cdrag,Cdrag0,ALFA,BETA,KAPPA,OMEGA,Acr
!
      CHARACTER*16 NOMX,NOMY
      LOGICAL DEJALU,OKX,OKY,OKC
      DATA DEJALU /.FALSE./
      REAL, ALLOCATABLE :: W(:)
      SAVE W
!
      INTRINSIC SQRT,MAX,ACOS
!
!-----------------------------------------------------------------------
!  EXTRACTS X COORDINATES, NUMBER OF POINTS P1
!                          AND P1 ELEMENT OF THE MESH
!-----------------------------------------------------------------------
!
      IELM1 = MESH%X%ELM
      NPOIN = MESH%NPOIN
!
!-----------------------------------------------------------------------
!  INITIALISES
!-----------------------------------------------------------------------
!
      CALL CPSTVC(UN,FU)
      CALL CPSTVC(VN,FV)
      CALL OS( 'X=0     ' , X=FU )
      CALL OS( 'X=0     ' , X=FV )
!
!-----------------------------------------------------------------------
!
!  COMPUTATION WITH WIND
!  ----------------
!
!                               1                         2     2
!              FU           =  --- * F    * U    * SQRT( U   + V    )
!                VENT           H     AIR    AIR          AIR   AIR
!
!
!                               1                         2     2
!              FV           =  --- * F    * V    * SQRT( U   + V    )
!                VENT           H     AIR    AIR          AIR   AIR
!
!
!  ----------------
!
!  COEFFICIENT OF WIND INFLUENCE
!
!  ----------------
!  ----------------
!
!  REDUCED DRAG COEFFICIENT FOR HIGH WIND SPEEDS 
!
!  ----------------     

      IF(VENT) THEN
!
!  TEMPORARY TREATMENT OF TIDAL FLATS
!  THE WIND EFFECT IS ONLY CONSIDERED IF THE WATER DEPTH IS
!  GREATER THAN 1 M.
!
!  ASSUMES HERE THAT THE WIND IS GIVEN IN P1
!
        DO N=1,NPOIN
          IF (HN%R(N).GT.HWIND) THEN
            WD = SQRT( WINDX%R(N)**2 + WINDY%R(N)**2 )
!
!  ----------------
!
!  A) Flather’s FORMULATION (1976)
!
!  ----------------       

!      IF (WD.LE.5)                 Cdrag=0.565*1.D-3
!      IF (WD.GT.5.AND.WD.LT.19.22) Cdrag=(-0.12+0.137*WD)*1.D-3
!      IF (WD.GE.19.22)             Cdrag=2.513*1.D-3

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'Cdrag FL' , Cdrag

!
!  ----------------
!
!  B) Charnock’s FORMULATION (1955):
!  Cdrag=[κ/ln(g*z/(α*Cdrag*W^2))]^2
!
!  ----------------
!
!  CONSTATNT VALUE FOR α (ALFA)
!  0.011 < α < 0.035 
!  ----------------               
!

!      ALFA=0.011

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'ALFA FL' , ALFA

!
!  ----------------   
!
!  1st GUESS FOR Cdrag IS DEFAULT VALUE

!      Cdrag0=1.D-3
!      KAPPA=0.4

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'DEFAULT Cdrag' , Cdrag0

!  CALCULATING Cdrag BY ITERATION

!      DO I=1,100

!      Cdrag=(KAPPA/(log(GRAV*10/(ALFA*Cdrag0*WD**2))))**2

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'Cdrag FA' ,I, Cdrag

!      IF (ABS(Cdrag-Cdrag0).LE.1.D-10) EXIT 

!      Cdrag0=Cdrag

!      END DO

!
!  ----------------
!
!  C) MODIFIED VERSION OF Charnock’s FORMULATION (1955):
!  Cdrag=[κ/ln(g*z/(α*Cdrag*W^2))]^2
!
!  ----------------
!
!  1. Fairall’s FORMULATION (2003) FOR α (ALFA)
!
!  ----------------               
!

!      IF (WD.LE.10)                ALFA=0.011
!      IF (WD.GT.10.AND.WD.LT.18)   ALFA=0.011+0.000875*(WD-10.)
!      IF (WD.GE.18)                ALFA=0.018

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'ALFA FL' , ALFA

!
!  ----------------   
!
!  1st GUESS FOR Cdrag IS DEFAULT VALUE

!      Cdrag0=1.D-3
!      KAPPA=0.4

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'DEFAULT Cdrag' , Cdrag0

!  CALCULATING Cdrag BY ITERATION

!      DO I=1,100

!      Cdrag=(KAPPA/(log(GRAV*10/(ALFA*Cdrag0*WD**2))))**2

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'Cdrag FA' ,I, Cdrag

!      IF (ABS(Cdrag-Cdrag0).LE.1.D-10) EXIT 

!      Cdrag0=Cdrag

!      END DO

!
!  ----------------    
!
!  2. Makin’s FORMULATION (2005) FOR α (ALFA)
!   
!  ----------------          
!
      ALFA=0.0075+0.02*TANH(0.075*WD-0.75)
      IF (TANH(0.075*WD-0.75).LE.0.) ALFA=0.0075

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'ALFA MA' , ALFA

!
!  ----------------   
!
!  1st GUESS FOR Cdrag IS DEFAULT VALUE

      Cdrag0=1.D-3
      KAPPA=0.4
      Acr=0.64

!      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) 'DEFAULT Cdrag' , Cdrag0

!  CALCULATING Cdrag BY ITERATION

      DO I=1,100

!  ω (OMEGA) THE SPRAY IMPACT

      OMEGA=Acr/(KAPPA*WD*Cdrag0**0.5)

!  WD < 33 m/s => ω > 1 => NO SPRAY IMPACT

      BETA=ALFA

      IF (OMEGA.LT.1.0) THEN

      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) I,'ALFA MA' , ALFA

      BETA=(10.**(1.-1./OMEGA))*(ALFA**(1./OMEGA))

      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) I,'BETA SPRAY ω' , BETA

      END IF 

      Cdrag=(KAPPA/(log(GRAV*10/(BETA*Cdrag0*WD**2))))**2

      IF (ABS(Cdrag-Cdrag0).LE.1.D-10) EXIT 

      IF(LT.EQ.1.AND.N.EQ.NPOIN) WRITE(LU,*) I,'Cdrag Makin', Cdrag

      Cdrag0=Cdrag

      END DO

!
!  ----------------   
!
!  CONSTANT VALUE FOR COEFFICIENT OF WIND INFLUENCE 
!
!  DEFAULT VALUE FOR COEFFICIENT OF WIND INFLUENCE

!            FAIR=Cdrag*DENSCOEF=3.178206E-6
!
!            FU%R(N) = FU%R(N) + FAIR * WINDX%R(N) * WD / HN%R(N)
!            FV%R(N) = FV%R(N) + FAIR * WINDY%R(N) * WD / HN%R(N)
!
!  ----------------   
!
!  VARIABLE VALUE FOR COEFFICIENT OF WIND INFLUENCE 
!            FAIR=Cdrag*DENSCOEF
!
!  DEFAULT AIR DENSITY ρair ~ 1.29 kg/m3 
!  DEFAULT WATER DENSITY ρwater ~ 1030 kg/m3

      DENSCOEF=1.29/1030.
!
!  ----------------   
!
            FU%R(N) = FU%R(N)+Cdrag*DENSCOEF*WINDX%R(N)*WD/HN%R(N)
            FV%R(N) = FV%R(N)+Cdrag*DENSCOEF*WINDY%R(N)*WD/HN%R(N)
          ENDIF
        ENDDO
!
      ENDIF
!
!***********************************************************************
!
!     * WITH CORIOLIS FORCE
!       --------------------------
!
!                FU           =  + FCOR * V
!                  CORIOLIS
!
!                FV           =  - FCOR * U
!                  CORIOLIS
!
      IF(CORIOL) THEN
!
      PI = ACOS(-1.D0)
!
        IF(SPHERI) THEN
!
          WROT = 2 * PI / 86164.D0
          DO I=1,NPOIN
!           FORMULATION INDEPENDENT OF THE DIRECTION OF NORTH
            FU%R(I) = FU%R(I) + VN%R(I) * 2 * WROT * SINLAT%R(I)
            FV%R(I) = FV%R(I) - UN%R(I) * 2 * WROT * SINLAT%R(I)
          ENDDO
!
!         TAKES THE TIDAL FORCE INTO ACCOUNT
!
          IF(MAREE) THEN
            CALL MARAST(MARDAT,MARTIM,PHI0,NPOIN,AT,
     &                  FU%R,FV%R,MESH%X%R,SINLAT%R,COSLAT%R,GRAV)
          ENDIF
!
          IF(LT.EQ.1) THEN
            IF(LNG.EQ.1) WRITE(LU,11)
            IF(LNG.EQ.2) WRITE(LU,12)
          ENDIF
11        FORMAT(1X,'PROSOU : EN COORDONNEES SHERIQUES, LE',/,
     &           1X,'         COEFFICIENT DE CORIOLIS EST',/,
     &           1X,'         CALCULE EN FONCTION DE LA LATITUDE.',/,
     &           1X,'         LE MOT-CLE ''COEFFICIENT DE CORIOLIS''',/,
     &           1X,'         N''EST DONC PAS PRIS EN COMPTE.')
12        FORMAT(1X,'PROSOU : IN SPHERICAL COORDINATES, THE CORIOLIS',/,
     &           1X,'         PARAMETER DEPENDS ON THE LATITUDE.',/,
     &           1X,'         THE KEY WORD ''CORIOLIS COEFFICIENT''',/,
     &           1X,'         IS CONSEQUENTLY IGNORED.')
!
        ELSE
!
          CALL OS( 'X=X+CY  ' , FU , VN , VN ,  FCOR )
          CALL OS( 'X=X+CY  ' , FV , UN , UN , -FCOR )
!
          IF(LT.EQ.1) THEN
            IF(LNG.EQ.1) WRITE(LU,21)
            IF(LNG.EQ.2) WRITE(LU,22)
          ENDIF
21        FORMAT(1X,'PROSOU : EN COORDONNEES CARTESIENNES, LE',/,
     &           1X,'         COEFFICIENT DE CORIOLIS EST LU DANS LE',/,
     &           1X,'         FICHIER DES PARAMETRES ET CORRESPOND',/,
     &           1X,'         AU MOT-CLE ''COEFFICIENT DE CORIOLIS''',/,
     &           1X,'         IL EST ALORS CONSTANT EN ESPACE')
22        FORMAT(1X,'PROSOU : IN CARTESIAN COORDINATES, THE CORIOLIS',/,
     &           1X,'         PARAMETER IS READ IN THE STEERING FILE',/,
     &           1X,'         IT IS THE KEY WORD ''CORIOLIS',/,
     &           1X,'         COEFFICIENT'', IT IS UNIFORM IN SPACE')
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!  THE SECOND MEMBERS ARE PROPERLY DISCRETISED
!
      IELMU=UN%ELM
!
      IF(IELMU.NE.IELM1) THEN
        CALL CHGDIS(FU,IELM1,IELMU,MESH)
        CALL CHGDIS(FV,IELM1,IELMU,MESH)
      ENDIF
!
!-----------------------------------------------------------------------
!
      IELMH=HN%ELM
      CALL CPSTVC(HN,SMH)
      YASMH=.FALSE.
      CALL OS('X=0     ',X=SMH)
!
!     RAIN-EVAPORATION
!
      IF(RAIN) THEN
        RAIN_MPS=RAIN_MMPD/86400000.D0
        SURDT=1.D0/DT
        IF(BANDEC) THEN
!         EVAPORATION (TENTATIVELY...) LIMITED BY AVAILABLE WATER
          DO I=1,NPOIN
            PLUIE%R(I)=MAX(RAIN_MPS,-MAX(HN%R(I),0.D0)*SURDT)
          ENDDO
        ELSE
          CALL OS('X=C     ',X=PLUIE,C=RAIN_MPS)
        ENDIF 
      ENDIF
!
!     SOURCES
!
      IF(NREJET.GT.0) THEN
!
        YASMH = .TRUE.
!
!       SOURCE TERMS IN THE CONTINUITY EQUATION
!       BEWARE, SMH IS ALSO USED FOR TRACER
!
        DO I = 1 , NREJET
          IR = ISCE(I)
!         THE TEST IS USEFUL IN PARALLEL MODE, WHEN THE POINT SOURCE
!         IS NOT IN THE SUB-DOMAIN
          IF(IR.GT.0) THEN
            IF(OPTSOU.EQ.1) THEN
!             "NORMAL" VERSION
              SMH%R(IR)=SMH%R(IR)+DSCE(I)*UNSV2D%R(IR)
            ELSE
!             "DIRAC" VERSION
              SMH%R(IR)=SMH%R(IR)+DSCE(I)
            ENDIF
          ENDIF
        ENDDO
!
!       SOURCE TERMS IN THE MOMENTUM EQUATIONS
!       EXPLICIT TREATMENT OF MOMENTUM CONTRIBUTIONS TO THE SOURCES
!
        IF(NREJEU.GT.0) THEN
          DO I = 1 , NREJEU
            IR = ISCE(I)
!           THE TEST IS USEFUL IN PARALLEL MODE, WHEN THE POINT SOURCE
!           IS NOT IN THE SUB-DOMAIN
            IF(IR.GT.0) THEN
!             MOMENTUM ADDED BY THE SOURCE
!      -      MOMENTUM TAKEN BY THE SOURCE
              FU%R(IR)=FU%R(IR) + (VUSCE(AT,I)-UN%R(IR))*
     &        DSCE(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
              FV%R(IR)=FV%R(IR) + (VVSCE(AT,I)-VN%R(IR))*
     &        DSCE(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
            ENDIF
          ENDDO
        ENDIF
!
      ENDIF
!
!     CULVERTS
!
      IF(NSIPH.GT.0) THEN
!
        YASMH = .TRUE.
!
        DO I = 1 , NSIPH
        IR = ENTSIP(I)
          IF(IR.GT.0) THEN
            IF(OPTSOU.EQ.1) THEN
!             "NORMAL" VERSION
              SMH%R(IR)=SMH%R(IR)-DSIP(I)*UNSV2D%R(IR)
            ELSE
!             "DIRAC" VERSION
              SMH%R(IR)=SMH%R(IR)-DSIP(I)
            ENDIF
            FU%R(IR) = FU%R(IR) - (USIP(I,1)-UN%R(IR))*
     &      DSIP(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
            FV%R(IR) = FV%R(IR) - (VSIP(I,1)-VN%R(IR))*
     &      DSIP(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
          ENDIF
          IR = SORSIP(I)
          IF(IR.GT.0) THEN
            IF(OPTSOU.EQ.1) THEN
!             "NORMAL" VERSION
              SMH%R(IR)=SMH%R(IR)+DSIP(I)*UNSV2D%R(IR)
            ELSE
!             "DIRAC" VERSION
              SMH%R(IR)=SMH%R(IR)+DSIP(I)
            ENDIF
            FU%R(IR) = FU%R(IR) + (USIP(I,2)-UN%R(IR))*
     &      DSIP(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
            FV%R(IR) = FV%R(IR) + (VSIP(I,2)-VN%R(IR))*
     &      DSIP(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
          ENDIF
        ENDDO
      ENDIF
!
!     TUBES OR BRIDGES
!
      IF(NBUSE.GT.0) THEN
!
        YASMH = .TRUE.
!
        DO I = 1 , NBUSE
          IR = ENTBUS(I)
          IF(IR.GT.0) THEN
            IF(OPTSOU.EQ.1) THEN
!             "NORMAL" VERSION
              SMH%R(IR)=SMH%R(IR)-DBUS(I)*UNSV2D%R(IR)
            ELSE
!             "DIRAC" VERSION
              SMH%R(IR)=SMH%R(IR)-DBUS(I)
            ENDIF
            FU%R(IR) = FU%R(IR) - (UBUS(I,1)-UN%R(IR))*
     &      DBUS(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
            FV%R(IR) = FV%R(IR) - (VBUS(I,1)-VN%R(IR))*
     &      DBUS(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
          ENDIF
          IR = SORBUS(I)
          IF(IR.GT.0) THEN
            IF(OPTSOU.EQ.1) THEN
!             "NORMAL" VERSION
              SMH%R(IR)=SMH%R(IR)+DBUS(I)*UNSV2D%R(IR)
            ELSE
!             "DIRAC" VERSION
              SMH%R(IR)=SMH%R(IR)+DBUS(I)
            ENDIF
            FU%R(IR) = FU%R(IR) + (UBUS(I,2)-UN%R(IR))*
     &      DBUS(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
            FV%R(IR) = FV%R(IR) + (VBUS(I,2)-VN%R(IR))*
     &      DBUS(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
          ENDIF
        ENDDO
      ENDIF
!
!     WEIRS (ONLY IF TYPSEUIL=2)
!
      IF(NWEIRS.GT.0.AND.TYPSEUIL.EQ.2) THEN
!
        YASMH = .TRUE.
!
         DO N=1,NWEIRS
            DO I=1,NPSING%I(N)
               IR=NDGA1%ADR(N)%P%I(I)
               IF(IR.GT.0) THEN
                 SMH%R(IR)= SMH%R(IR)+QWA%ADR(N)%P%R(I)*UNSV2D%R(IR)
! QUANTITY OF MOVEMENTS NOT TAKING IN ACCOUNT FOR THE MOMENT
! The following lines generate instability and crash
! Probably because we would like to impose velocities accross  solid boundaries!
!
!                 FU%R(IR) = FU%R(IR) + (UWEIRA%ADR(N)%P%R(I)-UN%R(IR))*
!     &              QWA%ADR(N)%P%R(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
!                 FV%R(IR) = FV%R(IR) + (VWEIRA%ADR(N)%P%R(I)-VN%R(IR))*
!     &              QWA%ADR(N)%P%R(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
               ENDIF
               IR=NDGB1%ADR(N)%P%I(I)
               IF(IR.GT.0) THEN
                 SMH%R(IR)=SMH%R(IR)+QWB%ADR(N)%P%R(I)*UNSV2D%R(IR)
!                 FU%R(IR) = FU%R(IR) + (UWEIRB%ADR(N)%P%R(I)-UN%R(IR))*
!     &              QWB%ADR(N)%P%R(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
!                 FV%R(IR) = FV%R(IR) + (VWEIRB%ADR(N)%P%R(I)-VN%R(IR))*
!     &              QWB%ADR(N)%P%R(I)*UNSV2D%R(IR)/MAX(HN%R(IR),0.1D0)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!***********************************************************************
!
!     * WITH WAVE DRIVEN CURRENTS
!       -------------------------------------
!
!                FU        =  FXWAVE
!                  COUROU
!
!                FV        =  FYWAVE
!                  COUROU
!
!       FXWAVE AND FYWAVE ARE TAKEN IN A RESULTS FILE FROM
!       ARTEMIS OR TOMAWAC
!
!       BEWARE   : 1. MESHES MUST BE THE SAME
!       ---------
!
!                  2. STATIONARY FORCING
!
      IF(COUROU) THEN
!
!        WITH NO COUPLING, TAKING THE WAVE STRESSES ONCE FOR ALL
!        IN A BINARY DATA FILE
!
         IF(.NOT.DEJALU.AND..NOT.INCLUS(COUPLING,'TOMAWAC')) THEN
!
            ALLOCATE(W(NPOIN),STAT=ERR)
            IF(ERR.NE.0) THEN
              IF(LNG.EQ.1) THEN
                WRITE(LU,*) 'ERREUR D''ALLOCATION DE W DANS PROSOU'
                WRITE(LU,*) 'CODE ERREUR ',ERR
                WRITE(LU,*) 'NOMBRE DE POINTS : ',NPOIN
              ENDIF
              IF(LNG.EQ.2) THEN
                WRITE(LU,*) 'MEMORY ALLOCATION ERROR OF W IN PROSOU'
                WRITE(LU,*) 'ERROR CODE ',ERR
                WRITE(LU,*) 'NUMBER OF POINTS: ',NPOIN
              ENDIF
            ENDIF
!
!           NBI1 : BINARY DATA FILE 1
            NOMX='FORCE FX        '
            NOMY='FORCE FY        '
            CALL FIND_IN_SEL(FXWAVE,NOMX,T2D_FILES(T2DBI1)%LU,
     &                       T2D_FILES(T2DBI1)%FMT,W,OKX,NPTH,NP,ATH)
            CALL FIND_IN_SEL(FYWAVE,NOMY,T2D_FILES(T2DBI1)%LU,
     &                       T2D_FILES(T2DBI1)%FMT,W,OKY,NPTH,NP,ATH)
            IF(.NOT.OKX.OR..NOT.OKY) THEN
!             SECOND TRY (OLD VERSIONS OF ARTEMIS OR TOMAWAC)
              NOMX='FORCE_FX        '
              NOMY='FORCE_FY        '
              CALL FIND_IN_SEL(FXWAVE,NOMX,T2D_FILES(T2DBI1)%LU,
     &                         T2D_FILES(T2DBI1)%FMT,W,OKX,NPTH,NP,ATH)
              CALL FIND_IN_SEL(FYWAVE,NOMY,T2D_FILES(T2DBI1)%LU,
     &                         T2D_FILES(T2DBI1)%FMT,W,OKY,NPTH,NP,ATH)
            ENDIF
!           CLANDESTINE VARIABLES FROM TOMAWAC TO SISYPHE
            IF(NVARCL.GT.0) THEN
              DO I=1,NVARCL
              CALL FIND_IN_SEL(VARCL%ADR(I)%P,VARCLA(I)(1:16),
     &                         T2D_FILES(T2DBI1)%LU,
     &                         T2D_FILES(T2DBI1)%FMT,
     &                         W,OKC,NPTH,NP,ATH)
              IF(.NOT.OKC) THEN
                IF(LNG.EQ.1) WRITE(LU,7) VARCLA(I)(1:16)
                IF(LNG.EQ.2) WRITE(LU,8) VARCLA(I)(1:16)
7             FORMAT(1X,'PROSOU : VARIABLE CLANDESTINE :',/,1X,A16,/,1X,
     &                  '         NON TROUVEE',/,1X,
     &                  '         DANS LE FICHIER DE HOULE')
8             FORMAT(1X,'PROSOU : CLANDESTINE VARIABLE:',/,1X,A16,/,1X,
     &                  '         NOT FOUND',/,1X,
     &                  '         IN THE WAVE RESULTS FILE')
              CALL PLANTE(1)
              STOP
              ENDIF
              ENDDO
            ENDIF
!
            IF(.NOT.OKX.OR..NOT.OKY) THEN
              IF(LNG.EQ.1) WRITE(LU,5)
              IF(LNG.EQ.2) WRITE(LU,6)
5             FORMAT(1X,'PROSOU : FORCE FX OU FY NON TROUVES',/,1X,
     &                  '         DANS LE FICHIER DE HOULE')
6             FORMAT(1X,'PROSOU: FORCE FX OR FY NOT FOUND',/,1X,
     &                  '         IN THE WAVE RESULTS FILE')
              CALL PLANTE(1)
              STOP
            ENDIF
            IF(NP.NE.NPOIN) THEN
              IF(LNG.EQ.1) WRITE(LU,95)
              IF(LNG.EQ.2) WRITE(LU,96)
 95           FORMAT(1X,'PROSOU : SIMULATION DES COURANTS DE HOULE.',/,
     &               1X,'LES MAILLAGES HOULE ET COURANTS SONT ',/,
     &               1X,'DIFFERENTS : PAS POSSIBLE POUR LE MOMENT.')
 96           FORMAT(1X,'PROSOU: WAVE DRIVEN CURRENTS MODELLING.',/,
     &               1X,'WAVE AND CURRENT MODELS MESHES ARE ',/,
     &               1X,'DIFFERENT : NOT POSSIBLE AT THE MOMENT.')
!
              CALL PLANTE(1)
              STOP
            ENDIF
!           WRITES OUT TO THE LISTING
            IF(LNG.EQ.1) WRITE(LU,115) ATH
            IF(LNG.EQ.2) WRITE(LU,116) ATH
115         FORMAT(1X,/,1X,'PROSOU : COURANTS DE HOULE',/,
     &                  1X,'         LECTURE AU TEMPS ',F10.3,/)
116         FORMAT(1X,/,1X,'PROSOU: WAVE DRIVEN CURRENTS MODELLING',/,
     &                  1X,'         READING FILE AT TIME ',F10.3,/)
            IF(IELMU.NE.IELM1) THEN
              CALL CHGDIS(FXWAVE,IELM1,IELMU,MESH)
              CALL CHGDIS(FYWAVE,IELM1,IELMU,MESH)
            ENDIF
            DEJALU = .TRUE.
!
         ENDIF
!
!        ADDS INTO FU AND FV
!
         IF(INCLUS(COUPLING,'TOMAWAC')) THEN
           IF(IELMU.NE.IELM1) THEN
             CALL CHGDIS(FXWAVE,IELM1,IELMU,MESH)
             CALL CHGDIS(FYWAVE,IELM1,IELMU,MESH)
           ENDIF
         ENDIF
         CALL OS('X=X+Y   ',X=FU,Y=FXWAVE)
         CALL OS('X=X+Y   ',X=FV,Y=FYWAVE)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!     TAKES SEEPAGE IN THE SOIL INTO ACCOUNT
!     COMMUNICATES WITH ESTEL-3D
!
!     GETS SOURCE TERM FROM ESTEL-3D TO ACCOUNT FOR SEEPAGE
!     CALLS THE INFILTRATION ROUTINE
!
      CALL INFILTRATION_GET(SMH%R,UNSV2D%R,YASMH)
!
!-----------------------------------------------------------------------
!
      RETURN
      END
!                    ***************
                     SUBROUTINE INTERPOLATION
!                    ***************
!
!         BILINEAR INTERPOLATION 
!
!     &(X,Y,NPOIN,XRELV,YRELV,NP,L1,L2,L3,L4,NP1,NP2,NP3,NP4)
!
!         BARYCENTRIC INTERPOLATION     
!
     &(X,Y,NPOIN,XRELV,YRELV,NP,L1,L2,NP1,NP2,NP3)
!
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    CREATING A TABLE OF WEIGHTING FACTOR TO ACCELERATE THE 
!         INTERPOLATION OF DATA FROM STRUCTURAL DOMAIN (NetCDF) 
!         ON THE COMPUTATIONAL MESH NODES.
!
!history  EHSAN SARHADI ZADEH (AnteaGroup, BELGIUM)
!+        23/04/14
!+        V6P3
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| L1,L2,L3,L4    |<--| WEIGHTING FACTOR FOR INTERPOLATION
!| NP             |-->| NUMBER OF POINTS IN THE INPUT FILE
!| NP1,NP2,NP3,NP4|<--| ADDRESS OF INPUT POINTS
!| NPOIN          |-->| NUMBER OF POINTS IN THE COMPUTATIONAL MESH
!| X,Y            |-->| MESH COORDINATES
!| XRELV          |-->| ABCISSAE OF INPUT FILE
!| YRELV          |-->| ORDINATES OF INPUT FILE
!| OK             |-->| IF YES, NO CROSSING OF BOUNDARY
!|                |   | IF NO, POINT OUTSIDE THE DOMAIN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
      USE BIEF, EX_FASP => FASP
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPOIN,NP
      DOUBLE PRECISION, INTENT(IN)  :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN)  :: XRELV(NP),YRELV(NP)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER N,I
!
      DOUBLE PRECISION DIST1,DIST2,DIST3,DIST4
      DOUBLE PRECISION DIFX,DIFY,DIST,X1,Y1,X2,Y2,X3,Y3,X4,Y4

      INTEGER, INTENT(OUT)          :: NP1(NPOIN),NP2(NPOIN),NP3(NPOIN)
      
!         BILINEAR INTERPOLATION       
!      DOUBLE PRECISION, INTENT(OUT) :: L1(NPOIN),L2(NPOIN)
!      DOUBLE PRECISION, INTENT(OUT) :: L3(NPOIN),L4(NPOIN)
!      INTEGER, INTENT(OUT)          :: NP4(NPOIN)

!         BARYCENTRIC INTERPOLATION      
      DOUBLE PRECISION, INTENT(OUT) :: L1(NPOIN),L2(NPOIN)
      INTEGER, DIMENSION(NPOIN)     :: NP4
      DOUBLE PRECISION, DIMENSION(NPOIN)     :: L3
!
      LOGICAL OK1,OK2,OK3,OK4
!
!-----------------------------------------------------------------------
!
!  LOOP ON THE MESH NODES:
!
      DO I = 1 , NPOIN
!
!     INTERPOLATES THE BOTTOM FROM 4 QUADRANTS
!
! ---->  INITIALISES:
!
!     min and max of MERCATOR coordinates
!     φmax = ±85.05113° 
!
      DIST1=625.D12
      DIST2=625.D12
      DIST3=625.D12
      DIST4=625.D12
!
      OK1 = .FALSE.
      OK2 = .FALSE.
      OK3 = .FALSE.
      OK4 = .FALSE.
!
! --------->  LOOP ON THE SET OF POINTS (THERE ARE NP):
!
       DO N = 1 , NP
           DIFX = XRELV(N)-X(I)
           DIFY = YRELV(N)-Y(I)
           DIST = DIFX*DIFX + DIFY*DIFY
!
             IF ( DIST.LT.1.D-6 ) DIST=1.D-6
! ->QUADRANT 1 :
               IF( DIFX.LE.0.D0.AND.DIFY.LE.0.D0) THEN
                 IF(DIST.LE.DIST1)THEN
                      X1=XRELV(N)
                      Y1=YRELV(N)
                      NP1(I)=N
                      DIST1=DIST
                      OK1 = .TRUE.
                 ENDIF
! ->QUADRANT 2 :
              ELSE IF( DIFX.GE.0.D0.AND.DIFY.LE.0.D0) THEN
                 IF(DIST.LE.DIST2)THEN
                      X2=XRELV(N)
                      Y2=YRELV(N)
                      NP2(I)=N
                      DIST2=DIST
                      OK2 = .TRUE.
                 ENDIF
! ->QUADRANT 3 :
              ELSE IF( DIFX.GE.0.D0.AND.DIFY.GE.0.D0) THEN
                 IF(DIST.LE.DIST3)THEN
                      X3=XRELV(N)
                      Y3=YRELV(N)
                      NP3(I)=N
                      DIST3=DIST
                      OK3 = .TRUE.
                 ENDIF
! ->QUADRANT 4 :
              ELSE IF( DIFX.LE.0.D0.AND.DIFY.GE.0.D0) THEN
                 IF(DIST.LE.DIST4)THEN
                      X4=XRELV(N)
                      Y4=YRELV(N)
                      NP4(I)=N
                      DIST4=DIST
                      OK4 = .TRUE.
                 ENDIF
              ENDIF
       END DO
!
! --------->  END OF LOOP ON THE SET OF POINTS

         L1(I)=0
         L2(I)=0
!         L3(I)=0
!         L4(I)=0
!
!-----------------------------------------------------------------------
!         BILINEAR INTERPOLATION 
!         4 CORRESPONDING NODES
!-----------------------------------------------------------------------
!
!        L1(I)=((X2-X(I))*(Y3-Y(I)))/((X2-X1)*(Y3-Y1))
!        L2(I)=((X(I)-X1)*(Y3-Y(I)))/((X2-X1)*(Y3-Y1))
!        L3(I)=((X2-X(I))*(Y(I)-Y1))/((X2-X1)*(Y3-Y1))
!        L4(I)=((X(I)-X1)*(Y(I)-Y1))/((X2-X1)*(Y3-Y1))
!
!-----------------------------------------------------------------------
!         BARYCENTRIC INTERPOLATION 
!         3 CORRESPONDING NODES
!-----------------------------------------------------------------------
!
!         ELIMINATE THE FARTHEST POINT
          
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST3) X3=X4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST3) Y3=Y4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST3) NP3(I)=NP4(I)
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST2) X2=X4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST2) Y2=Y4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST2) NP2(I)=NP4(I)
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST1) X1=X4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST1) Y1=Y4
        IF (MAX(DIST4,DIST3,DIST2,DIST1)==DIST1) NP1(I)=NP4(I)

!-----------------------------------------------------------------------
        L1(I)=((Y2-Y3)*(X(I)-X3)+(X3-X2)*(Y(I)-Y3))/ 
     &        ((Y2-Y3)*(X1-X3)+(X3-X2)*(Y1-Y3))
        L2(I)=((Y3-Y1)*(X(I)-X3)+(X1-X3)*(Y(I)-Y3))/ 
     &        ((Y2-Y3)*(X1-X3)+(X3-X2)*(Y1-Y3))

!         CHECK THE POINTS TO BE INSIDE A TRIANGLE 
!-----------------------------------------------------------------------
!
!      
      L3(I)=1-L1(I)-L2(I)

      IF (L1(I).LE.0.OR.L1(I).GE.1) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'INTERPOLATION'
      IF(LNG.EQ.1) WRITE(LU,*) 'CENTREPOINT EST EXTÉRIEUR AU TRIANGLE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CENTERPOINT IS OUTSIDE THE TRIANGLE'
      WRITE(LU,*) 'POINT: ',I
      WRITE(LU,*) L1(I),L2(I),L3(I),'L1,L2,L3 ≠ 0.0 & 1.0'
      ENDIF

      IF (L2(I).LE.0.OR.L2(I).GE.1) THEN

      WRITE(LU,*) ' '
      WRITE(LU,*) 'INTERPOLATION'
      IF(LNG.EQ.1) WRITE(LU,*) 'CENTREPOINT EST EXTÉRIEUR AU TRIANGLE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CENTERPOINT IS OUTSIDE THE TRIANGLE'
      WRITE(LU,*) 'POINT: ',I
      WRITE(LU,*) L1(I),L2(I),L3(I),'L1,L2,L3 ≠ 0.0 & 1.0'

      ENDIF

!         L3 WOULD BE CALCULATED IN THE METEO CODES
!
!        L3(I)=1-L1(I)-L2(I)
!        L4(I)=0.0    

!-----------------------------------------------------------------------
!
      END DO
!
!-----------------------------------------------------------------------
!
      RETURN
      END

!                *************************
                 SUBROUTINE DATE_MJD
!                *************************
!
     &( MM,ID,IYYY,STARTTIME )
!
!***********************************************************************
! TELEMAC2D   V6P2                                   16/01/2012
!***********************************************************************
!
!brief    CONVERTS DATE TO MJD (MODIFIED JULIAN DAYS)
!+  INPUT:  ID - DAY, MM - MONTH, IYYY - YEAR
!+  OUTPUT: MJD > 0 - MODIFIED JULIAN DAYS
!+  DATE >= 11.17.1858 CORRESPONDS TO MJD = 0
!
!history  EHSAN SARHADI ZADEH (AnteaGroup, BELGIUM)
!+        23/04/2014
!+        V6P3
!+   Modification for comparing simulation date with netCDF date
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|  ID            |-->| DAY
!|  IYYY          |-->| YEAR
!|  MM            |-->| MONTH
!|  STARTTIME     |<--| SIMULATION TIME IN MINUTES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: MM,ID,IYYY
      INTEGER, INTENT(OUT) :: STARTTIME
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER DPM(12),DAYS,I,NLEAP,K
      DATA DPM/31,28,31,30,31,30,31,31,30,31,30,31/
!
!-----------------------------------------------------------------------
!
      STARTTIME = 0
!     NO EARLIER DATES THAN NOVEMBER 17TH 1858
      IF(     IYYY.LT.1858.OR.(IYYY.EQ.1858.AND.MM.LT.11)
     &   .OR.(IYYY.EQ.1858.AND.MM.EQ.11.AND.ID.LT.17) ) THEN
         IF(LNG.EQ.1) WRITE(LU,*) 'PAS DE DATES ANTERIEURES ' //
     &                'AU 17 NOVEMBRE 1858 NE SONT PERMISES'
         IF(LNG.EQ.2) WRITE(LU,*) 'NO EARLIER DATES ' //
     &                'THAN NOVEMBER 17TH 1858 ARE ALLOWED'
         CALL PLANTE(1)
         STOP
      ENDIF
!
      DAYS = 0
      DO I = 1,MM-1
         DAYS = DAYS+DPM(I)
         IF( I.EQ.2.AND.INT(IYYY/4)*4.EQ.IYYY ) DAYS = DAYS+1
      ENDDO
!     321TH DAY CORRESPONDS TO NOVEMBER 17TH FOR A NON LEAP YEAR
      DAYS = DAYS+ID-321

!     LEAP DAY CORRECTION
      DO K = 1900,IYYY,100
         IF( K.EQ.IYYY.AND.MM.GT.2 ) DAYS = DAYS-1
      ENDDO
      DO K = 2000,IYYY,400
         IF( K.EQ.IYYY.AND.MM.GT.2 ) DAYS = DAYS+1
      ENDDO
!     EACH 4TH YEAR IS LEAP YEAR
      NLEAP = INT(REAL(IYYY-1-1860)*0.25)
      IF( IYYY.GT.1860 ) NLEAP = NLEAP+1
!     EXCEPT
      DO K = 1900,IYYY-1,100
        IF( K.LT.IYYY ) NLEAP = NLEAP-1
!     THE FOLLOWING LINE IS USELESS AS K.GE.IYYY-1
!       IF( K.EQ.IYYY.AND.MM.GT.2 ) DAYS = DAYS-1
      ENDDO
!     BUT EACH IN THE ROW 2000:400:... IS LEAP YEAR AGAIN
      DO K = 2000,IYYY-1,400
        IF( K.LT.IYYY ) NLEAP = NLEAP+1
!     THE FOLLOWING LINE IS USELESS AS K.GE.IYYY-1
!       IF( K.EQ.IYYY.AND.MM.GT.2 ) DAYS = DAYS+1
      ENDDO
      STARTTIME = 365*(IYYY-1858)+NLEAP+DAYS
!
!-----------------------------------------------------------------------
!
      RETURN
      END


MODULE DWM_Wake_Sub

USE   NWTC_Library
USE   InflowWind
USE   Wind, ONLY: RHO
USE   Element
USE   BLADE

USE   DWM_Wake_main_Data

IMPLICIT NONE

    !......... Public Subroutines .................................................

PUBLIC :: CalVelScale
PUBLIC :: turbine_average_velocity
PUBLIC :: pass_velocity
PUBLIC :: filter_average_induction_factor
PUBLIC :: calculate_mean_u
PUBLIC :: calculate_element_area
PUBLIC :: calculate_induction_factor
PUBLIC :: CalAtmShear
PUBLIC :: get_initial_condition
PUBLIC :: calculate_wake
PUBLIC :: Gauss
PUBLIC :: TI_downstream_total
PUBLIC :: smallscale_TI
PUBLIC :: shear_correction
PUBLIC :: AD_GetUndisturbedWind
PUBLIC :: filter_velocity
PUBLIC :: Get_wake_center
PUBLIC :: shifted_velocity
PUBLIC :: read_parameter_file
PUBLIC :: read_upwind_result_file
PUBLIC :: write_result_file
PUBLIC :: turbsim_mean_velocity
PUBLIC :: read_turbine_position
PUBLIC :: smooth_out_wake
PUBLIC :: smooth_wake_shifted_velocity
!PUBLIC :: write_wake_center
PUBLIC :: rename_FAST_output
PUBLIC :: advection_time
PUBLIC :: DWM_init
PUBLIC :: Hub_height_mean_velocity
PUBLIC :: collect_velocity

    CONTAINS
!----------------------------------------------------------------------------------
SUBROUTINE CalVelScale(AtmU,AtmV)
!..................................................................................
! This routine is to calculat the atmospheric length scale before 
! introducing the TI term (which will be used later)
!..................................................................................
    USE DWM_Wake_main_Data, ONLY:AtmUscale   

    REAL               ::    AtmU             ! atmospheric U velocity
    REAL               ::    AtmV             ! atmospheric V velocity
    INTEGER,SAVE       ::    counter = 0
    REAL,SAVE          ::    Denominator = 0
    REAL,SAVE          ::    Numerator = 0
    
    counter = counter + 1
    
    Denominator =  (Denominator * (counter-1) + AtmU * AtmU) / counter
    Numerator   =  (Numerator   * (counter-1) + AtmU * AtmV) / counter
    
    AtmUscale = Numerator / Denominator
    
        
END SUBROUTINE CalVelScale

!----------------------------------------------------------------------------------
SUBROUTINE turbine_average_velocity( single_velocity, blade_num, element, average_velocity_array_local )
!..................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! To calculate the average of the wind speed of a specific blade ring
! the outpout is the average_velocity_array
!..................................................................................
    
    REAL                ::   single_velocity
    REAL,ALLOCATABLE    ::   average_velocity_array_local(:)
    INTEGER             ::   element
    INTEGER             ::   blade_num
    INTEGER             ::   i
    INTEGER, SAVE       ::   time_step_velocity = -1
    INTEGER,ALLOCATABLE,SAVE   :: time_step_velocity_array(:) ! counter of each section of the blade
    
    time_step_velocity = time_step_velocity +1
    
    IF ( time_step_velocity==0) THEN
       ALLOCATE ( average_velocity_array_local(NElm) )
       ALLOCATE ( time_step_velocity_array(NElm) )
       !ALLOCATE ( average_velocity_array_temp(NumElOut) )
       average_velocity_array_local(:)       = 0
       time_step_velocity_array(:)           = 0
       !average_velocity_array_temp           = 0
       average_velocity_array_local(element) = single_velocity
       time_step_velocity_array(element)     = 1
       
    ELSE IF (time_step_velocity > 0) THEN
       DO i=1,NElm
          IF ( element == i) THEN
             time_step_velocity_array(element)     = time_step_velocity_array(element)+1
             average_velocity_array_local(element) = (average_velocity_array_local(element)*( (time_step_velocity_array(element)-1) ) + single_velocity) &
                                               / (time_step_velocity_array(element))
             IF ( i == NElm .and. blade_num == NB) THEN
                CALL pass_velocity(average_velocity_array_local)
                !average_velocity_array_temp(:) = average_velocity_array_local(:)              ! the average velocity of a whole time step
                time_step_velocity = -1
                IF (ALLOCATED( average_velocity_array_local ))      DEALLOCATE ( average_velocity_array_local )
                IF (ALLOCATED( time_step_velocity_array ))          DEALLOCATE ( time_step_velocity_array )
                !IF (ALLOCATED( average_velocity_array_temp ))       DEALLOCATE ( average_velocity_array_temp )
             END IF                
          END IF
       END DO
    END IF
    
END SUBROUTINE turbine_average_velocity
!----------------------------------------------------------------------------------
SUBROUTINE pass_velocity(array_velocity)
!..................................................................................
! 
! 
!  
!..................................................................................
    USE  turbine_average_velocity_data,  ONLY: average_velocity_array_temp
    
    
    REAL                  ::  array_velocity(NElm)
    INTEGER,SAVE          ::  time_step_pass_velocity = -1
    
    time_step_pass_velocity = time_step_pass_velocity+1
    
    IF (time_step_pass_velocity==0) THEN
       ALLOCATE (average_velocity_array_temp(NElm))
       average_velocity_array_temp(:) = array_velocity(:)
    ELSE IF(time_step_pass_velocity>0) THEN
       average_velocity_array_temp(:) = array_velocity(:)
    END IF
        
END SUBROUTINE pass_velocity
!----------------------------------------------------------------------------------
SUBROUTINE filter_average_induction_factor( thrust_force, num_of_element, dr_blade )
!..................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! The output of the subroutine is induction factor at this time step 
! and the average induction factor through all the time steps which have been simulated
!.................................................................................. 
    USE turbine_average_velocity_data
    
    
    INTEGER                      ::  num_of_element                   ! The number of the nodes in the blade
    REAL                         ::  thrust_force ( num_of_element,NB )  ! Thrust force at each node    N/m
    INTEGER, SAVE                ::  time_step_force = -1                   ! The time step (save attribute) of the FAST simulation
    REAL                         ::  thrust_coefficient ( num_of_element )
    REAL                         ::  average_induction_factor   ( num_of_element )
    REAL                         ::  induction_factor_local_temp ( num_of_element )
    REAL                         ::  dr_blade (num_of_element)
    INTEGER                      ::  I,J
    
    
    time_step_force = time_step_force +1
          
    IF ( time_step_force==0) THEN
       ALLOCATE (induction_factor ( num_of_element ))
       ALLOCATE (average_velocity_array ( num_of_element ))
       ALLOCATE (turbine_thrust_force (num_of_element ))
       ALLOCATE (swept_area (num_of_element ))
       
       turbine_thrust_force = 0
       DO I = 1,num_of_element
           DO J = 1,NB
               turbine_thrust_force (I) = turbine_thrust_force (I) + thrust_force (I,J)      !dFn/m
           END DO
       END DO
       
       DO I = 1,num_of_element
           turbine_thrust_force (I) = turbine_thrust_force(I)  !! * dr_blade(I)           ! integrate dFn through blade
       END DO
                
       !turbine_thrust_force (:) = NB * thrust_force(:)
       CALL calculate_element_area ( R, NElm, RELM(:), swept_area )
       CALL calculate_induction_factor ( turbine_thrust_force , swept_area , num_of_element, average_velocity_array_temp, induction_factor_local_temp )
       induction_factor = induction_factor_local_temp
       average_velocity_array = average_velocity_array_temp
       
    ELSE IF ( time_step_force>0) THEN
       turbine_thrust_force = 0
       DO I = 1,num_of_element
           DO J = 1,NB
               turbine_thrust_force (I) = turbine_thrust_force (I) + thrust_force (I,J)  !dFn/m
           END DO
       END DO
       
       DO I = 1,num_of_element
           turbine_thrust_force (I) = turbine_thrust_force(I)   !!* dr_blade(I)           ! integrate dFn through blade
       END DO
       
       !turbine_thrust_force (:) = NB * thrust_force(:)
       
       CALL calculate_induction_factor ( turbine_thrust_force , swept_area , num_of_element, average_velocity_array_temp, induction_factor_local_temp )
       induction_factor = ( induction_factor(:) * time_step_force + induction_factor_local_temp(:) ) / ( time_step_force+1 )
       average_velocity_array = ( average_velocity_array(:) * time_step_force + average_velocity_array_temp(:) ) / ( time_step_force+1 )
    END IF
    
END SUBROUTINE filter_average_induction_factor

!----------------------------------------------------------------------------------
SUBROUTINE calculate_mean_u( num_element,r_t,turbine_mean_velocity,TI_normalization )
!..................................................................................
! This routine is called to calculate the mean velocity and the TI of the turbine
! Using weighting method according to the blade ring area 
!.................................................................................. 
    USE turbine_average_velocity_data, ONLY : average_velocity_array
    USE read_turbine_position_data   , ONLY : SimulationOrder_index,upwindturbine_number
    USE parameter_file_data             , ONLY : TI_amb,Uambient,TI_wake,Uambient
    USE weighting_method             , ONLY : sweptarea,weighting_denominator
    USE read_upwind_result_file_data , only : upwind_TI,upwind_meanU
    
    INTEGER           ::    num_element
    INTEGER           ::    i
    REAL              ::    r_t ( num_element )              ! The distance from the node to the hub
    REAL              ::    turbine_mean_velocity            ! turbine mean velocity
    REAL              ::    node_radius    ( num_element )
    REAL              ::    element_length ( num_element )
    REAL              ::    TI_normalization

    

    ALLOCATE (sweptarea(num_element))
    
    weighting_denominator = 0
    turbine_mean_velocity = 0
    
    element_length (num_element) = 2.0*( R - r_t(num_element) )
    DO i=num_element-1,1,(-1)
       element_length(i)= 2.0*( r_t(i+1)-r_t(i) ) - element_length (i+1)
    END DO

    node_radius ( num_element ) = R - element_length (num_element)
    DO i=num_element-1,1,(-1)
       node_radius (i) = node_radius (i+1) - element_length (i)
    END DO

    DO i=1, num_element-1,1
       sweptarea (i) = Pi * (node_radius (i+1) **2 - node_radius (i) **2)
    END DO
    sweptarea (num_element) = Pi * (R**2-node_radius (num_element)**2) ! ring area
    
    DO i=1,num_element
       weighting_denominator = weighting_denominator + sweptarea (i)   ! denominator
    END DO
    
    ! calculate the mean velocity of the turbine using weighting method
    DO i=1,num_element
       turbine_mean_velocity = turbine_mean_velocity + sweptarea (i) / weighting_denominator * average_velocity_array(i)       
    END DO
    
    IF (SimulationOrder_index == 1 .OR. SimulationOrder_index == 0) THEN
       TI_normalization = TI_amb
    ELSE
        IF(upwindturbine_number /= 0) THEN           ! superimpose the TI from upstream wakes
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            ! superimpose the TI from several upstream turbines
            TI_normalization  =  0
            DO I = 1,upwindturbine_number               
               TI_normalization = (TI_normalization**2 + upwind_TI(I)**2)**0.5
            END DO
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           TI_normalization = upwind_TI(1)          ! only take the TI effect from the closest upstream turbine
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           TI_normalization = TI_normalization/(turbine_mean_velocity/Uambient)
           !TI_normalization = TI_normalization/(turbine_mean_velocity/upwind_meanU(1))
        ELSE
           TI_normalization = TI_amb
        END IF
    END IF
    
    !TI_normalization = TI_amb/(turbine_mean_velocity/Uambient)   ! 8/5/2014 test
    
    !OPEN (unit=25,file="D:\5MW_simulation\after_release\results\TI.txt")
    
    !OPEN (unit=25,file="DWM\results\TI.txt")
    !WRITE (25,'(f13.7)'), TI_normalization 
     
END SUBROUTINE calculate_mean_u

!----------------------------------------------------------------------------------------
SUBROUTINE calculate_element_area (blade_radius, num_element, r_t, sweptarea)
!..................................................................................
! This routine is called when the Aerodyn simuilation finishes.
! This routine is called to calculate the swept area of each blade section.
! The output of the subroutine is swept_area (:). 
! Which is the swept area of of each blade element.
!..................................................................................
    IMPLICIT                                NONE
    
    
    INTEGER         :: num_element
    INTEGER         :: i
    INTEGER         :: j
    INTEGER         :: k
    REAL            :: blade_radius
    REAL            :: r_t ( num_element )              ! The distance from the node to the hub
    REAL            :: node_radius    ( num_element )
    REAL            :: element_length ( num_element )
    INTEGER         :: time_step_area
    REAL            :: sweptarea(num_element)

    element_length (num_element) = 2.0*( blade_radius - r_t(num_element) )
    DO i=num_element-1,1,(-1)
       element_length(i)= 2.0*( r_t(i+1)-r_t(i) ) - element_length (i+1)
    END DO

    node_radius ( num_element ) = blade_radius - element_length (num_element)
    DO k=num_element-1,1,(-1)
       node_radius (k) = node_radius (k+1) - element_length (k)
    END DO

    DO j=1, num_element-1,1
       sweptarea (j) = Pi * (node_radius (j+1) **2 - node_radius (j) **2)
    END DO
    sweptarea (num_element) = Pi * (blade_radius**2-node_radius (num_element)**2)
    
    
END SUBROUTINE calculate_element_area

!----------------------------------------------------------------------------------
SUBROUTINE calculate_induction_factor ( normal_force , element_swept_area , num_element, FFWS_array, induction_factor_local )
!..................................................................................
! This routine is called to calculate thrust coefficient then the induction factor using local thrust force.
! The output of the subroutine is induction_factor (:). 
! Which is the induction factor of each blade node.
!..................................................................................
    IMPLICIT                                NONE
    
    
    INTEGER   ::   num_element
    INTEGER   ::   i
    INTEGER   ::   j
    INTEGER   ::   ErrStat
    REAL      ::   temp_wind(3)
    REAL      ::   normal_force ( num_element )
    REAL      ::   element_swept_area ( num_element )
    REAL      ::   thrust_coefficient ( num_element )
    REAL      ::   induction_factor_local   ( num_element )
    REAL      ::   FFWS_array( num_element)
    REAL      ::   Ct_1
    REAL      ::   a_t
    REAL      ::   Ct_critical

    !Calculate thrust coefficient
    DO i=1,num_element
      thrust_coefficient (i) = normal_force(i)/(0.5* RHO* element_swept_area(i)* FFWS_array(i)**2)   ! previously use VROTORX
    END DO
    
    ! Then calculate the induction factor by solving 4a(1-a)= Ct
    ! Applying the Glauert empirical Ct modification (10.7.2013)
      
    Ct_1       = 1.816
    a_t        = 1 - 0.5*SQRT(Ct_1)
    Ct_critical = 4*a_t*(1-a_t)
    
    DO j=1,num_element
         IF (thrust_coefficient(j)<=Ct_critical) THEN
            induction_factor_local (j) = (-4 + (16-16*thrust_coefficient(j))**(0.5))/(2*(-4))
         ELSE
            induction_factor_local (j) = 1 - (thrust_coefficient(j)-Ct_1) / (-4*(SQRT(Ct_1)-1))
         END IF
    END DO

END SUBROUTINE calculate_induction_factor

!----------------------------------------------------------------------------------
SUBROUTINE CalAtmShear ( Uscale, DuDzAbl, TIAmb)
!..................................................................................
! This routine is called to calculate the atmospheric shear gradient du/dz
! The output is the DuDzAbl
!..................................................................................
    USE BLADE,   ONLY : R
    
    REAL        ::        Uscale
    REAL        ::        DuDzAbl
    REAL        ::        TIAmb
    REAL        ::        Von_Karman
    REAL        ::        Height
    REAL        ::        Ustar_scale
    REAL        ::        L_scale
    
    Ustar_scale = ( TIAmb * TIAmb * ABS(Uscale) )**0.5
    
    Von_Karman = 0.41         ! Von Karman constant
    Height     = 100
    
    L_scale = Von_Karman * Height / R
    
    !DuDzAbl = Ustar_scale / L_scale
    DuDzAbl = 0.5
   
END SUBROUTINE CalAtmShear
!----------------------------------------------------------------------------------
SUBROUTINE get_initial_condition( induc_array, r_t, element_num, r_w, U_w )
!..................................................................................
! This routine is called at the end of the subroutine calculate_initial_condition.
! This routine is called to calculate the initial condition of the DWM model.
! The output of the subroutine is r_wake (:) and U_wake (:).
! Which are the the scaled rotor radius and the scaled velocity at the rotor.
!..................................................................................
    USE     read_turbine_position_data,   ONLY: SimulationOrder_index
    USE     parameter_file_data,          ONLY: smoothed_wake,Uambient
    USE     DWM_Wake_main_Data,           ONLY: Mean_FFWS,avg_ct
    USE     read_turbine_position_data,   ONLY: upwindturbine_number
    USE     read_upwind_result_file_data, ONLY: upwind_smoothWake
    USE     weighting_method,             ONLY: sweptarea,weighting_denominator
    USE     InitCond,                     ONLY: NacYaw
    
    INTEGER           :: element_num
    INTEGER           :: i,J
    REAL              :: induc_array(element_num)
    REAL              :: r_t(element_num)
    REAL              :: dA (element_num-1)
    REAL              :: a_cellC (element_num-1)
    REAL              :: mean_a
    REAL              :: f_w
    REAL              :: fU !fU factor (realised induction for wake depth) {0-1}
    REAL              :: fR !fR factor (realised expansion for wake width) {0-1} 
    REAL, ALLOCATABLE :: r_w(:)
    REAL, ALLOCATABLE :: U_w(:)
    !REAL              :: avg_induction_factor   ! average induction factor of the rotor plane
    !REAL              :: avg_ct
    

    ALLOCATE       (r_w(element_num))
    ALLOCATE       (U_w(element_num))
    fU             = 1.10           !1.10 is not working when induction factor is reaching 0.5: (1-a*(fu+1))<1     !!!! 0.92
    fR             = 0.98
    
         !open (unit=25,file="DWM\results\inducion_factor.txt")
         !write (25,'(f13.7)'), induc_array(:)
    
   !-------------------------------------------------------------------------------------------------
   ! apply the smoothed wake profile as the input wind profile for downstream turbine
   !-------------------------------------------------------------------------------------------------
    !IF (turbine_index/=1 .and. smooth_flag==1 ) THEN
      ! DO i=1,element_num,1
         ! !induc_array(i) = ( 1 - (1-2*induc_array(i))*smoothed_wake(i) )/2
       !   induc_array(i) = ( 1 - (1-2*induc_array(i))*smoothed_wake(i)*(Mean_FFWS/Uambient) )/2
       !END DO
    !ELSE
     !  induc_array = induc_array
    !END IF
    
    !IF (turbine_index/=1) THEN
     ! DO i=1,element_num,1
      !   induc_array(i) = (1-Mean_FFWS*(1-2*induc_array(i))/Uambient)/2
      !END DO
    !ELSE
      !induc_array = induc_array
    !END IF

   !-------------------------------------------------------------------------------------------------
   ! calculate the boundary condition
   !-------------------------------------------------------------------------------------------------
   
    !===Initial Condition===
    DO i=1,element_num-1
       dA (i) =( (r_t(i+1)/R)**2 - (r_t(i)/R)**2 )*Pi
    END DO
    !== simple approximation of cell center value
    DO i=1,element_num-1
       a_cellC (i) = ( induc_array(i+1) +induc_array(i) ) /2
    END DO

    !== Boundary Condition r_w
    mean_a = 0.0
    DO i=1,element_num-1
      mean_a = ( a_cellC (i) * dA (i) ) /Pi + mean_a
    END DO

    !== Uniform expansion
    f_w = ( (1-mean_a) / (1- ((1+fR) * mean_a)) )**0.5
    r_w = r_t /R * f_w

    !== Boundary velocity
    !U_w = 1 - (induc_array * (1+fU))
    
    ! superimpose the smoothed wake from upwind turbines
    IF (SimulationOrder_index == 1 .OR. SimulationOrder_index == 0) THEN
        U_w = 1 - (induc_array * (1+fU))
    ELSEIF(SimulationOrder_index > 1) THEN
        IF (upwindturbine_number == 0) THEN
            U_w = 1 - (induc_array * (1+fU))
        ELSEIF (upwindturbine_number > 0) THEN
            !!ALLOCATE (smoothed_wake(element_num))
            !!smoothed_wake = 1
            !!DO I = 1,upwindturbine_number
              !!  DO J = 1,element_num
                !!    smoothed_wake(J) = 1- ( (1-smoothed_wake(J))**2 + (1-upwind_smoothWake(I,J))**2 )**0.5
                !!END DO
            !!END DO
            
            DO I = 1,element_num
                !U_w(i) = (Uambient/Mean_FFWS) * upwind_smoothWake(1,I)*(1 - (induc_array(i) * (1+fU)))
                U_w(i) = (Mean_FFWS/Uambient) * (1 - (induc_array(i) * (1+fU)))
            END DO
        END IF
    END IF
    
    DO i=1,element_num,1           ! modification for low wind speed, high thrust situation
        IF (U_w(i) < 0.0) THEN
            U_w(i) = 0.01
        END IF
    END DO
    
    
    !! IF (turbine_index/=1 .and. smooth_flag==1 ) THEN
       !!  DO i=1,element_num,1
          !! !U_w = (Mean_FFWS/Uambient) * smoothed_wake(i)*(1 - (induc_array * (1+fU)))
           !! U_w = (Uambient/Mean_FFWS) * smoothed_wake(i)*(1 - (induc_array * (1+fU)))
           !! !U_w =  smoothed_wake(i)*(1 - (induc_array * (1+fU)))
           !! !U_w = 1 - (induc_array * (1+fU))
        !! END DO
    !! ELSE
        !! U_w = 1 - (induc_array * (1+fU))
    !! END IF
    
    ! calculate the average induction factor of the rotor plane   
    !avg_induction_factor = 0
    
    !DO i=1,element_num
     !  avg_induction_factor = avg_induction_factor + sweptarea (i) / weighting_denominator * induc_array(i)       
    !END DO
    
    !skew_angle = 0.60*avg_induction_factor*NacYaw *(-1)    ! minus sign means different direction
    
    
    ! calculate the average thrust coefficient and the ct_tilde
    
    avg_ct = 0
    
    DO i=1,element_num
       avg_ct = avg_ct + sweptarea (i) / weighting_denominator * ( 4*induc_array(i)*(1-induc_array(i)) )       
    END DO
    
    !ct_tilde  = 0.5*COS(NacYaw)**2*SIN(NacYaw)*avg_ct
    
    ct_tilde  = avg_ct
    
    ! Calculate the averaged induction factor  6.24.2014
    
    MeanInd = 0
    
    DO i=1,element_num
       MeanInd = MeanInd + sweptarea (i) / weighting_denominator * induc_array(i)       
    END DO
    
    ! for advection speed test
         OPEN (unit=25,file='DWM-results\average_Ct.txt')
         WRITE (25,'(f13.7)'), avg_ct
         WRITE (25,'(f13.7)'), sweptarea
         WRITE (25,'(f13.7)'), weighting_denominator

END SUBROUTINE get_initial_condition

!-------------------------------------------------------------------------
SUBROUTINE calculate_wake(r_w, U_w, element_num, U, b)
!..................................................................................
! This routine is the main routine to calculate the wake
! This routine is called after receiving the scaled rotor radius and the scaled velocity at the rotor
! The output of this routine is the wake velocity which is "U"
!   and the wake width which is the "b"
!..................................................................................
    USE   DWM_Wake_Deficit_Data
    USE   parameter_file_data, ONLY: p_p_r, r_domain, Uambient, x_domain, hub_height, TI_amb
    USE   read_turbine_position_data, ONLY:upwindturbine_number
    
    IMPLICIT                            NONE
    
    INTEGER              ::    element_num
    REAL                 ::    r_w (element_num)        ! scaled rotor radius r_w
    REAL                 ::    U_w (element_num)        ! scaled velocity  U_w
    REAL,ALLOCATABLE     ::    U(:,:)
    INTEGER,ALLOCATABLE     ::    b(:)
    
    REAL   ::  mtemp
    REAL   ::  ntemp
    REAL   ::  xtemp
    REAL   ::  ytemp
    
    real   ::  temp_test1
    real   ::  temp_test2
   
    
    
    
    ppR       = p_p_r
    Domain_R  = r_domain       !10.0  domain size in R [R]
    Domain_X  = x_domain      !42.0  domain size in X [R]
    filter1   = (/0.0,   4.0 /)
    filter2   = (/0.035, 0.35/)
    k1        = 0.0919 !0.1519 !0.0919   ! 7.31.2014 test
    k2        = 0.0178 !0.0108 !0.0178
    R_WTG     = 1.0
    U0        = 1.0
    D_WTG     = 2.0
    R_length  = Domain_R
    X_length  = Domain_X
    
    TI_original = TI*(Mean_FFWS/Uambient)    ! calculate the TI if under ambient wind speed
    !TI_original = TI
    
    Print*, TI_original
    
    np_r       =  ppR          ! per R   ie. R resolution is 50
    np_x       =  ppR          ! per D   ie. X resolution is 50
    delrad     =  R_WTG/np_r   ! dr
    delaxi     =  D_WTG/np_x   ! dx
    n_x_vector =  floor((X_length)/D_WTG*np_x)      ! number of point in equally spaced array x_vector
    n_r_vector =  floor((R_length)/R_WTG*np_r)      ! number of point in equally spaced array r_vector

    ! create coordinate vectors
    ALLOCATE    (x_vector(n_x_vector))                             
    ALLOCATE    (r_vector(n_r_vector))
    ! similar to linspace function
    x_vector   = ( (X_length-delaxi)/(n_x_vector-1 ) )*[(i,i=1,n_x_vector)]+(0-( (X_length-delaxi)/(n_x_vector-1) ))        
    r_vector   = ( (R_length-delrad)/(n_r_vector-1 ) )*[(i,i=1,n_r_vector)]+(0-( (R_length-delrad)/(n_r_vector-1) )) 
    
    ! Create the F1 filter
    
    CALL create_F1_filter (F1_vector, filter1, length_F1_vector,np_x,X_length)
         !OPEN (unit=25,file="DWM\results\F1_filter.txt")
         !WRITE (25,'(f13.6)'), F1_vector(:)
         !CLOSE(25)
         
    CALL create_F2_filter (F2_vector, filter2, np_x, length_F1_vector)
    

    !CALL create_filter_vector ( filter1, F1_vector )
    
    !CALL create_filter_vector ( filter2, F2_vector )
    
         !OPEN (unit=25,file="DWM\results\F2_filter.txt")
         !WRITE (25,'(f13.6)'), F2_vector(:)
         !CLOSE(25)
  
     ! Initiate the U, V, visc, TI_add, Turb_stress matrices
    ALLOCATE (V               (n_x_vector,n_r_vector))
    ALLOCATE (U               (n_x_vector,n_r_vector))    !axial velocity matrix
    ALLOCATE (visc            (n_x_vector,n_r_vector))
    ALLOCATE (visc_DWM        (n_x_vector,n_r_vector))
    ALLOCATE (du_dr_DWM       (n_x_vector,n_r_vector))
    ALLOCATE (du_dr_total     (n_x_vector,n_r_vector))
    ALLOCATE (Turb_Stress_DWM (n_x_vector,n_r_vector))
    !ALLOCATE (TI_DWM          (n_x_vector,n_r_vector))
    !ALLOCATE (U_face          (n_x_vector,n_r_vector-1))
    !ALLOCATE (VOL_x_jhigh     (n_x_vector,n_r_vector-1))
    !ALLOCATE (VOL_x_jlow      (n_x_vector,n_r_vector-1))
    !ALLOCATE (VOL_r_ihigh     (n_x_vector,n_r_vector-1))
    !ALLOCATE (VOL_r_ilow      (n_x_vector,n_r_vector-1))
    ALLOCATE (visc_wake         (n_x_vector,n_r_vector))
    ALLOCATE (visc_wake1        (n_x_vector,n_r_vector))
    ALLOCATE (visc_wake2        (n_x_vector,n_r_vector))
    ALLOCATE (alfa_1            (n_r_vector           ))
    ALLOCATE (alfa_2            (n_r_vector           ))
    ALLOCATE (du_dr_tot         (n_x_vector,n_r_vector))
    ALLOCATE (shear_flag        (n_r_vector           )) 
    ALLOCATE (One_div_du_dr_DWM (n_x_vector,n_r_vector))
    ALLOCATE (visc_fac          (n_r_vector           ))
    
    ALLOCATE (main_diagonal     (n_r_vector           ))
    ALLOCATE (sub_diagonal      (n_r_vector           ))
    ALLOCATE (sup_diagonal      (n_r_vector           ))
    

    V                  = 0
    U                  = 0
    visc               = 0
    du_dr_DWM          = 0
    Turb_Stress_DWM    = 0
    !TI_DWM             = 0
    !U_face             = 0
    !VOL_x_jhigh        = 0
    !VOL_x_jlow         = 0
    !VOL_r_ihigh        = 0
    !VOL_r_ilow         = 0

    !%%%% BOUNDARY CONDITIONS
    ! ROTOR PLANE
    ALLOCATE (bin_filter(n_r_vector))
    DO i=1,n_r_vector
      IF (MAXVAL(r_w)>r_vector(i)) THEN
      bin_filter(i) = 1
      ELSE
      bin_filter(i) = 0
      END IF
    END DO


    n_xi=floor(sum(bin_filter))
    ALLOCATE (xi(n_xi))
    xi=r_vector(1:n_xi)*bin_filter(1:n_xi)

    ALLOCATE (U_tmp_1(n_xi))
    ILo = 1
    DO i=1,n_xi
       U_tmp_1(i) = InterpBin( xi(i), r_w, U_w, ILo, size(r_w))
    END DO

    n_U_tmp_2 = size(r_vector)-n_xi 
    ALLOCATE (U_tmp_2(n_U_tmp_2))
    U_tmp_2 = U0

    ALLOCATE (U_tmp(n_xi +n_U_tmp_2))    !  =n_r_vector
    U_tmp = (/U_tmp_1, U_tmp_2/)
    U (1,:) = U_tmp(1:n_r_vector)

    ! Centerline
    V (1,:) = 0
    

       
    !%%%% SOLVING FLOW FIELD
    ALLOCATE  (b(n_x_vector))
    ALLOCATE  (counter(n_x_vector))
              counter=1
    ALLOCATE (AugMat (n_r_vector,n_r_vector+1) )
    ALLOCATE (Soln   (n_r_vector)              )

    ALLOCATE (mat(n_r_vector,n_r_vector))
    ALLOCATE (RHS(n_r_vector))

    n_r_vec_DWM = floor(ppR*Domain_R)
    ALLOCATE ( r_vec_DWM (n_r_vec_DWM) )
    ALLOCATE ( dA_DWM (n_r_vec_DWM-1)  )
    
    !%%%%%%%%%%%%%%%  Atmospheric stability effects  %%%%%%%%%%%%%%%%%%
    L_ABL_vector         = (/26.5352,      34.026,      40.7458/)
    UW_UU_vector         = (/-0.27359,    -0.27887,    -0.27935/)
    L_DEF_vector         = (/11.065,       12.9746,     14.4395/)
    UU_DEF_UU_ABL_vector = (/0.63044,      0.57982,      0.5287/)
    UW_DEF_UU_DEF_vector = (/-0.27341,    -0.25684,    -0.24217/)
    
    !%%%%% interpolation wrt to the hub height
    x_ary = (/40,100,160/)
    L_ABL         = InterpBin( hub_height, x_ary, L_ABL_vector, ILo, size(x_ary) ) !int_Lww(k3) i.e (Integral length scale (in vertical directions), from ww(k3)) 
    UW_UU         = InterpBin( hub_height, x_ary, UW_UU_vector, ILo, size(x_ary) ) ! ratio of UW and UU stresses for whole spectra   
    L_DEF         = InterpBin( hub_height, x_ary, L_DEF_vector, ILo, size(x_ary) ) ! Integral length scale (in vertical directions), Meandering length scale subtracted, from ww(k3)   
    UU_DEF_UU_ABL = InterpBin( hub_height, x_ary, UU_DEF_UU_ABL_vector, ILo, size(x_ary) ) ! Part of normal stress in the deficit module   
    UW_DEF_UU_DEF = InterpBin( hub_height, x_ary, UW_DEF_UU_DEF_vector, ILo, size(x_ary) ) ! ratio of UW and UU stresses for spectra in deficit scales
    
    !%%%%% normalized by the fixed rotor R
    Rotor_fixed_R  = 40 !ATMOSTAB ANALYSIS IS CARRIED OUT OVER R = 40m, which should be used to normalize the length scales
    l_star_ABL     = L_ABL / Rotor_fixed_R;
    l_star_DEF     = L_DEF / Rotor_fixed_R;
    
    !%%%%% Normalize UU_160m to neutral condition
    UU_DEF_UU_ABL_fac = InterpBin( hub_height, x_ary, (/0.63044, 0.57982, 0.5287/), ILo, size(x_ary) )
    UU_DEF_UU_ABL     = UU_DEF_UU_ABL / UU_DEF_UU_ABL_fac
    
    !%%%%% CALCULATE u* according to:
    ! 1. u* ~= (mean(u'w')^2 )^0.25 
    ! 2. {mean(u'w') = mean(u'u')*Cuw_uu} 
    ! 3. {u' ~= TI (in normalized form)}
    ! => u* ~= ((TI^2 * Cuw_uu )^2)^0.25 
    u_star_ABL     = ( (  (TI_amb/100)**2                      * ABS(UW_UU)         )**2 )**0.25
    u_star_DEF     = ( (  (TI_original/100)**2 * UU_DEF_UU_ABL * abs(UW_DEF_UU_DEF) )**2 )**0.25;
    
    Shear_add_du_dz = u_star_ABL / l_star_ABL
    
    !k_wiener = du_dz_ABL + delrad**2

    DO j=2, n_x_vector, 1  ! start from the plane next to the rotor plane

    !==== Calculating wake width "b" where 95% of the deficit is captured
      dr_DWM = 1.0/ppR
  
       DO i=1,n_r_vec_DWM    ! build r_vec_DWM
          r_vec_DWM (i) = dr_DWM/2 + (i-1)*dr_DWM
       END DO
   
       DO i=1,n_r_vec_DWM-1    ! build dA_DWM
          dA_DWM (i) = Pi * ( r_vec_DWM (i+1)**2 - r_vec_DWM (i)**2 )
       END DO
   
       Def_DWM = 0                             ! Calculate Def_DWM 
       DO i=1, n_r_vec_DWM-1
          Def_DWM = (1-U (j-1,i+1)) * dA_DWM (i) + Def_DWM
       END DO
   
       Def_DWM_mixL = 0                        ! Calculate the wake width "b"
       DO i = 2,ppR
          Def_DWM_mixL = ( 1- U(j-1,i) ) * dA_DWM(i) + Def_DWM_mixL
       END DO
       DO b_counter = ppR+1, n_r_vec_DWM-1          
          Def_DWM_mixL = ( 1- U(j-1,b_counter) ) * dA_DWM(b_counter) + Def_DWM_mixL
          IF ( Def_DWM_mixL > Def_DWM * 0.999 ) THEN
             EXIT
          ELSE IF (b_counter == n_r_vec_DWM-1) THEN
             EXIT
          END IF
       END DO
       b(j-1) = b_counter
 
       ! %%%%% Calculate eddy viscosity
       ! Include blend between original Prandtl model and Ainslie to avoid issues when wake turbulence goes to 0.
       ! The largest eddy viscosity at each point is applied.
       
       ! Calculate mean flow gradient - du/dr is created with CDS (apart from 1st and last point)
       du_dr_DWM(j-1,1)                  = (U(j-1,2) - U(j-1,1))/delrad
       du_dr_DWM(j-1,2:R_length*np_r-1)  = (U(j-1,3:(R_length*np_r-1)+1) - U(j-1,1:(R_length*np_r-1)-1))/(2*delrad)
       du_dr_DWM(j-1,R_length*np_r)      = (U(j-1,R_length*np_r) - U(j-1,R_length*np_r-1))/delrad
       
       ! %%% Blend of mixL and Ainslie eddy visc
       DO I = 1,n_r_vector
           visc_wake1(j-1,I)     = F2_vector(j-1)* k2 *( r_vector(b(j-1))/R_WTG )**2 * ABS(du_dr_DWM(j-1,I));
           visc_wake2(j-1,I)     = F2_vector(j-1)* k2 *( r_vector(b(j-1))/R_WTG )    * ( 1 - min_of_array( U(j-1,:),SIZE(U(j-1,:)) ) );
           visc_wake (j-1,I)     = max_of_TwoNum( visc_wake1(j-1,I),visc_wake2(j-1,I) );
       END DO
       
       ! %%% Atmospheric eddy visc as u*l*, yields total eddy viscosity
       visc_norm_factor = 6.3918
       DO I = 1,n_r_vector
           visc(j-1,I)           = F1_vector(j-1) * k1 * visc_norm_factor * u_star_DEF * l_star_DEF + visc_wake(j-1,I);
       END DO
       
       ! %%%%% Include contribution from atmospheric boundary layer on DWM
       ! % 1. Calculate the azimuthally averaged local gradient (du/dr tot) acting of the eddy viscosity as a combination of du/dr in the DWM model and du/dz from ABL
       ! % 2. The du/dr contribution is constant in azimuthal direction. The du/dz part is assumed linear, which gives a sinus curve in a du/dr system
       
       !% Calculate total mean flow gradient - adds shear contribution via
       !% sin function. This gets the stresses right, but sign is wrong in
       !% regions where du/dr_DWM - sign of du/dz_ABL is negative 
       
       DO I = 1,n_r_vector
          !alfa_1(I)      = ASIN(ABS(du_dr_DWM(j-1,I)) / Shear_add_du_dz)
          !alfa_2(I)      = Pi - alfa_1(I)
          
          ! % condition for added shear gradient (if du/dr_DWM >= du/dz_ABL there are no contribution)
          IF ( ABS(du_dr_DWM(j-1,I)) < Shear_add_du_dz ) THEN
              shear_flag(I) = 1
              alfa_1(I)      = ASIN(ABS(du_dr_DWM(j-1,I)) / Shear_add_du_dz)
              alfa_2(I)      = Pi - alfa_1(I)
          ELSE
              shear_flag(I) = 0
              alfa_1(I)      = 0
              alfa_2(I)      = 0
          END IF
          
          
          du_dr_tot(j-1,I) = (  ABS(du_dr_DWM(j-1,I))*2*Pi + shear_flag(I)*2*&
                                ( Shear_add_du_dz*2*COS(alfa_1(I))-ABS(du_dr_DWM(j-1,I))*(alfa_2(I) - alfa_1(I))) )/(2*Pi)
          temp_test1 = ABS(du_dr_DWM(j-1,I))
          temp_test2 = du_dr_tot(j-1,I)
       END DO
       
       ! %%% Use "wiener filter" for numerical stability:  1/f(x) ~= f(x) / (f(x)^2 + k)
       k_wiener                 = 2*Shear_add_du_dz * delrad**2;
       DO I = 1,n_r_vector
           One_div_du_dr_DWM(j-1,I) = du_dr_DWM(j-1,I) / (du_dr_DWM(j-1,I)**2 + k_wiener)
           visc_fac(I)              = max_of_TwoNum(1.0, (du_dr_tot(j-1,I) * ABS(One_div_du_dr_DWM(j-1,I))))
           visc(j-1,I)              = visc(j-1,I) * visc_fac(I)
       END DO
       
          
       !IF (upwindturbine_number > 0) THEN
           !visc = visc*0.9
       !END IF
       

       mat=0
       ! ====SHORT INSTRUCIONS TO SOLVE RUTINE:
       ! The terms LHS and RHS (left/right hand side) refers to the terms of
       ! the coefficient matrix developed to solve the then shear layer
       ! approximation of NS. The numbers indicate the position in the equation,
       ! ex LHS21 is the 2nd part of the 1st term on the left side in eq.2.8,
       ! see document "Numerical implementation of DWM deficit module" for details.
   
       ! Input BC for wake center
   
       LHS2       = U(j-1,1)/delaxi     + (2*visc(j-1,1)/(delrad**2))
       LHS3       = -(2*visc(j-1,1) /(delrad**2))
       RHS(1)     = (U(j-1,1)**2    / delaxi)
       mat(1,1)   = LHS2
       mat(2,1)   = LHS3
   
       ! Calculation of U for the wake body
       DO i=2,(n_r_vector-1),1      ! starts from the point next to the hub center
          LHS11             = -V(j-1,i)      / (2*delrad)
          LHS21             = visc(j-1,i)    / (2*r_vector(i)*delrad)
          LHS31             = -visc(j-1,i)   / (delrad**2)
          LHS41             = (visc(j-1,i+1) - visc(j-1,i-1))  / (2*delrad)**2;  ! new term due to d(nu_t)/dr dependence
          LHS12             = U(j-1,i)       / (delaxi)
          LHS22             = 2*visc(j-1,i)  / (delrad**2)
          LHS13             = V(j-1,i)       / (2*delrad)
          LHS23             = -visc(j-1,i)   / (2*r_vector(i)*delrad)
          LHS33             = -visc(j-1,i)   / (delrad**2)
          LHS43             = -(visc(j-1,i+1) - visc(j-1,i-1))  / (2*delrad)**2; ! new term due to d(nu_t)/dr dependence
          LHS1              = LHS11 + LHS21 + LHS31 + LHS41
          LHS2              = LHS12 + LHS22
          LHS3              = LHS13 + LHS23 + LHS33 + LHS43
          RHS(i)            = (U(j-1,i)**2  / delaxi)
          ! Build the matrix for X =A/B
          mat(i-1,i) = LHS1
          mat(i  ,i) = LHS2
          mat(i+1,i) = LHS3

       END DO
   
       ! Input BC for wake edge
       LHS1                     = 0
       LHS2                     = 1/delaxi
       RHS(R_length*np_r)       = (U(j-1,R_length*np_r)/ delaxi)
       mat(R_length*np_r-1,    R_length*np_r)     = LHS1
       mat(R_length*np_r  ,    R_length*np_r)     = LHS2
   
       ! Solve for the U
       ! Use Gauss-Jordan elimination
       AugMat (1:n_r_vector, 1:n_r_vector) = TRANSPOSE(mat)
       AugMat (:           , n_r_vector+1) = RHS
       NumEqu                              = n_r_vector
       !CALL Gauss(AugMat, NumEqu, Soln)
       !U(j,:)=Soln
       
       
       ! === USE Thomas Algorithm to solve the matrix ====  6.30.2014
       main_diagonal (1) = AugMat(1,1)
       sub_diagonal  (1) = 0               ! means it is the diagonal below the main diagonal
       sup_diagonal  (1) = AugMat(1,2)     ! means it is the diagonal above the main diagonal
      
       DO I = 2,n_r_vector-1
           main_diagonal (I) = AugMat(I,I)
           sub_diagonal  (I) = AugMat(I,I-1)
           sup_diagonal  (I) = AugMat(I,I+1)
       END DO
       
       main_diagonal (n_r_vector) = AugMat(n_r_vector, n_r_vector)
       sub_diagonal  (n_r_vector) = AugMat(n_r_vector, n_r_vector-1)
       sup_diagonal  (n_r_vector) = 0
       
       CALL Thomas_diagonal (sub_diagonal, main_diagonal, sup_diagonal, RHS, Soln, NumEqu)
       U(j,:)=Soln

       ! === Solve for V
   
       DO i = 1,R_length*np_r-1,1
         V(j,i+1) = (r_vector(i) / r_vector(i+1)) * V(j,i) -(delrad/(2*delaxi))*( (U(j,i+1) - U(j-1,i+1)) + &
                (r_vector(i) / r_vector(i+1)) * ((U(j,i) - U(j-1,i))) )
       END DO
   
       ! POST PROCESSING SIGNAL: Turbulent stress
       DO i=1,n_r_vector,1
         !Turb_Stress_DWM(j-1,i) = visc_DWM(j-1,i) * du_dr_total(j-1,i)
         Turb_Stress_DWM(j-1,i) = visc(j-1,i) * du_dr_DWM(j-1,i)
       END DO
   
       ! Control calculatoins of mass flux over cells
  
       !!DO i=1,n_r_vector-1,1
          !! VOL_x_jhigh(j-1,i)  = (Pi/3) *((U(j,i) *(r_vector(i+1)**3 - (3 *r_vector(i)**2 *r_vector(i+1)) + 2 *r_vector(i)**3)) +&
          !!                    (U(j,i+1) *(r_vector(i)**3 - (3*r_vector(i+1)**2 *r_vector(i)) + 2 *r_vector(i+1)**3)))/ delrad
          !! VOL_x_jlow(j-1,i)   = (Pi/3) *((U(j-1,i) *(r_vector(i+1)**3 - (3*r_vector(i)**2 *r_vector(i+1)) + 2 *r_vector(i)**3)) +&
            !!                   (U(j-1,i+1) *(r_vector(i)**3   - (3*r_vector(i+1)**2 *r_vector(i)) + 2 *r_vector(i+1)**3)))/ delrad
          !! VOL_r_ilow(j-1,i)   =  Pi * r_vector(i) * (V(j-1,i)+V(j,i)) * delaxi
      
          !! V(j,i+1)            = ((VOL_x_jlow(j-1,i) - VOL_x_jhigh(j-1,i) + VOL_r_ilow(j-1,i)) / (Pi*(r_vector(i+1)) *delaxi)) - V(j-1,i+1)     !! changed to version 2 2012/4/11
          !! VOL_r_ihigh(j-1,i)  = Pi * r_vector(i+1) * (V(j-1,i+1)+V(j,i+1)) * delaxi
          !! specificly U_face for mass flow and momentum calculations
          !! U_face(j-1,i) = VOL_x_jlow(j-1,i)  / (Pi*((r_vector(i+1)**2)-(r_vector(i))**2))
       !!END DO
   
    END DO

    b(n_x_vector) = b(n_x_vector-1)
    
    IF (ALLOCATED( V ))                  DEALLOCATE ( V )
    IF (ALLOCATED( visc ))               DEALLOCATE ( visc )
    IF (ALLOCATED( du_dr_DWM ))          DEALLOCATE ( du_dr_DWM )
    !IF (ALLOCATED( Turb_Stress_DWM ))    DEALLOCATE ( Turb_Stress_DWM )
    !IF (ALLOCATED( TI_DWM ))             DEALLOCATE ( TI_DWM )
    !IF (ALLOCATED( U_face ))             DEALLOCATE ( U_face )
    !IF (ALLOCATED( VOL_x_jhigh ))        DEALLOCATE ( VOL_x_jhigh )
    !IF (ALLOCATED( VOL_x_jlow ))         DEALLOCATE ( VOL_x_jlow )
    !IF (ALLOCATED( VOL_r_ihigh ))        DEALLOCATE ( VOL_r_ihigh )
    !IF (ALLOCATED( VOL_r_ilow ))         DEALLOCATE ( VOL_r_ilow )
    IF (ALLOCATED( r_vec_DWM ))          DEALLOCATE ( r_vec_DWM )
    IF (ALLOCATED( dA_DWM ))             DEALLOCATE ( dA_DWM )
    IF (ALLOCATED( bin_filter ))         DEALLOCATE ( bin_filter )
    IF (ALLOCATED( xi ))                 DEALLOCATE ( xi )
    IF (ALLOCATED( U_tmp_1 ))            DEALLOCATE ( U_tmp_1 )
    IF (ALLOCATED( U_tmp_2 ))            DEALLOCATE ( U_tmp_2 )
    IF (ALLOCATED( U_tmp ))              DEALLOCATE ( U_tmp )
    IF (ALLOCATED( mat ))                DEALLOCATE ( mat )
    IF (ALLOCATED( U_tmp ))              DEALLOCATE ( U_tmp )
    IF (ALLOCATED( RHS ))                DEALLOCATE ( RHS )
    IF (ALLOCATED( Soln ))               DEALLOCATE ( Soln )
    IF (ALLOCATED( AugMat ))             DEALLOCATE ( AugMat )
    
    !OPEN(unit = 10, status='replace',file='sizeof_Uvelocity_2nd.bin',form='unformatted')   ! create sizeof_Uvelocity_2nd.bin, or overwrite an existing on
    !WRITE(10)   ppR,Domain_R                                                               ! write the length of the velocity vector                                                                                                                                                                      
    !CLOSE(10)
    
    !OPEN(unit = 10, status='replace',file='DWM\results\Uvelocity.bin',form='unformatted')          
    !WRITE(10)   U(floor(spacing_turbine * ppR)+1,:)  ! write the wind data of the plane where the downstream turbine locates                                                                                                                                                                                                                                                                                                    
    !CLOSE(10)
    
    !OPEN (unit=25,file="DWM\results\wake_width.txt")
    !WRITE (25,'(I5)'), b(:)
    !close(25)

    
END SUBROUTINE calculate_wake

!---------------------------------------------------------------------------------------------
SUBROUTINE Thomas_diagonal (lowerDia, mainDia, upperDia, RightHS, SolnVec, NumEq)
!.............................................................................................
! This function returns the F1 filter function
!.............................................................................................

    INTEGER  ::   NumEq
    
    REAL     ::   lowerDia(NumEq)
    REAL     ::   mainDia(NumEq)
    REAL     ::   upperDia(NumEq)
    REAL     ::   RightHS(NumEq)
    REAL     ::   SolnVec(NumEq)
    
    REAL     ::   cp_vec(NumEq)
    REAL     ::   dp_vec(NumEq)
    REAL     ::   temp
    INTEGER  ::   I
    
    ! initialize c-prime and d-prime
    cp_vec(1) = upperDia(1) / mainDia(1)
    dp_vec(1) = RightHS(1)  / mainDia(1)
    
    ! solve for vectors c-prime and d-prime
    DO I = 2,NumEq
        temp = mainDia(i) - cp_vec(i-1)*lowerDia(i)
        cp_vec(i) = upperDia(i)/temp
        dp_vec(i) = (RightHS(i)-dp_vec(i-1)*lowerDia(i))/temp
    END DO
    
    ! initialized SolnVec
    SolnVec(NumEq) = dp_vec(NumEq)
    
    ! solve for x from the vectors c-prime and d-prime
    DO I = NumEq-1, 1,-1
        SolnVec(i) = dp_vec(i) - cp_vec(i)*SolnVec(i+1)
    END DO
    
END SUBROUTINE Thomas_diagonal

!---------------------------------------------------------------------------------------------
SUBROUTINE create_F1_filter (F1_vector, filter1, length_F1_vector,np_x,X_length)
!.............................................................................................
! This function returns the F1 filter function
!.............................................................................................
    REAL,ALLOCATABLE   ::     F1_vector(:)
    REAL               ::     filter1(2)
    INTEGER            ::     length_F1_vector
    INTEGER            ::     np_x
    REAL               ::     X_length
    
    INTEGER            ::     length_F1_vector_1
    INTEGER            ::     length_F1_vector_2
    REAL,ALLOCATABLE   ::     F1_vector_1(:)
    REAL,ALLOCATABLE   ::     F1_vector_2(:)
    INTEGER            ::     I
    
    length_F1_vector_1 = floor(filter1(2)*np_x/2)
    length_F1_vector_2 = floor(X_length*np_x/2)
    length_F1_vector   = length_F1_vector_1 + length_F1_vector_2
    ALLOCATE    (F1_vector_1(length_F1_vector_1))
    ALLOCATE    (F1_vector_2(length_F1_vector_2))
    ALLOCATE    (F1_vector  (length_F1_vector  ))
    
    F1_vector_1 = ( (1-filter1(1))   /(length_F1_vector_1-1 ) )*[(i,i=1,length_F1_vector_1)]+(0-( (1-filter1(1) )/(length_F1_vector_1-1 ) ))
    !r_vector    = ( (R_length-delrad)/(n_r_vector-1 ) )*[(i,i=1,n_r_vector)]+(0-( (R_length-delrad)/(n_r_vector-1) ))
    F1_vector_2 = 1
    
    F1_vector   = (/F1_vector_1,F1_vector_2/)
  
END SUBROUTINE create_F1_filter

!---------------------------------------------------------------------------------------------
SUBROUTINE create_F2_filter (F2_vector, filter2, np_x, length_F1_vector)
!.............................................................................................
! This function returns the F2 filter function
!.............................................................................................
    REAL,ALLOCATABLE   ::    F2_vector(:)
    REAL               ::    filter2(2)
    INTEGER            ::    np_x
    INTEGER            ::    length_F1_vector
    
    REAL,ALLOCATABLE   ::    F2_vector_x(:)
    REAL,ALLOCATABLE   ::    F2_vector_1(:)
    REAL,ALLOCATABLE   ::    F2_vector_2(:)
    INTEGER            ::    length_F2_vector_x
    INTEGER            ::    length_F2_vector_1
    INTEGER            ::    length_F2_vector_2
    INTEGER            ::    length_F2_vector
    INTEGER            ::    I
    
    length_F2_vector_x = floor(( REAL(length_F1_vector) * (1/REAL(np_x)) - (2+1/REAL(np_x)) ) / (1/REAL(np_x)) + 1)
    length_F2_vector_1 = 2*np_x
    length_F2_vector_2 = length_F2_vector_x
    length_F2_vector   = length_F2_vector_1 + length_F2_vector_2
    
    ALLOCATE ( F2_vector_x(length_F2_vector_x) )
    ALLOCATE ( F2_vector_1(length_F2_vector_1) )
    ALLOCATE ( F2_vector_2(length_F2_vector_2) )
    ALLOCATE ( F2_vector  (length_F2_vector  ) )
    
    F2_vector_x = ( (length_F1_vector * (1/REAL(np_x)) - (2+1/REAL(np_x)))   /(length_F2_vector_x-1 ) )*[(i,i=1,length_F2_vector_x)]+(2+1/REAL(np_x)-( (length_F1_vector * (1/REAL(np_x)) - (2+1/REAL(np_x)))   /(length_F2_vector_x-1 ) ))
    F2_vector_1 = filter2(1)
    
    DO I = 1,length_F2_vector_2
        F2_vector_2(I) = 1-(1-filter2(1))*EXP(-filter2(2)*(F2_vector_x(I)-2))
    END DO
    
    F2_vector   = (/F2_vector_1,F2_vector_2/)
  
END SUBROUTINE create_F2_filter

!----------------------------------------------------------------------------------
SUBROUTINE Gauss( AugMatIn, NumEq, SolnVec )
!..................................................................................
! This routine uses the Gauss-Jordan elimination method for the solution of
!   a given set of simultaneous linear equations.
! NOTE: this routine works if no pivot points are zero and you don't want 
!   the eschelon or reduced eschelon form of the augmented matrix. The form of
!   the original augmented matrix IS preserved in this call.
!..................................................................................
    IMPLICIT                        NONE


       ! Passed variables:

    INTEGER(4), INTENT(IN )      :: NumEq                                           ! Number of equations in augmented matrix.

    REAL(ReKi), INTENT(IN )      :: AugMatIn (NumEq,( NumEq + 1 ))                  ! Augmented matrix passed into this subroutine.
    REAL(ReKi), INTENT(OUT)      :: SolnVec  (NumEq)                                ! Solution vector.


       ! Local variables:

    REAL(ReKi)                   :: AugMat   (NumEq,( NumEq + 1 ))                  ! A copy of the augmented matrix.

    INTEGER(4)                   :: I                                               ! Steps through columns
    INTEGER(4)                   :: J                                               ! Steps through rows
    INTEGER(4)                   :: L                                               ! Steps through rows
    INTEGER(4)                   :: NAug                                            ! Column dimension of augmented matrix



       ! Transfer the data from AugMatIn to AugMat:

    AugMat = AugMatIn


       ! Find the column dimension of the augmented matrix:

    NAug = NumEq + 1


       ! Perform Gauss-Jordan elimination and store the solution vector
       !   in the last column of the augmented matrix:

    DO L = 1,NumEq             ! Loop through all rows
       DO I = ( L + 1 ), NAug  ! Loop through all columns above current row number
          AugMat(L,I) = AugMat(L,I) / AugMat(L,L)
          DO J = 1,NumEq       ! Loop through all rows except L
             IF ( J /= L )  AugMat(J,I) = AugMat(J,I) - ( AugMat(J,L)*AugMat(L,I) )
          ENDDO                ! J - All rows except L
       ENDDO                   ! I - All columns above current row number
    ENDDO                      ! L - All rows


       ! Transfer the solution vector from AugMat() to SolnVec():

    SolnVec = AugMat(:,NAug)



    RETURN
END SUBROUTINE Gauss

!---------------------------------------------------------------------------------------------
FUNCTION shear_correction (du_dr_dwm,du_dz)
!.............................................................................................
! This function returns the shear correction factor A1 and A2
!.............................................................................................
    
    USE shear_correction_data
    
    REAL   ::    shear_correction
    REAL   ::    du_dr_dwm
    REAL   ::    du_dz             ! du_dz_abl
    
    alpha_1 = ASIN(du_dr_dwm/du_dz)
    alpha_2 = Pi/2 - alpha_1
    temp_integration = 0
    
    temp_n = 100
    ALLOCATE ( alpha_array (temp_n) )
    
    alpha_array = ((alpha_2-alpha_1)/(temp_n-1))*[(i,i=1,temp_n)]+(alpha_1-((alpha_2-alpha_1)/(temp_n-1)))
    delta_alpha = (alpha_2-alpha_1)/(temp_n-1)
    
    DO i = 1,temp_n
        temp_integration = du_dz * SIN(alpha_array(i)) * delta_alpha + temp_integration
    END DO
    
    shear_correction = temp_integration - (alpha_2 - alpha_1) * du_dr_dwm

    
    IF (ALLOCATED( alpha_array ))                DEALLOCATE ( alpha_array )
    
END FUNCTION shear_correction

!---------------------------------------------------------------------------------------------
FUNCTION AD_GetUndisturbedWind ( Time, InputPosition, ErrStat)
!.............................................................................................
! This function returns the U-V-W wind speeds at the specified time and X-Y-Z location
!.............................................................................................

       ! Passed variables
    USE  SharedInflowDefns
    
    REAL(ReKi), INTENT(IN)   :: Time
    REAL(ReKi), INTENT(IN)   :: InputPosition(3)
    INTEGER,    INTENT(OUT)  :: ErrStat

       ! function definition

    REAL(ReKi)               :: AD_GetUndisturbedWind(3)

       ! local variables

    TYPE(InflIntrpOut)       :: InflowVel


    !-------------------------------------------------------------------------------------------------
    ! get the wind speed (the wind inflow module will check that it's initialized)
    !-------------------------------------------------------------------------------------------------

    InflowVel = WindInf_GetVelocity( Time, InputPosition, ErrStat)

    IF (ErrStat /=0) CALL ProgWarn( ' Error getting velocity in AeroDyn/AD_GetUndisturbedWind().' )

    AD_GetUndisturbedWind(:) = InflowVel%Velocity(:)

END FUNCTION AD_GetUndisturbedWind

!-------------------------------------------------------------------------------
FUNCTION filter_velocity (timestep,y_0,z_0,wake_radius)
!...............................................................................
! This function is called to calculate the filtered wake velocity
! The filter is a low pass filter
! The output is the filtered wake velocity at a certain wake center
!...............................................................................
    USE DWM_Wake_Deficit_Data, ONLY: ppR
    USE filter_velocity_data
    USE parameter_file_data, ONLY: WFLowerBd
    

    REAL         ::  timestep                     ! upper limit = usable time + grid width /  mean wind speed              %%% will change wrt wind speed  
                                                  ! = second / 0.05    timestep >= 1
    REAL         ::  y_0                          ! wake center point
    REAL         ::  z_0
    INTEGER      ::  wake_radius                  ! b(:) in cal_mixl
    REAL         ::  filter_velocity (3)          ! only v,w components


    !Print*, y_0
    
    temp_filter_velocity = 0.0
    number_counter       = 0
    !radius_length = NINT( 2*wake_radius/ppR*R )             ! R(m): turbine radius
    
    radius_length = NINT(2*50/ppR*R )                  ! 2D filter size (2*50= 1D)

    DO y_axis = y_0-radius_length,y_0+radius_length,1
        !IF (y_axis > WFLowerBd) THEN
            
      DO z_axis = z_0-radius_length,z_0+radius_length,1
        IF ( z_axis > WFLowerBd )  THEN              !(make sure the circle does not exceed wind field)      1.12.2015
          IF ( ((y_axis-y_0)**2+(z_axis-z_0)**2)**0.5 <= radius_length )  THEN
                                
          temp_filter_velocity(:) = temp_filter_velocity(:)  +  AD_GetUndisturbedWind ( (REAL(timestep,ReKi)), (/0.0,&             
                                                                                       REAL(y_axis,ReKi),REAL(z_axis,ReKi)/), ErrStat)
                                                        !   AD_GetUndisturbedWind ( (REAL(timestep,ReKi)/20.0)-315.0, (/0.0,&             
                                                                                   !REAL(y_axis,ReKi),REAL(z_axis,ReKi)/), ErrStat)
           number_counter = number_counter + 1
          END IF
        END IF
      END DO
        !END IF        
    END DO

    filter_velocity (1) = temp_filter_velocity(1) / number_counter
       
    filter_velocity (2) = temp_filter_velocity(2) / number_counter                 ! Filtered V velocity in the certain radius circle 
    
    filter_velocity (3) = temp_filter_velocity(3) / number_counter                 ! Filtered W velocity in the certain radius circle

END FUNCTION filter_velocity 

!---------------------------------------------------------------------------------
SUBROUTINE Get_wake_center ( wakewidth, wake_center )
!................................................................................
! This routine is called to calculate the wake center of a specific release time 
!   and flying time wind plane.
! The wake center is passed to the filter to calculate the averaged wind velocity for
!   the downstream turbine.
!.................................................................................
    USE DWM_Wake_Deficit_Data, ONLY : n_x_vector,ppR 
    USE meandering_data 
    USE parameter_file_data, ONLY : WakePosition_1, WakePosition_2, hub_height, Uambient
    USE InitCond,                     ONLY: NacYaw
    USE Ptfm_pitch_data,    ONLY: colNum, ptfm_pitch 
   
       ! local variables

    REAL, ALLOCATABLE     ::   wake_center (:,:,:)
    INTEGER, ALLOCATABLE    ::   wakewidth(:)
    REAL       ::      Modified_U
    REAL       ::      x_step
    
    REAL       ::      art_scale_factor
    
    art_scale_factor = 1
    
    !-------------------------------------------------------------
    !!n_x_vector = 1700
    !!ALLOCATE  (wakewidth(n_x_vector)) 
    !!wakewidth = 50
    !!ppR = 50
    !---------------------------------------------------------------
    
    U_factor  =  1.00

    Modified_U = Mean_FFWS * U_factor
    
    DWM_time_step = (2*R/ppR)/Modified_U          ! resolution (126m/50) / wind speed (8m/s) => make sure there is always a wake width at every time step
                                                  ! D/(DWM_time_step*Mean_FFWS)= 50 which is the X resolution

    U_Scale_Factor =  Modified_U / (Uambient*U_factor)       ! modify the wake displacement error caused by the change of Mean_FFWS 
    
    U_Scale_Factor = 1   ! 2015.7.15
    
    simulation_time_length = WakePosition_1    !80   in reality, 80*scale_factor*DWM_time_step
                                 ! from 1 to 800 (scale_factor : 800/80=10)         to 16D (16*50)
    moving_time       = WakePosition_2         !50   from 0 to 49  0: wind turbine plane
                                 !               ppR/scale_factor = 5 presents 1D
                                  
                
    release_time      = simulation_time_length       
    flying_time       = moving_time       
    scale_factor      = 10       ! to decrease the calculation time
    ALLOCATE (wake_center (release_time,flying_time+1,3) )
                                 ! ex. @8D: (1~release_time,8*[ppR/scale_factor]+1,:)

    DO release_time = 1,simulation_time_length,1               ! wake center position at turbine plane
    wake_center (release_time,1,1) = 0
    wake_center (release_time,1,2) = 0
    wake_center (release_time,1,3) = REAL(hub_height)
    END DO

    ! get the initial wake center position of each cross scetion  (from the velocity at the turbine plane * dt)
    x_step = Modified_U * (DWM_time_step*scale_factor)
    DO release_time=1,simulation_time_length,1                        
       wake_center (release_time,2,1) = Modified_U * (DWM_time_step*scale_factor) +0             


       temp_center_wake (:) = AD_GetUndisturbedWind ( (REAL(((release_time-1)+1)*DWM_time_step*scale_factor,ReKi)), (/0.0,&                  
                                                   REAL(0,ReKi),REAL(hub_height,ReKi)/), ErrStat)                                ! get the velocity at the turbine plane                                     
       wake_center (release_time,2,2) = temp_center_wake (2) * (DWM_time_step*scale_factor) * U_Scale_Factor + wake_center (release_time,1,2) + & 
                                        local_skew_angle(NacYaw, ct_tilde, wake_center (release_time,2,1), NINT(ppR), ppR) * x_step !+ &
                                        !skew_lateral_offset(NacYaw, ct_tilde, wake_center (release_time,2,1), 2*RotorR) + &
                                        !rotation_lateral_offset( wake_center (release_time,2,1) )
       
       wake_center (release_time,2,3) = temp_center_wake (3) * (DWM_time_step*scale_factor) * U_Scale_Factor + wake_center (release_time,1,3) + &
                                        local_ptfm_pitch_angle(ptfm_pitch(release_time,1), ct_tilde, wake_center (release_time,2,1), NINT(ppR), ppR) * x_step
    END DO


    DO flying_time = 2,moving_time,1
       DO release_time = 1,simulation_time_length,1
       wake_center (release_time,flying_time+1,1) = wake_center (release_time,flying_time+1-1,1) + Modified_U * (DWM_time_step*scale_factor)
   
       temp_velocity(:) = filter_velocity (REAL(((release_time-1)+1)*DWM_time_step*scale_factor,ReKi), wake_center (release_time,flying_time+1-1,2), &
                                          wake_center (release_time,flying_time+1-1,3), wakewidth((flying_time-1)*scale_factor) )
   
       wake_center (release_time,flying_time+1,2) = art_scale_factor * 1.00 * temp_velocity (2) * (DWM_time_step*scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,2) + &
                local_skew_angle(NacYaw, ct_tilde, wake_center (release_time,flying_time,1), wakewidth((flying_time-1)*scale_factor), ppR) * x_step  !+ &
                                                    !rotation_lateral_offset( wake_center (release_time,flying_time+1,1) )                            - &
                                                    !rotation_lateral_offset( wake_center (release_time,flying_time,  1) )
       
       wake_center (release_time,flying_time+1,3) = art_scale_factor * 1.00 * temp_velocity (3) * (DWM_time_step*scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,3) + &
                local_ptfm_pitch_angle(ptfm_pitch(release_time,1), ct_tilde, wake_center (release_time,flying_time,1), wakewidth((flying_time-1)*scale_factor), ppR) * x_step
   
   
       END DO
    END DO
    
    
    !OPEN(unit = 10, status='replace',file='sizeof_meandered_wake_center.bin',form='unformatted')    
    !WRITE(10)   simulation_time_length, moving_time                                                                                                                                                                                                                                                                                                              
    !CLOSE(10)
    
    !OPEN(unit = 10, status='replace',file='DWM\results\meandered_wake_center.bin',form='unformatted')    
    !WRITE(10)   wake_center(:,:,:)                                                ! write the downstream meandered wake center                                                                                                                                                                      
    !CLOSE(10)
    
    !OPEN(unit = 10, status='replace',file='DWM\results\upstream_meanU.bin',form='unformatted')    
    !WRITE(10)   Mean_FFWS                                                             ! write the Mean_FFWS                                                                                                                                                                       
    !CLOSE(10)
        
END SUBROUTINE Get_wake_center

!---------------------------------------------------------------------------------
SUBROUTINE Get_wake_center_RW (wake_center)
!................................................................................
! This routine is called to calculate the wake center by using RANDOM WALK MODEL
! Have the same output format as the original Get_wake_center model
!.................................................................................
    USE meandering_data 
    USE parameter_file_data, ONLY : WakePosition_1, WakePosition_2, hub_height
    USE RW_data
    
    REAL, ALLOCATABLE     ::   wake_center (:,:,:)   ! 1st dimension: at each second, 2nd dimension: WakePosition_2+1, 3rd dimension: x(fake), y, and z
    INTEGER :: maxVal = 0
    INTEGER :: minVal = 0
    
    
    use_nth   = 40
    npoint    = WakePosition_1 * use_nth
    mean      = 0.0
    sd        = 240.00
    which     = 1
    step_size = sd/80
    scale_factor = 10
    
    IF (.NOT. ALLOCATED(arr))           ALLOCATE (arr(npoint))
    IF (.NOT. ALLOCATED(wake_center))   ALLOCATE (wake_center(WakePosition_1, WakePosition_2+1, 3))
    
    arr(1)    = 0.0
    total     = 0
    CALL random_seed()
    
    DO i=2,npoint
        cur = arr(i-1)
        CALL cdfnor(which, cdf, cdf_right, cur, mean, sd, status, bound) ! calculate the CDF
        CALL random_number (randNum) 
        IF (randNum <= cdf)   THEN   ! will move to left
            CALL random_number (randNum)
            IF (randNum<=0.1)   THEN
                step = 0.0
            ELSE
                step = step_size
            END IF
            arr(i) = arr(i-1) - step
        ELSE
            CALL random_number (randNum)
            IF (randNum<=0.1)   THEN
                step = 0.0
            ELSE
                step = step_size
            END IF
            arr(i) = arr(i-1) + step
        END IF
    END DO
    
    ! thinning process
    DO i = 1, WakePosition_1
        DO j = 1, WakePosition_2+1
            wake_center(i,j,1) = 1
            wake_center(i,j,2) = arr((i-1)*use_nth+1) 
            wake_center(i,j,3) = arr((i-1)*use_nth+1) + REAL(hub_height)
            IF (wake_center(i,j,3) < 5) THEN
                wake_center(i,j,3) = 5
            END IF
            
            IF (wake_center(i,j,2) > maxVal) THEN
                maxVal = wake_center(i,j,2)
            END IF
            IF (wake_center(i,j,2) < minVal) THEN
                minVal = wake_center(i,j,2)
            END IF
        END DO
    END DO
    
    res_mean = SUM(wake_center(1:WakePosition_1,1,2)) / WakePosition_1
    res_std = SQRT (SUM((wake_center(1:WakePosition_1,1,2)-res_mean)**2) / WakePosition_1)

END SUBROUTINE Get_wake_center_RW 
    
!----------------------------------------------------------------------------------
FUNCTION TI_downstream_total (spacing,angle,velocity_matrix)   ! name should be calculate_TI_downstream
!..................................................................................
!  This subroutine is called to calculate the TI of the wake deficit
!  The method is by the paper of Rolf-Erik
!  The output is TI_downstream_matrix which is the TI for each computating node in the DWM domain
!..................................................................................
    USE  DWM_Wake_Deficit_Data,             ONLY: Turb_Stress_DWM, n_x_vector, n_r_vector, ppR, TI_original
    USE  TI_downstream_data
    USE  DWM_Wake_main_Data,              ONLY: wake_position
    USE  meandering_data,                 ONLY: moving_time, scale_factor
    USE  parameter_file_data,                ONLY: TI_amb, hub_height,Uambient,WakePosition_1,Uambient, next_hub_height
    USE  smooth_out_wake_data,            ONLY: length_velocity_array
    USE  DWM_Wake_main_Data,              ONLY: wake_u,wake_position
    
    REAL             ::       TI_downstream_total     ! TI of a downstream turbine
    REAL             ::       spacing                 ! the spacing between the downwind turbine and this upwind turbine
    REAL             ::       angle                   ! the angle between the downwind turbine and the line connecting this upwind turbine and the wind direction
    REAL             ::       velocity_matrix(:,:)    ! the velocity matrix at the certain downswind turbine
    
    REAL             ::       c_uw
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI at each node at the downstream turbine plane from the wake deficit calculation
   !------------------------------------------------------------------------------------------------- 
    IF (ALLOCATED( TI_downstream_matrix ))                DEALLOCATE ( TI_downstream_matrix )
    ALLOCATE (TI_downstream_matrix(n_x_vector,n_r_vector))
    
    c_uw = (0.7550 - TI_original/100 *1.75) / 2
    
    DO i=1,n_x_vector
       DO j=1,n_r_vector
          TI_downstream_matrix(i,j) = abs( ( 1/c_uw *Turb_Stress_DWM(i,j)   ) )**0.5
          !TI_downstream_matrix(i,j) = abs( ( 10/3*Turb_Stress_DWM(i,j)   ) )**0.5
       END DO
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI of the downstream turbine grid considering the meandering effect
   !------------------------------------------------------------------------------------------------- 
    cross_plane_position_TI = ANINT( ppR*spacing+1 )
                     !cross_plane_position_TI = ANINT( ppR*4.4+1 )
    cross_plane_position_ds = ANINT( (ppR/scale_factor)*spacing+1 )  ! the moving time index of the cross plane in the wake_position(:,:,:)
                     !cross_plane_position_ds = ANINT( (ppR/scale_factor)*4.4+1 )
    
    Rscale = 2 !1.3  ! 7.31.2014 test
    !HubHt = hub_height
    HubHt = next_hub_height    ! HubHt is the hub height of the downwind turbine
    counter1 = 0
    counter2 = 0
    TI_accumulation = 0
    TI_apprant_accumulation = 0
    
    DO i=1,WakePosition_1,1
       DO j=ANINT(HubHt-Rscale*R),ANINT(HubHt+Rscale*R),1      ! Z direction
          DO k=1,2*Rscale*R+1,1                         ! Y direction
             y_axis_turbine = k-(Rscale*R+1) + 2*R*spacing*TAN(angle*Pi/180) 
             z_axis_turbine = j
             
             wake_center_y=wake_position(i,cross_plane_position_ds,2) !+ 2*R*spacing*TAN(skew_angle)         ! shift effect
             wake_center_z=wake_position(i,cross_plane_position_ds,3) !- REAL(TurbRefHt-hub_height)
             
             distance = ( (y_axis_turbine-wake_center_y)**2 + (z_axis_turbine-wake_center_z)**2)**0.5
             
             distance_index = FLOOR(distance/(R/ppR)) + 1
              
             TI_node_temp = TI_downstream_matrix( cross_plane_position_TI,distance_index )
             
             IF ( TI_node_temp > (TI/100*(Mean_FFWS/Uambient)) ) THEN
                 TI_node = TI_node_temp
             ELSE
                 TI_node = TI/100*(Mean_FFWS/Uambient)
             END IF
             
             TI_accumulation = TI_accumulation + TI_node
             
             !TI_accumulation = TI_accumulation + TI_node_temp       ! 7.31.2014 test
             
             counter1 = counter1+1
             
          END DO
       END DO
    END DO
    
    TI_average = TI_accumulation / REAL(counter1) 
   
   !-------------------------------------------------------------------------------------------------
   ! calculate the apprant TI due to the meandering
   !------------------------------------------------------------------------------------------------- 
    zero_spacing = 0
    initial_timestep = 1
    !ALLOCATE (velocity_matrix(2*length_velocity_array,2*length_velocity_array))
    
    DO i=1,2*length_velocity_array,1                                                        ! velocity_matrix (i,:)
      DO j=1,2*length_velocity_array,1                                                      ! velocity_matrix (:,j)
          y = (0-length_velocity_array) + (i-1)                                              ! y coordinate
          z = (hub_height-length_velocity_array) + (j-1)                                     ! z coordinate 
          TI_apprant_accumulation = TI_apprant_accumulation + &
                                    !(smooth_wake_shifted_velocity(y,z,wake_u(floor(spacing_turbine * ppR)+1,:), wake_position(:,:,:),spacing_turbine,k) - &
                                    !velocity_matrix(i,j))**2
                                    (velocity_matrix(i,j) - &
                                    smooth_wake_shifted_velocity(y,z,wake_u(floor(spacing * ppR)+1,:), wake_position(:,:,:),zero_spacing,initial_timestep,angle))**2
                                           !smooth_wake_shifted_velocity(y,z,wake_u(floor(4.4 * ppR)+1,:), wake_position(:,:,:),zero_spacing,initial_timestep))**2
          counter2=counter2+1
      END DO
    END DO
    
    TI_apprant = ((TI_apprant_accumulation / REAL(counter2))**0.5)
  
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the total TI
   !------------------------------------------------------------------------------------------------- 
    !!!!TI_downstream_total = (TI_average**2 + TI_apprant**2)**0.5*100
    
    TI_downstream_total = TI_average * 100     !9.30.2014
     
    
    !OPEN(unit = 10, status='replace',file='DWM\results\Downstream_TI_b4normalization.bin',form='unformatted')          
    !WRITE(10)   TI_total                                                                                                                                                                                                                                                                                                         
    !CLOSE(10)
    
    !open (unit=25,file="D:\5MW_simulation\after_release\results\TI.txt")
    !write (25,'(f13.7)'), TI_total
    
    !print*, TI_downstream_matrix(300,120)
    
END FUNCTION TI_downstream_total

!------------------------------------------------------------------------------------------------ 
FUNCTION smallscale_TI (spacing,angle,velocity_matrix)
!............................................................................
! This routine is called to obtain the small scale TI due to the wake deficit
! 
!............................................................................
    USE  DWM_Wake_Deficit_Data,             ONLY: Turb_Stress_DWM, n_x_vector, n_r_vector, ppR, TI_original
    USE  TI_downstream_data
    USE  DWM_Wake_main_Data,              ONLY: wake_position
    USE  meandering_data,                 ONLY: moving_time, scale_factor
    USE  parameter_file_data,                ONLY: TI_amb, hub_height,Uambient,WakePosition_1,Uambient, next_hub_height
    USE  DWM_Wake_main_Data,              ONLY: wake_u,wake_position
    
    REAL             ::       smallscale_TI           ! small scale TI
    REAL             ::       spacing                 ! the spacing between the downwind turbine and this upwind turbine
    REAL             ::       angle                   ! the angle between the downwind turbine and the line connecting this upwind turbine and the wind direction
    REAL             ::       velocity_matrix(:,:)    ! the velocity matrix at the certain downswind turbine
    
    REAL             ::       c_uw
    REAL             ::       curTI
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI at each node at the downstream turbine plane from the wake deficit calculation
   !------------------------------------------------------------------------------------------------- 
    IF (ALLOCATED( TI_downstream_matrix ))                DEALLOCATE ( TI_downstream_matrix )
    ALLOCATE (TI_downstream_matrix(n_x_vector,n_r_vector))
    
    c_uw = (0.7550 - TI_original/100 *1.75) / 2
    
    DO i=1,n_x_vector
       DO j=1,n_r_vector
          TI_downstream_matrix(i,j) = abs( ( 1/c_uw *Turb_Stress_DWM(i,j)   ) )**0.5 
          !TI_downstream_matrix(i,j) = abs( ( 10/3*Turb_Stress_DWM(i,j)   ) )**0.5
       END DO
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI of the downstream turbine grid considering the meandering effect
   !------------------------------------------------------------------------------------------------- 
    cross_plane_position_TI = ANINT( ppR*spacing+1 )
                     !cross_plane_position_TI = ANINT( ppR*4.4+1 )
    cross_plane_position_ds = ANINT( (ppR/scale_factor)*spacing+1 )  ! the moving time index of the cross plane in the wake_position(:,:,:)
                     !cross_plane_position_ds = ANINT( (ppR/scale_factor)*4.4+1 )
    
    Rscale = 2  !1.3   ! 7.31.2014 test
    !HubHt = hub_height
    HubHt = next_hub_height
    counter1 = 0
    counter2 = 0
    TI_accumulation = 0
    TI_apprant_accumulation = 0
    curTI  = TI/100*(Mean_FFWS/Uambient)
    
    DO i=1,WakePosition_1,1
       DO j=ANINT(HubHt-Rscale*R),ANINT(HubHt+Rscale*R),1      ! Z direction
          DO k=1,2*Rscale*R+1,1                         ! Y direction
             y_axis_turbine = k-(Rscale*R+1) + 2*R*spacing*TAN(angle*Pi/180) 
             z_axis_turbine = j
             
             wake_center_y=wake_position(i,cross_plane_position_ds,2) !+ 2*R*spacing*TAN(skew_angle)         ! shift effect
             wake_center_z=wake_position(i,cross_plane_position_ds,3) !- REAL(TurbRefHt-hub_height)
             
             distance = ( (y_axis_turbine-wake_center_y)**2 + (z_axis_turbine-wake_center_z)**2)**0.5
             
             distance_index = FLOOR(distance/(R/ppR)) + 1
             
             TI_node_temp = TI_downstream_matrix( cross_plane_position_TI,distance_index )
             
             
             IF ( TI_node_temp > curTI ) THEN
                 TI_accumulation = TI_accumulation + TI_node_temp
             ELSE
                 TI_accumulation = TI_accumulation + curTI
             END IF
             
             !TI_accumulation = TI_accumulation + MAX(TI_node_temp, TI/100*(Mean_FFWS/Uambient))
             counter1 = counter1+1
             
          END DO
       END DO
    END DO
    
    TI_average = TI_accumulation / REAL(counter1)
    
    smallscale_TI = TI_average * 100
    
END FUNCTION smallscale_TI

!------------------------------------------------------------------------------------------------ 
FUNCTION shifted_velocity( y, z, upwind_mean_u, Uwake, WakeCenter,spacing,angle)
!............................................................................
! This routine is called to get the DWM wake velocity at a certain point in the downstream turbine plane
! Consideirng the meandered wake center
! Uwake(:) is the axial velocity of the wake at the downstream turbine plane
! WakeCenter(:,:,:) is the wake center (y,z) at the downstream turbine plane 
!............................................................................
 
    USE    SimCont,                       ONLY: ZTime
    USE    parameter_file_data,              ONLY: p_p_r, Wind_file_Mean_u,hub_height, ranW, WakePosition_1
    USE    meandering_data,               ONLY: U_factor
    
    USE    parameter_file_data,           ONLY: hub_height
    USE    TurbConf,                      ONLY: TowerHt
    USE    RW_data,                       ONLY: resTemp, curTime

    REAL       ::   y,z                           ! point location on the y,z axis
    REAL       ::   Uwake(:)                      ! axial velocity of the wake at the downstream turbine plane
    REAL       ::   upwind_mean_u                 ! the mean velocity of the turbine UPstream
    REAL       ::   WakeCenter(:,:,:)             ! wake_center
    REAL       ::   spacing                       ! the distance from the downstream turbine to the upstream turbine
    REAL       ::   angle                         ! the angle between the investigated turbine and the line connecting the upwind turbine and wind origin
    
    REAL       ::   shifted_velocity                   ! the output
                                                  !   the velocity at a certain point
    
    REAL       ::   distance                      ! the distance from the point to the meandered wake center
    REAL       ::   y0                            ! wake center position on y axis
    REAL       ::   z0                            ! wake center position on z axis
    REAL       ::   unit                          ! single unit length  R/ppR
    REAL       ::   scale_factor
    INTEGER    ::   p1
    INTEGER    ::   p2
    INTEGER    ::   time_position                 ! to define which plane's wake center is used
    REAL       ::   Yshifted
    REAL       ::   Zshifted
    

    INTEGER    ::   ILo
    
    
    !z = (hub_height - TowerHt) + z                ! shift the z coordinate based on the new hub height
    
    
    !ALLOCATE  (Uwake(NINT( ppR*Rdomain ))) ! the axis symmetrical velocity
    !ALLOCATE  (WakeCenter( size_of_WakeCenter1,size_of_WakeCenter2,3 ))
    
    IF (ranW /= 1) THEN
        scale_factor = 10
    
        time_position = floor(ZTime/( (2*R/p_p_r/upwind_mean_u/1.00)*scale_factor ))+1  ! ZTime/(DWM_time_step*scale_factor)
                !OPEN (unit=25,file='DWM_WIND_FARM\results\time_position.txt',POSITION='APPEND')
                !OPEN (unit=25,file="DWM\results\lateral_wake_center.txt",POSITION='APPEND')
                !WRITE (25,'(I5)'), time_position
    
    
        y0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,2) !+ 2*R*spacing*TAN(skew_angle)
        z0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,3) !- REAL(TurbRefHt-hub_height)
    ELSE
        ILo = 1
        z0 = InterpBin( ZTime, curTime, resTemp, ILo, WakePosition_1)       ! vertical direction
        IF (z0 < 5) THEN
            z0 = 5
        END IF
        y0 = z0-hub_height                                                  ! horizontal direction
    END IF
                      
    
    Yshifted = y + 2*R*spacing*TAN(angle*Pi/180)
    Zshifted = z
    
    distance = ( (Yshifted-y0)**2 + (Zshifted-z0)**2 )**(0.5)
    unit=R/p_p_r
    
    p1=FLOOR(distance/unit)
    p2=p1+1
    IF (p1>0) THEN
       shifted_velocity = Uwake(p1)+( Uwake(p2)-Uwake(p1) )*( (distance/unit)-p1 )    ! Weighting method
    ELSE
       shifted_velocity = Uwake(p2)
    END IF
    
    !shifted_velocity = shifted_velocity * upwind_mean_u
       
    
    !shifted_velocity = shifted_velocity*Wind_file_Mean_u                   !            ! real velocity in m/s
    
    
   ! IF (ALLOCATED( Uwake ))                  DEALLOCATE ( Uwake )
   ! IF (ALLOCATED( WakeCenter ))             DEALLOCATE ( WakeCenter )


END FUNCTION shifted_velocity
!------------------------------------------------------------------------------------------------ 
SUBROUTINE read_parameter_file()
!............................................................................
! This routine is called to read the parameter file from the DWM simulation of upstream turbine
! read wake velocity @ the downstream turbine from the upstream wake
! read the meandered wake center
! read the mean velocity of the upstream turbine
!............................................................................
    USE   parameter_file_data
    USE   TurbConf
    !USE   read_turbine_position_data,   ONLY: SimulationOrder_index 

    OPEN(unit = 10, status='old',file='DWM-driver\DWM_parameter.bin',form='unformatted')  ! open an existing file
    READ(10) hub_height, next_hub_height, NumWT, Uambient, TI_amb, r_domain, x_domain, p_p_r, WakePosition_1, WakePosition_2,WFLowerBd, Winddir,Tinfluencer, ranW
    CLOSE(10) ! close the file
    
    !hub_height = TowerHt
    RotorR     = TipRad
    
    print*,NumWT
    print*,hub_height, next_hub_height
           
    !IF (SimulationOrder_index > 1) THEN
      !OPEN(unit = 10, status='old',file='DWM_WIND_FARM\results\Mean_U_Turbine_1.txt',form='formatted')
      !READ(10,'(f13.7)'), Wind_file_Mean_u
      !CLOSE(10)
      
      !!ALLOCATE (WakePosition(WakePosition_1,WakePosition_2,3))
      !!OPEN(unit = 10, status='old',file='DWM\results\meandered_wake_center.bin',form='unformatted')  
      !!READ(10) WakePosition 
      !!CLOSE(10)
      
      !!ALLOCATE (velocityU(floor(p_p_r*r_domain)))
      !!OPEN(unit = 10, status='old',file='DWM\results\Uvelocity.bin',form='unformatted')  
      !!READ(10) velocityU
      !!CLOSE(10)
      
      !!OPEN(unit = 10, status='old',file='DWM\results\upstream_meanU.bin',form='unformatted')  
      !!READ(10) U_mean 
      !!CLOSE(10)
      
      !!OPEN(unit = 10, status='old',file='DWM\results\Downstream_TI_b4normalization.bin',form='unformatted')          
      !!READ(10) TI_wake                                                                                                                                                                                                                                                                                                         
      !!CLOSE(10)
      
      !!ALLOCATE ( smoothed_wake (NELM) )
      !!OPEN(unit = 10, status='old',file='DWM\results\smoothed_vlocity.bin',form='unformatted')          
      !!READ(10) smoothed_wake                                                                                                                                                                                                                                                                                                         
      !!CLOSE(10)
      

    !END IF
    
END SUBROUTINE read_parameter_file

!------------------------------------------------------------------------------------------------ 
SUBROUTINE read_upwind_result_file()
!............................................................................
! This routine is called to read the results from the DWM simulation of upwind turbines
! and to generate the output variables
!............................................................................
    USE read_turbine_position_data, ONLY:upwindturbine_number,upwind_turbine_index,WT_index,downwindturbine_number,SimulationOrder_index,Turbine_sort_order
    USE parameter_file_data,        ONLY:p_p_r,r_domain,WakePosition_1,WakePosition_2,Wind_file_Mean_u, ranW
    USE read_upwind_result_file_data
    USE Element,                    ONLY:NELM
    USE smooth_out_wake_data,       ONLY:length_velocity_array
    USE RW_data,                    ONLY:resTemp, curTime
    
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    CHARACTER(LEN=3)  :: upwind_turbine_index_character
    CHARACTER(LEN=80) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin,filename_TI_bin,filename_smoothWake_bin,filename_smallTI_bin
    CHARACTER(LEN=80) :: filename_meanU_txt
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'
    CHARACTER(LEN=3)  :: Tiprefix_bin       = 'TI_'
    CHARACTER(LEN=8)  :: smallTIprefix_bin  = 'SmallTI_'
    CHARACTER(LEN=11) :: SmoothWprefix_bin  = 'Smoothwake_'
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results\'
    CHARACTER(LEN=4)  :: connectionprefix   = '_to_'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    INTEGER           :: I
    CHARACTER(LEN=3)  :: turbine_sort_order_char
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! read the wind file mean velocity at the turbine plane from the very first turbine
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    IF (SimulationOrder_index > 1) THEN
        IF (Turbine_sort_order(1) <= 9) THEN
            write(turbine_sort_order_char,'(i1)') Turbine_sort_order(1)
        ELSEIF (Turbine_sort_order(1) <= 99) THEN
            write(turbine_sort_order_char,'(i2)') Turbine_sort_order(1)
        ELSE
            write(turbine_sort_order_char,'(i3)') Turbine_sort_order(1)
        END IF
        
        filename_meanU_txt = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(turbine_sort_order_char)//".txt"
        OPEN(unit = 10, status='old',file=filename_meanU_txt,form='formatted')
        READ(10,'(f13.7)'), Wind_file_Mean_u
        CLOSE(10)
    END IF
    
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! read the upwind results if have any
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    IF (upwindturbine_number > 0) THEN
        ALLOCATE (upwind_U          (upwindturbine_number,floor(p_p_r*r_domain)))                   ! declare the input
        ALLOCATE (upwind_wakecenter (upwindturbine_number,WakePosition_1,WakePosition_2,3))
        ALLOCATE (upwind_meanU      (upwindturbine_number))
        ALLOCATE (upwind_TI         (upwindturbine_number))
        ALLOCATE (upwind_small_TI   (upwindturbine_number))
        ALLOCATE (upwind_smoothWake (upwindturbine_number,NELM))
        ALLOCATE (velocity_aerodyn  (upwindturbine_number))                       ! the temp velocity used by the aerodyn
        
        ! transfer the turbine index from integer to character
        IF (WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') WT_index
        ELSEIF (WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') WT_index
        END IF
        
        DO I = 1,upwindturbine_number
            
           IF (upwind_turbine_index(I) <= 9) THEN
              write(upwind_turbine_index_character,'(i1)') upwind_turbine_index(I)
           ELSEIF (upwind_turbine_index(I) <= 99) THEN
              write(upwind_turbine_index_character,'(i2)') upwind_turbine_index(I)
           ELSE
              write(upwind_turbine_index_character,'(i3)') upwind_turbine_index(I)
           END IF
           
           ! obtain the coresponded profile name
            
           filename_u_bin          = trim(Prefix)//trim(Uprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"         ! the file name needs to be read

           filename_TI_bin         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           filename_smallTI_bin    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           
           filename_smoothWake_bin = trim(Prefix)//trim(SmoothWprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)//".bin"
           
           filename_meanU_bin      = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)//".bin"
           
           ! open the file and read
           OPEN(unit = 10, status='old',file=filename_u_bin,form='unformatted')  
           READ(10) upwind_U(I,:)          
           CLOSE(10)
           
           OPEN(unit = 10, status='old',file=filename_TI_bin,form='unformatted')  
           READ(10) upwind_TI(I) 
           CLOSE(10)
           
           OPEN(unit = 10, status='old',file=filename_smallTI_bin,form='unformatted')  
           READ(10) upwind_small_TI(I) 
           CLOSE(10)
           
           OPEN(unit = 10, status='old',file=filename_smoothWake_bin,form='unformatted')  
           READ(10) upwind_smoothWake(I,:) 
           CLOSE(10)
           
           OPEN(unit = 10, status='old',file=filename_wakecenter_bin,form='unformatted')  
           READ(10) upwind_wakecenter(I,:,:,:) 
           CLOSE(10)
           
           OPEN(unit = 10, status='old',file=filename_meanU_bin,form='unformatted')  
           READ(10) upwind_meanU(I) 
           CLOSE(10)
        
        END DO       
    END IF
   
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! declare the downwind output variables
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    IF (SimulationOrder_index > 0) THEN                ! not for the base 0 turbine
        IF (downwindturbine_number /= 0 ) THEN         ! not for the turbines that don't have a downwind turbine
            length_velocity_array = NINT(1.2*R)
            ALLOCATE ( TI_downstream             (downwindturbine_number                                                ) )
            ALLOCATE ( small_scale_TI_downstream (downwindturbine_number                                                ) )
            ALLOCATE ( smoothed_velocity_array   (downwindturbine_number,NELM                                           ) )
            ALLOCATE ( vel_matrix                (downwindturbine_number,2*length_velocity_array,2*length_velocity_array) ) 
        END IF
    END IF
    
    
    ! If random walk model is used, build the array used in the shifted_velocity subroutine
    IF (upwindturbine_number > 0) THEN
        IF (ranW == 1) THEN
            IF (.NOT. ALLOCATED(curTime)) ALLOCATE (curTime(WakePosition_1))    ! time array 0:1:....
            DO i = 1, WakePosition_1
                curTime(i) = i-1
            END DO
        
            IF (.NOT. ALLOCATED(resTemp))  ALLOCATE (resTemp(WakePosition_1))   ! wake center array corresponds to the time array
            DO i = 1, WakePosition_1
                resTemp(i) = upwind_wakecenter(1, i, WakePosition_2, 3)                      ! vertical direction (z), y = z-hubHt 
            END DO
        END IF
    END IF
    
END SUBROUTINE read_upwind_result_file

!------------------------------------------------------------------------------------------------ 
SUBROUTINE write_result_file()
!............................................................................
! This routine is called to write the results from the DWM simulation of this turbine
!............................................................................
    USE read_turbine_position_data, ONLY: WT_index,downwindturbine_number,downwind_turbine_index,downwind_turbine_projected_distance,SimulationOrder_index
    USE DWM_Wake_main_Data,         ONLY: TI,Mean_FFWS,induction_factor,wake_u,wake_position,wake_width,travel_time, avg_ct, Agenpwr, MeanInd, DWM_pitch
    USE parameter_file_data,        ONLY: p_p_r
    USE read_upwind_result_file_data, ONLY:TI_downstream,smoothed_velocity_array,small_scale_TI_downstream
    
    USE DWM_Wake_Deficit_Data,      ONLY: n_x_vector, n_r_vector
    

    
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    CHARACTER(LEN=3)  :: downwind_turbine_index_character
    CHARACTER(LEN=80) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin,filename_TI_bin,filename_smallTI_bin,filename_smoothWake_bin,filename_wakewidth_bin,filename_wake_bin
    CHARACTER(LEN=80) :: filename_TI_txt,filename_meanU_txt,filename_induction_txt,filename_wake_txt,filename_wakecenter_txt,filename_advection_time_txt,filename_Ct_txt,filename_meanpower_txt,filename_meanInd_txt
    CHARACTER(LEN=80) :: filename_pitch_txt
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'
    CHARACTER(LEN=3)  :: Tiprefix_bin       = 'TI_'
    CHARACTER(LEN=8)  :: smallTIprefix_bin  = 'SmallTI_'
    CHARACTER(LEN=14) :: ATprefix_bin       = 'AdvectionTime_'
    CHARACTER(LEN=11) :: SmoothWprefix_bin  = 'Smoothwake_'
    CHARACTER(LEN=10) :: InductionPrefix    = 'Induction_'
    CHARACTER(LEN=6)  :: Wakeprefix         = 'WakeU_'
    CHARACTER(LEN=11) :: WWprefix_bin       = 'Wake_width_'
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results\'
    CHARACTER(LEN=4)  :: connectionprefix   = '_to_'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    CHARACTER(LEN=7)  :: MeanCtPrefix       = 'MeanCt_'
    CHARACTER(LEN=10) :: Powerprefix        = 'Meanpower_'
    CHARACTER(LEN=12) :: MeanAprefix        = 'AveragedInd_'
    CHARACTER(LEN=11) :: PitchPrefix        = 'PitchAngle_'
    INTEGER           :: I,RESULT
    
    CHARACTER(LEN=80) :: filename_TI_to_txt
    CHARACTER(LEN=80) :: filename_SmallTI_to_txt
    
    CHARACTER(LEN=10) :: fmt_1
    CHARACTER(LEN=20) :: fmt
    
    IF ( SimulationOrder_index > 0 ) THEN            ! exclude the first turbine
        
        IF (WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') WT_index
        ELSEIF (WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') WT_index
        END IF
    
        ! Write the TI of this turbine
        filename_TI_txt    = trim(Prefix)//trim(Tiprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
    
          OPEN (unit=25,file=filename_TI_txt)
          WRITE (25,'(f13.7)'), TI 
          CLOSE (25)
    
        ! Write the mean velocity of this turbine
        filename_meanU_bin = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        filename_meanU_txt = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
    
          OPEN (unit=25,file=filename_meanU_txt)
          WRITE (25,'(f13.7)'), Mean_FFWS
          CLOSE (25)
    
          OPEN(unit = 10, status='replace',file=filename_meanU_bin,form='unformatted')    
          WRITE(10)   Mean_FFWS                                                                                                                                                                                                                                 
          CLOSE(10)
    
        ! Write the induction factor of this turbine
        filename_induction_txt = trim(Prefix)//trim(InductionPrefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
    
          OPEN (unit=25,file=filename_induction_txt)
          WRITE (25,'(f13.7)'), induction_factor(:)
          CLOSE (25)
          
        ! Write the pitch angle of this turbine
        filename_pitch_txt = trim(Prefix)//trim(PitchPrefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
          OPEN (unit=25,file=filename_pitch_txt,POSITION = 'APPEND')
          WRITE (25,'(f13.7)'), DWM_pitch
          CLOSE (25)
          
        ! Write the mean Ct of this turbine
        filename_Ct_txt = trim(Prefix)//trim(MeanCtPrefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
          
          OPEN (unit=25,file=filename_Ct_txt)
          WRITE (25,'(f13.7)'), avg_ct
          CLOSE (25)
          
        ! Write the averaged SD power of this turbine
        filename_meanpower_txt = trim(Prefix)//trim(Powerprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
        
          OPEN (unit=25,file=filename_meanpower_txt,POSITION = 'APPEND')
          WRITE (25,'(f13.7)'), Agenpwr
          CLOSE (25)
          
        ! Write the averaged induction factor of this turbine
        filename_meanInd_txt = trim(Prefix)//trim(MeanAprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
        
          OPEN (unit=25,file=filename_meanInd_txt,POSITION = 'APPEND')
          WRITE (25,'(f13.7)'), MeanInd
          CLOSE (25)
        
    
        ! Write the wake deficit profile of this turbine
        filename_wake_txt = trim(Prefix)//trim(Wakeprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
        filename_wake_bin = trim(Prefix)//trim(Wakeprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        
        OPEN (unit=25,file=filename_wake_txt)
          !---- format U ------- 
          !WRITE(fmt_1,'(i4)' ) n_x_vector       
          !fmt = '('// trim(fmt_1) // "(f12.7))"
          
          !DO I = 1,n_r_vector
           !   WRITE(25,fmt),wake_u(1:n_x_vector,I)
          !END DO
          !----------
          
          WRITE (25,'(f9.7)'), wake_u(:,:)
          CLOSE (25)
          
          OPEN(unit = 10, status='replace',file=filename_wake_bin,form='unformatted')    
          WRITE(10)   wake_u(:,:)                                                                                                                                                                                              
          CLOSE(10)
    
        ! Write the meandered wake center result of this turbine
        filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        filename_wakecenter_txt = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
    
          OPEN (unit=25,file=filename_wakecenter_txt)
          WRITE (25,'(f13.7)'), wake_position(:,:,:)
          CLOSE (25)
    
          OPEN(unit = 10, status='replace',file=filename_wakecenter_bin,form='unformatted')    
          WRITE(10)   wake_position(:,:,:)                                                                                                                                                                                              
          CLOSE(10)
    
        ! Write the downstream turbine customized output files
        IF (downwindturbine_number /= 0) THEN
        
            DO I = 1,downwindturbine_number
               IF (downwind_turbine_index(I) <= 9) THEN
                  write(downwind_turbine_index_character,'(i1)') downwind_turbine_index(I)
               ELSEIF (downwind_turbine_index(I) <= 99) THEN
                  write(downwind_turbine_index_character,'(i2)') downwind_turbine_index(I)
               ELSE
                  write(downwind_turbine_index_character,'(i3)') downwind_turbine_index(I)
               END IF

               filename_u_bin          = trim(Prefix)//trim(Uprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               filename_TI_bin         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               filename_smallTI_bin    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"             
               filename_smoothWake_bin = trim(Prefix)//trim(SmoothWprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               
               filename_TI_to_txt         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".txt"
               filename_SmallTI_to_txt    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".txt"
            
           
               ! Write the wake velocity at the certain downstream turbine plane
               
                ! transfer the U to real velocity, not normalized
                
                OPEN(unit = 10, status='replace',file=filename_u_bin,form='unformatted')          
                WRITE(10)   wake_u(floor(downwind_turbine_projected_distance(I) * p_p_r)+1,:)               
                CLOSE(10)
            
               ! Write the TI at the certain downstream turbine plane
                OPEN(unit = 10, status='replace',file=filename_TI_bin,form='unformatted')          
                WRITE(10)   TI_downstream (I)                                                                                                                                                                                                                                                                                                        
                CLOSE(10)
                
                OPEN(unit = 10, status='replace',file=filename_smallTI_bin,form='unformatted')          
                WRITE(10)   small_scale_TI_downstream (I)                                                                                                                                                                                                                                                                                                        
                CLOSE(10)
                
                OPEN (unit=25,file=filename_TI_to_txt)
                WRITE (25,'(f14.7)'), TI_downstream (I)
                
                OPEN (unit=25,file=filename_smallTI_to_txt)
                WRITE (25,'(f14.7)'), small_scale_TI_downstream (I)
            
               ! Write the smoothed wake profile at the certain downstream turbine plane
                OPEN(unit = 10, status='replace',file=filename_smoothWake_bin,form='unformatted')          
                WRITE(10)   smoothed_velocity_array(I,:)          ! write the wind data of the plane where the downstream turbine locates,                                                                                                                                                                                                                                                                
                CLOSE(10)
            END DO
        END IF
    END IF
    
        ! Write the meandered wake center and the wake width result from the 0 base wind turbine
    IF (SimulationOrder_index == 0) THEN
        
            filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//"0"//".bin"
            filename_wakewidth_bin  = trim(Prefix)//trim(WWprefix_bin)//trim(Turbineprefix)//"0"//".bin"
            
            OPEN(unit = 10, status='replace',file=filename_wakecenter_bin,form='unformatted')    
            WRITE(10)   wake_position(:,:,:)                                                                                                                                                                                              
            CLOSE(10)
            
            OPEN(unit = 10, status='replace',file=filename_wakewidth_bin,form='unformatted')    
            WRITE(10)   wake_width(:)                                                                                                                                                                                              
            CLOSE(10)
            
    END IF

    ! Write the wake advection time
    IF (SimulationOrder_index /= 0) THEN
        IF (downwindturbine_number > 0) THEN
            
            DO I = 1,downwindturbine_number
               IF (downwind_turbine_index(I) <= 9) THEN
                  write(downwind_turbine_index_character,'(i1)') downwind_turbine_index(I)
               ELSEIF (downwind_turbine_index(I) <= 99) THEN
                  write(downwind_turbine_index_character,'(i2)') downwind_turbine_index(I)
               ELSE
                  write(downwind_turbine_index_character,'(i3)') downwind_turbine_index(I)
               END IF
            
               filename_advection_time_txt = trim(Prefix)//trim(ATprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                            //trim(connectionprefix)//trim(downwind_turbine_index_character)//".txt"
            
               OPEN (unit=25,file=filename_advection_time_txt)
               WRITE (25,*), travel_time(1)
               WRITE (25,*), travel_time(2)
               WRITE (25,*), travel_time(3)
               WRITE (25,*), travel_time(4)
            END DO
        END IF
    END IF
       
END SUBROUTINE write_result_file

!----------------------------------------------------------------------------------
SUBROUTINE rename_FAST_output()
!............................................................................
! This routine is called to rename the fast output
!............................................................................
    USE IFPORT
    USE DFLIB
    USE read_turbine_position_data, ONLY: WT_index, SimulationOrder_index
    
    CHARACTER(LEN=80) :: filename_FastOutput,filename_FastElm
    CHARACTER(LEN=11) :: Fastprefix         = 'FastOutput_' ! Fast output file
    CHARACTER(LEN=8)  :: FastElmprefix      = 'FastElm_'    ! Fast Elm output file
    INTEGER           :: RESULT
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results\'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    
    IF ( SimulationOrder_index > 0 ) THEN            ! exclude the first turbine
        
        IF (WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') WT_index
        ELSEIF (WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') WT_index
        END IF
    
            ! Rename the FAST output wrt the turbine index
        filename_FastOutput = trim(Prefix)//trim(Fastprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".out"
        filename_FastElm    = trim(Prefix)//trim(FastElmprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".elm"
    
        !RESULT =rename('V80_2MW.out',filename_FastOutput)          
        !RESULT =rename('V80_2MW.elm',filename_FastElm   )
        RESULT =rename('SWT-2.3.out',filename_FastOutput)          
        RESULT =rename('SWT-2.3.elm',filename_FastElm   )
    END IF

END SUBROUTINE rename_FAST_output

!----------------------------------------------------------------------------------
FUNCTION turbsim_mean_velocity()
!............................................................................
! This routine is called to get the mean velocity of the turbsim wind file
!............................................................................
    USE parameter_file_data,  ONLY : Wind_file_Mean_u, Uambient
    INTEGER   ::   ErrStat
    REAL      ::   temp_wind(3)
    REAL      ::   turbsim_mean_velocity
    

    !temp_wind(:)           = WindInf_ADhack_diskVel( 0.0, (/0.0,0.0,90.0/), ErrStat )
    !turbsim_mean_velocity  = temp_wind(1)
    turbsim_mean_velocity = Wind_file_Mean_u
    
END FUNCTION turbsim_mean_velocity

!----------------------------------------------------------------------------------
SUBROUTINE read_turbine_position()
!..................................................................................
! This routine is called at the first of the FAST.
! To decide the position of the turbine in a row
!   if it is the first turbine or not
!..................................................................................
    
    USE read_turbine_position_data
    USE parameter_file_data,           ONLY: NumWT,Tinfluencer
    
    INTEGER   ::    I,J
    
    N_ARGU = IARGC()
    CALL GETARG(1,SimulationOrder_index_char)
    print*,SimulationOrder_index_char
    READ(SimulationOrder_index_char,*) SimulationOrder_index
    print*,SimulationOrder_index
    
    
    IF ( SimulationOrder_index /= 0 ) THEN              ! exclude the base turbine
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the wind turbine index
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 
         
        ALLOCATE (Turbine_sort_order(NumWT))
        OPEN(unit = 10, status='old',file='DWM-results\wind_farm_turbine_sort.bin',form='unformatted')  ! open an existing file
        READ(10)  Turbine_sort_order 
        CLOSE(10) ! close the file
    

        WT_index = Turbine_sort_order(SimulationOrder_index)

        
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         !''''''''''''''''''''''''''''''''''''UPWIND  DIRECTION''''''''''''''''''''''''''''''''''''''
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the index of upwind turbines that affecting this turbine, and the distance/angle
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        IF (SimulationOrder_index > 0) THEN
            ALLOCATE (TurbineInfluenceData(NumWT,Tinfluencer ) )
        
            OPEN(unit = 10, status='old',file='DWM-results\turbine_interaction.bin',form='unformatted')  
            READ(10)  TurbineInfluenceData
            CLOSE(10) 
        
            ! obtain the upwind turbine index
            ALLOCATE (upwind_turbine_index(Tinfluencer))
            DO I = 1,Tinfluencer
                upwind_turbine_index(I) = TurbineInfluenceData(WT_index,I)
            END DO
        
            ! calculate the number of upwind turbines affecting the downwind turbine
            upwindturbine_number = 0
            DO I = 1,Tinfluencer
                IF (upwind_turbine_index(I) /=0) THEN
                upwindturbine_number = upwindturbine_number + 1
                END IF
            END DO
        
            ! obtain the upwind turbine coordinates
            ALLOCATE (wind_farm_Xcoor (NumWT))
            ALLOCATE (wind_farm_Ycoor (NumWT))
        
            OPEN(unit = 10, status='old',file='DWM-driver\wind_farm_coordinate.bin',form='unformatted')  
            READ(10)  wind_farm_Xcoor,wind_farm_Ycoor
            CLOSE(10)
        
        
            IF (upwindturbine_number /= 0) THEN
                ALLOCATE(upwind_turbine_Xcoor(upwindturbine_number))
                ALLOCATE(upwind_turbine_Ycoor(upwindturbine_number))
                DO I = 1,upwindturbine_number
                    upwind_turbine_Xcoor(I) = wind_farm_Xcoor(upwind_turbine_index(I))
                    upwind_turbine_Ycoor(I) = wind_farm_Ycoor(upwind_turbine_index(I))
                END DO
            END IF
        
            ! obtain the distance beween the upwind turbine and this turbine
            ALLOCATE (turbine_windorigin_length (NumWT      ))
        
            OPEN(unit = 10, status='old',file='DWM-results\turbine_distance.bin',form='unformatted')  
            READ(10)  turbine_windorigin_length
            CLOSE(10)
        
        
            IF (upwindturbine_number /= 0) THEN
                ALLOCATE (upwind_turbine_projected_distance(upwindturbine_number))
                DO I = 1,upwindturbine_number
                    upwind_turbine_projected_distance(I) = turbine_windorigin_length(WT_index) - turbine_windorigin_length(upwind_turbine_index(I))
                END DO
            END IF
        
        
            ! obtain the angle beween the line connecting the upwind turbine and this turbine and the wind direction vector
            ALLOCATE (turbine_angle(NumWT,NumWT))
        
            OPEN(unit = 10, status='old',file='DWM-results\turbine_angles.bin',form='unformatted')  
            READ(10)  turbine_angle
            CLOSE(10)
        
            IF (upwindturbine_number /= 0) THEN
                ALLOCATE (upwind_align_angle(upwindturbine_number))
                DO I = 1,upwindturbine_number
                       upwind_align_angle(I) = turbine_angle(upwind_turbine_index(I),WT_index)
                END DO
            END IF 
        
        END IF
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         !''''''''''''''''''''''''''''''''''DOWNWIND  DIRECTION''''''''''''''''''''''''''''''''''''''
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the index of downwind turbines that being affected by this turbine, and the distance/angle
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 
    
        IF (SimulationOrder_index/=NumWT) THEN    ! WHEN index = 0, there will not such files (which are generated by the 0 turbine)
        
           ! obtain the downwind turbine index and calculate the downwind turbine numbers
            downwindturbine_number = 0
            ALLOCATE (downwind_turbine_index(NumWT-1))
            DO I = 1,Tinfluencer
                DO J = 1,NumWT
                    IF (TurbineInfluenceData(J,I) == WT_index) THEN
                        downwindturbine_number = downwindturbine_number + 1
                        downwind_turbine_index(downwindturbine_number) = J
                    END IF
                END DO
            END DO
        
            ! obtain the downwind turbine coordinates
            IF (downwindturbine_number /= 0) THEN
                ALLOCATE(downwind_turbine_Xcoor(downwindturbine_number))
                ALLOCATE(downwind_turbine_Ycoor(downwindturbine_number))
                DO I = 1,downwindturbine_number
                    downwind_turbine_Xcoor(I) = wind_farm_Xcoor(downwind_turbine_index(I))
                    downwind_turbine_Ycoor(I) = wind_farm_Ycoor(downwind_turbine_index(I))
                END DO
            END IF
        
            ! obtain the distance beween the upwind turbine and this turbine
            IF (downwindturbine_number/=0) THEN
                ALLOCATE (downwind_turbine_projected_distance(downwindturbine_number))
                DO I = 1,downwindturbine_number
                   downwind_turbine_projected_distance(I) = turbine_windorigin_length(downwind_turbine_index(I)) - turbine_windorigin_length(WT_index)
                END DO
            END IF
        
            ! obtain the angle beween the line connecting the downwind turbine and this turbine and the wind direction vector
            IF (downwindturbine_number/=0) THEN
                ALLOCATE (downwind_align_angle(downwindturbine_number))
                DO I = 1,downwindturbine_number
                    downwind_align_angle(I) = turbine_angle(WT_index,downwind_turbine_index(I))
                END DO
            END IF
        
        END IF
    END IF   
              
END SUBROUTINE read_turbine_position

!----------------------------------------------------------------------------------
SUBROUTINE Turbulence_kaimal_spectrum(intensity,frequency_range,KS)
!..................................................................................
!  This routine is used to get the kaimal spectrum of a specific TI level turbulence
!  The output is the frequency_range and the KS (spectrum)
!  Which are both one dimension arrays.  
!..................................................................................
    USE Turbulence_KS
    USE parameter_file_data,   ONLY: Uambient
    
    REAL                ::   intensity
    REAL,ALLOCATABLE    ::   frequency_range(:)     
    REAL,ALLOCATABLE    ::   KS(:)
    
    fs = 1000                       ! sample frequency
    low_f = 0                       ! lower bound of frequency range
    high_f = 100                    ! upper bound of frequency range
    temp_n = floor(fs*(high_f-low_f))+1
    
    ALLOCATE (frequency_range(temp_n))
    frequency_range = ((high_f-low_f)/(temp_n-1))*[(i,i=1,temp_n)]+(low_f-((high_f-low_f)/(temp_n)))
    
    STD = intensity*Uambient
    lk_facor = 8.1*0.7*30
    
    DO i=1,temp_n
       KS(i)= ( 4*STD**2*lk_facor/Uambient ) / (1+6*frequency_range(i)*lk_facor/Uambient)**(5/3)
    END DO
    
END SUBROUTINE Turbulence_kaimal_spectrum

!----------------------------------------------------------------------------------
SUBROUTINE shinozuka(frequency_scale,Kspectrum,u_syn)
!..................................................................................
! This routine is called to calculate the turbulence wind speed
! Using the Kaimal spectrum abtained previously
! The output is the wind speed time series u_syn 
!..................................................................................
    USE shinozuka_data
    USE Turbulence_KS,    ONLY:  low_f, high_f, temp_n
    
    REAL,ALLOCATABLE   ::  u_syn(:) 
    REAL,ALLOCATABLE   ::  frequency_scale(:)
    REAL,ALLOCATABLE   ::  Kspectrum(:)
    
    ALLOCATE (frequency_scale(temp_n) )
    ALLOCATE (Kspectrum      (temp_n) )
    
    num_points = 10000      ! the total point number of the time series
    dt = 1/(high_f-low_f)
    t_min = 0
    t_max = (num_points-1)*dt
    df = 1/(dt*num_points)  ! the sample frequency of the frequency in the spectrum
    
    ALLOCATE (f_syn(num_points-1))
    ALLOCATE (t_syn(num_points))   ! every time point must have a velocity
    ALLOCATE (u_syn(num_points))
    ALLOCATE (phi(num_points))
    ALLOCATE (p_k(num_points-1))
    ALLOCATE (a_k(num_points-1))
    
    ! generate frequency series of the spectrum (discretization)
    f_syn(1) = low_f + df/2     
    DO i=2,num_points-1
       f_syn(i) = f_syn(i-1)+df
    END DO
    
    ! generate time series
    ! t_syn=((t_max-t_min)/(num_points-1))*[(i,i=1,num_points)]+(t_min-((t_max-t_min)/(num_points-1)))
    t_syn(1) = t_min
    DO i=2,num_points
       t_syn(i)=t_syn(i-1)+dt
    END DO
    
    ! generate random phase angle series
    CALL RANDOM_SEED()              
    DO i=1,num_points-1
       CALL RANDOM_NUMBER(phi(i))
       phi(i)=phi(i)*2*Pi
    END DO
    
    ! generate random phase angle and parameter Ak
    ILo=1
    DO i=1,num_points-1
       p_k(i) = InterpBin( f_syn(i), frequency_scale, Kspectrum, ILo, temp_n)
       a_k(i) = (p_k(i)*df)**0.5
    END DO
    
    ! generate wind speed time series
    u_syn = 0                ! initial value
    DO i=1,num_points        ! wind speed time series u_syn
       DO j=1,num_points-1   ! the number of the points in the frequency spectrum
       u_syn(i)=u_syn(i)+(2**0.5)*a_k(j)*cos(2*Pi*f_syn(j)*t_syn(i)+phi(j))
       END DO
    END DO
    
END SUBROUTINE shinozuka

!----------------------------------------------------------------------------------
SUBROUTINE smooth_out_wake(Uvelocity,Uwake_center,wake_array,spacing,angle,velocity_matrix)
!..................................................................................
! This routine is called to fillter out the smoothed out upstream wake profile
! The output is the wake_array
! Which is the axisymmetrical wake velocity profile
!..................................................................................
    USE DWM_Wake_Deficit_Data,    ONLY: n_x_vector, n_r_vector, ppR
    USE parameter_file_data,       ONLY: hub_height, WakePosition_1
    USE smooth_out_wake_data
    
    REAL,ALLOCATABLE   ::   Uvelocity(:,:)
    REAL,ALLOCATABLE   ::   Uwake_center(:,:,:)
    REAL               ::   wake_array(:)
    REAL               ::   spacing   ! the spacing between two turbines
    REAL               ::   angle     ! the angle between the downwind turbine and the line conneting the upwind(investigated) turbine and the wind origin
    REAL               ::   velocity_matrix(:,:)      ! the velocity matrix that store the velocity of the downstream turbine plane
    
    !length_velocity_array = NINT(1.2*R)
    
    !ALLOCATE   ( Uvelocity(n_x_vector,n_r_vector) )
    !ALLOCATE   ( Uwake_center(release_time,flying_time+1,3) )
    IF (ALLOCATED( velocity_array ))                DEALLOCATE ( velocity_array )
    IF (ALLOCATED( counter_array ))                 DEALLOCATE ( counter_array )
    ALLOCATE   ( velocity_array (length_velocity_array) )
    !ALLOCATE   ( velocity_matrix (2*length_velocity_array,2*length_velocity_array) )  ! twice the size of the velocity array
    ALLOCATE   ( counter_array (length_velocity_array) )
    !ALLOCATE   ( wake_array (NELM) )
    
    velocity_array=0
    velocity_matrix = 0
    counter=0
    counter_array=0
    
   !-------------------------------------------------------------------------------------------------
   ! get the time averaged velocity matrix
   !-------------------------------------------------------------------------------------------------
    DO i=1,WakePosition_1,1
       DO j=1,2*length_velocity_array,1                        ! y axis
          DO k=1,2*length_velocity_array,1                     ! z axis
             y = (0-length_velocity_array) + (j-1)             ! y coordinate
             z = (hub_height-length_velocity_array) + (k-1)    ! z coordinate
             velocity_matrix(j,k) = velocity_matrix(j,k) + smooth_wake_shifted_velocity( y, z, Uvelocity(floor(spacing * ppR)+1,:), Uwake_center(:,:,:),spacing,i,angle)
                       !velocity_matrix(j,k) = velocity_matrix(j,k) + smooth_wake_shifted_velocity( y, z, Uvelocity(floor(4.4 * ppR)+1,:), Uwake_center(:,:,:),4.4,i)
          END DO
       END DO
       counter = counter+1
    END DO
    
    velocity_matrix = velocity_matrix / counter                

   !-------------------------------------------------------------------------------------------------
   ! get the time averaged axisymmetrical velocity array
   !-------------------------------------------------------------------------------------------------
    DO i=1,length_velocity_array,1                          ! velocity array
       DO j=1,2*length_velocity_array,1                       ! velocity_matrix (j,:)
          DO k=1,2*length_velocity_array,1                    ! velocity_matrix (:,k)
             y = (0-length_velocity_array) + (j-1)            ! y coordinate
             z = (hub_height-length_velocity_array) + (k-1)   ! z coordinate
             IF ( ((y-0)**2+(z-hub_height)**2)**0.5>(i-1) .and. ((y-0)**2+(z-hub_height)**2)**0.5<=i) THEN        
                velocity_array(i) = velocity_array(i)+velocity_matrix(j,k)
                counter_array(i) = counter_array(i)+1
             END IF
          END DO
       END DO
    END DO
    
    DO i=1,length_velocity_array,1                           
       velocity_array(i) = velocity_array(i) / counter_array(i)
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! get the wake array at the RELM node point
   !-------------------------------------------------------------------------------------------------
    
    DO i=1,NELM,1
       low  = FLOOR(RELM(i))
       high = low+1
       wake_array(i) = velocity_array(low) + ( velocity_array(high)-velocity_array(low) )*(RELM(i)-low)  ! Weighting method
    END DO
    
    !OPEN(unit = 10, status='replace',file='DWM\results\smoothed_vlocity.bin',form='unformatted')          
    !WRITE(10)   wake_array(:)                                          ! write the wind data of the plane where the downstream turbine locates,                                                                                                                                                                                                                                                                
    !CLOSE(10)
    
    !open (unit=25,file="D:\5MW_simulation\after_release\results\Uarray.txt")
    !write (25,'(f13.7)'), velocity_array(:)
    
    !open (unit=25,file="D:\5MW_simulation\after_release\results\wakearray.txt")
    !write (25,'(f13.7)'), wake_array(:)
                
END SUBROUTINE smooth_out_wake

!------------------------------------------------------------------------------------------------ 
FUNCTION smooth_wake_shifted_velocity( y, z, Uwake, WakeCenter,spacing,time_position,angle)
!............................................................................
! This routine is called to get the DWM wake velocity at a certain point in the downstream turbine plane
! Consideirng the meandered wake center
! Uwake(:) is the axial velocity of the wake at the downstream turbine plane
! WakeCenter(:,:,:) is the wake center 
! Used to calculate the smoothed out wake profile
! (y,z) at the downstream turbine plane 
!............................................................................
 
    USE    parameter_file_data,              ONLY: p_p_r, hub_height
    USE    smooth_wake_shifted_velocity_data
    USE    meandering_data,               ONLY: scale_factor

    REAL       ::   y,z                           ! point location on the y,z axis
    REAL       ::   Uwake(:)                      ! axial velocity of the wake at the downstream turbine plane
    REAL       ::   WakeCenter(:,:,:)             ! wake_center
    REAL       ::   spacing                       ! the distance from the downstream turbine to the upstream turbine
    INTEGER    ::   time_position                 ! to define which plane's wake center is used
    REAL       ::   angle                         ! the angle between the downwind turbine and the line conneting the upwind(investigated) turbine and the wind origin
    
    REAL       ::   smooth_wake_shifted_velocity  ! the output
                                                  !   the velocity at a certain point
    REAL       ::   Yshifted
    REAL       ::   Zshifted
    
    y0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,2) !+ 2*R*spacing*TAN(skew_angle) 
    z0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,3) !- REAL(TurbRefHt-hub_height)
    
    Yshifted = y + 2*R*spacing*TAN(angle*Pi/180) 
    Zshifted = z
    
    distance = ( (Yshifted-y0)**2 + (Zshifted-z0)**2 )**(0.5)
    unit=R/p_p_r
    
    p1=FLOOR(distance/unit)
    p2=p1+1
    IF (p1>0) THEN
       smooth_wake_shifted_velocity = Uwake(p1)+( Uwake(p2)-Uwake(p1) )*( (distance/unit)-p1 )   ! Weighting method
    ELSE
       smooth_wake_shifted_velocity = Uwake(p2)
    END IF
    
    
END FUNCTION smooth_wake_shifted_velocity

!------------------------------------------------------------------------------------------------ 
SUBROUTINE advection_time(Uvelocity,WakeWidth,spacing,AdvectionT_array)
!............................................................................
! This routine is called to ontain the wake advection time from this turbine
! to the first downwind turbine ---- AdvectionT
!............................................................................
    USE   parameter_file_data,              ONLY: p_p_r,RotorR,Wind_file_Mean_u
    USE   DWM_Wake_main_Data,               ONLY: Mean_FFWS
    USE   read_turbine_position_data,       ONLY: SimulationOrder_index
    
    REAL,ALLOCATABLE      ::    Uvelocity(:,:)
    INTEGER, ALLOCATABLE  ::    WakeWidth(:)
    REAL,ALLOCATABLE      ::    AdvectionT_array(:)
    REAL                  ::    spacing   

  
    INTEGER            ::    D_WTG
    INTEGER            ::    R_WTG
    INTEGER            ::    I,J
    INTEGER            ::    downwind_turbine_location
    INTEGER            ::    critical_r_length
    REAL               ::    delaxi
    REAL               ::    delrad
    REAL               ::    critical_r
    REAL               ::    local_mean_axial_u
    REAL               ::    local_r_length
    REAL               ::    local_weighting_velocity
    REAL               ::    velocity_sum
    REAL               ::    windU
    
    INTEGER            ::    crossplane1,half_way,crossplane3
    
    IF (ALLOCATED( AdvectionT_array ))      DEALLOCATE ( AdvectionT_array )
    ALLOCATE (AdvectionT_array(4))
    AdvectionT_array = 0
      
    D_WTG    =  2
    R_WTG    =  1
    delaxi   =  D_WTG/(p_p_r*1.0) * 2 * RotorR   ! dx in meters
    delrad   =  R_WTG/(p_p_r*1.0) * RotorR       ! dr in meters
    
    downwind_turbine_location = NINT(spacing * 2 * R / delaxi)
    crossplane1 = NINT(spacing * 2 * R / delaxi / 4)      ! 1/4
    half_way    = NINT(spacing * 2 * R / delaxi / 4*2)    ! 2/4
    crossplane3 = NINT(spacing * 2 * R / delaxi / 4*3)    ! 3/4
    
   
    critical_r        = 2*R
    critical_r_length = NINT(critical_r/delrad)
    
    IF (SimulationOrder_index > 1) THEN
        windU = Wind_file_Mean_u      ! the spacial and time averaged wind velocity for the first turbine
    ELSE
        windU = Mean_FFWS
    END IF
    
    DO J = 1,downwind_turbine_location
        local_r_length    = 0
        velocity_sum      = 0
        
        DO I = 1,WakeWidth(J)
           local_weighting_velocity = Uvelocity(J,I) * Pi *( (local_r_length+delrad)**2 - local_r_length**2 )
           velocity_sum   = velocity_sum + local_weighting_velocity
           local_r_length = local_r_length + delrad
        END DO
        
        local_mean_axial_u = velocity_sum / ( Pi * local_r_length**2 )
        
        AdvectionT_array(4) = AdvectionT_array(4) + delaxi/ (local_mean_axial_u*windU)
    END DO
    
    
    DO J = 1,half_way
        local_r_length    = 0
        velocity_sum      = 0
        
        DO I = 1,WakeWidth(J)
           local_weighting_velocity = Uvelocity(J,I) * Pi *( (local_r_length+delrad)**2 - local_r_length**2 )
           velocity_sum   = velocity_sum + local_weighting_velocity
           local_r_length = local_r_length + delrad
        END DO
        
        local_mean_axial_u = velocity_sum / ( Pi * local_r_length**2 )
        
        AdvectionT_array(2) = AdvectionT_array(2) + delaxi/ (local_mean_axial_u*windU)
    END DO
    
    DO J = 1,crossplane1
        local_r_length    = 0
        velocity_sum      = 0
        
        DO I = 1,WakeWidth(J)
           local_weighting_velocity = Uvelocity(J,I) * Pi *( (local_r_length+delrad)**2 - local_r_length**2 )
           velocity_sum   = velocity_sum + local_weighting_velocity
           local_r_length = local_r_length + delrad
        END DO
        
        local_mean_axial_u = velocity_sum / ( Pi * local_r_length**2 )
        
        AdvectionT_array(1) = AdvectionT_array(1) + delaxi/ (local_mean_axial_u*windU)
    END DO
    
    DO J = 1,crossplane3
        local_r_length    = 0
        velocity_sum      = 0
        
        DO I = 1,WakeWidth(J)
           local_weighting_velocity = Uvelocity(J,I) * Pi *( (local_r_length+delrad)**2 - local_r_length**2 )
           velocity_sum   = velocity_sum + local_weighting_velocity
           local_r_length = local_r_length + delrad
        END DO
        
        local_mean_axial_u = velocity_sum / ( Pi * local_r_length**2 )
        
        AdvectionT_array(3) = AdvectionT_array(3) + delaxi/ (local_mean_axial_u*windU)
    END DO
        
END SUBROUTINE advection_time


!!------------------------------------------------------------------------------------------------ 
!SUBROUTINE write_wake_center(local_time,Uwake,mean_u, WakeCenter,spacing)
!............................................................................
! This routine is called to write the wake center at local time
! 
!............................................................................
 !   USE    parameter_file_data,              ONLY: p_p_r, Wind_file_Mean_u,hub_height,TurbRefHt,spacing_turbine,AlignAngle
    
 !   REAL       ::   local_time
 !   REAL       ::   Uwake(:)                      ! axial velocity of the wake at the downstream turbine plane
 !   REAL       ::   mean_u                        ! the mean velocity of the turbine UPstream
 !   REAL       ::   WakeCenter(:,:,:)             ! wake_center
 !   REAL       ::   spacing                       ! the distance from the downstream turbine to the upstream turbine
 !   
 !   REAL       ::   y0                            ! wake center position on y axis
 !   REAL       ::   z0                            ! wake center position on z axis
 !   REAL       ::   scale_factor
 !   INTEGER    ::   time_position                 ! to define which plane's wake center is used
    
    
    !ALLOCATE  (Uwake(NINT( ppR*Rdomain ))) ! the axis symmetrical velocity
    !ALLOCATE  (WakeCenter( size_of_WakeCenter1,size_of_WakeCenter2,3 ))
    
  !  scale_factor = 10
    
  !  time_position = floor(local_time/( (2*R/p_p_r/mean_u)*scale_factor ))+1  ! ZTime/(DWM_time_step*scale_factor)
    
    
    
  !  y0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,2) - 2*R*spacing_turbine*TAN(AlignAngle*Pi/180)
  !  z0 = WakeCenter(time_position,FLOOR(spacing*p_p_r/scale_factor)+1,3) - REAL(TurbRefHt-hub_height)
    
  !  ! write the meandered wake center as a time series
  !            OPEN (unit=25,file="DWM\results\lateral_wake_center.txt",POSITION='APPEND')
  !            WRITE (25,'(f13.7)'), y0
              
  !            OPEN (unit=25,file="DWM\results\vertical_wake_center.txt",POSITION='APPEND')
  !            WRITE (25,'(f13.7)'), z0
              
!END SUBROUTINE write_wake_center

!------------------------------------------------------------------------------------------------ 
FUNCTION min_of_array(ary, ary_length)
!............................................................................
! This routine is called to return the minmum value in an array
!............................................................................
    
    INTEGER ::    ary_length
    REAL    ::    ary(ary_length)
    REAL    ::    min_of_array
    INTEGER ::    I
    
    min_of_array = 10000.00
    
    DO I = 1,ary_length
        IF (ary(I) < min_of_array) THEN
            min_of_array = ary(I)
        END IF
    END DO

END FUNCTION min_of_array

!------------------------------------------------------------------------------------------------ 
FUNCTION max_of_array(ary, ary_length)
!............................................................................
! This routine is called to return the maximum value in an array
!............................................................................
    
    INTEGER ::    ary_length
    REAL    ::    ary(ary_length)
    REAL    ::    max_of_array
    INTEGER ::    I
    
    max_of_array = -10000.00
    
    DO I = 1,ary_length
        IF (ary(I) > max_of_array) THEN
            max_of_array = ary(I)
        END IF
    END DO

END FUNCTION max_of_array

!------------------------------------------------------------------------------------------------ 
FUNCTION max_of_TwoNum(Num1, Num2)
!............................................................................
! This routine is called to return the maximum value of two numbers
!............................................................................

    REAL   ::    Num1
    REAL   ::    Num2
    REAL   ::    max_of_TwoNum
    
    IF (Num1 > Num2 .OR. Num1 == Num2) THEN
        max_of_TwoNum = Num1
    ELSE
        max_of_TwoNum = Num2
    END IF
    
END FUNCTION max_of_TwoNum

!------------------------------------------------------------------------------------------------ 
FUNCTION rotation_lateral_offset(x_spacing)
!............................................................................
! This routine is called to return the wake lateral offset due to the turbine rotation
! (assume the rotor spins clockwise) --> always shift to right (negative)
!............................................................................

    REAL      ::       x_spacing
    REAL      ::       rotation_lateral_offset
    
    ! parameters
    REAL      ::       ad = -4.5
    REAL      ::       bd = -0.01
    
    rotation_lateral_offset = ad + bd*x_spacing
 
END FUNCTION rotation_lateral_offset
    
!------------------------------------------------------------------------------------------------ 
FUNCTION local_skew_angle(yaw_angle, tilde_ct, x_spacing, wake_width, ppr)
!............................................................................
! This routine is called to return the local skew angle at a certain downstream location
!............................................................................
    
    REAL     ::     yaw_angle
    REAL     ::     tilde_ct
    REAL     ::     x_spacing
    INTEGER  ::     wake_width
    REAL     ::     ppr
    REAL     ::     local_skew_angle
    
    IF ( ABS(yaw_angle) > 0.000001 ) THEN
        local_skew_angle = (ppr/wake_width)**2 *COS(yaw_angle)**2 *SIN(yaw_angle) *tilde_ct/2
        local_skew_angle = -local_skew_angle        ! the direction (positive or negative) is opposite to the turbine yaw angle
    ELSE
        local_skew_angle = 0.0
    END IF
    
    local_skew_angle = TAN(local_skew_angle)
 
END FUNCTION local_skew_angle

!------------------------------------------------------------------------------------------------ 
FUNCTION local_ptfm_pitch_angle(yaw_angle, tilde_ct, x_spacing, wake_width, ppr)
!............................................................................
! This routine is called to return the local skew angle at a certain downstream location
!............................................................................
    
    REAL     ::     yaw_angle
    REAL     ::     tilde_ct
    REAL     ::     x_spacing
    INTEGER  ::     wake_width
    REAL     ::     ppr
    REAL     ::     local_ptfm_pitch_angle
    
    IF ( ABS(yaw_angle) > 0.000001 ) THEN
        local_ptfm_pitch_angle = (ppr/wake_width)**2 *COS(yaw_angle)**2 *SIN(yaw_angle) *tilde_ct/2
    ELSE
        local_ptfm_pitch_angle = 0.0
    END IF
    
    local_ptfm_pitch_angle = TAN(local_ptfm_pitch_angle)
 
END FUNCTION local_ptfm_pitch_angle

!------------------------------------------------------------------------------------------------ 
SUBROUTINE calculate_mean_genpwr (total_power, total_timestep)
!------------------------------------------------------------------------------------------------
! This routine is used to calculate the average time step power from SD subroutine
!------------------------------------------------------------------------------------------------
    USE DWM_Wake_main_Data, ONLY : Agenpwr
    
    REAL       ::     total_power
    INTEGER    ::     total_timestep
    
    Agenpwr = total_power / total_timestep
    
END SUBROUTINE  calculate_mean_genpwr

!------------------------------------------------------------------------------------------------ 
SUBROUTINE DWM_init()
!------------------------------------------------------------------------------------------------
! 
!------------------------------------------------------------------------------------------------
    USE DWM_Wake_main_Data, ONLY: DFN_DWM

    ALLOCATE (DFN_DWM(NELM,NB))

END SUBROUTINE DWM_init

!------------------------------------------------------------------------------------------------ 
SUBROUTINE Hub_height_mean_velocity()
!------------------------------------------------------------------------------------------------
!  This subroutine is used to return the velocity at the hub height
!------------------------------------------------------------------------------------------------

    USE DWM_Wake_main_Data,    ONLY: U_hubht
    USE parameter_file_data,   ONLY: hub_height
    USE    TurbConf,                      ONLY: TowerHt
    
    INTEGER      ::  ErrStat
    REAL         ::  Velocity_hub(3)
    
    Velocity_hub = WindInf_GetMean(0.0, 500.0, 0.01, (/0.0,0.0,REAL(TowerHt,ReKi)/),  ErrStat )
    
    U_hubht = Velocity_hub(1)
    
END SUBROUTINE Hub_height_mean_velocity

!------------------------------------------------------------------------------------------------ 
SUBROUTINE collect_velocity(axial_u)
!------------------------------------------------------------------------------------------------
!  This subroutine is used to collect the axial velocity
!------------------------------------------------------------------------------------------------
    
    USE filter_velocity_data,  ONLY : total_velocity, total_velocity_counter
    
    REAL   ::   axial_u
    
    total_velocity = total_velocity + axial_u
    total_velocity_counter = total_velocity_counter + 1

END SUBROUTINE collect_velocity

!------------------------------------------------------------------------------------------------ 
SUBROUTINE Get_ptfmPitch()
!------------------------------------------------------------------------------------------------
!  This subroutine is used to collect the axial velocity
!------------------------------------------------------------------------------------------------
    
    USE Ptfm_pitch_data,    ONLY: colNum, ptfm_pitch
    
    Integer :: use_ptfm_pitch
    Integer :: I
    REAL    :: curPitch
    
    use_ptfm_pitch = 0
    
    OPEN (unit=99, file='C:\Users\Administrator\Dropbox\Research\DATA\OC3_pitch\ptfmPitch_8.txt', status='old', action='read')
    READ (99, *) colNum
    
    ALLOCATE (ptfm_pitch(colNum,1))
    
    if (use_ptfm_pitch == 1) THEN
        DO I = 1, colNum
            READ(99,*) curPitch
            ptfm_pitch(I,1) = curPitch * Pi / 180 
        END DO
    else
        DO I = 1, colNum
            ptfm_pitch(I,1) = 0
        END DO
    end IF
    
    CLOSE(99)

END SUBROUTINE Get_ptfmPitch


END MODULE DWM_Wake_Sub
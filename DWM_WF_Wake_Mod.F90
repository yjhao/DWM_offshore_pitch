! DWM_Wake Modules
!=======================================================================
MODULE DWM_Wake_main_Data       
!...............................................................................
!       public variables
!...............................................................................

    REAL, ALLOCATABLE        ::   turbine_thrust_force (:)
    REAL, ALLOCATABLE        ::   induction_factor     (:)
    REAL, ALLOCATABLE        ::   r_initial            (:)     ! scaled rotor radius 
    REAL, ALLOCATABLE        ::   U_initial            (:)     ! scaled velocity at the rotor
    REAL, ALLOCATABLE        ::   Mean_FFWS_array      (:)     ! Mean velocity of each section on the blade
    REAL                     ::   Mean_FFWS                    ! Mean (total) wind speed at the hub height
    REAL                     ::   TI                           ! the turbulence intensity of the turbine
    REAL                     ::   TI_downstream                ! the TI of a downstream turbine before normalization

!-----------------------------------------
    ! Output Data
!-----------------------------------------

    INTEGER,ALLOCATABLE      ::    wake_width          (:)     ! wake width
    REAL,ALLOCATABLE         ::    wake_u              (:,:)   ! wake velocity
    REAL,ALLOCATABLE         ::    wake_position       (:,:,:) ! wake center position
    REAL,ALLOCATABLE         ::    smoothed_velocity_array (:) ! smoothed out upstream axisymetric wake profile
    !REAL                     ::    DSTU                        ! downstream turbine wind speed
    
    REAL    ::      AtmUscale    ! atmospheric velocity scale before introducing TI
    REAL    ::      du_dz_ABL    ! atmosperic shear gradient
    
    !REAL    ::      skew_angle   ! the difference between the skew angle and the yaw angle
    REAL    ::      ct_tilde
    REAL    ::      avg_ct        ! average thrust coefficient
    
    REAL,ALLOCATABLE    ::      travel_time(:)   ! wake advection time
    INTEGER             ::      AP_timestep = 0  ! time step of average power
    REAL                ::      Totalgenpwr = 0  ! total generator power
    REAL                ::      Agenpwr          ! average generator power
    REAL                ::      MeanInd          ! averaged induction factor
    
    REAL                ::      DWM_pitch
    REAL,ALLOCATABLE    ::      DFN_DWM(:,:)
    
    REAL                ::      U_hubht



END MODULE DWM_Wake_main_Data


!===============================================================================
MODULE turbine_average_velocity_data
!...............................................................................
!       
!...............................................................................
    REAL,ALLOCATABLE    ::   average_velocity_array_temp(:)   ! the average velocity of the whole blade sections 
                                                              !   in a specific time step
    REAL,ALLOCATABLE    ::   average_velocity_array (:)
    
    REAL,ALLOCATABLE    ::   swept_area (:)

END MODULE turbine_average_velocity_data

!=======================================================================
MODULE DWM_Wake_Deficit_Data          
!...............................................................................
! The DWM wake deficit subroutine data   
!...............................................................................

REAL, DIMENSION(2)     ::   filter1
REAL, DIMENSION(2)     ::   filter2
REAL, ALLOCATABLE       ::   F1_vector (:)
REAL, ALLOCATABLE       ::   F2_vector (:)
REAL                    ::   ppR                           ! Point_per_R_resoulution
REAL                    ::   Domain_R                      ! Domain_size_in_radial_direction
REAL                    ::   Domain_X                      ! Domain_size_in_flow_direction
REAL                    ::   TI_original                ! Turbulence_intensity normalized back to ambient wind speed
REAL                    ::   k1                            ! Amb turb. coeff.
REAL                    ::   k2                            ! Shear layer coeff.



INTEGER                 ::   n_x_vector
INTEGER                 ::   n_r_vector 

!%%%%% Rolf modification
INTEGER                 ::   length_F1_vector
REAL                    ::   L_ABL_vector(3)
REAL                    ::   UW_UU_vector(3)
REAL                    ::   L_DEF_vector(3)
REAL                    ::   UU_DEF_UU_ABL_vector(3)
REAL                    ::   UW_DEF_UU_DEF_vector(3)
REAL                    ::   x_ary(3)
REAL                    ::   L_ABL
REAL                    ::   UW_UU
REAL                    ::   L_DEF
REAL                    ::   UU_DEF_UU_ABL
REAL                    ::   UW_DEF_UU_DEF
REAL                    ::   Rotor_fixed_R
REAL                    ::   l_star_ABL
REAL                    ::   l_star_DEF
REAL                    ::   UU_DEF_UU_ABL_fac
REAL                    ::   u_star_ABL
REAL                    ::   u_star_DEF
REAL                    ::   Shear_add_du_dz
REAL,ALLOCATABLE        ::   visc_wake(:,:)
REAL,ALLOCATABLE        ::   visc_wake1(:,:) 
REAL,ALLOCATABLE        ::   visc_wake2(:,:) 
REAL                    ::   visc_norm_factor
REAL,ALLOCATABLE        ::   alfa_1(:)
REAL,ALLOCATABLE        ::   alfa_2(:)
REAL,ALLOCATABLE        ::   du_dr_tot(:,:)
INTEGER,ALLOCATABLE     ::   shear_flag(:)
REAL,ALLOCATABLE        ::   One_div_du_dr_DWM(:,:)
REAL,ALLOCATABLE        ::   visc_fac(:)

REAL                    ::    R_WTG                        ! normalized radius
REAL                    ::    U0                           ! normalized wind speed
REAL                    ::    D_WTG                        ! normalized diameter
REAL                    ::    R_length                     ! normalized length in radial direction
REAL                    ::    X_length                     ! normalized length in axial direction
INTEGER                 ::    np_r                         ! point per radial distance
INTEGER                 ::    np_x                         ! point per axial distance
REAL                    ::    delrad                       ! delta r
REAL                    ::    delaxi                       ! delta x

REAL, ALLOCATABLE       ::    x_vector(:)
REAL, ALLOCATABLE       ::    r_vector(:)

REAL, ALLOCATABLE       ::    V(:,:)
REAL, ALLOCATABLE       ::    visc(:,:)
REAL, ALLOCATABLE       ::    visc_DWM(:,:)
REAL, ALLOCATABLE       ::    du_dr_DWM(:,:)
REAL, ALLOCATABLE       ::    du_dr_total(:,:)
REAL, ALLOCATABLE       ::    Turb_Stress_DWM(:,:)
!REAL, ALLOCATABLE       ::    TI_DWM(:,:)
!REAL, ALLOCATABLE       ::    U_face(:,:)
!REAL, ALLOCATABLE       ::    VOL_x_jhigh(:,:)
!REAL, ALLOCATABLE       ::    VOL_x_jlow (:,:)
!REAL, ALLOCATABLE       ::    VOL_r_ihigh(:,:)
!REAL, ALLOCATABLE       ::    VOL_r_ilow (:,:)
REAL, ALLOCATABLE       ::    r_vec_DWM (:)
REAL, ALLOCATABLE       ::    dA_DWM (:)

INTEGER                 ::    n_r_vec_DWM
INTEGER                 ::    b_loop
INTEGER                 ::    b_counter
REAL                    ::    dr_DWM
REAL                    ::    Def_DWM
REAL                    ::    Def_DWM_mixL
REAL                    ::    A_total
REAL                    ::    k_wiener

INTEGER, ALLOCATABLE    ::    counter(:)
INTEGER                 ::    i
INTEGER                 ::    j
INTEGER                 ::    k
INTEGER                 ::    ILo
INTEGER                 ::    NumEqu
INTEGER                 ::    n_xi
INTEGER                 ::    n_U_tmp_2

REAL, ALLOCATABLE       ::    bin_filter(:)
REAL, ALLOCATABLE       ::    xi(:)
REAL, ALLOCATABLE       ::    U_tmp_1(:)
REAL, ALLOCATABLE       ::    U_tmp_2(:)
REAL, ALLOCATABLE       ::    U_tmp(:)
REAL, ALLOCATABLE       ::    mat(:,:)
REAL, ALLOCATABLE       ::    RHS(:)
REAL, ALLOCATABLE       ::    Soln(:)
REAL, ALLOCATABLE       ::    AugMat(:,:)

REAL                    ::    LHS1
REAL                    ::    LHS2
REAL                    ::    LHS3
REAL                    ::    LHS11
REAL                    ::    LHS12
REAL                    ::    LHS13
REAL                    ::    LHS21
REAL                    ::    LHS22
REAL                    ::    LHS23
REAL                    ::    LHS31
REAL                    ::    LHS41
REAL                    ::    LHS32
REAL                    ::    LHS33
REAL                    ::    LHS43

REAL,ALLOCATABLE        ::    main_diagonal(:)
REAL,ALLOCATABLE        ::    sub_diagonal(:)
REAL,ALLOCATABLE        ::    sup_diagonal(:)


END MODULE DWM_Wake_Deficit_Data

!===============================================================================
MODULE filter_vector_data
!...............................................................................
! The filter subroutine data
!...............................................................................

REAL, ALLOCATABLE       ::    filter_vector_1(:)
REAL, ALLOCATABLE       ::    filter_1_vector(:)
REAL, ALLOCATABLE       ::    filter_2_vector(:)
REAL, ALLOCATABLE       ::    filter_3_vector(:)
REAL, ALLOCATABLE       ::    filter_4_vector(:)
REAL, ALLOCATABLE       ::    filter_5_vector(:)
REAL, ALLOCATABLE       ::    filter_6_vector(:)
REAL, ALLOCATABLE       ::    filter_vector_fill(:)
INTEGER                 ::    i
INTEGER                 ::    n_filter_1_vector
INTEGER                 ::    n_filter_2_vector
INTEGER                 ::    n_filter_3_vector
INTEGER                 ::    n_filter_4_vector
INTEGER                 ::    n_filter_5_vector
INTEGER                 ::    n_filter_6_vector
INTEGER                 ::    n_filter_vector_fill



END MODULE filter_vector_data

!===============================================================================
MODULE main_wake_calculation_data
!...............................................................................
! The DWM wake deficit calculation data
!...............................................................................

REAL, ALLOCATABLE       ::    V(:,:)
REAL, ALLOCATABLE       ::    visc(:,:)
REAL, ALLOCATABLE       ::    du_dr_DWM(:,:)
REAL, ALLOCATABLE       ::    Turb_Stress_DWM(:,:)
REAL, ALLOCATABLE       ::    TI_DWM(:,:)
REAL, ALLOCATABLE       ::    U_face(:,:)
REAL, ALLOCATABLE       ::    VOL_x_jhigh(:,:)
REAL, ALLOCATABLE       ::    VOL_x_jlow (:,:)
REAL, ALLOCATABLE       ::    VOL_r_ihigh(:,:)
REAL, ALLOCATABLE       ::    VOL_r_ilow (:,:)
REAL, ALLOCATABLE       ::    r_vec_DWM (:)
REAL, ALLOCATABLE       ::    dA_DWM (:)

INTEGER                 ::    n_r_vec_DWM
INTEGER                 ::    b_loop
INTEGER                 ::    b_counter
REAL                    ::    dr_DWM
REAL                    ::    Def_DWM
REAL                    ::    Def_DWM_mixL

INTEGER, ALLOCATABLE    ::    counter(:)
INTEGER                 ::    i
INTEGER                 ::    j
INTEGER                 ::    k
INTEGER                 ::    ILo
INTEGER                 ::    NumEqu
INTEGER                 ::    n_xi
INTEGER                 ::    n_U_tmp_2

REAL, ALLOCATABLE       ::    bin_filter(:)
REAL, ALLOCATABLE       ::    xi(:)
REAL, ALLOCATABLE       ::    U_tmp_1(:)
REAL, ALLOCATABLE       ::    U_tmp_2(:)
REAL, ALLOCATABLE       ::    U_tmp(:)
REAL, ALLOCATABLE       ::    mat(:,:)
REAL, ALLOCATABLE       ::    RHS(:)
REAL, ALLOCATABLE       ::    Soln(:)
REAL, ALLOCATABLE       ::    AugMat(:,:)
REAL                    ::    LHS1
REAL                    ::    LHS2
REAL                    ::    LHS3
REAL                    ::    LHS11
REAL                    ::    LHS12
REAL                    ::    LHS13
REAL                    ::    LHS21
REAL                    ::    LHS22
REAL                    ::    LHS23
REAL                    ::    LHS31
REAL                    ::    LHS41
REAL                    ::    LHS32
REAL                    ::    LHS33
REAL                    ::    LHS43

    END MODULE main_wake_calculation_data

!===============================================================================
MODULE shear_correction_data
!...............................................................................
! The shear correction data
!...............................................................................

REAL                    ::  alpha_1                       
REAL                    ::  alpha_2
INTEGER                 ::  temp_n
REAL                    ::  temp_integration
REAL                    ::  correction_factor
INTEGER                 ::  i
REAL,ALLOCATABLE        ::  alpha_array(:)
REAL                    ::  delta_alpha

END MODULE shear_correction_data
    
!===============================================================================
MODULE filter_velocity_data
!...............................................................................
! The filter velocity data
!...............................................................................

INTEGER      ::  number_counter               ! counter : how many points are in the circle
INTEGER      ::  radius_length                ! wake radius (meters)
INTEGER      ::  ErrStat
REAL         ::  y_axis                       
REAL         ::  z_axis
REAL         ::  temp_filter_velocity (3)     ! interpolation function, has u,v,w three components

REAL         ::  total_velocity = 0
INTEGER      ::  total_velocity_counter = 0


END MODULE filter_velocity_data

!===============================================================================
MODULE meandering_data
!...............................................................................
! The meandering subroutine data
!...............................................................................

INTEGER                 ::    scale_factor
INTEGER                 ::    ErrStat
INTEGER                 ::    release_time
INTEGER                 ::    flying_time
INTEGER                 ::    simulation_time_length
INTEGER                 ::    moving_time
REAL                    ::    DWM_time_step
REAL                    ::    temp_center_wake (3)
REAL                    ::    temp_velocity (3)
REAL                    ::    U_Scale_Factor
REAL                    ::    U_factor

END MODULE meandering_data

!===============================================================================
MODULE parameter_file_data
!...............................................................................
! The reading binary file subroutine data      
!...............................................................................

                        
REAL, ALLOCATABLE :: velocityU(:)                         ! the wake velocity profile @ the downstream turbine plane
REAL, ALLOCATABLE :: smoothed_wake(:)
REAL, ALLOCATABLE :: WakePosition(:,:,:)                  ! meandered wake center
INTEGER           :: WakePosition_1, WakePosition_2       ! size of the WakePosition
INTEGER           :: smooth_flag                          ! Whether or not use the smoothed out upstream wake profile (1-yes, 0-no)
INTEGER           :: p_p_r
INTEGER           :: NumWT                                ! The number of wind turbines in the wind farm
INTEGER           :: Tinfluencer
REAL              :: RotorR                               ! Rotor radius
!REAL              :: U_mean                               ! mean velocity of upstream turbine
REAL              :: r_domain
REAL              :: x_domain
REAL              :: Uambient                             ! The ambient wind velocity
REAL              :: TI_amb                               ! Ambient turbulence intensity (%)
REAL              :: TI_wake
REAL              :: hub_height
!REAL              :: spacing_turbine                         
REAL              :: length_velocityU
REAL              :: TurbRefHt                            ! The Turbsim wind file reference height
!REAL              :: AlignAngle                           ! The wind alignment angle along the row of turbine (left-positive,right-negative)
REAL              :: Wind_file_Mean_u                     ! The mean velocity of the first turbine
REAL              :: Winddir
REAL              :: WFLowerBd

REAL              :: next_hub_height                      ! downwind turbine hub height
INTEGER           :: ranW                                 ! Whether or not use random walk model


END MODULE parameter_file_data

!===============================================================================
MODULE read_turbine_position_data
!...............................................................................
! The reading turbine index subroutine data     
!...............................................................................
    
    CHARACTER(LEN=20)   ::  SimulationOrder_index_char
    INTEGER             ::  SimulationOrder_index
    INTEGER             ::  N_ARGU
    INTEGER,ALLOCATABLE ::  Turbine_sort_order(:)
    INTEGER             ::  WT_index                      ! wind turbine index in the wind farm
    INTEGER,ALLOCATABLE ::  TurbineInfluenceData(:,:)
    INTEGER,ALLOCATABLE ::  upwind_turbine_index(:)       ! the upwind turbines that affecting this turbine
    INTEGER,ALLOCATABLE ::  downwind_turbine_index(:)
    INTEGER             ::  upwindturbine_number          ! the number of upwind turbines affecting the downwind turbine
    INTEGER             ::  downwindturbine_number
    REAL,ALLOCATABLE    ::  turbine_windorigin_length(:)  
    REAL,ALLOCATABLE    ::  upwind_turbine_projected_distance(:) ! the projected distacne between two turbines
    REAL,ALLOCATABLE    ::  downwind_turbine_projected_distance(:)
    REAL,ALLOCATABLE    ::  turbine_angle(:,:)
    REAL,ALLOCATABLE    ::  upwind_align_angle(:)         ! the angle beween the line connecting the upwind turbine and this turbine and the wind direction vector
    REAL,ALLOCATABLE    ::  downwind_align_angle(:)
    REAL,ALLOCATABLE    ::  upwind_turbine_Xcoor(:)       ! the coordinate of the upwind turbine which affects this investigated turbine
    REAL,ALLOCATABLE    ::  upwind_turbine_Ycoor(:)
    REAL,ALLOCATABLE    ::  wind_farm_Xcoor(:)            ! the coordinates of all the turbines in the wind farm
    REAL,ALLOCATABLE    ::  wind_farm_Ycoor(:)
    REAL,ALLOCATABLE    ::  downwind_turbine_Xcoor(:)       ! the coordinate of the downwind turbine which is affected by this investigated turbine
    REAL,ALLOCATABLE    ::  downwind_turbine_Ycoor(:)
    
END MODULE read_turbine_position_data

!===============================================================================
MODULE weighting_method
!...............................................................................
! The weighting method subroutine data      
!...............................................................................
    
    REAL, ALLOCATABLE  ::  sweptarea(:)
    REAL               ::  weighting_denominator
    
    
END MODULE weighting_method
    
!===============================================================================
MODULE read_upwind_result_file_data
!...............................................................................
! The read upwind result file      
!...............................................................................
    REAL,ALLOCATABLE  ::  upwind_U(:,:)
    REAL,ALLOCATABLE  ::  upwind_wakecenter(:,:,:,:)
    REAL,ALLOCATABLE  ::  upwind_meanU(:)
    REAL,ALLOCATABLE  ::  upwind_TI(:)
    REAL,ALLOCATABLE  ::  upwind_small_TI(:)
    REAL,ALLOCATABLE  ::  upwind_smoothWake(:,:)
    REAL,ALLOCATABLE  ::  velocity_aerodyn(:)
    REAL,ALLOCATABLE  ::  TI_downstream(:)
    REAL,ALLOCATABLE  ::  small_scale_TI_downstream(:)
    REAL,ALLOCATABLE  ::  smoothed_velocity_array(:,:)
    REAL,ALLOCATABLE  ::  vel_matrix(:,:,:)             ! The smoothed out wake velocity matrix for n downwind turbine
    
END MODULE read_upwind_result_file_data

!===============================================================================
MODULE TI_downstream_data
!...............................................................................
! The TI calculation for downstream turbine subroutine data      
!...............................................................................
    
    REAL, ALLOCATABLE  ::  TI_downstream_matrix(:,:)
    INTEGER            ::  i,j,k
    INTEGER            ::  cross_plane_position_ds   ! the cross plane position which to be investigated in term of the flying time
    INTEGER            ::  cross_plane_position_TI   ! the cross plane position which to be investigated in term of the n_x_vector
    INTEGER            ::  distance_index            ! the index of the distance in the TI axisymmetric array
    INTEGER            ::  counter1
    INTEGER            ::  counter2
    INTEGER            ::  initial_timestep
    REAL               ::  y_axis_turbine
    REAL               ::  z_axis_turbine
    REAL               ::  distance                  ! the distance between one point to the meandered wake center
    REAL               ::  TI_downstream_node        ! the TI at a specfic point in the inbestigated cross plane
    REAL               ::  TI_node_temp
    REAL               ::  TI_node
    REAL(kind=8)       ::  TI_accumulation
    REAL(kind=8)       ::  TI_apprant_accumulation
    REAl               ::  TI_average                ! THE AVERAGE TI OF THE CROSS PLANE
    !REAL               ::  TI_total                  ! The Total TI
    REAL               ::  TI_apprant                ! The TI due to the meadering
    REAL               ::  HubHt 
    REAL               ::  wake_center_y
    REAL               ::  wake_center_z
    REAL               ::  Rscale
    REAL               ::  y
    REAL               ::  z  
    REAL               ::  zero_spacing
    REAL               ::  temp1,temp2,temp3   
    
    
END MODULE TI_downstream_data

!===============================================================================
MODULE Turbulence_KS
!...............................................................................
! The Kaimal spectrum subroutine data     
!...............................................................................
    
    INTEGER            ::   fs       ! sample frequency
    INTEGER            ::   temp_n
    INTEGER            ::   i
    REAL               ::   low_f    ! lower bound of frequency range
    REAL               ::   high_f   ! upper bound of frequency range
    REAL               ::   lk_facor ! turbulence length-scale
    REAL               ::   STD      ! standard deviation of the turbulence
    
    
END MODULE Turbulence_KS

!===============================================================================
MODULE shinozuka_data
!...............................................................................
! The shinozuka subroutine data      
!...............................................................................
    
    REAL,ALLOCATABLE  ::   f_syn(:)       ! frequency series
    REAL,ALLOCATABLE  ::   t_syn(:)       ! time series
    REAL,ALLOCATABLE  ::   phi(:)         ! random phase angle
    REAL,ALLOCATABLE  ::   p_k(:)         
    REAL,ALLOCATABLE  ::   a_k(:)
    INTEGER           ::   num_points     ! total number of points
    INTEGER           ::   ILo
    INTEGER           ::   i,j
    REAL              ::   dt             ! time step
    REAL              ::   t_min
    REAL              ::   t_max
    REAL              ::   df             ! frequency step

    
END MODULE shinozuka_data

!===============================================================================
MODULE smooth_out_wake_data
!...............................................................................
! The smooth_out_wake subroutine data     
!...............................................................................
    
    INTEGER,ALLOCATABLE  ::   counter_array(:)
    !REAL, ALLOCATABLE    ::   velocity_matrix(:,:)      ! the velocity matrix that store the velocity of the downstream turbine plane
    REAL,ALLOCATABLE     ::   velocity_array(:)
    INTEGER              ::   length_velocity_array     ! the length of velocity_array
    INTEGER              ::   low
    INTEGER              ::   high
    INTEGER              ::   i,j,k,n
    INTEGER              ::   counter
    REAL                 ::   y                         ! y coordinate
    REAL                 ::   z                         ! z coordinate
 
END MODULE smooth_out_wake_data

!===============================================================================
MODULE smooth_wake_shifted_velocity_data
!...............................................................................
! The smooth_wake_shifted_velocity function data     
!...............................................................................
    
    INTEGER    ::   p1
    INTEGER    ::   p2
    REAL       ::   distance                      ! the distance from the point to the meandered wake center
    REAL       ::   y0                            ! wake center position on y axis
    REAL       ::   z0                            ! wake center position on z axis
    REAL       ::   unit                          ! single unit length  R/ppR
  
END MODULE smooth_wake_shifted_velocity_data

    
!===============================================================================
MODULE RW_data
!...............................................................................
! Randow walk model subroutine data      
!...............................................................................

    REAL ( kind = 8 ) cdf    ! left region
    REAL ( kind = 8 ) cdf_right    ! right region
    REAL ( kind = 8 ) pdf
    REAL ( kind = 8 ) x
    REAL ( kind = 8 ) mean
    REAL ( kind = 8 ) bound
    REAL ( kind = 8 ) sd
    REAL ( kind = 8 ) cur
    REAL randNum
    REAL step
    REAL step_size  
    REAL(kind = 8), ALLOCATABLE :: res(:)
    REAL(kind = 8), ALLOCATABLE :: arr(:)
    REAL, ALLOCATABLE :: resTemp(:)     ! the wake center of the UPWIND turbine
    REAL, ALLOCATABLE :: curTime(:)
    INTEGER ( kind = 4 ) status
    INTEGER ( kind = 4 ) which
    INTEGER :: npoint
    INTEGER :: i
    INTEGER :: use_nth
    INTEGER :: size_use
    REAL  :: total
    REAL  :: res_std
    REAL  :: res_mean
    REAL  :: interpLoc
    REAL  :: ILo
    REAL  :: IHi
    INTEGER :: totalTime
    INTEGER :: J

END MODULE RW_data
MODULE DWM_main

USE   DWM_Wake_Sub

IMPLICIT NONE

    !......... Public Subroutines .................................................
    
PUBLIC :: calculate_DWM

CONTAINS
!----------------------------------------------------------------------------------------
SUBROUTINE calculate_DWM ( blade_radius, num_element, r_t  )
!..................................................................................
! This routine is called when the DWM model is asked to run
! There are two sub-modules in this routine
!----------------------------------------------------------
! The first sub-module is to calculate the wake velocity
!    which is the calculate_wake()
!    The output are wake_u, wake_width
!----------------------------------------------------------
! The second sub-module is to calculate the wake center position in the space
!    which is the Get_wake_center()
!    The output is wake_position
!.................................................................................. 
    
    USE   parameter_file_data,          ONLY: p_p_r, ranW
    USE   read_turbine_position_data,   ONLY: downwindturbine_number,downwind_turbine_projected_distance,downwind_align_angle
    USE   read_upwind_result_file_data, ONLY: TI_downstream, smoothed_velocity_array, vel_matrix, small_scale_TI_downstream
    USE   filter_velocity_data,         ONLY: total_velocity, total_velocity_counter
    
    INTEGER         ::           num_element
    REAL            ::           blade_radius
    REAL            ::           r_t ( num_element )
    INTEGER         ::           I

    
    !U_hubht = total_velocity/total_velocity_counter
    
    !Print*, U_hubht
  !----------------------------------------------------------------------------------------------------
      
    CALL calculate_mean_u( NElm, RELM(:), Mean_FFWS, TI)
    
         !OPEN (unit=25,file="DWM\results\TI.txt")
         !WRITE (25,'(f13.7)'), TI 
         
         !OPEN (unit=25,file="DWM\results\Mean_FFWS.txt")
         !WRITE (25,'(f13.7)'), Mean_FFWS
         !PRINT*, Mean_FFWS
         
         !OPEN(unit = 10, status='replace',file='DWM\results\upstream_meanU.bin',form='unformatted')    
         !WRITE(10)   Mean_FFWS                                                                                                                                                                                                                                 
         !CLOSE(10)
         
  !----------------------------------------------------------------------------------------------------
    CALL CalAtmShear( AtmUscale, du_dz_ABL, TI)
    
    
  !---------------------------------------------------------------------------------------------------- 
    
    CALL get_initial_condition ( induction_factor, RELM(:), NElm, r_initial, U_initial )
    
         !OPEN (unit=25,file="DWM\results\inducion_factor.txt")
         !WRITE (25,'(f13.7)'), induction_factor(:)
         
         
  !----------------------------------------------------------------------------------------------------
         
    CALL calculate_mean_genpwr ( Totalgenpwr, AP_timestep )
         
  !----------------------------------------------------------------------------------------------------
         
    CALL calculate_wake ( r_initial, U_initial, NElm, wake_u, wake_width )
    
         !OPEN (unit=25,file="DWM\results\U_wake.txt")
         !WRITE (25,'(f9.7)'), wake_u(:,:)
         
         !OPEN(unit = 10, status='replace',file='DWM\results\Uvelocity.bin',form='unformatted')          
         !WRITE(10)   wake_u(floor(spacing_turbine * p_p_r)+1,:)  ! write the wind data of the plane where the downstream turbine locates                                                                                                                                                                                                                                                                                                        
         !CLOSE(10)
         
         !!!OPEN (unit=25,file="DWM\results\wake_width.txt")
         !!!WRITE (25,'(I5)'), wake_width(:)
         !!!CLOSE(25)
         
         !IF (SimulationOrder_index==1) THEN
             !OPEN(unit = 10, status='replace',file='DWM\results\wake_width.bin',form='unformatted')    
             !WRITE(10)   wake_width(:)                        ! write the downstream wake width                                                                                                                                                                      
             !CLOSE(10)
         !END IF
    ! debug test
    
    !CALL Wake_time_advance (wake_u, wake_width)
  !-----------------------------------------------------------------------------------------------------    
    IF (downwindturbine_number >0 ) THEN
        CALL advection_time(wake_u,wake_width,downwind_turbine_projected_distance(1),travel_time)
    END IF
                
  !-----------------------------------------------------------------------------------------------------
    IF (ranW==0)   THEN     
        CALL Get_wake_center ( wake_width, wake_position )
    ELSE IF (ranW == 1)   THEN
        CALL Get_wake_center_RW (wake_position)
    END IF
    
         !OPEN(unit = 10, status='replace',file='DWM\results\meandered_wake_center.bin',form='unformatted')    
         !WRITE(10)   wake_position(:,:,:)                        ! write the downstream meandered wake center                                                                                                                                                                      
         !CLOSE(10)
         
         !OPEN (unit=25,file="DWM\results\wake_center.txt")
         !WRITE (25,'(f13.7)'), wake_position(:,:,:)
  !-----------------------------------------------------------------------------------------------------
    IF (downwindturbine_number >0 ) THEN
        DO I = 1,downwindturbine_number
           CALL smooth_out_wake(wake_u,wake_position,smoothed_velocity_array(I,:),downwind_turbine_projected_distance(I),downwind_align_angle(I),vel_matrix(I,:,:))
           TI_downstream (I)             = TI_downstream_total (downwind_turbine_projected_distance(I),downwind_align_angle(I),vel_matrix(I,:,:))
           small_scale_TI_downstream (I) = smallscale_TI (downwind_turbine_projected_distance(I),downwind_align_angle(I),vel_matrix(I,:,:))
        END DO
    END IF
    
    
         !OPEN(unit = 10, status='replace',file='DWM\results\smoothed_vlocity.bin',form='unformatted')          
         !WRITE(10)   smoothed_velocity_array(:)          ! write the wind data of the plane where the downstream turbine locates,                                                                                                                                                                                                                                                                
         !CLOSE(10)
        
  !------------------------------------------------------------------------------------------------------
    
    !CALL calculate_TI_node_downstream(TI_downstream)
    
         !OPEN(unit = 10, status='replace',file='DWM\results\Downstream_TI_b4normalization.bin',form='unformatted')          
         !WRITE(10)   TI_downstream                                                                                                                                                                                                                                                                                                         
         !CLOSE(10)
         
    CALL write_result_file()

         

         
    
END SUBROUTINE calculate_DWM

END MODULE DWM_main
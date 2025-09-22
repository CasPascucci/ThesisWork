    SUBROUTINE Targeting_Beyond(R_I, V_I, gamma, kr, g0, R_END, V_END, AT_END, &
                    Delta_t, t_go_true, r_f, V_f, a_f, t_go_beyond)

!==================================================================================
!   Find the targeting condition for Fractional-Polynominal powered descent guidance (FP2DG)
!   after a period of time Delta t beyond  the termination point of the current phase, 
!   so that the trajectory achieves the required final condition at the termination where
!   t_go = Delta t>0 (as opposed to at t_go =0 which will cause a singlarity in the FP2DG law).
!   Delta t can have a wide range of values. It is suggested that Delta=10 (s) as a default. 
!       See P. Lu, "The Theory of Fractional-Polynominal Powered Descent Guidance", to be presented
!       at 2020 GN&C Conference, and to appear in J. Guidance, Control, and Dynamics. 
!
!       This subroutine assumes that all the vectors (inputs and outputs) are expressed in the same
!       Cartesian frame (henceforth referred to as the "guidance frame") 
! Inputs: 
!   R_I(3)              Current position vector in the guidance frame                           (m)
!   V_I(3)              Current velocity vector in the guidance frame                           (m/s)
!   gamma               A parameter for FP2DG law (gamma>0)
!   kr                  Another parameter for FP2DG law (kr>=2(gamma+2))
!   g0                  Gravitational acceleration at the equator of the planet                 (m/s^2) 
!   R_END(3)            The required final position vector in the guidance frame                (m)
!   V_END(3)            The required final velocity vector in the guidance frame                (m/s)
!   AT_END(3)           The required final thrust acceleration vector in the guidance frame     (m/s^2)
!   Delta_t             The time increment beyond the landing site (a constant)                 (s)
!   t_go_true           Time-to-go to the required targeting condition                          (s)
!                       
!
! Outputs:
!   r_f(3)              The targeting position vector in the guidance frame used in FP2DG law   (m)
!   V_f(3)              The targeting velocity vector in the guidance frame used in FP2DG law   (m/s)                   
!   a_f(3)              The targeting acceleration vector in the guidance frame used in FP2DG law(m/s^2)
!   t_go_beyond         The time-to-go to the above output targeting condition                  (s)
!
!  Ping Lu
!  July 17, 2019
!  Oct. 5, 2019
!=================================================================================
    IMPLICIT NONE
	REAL(8), DIMENSION (3), INTENT(IN) :: R_I	
    REAL(8), DIMENSION (3), INTENT(IN) :: V_I
    REAL(8), INTENT(IN) :: gamma
    REAL(8), INTENT(IN) :: kr
    REAL(8), INTENT(IN) :: g0
	REAL(8), DIMENSION (3), INTENT(IN) :: R_END	
    REAL(8), DIMENSION (3), INTENT(IN) :: V_END
	REAL(8), DIMENSION (3), INTENT(IN) :: AT_END
    REAL(8), INTENT(IN) :: Delta_t
    REAL(8), INTENT(IN) :: t_go_true
    
    REAL(8), DIMENSION (3), INTENT(OUT) :: r_f	
    REAL(8), DIMENSION (3), INTENT(OUT) :: V_f
	REAL(8), DIMENSION (3), INTENT(OUT) :: a_f
    REAL(8), INTENT(OUT) :: t_go_beyond    
    
    REAL(8) :: phi1, phi2, phi1_bar, phi2_bar, phi1_hat, phi2_hat
    REAL(8) :: Delta, gamma1, gamma2,r, t_go
    REAL(8) :: k1r, k1v, k1a, k2r, k2v, k2a
    REAL(8), DIMENSION (3) :: g, b1, b2, b3, d1, d2
    REAL(8), DIMENSION (3,3) :: a
    REAL(8), DIMENSION (9,9) :: M
    REAL(8), DIMENSION (9) :: b,x,x_star
    REAL(8), DIMENSION (6,6) :: MM
    REAL(8), DIMENSION (6) :: bb
    INTEGER :: I, J, N, mp

    ! Check the input parameters
    IF(gamma <0.0D0) THEN
        print*, "Warning: Incorrect value for gamma in subroutine Targeting_Beyond"
        STOP 1203
    ELSE IF (kr<2.0D0*(gamma+2.0D0)) THEN
        print*, "Warning: Incorrect value for kr in subroutine Targeting_Beyond"
        STOP 1204
    END IF
    
     t_go_beyond = t_go_true + Delta_t
     t_go = t_go_beyond
     
    IF(Delta_t<1.0D-15) THEN  ! No beyond-termination targeting is done
        t_go_beyond = t_go_true
        r_f     = R_END
        V_f     = V_END
        a_f     = AT_END
        RETURN
    END IF
    
    r   = Dsqrt(R_I(1)**2+R_I(2)**2+R_I(3)**2)      ! radius
    DO I = 1, 3
        g(I) = -R_I(I)*g0/r                         ! gravity vector in gudiance frame
    END DO
    gamma1      = gamma
    gamma2      = kr/(gamma+2.0D0) -2.0D0

    IF(DABS(gamma1-gamma2)<1.0D-15) THEN    ! The case where kr=(gamma+2)^2
        t_go_beyond = t_go_true             ! Do not do beyond-termination targeting
        r_f     = R_END
        V_f     = V_END
        a_f     = AT_END
        RETURN
    END IF
    
    phi1        = t_go**gamma1
    phi2        = t_go**gamma2
    IF(gamma2<=1.0D-15 .AND. t_go<1.D-15) THEN
        phi2 = 0.0D0
    END IF
    phi1_bar    = -(1.0D0/(gamma1+1.0D0))*t_go**(gamma1+1.0D0)
    phi2_bar    = -(1.0D0/(gamma2+1.0D0))*t_go**(gamma2+1.0D0)
    phi1_hat    = (1.0D0/((gamma1+1.0D0)*(gamma1+2.0D0)))*t_go**(gamma1+2.0D0)
    phi2_hat    = (1.0D0/((gamma2+1.0D0)*(gamma2+2.0D0)))*t_go**(gamma2+2.0D0)
    Delta       = phi1_hat*phi2_bar-phi2_hat*phi1_bar
    
    k1r =-phi2_bar/Delta
    k1v =(phi2_bar*t_go+phi2_hat)/Delta
    k1a =-(0.5D0*t_go*phi2_bar+phi2_hat)*t_go/Delta
    d1  =(1.0D0/Delta)*(phi2_bar*R_I-0.5D0*phi2_bar*t_go**2*g-phi2_hat*V_I-phi2_hat*t_go*g)
    k2r =phi1_bar/Delta
    k2v =-(phi1_bar*t_go+phi1_hat)/Delta
    k2a =(0.5D0*t_go*phi1_bar+phi1_hat)*t_go/Delta
    d2  =-(1.0D0/Delta)*(phi1_bar*R_I-0.5D0*phi1_bar*t_go**2*g-phi1_hat*V_I-phi1_hat*t_go*g)

! Note the following quantities are re-defined with different time (Delta_t, not t_go)
    phi1        = Delta_t**gamma1
    phi2        = Delta_t**gamma2
    IF(gamma2<=1.0D-14 .AND. Delta_t<1.D-14) THEN
        phi2 = 0.0D0
    END IF
    phi1_bar    = -(1.0D0/(gamma1+1.0D0))*Delta_t**(gamma1+1.0D0)
    phi2_bar    = -(1.0D0/(gamma2+1.0D0))*Delta_t**(gamma2+1.0D0)
    phi1_hat    = (1.0D0/((gamma1+1.0D0)*(gamma1+2.0D0)))*Delta_t**(gamma1+2.0D0)
    phi2_hat    = (1.0D0/((gamma2+1.0D0)*(gamma2+2.0D0)))*Delta_t**(gamma2+2.0D0)
    
    a(1,1)  = k1r*phi1_hat+k2r*phi2_hat+1.0D0
    a(1,2)  = k1v*phi1_hat+k2v*phi2_hat-Delta_t
    a(1,3)  = k1a*phi1_hat+k2a*phi2_hat+0.5D0*Delta_t**2
    b1      = phi1_hat*d1+phi2_hat*d2+0.5D0*Delta_t**2*g
    a(2,1)  = k1r*phi1_bar+k2r*phi2_bar
    a(2,2)  = k1v*phi1_bar+k2v*phi2_bar+1.0D0
    a(2,3)  = k1a*phi1_bar+k2a*phi2_bar-Delta_t
    b2      = phi1_bar*d1+phi2_bar*d2-Delta_t*g
    a(3,1)  = k1r*phi1+k2r*phi2
    a(3,2)  = k1v*phi1+k2v*phi2
    a(3,3)  = k1a*phi1+k2a*phi2+1.0D0
    b3      = phi1*d1+phi2*d2
    ! Define the vector b in Mx=b   
    DO I= 1, 3
        x_star(I)   = R_END(I)
        x_star(I+3) = V_END(I)
        x_star(I+6) = AT_END(I)
        b(I)        = b1(I)
        b(I+3)      = b2(I)
        b(I+6)      = b3(I)
    END DO
    b   = x_star -b
    ! Define the matrix in Mx=b
    M=0.0D0
    DO I =1, 3
        M(I, I)     = a(1,1)
        
        M(I, I+3)   = a(1,2)
        M(I, I+6)   = a(1,3)
        M(I+3,I)    = a(2,1)
        M(I+3,I+3)  = a(2,2)
        M(I+3,I+6)  = a(2,3)
        M(I+6,I)    = a(3,1)
        M(I+6,I+3)  = a(3,2)
        M(I+6,I+6)  = a(3,3)
    END DO
    
    ! Define the dimension of the linear algebraic system
    IF(gamma2 >0.0D0) THEN
        N   = 9
        mp  = 1
        CALL gaussj(M,N,N,b,mp,mp)
        ! The solution to Mx=b is now contained in b; assign the outputs:
        DO I = 1, 3
            r_f(I)  = b(I)
            V_f(I)  = b(I+3)
            a_f(I)  = b(I+6)
        END DO
    ELSE            ! Degenerate cases where either gamma1=0 or gamma1=gamma
        DO I=1,6
            DO J=1, 6
                MM(I,J)= M(I,J)
            END DO
            bb(I) = b(I)
        END DO
        N   = 6
        mp  = 1
        CALL gaussj(MM,N,N,bb,mp,mp)
        ! The solution to MMx=bb is now contained in b; assign the outputs:
        DO I = 1, 3
            r_f(I)  = bb(I)
            V_f(I)  = bb(I+3)
            a_f(I)  = AT_END(I)    ! Just to assign some value
        END DO
    END IF
    
    RETURN
    END SUBROUTINE Targeting_Beyond
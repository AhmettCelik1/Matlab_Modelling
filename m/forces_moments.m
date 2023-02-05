% forces_moments.m
%   Computes the forces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%

function out = forces_moments(x, delta, wind, P, MAV)

    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis
    
    % compute wind data in NED
    w_n = w_ns - u_wg;
    w_e = w_es - v_wg;
    w_d = w_ds - w_wg;
    
    % compute air data
    Va = sqrt((u - w_n)*(u - w_n) + (v - w_e)*(v - w_e) + (w - w_d)*(w - w_d));
    alpha = atan((w-w_d)/(u-w_n));
    beta = asin((v-w_e)/Va);
    
    C_x_alpha = -0.030*cos(alpha) + MAV.C_L_alpha*sin(alpha);
    C_x_q_alpha = -MAV.C_D_q*cos(alpha) + MAV.C_L_q*sin(alpha);
    C_x_delta_e_alpha = -MAV.C_D_delta_e*cos(alpha) + MAV.C_L_delta_e*sin(alpha);
    C_z_alpha = -0.030*sin(alpha) - MAV.C_L_alpha*cos(alpha);
    C_z_q_alpha = -MAV.C_D_q*sin(alpha) - MAV.C_L_q*cos(alpha);  
    C_z_delta_e_alpha = -MAV.C_D_delta_e*sin(alpha) -  MAV.C_L_delta_e*cos(alpha);
    n = 10;

    J = Va / (n * MAV.D_prop);
    C_T_J = MAV.C_T2 * J * J + MAV.C_T1 * J + MAV.C_T0;
    C_T_Q = MAV.C_Q2 * J * J + MAV.C_Q1 * J + MAV.C_Q0;
    T_p = MAV.rho * n * n * MAV.D_prop^4 * C_T_J; 
    Q_p = MAV.rho * n * n * MAV.D_prop^5 * C_Q_J;
    % compute external forces and torques on aircraft
    Force(1) =  -m*g*sin(theta) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(C_x_alpha + C_x_q_alpha * MAV.c * q/(2*Va)) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(C_x_delta_e_alpha*delta_e) + (T_p*delta_t*Va);%tp sıkıntı
    Force(2) =  m*g*cos(theta)*sin(psi) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.C_Y_0 + MAV.C_Y_beta*beta + MAV.C_Y_p*MAV.b*p/(2*Va) + MAV.C_Y_r*MAV.b*r/(2*Va)) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.C_Y_delta_a*delta_a + MAV.C_Y_delta_r*delta_r); 
    Force(3) =  m*g*cos(theta)*cos(psi) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(C_z_alpha + C_z_q_alpha * MAV.c * q/(2*Va) ) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(C_z_delta_e_alpha*delta_e);
    
    Torque(1) = (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.b*(MAV.C_L_0 + MAV.C_ell_beta*beta + MAV.C_ell_p*MAV.b*p/(2*Va) + MAV.C_ell_r*MAV.b*r/(2*Va))) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(b*(MAV.C_ell_delta_a*delta_a + MAV.C_ell_delta_r*delta_r) +  + Q_p );%Qp eksik
    Torque(2) = (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.c*(MAV.C_m_0 + MAV.C_m_alpha*alpha + MAV.C_m_q*MAV.c*q/(2*Va))) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.c*(MAV.C_m_delta_e*delta_e));   
    Torque(3) = (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.b*(MAV.C_n_0 + MAV.C_n_beta*beta + MAV.C_n_p*MAV.b*p/(2*Va) + MAV.C_n_r*MAV.b*r/(2*Va) )) + (1/2)*MAV.rho*Va*Va*MAV.S_wing*(MAV.b*(MAV.C_n_delta_a*delta_a + MAV.C_n_delta_r*delta_r));
   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
end




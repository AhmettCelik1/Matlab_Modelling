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

function out = forces_moments(x, delta, wind, P)

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
    %w_n = 0;
    %w_e = 0;
    %w_d = 0;
    
    % compute air data
    %Va = 0;
    %alpha = 0;
    %beta = 0;
    
    % compute external forces and torques on aircraft
    %Force(1) =  0;
    %Force(2) =  0;
    %Force(3) =  0;
    
    %Torque(1) = 0;
    %Torque(2) = 0;   
    %Torque(3) = 0;
  % compute wind data in NED
    w_n = w_ns - u_wg;
    w_e = w_es - v_wg;
    w_d = w_ds - w_wg;
    
    % compute air data
    Va = sqrt((u - w_n)*(u - w_n) + (v - w_e)*(v - w_e) + (w - w_d)*(w - w_d));
    alpha = atan((w-w_d)/(u-w_n));
    beta = asin((v-w_e)/Va);
    
    C_x_alpha = -0.030*cos(alpha) + 5.61*sin(alpha);
    C_x_q_alpha = -0.0*cos(alpha) + 7.95*sin(alpha);
    C_x_delta_e_alpha = -0.0135*cos(alpha) + 0.13*sin(alpha);
    C_z_alpha = -0.030*sin(alpha) - 5.61*cos(alpha);
    C_z_q_alpha = -0.0*sin(alpha) - 7.95*cos(alpha);  
    C_z_delta_e_alpha = -0.0135*sin(alpha) -  0.13*cos(alpha);
    n = 10;

    J = Va / (n * 0.508);
    C_T_J = -0.1079 * J * J + -0.06044 * J + 0.09357;
    C_Q_J = -0.01664 * J * J + 0.004970 * J + 0.005230;
    T_p = 1.2682 * n * n * 0.508^4 * C_T_J; 
    Q_p = 1.2682 * n * n * 0.508^5 * C_Q_J;
    % compute external forces and torques on aircraft
    Force(1) =  -11.0*9.81*sin(theta) + (1/2)*1.2682*Va*Va*0.55*(C_x_alpha + C_x_q_alpha * 0.19 * q/(2*Va)) + (1/2)*1.2682*Va*Va*0.55*(C_x_delta_e_alpha*delta_e) + (T_p*delta_t*Va);%tp sıkıntı
    Force(2) = 11.0*9.81*cos(theta)*sin(psi) + (1/2)*1.2682*Va*Va*0.55*(0.0 + -0.98*beta + 0.0*2.90*p/(2*Va) + 0.0*2.90*r/(2*Va)) + (1/2)*1.2682*Va*Va*0.55*(0.075*delta_a + 0.19*delta_r); 
    Force(3) = 11.0*9.81*cos(theta)*cos(psi) + (1/2)*1.2682*Va*Va*0.55*(C_z_alpha + C_z_q_alpha * 0.19 * q/(2*Va) ) + (1/2)*1.2682*Va*Va*0.55*(C_z_delta_e_alpha*delta_e);
    
    Torque(1) = (1/2)*1.2682*Va*Va*0.55*(2.90*(0.23 + -0.13*beta + -0.51*2.90*p/(2*Va) + 0.25*2.90*r/(2*Va))) + (1/2)*1.2682*Va*Va*0.55*(2.90*(0.17*delta_a + 0.0024*delta_r) +  + Q_p );%Qp eksik
    Torque(2) = (1/2)*1.2682*Va*Va*0.55*(0.19*(0.0135 + -2.74*alpha + -38.21*0.19*q/(2*Va))) + (1/2)*1.2682*Va*Va*0.55*(0.19*(-0.99*delta_e));   
    Torque(3) = (1/2)*1.2682*Va*Va*0.55*(2.90*(0.0 + 0.073*beta + -0.069*2.90*p/(2*Va) + -0.095*2.90*r/(2*Va) )) + (1/2)*1.2682*Va*Va*0.55*(2.90*(-0.011*delta_a + -0.069*delta_r));
   

   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
    fprintf("VAA");
    disp(Va);
end




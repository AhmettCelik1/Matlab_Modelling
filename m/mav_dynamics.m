function [sys,x0,str,ts,simStateCompliance] = mav_dynamics(t,x,u,flag,MAV)
%function [sys,x0,str,ts,simStateCompliance] = mav_dynamics(t,x,u,flag,MAV)
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(MAV);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=mdlDerivatives(t,x,u,MAV);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(MAV)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 13;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [...
    MAV.pn0;...
    MAV.pe0;...
    MAV.pd0;...
    MAV.u0;...
    MAV.v0;...
    MAV.w0;...
    MAV.e0;...
    MAV.e1;...
    MAV.e2;...
    MAV.e3;...
    MAV.p0;...
    MAV.q0;...
    MAV.r0;...
    ];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,uu, MAV)

    pn    = x(1);
    pe    = x(2);
    pd    = x(3);
    u     = x(4);
    v     = x(5);
    w     = x(6);
    e0    = x(7);
    e1    = x(8);
    e2    = x(9);
    e3    = x(10);
    p     = x(11);
    q     = x(12);
    r     = x(13);
    fx    = uu(1);
    fy    = uu(2);
    fz    = uu(3);
    ell   = uu(4);
    m     = uu(5);
    n     = uu(6);
    
    pndot = (e1*e1 + e0*e0- e2*e2- e3*e3)*u + 2*v*(e1*e2 - e3*e0) + 2*w*(e1*e3 + e2*e0);
    pedot = 2*u*(e1*e2 + e3*e0) + v*(e2*e2 + e0*e0 - e1*e1 - e3*e3) + 2*w*(e2*e3 - e1*e0);
    pddot = 2*u*(e1*e3 - e2*e0) + 2*v*(e2*e3 + e1*e0) + w*(e3*e3 + e0*e0 - e1*e1 - e2*e2); 
    udot = r*v - q*w + fx/(m+0.0000001);  
    vdot = p*w - r*u + fy/(m+0.0000001);
    wdot = q*u - p*v + fz/(m+0.0000001);
    e0dot = (-p*e1 -q*e2 -r*e3)/2;
    e1dot = (p*e0 + r*e2 -q*3)/2;
    e2dot = (q*e0 -r*e1 +p*e3)/2;
    e3dot = (r*e0 +q*e1 -p*e2)/2;
    pdot = MAV.Gamma1*p*q - MAV.Gamma2*q*r + MAV.Gamma3*ell + MAV.Gamma4*n;
    qdot = MAV.Gamma5*p*r - MAV.Gamma6*(p*p-r*r) + m*(1/MAV.Jy);
    rdot = MAV.Gamma7*p*q - MAV.Gamma1*q*r + MAV.Gamma4*ell + MAV.Gamma8*n;

sys = [pndot; pedot; pddot; udot; vdot; wdot; e0dot; e1dot; e2dot; e3dot; pdot; qdot; rdot];
sys(isnan(sys))=0;

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [ ];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x)
   
    quaternion= [x(7) x(8) x(9) x(10) ]
    eulXYZ = quat2eul(quaternion,'XYZ')
    phi = eulXYZ(1)
    theta = eulXYZ(2)
    psi = eulXYZ(3)

    y = [x(1),
    x(2),
    x(3),
    x(4),
    x(5),
    x(6),
    phi,%phi
    theta,%theta
    psi,%psi
    x(11),
    x(12),
    x(13),
  ];
%phi theta psi    
sys = y;

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate

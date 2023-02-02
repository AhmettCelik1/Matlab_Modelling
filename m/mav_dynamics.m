function [sys,x0,str,ts,simStateCompliance] = mav_dynamics(t,x,u,flag,MAV)

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
    MAV.Count =  MAV.Count + 1;
    %fprintf('The following ID does not have an assigned stress value: %d\n',MAV.Count)
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
    %disp(x,uu)
    fprintf('pn value: %d\n',pn)
    fprintf('fx value: %d\n',fx)
    %fprintf('pe value: %d\n',pe)
    %fprintf('pd value: %d\n',pd)
    %fprintf('u value: %d\n',u)
    %fprintf('v value: %d\n',v)
    %fprintf('w value: %d\n',w)


    pndot = (e1*e1 + e0*e0- e2*e2- e3*e3)*u + 2*v*(e1*e2 - e3*e0) + 2*w*(e1*e3 + e2*e0);
    %fprintf('pndot value: %d\n',pndot)

    pedot = 2*u*(e1*e2 + e3*e0) + v*(e2*e2 + e0*e0 - e1*e1 - e3*e3) + 2*w*(e2*e3 - e1*e0);
    %fprintf('pedot value: %d\n',pedot)
    pddot = 2*u*(e1*e3 - e2*e0) + 2*v*(e2*e3 + e1*e0) + w*(e3*e3 + e0*e0 - e1*e1 - e2*e2); 
    %fprintf('pddot value: %d\n',pddot)
    udot = r*v - q*w + fx/(m+0.0000001);
    fprintf('m value: %d\n',m)
    %fprintf("r*v: %d\n",r*v)
    %fprintf("q*w: %d\n",q*w)

    %fprintf("fx/m: %d\n",fx/m)
    vdot = p*w - r*u + fy/(m+0.0000001);
    %fprintf('vdot value: %d\n',vdot)
    wdot = q*u - p*v + fz/(m+0.0000001);
    %fprintf('wdot value: %d\n',wdot)
    %wdot = 0;
    %fprintf('wdot value: %d\n',wdot)
    %udot = 0;
    %fprintf('udot value: %d\n',udot)
    %vdot = 0;
    %fprintf('vdot value: %d\n',vdot)
    e0dot = (-p*e1 -q*e2 -r*e3)/2;
    %fprintf('e0dot value: %d\n',e0dot)
    e1dot = (p*e0 + r*e2 -q*3)/2;
    %fprintf('e1dot value: %d\n',e1dot)
    e2dot = (q*e0 -r*e1 +p*e3)/2;
    %fprintf('e2dot value: %d\n',e2dot)
    e3dot = (r*e0 +q*e1 -p*e2)/2;
    %fprintf('e3dot value: %d\n',e3dot)
    pdot = MAV.Gamma1*p*q - MAV.Gamma2*q*r + MAV.Gamma3*ell + MAV.Gamma4*n;
    %fprintf('pdot value: %d\n',pdot)
    qdot = MAV.Gamma5*p*r - MAV.Gamma6*(p*p-r*r) + m*(1/MAV.Jy);
    %fprintf('qdot value: %d\n',qdot)
    rdot = MAV.Gamma7*p*q - MAV.Gamma1*q*r + MAV.Gamma4*ell + MAV.Gamma8*n;
    %fprintf('rdot value: %d\n',rdot)

    %pn, pn_hat, pn_c
    %fprintf('pn: %d\n',uu(1))
    %fprintf('pn_hat: %d\n',uu(31))
    %fprintf('pn_c: %d\n',uu(19))
sys = [pndot; pedot; pddot; udot; vdot; wdot; e0dot; e1dot; e2dot; e3dot; pdot; qdot; rdot];
sys(isnan(sys))=0;
disp(sys)
%disp(sys)
% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x)
    y = [...
        ];
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

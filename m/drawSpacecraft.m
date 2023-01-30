function drawSpacecraft(uu)

    % process inputs to function
    pn       = uu(1);       % inertial North position     
    pe       = uu(2);       % inertial East position
    pd       = uu(3);           
    u        = uu(4);       
    v        = uu(5);       
    w        = uu(6);       
    phi      = uu(7);       % roll angle         
    theta    = uu(8);       % pitch angle     
    psi      = uu(9);       % yaw angle     
    p        = uu(10);       % roll rate
    q        = uu(11);       % pitch rate     
    r        = uu(12);       % yaw rate    
    t        = uu(13);       % time

    % define persistent variables 
    persistent spacecraft_handle;
    persistent Vertices
    persistent Faces
    persistent facecolors
    
    % first time function is called, initialize plot and persistent vars
    if t==0
        figure(1), clf
        [Vertices, Faces, facecolors] = defineSpacecraftBody;
        spacecraft_handle = drawSpacecraftBody(Vertices,Faces,facecolors,...
                                               pn,pe,pd,phi,theta,psi,...
                                               [],'normal');
        title('Spacecraft')
        xlabel('East')
        ylabel('North')
        zlabel('-Down')
        view(32,47)  % set the vieew angle for figure
        axis([-10,10,-10,10,-10,10]);
        hold on
        
    % at every other time step, redraw base and rod
    else 
        drawSpacecraftBody(Vertices,Faces,facecolors,...
                           pn,pe,pd,phi,theta,psi,...
                           spacecraft_handle);
    end
end

  
%=======================================================================
% drawSpacecraft
% return handle if 3rd argument is empty, otherwise use 3rd arg as handle
%=======================================================================
%
function handle = drawSpacecraftBody(V,F,patchcolors,...
                                     pn,pe,pd,phi,theta,psi,...
                                     handle,mode)
  V = rotate(V', phi, theta, psi)';  % rotate spacecraft
  V = translate(V', pn, pe, pd)';  % translate spacecraft
  % transform vertices from NED to XYZ (for matlab rendering)
  R = [...
      0, 1, 0;...
      1, 0, 0;...
      0, 0, -1;...
      ];
  V = V*R;
  
  if isempty(handle)
  handle = patch('Vertices', V, 'Faces', F,...
                 'FaceVertexCData',patchcolors,...
                 'FaceColor','flat',...
                 'EraseMode', mode);
  else
    set(handle,'Vertices',V,'Faces',F);
    drawnow
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
function XYZ=rotate(XYZ,phi,theta,psi)
  % define rotation matrix
  R_roll = [...
          1, 0, 0;...
          0, cos(phi), -sin(phi);...
          0, sin(phi), cos(phi)];
  R_pitch = [...
          cos(theta), 0, sin(theta);...
          0, 1, 0;...
          -sin(theta), 0, cos(theta)];
  R_yaw = [...
          cos(psi), -sin(psi), 0;...
          sin(psi), cos(psi), 0;...
          0, 0, 1];
  R = R_roll*R_pitch*R_yaw;
  % rotate vertices
  XYZ = R*XYZ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate vertices by pn, pe, pd
function XYZ = translate(XYZ,pn,pe,pd)
  XYZ = XYZ + repmat([pn;pe;pd],1,size(XYZ,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define spacecraft vertices and faces
function [V,F,colors] = defineSpacecraftBody()
    % Define the vertices (physical location of vertices
    V = [...
        1    1    0;... % point 1
        1   -1    0;... % point 2
        -5   -1    -1;... % point 3
        -1    1    0;... % point 4
        1    1   -2;... % point 5
        1   -1   -2;... % point 6
        -1   -1   -2;... % point 7
        -1    1   -2;... % point 8
        1.5  1.5  0;... % point 9
        1.5 -1.5  0;... % point 10
        -1.5 -1.5  0;... % point 11
        -1.5  1.5  0;... % point 12
        2    0    -1;... % point 13
        0 -8 10;...%point14
        -5 1 -1;... %point15
        -7 0 -1;...%point 16
        -1 -3 -1;...%point17
        -1 3 -1;...%point18
        -3 3 -1;...%point19
        -3 -3 -1;...%point20 
       
        -7 -2 -1;...%point21
        -7 2 -1;...%point22
        -5 2 -1;...%point23
        -5 -2 -1;...%point24 
                
        -7, 0, -4;...%point25
         -7 0 -1;...%point 26
         -5,0, -1;...%point 27
    ];

    % define faces as a list of vertices numbered above
    F = [...
        1, 2,  6,  5, 13,13;...  % front
        14, 14,  14,  14, 14,14;...  % back
        1, 5,  16,  16, 16,16;...  % right 
        2, 6,  16,  16, 16,16;...  % left
        5, 6, 16, 16, 16, 16;...  % top
        1, 2, 16, 16, 16,16;... % bottom
        17, 18, 19, 20, 20,20;...%wing
        21, 22, 23, 24, 24, 24;...%tailiş
        26, 27, 25, 25, 25, 25;...%rt
        ];

    % define colors for each face    
    myred = [1, 0, 0];
    mygreen = [0, 1, 0];
    myblue = [0, 0, 1];
    myyellow = [1, 1, 0];
    mycyan = [0, 1, 1];

    colors = [...
        myblue;...    % front
        myblue;...  % back
        myblue;...   % right
        myblue;... % left
        myblue;...   % top
        myblue;...   % bottom
        mygreen;...
        mygreen;...
        mygreen;...
        ];
end
  

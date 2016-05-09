% Generates a Cart3D mesh and input.cntl file and calls a script that runs
%   Cart3D
%
%       M - freestream Mach number
%       alpha - freestream angle of attack
%       x - vector of x/c points
%       top - vector of points on the airfoil upper surface y/c
%       bot - vector of points on the airfoil lower surface y/c
%       foilname - mesh file name desired
%
%       Cl - lift coefficient
%       Cd - drag coefficient
%               Cd=-1 if Cart3D failed
function [Cl_Cart,Cd_Cart]=CartMesh(M,alpha,Xin,top,bot,foilname)
% Mesh a 2.5D wing % this is the newer one
span=1;
chord=1;
comp=1; % component number in mesh (1 is default)
compute_norms=false;

thick=top-bot;
if(any(thick<0))
    Cl_Cart=0;
    Cd_Cart=-1;
    disp('Cart Case Failed-thickness')
else
    % % For a NACA 4-Series Airfoil
    % foilname='0012';
    % numpts=50; for NACA 4-series
    % numsec=10; % number spanwise sections
    % mc=str2double(foilname(1));
    % pc=str2double(foilname(2));
    % tc=str2double(foilname(3:4));
    % [x,z]=NACA4_2(mc,pc,tc,numpts);

    % to use top, bot, and x
    x=[fliplr(Xin),Xin(2:end)];
    z=[fliplr(top),bot(2:end)];
    z(length(z))=(top(end)*0.9+0.1*bot(length(bot)-1));
    numpts=length(x);
    numsec=floor(numpts/10); % number spanwise sections

    nodes=zeros(numpts*numsec+1,3);
    ys=linspace(-span/2,span/2,numsec+1);
    for i=1:numsec+1
        % standard direction:
        nodes((i-1)*numpts+1:i*numpts,1)=x*chord;
        nodes((i-1)*numpts+1:i*numpts,2)=ys(i);
        nodes((i-1)*numpts+1:i*numpts,3)=z*chord;
        % Cart3D directions
        nodes((i-1)*numpts+1:i*numpts,1)=x*chord;
        nodes((i-1)*numpts+1:i*numpts,2)=z*chord;
        nodes((i-1)*numpts+1:i*numpts,3)=ys(i);
    end
    quads=zeros(numpts*numsec,4);
    count=1;
    sec=1;
    for i=1:numpts*numsec
        if(mod(i,numpts)==0)
            quads(count,1)=count;
            quads(count,2)=count-numpts+1;
            quads(count,3)=count+1;
            quads(count,4)=count+numpts;
        else
            quads(count,1)=count;
            quads(count,2)=count+1;
            quads(count,3)=count+numpts+1;
            quads(count,4)=count+numpts;
        end
        count=count+1;
    end
    %Plot quads if desired
    % plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.r')
    % hold on;
    % for i=1:length(quads(:,1))
    %     xpts=nodes([quads(i,:),quads(i,1)],1);
    %     ypts=nodes([quads(i,:),quads(i,1)],2);
    %     zpts=nodes([quads(i,:),quads(i,1)],3);
    %     plot3(xpts,ypts,zpts)
    % end
    % hold off;
    % figure
    tris=zeros(numpts*numsec*2+2*(numpts-2),3);
    for i=1:length(quads(:,1));
        % ordering is important...
        % must be counterclockwise from flow side% ncpts=5;
    % FoilGen=1;
    % DVs=[1,zeros(1,ncpts)+0.025,zeros(1,ncpts)-0.025];
    % [x0,top0,bot0]=GenFoilSpline(DVs,false);
    % 
    % LB=[-5,zeros(1,ncpts),zeros(1,ncpts)-.03];
    % UB=[5,zeros(1,ncpts)+.03,zeros(1,ncpts)];
    % 
    % nvar=2*ncpts+1;
    % ndiv=10;
    % % NEED TO VIOLATE FEWER HIGH-FIDELITY CHOICES...
    % LHS=lhsdesign(ndiv,nvar);
    % points=zeros(ndiv,nvar);
    % f1=zeros(ndiv,1);
    % f2=zeros(ndiv,1);
    % for i=1:ndiv
    %     points(i,:)=LB+(UB-LB).*LHS(i,:);
    %     f1(i)=analyze(points(i,:),fid(1),M,gamma,FoilGen);
    %     f2(i)=analyze(points(i,:),fid(2),M,gamma,FoilGen);
    %     % this maybe should change...
    %     if(f2(i)>=100)
    %         f2(i)=f1(i);
    %     end
    % end
    % fvals=f2-f1;

        tris(2*i-1,:)=quads(i,3:-1:1);
        tris(2*i,:)=quads(i,[4,3,1]);
    end
    % mesh the first side:
    count=2*length(quads(:,1));
    tris(count+1,:)=[(numpts-1)/2,(numpts-1)/2+1,(numpts-1)/2+2];
    %tris(count+2,:)=[1,2,numpts];
    count=count+2;
    for i=1:(numpts-1)/2-1
        tris(count,:)=[i,i+1,numpts-i+1];
        tris(count+1,:)=[numpts-i,numpts-i+1,i+1];
        count=count+2;
    end
    %mesh the back side:
    tris(count,:)=[(numpts-1)/2+2,(numpts-1)/2+1,(numpts-1)/2]+numpts*numsec;
    %tris(count+1,:)=[numpts,2,1]+numpts*numsec;
    count=count+1;
    for i=1:(numpts-1)/2-1
        tris(count,:)=[numpts-i+1,i+1,i]+numpts*numsec;
        tris(count+1,:)=[i+1,numpts-i+1,numpts-i]+numpts*numsec;
    %     tris(count,:)=[numpts-i+1,i+1,i]+numpts*numsec;
    %     tris(count+1,:)=[i,numpts-i+2,numpts-i+1]+numpts*numsec;
        count=count+2;
    end
    % figure
    % plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.r')
    % hold on;
    % % Plot Triangles
    % for i=1:length(tris(:,1))
    %     xpts=nodes([tris(i,:),tris(i,1)],1);
    %     ypts=nodes([tris(i,:),tris(i,1)],2);
    %     zpts=nodes([tris(i,:),tris(i,1)],3);
    %     plot3(xpts,ypts,zpts,'g');
    % end
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')

    % Switch orientation of triangles:
    % tri2=tris;
    % tris(:,1)=tri2(:,3);
    % tris(:,3)=tri2(:,1);
    if(compute_norms)
        % Compute Surface Normals:
        normals=zeros(size(tris));
        centroids=zeros(size(tris));
        tri_nodes=zeros(3,3);
        for i=1:length(tris(:,1))
            tri_nodes(1,:)=nodes(tris(i,1),:);
            tri_nodes(2,:)=nodes(tris(i,2),:);
            tri_nodes(3,:)=nodes(tris(i,3),:);
            vec1=tri_nodes(2,:)-tri_nodes(1,:);
            vec2=tri_nodes(3,:)-tri_nodes(2,:);
            normals(i,:)=cross(vec1,vec2);
            normals(i,:)=normals(i,:)/norm(normals(i,:),2);
            centroids(i,:)=mean(tri_nodes,1);
        end
        h=quiver3(centroids(:,1),centroids(:,2),centroids(:,3),normals(:,1),normals(:,2),normals(:,3),'c');
        hold off;
    end
    % for i=1:length(nodes(:,1))
    %     text(nodes(i,1),nodes(i,2),nodes(i,3),num2str(i));
    % end
    %axis equal;
    fid=fopen([foilname,'.tri'],'w');
    fprintf(fid,'%d %d\n',length(nodes(:,1)),length(tris(:,1)));
    for i=1:length(nodes(:,1))
        fprintf(fid,'%10.6f %10.6f %10.6f\n',nodes(i,1),nodes(i,2),nodes(i,3));
    end
    for i=1:length(tris(:,1))
        fprintf(fid,'%d %d %d\n',tris(i,1),tris(i,2),tris(i,3));
    end
    for i=1:length(tris(:,1))
        fprintf(fid,'%d\n',comp);
    end
    fclose(fid);
    fid=fopen('input.cntl','w');
    fprintf(fid,'#    +--------------------------------------------------------+\n');
    fprintf(fid,'#    |       Steering and Control file for "flowCart"         |\n');
    fprintf(fid,'#    |           3D Cut-Cell Cartesian Flow Solver            |\n');
    fprintf(fid,'#    +--------------------------------------------------------+\n');
    fprintf(fid,'#   \n');
    fprintf(fid,'#      NOTE:  o Start Comments in this file with the "#" character\n');
    fprintf(fid,'#             o Blocks can come in any order\n');
    fprintf(fid,'#             o info within blocks can come in any order\n');
    fprintf(fid,'#\n');
    fprintf(fid,'#  Matlab Automated submission\n');
    fprintf(fid,'#\n');
    fprintf(fid,'$__Case_Information:          # ...Specify Free Stream Quantities\n');
    fprintf(fid,'Mach     %3.2f   #  (double)\n',M);
    fprintf(fid,'alpha    %3.2f   #  (double) - angle of attack\n',alpha);
    fprintf(fid,'beta     0.0   #  (double) - sideslip angle\n');
    fprintf(fid,'\n');
    fprintf(fid,'$__File_Name_Information:\n');
    fprintf(fid,'MeshInfo         Mesh.c3d.Info   # Mesh info file (usually Mesh.c3d.Info)\n');
    fprintf(fid,'MeshFile           Mesh.mg.c3d   # Mesh file\n');
    fprintf(fid,'\n');
    fprintf(fid,'# --NOTE:   surface triangulation specified in "MeshInfo" file ------\n');
    fprintf(fid,'\n');
    fprintf(fid,'$__Solver_Control_Information:\n');
    fprintf(fid,'#   Runge-Kutta Stage Coefficients\n');
    fprintf(fid,'#   stageCoef    GradEval  ->to run 1st order, set GradEval to 0 in all stages\n');
    fprintf(fid,'#    --------    -------\n');
    fprintf(fid,'RK        0.0695      1  #         van Leer 5-stage \n');
    fprintf(fid,'RK        0.1602      0  #         "optimally damped 2nd order scheme"\n');
    fprintf(fid,'RK        0.2898      0  #          AIAA 89-1933-CP (CFLopt = 2.5 1st order)\n');
    fprintf(fid,'RK        0.5060      0  #                          (CFLopt = ~1.2 2nd order)\n');
    fprintf(fid,'RK        1.0         0  #\n');
    fprintf(fid,'                         #                          (CFLopt = 0.694)\n');
    fprintf(fid,'                       # NOTE: GradEval = 0 = no new evaluation at this stage, \n');
    fprintf(fid,'                       #       GradEval = 1 = Yes, re-evaluate at this stage\n');
    fprintf(fid,'CFL           1.2 # CFL number \n');
    fprintf(fid,'Limiter       0   # (int) default is 2, organized in order of increasing \n');
    fprintf(fid,'                  #       dissipation.\n');
    fprintf(fid,'                  #         Limiter Type: 0 = no Limiter\n');
    fprintf(fid,'                  #                       1 = Barth-Jespersen\n');
    fprintf(fid,'                  #                       2 = van Leer\n');
    fprintf(fid,'                  #                       3 = sin limiter\n');
    fprintf(fid,'                  #                       4 = van Albada\n');
    fprintf(fid,'                  #                       5 = MinMod\n');
    fprintf(fid,'                  #\n');
    fprintf(fid,'FluxFun       0   # (int) - Flux Function:   0 = van Leer\n');
    fprintf(fid,'                  #                          1 = Colella 1998\n');
    fprintf(fid,'                  #                          2 = HLLC (alpha test)\n');
    fprintf(fid,'\n');
    fprintf(fid,'Precon        0   # (int) - Preconditioning: 0 = scalar timestep\n');
    fprintf(fid,'wallBCtype    0   # Cut-Cell Boundary Condition type   0 = Agglomerated Normals\n');
    fprintf(fid,'                  #                                    1 = SubCell Resolution\n');
    fprintf(fid,'nMGlev        1   # (int) - Number of Multi-Grid levels  (1 = single grid)\n');
    fprintf(fid,'MG_cycleType  2   # (int) - MultiGrid cycletype: 1 = "V-cycle", 2 = "W-cycle"\n');
    fprintf(fid,'                  # "sawtooth" cycle is: nPre = 1, nPost = 0\n');
    fprintf(fid,'MG_nPre       1   # (int) - no of pre-smoothing  passes in multigrid\n');
    fprintf(fid,'MG_nPost      1   # (int) - no of post-smoothing passes in multigrid\n');
    fprintf(fid,'\n');
    fprintf(fid,'    \n');
    fprintf(fid,'$__Boundary_Conditions: # BC types: 0 = FAR FIELD \n');
    fprintf(fid,'                        #           1 = SYMMETRY\n');
    fprintf(fid,'                        #           2 = INFLOW  (specify all)\n');
    fprintf(fid,'                        #           3 = OUTFLOW (simple extrap)\n');
    fprintf(fid,'Dir_Lo_Hi     0   0 0   # (int) (0/1/2) direction  (int) Low BC   (int) Hi BC\n');
    fprintf(fid,'Dir_Lo_Hi     1   0 0   # (int) (0/1/2) direction  (int) Low BC   (int) Hi BC\n');
    fprintf(fid,'Dir_Lo_Hi     2   1 1   # (int) (0/1/2) direction  (int) Low BC   (int) Hi BC\n');
    fprintf(fid,'\n');
    fprintf(fid,'$__Convergence_History_reporting:\n');
    fprintf(fid,'iForce     1   # (int) - Report residual information every iSkip cycles.\n');
    fprintf(fid,'iHist      1   # (int) - Update "HistoryFile" every iHist cycles.\n');
    fprintf(fid,'nOrders   12   # (int) - Num of orders of Magnitude reduction in residual.\n');
    fprintf(fid,'    \n');
    fprintf(fid,'$__Partition_Information:  \n');
    fprintf(fid,'nPart   1      # (int) - Number of SubDomains to partition into: \n');
    fprintf(fid,'type    1      # (int) - Type of partitioning: 1 = SpaceFillingCurve\n');
    fprintf(fid,'\n');
    fprintf(fid,'$__Post_Processing:\n');
    fprintf(fid,'#                                   Pretty printed cutting planes\n');
    fprintf(fid,'Zslices   0.0001\n');
    fprintf(fid,'# lineSensor LINE1   0.  1.023  0.2331  3.1 1.023  0.2331 \n');
    fprintf(fid,'\n');
    fprintf(fid,'$__Force_Moment_Processing:\n');
    fprintf(fid,'# \n');
    fprintf(fid,'# ... Axis definitions (with respect to body axis directions (Xb,Yb,Zb)\n');
    fprintf(fid,'#                       w/ usual stability and control orientation)\n');
    fprintf(fid,'#  see:\n');
    fprintf(fid,'#  http://people.nas.nasa.gov/~aftosmis/cart3d/clic/html/clic_doc.html#frames\n');
    fprintf(fid,'Model_X_axis  -Xb \n');
    fprintf(fid,'Model_Y_axis  -Zb \n');
    fprintf(fid,'Model_Z_axis  -Yb \n');
    fprintf(fid,'\n');
    fprintf(fid,'# ... reference area and length specifictions\n');
    fprintf(fid,'\n');
    fprintf(fid,'#$_Reference_Area %%f  compNumberList or compNameList\n');
    fprintf(fid,'Reference_Area    1.0  all\n');
    fprintf(fid,'\n');
    fprintf(fid,'#$_Reference_Length %%f compNumberList or compNameList\n');
    fprintf(fid,'Reference_Length  1.  all \n');
    fprintf(fid,'\n');
    fprintf(fid,'# ... Force and Moment Info\n');
    fprintf(fid,'Force        entire\n');
    fprintf(fid,'Moment_Point 0.25 0. 0. entire\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'$__Design_Info:\n');
    fprintf(fid,'\n');
    fprintf(fid,'# Objective Function: SUM of functionals (J)\n');
    fprintf(fid,'# J = 0 -> W(P-T)^N\n');
    fprintf(fid,'# J = 1 -> W(1-P/T)^N\n');
    fprintf(fid,'\n');
    fprintf(fid,'# Ref. Frame = 0 Aerodynamic Frame\n');
    fprintf(fid,'#            = 1 Aircraft (Body) Frame\n');
    fprintf(fid,'\n');
    fprintf(fid,'# L/D -> SIGN(CL)*ABS(CL)^A/(CD+Bias) in Aero Frame\n');
    fprintf(fid,'#     -> SIGN(CN)*ABS(CN)^A/(CA+Bias) in Body Frame\n');
    fprintf(fid,'# Format: \n');
    fprintf(fid,'#\n');
    % *********************Use these lines to refine based on L/D
    % fprintf(fid,'#      Name   Frame   J     N     A     Bias  Target  Weight  Bound  GMP_Comp\n');
    % fprintf(fid,'#    (String) (0,1) (0,1) (int) (dble) (dble) (dble)  (dble)   (0)\n');
    % fprintf(fid,'# ---------------------------------------------------------------------------\n');
    % fprintf(fid,'optLD  LoD      0     0     1     1.     0.     0.     1.       0    entire\n');
    % **********************Use these lines to refine  based on L+D
    fprintf(fid,'# Force Codes: CD=0 Cy=1, CL=2 in Aero Frame\n');
    fprintf(fid,'# Force Codes: CA=0 CY=1 CN=2 in Aircraft Frame\n');
    fprintf(fid,'# Format:\n');
    fprintf(fid,'#      Name    Force     Frame  J      N     Target    Weight    Bound   GMP Comp\n');
    fprintf(fid,'#    (String)  (0,1,2)  (0,1) (0,1,2) (int)   (dble)  (dble)   (-1,0,1)\n');
    fprintf(fid,'#--------------------------------------------------------------------------------\n');
    fprintf(fid,'optForce  CD     0       0    0        1       0.       1.0       0      entire\n');
    fprintf(fid,'optForce  CL     2       0    0        1       0.       0.5      0       entire\n');
    %*********************************
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fclose(fid);
    [r,s]=system('./cart');
    fid=fopen('loadsCC.dat','r');
    if(fid>0)
        for i=1:8
            fgetl(fid);
        end
        str=fgetl(fid);
        Cd_Cart=sscanf(str,'%*s %*s %*s %*s %f');
        str=fgetl(fid);
        sideF=sscanf(str,'%*s %*s %*s %*s %f');
        str=fgetl(fid);
        Cl_Cart=sscanf(str,'%*s %*s %*s %*s %f');

        fclose(fid);
        if(isnan(Cd_Cart) || isnan(Cl_Cart))
            % cart failed
            Cl_Cart=0;
            Cd_Cart=-1;
            disp('Cart Case Failed')
        else
            disp('Cart Case Worked!')
        end
    else
        % Cart3D Failed
        Cl_Cart=0;
        Cd_Cart=-1;
        disp('Cart Case Failed')
    end
end
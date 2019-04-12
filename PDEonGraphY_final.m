% PDE model on a Y shaped graph  
%                 _________b2 
%                /  e2 
%      e1     a2/          
% a1--------b1               
%             a3\
%                \_________   
%                   e3     b3 

function PDEonGraphY_final( ) 
global dx p Tstep int_dx x N1 

    N = 100; %default('Number of grid points (Default: 100)',100); 
    T = 15; %default('Final time (Default: 4)',4); 
    
%% Space discretization, FDM 
    L1min = 0; L1max = 1; 
    N1 = N+1; 
    x   = linspace(L1min, L1max, N1(1))';  dx = (L1max-L1min)/N(1); 
    int_dx = ones(N(1)+1,1)*dx; int_dx(1) = int_dx(1)/2; int_dx(end) = int_dx(end)/2; 
    FDM_1D( N(1) ); 
    
%% IC 
    U = IC;  
    
%% Parameter 
    cD  = default('Diffusion coefficient D (Default: 0.01)', 0.01);
    cW  = [1, 0.5, 0.5]; % 'vertical fluctuation on each edges' 
    cA  = default('Advection coefficient C (Default: 0.2)', 0.2); 
    p   = default('asymmetry parameter p (Default: 0.5)', 0.5); 
    
    nBC = 2; %default('Boundary Condition 1/2/3 (Default: 2)', 2);
    set_parameters( cD, cA, p ); 
    
%% Time integration     
    dt=0.001; %min( dx/max(max(abs(U)))/2, dx^2/max(cD)/2 );   % Time step     
    t = 0;     Tstep = T/5;      % parameters for time dependent forcing/ plot 
    %% plot initial     
    nfig = floor( rand(1)*1000 ); disp( strcat( '......Time marching.... Results in Figure ', int2str(nfig) ) ); 
    plot_graph3d( U, nfig, 1 ); 
    plot_graph2d( U, nfig+1 );       
    
    Utot = zeros( T/Tstep, 1 ); 
    for n = 1:3; Utot(1,n) = U(:,n)' * int_dx; end 
    
    Tstart = tic;
    for nt = 1:floor( (T+eps)/Tstep )
      for time = 1:(Tstep/dt) 
        
        %% RK3 TVD 
        RKdu = Compute_du( U ); 
        u = U + RKdu * dt;
        u = BC( u, t, cD, cA, cW );   
        
        RKdu = Compute_du( u ); 
        u = 0.75*U + 0.25*u + RKdu*0.25*dt; 
        u = BC( u, t+0.5*dt, cD, cA, cW );   

        RKdu = Compute_du( u ); 
        U = U/3 + u*2/3 + RKdu*2/3*dt; 
        U = BC( U, t+dt, cD, cA, cW ); 
        
        t = t + dt; 
        
        % postprecess positivity 
        % U( U < 0 ) = 0 ;  
        
      end
      
      %  dt=min( dt, min( dx/max(max(abs(U)))/2, dx^2/max(cD)/2 ) );
      if( min( dx/max(max(abs(U)))/2, dx^2/max(cD)/2 ) < dt ) 
          disp('reduce time step' ); 
      end 
      %%%%% print results 
          plot_graph3d( U, nfig, nt+1 ); 
          plot_graph2d( U, nfig+1); 
      %%%%% Compute Tot cell 
      for n = 1:3; Utot(nt+1,n) = U(:,n)' * int_dx; end 
      
    end
    toc(Tstart) 

    clear Dx Dxx x dx cDiff cAdv cAdvM N1 p int_dx  
    
end 

function du = Compute_du( U, time ) 
global Dxx Dx cDiff cAdvM 

    %% PDF on graph 
    du = Dxx * U; 
    for n = 1:3 
        du(:,n) = cDiff(n) * du(:,n); 
    end 
    
    du = du - Dx * (cAdvM.* U); 
    
end

function set_parameters( cD, cA, p ) 
global cDiff cAdvM x 

  cDiff = cD*ones(3, 1 ); 

  cAdvM = ones( length(x), 1 )* cA;  
  cAdvM(:,2) = (1-x.^2) * cA *(1-p);   
  cAdvM(:,3) = (1-x.^2) * cA *p; 

end 

function  U = BC( U, t, cD, cA, cW ) 
global p dt 
%% 4th order BC 
%%%% Neunman zero condition at b2, b3 
    U(end,2) = (48*U(end-1,2)-36*U(end-2,2)+16*U(end-3,2)-3*U(end-4,2))/25; 
    U(end,3) = (48*U(end-1,3)-36*U(end-2,3)+16*U(end-3,3)-3*U(end-4,3))/25; 
    
%%%% Dirichlet boundry condition at a1 
    U( 1, 1) = Normal( 0, 0.2, -cA*t )* sqrt( 2.0*pi ) * 0.2; 

% gluing and continuity condition 
    U(end,1) = ( ( 48*U(2,2)-36*U(3,2)+16*U(4,2)-3*U(5,2) )*cW(2) ... 
                +( 48*U(2,3)-36*U(3,3)+16*U(4,3)-3*U(5,3) )*cW(3) ... 
                -(-48*U(end-1,1)+36*U(end-2,1)-16*U(end-3,1)+3*U(end-4,1)) )/(25*sum(cW)); 
    U(1,2) = U(end,1); U(1,3) = U(end,1); 
    
end 

function U = IC  
global x N1 

    U = zeros( N1, 3 ); 
    U(:,1) = Normal( -1, 0.2, -x(end:-1:1) )* sqrt( 2.0*pi ) * 0.2; 
    U(:,2) = Normal( -1, 0.2, x )* sqrt( 2.0*pi ) * 0.2; 
    U(:,3) = U(:,2); 
    
end 

function plot_graph3d( U, nfig, nfigsub )
global x N1 Tstep 
    if( nargin == 1 );         figure;        nfigsub = 1; 
    elseif( nargin == 2);      figure(nfig);  nfigsub = 1; 
    else ;      figure(nfig); 
    end 
    hold on; 
    subplot( 1, 6, nfigsub ); hold on; 
    plot3( -x(end:-1:1), zeros(N1,1), U(:,1), 'black-', 'linewidth', 2 ); 
    plot3( x, sqrt(x), U(:,2), 'b--', 'linewidth', 2  ); 
    plot3( x, -sqrt(x), U(:,3), 'r-.', 'linewidth', 2  ); hold off; 
    title( strcat('t=',num2str((nfigsub-1)*Tstep)), 'Fontsize', 14);
    view( [-17 6] ); box on; grid on; hold off; 
    set(findall(gcf,'-property','FontSize'),'FontSize',14);    
    
end 

function plot_graph2d( U, nfig )
global x  
    if( nargin == 1 );          figure;   
    elseif( nargin == 2 );      figure(nfig); 
    end 
    hold on; plot( -x(end:-1:1), U(:,1), 'black-', 'linewidth', 1 ); 
    plot( x,  U(:,2), 'b--', 'linewidth', 1  ); 
    plot( x,  U(:,3), 'r-.', 'linewidth', 1  ); 
    xlabel('$x$','Fontsize', 14,'Interpreter','Latex');
    ylabel('$y$','Fontsize', 14,'Interpreter','Latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);    
    
end 


function PDFonGraph12_StringRay_final(  )
global x dx dt N1 int_dx vertex edge BC nEdge printon T 
global cAdv cDiff 
global AMLTime AMLOn 

    N = 100; %default('Number of grid points (Default: 100)',100); 
    T = default('Final time (Default: 6)',6); 
    dt = 1*0.1^3; 
        
    printon = 1; 
    
    AMLOn = default('AML condition: Yes=1 / No=0 (Default: 0)', 0);
    if( AMLOn )         
        AMLTime = default('AML starts at: (Default: 4)', 4);
    else; AMLOn = 0; 
    end 
    
    
%% FDM 
    L1min = 0; L1max = 1; 
    N1 = N+1; 
    x   = linspace(L1min, L1max, N1(1))';  dx = (L1max-L1min)/N(1); 
    int_dx = ones(N1(1),1)*dx; int_dx(1) = int_dx(1)/2; int_dx(end) = int_dx(end)/2; 
    FDM_1D( N(1) ); 
    
%% IC & parameters 
   [vertex, edge, BC, nEdge] = edge_graph( );  

    U = IC( ); 
    
    nfig = floor( rand(1)*1000 ); disp( strcat( '.... Results in Figure : ', int2str(nfig), '....' ) ); 
    
    if( printon ) 
      plot_graph3d( U, nfig, [2,5,1] );   title('t=0')
      %%%% print log scale 
%       tmp = log10(U); tmp( tmp <= -10 ) = -10; 
%       plot_graph3d( tmp, 10*nfig, [2,5,1] );  %caxis( [-2 1] )      
      plot_graph2d( U, nfig+1, [2,5,1], 1  );   title('t=0') 
    end 
    
    nTotCell(1) = int_dx'*U*ones(nEdge,1); 
    nTotEdge(:,1) = int_dx'*U; 
    Usave_h(:,1) = plot_graph1d_half( U );   
    
    
%% Time integration 
    t = 0;  Tstep = 0.1;      Tstart = tic; 
    if( dx/max(cAdv)/2 < dt || dx^2/max(cDiff)/2 < dt ) 
        disp( 'reduce dt' ); 
    end 
    nBC = 0; 
    
    for nt = 1:(T/Tstep) 
      for time = 1:(Tstep/dt)  
        %% RK3 TVD 
        RKdu = Compute_du( U, t ); 
        u = U + RKdu * dt;
        u = Compute_BC( u, nBC+time ); %t ); % 
        
        RKdu = Compute_du( u, t+0.5*dt ); 
        u = 0.75*U + 0.25*u + RKdu*0.25*dt; 
        u = Compute_BC( u, nBC+time );  %t+0.5*dt ); %
        
        RKdu = Compute_du( u, t+dt ); 
        U = U/3 + u*2/3 + RKdu*2/3*dt; 
        U = Compute_BC( U, nBC+time ); %t+dt ); % 
        
        U( U < 0 ) = 0 ; 
        t = t + dt; 
        
      end 
      
      nTotCell(nt+1) = int_dx'*U*ones(nEdge,1);
      nTotEdge(:,nt+1) = int_dx'*U; 
      Usave_h(:,nt+1) = plot_graph1d_half( U );          
      
      if( mod(Tstep*nt, 1 ) == 0 )
        
        if( printon ) 
        plot_graph3d( U, nfig, [2,5, min(Tstep*nt+1,10)] );   title(strcat('t=',num2str(Tstep*nt))); if(max(U(:))>10^6); caxis( [0 10^6] ); end 
%         tmp = log10(U); tmp( tmp <= -10 ) = -10; 
%         plot_graph3d( tmp, 10*nfig, [2,5, min(Tstep*nt+1,10)] ); %caxis( [-2 1] )
        plot_graph2d( U, nfig+1, [2,5, min(Tstep*nt+1,10)], 1 );   title(strcat('t=',num2str(Tstep*nt)));  
        end 
      end 
      
      if(  AMLOn~=0 && Tstep*nt == AMLTime ) 
          UpdateC( ); 
          
      end; 
      
      nBC = min( length(BC.Uinit)-time, floor(nt*(Tstep/(dt-eps))) ); 
      
    end 
    toc(Tstart) 

    if( printon ) 
      figure(nfig+3); hold on; plot( [0:Tstep:T], nTotCell/nTotCell(1)    );    
      xlabel('t'); ylabel('Total number of cells'); set(gca,'yscale','log' ); if(max(U(:))>10^6); ylim( [1 10^6] ); else; ylim( [1 10^4] ); end 
      figure(nfig+2); plot( [0:Tstep:T], Usave_h'/nTotCell(1) ); 
      xlabel('t'); ylabel('Cell number in each cluster'); set(gca,'yscale','log' ); if(max(U(:))>10^6); ylim( [.01 10^6] ); else; ylim( [.01 10^4] ); end
    end
    
 
end 

function du = Compute_du( U, time ) 
global N1 Dxx Dx cDiffM cAdvM cRctM 

    %% PDF on graph 
    du = (Dxx * U)*cDiffM;     
    du = du - Dx * (cAdvM.* U); 
    du = du + cRctM .* U; 


end

function  U = Compute_BC( U, time ) 
global BC edge cDiff 

    for nn = BC.noout 
      ind =  BC.indin(nn,1:BC.Nindin(nn));   %find( edge(:,2) == nn ); 
        tmp = sum( (48*U(end-1,ind)-36*U(end-2,ind)+16*U(end-3,ind)-3*U(end-4,ind)) .*BC.ratein(ind)' , 2 ); 
        tmp = tmp / 25; %/sum(BC.ratein(ind)); 
        U(end,ind) = tmp;       
    end 
    for nn = setdiff( 9:12, BC.noout )   
      ind =  BC.indin(nn,1:BC.Nindin(nn));   %find( edge(:,2) == nn ); 
        tmp = sum( (48*U(end-1,ind)-36*U(end-2,ind)+16*U(end-3,ind)-3*U(end-4,ind)) .*BC.rateinprev(ind)' , 2 ); 
        tmp = tmp / 25;
        U(end,ind) = tmp;       
    end     
    nn = 10; ind =  BC.indin(nn,1); 
    U(1,49:51) = U(end,ind); 
    
    for nn = BC.noin 
      ind =  BC.indout(nn,1:BC.Nindout(nn));        
      tmp = BC.Uinit(nn,time); 
      U(1,ind) = tmp; 
    end 
    for nn = 1:3 
      U(end,nn) = BC.Uinit(edge(nn,2),time); 
    end 
    
    for nn = BC.inNout 
      indout = BC.indout(nn,1:BC.Nindout(nn));        
      indin  = BC.indin( nn,1:BC.Nindin(nn)); 
      
      % Diff 
      tmp = sum( (48*U(2,indout)-36*U(3,indout)+16*U(4,indout)-3*U(5,indout)) .*cDiff(indout)' , 2 ) ;  
      tmp = tmp - sum( (-48*U(end-1,indin)+36*U(end-2,indin)-16*U(end-3,indin)+3*U(end-4,indin)) .*cDiff(indin)', 2 ) ; 
      tmp = tmp / 25 / (sum(cDiff(indout))+sum(cDiff(indin))) ; 
            
      U(1,indout) = tmp; U(end,indin) = tmp; 
    end 
    
end 

function U = IC( ) 
global x N1 edge nEdge BC int_dx 
    nIC = [24    66   155 0 0 0 0 0 0 0 0 0]; %/ 245; 
    pp  = ones(12,1) * 0.1; 
    
    U = zeros( N1, nEdge ); 
    for n = 1:3 
        ind = find( edge(:,1) == n );
        for m = 1:length(ind) 
            sdtmp = pp(n); 
            Uinit1 = Normal( 0, sdtmp, x )*(sqrt(2*pi)*sdtmp); 
            sdtmp = pp(edge(ind(m),2)); 
            Uinit2 = Normal( 1, sdtmp, x )*(sqrt(2*pi)*sdtmp); 
            coeff = [Uinit1(1),Uinit1(end); Uinit2(1),Uinit2(end)]\ [nIC(edge(ind(m),1));nIC(edge(ind(m),2))]; 
            if( sum(coeff<0) ); 
              coeff(coeff<0) = 0; 
            end 
            U(:, ind(m)) = coeff(1)*Uinit1 + coeff(2)*Uinit2;           
%           end 
        end 
    end
    tmp = plot_graph1d_half( U );
    U = U / sum(tmp) * 245; 

    U( U < eps ) = 0; 
    
end 

function plot_graph2d( U, nfig, nfigsub, ncase ) 
global x N1 vertex edge colororder 
    if( nargin == 1 );          figure; 
    elseif( nargin == 2 );      figure(nfig); 
    else; figure(nfig);  subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); 
    end 
    if( nargin ~= 4);  ncase = 1; end 
    
    subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); 
    if( ncase == 1 ) 
    tmp = []; for n = 1:12; ind = find(edge(:,1)==n); 
        ntmp(n) = length(ind); 
        tmp = [tmp; ind ]; end; 
    imagesc( U(:,tmp(end:-1:1))' ); 
 
    tt = cumsum( ntmp ); [tt,ind] = unique( tt ); 
    set(gca,'YTick', tt, 'YTicklabel', ind(end:-1:1) ); 
    set(gca,'XTick', [1:(N1-1)/5:N1], 'XTicklabel', [0:.2:1] )
    
    elseif( ncase == 2 ) 
    subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); 
    tmp = []; for n = 1:12; ind = find(edge(:,2)==n);
        ntmp(n) = length(ind);         
        tmp = [tmp; ind ]; end;
    imagesc( U(:,tmp(end:-1:1))' ); 
 
    tt = cumsum( ntmp ); [tt,ind] = unique( tt ); 
    set(gca,'YTick', tt, 'YTicklabel', ind(end:-1:1) ); set( gca, 'YAxisLocation','right');
    set(gca,'XTick', [1:(N1-1)/5:N1], 'XTicklabel', [0:.2:1] )
    end 
end


function plot_graph3d( U, nfig, nfigsub ) 
global x N1 vertex edge colororder 
    if( nargin == 1 );          figure; 
    elseif( nargin == 2 );      figure(nfig); 
    else; figure(nfig);  subplot( nfigsub(1), nfigsub(2), nfigsub(3) ); 
    end 
    
    load( 'edgeline.mat' )  
    hold on;  axis( [-0.6 0.6 -0.6 0.7] ); 
    set(gca, 'color', 'black' ); 

    for m = 1:size(edge, 1)  
      surface([edgex{m};edgex{m}],[edgey{m};edgey{m}],[U(:,m)';U(:,m)'],[U(:,m)';U(:,m)'],'facecol','no','edgecol','interp','linew',1);  
    end 
    
%     load( 'colororder.mat' ) 

    xx =[    0.4600    0.4000    0.3800    0.0216    0.0700    0.0216 ... 
       -0.0566   -0.0566   -0.0213    -0.3059  0.3023   -0.4697  ];
    yy = [         0   -0.0700    0.0500   -0.0666   -0.0000    0.0666 ... 
        0.0411   -0.0411   -0.5359    -0.2311    0.5301  0.0104  ];
    for n = 0:2; xx(n+1) = 0.07*cos(2*pi/3*n) + 0.4; yy(n+1) = 0.07*sin(2*pi/3*n); end 
    for n = 0:4; xx(n+4) = 0.1*cos(2*pi/5*n); yy(n+4) = 0.1*sin(2*pi/5*n); end 


    plot3( xx, yy, -10*ones(12,1), 'o', 'color', [.5 .5 .5], 'markerfacecolor', [.5 .5 .5] ) 
    set(gca,'XTick', [], 'XTicklabel', [] )
    set(gca,'YTick', [], 'YTicklabel', [] )

    view( [0 90] );  box on;  %grid on; 

    
end 


function UU = plot_graph1d_half( U, nfig ) 
global x N1 edge nEdge int_dx 
    Nf = (N1+1)/2; 
    
    UU = zeros( 12, 1 ); 
    for n = 1:nEdge 
        UU(edge(n,1)) = UU(edge(n,1)) + sum(U(1:Nf, n).*int_dx(1:Nf)) ; 
        UU(edge(n,2)) = UU(edge(n,2)) + sum(U(Nf:end, n).*int_dx(Nf:end)) ; 
    end    
    
    if( nargin > 1 )
      figure(nfig); hold on;  plot( UU, '-x', 'linewidth', 1 )
    end     
end 

function [vertex, edge, BC, nEdge] = edge_graph( Cout, Cin ) 
global cAdv cDiff dt cDiffM cAdvM cRct cRctM N1 x T  
%     plotting index 
    vertex = [    3.8776e-01   2.6652e-01  -1.6446e-01
   3.4409e-01   2.4803e-01  -1.5479e-01
   3.0300e-01   1.6100e-01  -1.5800e-01
   3.3590e-02  -5.2962e-02  -9.3491e-02
   1.1397e-01  -3.6054e-02  -1.1687e-01
   1.2959e-01   1.1360e-01  -1.0723e-01
  -5.5520e-02   6.2699e-02  -2.0664e-04
  -9.4601e-02  -1.6594e-01  -6.6782e-02
  -1.9790e-02  -2.4617e-01  -1.0378e-01
  -2.3224e-01   7.4680e-02   3.4614e-02
   1.0773e-01   1.9809e-02   4.3378e-01
  -4.5300e-01   7.1000e-02   5.5000e-02]; 
   
% mean distance  1.8123e-02,  4.2047e-02 
    edge = [1 2; 2 3; 3 1]; % [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; 
    ntmp = 3; 
    for n = 1:3 
        edge((ntmp+1):(ntmp+5),1) = n;  edge((ntmp+1):(ntmp+5),2) = 4:8; 
        ntmp = ntmp + 5; 
    end 
    edge = [edge; 4 5; 5 6; 6 7; 7 8; 8 4; 4 6; 5 7; 6 8; 7 4; 8 5]; 
    ntmp = 28; 
    for n = 4:8 
        edge((ntmp+1):(ntmp+4),1) = n;  edge((ntmp+1):(ntmp+4),2) = 9:12; 
        ntmp = ntmp + 4; 
    end 
    edge = [edge; 10 9; 10 11; 10 12]; 
    nEdge = size( edge, 1 ); 
    
    cDiff = 0.1^2.5 * ones( nEdge, 1 ); 
    cDiff( [1:3, 19:28, 50:51 ] )  = 0.1^1.5; 
    
    cDiffM = diag( cDiff ); cDiffM = sparse( cDiffM ); 
    
    nIC = [24    66   155]/ sum([24    66   155]);   
    p1  = [236   36    27   11   60]/sum([236   36    27   11   60]);  
    p2  = [192  223   227   54]/sum([192   223    227 54]);  
    
%   [1-3]-[4-8] [4-8]-9  [4-8]-10   [4-8]-11   [4-8]-12 
% from.  for n = 1:12; for nn = 1:12; edgeLen(n,nn) = norm( vertex(n,:)- vertex(nn,:) ); end; end

    BC.cind = [4:18, 29:48]; 

    if( nargin == 0 ) 
        Cout = zeros( 12 ); 
        for n = 1:3;  Cout( n, 4:8  ) = p1;     end 
        for n = 4:8;  Cout( n, 9:12 ) = p2;     end 
        Cin = zeros( 12 ); 
        for n = 4:8;   Cin( 1:3, n ) = nIC; end 
        for n = 9:12;  Cin( 4:8, n ) = p1;  end 
    end 
    
    for n = 1:nEdge 
        BC.rateout(n,1) = Cout( edge(n,1),edge(n,2) ); 
        BC.ratein( n,1) = Cin( edge(n,1),edge(n,2) ); 
    end     
    
    BC.Nindout = zeros( 1, 12 ); BC.Nindin = zeros( 1, 12 ); 
    BC.inNout = [4:8]; 
    for nn = BC.inNout  
        ind = find( edge(:,1) == nn ); % ind = setdiff( ind, 19:28 ); 
        BC.Nindout(nn) = length(ind); 
        BC.indout(nn,1:BC.Nindout(nn)) =  ind; 
        
        ind = find( edge(:,2) == nn ); % ind = setdiff( ind, 19:28 );  
        BC.Nindin(nn) = length(ind); 
        BC.indin(nn,1:BC.Nindin(nn)) =  ind; 
    end    
    BC.noin = [1:3]; 
    for nn = BC.noin 
        ind = find( edge(:,1) == nn ); % ind = setdiff( ind, 1:3 ); 
        BC.Nindout(nn) = length(ind); 
        BC.indout(nn,1:BC.Nindout(nn)) =  ind; 
    end 
    BC.noout = [9:12]; 
    for nn = BC.noout     
        ind = find( edge(:,2) == nn ); 
        BC.Nindin(nn) = length(ind); 
        BC.indin(nn,1:BC.Nindin(nn)) =  ind; 
    end 

    load( '180310_PDFonGraph12_coeff.mat' ); 
%     global cAdv0 cRct0 
    cAdv = cAdv0; 
    
    cAdvM = zeros( N1, nEdge ); 
    for n = 1:nEdge  % [BC.inNout,BC.noin]  
        cAdvM(:,n) = linspace( cAdv(edge(n,1))* BC.rateout(n), cAdv(edge(n,2))*BC.ratein(n), N1 ); 
    end 
    for nn = BC.noout 
      ind = find( edge(:,2) == nn ); 
      for m = 1:length(ind) 
        n = ind(m); 
        cAdvM(:,n) = (1-x.^2) .* linspace( cAdv(edge(n,1))* BC.rateout(n), cAdv(edge(n,2))*BC.ratein(n), N1 )'; 
      end 
    end 
    cAdvM = sparse(cAdvM); 
    
    cRct = cRct0*ones( 12, 1 ); cRct(1:3) = cRct(1:3)/100; 
    cRctM = zeros( N1, nEdge ); 
    for n = 1:nEdge 
        cRctM(:,n) = linspace( cRct(edge(n,1)) , cRct(edge(n,2)) , N1 ); 
    end 
    
    BC.Diclt(1) = 3;     BC.Diclt(2) = 1;     BC.Diclt(3) = 2; 
    
    BC.Uinit(1,:) = 24 * Normal( 0, 0.2, cAdv(1)*[0:dt:2] )*(sqrt(2*pi)*0.2); 
    BC.Uinit(2,:) = 66 * Normal( 0, 0.2, cAdv(2)*[0:dt:2] )*(sqrt(2*pi)*0.2); 
    BC.Uinit(3,:) = 155 * Normal( 0, 0.2, cAdv(3)*[0:dt:2] )*(sqrt(2*pi)*0.2); 
    
end 

function UpdateC( ) 
global vertex edge BC nEdge cRctM N1 cDiff cAdvM cDiffM x AMLTime 
 
    rateinprev = BC.ratein; 
    
    nIC = [24    66   155]/ 245; 
    p1 =  [236    36    27    11    60 ]/370; 
    p2 = [192   223    0 54]/sum([192   223    0 54]); 
    
    Cout = zeros( 12 ); 
    for n = 1:3;   Cout( n, 4:8  ) = p1;     end 
    for n = 4:8;   Cout( n, 9:12 ) = p2;     end 
    Cin = zeros( 12 ); 
    for n = 4:8;   Cin( 1:3, n ) = nIC; end 
    for n = [9:10,12];  Cin( 4:8, n ) = p1;  end 
    
    for n = 1:nEdge 
        BC.rateout(n,1) = Cout( edge(n,1),edge(n,2) ); 
        BC.ratein( n,1) = Cin( edge(n,1),edge(n,2) ); 
    end 
    
    ind = find( edge(:,2) == 11 ); 
    for m = 1:length(ind) 
      cAdvM(:,ind(m)) = 0; 
    end 
    cAdvM = sparse(cAdvM);     
    
    cRctM = zeros( N1, nEdge ); 
    cRct = 0.6482 * ones( 1, 12 ); 
    cRct(11) = cRct(11)*10; 
    cRct([3,5]) = cRct([3,5])*2.5; 
    
    Nf = (N1+1)/2; 
    for n = 1:size( edge, 1 ) 
        cRctM(:,n) =  linspace( cRct(edge(n,1)), cRct(edge(n,2)), N1 ); 
    end 
    
    BC.rateinprev = rateinprev; 
    BC.noout = setdiff( BC.noout, 11 ); 
    
    ind = find( edge(:,2) == 11 ); 
%     cDiff( ind ) = cDiff( ind ) / 10 ; 
%     cDiffM = diag( cDiff ); cDiffM = sparse( cDiffM ); 
    
    
    disp( strcat( 'AML condition at t=', int2str(AMLTime) ) ) 

end 

function sensitivity_system(ncase, np) 
global cAdv cDiff cDiffM cAdvM cRct cRctM 

    if( ncase == 1 ) 
        cDiffM = cDiffM*(1+np); 
    elseif( ncase == 2 ) 
        cAdvM = cAdvM*(1+np);    
    elseif( ncase == 3 ) 
        cRctM = cRctM*(1+np);        
    end        

end 



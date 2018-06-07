
function Secant4cAdv()
global cTarget printon cAdv0 cRct0 

cTarget = [   24   66  155    0    0    0    0    0    0    0    0   0; 
               0    0    0  236   36   27   11   60    0    0    0   0; 
               0    0    1    0    0    0    0    0  192  223  227  54]; 
cTargetsum = sum( cTarget(2,4:8) ); 
cRct0 =  0.8096; 

itrMax = 10; 
cAdv  = 1*ones( 12, itrMax ); 
Usave = zeros( 12, itrMax ); 
fsave = zeros( 12, itrMax ); 

% cAdv(4:8,1) = 0.5;   cAdv(4:8,2) = 1; 
% 
% for nitr = 1:2 
%     cAdv0 = cAdv(:,nitr); 
%    [~, Usave(:,nitr)] = PDFonGraph12_StringRay_final( 2, 0, [], 90 );  
%    [fsave(:,nitr), error( nitr )] = myfunc( Usave(:,nitr), 2 ); 
%    
%     cRct0 = cRct0 + log( cTargetsum / sum( Usave(4:8,nitr) ) );    
% end 
% 
% while( error( nitr )  > 0.1^2 && nitr <= itrMax ) 
%     nitr = nitr + 1; 
%     
%     tmp = secantUpdate( cAdv(:,(nitr-2):(nitr-1)), fsave(:,(nitr-2):(nitr-1)) ); 
%     cAdv(4:8,nitr) = tmp(4:8); 
%     cAdv0 = cAdv(:,nitr); 
%     
%    [~, Usave(:,nitr) ] = PDFonGraph12_StringRay_final( 2, 0, [], 90 );      
%     cRct0 = cRct0 + log( cTargetsum / sum( Usave(4:8,nitr) ) ); 
%     cRct(nitr) = cRct0; 
%     
%    [fsave(:,nitr), error( nitr ) ] = myfunc( Usave(:,nitr), 2 ); 
%     disp( error(nitr) );  
% end 
% 
% [~,ind] = min( error ); 
%  cAdv0 = cAdv(:, ind);  cRct0 = cRct(ind); 
%  
% [U, Usave, nTotCell] = PDFonGraph12_StringRay_final( 2, 0, [], 90 ); 
% save( '180310_UatT1.mat', 'U', 'cAdv0', 'cRct0' ); 

load( '180310_UatT1.mat' ) 
cTargetsum = sum( cTarget(3,9:12) ); 
cRct0 =  0.4;  

cAdv = cAdv0*ones(1,itrMax); cAdv(4:end,1) = cAdv(4:end,1)+0.4;  cAdv(4:end,2) = cAdv(4:end,2)+0.5;  cAdv(4:end,3) = cAdv(4:end,3)+2; 
for nitr = 1:2  
    cAdv0 = cAdv(:,nitr); 
   [~, Usave(:,nitr)] = PDFonGraph12_StringRay_final( 2, 0, U, 90 ); 
   [fsave(:,nitr), error( nitr )] = myfunc( Usave(:,nitr), 3 ); 
    figure(100414); hold on; plot( Usave(:,nitr) )
    disp( error(nitr) ); 
    
end 

[~,ind] = min( error(1:2) ); cAdv(:,3:end) = cAdv(:,ind)*ones(1,itrMax-2); 

while( error( nitr )  > 0.1^2 && nitr <= itrMax )     
    nitr = nitr + 1; 
    
    tmp = secantUpdate( cAdv(:,(nitr-2):(nitr-1)), fsave(:,(nitr-2):(nitr-1)) ); 
    cAdv(9:12,nitr) = tmp(9:12); 
    cAdv0 = cAdv(:,nitr); 
   [~, Usave(:,nitr)] = PDFonGraph12_StringRay_final( 2, 0, U, 90 );     
   [fsave(:,nitr), error( nitr )] = myfunc( Usave(:,nitr), 3 ); 
    figure(100414); hold on; plot( Usave(:,nitr) )
    disp( error(nitr) ); 
end 
 
[~,ind] = min( error ); 
 cAdv0 = cAdv(:, ind); 
 
 
cAdv1 = cAdv0; cRct1 = cRct0;  
load( '180310_UatT1.mat' ) 
save( '180310_PDFonGraph12_coeff.mat', 'cAdv0', 'cRct0', 'cAdv1', 'cRct1' ); 

end 

function cAdvUpdate = secantUpdate( cAdv, fsave, ncase ) 
    cAdvUpdate = cAdv(:,2) - fsave(:,2)./ (fsave(:,2)-fsave(:,1)) .* (cAdv(:,2)-cAdv(:,1) ); 
%     for n = union( find(cAdvUpdate<0), find(cAdvUpdate>5)  ) 
%         cAdvUpdate(n) = cAdv(n,2);           
%     end 
    for n = find(cAdvUpdate<0)
        cAdvUpdate(n) = 0; % cAdv(n,2);           
    end 
    for n = find(cAdvUpdate>5) 
        cAdvUpdate(n) = 5;           
    end     
    for n = find( fsave(:,2) == fsave(:,1) ) 
      if( abs(fsave(n,2)) < 0.1^2 ) 
          cAdvUpdate(n) = cAdv(n,2); 
      else 
          cAdvUpdate(n) = cAdv(n,2)-0.1; 
      end 
    end 
end 

function [fval, error, error_] = myfunc( Usave, ncase )
global cTarget 

fval = zeros( 12, 1 ); 
if( ncase == 2 ) 
    fval(4:8) = Usave(4:8)/sum(Usave(4:8)) - cTarget(ncase,4:8)'/sum(cTarget(ncase,4:8)); 
%     fval(4:8) = Usave(4:8) - cTarget(ncase,4:8)'; 
%     fval(1:8) = Usave(1:8)/sum(Usave(1:8)) - cTarget(ncase,1:8)'/sum(cTarget(ncase,1:8)); 
    error = norm( fval ); 
    error_ = max( abs(fval) ); 
elseif( ncase == 3 ) 
%     fval(4:12) = Usave(4:12)/sum(Usave(4:12)) - cTarget(ncase,4:12)'/sum(cTarget(ncase,4:12)); 
    fval(9:12) = Usave(9:12) - cTarget(ncase,9:12)'; 
    error = norm( fval(9:12) )/norm(cTarget(ncase,9:12)); 
    error_ = max( abs(fval(9:12)) ); 
end     
end 

    


%%**************************************************************************
%% Compute pair-wise distances within R limit between known-unknown and
%% unknown-unknown pairs.
%%
%% Dall = randistance(P0,PP,Radius,nf,noisetype);
%%
%% P0 : anchor positions
%%      nfix being the number of anchors
%% PP : sensor positions 
%%      npts being the number of unknown sensors
%% nf : noise factor 
%% noisetype : 1 - Normal :  dnoisy = dactual + N(0,1)*nf (Default)
%%             2 - Multiplicative Normal : dnoisy= dactual*(1 + N(0,1)*nf)
%%             3 - Log Normal : dnoisy= dactual * 10^(N(0,1)*nf)
%% Dall = [DD, D0]
%% DD = [sensor-sensor distance];
%% D0 = [anchor-sensor distance];
%%*************************************************************************
function [Dall] = randistance(P0,PP,Radius)
  [dim,npts] = size(PP);
  nfix = size(P0,2);
  D0 = sparse(npts,nfix);
  DD = sparse(npts,npts); 
%%
  for j = 1:npts    
     if (nfix > 0)
        tmp = PP(:,j)*ones(1,nfix)-P0;
        rr = sqrt(sum(tmp.*tmp));      
        idx = find(rr < Radius);
        rr = rr(idx); 
        if (~isempty(idx))
            D0(j,idx) = rr';
        end
     end        
     if (j > 1)
        tmp = PP(:,j)*ones(1,j-1) - PP(:,1:j-1);
        rr = sqrt(sum(tmp.*tmp));
        idx = find(rr < Radius);
        rr = rr(idx); 
        if (~isempty(idx))
            DD(idx,j) = rr';
        end
     end
  end
  DD = triu(DD,1) + triu(DD,1)';
  Dall = [DD, D0];
%%*************************************************************************

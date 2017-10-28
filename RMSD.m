%% Calculate RMSD
%%
function rmsd = RMSD(X0,PP)
    sum=0;
    for i=1:100
    add=(norm(X0(:,i)-PP(:,i)))^2;
    sum=sum+add;
    end
    rmsd=1/sqrt(100)*sqrt(sum);
end
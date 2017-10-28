%% Calculate Euclidean Distance
%%
function ans = distance(x,y)
ans = 0;
for i = 1:length(x)
    ans = ans + (x(i)-y(i))*(x(i)-y(i));
end
end
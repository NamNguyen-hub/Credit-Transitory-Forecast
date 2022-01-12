function y = transform_data(data,transform_code)
%Benjamin Wong
%RBNZ
%August 2016
%
%Performs transformation of time series according to transformation code
%
%% CODE
% 1         No transformation
% 2         100 X difference in logs
% 3         Percent Change (if diff log approximation is poor)

%%
if transform_code == 1
    y = data;
elseif transform_code == 2
    y = [NaN; 100*diff(log(data))];
elseif transform_code == 3
    y = [NaN; 100*((data(2:end)./data(1:end-1))-1)];
end

end


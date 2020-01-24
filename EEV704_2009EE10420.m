function [action,EntryNo] =EEV704_2009EE10420(Data)
persistent Data_Seen i
persistent long short filter_iter


EntryNo = 'EEV_2009EE10420';
% replace EEV_XXX with your complete entry number
if(isempty(Data_Seen))
 i = 1;
 long = 0;
 short = 0;
 filter_iter = 1;
else
    filter_iter = filter_iter + 1;
    
    if(rem(filter_iter,10)~=1)
        action = 0;
        return;
    else
        i = i+1;
    end
end
Data_Seen(i,:) = Data; 

Close = Data(:,6);

if(length(Data_Seen)>50)
    %This value ‘50’ implies that it waits for 50 data points 
       
    MovAvgS = wma(Data_Seen(:,6),4);
    MovAvgL = wma(Data_Seen(:,6),50);
    RSI = rsi(Data_Seen(:,6),25);
           
    if( (MovAvgS(i) > MovAvgL(i)) && RSI(i) < 50)
        if (long == 0 && short == 0)
            action = 1;
            long = 1;
        elseif (short == 1 && long == 0)
            action = -11;
            short = 0;
        else
            action = 0;
        end
    elseif(MovAvgS(i) < MovAvgL(i))
        if (long == 1 && short == 0)
            action = 11;
            long = 0;
        elseif (short == 0 && long == 0)
            action = -1;
            short = 1;
        else
            action = 0;
        end
            
    else
        action = 0;
    end
   
else
    action = 0;
end
end

function out = rsi(data,period)
% Function to calculate the Relative Strength Index of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the length of the Wilder Smoothing window.
%
% Example:
% out = rsi(data,period)
%

% Error check
if nargin ~=2 
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The first input to ',mfilename,' must be a vector.']);
end
if numel(period) ~= 1
    error('The Wilder smoothing period must be a scalar.');
end
if length(data) < period+1
    error('The data set must be at least 1 element longer than the requested RSI period.');
end

% calculate the up and down data changes
dd = diff(data);
uc = dd;
uc(uc<0)=0;
dc = dd;
dc(dc>0)=0;
dc = -dc;
% perform Wilder Smoothing
wuc = wildersmoothing(uc,period);
wdc = wildersmoothing(dc,period);
% calculate the RSI (taking account of nan's in wuc and wdc)
out = [nan*ones(period-1,1); 100-(100./(1+wuc(period-1:end)./wdc(period-1:end)))];
end

function [Exposure_Net , No_of_Shares] =sell(Price,Exposure,Transaction_Cost)


No_of_Shares = floor(Exposure/Price) ;
Exposure_Net = Exposure + No_of_Shares*Price - Transaction_Cost ;

end

function out = sma(data,period)
% Function to calculate the simple moving average of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the number of periods over which to calculate the average
%
% Example:
% out = sma(data,period)
%

% Error check
if nargin ~= 2
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end
if (numel(period) ~= 1) || (mod(period,1)~=0)
    error('The period must be a scalar integer.');
end
if length(data) < period
    error('The length of the data must be at least the specified ''period''.');
end

% calculate the SMA
out = filter(ones(1,period),period,data);
out(1:period-1) = nan; % these are just the filter buffer filling
end

function out = wildersmoothing(data,period)
% Function to perform Wilder Smoothing on a data set
% 'data' is the vector of values to be smoothed.  The first element is assumed to be
% the oldest data.
% 'period' is the length of the smoothing window.
%
% Example:
% out = wildersmoothing(data,period)
%

% Error check
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The first input to ',mfilename,' must be a vector.']);
end
if numel(period) ~= 1
    error('The Wilder smoothing period must be a scalar.');
end

% perform the filtering
ld = length(data);
if ld < period
    error('The data vector must be at least as long as the required period of smoothing.');
elseif ld == period
    out = mean(data);
else
    out = nan*ones(size(data));
    out(period:end) = filter(1/period,[1 -(period-1)/period],data(period:end),sum(data(1:(period-1)))/period);
end
end

function out = wma(data,period)
% Function to calculate the weighted moving average of a data set
% The weight is based on the number of days in the moving average
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the number of periods over which to calculate the average
%
% Example:
% wma(data,period)
%

% Error check
if nargin ~= 2
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end
if (numel(period) ~= 1) || (mod(period,1)~=0)
    error('The period must be a scalar integer.');
end
if length(data) < period
    error('The length of the data must be at least the specified ''period''.');
end

% calculate the WMA
den = sum(1:period);
out = filter(period:-1:1,den,data);
out(1:period-1) = nan; % these are just the filter buffer filling

end

%Peaks detect functions
function peaks=mypeaks(data,w1,w2,sig1)
  %if the data is not nan add it together
  presum = zeros(1+length(data), 1);
  for i = 1:length(data)
    if ~isnan(data(i))
      presum(i+1) = presum(i)+data(i);
    else
      presum(i+1) = presum(i);
    end
  end

  %define the window and calculate the mean value of the window
  w = floor(round(length(data)/w1)/w2);
  mm = zeros(length(data), 1);
  for i = 1:length(data)
    s = max(i-w,1);
    t = min(i+w+1,length(data)+1);
    c = 0;
    for j = s:t-1
      if ~isnan(data(j))
        c = c+1;
      end
    end
    mm(i) = (presum(t)-presum(s))/c;
  end
  %calculate RMSE
  c = 0;
  dx = 0;
  for i = 1:length(data)
    if ~isnan(data(i))
      dx = dx + (data(i)-mm(i))*(data(i)-mm(i));
      c = c+1;
    end
  end
  dx = dx/c;
  sig = dx^(1/2);
  %find the peaks and return the values
  peaks = zeros(length(data),1);
  index = 1;
  for i = 1:length(data)
    if data(i) > mm(i)+sig/sig1 && (i==1 || data(i-1) < data(i)) && (i==length(data) || data(i) >= data(i+1))
      peaks(index) = i;
      index = index+1;
    end
  end
  
  peaks(index:end) = [];
  
end
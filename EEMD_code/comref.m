%This function is to calculate the respiration rate (bpm)
function complot = comref(data,t,refx)

    %tpeaks = t(data);
    tpeaks = data;
    complot = zeros(1,length(refx));
    for i = 1:length(refx)
   
        c = 0;
        for j = 1:length(tpeaks)
            if tpeaks(j)<=refx(i)
               c = c + 1; 
            end
        end
        complot(i) = (c/refx(i))*60;
        
    end

    complot(i+1:end) = [];
    
end
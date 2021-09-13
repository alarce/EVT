function [prob, EVTProb] = EVTfromPDF(data, threshold, pdf_x, pdf_y,tail, target)
%Extreme value theory -> Generalized pareto distribution

%TODO apply EVT to a PDFintead of raw data, is only supported for right
%tail, todo left tail
%% Argument validation
%input dataset must be a row vector
if (size(data,2) == 1) %column vector
    data = transpose(data);  %transpose
elseif(size(data,1) == 1) %row vector
    %do nothing
end
%check that tail is either 'right' or 'left'
%check that for right tail target is higher than threshold
%check that for left tail target is lower than threshold
%% Initialization
prob = 0;
EVTProb = 0;
    cdf_data(1) = 0;
for i = 2:size(pdf_x,2)
    cdf_data(i) = trapz(pdf_x(1,1:i), pdf_y(1,1:i));
end 
    index = find(pdf_x >= threshold);
    percentile = 100 * (cdf_data(index(1)));

%% logic
if(strcmp(tail, 'left'))  %left tail
    %threshold = percentile10;
    percentageTail = percentile/100; %probability mass in the tail we are modelling
    numExceedance = sum(data <= threshold);
    exceedanceIndex = data <= threshold;
    exceedance = threshold - data(exceedanceIndex);
    exceedance(exceedance == 0) = 0.001;
    paramEsts = gpfit(exceedance);
    kHat      = paramEsts(1);   % Tail index parameter
    sigmaHat  = paramEsts(2);   % Scale parameter
    
  
    
    %Probability: LSAEB-> collision || Vision-> N/A (early auto rear view) 
    numCollisions = sum(data <= target);
    rawcollisionProb = (numCollisions / size(data,2) ) * 100; %classical definition of probability
    rawcollisionProb = round(rawcollisionProb, 2);
    
    exceedanceTarget = threshold - target; %above this exceedance it is considered collision, conversion of the target value to the exceedance scale
    probCollisionInTail = 1 - gpcdf(exceedanceTarget ,kHat,sigmaHat);
    probCollisionLSAEB  = ( percentageTail * probCollisionInTail)*100;
    probCollisionLSAEB  = round(probCollisionLSAEB, 2);
    
    %{
    figure()
    ygrid = linspace(0,2*max(exceedance),100);
    line(ygrid,gppdf(ygrid,kHat,sigmaHat));
    hold on
    xlim([0, 2*max(exceedance)]);
    xlabel('Exceedance');
    ylabel('Probability Density');
    title('Extreme value theory (EVT-GPD): left tail');
    area(ygrid(ygrid >= exceedanceTarget), gppdf(ygrid(ygrid>=exceedanceTarget),kHat,sigmaHat));
    text = strcat('Prob collision = Prob(stopDist < 0 cm) = ', {' '},num2str(probCollisionLSAEB),'%');
    legend('Generalized Pareto Distribtion for left tail',text{1})
    
    dim = [.2 .5 .3 .3];
    str = strcat('raw Probability of collision: ', {' '}, num2str(rawcollisionProb), '%');
    annotation('textbox', dim, 'string', str{1}, 'FitBoxToText', 'on');
    %}
    
    prob = rawcollisionProb;
    EVTProb = probCollisionLSAEB;
    
else  %right tail
    %threshold = percentile90;
    percentageTail = 1 - (percentile/100); %probability mass in the tail we are modelling 
    numExceedance = sum(data >= threshold);
    exceedanceIndex = data >= threshold;
    exceedance = data(exceedanceIndex) - threshold;
    exceedance(exceedance == 0) = 0.001;
    paramEsts = gpfit(exceedance);
    kHat      = paramEsts(1);   % Tail index parameter
    sigmaHat  = paramEsts(2);   % Scale parameter
    
    
    %Probability: LSAEB-> nuisance stop || Vision: late automatic rear view
    %nuisanceThreshold = 90; %>90 cm  -> Target
    rawProbNuisanceLSAEB = (sum(data >= target) ./ numel(data))*100;
    rawProbNuisanceLSAEB = round(rawProbNuisanceLSAEB, 2);
    exceedanceTarget = target - threshold; %above this exceedance it is considered nuisance, conversion of the target value to the exceedance scale
    probNuisanceInTail = 1 - gpcdf(exceedanceTarget ,kHat,sigmaHat);
    probNuisanceLSAEB  = ( percentageTail * probNuisanceInTail)*100;
    probNuisanceLSAEB  = round(probNuisanceLSAEB, 2);
    
    %{
    figure()
    ygrid = linspace(0,2*max(exceedance),100);
    line(ygrid,gppdf(ygrid,kHat,sigmaHat));
    hold on
    xlim([0, 2*max(exceedance)]);
    xlabel('Exceedance');
    ylabel('Probability');
    title('Extreme value theory (EVT-GPD): right tail');
    area(ygrid(ygrid >= exceedanceTarget), gppdf(ygrid(ygrid>=exceedanceTarget),kHat,sigmaHat));
    text = strcat('Prob nuisance stop = Prob(stopDist > 90 cm) = ', {' '},num2str(probNuisanceLSAEB),'%');
    legend('Generalized Pareto Distribtion for right tail',text{1})
    %sanity check total probability is one
    %integral = trapz(ygrid,gppdf(ygrid,kHat,sigmaHat)); 
    dim = [.2 .5 .3 .3];
    str = strcat('raw Probability of nuisance stop: ', {''}, num2str(rawProbNuisanceLSAEB), '%');
    annotation('textbox', dim, 'string', str{1}, 'FitBoxToText', 'on');
    %}
    
    prob = rawProbNuisanceLSAEB;
    EVTProb = probNuisanceLSAEB;
end
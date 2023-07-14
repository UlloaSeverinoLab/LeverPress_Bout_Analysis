% This code can recognize a Lever press bout based on the inter-press interval (IPI) 
% probability distribution values calculated on the IPI distribution through neuroexplorer
%The whole workspace get saved in a matlab file at the end

clear BoutBegin_meanIPI BoutEnd_meanIPI LeverPress_TS

filename = 'FR10d3_WT_BoutAnalysis'; %Select file name for saving variables at the end
numberOfMice = 17; %Set parameters
IPI_pDistribution = [];
for i = 1: numberOfMice
    mouse = num2str(i);
    Name = "FR10d3_" + mouse; %Change according to mouse name
    disp(Name)
    LeverPress_TS = eval(Name);
    ISI = diff(LeverPress_TS);
    edges=0:0.2:10;
    h=histogram(ISI,edges, 'Normalization', 'probability');
    IPI_pDistribution(:,i) = h.Values';
    close Figure 1
end
Average_IPIDistribution = round(mean(IPI_pDistribution,2),8);
Q = quantile(Average_IPIDistribution,0.03);
plot(edges(2:end),Average_IPIDistribution)
hold on
plot(edges(Average_IPIDistribution == Q), Q,'*r')

prompt = 'Do you want to use previous threshold values? Y/N';
str = input(prompt, 's');
if str == 'N'
        savedThreshold(1:numberOfMice,1) = edges(Average_IPIDistribution == Q); %Set parameters if it is the first time you run the analysis
        threshold =savedThreshold(1:numberOfMice,1);
end
c=0;
for j = 1:1:numberOfMice
    if j~=0 %if there is a missing mouse number
    mouse = j;
    fprintf ('mouse #%d\n', mouse);
    clear BoutBegin BoutEnd LeverPress_TS 
    mouse= num2str(mouse);
    Name = "FR10d3_" + mouse; %Change according to mouse name
    LeverPress_TS = eval(Name); 
    if str == 'Y'
        threshold =savedThreshold(j,1);
    end
      
    BoutBegin(1,1) = LeverPress_TS(1);
    
    for i = 1:length(LeverPress_TS)-1
        
        DT_LP = LeverPress_TS(i+1)-LeverPress_TS(i);
        if DT_LP>=(threshold) %InterBoutInterval threshold
            BoutBegin(i+1,1) = LeverPress_TS(i+1);
            BoutEnd((i+1)-1,1)= LeverPress_TS(i);
            BoutBegin( ~any(BoutBegin,2), : ) = [];  %rows
            BoutEnd( ~any(BoutEnd,2), : ) = [];  %rows
        else
        end
    end
    BoutEnd(end+1)=LeverPress_TS(end);
    
    Boutlength = [];
    BoutDuration = [];
    InterBoutInterval = [];
    Bout_IPI = [];
    for k = 1:length(BoutEnd)
        Boutlength(k,1)= sum(LeverPress_TS<=BoutEnd(k) & LeverPress_TS>=BoutBegin(k));
        LeverID = find (LeverPress_TS<=BoutEnd(k) &  LeverPress_TS>=BoutBegin(k));
        Bout_IPI(k,1) = mean(diff(LeverPress_TS(LeverID(1:end))));
    end
    c=c+1;
    InterBoutInterval = BoutBegin(2:end)-BoutEnd(1:end-1);
    BoutDuration = BoutEnd(:)-BoutBegin(:);
    BoutDuration( ~any(BoutDuration,2), : ) = [];%remove zeros
    AllMice_BoutBegin{1,c} = BoutBegin;
    AllMice_BoutEnd{1,c} = BoutEnd;
    AllMice_BoutDuration{1,c}=BoutDuration;
    AllMice_InterBoutInterval{1,c}=InterBoutInterval;
    AllMice_BoutLength{1,c}=Boutlength;
    AllMice_Bout_IPI{1,c}=Bout_IPI;
    meanBoutDuration(c,1) = mean(BoutDuration);
    meanInterBoutInterval(c,1) = mean(InterBoutInterval);
    meanBoutLength(c,1) = mean(Boutlength);
    meanBout_IPI(c,1) = mean(Bout_IPI, 'omitnan');
    
    figure(2)
    clf
    plot(LeverPress_TS, 1, 'ob', BoutBegin,1.1,'*r', BoutEnd,1.1,'*g')
    xlim([0 300])
    ylim([0 1.5])
  
    else
    end
end
save(filename)

%% ISIdistribution graphs
ISI = diff(LeverPress_TS);
edges=0:0.2:10;
figure
h=histogram(ISI,edges, 'Normalization', 'probability');

% Plot the curve alone
%
% I could have retrieved values and edges from h (they are stored in the
% Values and BinEdges properties respectively) but I wanted to show how
% to get this information without actually creating a plot
[values, edges] = histcounts(ISI,edges, 'Normalization', 'probability');
centers = (edges(1:end-1)+edges(2:end))/2;
figure
plot(centers, values, 'k-')
% Plot both, line superimposed on the histogram
figure
h = histogram(ISI,edges, 'Normalization', 'probability');
hold on
plot(centers, values, 'k-')

%%

BoutBegin_1nex=AllMice_BoutBegin{1,1};
BoutBegin_2nex=AllMice_BoutBegin{1,2};
BoutBegin_3nex=AllMice_BoutBegin{1,3};
BoutBegin_4nex=AllMice_BoutBegin{1,4};
BoutBegin_5nex=AllMice_BoutBegin{1,5};
BoutBegin_6nex=AllMice_BoutBegin{1,6};
BoutBegin_7nex=AllMice_BoutBegin{1,7};
BoutBegin_8nex=AllMice_BoutBegin{1,8};
BoutBegin_9nex=AllMice_BoutBegin{1,9};
BoutBegin_10nex=AllMice_BoutBegin{1,10};
BoutBegin_11nex=AllMice_BoutBegin{1,11};
BoutBegin_12nex=AllMice_BoutBegin{1,12};






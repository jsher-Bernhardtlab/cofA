%Code to make a demograph from POLAR data

function m=newdemograph(cellList)
n=0; per= 8; % lower percentile fluorescence removed
pixsize= 0.064; maxlength=6; minlength=1.5; %micron
maxsize=round(maxlength/pixsize);
minsize=round(minlength/pixsize);
data=zeros(1000,2); 
relintarray1=zeros(1000,maxsize); %matrix for signal1
relintarray2=zeros(1000,maxsize); %matrix for signal2
leftpole1=zeros(1000,13);
rightpole1=zeros(1000,13);
leftpole2=zeros(1000,13);
rightpole2=zeros(1000,13);
for frame=1:length(cellList)
    if ~isempty(cellList{frame})
        for Cell=1:length(cellList{frame})
            if isfield(cellList{frame}{Cell},'length')
                if cellList{frame}{Cell}.length>minsize &&...
                        cellList{frame}{Cell}.length<5/pixsize &&...
                        length(cellList{frame}{Cell}.steplength)>1.5/pixsize
                    
                    signal1=cellList{frame}{Cell}.signal1; %get signal1
                    signal2=cellList{frame}{Cell}.signal2; %get ftsZ signal
                    steparea=cellList{frame}{Cell}.steparea; %get steparea
                    steplength=cellList{frame}{Cell}.steplength; %get steplength
                    nsignal1=signal1./steparea; % normalize by steparea
                    
                    if length(signal1)== length(steparea) && n < 4001 % quality control
                        % switch pole if required to put lower signal pole
                        % at the left

                        if sum(nsignal1(2:5))>sum(nsignal1(end-4:end-1))
                            signal1= signal1(end:-1:1);
                            signal2= signal2(end:-1:1);
                            steparea= steparea(end:-1:1);
                            steplength= steplength(end:-1:1);
                            cellList{frame}{Cell}.signal1=signal1;
                            cellList{frame}{Cell}.signal2=signal2;
                            cellList{frame}{Cell}.steparea=steparea;
                            cellList{frame}{Cell}.lengthvector=cumsum(steplength);
                        end
                        n=n+1;
                        nsignal1=signal1./steparea;
                        relint1= nsignal1/sum(nsignal1);
                        nsignal2=signal2./steparea;
                        relint2= nsignal2/sum(nsignal2);
                        k = floor(cellList{frame}{Cell}.length/2);
                        temp = cellList{frame}{Cell}.lengthvector-cellList{frame}{Cell}.length/2;
                        interpint1 = interp1(temp(1:length(relint1)),relint1,-k:k,'linear','extrap');
                        relintarray1(n,floor(0.5*maxsize-k:0.5*maxsize+k))=interpint1;
                        interpint2 = interp1(temp(1:length(relint2)),relint2,-k:k,'linear','extrap');
                        relintarray2(n,floor(0.5*maxsize-k:0.5*maxsize+k))=interpint2;
                        data(n,:)=[sum(steplength), sum(signal2)];
                        leftpole1(n,:)=interpint1(1:13);
                        leftpole2(n,:)=interpint2(1:13);
                        rightpole1(n,:)=interpint1(end-12:end);
                        rightpole2(n,:)=interpint2(end-12:end);
                    end
                end
            end
        end
    end
end
if n<4000
    index=data(:,1)>0;
    data=data(index,:); %remove unfilled array
    relintarray1=relintarray1(1:n,:);
    relintarray2=relintarray2(1:n,:);
    leftpole1=leftpole1(index,:);
    leftpole2=leftpole2(index,:);
    rightpole1=rightpole1(index,:);
    rightpole2=rightpole2(index,:);
end
%figure, hist(data(:,2),30)
% sort by signal2, remove top "per" percent
[~, ind]=sort(data(:,2));
sigtresh=floor(per/100*n);
data1=data(ind,:);
data=data1(sigtresh:end,:);
relintarray1=relintarray1(ind,:);
relintarray1=relintarray1(sigtresh:end,:);
relintarray2=relintarray2(ind,:);
relintarray2=relintarray2(sigtresh:end,:);
%%%

%sample k data points
[numbr,~]=size(relintarray1);
k=500;
if numbr>500
    sampling = randsample(numbr,k);
    leftpole1=leftpole1(sampling ,:);
    leftpole2=leftpole2(sampling ,:);
    rightpole1=rightpole1(sampling ,:);
    rightpole2=rightpole2(sampling ,:);
    relintarray1=relintarray1(sampling,:);
    relintarray2=relintarray2(sampling,:);
    data=data(sampling,1);
end

% sort by cell length
[~, ind]=sort(data(:,1));
relintarray1=relintarray1(ind,:);
relintarray2=relintarray2(ind,:);

% demographs
figure, subplot(1,2,1)
%image(relintarray1,'CDataMapping','scaled')
imagesc(relintarray1,[0 0.08])
set(gca,'TickDir','out')
xticks([1 .5*maxsize maxsize])
yticks([])
xticklabels([0 0.5*maxlength maxlength]-0.5*maxlength)

subplot(1,2,2)
%image(relintarray2,'CDataMapping','scaled')
imagesc(relintarray2,[0 0.07])
set(gca,'TickDir','out')
xticks([1 .5*maxsize maxsize])
yticks([])
xticklabels([0 0.5*maxlength maxlength]-0.5*maxlength)
[m,~]=size(relintarray2);

% plot averaged pole profiles
[~,n]=size(leftpole1);
figure, subplot(2,2,1), 
yyaxis left
plot(1:n, nanmean(leftpole1), '-g','linewidth', 2), 
disp('bait left: '), sum(nanmean(leftpole1))
ylim([0 0.08])
xlim([1 13])
hold on, yyaxis right
plot(1:n, nanmean(leftpole2), '-r','linewidth', 2)
disp('prey left: '), sum(nanmean(leftpole2))
ylim([0 0.06])
yticks([0 0.03 0.06])
xlim([1 13])
set(gca,'xtick', [1 7 13], 'ticklength', get(gca,'ticklength')*4, 'linewidth', 1)

subplot(2,2,2)
yyaxis left
plot(1:n, nanmean(rightpole1), '-g','linewidth', 2), 
disp('bait right: '), sum(nanmean(rightpole1))
ylim([0 0.08])
xlim([1 13])
hold on, yyaxis right
plot(1:n, nanmean(rightpole2), '-r','linewidth', 2)
disp('prey right: '), sum(nanmean(rightpole2))
ylim([0 0.06])
yticks([0 0.03 0.06])
xlim([1 13])

set(gca,'xtick', [1 7 13], 'ticklength', get(gca,'ticklength')*4, 'linewidth', 1)
end

                    
                    
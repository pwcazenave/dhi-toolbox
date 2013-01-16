function [statRes,res,phaseShift]=xcorrBoundary(dataIn,dataPredIn,posIn,posPredIn,orientation,plotFigs)
% Do a 1D cross-correlation for every point along each boundary for the
% time specified in timespan. Also calculate the RMS of the residual.

% Orientation (north-south or west-east boundary)
if strcmpi(orientation,'north')==1 || strcmpi(orientation,'south')==1
    colPos=1;
    yaxisLabel='Longitude';
elseif strcmpi(orientation,'west')==1 || strcmpi(orientation,'east')==1
    colPos=2;
    yaxisLabel='Latitude';
else
    error('Unknown orietation. Please choose ''west'', ''east'', ''south'' or ''north''.')
end

% Specify a timespan
timespan=1:24*14; % fortnight

% Start with the western boundary
time=dataIn(timespan,1); % always the same
phaseShift=nan(size(dataIn(:,2:end),2),1);
res=nan(size(time,1),size(dataIn(:,2:end),2));
statRes=nan(size(dataIn(:,2:end),2),2);
for i=1:size(dataPredIn(:,2:end),2) % don't want to include the time column
    dataUehara=dataIn(timespan,i+1);
    dataMIKE=dataPredIn(timespan+1,i+1); % shift MIKE data due to midnight start
    [c,lags]=xcorr(dataUehara,dataMIKE);
    phaseShift(i)=lags(c==max(c));
    % Residual analysis
    res(:,i)=dataMIKE-dataUehara;
    statRes(i,1)=sqrt(nanmean(res(:,i).^2));
    statRes(i,2)=std(res(:,i));
end
clear dataUehara dataMike c lags


if strcmpi(plotFigs,'yes')==1

    % Let's have a look see...
    close all

    % Compare a single tidal curve with the calculated residuals
    subfig(2,2,1)
    pointVal=ceil(size(posIn,1)/2); % Hlafway through the domain of lat/long
    plot(time,res(:,pointVal),'r-x','LineWidth',2)
    hold on
    plot(time,dataIn(timespan,pointVal),'b-x','LineWidth',2)
    plot(time,dataPredIn(timespan+1,pointVal),'k-x','LineWidth',2)
    % plot(time,phaseShift,'c-x','LineWidth',2)
    xlabel('Time (days)')
    ylabel('Height (m)')
    legend('Residual','Uehara','MIKE21','Location','SouthOutside','Orientation','horizontal')
    legend('boxoff')

    % Plot the spatially varying rms values
    subfig(2,2,3)
    imagesc(time,posIn(:,colPos),res')
    set(gca,'ydir','normal')
    hold on
    % Add the profile location above
    plot(time,ones(size(time,1))*posIn(pointVal,colPos),'k:','LineWidth',3)
    % Add the time varying rms and standard deviations in statRes
    plot(10*(statRes(:,1)),posIn(:,colPos),'k.','MarkerSize',20)
    plot(10*(statRes(:,1)),posIn(:,colPos),'w.','MarkerSize',5)
    % plot(zeros(size(posIn(:,colPos))),posIn(:,colPos),'w:','Linewidth',3)
    plot(ones(size(posIn(:,colPos)))*2,posIn(:,colPos),'k','Linewidth',3)
    plot(ones(size(posIn(:,colPos)))*2,posIn(:,colPos),'w','Linewidth',1)
    plot(ones(size(posIn(:,colPos))),posIn(:,colPos),'k','Linewidth',3)
    plot(ones(size(posIn(:,colPos))),posIn(:,colPos),'w','Linewidth',1)
    % Label the rms axes
    text(2*1.1,max(posIn(:,colPos))*0.99,'0.2 m','Rotation',270,'Color','k','FontWeight','bold')
    text(1*1.2,max(posIn(:,colPos))*0.99,'0.1 m','Rotation',270,'Color','k','FontWeight','bold')
    colorbar
    set(get(colorbar,'YLabel'),'String','Residual (m)')
    xlabel('Time (days)')
    ylabel(yaxisLabel)

    % Uehara's data from the western boundary (Atlantic)
    subfig(2,2,[2,4])
    subplot(2,1,1)
    imagesc(dataIn(timespan,1),posIn(:,colPos),dataIn(timespan,2:end)')
    set(gca,'ydir','normal')
    xlabel('Time (days)')
    ylabel(yaxisLabel)
    colorbar
    set(get(colorbar,'YLabel'),'String','Height (m)')
    hold on
    % Add the profile location from Figure 1 (pointVal).
    plot(time,ones(size(time,1))*posIn(pointVal,colPos),'k:','LineWidth',3)
    % MIKE21 data from the western boundary
    subplot(2,1,2)
    imagesc(dataIn(timespan,1),posPredIn(:,colPos),dataPredIn(timespan,2:end)')
    set(gca,'ydir','normal')
    xlabel('Time (days)')
    ylabel(yaxisLabel)
    colorbar
    set(get(colorbar,'YLabel'),'String','Height (m)')
    hold on
    % Add the profile location above
    plot(time,ones(size(time,1))*posIn(pointVal,colPos),'k:','LineWidth',3)
end

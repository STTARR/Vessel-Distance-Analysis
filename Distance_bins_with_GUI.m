%{
    INPUT: Set of csvs files corresponding to one slide each
    with many entries (one per cell)
    
    Purpose: data agglomeration and data visualization
    Distance bining
    
    Creates density scatterplot comparing 2 markers (A and B), uses a
    automatically generated threshold to determine positive populations

    Change DistanceMarker to whatever you'd like to calculate the distance to.
    Will look through all cells and their assigned mean distance mask value.
    To remove cells that are not of interest, uncomment and modify the
    rows_to_remove variable with the filter.
%}


clear
clc
%define number of distance bins, and width of each histogram bin. Largest
%distance bin will be number_of_bins*width_of_bins
number_of_bins = 20;
width_of_bins = 10;
%define output subfolder names
subfolder_tables='tables';
subfolder_graphs='graphs';
%% INPUT

% Provide root folder (containing all csv files exported by Definiens)
root = 'C:\Users\Mark Zaidi\Desktop\and more desktop stuff\hypoxia\2018_Frontiers_DistanceAnalysis\DC 215 L1\ThesisDistanceMapStatistics\215 l1'
root=uigetdir(root);
csv_file=dir([root,'\*.csv']);

%Setup dynamic GUI by basing options off of a sample csv in the set
sample=readtable(fullfile(root, csv_file(1).name));

% Populating GUI with column names from first loaded csv
GUI=vessel_distance_GUI;
ColumnNames=sort(sample.Properties.VariableNames);
GUI.XMarkerDropDown.Items=ColumnNames;
GUI.YMarkerDropDown.Items=ColumnNames;
GUI.DistanceMarkerDropDown.Items=ColumnNames;
GUI.MarkersofInterestListBox.Items=ColumnNames;
GUI.XCoordinateDropDown.Items=ColumnNames;
GUI.YCoordinateDropDown.Items=ColumnNames;


%obtain user input in GUI
uiwait(GUI.UIFigure)
fprintf('RESUMING\n')


%MarkerB is the marker to be plotted on the X axis (Marker X), and MarkerA
%is the marker to be plotted on the Y axis (Marker Y)
%Populating variables with GUI selection data
MarkerA=GUI.YMarkerDropDown.Value;
MarkerB=GUI.XMarkerDropDown.Value;
DistanceMarker=GUI.DistanceMarkerDropDown.Value;
variablenames=GUI.MarkersofInterestListBox.Value;
%If the spatial coordinate data is present, populate the x and y
%coordinates of every cell
xcoord=[];
ycoord=[];
do_coord_scatter=0;
if GUI.CellcoordinatespresentindataCheckBox.Value==1
    xcoord=GUI.XCoordinateDropDown.Value;
    ycoord=GUI.YCoordinateDropDown.Value;
    do_coord_scatter=1;
end
%clean up table by identifying unique columns used in analysis
columnsused=[MarkerA,MarkerB,DistanceMarker,variablenames,xcoord,ycoord];
[~,idx]=unique(columnsused);
columnsused_clean=columnsused(idx);
closereq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output directories (tables and graphs subfolders)
if ~isdir(fullfile(root,subfolder_tables));
    mkdir(fullfile(root,subfolder_tables));
end
if ~isdir(fullfile(root,subfolder_graphs));
    mkdir(fullfile(root,subfolder_graphs));
end

% Find all csv files in the directory, create and save graphs for each of
% them

for j=1:length(csv_file)
    
    file_name=csv_file(j).name;
    warning off
    data=readtable([root '\' file_name]);
    
     % Filter out any 0's for marker intensities used in log scatterplot
     rows_to_remove2=data.(MarkerA)==0 | data.(MarkerB)==0;
     data(rows_to_remove2,:) = [];
    %clear out any columns not specified by GUI selections
    all_names=data.Properties.VariableNames;
    %Specify which rows to delete from deletion selection
    columns2del=ismember(all_names,columnsused_clean);
    all_names(columns2del) = [];
    
    %clear variables not used in the analysis
    for q=1:length(all_names)
        data.(all_names{q})=[];
    end
    warning on
    filename_no_ext=strsplit(file_name,'.');
    %Convert any numbers accidentally stored as strings to doubles
    %First, iterate over each column to see wether it is a cell or a
    %double. If it's a cell, then check each of the elements if it matches
    %the string 'undefined'
    
    %clear out any rows containing features labelled as undefined by
    %Definiens
    for p=1:length(data.Properties.VariableNames)
        if iscell(data.(data.Properties.VariableNames{p}))
                rows_to_remove=strcmp(data.(data.Properties.VariableNames{p}),'undefined');
                data(rows_to_remove,:) = [];
        end
        %If an element is not NaN but is numerical value stored as a string, convert it into a double    
        curr=str2double(data.(data.Properties.VariableNames{p}));
        viable4conversion=~isnan(curr);
        if viable4conversion==1
            data.(data.Properties.VariableNames{p})=str2double(data.(data.Properties.VariableNames{p}));
        end

    end
    
    %% Manual user interaction
    %calculate thresholds for x and y markers based on the mean plus one
    %standard deviation. Values below are considered negative, values above
    %are positive
    EF5Gate=mean(data.(MarkerA))+std(data.(MarkerA));
    EDUGate=mean(data.(MarkerB))+std(data.(MarkerB));
    

    means_of_interest=removevars(data,DistanceMarker);
    %means_of_interest = data;
%     means_of_interest.MeanPerfusedVesselDistanceMap = [];
    
    MOI_names=means_of_interest.Properties.VariableNames;
    %% Bin cells based on the DistanceMarker variable
    for i = 0:number_of_bins-1
        % Perform binning for the x and y axis markers, as well as the
        % total number of cells, for later use in plotting
        bin_values_mean_edu(i+1)=mean(data.(MarkerB)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins));
        bin_values_mean_ef5(i+1)=mean(data.(MarkerA)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins));
        bin_value_cell_counts(i+1)=numel(data.(MarkerB)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins));
        
        %select all data in distance bin i, and perform the same binning as
        %described above
        mask=(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins);
        
        for m=1:width(means_of_interest)
            bin_value_MOI(i+1,m)=mean(means_of_interest.(MOI_names{m})(mask));
        end
        % Percent of total cells in each bin that are greater than gate
        % in each bin, find the total number of cells that are
        % greater than the lower bin limit & less than or equal to the higher
        % bin limit & are greater than their respective gate value. Take that
        % number of cells and divide by the total cells in that bin
        bin_values_percent_EDU_positive(i+1)=numel(data.(MarkerB)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins&data.(MarkerB)>EDUGate))./bin_value_cell_counts(i+1);
        bin_values_percent_EF5_positive(i+1)=numel(data.(MarkerA)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins&data.(MarkerA)>EF5Gate))./bin_value_cell_counts(i+1);
        
        %Calculating percent positive cells using total cells in the entire image as the
        %denominator
        bin_values_percent_EDU_positive_from_total(i+1)=numel(data.(MarkerB)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins&data.(MarkerB)>EDUGate))./height(data);
        bin_values_percent_EF5_positive_from_total(i+1)=numel(data.(MarkerA)(data.(DistanceMarker)>i*width_of_bins&data.(DistanceMarker)<=i*width_of_bins+width_of_bins&data.(MarkerA)>EF5Gate))./height(data);
    end
    %Preparing and writing summary tables
    bin_cap=width_of_bins:width_of_bins:width_of_bins*number_of_bins;
    final_values=[bin_cap; bin_values_mean_edu; bin_values_mean_ef5; bin_value_cell_counts; bin_values_percent_EDU_positive; bin_values_percent_EF5_positive; bin_values_percent_EDU_positive_from_total; bin_values_percent_EF5_positive_from_total;]';
    variable_names = {'Distance_to_vessel' (MarkerB) (MarkerA) 'Number_of_cells' ['Percent_' (MarkerB) '_positive'] ['Percent_' (MarkerA) '_positive'] ['Percent_' (MarkerB) '_positive_from_total'] ['Percent_' (MarkerA) '_positive_from_total']};
    t = array2table(final_values,'VariableNames',variable_names)    
    writetable(t,fullfile(root, subfolder_tables,file_name))
    writetable(array2table([bin_cap' bin_value_MOI],'VariableNames',['Distance_bin' MOI_names]),fullfile(root, subfolder_tables,['all_means_',file_name]))
    
    %% calculations for dscatter later on
     MarkerAgate = EF5Gate;
    MarkerBgate = EDUGate;
    gatelength = 500;   % Large arbitrary value to make gate line visible
    MarkerBgatelog=min(data.(MarkerA)(data.(MarkerA)>0)):0.001*255:500*255;
    MarkerAgatelog=min(data.(MarkerB)(data.(MarkerB)>0)):0.001*255:500*255;
    MarkerAgatelog(2,:) = MarkerAgate;
    MarkerBgatelog(2,:) = MarkerBgate;
    % Define positive and negative populations
    posMarkerA = data(data.(MarkerA)>MarkerAgate,:);
    negMarkerA = data(data.(MarkerA)<MarkerAgate,:);
    doublepos = posMarkerA(posMarkerA.(MarkerB)>MarkerBgate,:);
    posMarkerB = data(data.(MarkerB)>MarkerBgate,:);
    dp = height(doublepos);
    posMarkerAonly = height(posMarkerA)-dp;
    posMarkerBonly = height(posMarkerB)-dp;
    nn = height(negMarkerA)-posMarkerBonly;
    %% subplotting 
    
    % Define plot layout and display
    figure('units','normalized','outerposition',[0 0 1 1])
    
    % Cell count plot    
    subplot(3,2,1)
      
    bar(bin_cap,bin_value_cell_counts)
    title('Cell counts') 
    
    % Marker A and B subplots
    % Define upper bound for x-axis
    x_upper =(number_of_bins +1)* width_of_bins;

    
        %check if coordinate variables are present in table. If they are, make
    %a scatter plot where posMarkerB cells are in red, and posMarkerA cells
    %are green is (marker B = marker X = red)
    if do_coord_scatter==1
        subplot(3,2,2)
        fprintf('hai')
        scatter(data.(xcoord),data.(ycoord),1,'b')
        
        hold on
        scatter(posMarkerA.(xcoord),posMarkerA.(ycoord),1,'g')
        scatter(posMarkerB.(xcoord),posMarkerB.(ycoord),1,'r')
        scatter(doublepos.(xcoord),doublepos.(ycoord),1,'y')
        legend({'--','Y+','X+','++'},'Interpreter','none','FontSize',8,'Location','Best')
        axis tight
        axis xy
        axis image
        hold off
    end
    %mean intensity of marker B distance bin histogram
    subplot(3,2,3)
    bar(bin_cap,bin_values_mean_edu,'r')
    %axis([0 x_upper 0 50])
    title(['Mean ' (MarkerB)], 'interpreter','none')
    
    %mean intensity of marker A distance bin histogram
    subplot(3,2,4)
    bar(bin_cap,bin_values_mean_ef5,'g')
    %axis([0 x_upper 0 50])
    title(['Mean ' (MarkerA)], 'interpreter','none')
    
    %percent of cells positive for marker B in each distrance bin
    subplot(3,2,5)
    bar(bin_cap,bin_values_percent_EDU_positive,'r')
    %axis([0 x_upper 0 1])
    title(['Percent ' (MarkerB) ' positive per bin'], 'interpreter','none')
    
    %percent of cells positive for marker A in each distrance bin
    subplot(3,2,6)
    bar(bin_cap,bin_values_percent_EF5_positive,'g')
    %axis([0 x_upper 0 1])
    title(['Percent ' (MarkerA) ' positive per bin'], 'interpreter','none')
    

  
    
    
    
    % Save output file
    filename_no_ext = strsplit(file_name,'.csv');
    output_fname = strcat(root,'\',subfolder_graphs,'\',filename_no_ext{1},'.png')
    fprintf('Saving results in directory :: %s \n', output_fname);
    saveas(gcf,[output_fname])
    %% Producing the EDU-EF5 density scatterplots
    %Make scatterplot as separate plot
    figure;
    hold on
    
    % Use dscatter to compute scatter plot with color as density
    %marker transparency scale according to # of points. 25000 is an
    %arbitrary value that determines whether transparency plotting is
    %necessary, or if the point count is low enough to warrant complete
    %opaqueness
    if height(data)<25000
        transparency=1;
    else
        transparency=5000/height(data);
        %method of linearly scaling transparency based on density. The
        %value of 5000 can be changed as long as the value of transparency
        %is between 0 and 1. Lower values equate to a greater transparency
        %of points
    end
    hold on
    
    %Using the prewritten dscatter function to generate density
    %scatterplots. Type "help dscatter" for details on syntax
    
    %Currently, the script will output a transparency scatter plot. To
    %switch to a density scatter plot, set perform_dscatter to 1
    perform_dscatter=0;
    if perform_dscatter
        dscat=dscatter(data.(MarkerB),data.(MarkerA),'log',1);
    else
        scatter(data.(MarkerB),data.(MarkerA),50,data.(DistanceMarker),'marker','.','MarkerEdgeAlpha',transparency)
    end
    contour_plot=dscatter(data.(MarkerB),data.(MarkerA),'log',1,'plottype','contour');

    colorbar
    ylim([min(data.(MarkerA)(:)) max(data.(MarkerA)(:))])
    xlim([min(data.(MarkerB)(:)) max(data.(MarkerB)(:))])
    %option of autoscaling data on a per-image basis, or manually setting
    %fixed scale bounds as seen below. Set fixed_scale to 0 for
    %automaticalyl generated scales
    fixed_scale=1;
    if fixed_scale
        ylim([1 150])
        xlim([10 200])
    end
    %set color scale bounds. if fixed_colorscale=1, will use manually
    %defined color axes
    fixed_colorscale=1;
    if fixed_colorscale
        caxis([0 200])
    end
    % Use log scale
    set(gca,'xscale','log');
    set(gca,'yscale','log');


    % Plot gate lines
   
    axis manual
    plot(MarkerAgatelog(1,:),MarkerAgatelog(2,:),'Color','g');
    plot(MarkerBgatelog(2,:),MarkerBgatelog(1,:),'Color','r');
    
    %label axes with gate values
    xlabel(strcat((MarkerB),' gate:',num2str(MarkerBgate)),'Interpreter','none');
    ylabel(strcat((MarkerA),' gate:',num2str(MarkerAgate)),'Interpreter','none');
    

    hold off
    %write figure
    saveas(gcf,strcat(root,'\',subfolder_graphs,'\',filename_no_ext{1},'_scatterplot.png'))
    close all
    

    
    end
    fprintf('Done!\n')

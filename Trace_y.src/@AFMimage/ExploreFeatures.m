function [h] = ExploreFeatures(img, targetFig)

%
% DESCRIPTION
% – Displays the image and displays information on features upon selection
% such as filament segment length and average height. 
% – Also displays detailed visualisation on double clicking a filament.
% – Part of Trace_y by WFX
%
% USAGE
% Standard method usage
% >> img.ExploreFeatures;
%
% INPUTS
% img  –  AFMimage object. Implicit with dot notation
% targetFig  –  Target figure window handle to use
%
% OUTPUTS
% The image is ploted in a figure with traced/segmented features 
% highlighted
%
% DEPENDENCIES
% – Uses Statistics and Machine Learning toolbox
% – Method for Trace_y's @AFMimage/ object and uses other methods such as 
% DispImage.m and DispFilament.m 
% – Used also by other routines
%
% AUTHORS
% Wei-Feng Xue
%
% HISTORY
% 2024.08  –  WFX Trace_y update to add some interactivity in data
% visualisation and exploration.
%



%
% Defaults
%
doubleClick_timing = 0.4;


% Init so data
filaments = img.features.filaments;
imgSize = img.pixelPerLine;
sel = 0;
sel_prev = 0;
tt_prev = tic;


% Init the image
if ~exist('targetFig', 'var') || isempty(targetFig)
    h_fig = img.DispImage;
    h_fig.OuterPosition(3:4) = [960 960];
else
    h_fig = img.DispImage(targetFig);
end
%h_fig = img.DispImage;

set(h_fig, 'NumberTitle', 'off', 'Name', 'Trace_y – ExploreFeatures');
hold('on');

% Find the image, should only be one in the DispImage figure.
h_img = findobj(h_fig, 'Type', 'Image');
h_img.ButtonDownFcn = @SelectFeature;

h_sel = [];
h_sel_e1 = [];
h_sel_e2 = [];
h_dt = [];

% Plot the selected filament with details
h_filaments = cell(numel(img.features.filaments), 1);
for aa = flip(1:numel(img.features.filaments))
    filament = img.features.filaments(aa);
    h_filaments{aa} = plot(filament.x, filament.y, 'c-', 'LineWidth', 1);
    h_filaments{aa}.ButtonDownFcn = @SelectFeature;
end



% Callback function for feature selection
    function [] = SelectFeature(~, evt)

        tt = tic;
        dtt = toc(tt_prev);

        xx = evt.IntersectionPoint(1);
        yy = evt.IntersectionPoint(2);

        % Find closest filament from click coordinates
        minDist = zeros(numel(filaments), 1);
        for aaa = 1:numel(filaments)
            allDist = pdist2([xx yy], [filaments(aaa).x filaments(aaa).y]);
            minDist(aaa) = min(allDist);
        end

        % Action if sufficiently close to a traced filament
        if min(minDist) < imgSize/50
            sel = find(minDist == min(minDist), 1);

            selected_filament = filaments(sel);

        else
            sel = 0;
        end

        if sel && sel ~= sel_prev

            % Disp info on the selected
            % Test code
            %fprintf('tt %g, dtt %g, sel, %g, sel_prev %g, select\n', tt, dtt, sel, sel_prev);

            delete(h_sel_e1);
            delete(h_sel_e2);
            delete(h_sel);
            h_sel_e1 = plot(selected_filament.x(1), selected_filament.y(1), ...
                'b+', 'MarkerSize', 10, 'LineWidth', 2);
            h_sel_e2 = plot(selected_filament.x(end), selected_filament.y(end), ...
                'r+', 'MarkerSize', 10, 'LineWidth', 2);
            h_sel = plot(selected_filament.x, selected_filament.y, 'r-', 'LineWidth', 2);
            h_sel.ButtonDownFcn = @SelectFeature;

            % Position and show a datatip popup
            delete(h_dt);
            if max(selected_filament.x(:))-min(selected_filament.x(:)) ...
                    > 0.38*(max(selected_filament.y(:))-min(selected_filament.y(:)))
                %h_dt = datatip(h_filaments{filament_sel}, max(selected_filament.x), ...
                %    selected_filament.y(selected_filament.x == max(selected_filament.x)));
                dt_x = max(selected_filament.x);
                dt_y = selected_filament.y(selected_filament.x == max(selected_filament.x));
            else
                %h_dt = datatip(h_filaments{filament_sel}, ...
                %    selected_filament.x(selected_filament.y == max(selected_filament.y)), ...
                %    max(selected_filament.y));
                dt_x = selected_filament.x(selected_filament.y == max(selected_filament.y));
                dt_y = max(selected_filament.y);
            end

            % Format the datatip. Setting value is a bit awkward as the
            % value variable must be the same in size as the data in the
            % object
            h_dt = plot(dt_x, dt_y, 'r.');
            %h_dtbox = datatip(h_dt, dt_x, dt_y);
            datatip(h_dt, dt_x, dt_y);
            %dtRows = [dataTipTextRow('Filament', 0*selected_filament.x+filament_sel)];
            %h_filaments{filament_sel}.DataTipTemplate.DataTipRows = dtRows;
            %set(h_filaments{filament_sel}.DataTipTemplate, ...
            %    'FontName', 'helvetica', 'FontSize', 14);
            dtRows = [...
                dataTipTextRow('Filament', sel); ...
                dataTipTextRow(sprintf('Length / %s:', selected_filament.xyUnit), selected_filament.lContour); ...
                dataTipTextRow(sprintf('Mean height / %s:', selected_filament.zUnit)', mean(selected_filament.z)); ...
                ];
            if ~isempty(selected_filament.helicalFilament3DModel)
                dtRows = [dtRows; ...
                    dataTipTextRow(sprintf('Periodicity / %s', selected_filament.xyUnit), selected_filament.helicalFilament3DModel.periodicity); ...
                    dataTipTextRow('Handedness', string(selected_filament.helicalFilament3DModel.handedness))];
            end
            h_dt.DataTipTemplate.DataTipRows = dtRows;
            set(h_dt.DataTipTemplate, 'FontName', 'helvetica', 'FontSize', 12);

            sel_prev = sel;
            tt_prev = tt;

        elseif sel && sel == sel_prev && dtt < doubleClick_timing

            % Double click displays detailed graphics on the selected in 
            % separate figs
            % Test code
            %fprintf('tt %g, dtt %g, sel, %g, sel_prev %g, double click\n', tt, dtt, sel, sel_prev);

            img.DispFilament(sel, 'nm');
            hh = findall(0, 'Type', 'figure', 'Number', 12);
            if ~isempty(hh)
                close(12);
            end
            if ~isempty(selected_filament.helicalFilament3DModel)
                img.DispFilament3DModel(sel);
            end

        elseif sel && sel == sel_prev && dtt > doubleClick_timing

            % Check if first click of a double click 
            % Test code
            %fprintf('tt %g, dtt %g, sel, %g, sel_prev %g, double click (?)\n', tt, dtt, sel, sel_prev);

            tt_prev = tt;

        elseif ~sel

            % Unselect if not clicking on anything
            % Test code
            %fprintf('tt %g, dtt %g, sel, %g, sel_prev %g, unselect\n', tt, dtt, sel, sel_prev);

            % Unselect
            delete(h_sel_e1);
            delete(h_sel_e2);
            delete(h_sel);
            delete(h_dt);
            h_sel = [];
            h_dt = [];
            sel = 0;
            sel_prev = sel;
            tt_prev = tt;

        end

    end



if nargout == 1
    h = h_fig;
end



end
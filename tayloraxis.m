function h = tayloraxis(varargin)
%TAYLORAXIS Create axis for Taylor diagram
%
% h = tayloraxis;
% h = tayloraxis(ax, ...)
%
% This function sets up a polar axis to use for a Tayor diagram (Taylor,
% 2001, J Geophys Res Atmos, 106:7183?7192).  These diagrams are useful to
% compactly show the skill of one or more datasets (typically models)
% relative to a reference dataset (typically observations). In this plot,
% the r coordinate corresponds to the standard deviation of a dataset and
% the theta parameter to the correlation between the dataset and a
% reference dataset; the centered RMSD is proportional to the
% distance between the plotted reference point and the plotted model
% points.
%
% To plot data points, use any polar-plotting function (polarplot,
% polarscatter, etc), with r = standard deviation and theta =
% real(acos(correlation)).  Similarly, correlation and standard deviation
% tick/axis lines can be altered from their default values by adjusting
% the R- and Theta- properties of the polar axis.   
%
% Input variables:
%
%   ax:         handle to axis to be turned into a Taylor diagram polar
%               axis. If not included, current axis will be used (or a new
%               axis created if none exist).    
%
% Optional input variables (passed as parameter/value pairs)
%
%   stdmax:     maximum r-axis limit (i.e. maximum standard deviation of
%               datasets to be plotted) 
%
%   npanel:     number of panels to show.  If npanel = 1, the theta limits
%               are set to [0 pi/2], corresponding to positive
%               correlations.  If npanel = 2, the theta limits are set to
%               [0 pi], allowing for negative correlations to be plotted.   
%
%   stdref:     reference value to use when plotting RMSD circles.  This
%               should correspond to the standard deviation of the
%               reference dataset.
%
%   rmsdtick:   scalar tick interval, or vector of RMSD tick values, to be
%               plotted as semi-circles extending from reference point.
%
% Output variables:
%
%   h:          1 x 1 structure of graphics handles:
%
%               ax:         PolarAxes object, polar axis handle
%
%               rmsdline:   n x 1 array of Line objects, RMSD contour lines

% Copyright 2017 Kelly Kearney

% Parse and check input

p = inputParser;
p.addOptional('ax', gca);
p.addParameter('stdmax', 1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('npanel', 1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('stdref', 1, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
p.addParameter('rmsdtick', NaN, @(x) validateattributes(x, {'numeric'}, {'vector', 'increasing'}));

p.parse(varargin{:});
Opt = p.Results;

if ~ismember(Opt.npanel, [1 2])
    error('npanel input must be either 1 or 2');
end
if Opt.npanel == 1
    thlim = [0 pi/2];
    cortk = [1 .99 .95 .9:-.1:0]';
else
    thlim = [0 pi];
    cortk = [1 .99 .95 .9:-.1:0 -.1:-.1:-.9 -.95 -.99 -1]';
end

% Start with polar axes

axes(Opt.ax);

thetatk = real(acos(cortk));
thetatklbl = strtrim(cellstr(num2str(cortk)));

polarplot(NaN, NaN);
h.ax = gca;

set(h.ax, 'ThetaAxisUnits', 'radians', ...
          'RLim', [0 Opt.stdmax], ...
          'ThetaLim', thlim, ...
          'ThetaTick', thetatk, ...
          'ThetaTickLabel', thetatklbl, ...
          'RColor', 'b');      
      
% Add RMSD circles

thrmsd = linspace(0, pi, 100);
if isscalar(Opt.rmsdtick)
    if isnan(Opt.rmsdtick)
        Opt.rmsdtick = mean(diff(get(h.ax, 'RTick')));
    end
    if Opt.npanel == 1
        rmsdmax = sqrt(Opt.stdref.^2 + Opt.stdmax.^2);
    else
        rmsdmax = Opt.stdref + Opt.stdmax;
    end
    rrmsd = 0:Opt.rmsdtick:rmsdmax;
else
    rrmsd = Opt.rmsdtick;
end
[th,rr] = ndgrid(thrmsd, rrmsd);
x = Opt.stdref + rr.*cos(th);
y = rr.*sin(th);
[th,rr] = cart2pol(x,y);


if Opt.npanel == 1
    isout = rr > h.ax.RLim(2) | th > pi/2;
else
    isout = rr > h.ax.RLim(2);
end
rr(isout) = NaN;

hold(h.ax, 'on');
h.rmsdline = polarplot(th,rr, 'color', [0.082353 0.6902 0.10196], 'linewidth', 0.1, 'linestyle', ':');





function cmap = f_getColorMap(varargin)


cmap = [ ...
    0.8863, 0.3373, 0.1294  %1 red   
    0.929 , 0.694 , 0.125   %2 yellow 
    0.2431, 0.6314, 0.3412  %3 green 
    0.2118, 0.5137, 0.7333  %4 blue
    0.000 , 0.447 , 0.741   %5 dark blue 
    0.4549, 0.4314, 0.6863  %6 lila
    0.494 , 0.184 , 0.556   %7 dark lila
    0.3882, 0.3882, 0.3882  %8 gray
    0.635 , 0.078 , 0.184   %9 brown
    1.000 , 1.000 , 1.000   %10 white
    1.000 , 0.000 , 1.000   %11 magenta 
    0.750 , 0.000 , 0.750
    0.750 , 0.750 , 0.000
    0.250 , 0.250 , 0.250 ];



if nargin ~= 0
    
    n = varargin{1};
    
    if n > size(cmap,1)
        
        if ~isempty(get(groot,'Children'))
            
            cmap = colormap(jet(n));
            
        else
            
            cmap = colormap(jet(n));
            
            close gcf;
            
        end
        
        % END IF ~isempty
        
    end
    
    % END IF n <= size(cmap,1)
    
end

% END IF nargin ~= 0



end
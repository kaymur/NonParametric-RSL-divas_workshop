function [p] = lim_indic(y,z,direction,unc,tri_lim)

%% p = 0 or is based on the triangle distribution
	defval('tri_lim',100);

    if direction == -1 % marine limiting
        prop = y + 2 * unc;
        if prop >= z
            p = (prop-z)/tri_lim;
            %p = (z-prop)/100;
            %p = 1;
        else
            p = 0;
        end
    elseif direction == 1
        prop = y - 2 * unc;
        if prop <= z % terrestrial limiting
            p = (z-prop)/tri_lim;
            %p = (prop-z)/100;
            %p=1;
        else
            p=0;
        end
    else
        p=0;
    end
end
function [] = MakeSpatFreqFig(rows)
%---------------------------------------------------------------------------------------------------------------------------------
% MakeSpatFreqFig: Compiles spatial frequency vs speed subplots into one figure
    % INPUTS:
        % -
    % OUTPUTS:
        % -
%---------------------------------------------------------------------------------------------------------------------------------
    pp = 1;
    while 1
        Analyze_Asymmetry_Control(rows,pp)
        x = input('Continue: ');
        switch x
            case 1
            case 0
                disp('DONE')
                break
        end
        pp = pp + 1;
    end
%---------------------------------------------------------------------------------------------------------------------------------   
end


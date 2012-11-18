function data = vcTestManual(fname);

    if nargin == 0
        fname = '';
    end

    % when does the vc test pulse start in each trace in samples
    par.protocolParams.stepStartInd = 1139;

    % how wide is the vc test pulse in samples
    par.protocolParams.stepWidth = 50000;

    % how high is the voltage step in mV
    par.protocolParams.voltageDelta = 5;

    % which channel contains the voltage
    par.chVCNum = 1;

    % the holding potential for this cell when the protocol begins
    % (i.e. the baseline current in the traces keeps the cell at
    % at this potential, used to estimate the natural resting potential)
    par.holding = -70;

    data = processVCTest(fname, par);

end

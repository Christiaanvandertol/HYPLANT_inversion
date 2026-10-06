function retval = isOctave
    %% https://docs.octave.org/v9.1.0/How-to-Distinguish-Between-Octave-and-Matlab.html
    persistent cacheval;  % speeds up repeated calls
    
    if isempty (cacheval)
    cacheval = (exist ("OCTAVE_VERSION", "builtin") > 0);
    end
    
    retval = cacheval;
end
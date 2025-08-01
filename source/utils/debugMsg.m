function debugMsg(msg)
%DEBUGHEADER Prints a formatted header message for debugging
    info = dbstack;
    if numel(info) > 1
        caller = info(2).name;
    else
        caller = info(1).name;
    end

    line = repmat('-',1,60);
    fprintf('\n%s\n[%s] %-50s\n%s\n', line, caller, msg, line);
end
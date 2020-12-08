function testSwitchString(txt)

switch txt
    case 'initialized'
        disp('initialized')
    case 'running'
        disp('running')
    case 'stopped'
        disp('stopped')
    case 'completed'
        disp('completed')
    otherwise
end
clear all;
clear classes;
close all;
beep off;

test_all_dir = dir('../src/*.m');
for test_all_ll=1:length(test_all_dir)
    
    test_all_classname = test_all_dir(test_all_ll).name(1:end-2);
    name               = ['./' test_all_classname '_TEST.m'];
    
    if exist(name,'file')
        run(name);
    else
        warning('TEST_ALL:NoScriptFound', 'No testscript found for class %s', test_all_classname);
    end
    
end

fprintf('\n');
fprintf('[0] All tests successful \n');
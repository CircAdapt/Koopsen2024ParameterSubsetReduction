%% CircAdapt()
% Makes C++ CircAdapt Object
% CircAdaptCompile() should run first to compile the code
function CA = CircAdapt()
    CA = mex_interface(str2fun('CircAdaptMatlab'));
end
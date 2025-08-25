function getfunc(funcname)
    file = which(funcname);
    if isempty(file)
        error([funcname, ' function not found. Please make sure it is in your MATLAB path.']);
    end
    dest = '/mnt/cup/people/bichanw/SpikeSorting/Codes/biliRRR/final_repos/GCSA/matlab/tools/';
    system(['cp ', file, ' ', dest]);
end
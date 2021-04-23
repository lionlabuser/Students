function diskspace
% Query to know how many free GB are left on the C:\
FileObj      = java.io.File('C:\');
free_bytes   = FileObj.getFreeSpace/1e+9;
total_bytes  = FileObj.getTotalSpace/1e+9;
usable_bytes = FileObj.getUsableSpace/1e+9;
fprintf('Disk space for C:\\ >> %3.1f GB free out of %3.1f GB total available (%3.1f%% used space)\n',free_bytes,total_bytes,(total_bytes-free_bytes)/total_bytes*100)
end
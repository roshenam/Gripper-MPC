function writeLineDat(fname_dat, ...
     A, B, K, R)

fid = fopen(fname_dat,'w+');

fwrite(fid, A, 'double');        % Matrix<Type, 8, 8> &A
fwrite(fid, B, 'double');        % Matrix<Type, 8, 3> &B
fwrite(fid, K, 'double');        % Matrix<Type, 3, 8> &K
fwrite(fid, R, 'double');        % Matrix<Type, 8, 6> &R

fclose(fid);

disp(['Parameters written to: ' fname_dat])

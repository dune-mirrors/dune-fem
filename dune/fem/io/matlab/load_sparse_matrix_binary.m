function A = load_sparse_matrix_binary(filename)
%function A = load_sparse_matrix(filename)
%
% load binary file, which represents a sparse matrix, e.g. generated
% by saveSparseMatrixBinary. Two magic numbers are read (double and int)
% to check the compatibility of the
% binary file. If binary format does not match, the reading must be 
% modified in this file.
%
% format:
% magic number "DSM" for Dune Sparse Matrix
% magic numbers: int 111, double 111
% number of rows, cols and maxnoonzero_per_row
% number of total_nonzeros
% for 0.. total_nonzeros-1 : triples (int r,int c,double v)   
% where r,c start from 0
% "EOF" as marker of EOF
  
% Bernard Haasdonk 15.12.2006

  fid = fopen(filename,'r');
  % if standard reading is not the correct format for a given binary
  % file, activate the following:
  %fid = fopen(filename,'r','ieee-be');
  
  magicstr = char(fread(fid,3,'char'))';
  if ~isequal(magicstr,'DSM')
    error('read magicstr doe not indicate Dune Sparse Matrix!');
  end;

  magicint = fread(fid,1,'int');
  magicdouble = fread(fid,1,'double');
  
  if (magicint~=111) | (magicdouble~=111.0)
    error(['magic numbers not read correctly, change the binary format in' ...
	   ' this reading routine!!']);
  end;
  
  nrows = fread(fid,1,'int');
  ncols = fread(fid,1,'int');
  nnonzeros = fread(fid,1,'int');
  ntotalnonzeros = fread(fid,1,'int');
  
  disp(['generating ',num2str(nrows),'x',num2str(ncols),...
	' sparse matrix with ',num2str(ntotalnonzeros),' totalnonzeros.']);
  A = sparse(nrows,ncols,nnonzeros);
  
  for i=1:ntotalnonzeros
    row = fread(fid,1,'int');
    col = fread(fid,1,'int');
    val = fread(fid,1,'double');
    A(row+1,col+1) = val;
  end;  
  
  eofstr = char(fread(fid,3,'char'))';
  if ~isequal(eofstr,'EOF')
    error('read eofstr does not indicate end of binary file!');
  end;

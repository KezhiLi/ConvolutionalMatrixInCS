function golay_2d=ext_golay_2D(img_size)
row=img_size(1);
col=img_size(2);

% Current program only works for even length sequence;
n0_row=floor(row/2);
n0_col=floor(col/2);


zad_row=zeros(row,1);
zad_col=zeros(col,1);

ind=1:n0_row;
ind_col=1:n0_col;

[ga_row,gb_row]=generate_golay(floor(log2(n0_row)));
[ga_col,gb_col]=generate_golay(floor(log2(n0_col)));

zad_row(1:n0_row)=ga_row;
zad_row(n0_row+1)=1;
zad_row(row:-1:n0_row+2)=zad_row(n0_row:-1:2);


zad_col(1:n0_col)=gb_col;
zad_col(n0_col+1)=1;
zad_col(col:-1:n0_col+2)=zad_col(n0_col:-1:2);

AA=ones(row,col);
golay_2d=diag(zad_row)*AA*diag(zad_col);
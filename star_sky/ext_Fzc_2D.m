function fzc_2d=ext_Fzc_2D(img_size,root_2D)
row=img_size(1);
col=img_size(2);

% Current program only works for even length sequence;
n0_row=floor(row/2);
n0_col=floor(col/2);
root_row=root_2D(1);
root_col=root_2D(2);

zad_row=zeros(row,1);
zad_col=zeros(col,1);

ind=1:n0_row;
ind_col=1:n0_col;

zad_row(1:n0_row)=exp(-j*pi*root_row*(ind-1).^2/row);
zad_row(n0_row+1)=1;
zad_row(row:-1:n0_row+2)=conj(zad_row(n0_row:-1:2));


zad_col(1:n0_col)=exp(-j*pi*root_col*(ind_col-1).^2/col);
zad_col(n0_col+1)=1;
zad_col(col:-1:n0_col+2)=conj(zad_col(n0_col:-1:2));

AA=ones(row,col);
fzc_2d=diag(zad_row)*AA*diag(conj(zad_col));
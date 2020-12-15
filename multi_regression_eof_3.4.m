X=load('[X+Y+]data.tsv');
nx=48; ny=12; nt=58;
Jan1983ssta=X(33,:);
Jan1983ssta_xy=reshape(Jan1983ssta,nx,ny);
Jan1989ssta=X(39,:);
Jan1989ssta_xy=reshape(Jan1989ssta,nx,ny);
Jan1998ssta=X(48,:);
Jan1998ssta_xy=reshape(Jan1998ssta,nx,ny);
Jan2002ssta=X(52,:);
Jan2002ssta_xy=reshape(Jan2002ssta,nx,ny);
long = 32.5:5:267.5;
lat = -27.5:5:27.5;
range = -30:30;
figure(1);
subplot(2,1,1);
[X1,Y1] = meshgrid(long,lat);
contourf(X1, Y1, Jan1983ssta_xy',range);
colorbar;
caxis([-5 5]);
title('Jan 1983 SSTa');
subplot(2,1,2)
contourf(X1, Y1, Jan1989ssta_xy',range)
colorbar
caxis([-5 5])
title('Jan 1989 SSTa')
figure(2)
subplot(2,1,1)
contourf(X1, Y1, Jan1998ssta_xy',range)
colorbar
caxis([-5 5])
title('Jan 1998 SSTa')
subplot(2,1,2)
contourf(X1, Y1, Jan2002ssta_xy',range)
colorbar
caxis([-5 5])
title('Jan 2002 SSTa')
iocean=find(X(1,:) ~= 1.00000000E30);
Xo=X(: , iocean);

Corr=corrcoef(Xo);
[EOFs,ve]=eig(Corr);
Var=real( diag(ve./trace(ve)) );
EOFLand=nan(nx*ny,2);
EOFLand(iocean,1)=EOFs(:,1);
EOFLand(iocean,2)=EOFs(:,2);
PCs=Xo*EOFs;
PC1 = PCs(:,1);
PC2 = PCs(:,2);
eof1_xy = reshape(EOFLand(:,1),nx,ny);
eof2_xy = reshape(EOFLand(:,2),nx,ny);
year = 1951:2008;
figure(3)
subplot(2,1,1)
[X1,Y1] = meshgrid(long,lat);
contourf(X1, Y1, eof1_xy')
caxis([-0.1 0.1])
colorbar
title('EOF1')
subplot(2,1,2)
plot(year,PC1,'-*')
axis([1951 2008 -30 30])
grid on
title('PC1')
figure(4)
subplot(2,1,1)
[X1,Y1] = meshgrid(long,lat);
contourf(X1, Y1, eof2_xy')
caxis([-0.1 0.1])
colorbar
title('EOF2')
subplot(2,1,2)
plot(year,PC2,'-*')
axis([1951 2008 -30 30])
grid on
title('PC2')

Cova = cov(Xo);
[EOFs1,ve] = eig(Cova);
Var=real( diag(ve./trace(ve)) );
EOFLand=nan(nx*ny,2);
EOFLand(iocean,1)=EOFs1(:,end);
EOFLand(iocean,2)=EOFs1(:,end-1);
PCs=Xo*EOFs1;
PC1 = PCs(:,end);
PC2 = PCs(:,end-1);
eof1_xy = reshape(EOFLand(:,1),nx,ny);
eof2_xy = reshape(EOFLand(:,2),nx,ny);
year = 1951:2008;
figure(3)
subplot(2,1,1)
[X1,Y1] = meshgrid(long,lat);
contourf(X1, Y1, eof1_xy')
caxis([-0.1 0.1])
colorbar
title('EOF1')
subplot(2,1,2)
plot(year,PC1,'-*')
axis([1951 2008 -30 30])
grid on
title('PC1')
figure(4)
subplot(2,1,1)
[X1,Y1] = meshgrid(long,lat);
contourf(X1, Y1, eof2_xy')
caxis([-0.1 0.1])
colorbar
title('EOF2')
subplot(2,1,2)
plot(year,PC2,'-*')
axis([1951 2008 -30 30])
grid on
title('PC2')
variance=var*100
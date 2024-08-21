function [mereni,dh,fi] = ul_0(nazev,m)

RAD=pi/200;
fid = fopen(nazev,'r');
mereni=fscanf(fid,'%d %f %f %f',[4,inf])';
fclose(fid);

fi = 0.00998*mereni(:,2)/1000.*sin(mereni(:,4)*RAD)*RAD;
dh=mereni(:,2).*sin(mereni(:,4)*RAD-fi)*m;
end


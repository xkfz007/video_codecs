a=dlmread('tmp.txt');
%%
frs=10;
for i=1:frs
    s=sum(a(i*99-98:i*99));
    disp(s);
end
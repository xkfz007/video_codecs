%this file is used for the icip paper experiment.
%here we present the performance of ssim, stvssim,ssim3d under the quality metric stvssim,ssim and psnr.
clear
clc

yuvname_seq2={'akiyo_qcif.yuv','bridge_qcif.yuv','bus_qcif.yuv','carphone_qcif.yuv',...
         'coastguard_qcif.yuv','container_qcif.yuv','football_qcif.yuv','foreman_qcif.yuv',...
         'grandma_qcif.yuv','hall_qcif.yuv','harbour_qcif.yuv','highway_qcif.yuv',...
         'ice_qcif.yuv','mobile_qcif.yuv','news_qcif.yuv','salesman_qcif.yuv',...
         'silent_qcif.yuv','soccer_qcif.yuv','stefan_qcif.yuv','suzie_qcif.yuv'};
yuvname=yuvname_seq2;
len=length(yuvname);
%here are the data of different rdo methods.
%these data files are the selected results, which are good for the icip2012 paper.
avgdata_all_mserdo;
%
qp=28:42;
%%
%psnr
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
	%mserdo 
    rate=[seqname '_rate'];
    eval([rate '=avg_mserdo(:,1)*30/49000;']);
    psnrY=[seqname '_psnrY'];
    eval([ psnrY '=avg_mserdo(:,2);']);
%   rate=avg_mserdo(:,1)*30/49000;
%   psnrY=avg_mserdo(:,2);
%   plot(rate,psnrY,'*k-');
   
%   pause
end
%%
%ssim
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
	%mserdo 
    rate=[seqname '_rate'];
    eval([rate '=avg_mserdo(:,1)*30/49000;']);
    ssimY=[seqname '_ssimY'];
    eval([ ssimY '=avg_mserdo(:,5);']);
   
%   pause
end
%%
%stvssim
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
	%mserdo 
    rate=[seqname '_rate'];
    eval([rate '=avg_mserdo(:,1)*30/49000;']);
    stvssimY=[seqname '_stvssimY'];
    eval([ stvssimY '=avg_mserdo(:,11);']);
   
%   pause
end

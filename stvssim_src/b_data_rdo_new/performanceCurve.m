%this file is used for the icip paper experiment.
%here we present the performance of ssim, stvssim,ssim3d under the quality metric stvssim,ssim and psnr.
clear
clc

yuvname_seq1={'akiyo_qcif.yuv','bridge_qcif.yuv','carphone_qcif.yuv',...
    'coastguard_qcif.yuv','container_qcif.yuv','football_qcif.yuv',...
    'foreman_qcif.yuv','grandma_qcif.yuv','hall_qcif.yuv',...
    'highway_qcif.yuv','mobile_qcif.yuv','news_qcif.yuv','salesman_qcif.yuv',...
    'silent_qcif.yuv','soccer_qcif.yuv','suzie_qcif.yuv'};
yuvname_seq2={'akiyo_qcif.yuv','bridge_qcif.yuv','bus_qcif.yuv','carphone_qcif.yuv',...
         'coastguard_qcif.yuv','container_qcif.yuv','football_qcif.yuv','foreman_qcif.yuv',...
         'grandma_qcif.yuv','hall_qcif.yuv','harbour_qcif.yuv','highway_qcif.yuv',...
         'ice_qcif.yuv','mobile_qcif.yuv','news_qcif.yuv','salesman_qcif.yuv',...
         'silent_qcif.yuv','soccer_qcif.yuv','stefan_qcif.yuv','suzie_qcif.yuv'};
seqList=[5 7 8 13 14 18 19];
%seqList=[2 6 8 9 15 16 ];
yuvname=yuvname_seq2(seqList);
len=length(yuvname);
%here are the data of different rdo methods.
%these data files are the selected results, which are good for the icip2012 paper.
avgdata_all_mserdo;
avgdata_all_ssimrdo_lam1_2_seq2;
%avgdata_all_stvssimrdo_lam2_3_seq1_psnr;
avgdata_all_stvssimrdo_lam3_15;
%avgdata_all_stvssimrdo_lam2_3_r26;
%
%qp=28:42;
qp=33:42;
range=6:15;
%%
%1. firestly, we use stvssim as the metric(y-axis)
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
    aa_ssimrdo=[seqname '_avg_ssimrdo'];
    eval(['avg_ssimrdo=' aa_ssimrdo ';']);
%    aa_ssim3drdo=[seqname '_avg_ssim3drdo'];
%   eval(['avg_ssim3drdo=' aa_ssim3drdo ';']);
    aa_stvssimrdo=[seqname '_avg_stvssimrdo'];
    eval(['avg_stvssimrdo=' aa_stvssimrdo ';']);
    figure('name',[seqname '_stvssim_rate'])
    grid on
    hold on 
	%mserdo 
    rate=avg_mserdo(range,10)*30/49000;
    stvssimY=avg_mserdo(range,7);
    plot(rate,stvssimY,'*k-');
	%ssimrdo
    rate=avg_ssimrdo(range,10)*30/49000;
    stvssimY=avg_ssimrdo(range,7);
    plot(rate,stvssimY,'*b-');
	%ssim3drdo
%   rate=avg_ssim3drdo(:,10);
%   stvssimY=avg_ssim3drdo(:,7);
%   plot(rate,stvssimY,'*y-');
    %stvssimrdo
    rate=avg_stvssimrdo(range,10)*30/49000;
    stvssimY=avg_stvssimrdo(range,7);
    plot(rate,stvssimY,'*r-');
   
    pause
end
%%
%2. secondly, ssim is used.
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
    aa_ssimrdo=[seqname '_avg_ssimrdo'];
    eval(['avg_ssimrdo=' aa_ssimrdo ';']);
%    aa_ssim3drdo=[seqname '_avg_ssim3drdo'];
%    eval(['avg_ssim3drdo=' aa_ssim3drdo ';']);
    aa_stvssimrdo=[seqname '_avg_stvssimrdo'];
    eval(['avg_stvssimrdo=' aa_stvssimrdo ';']);
    figure('name',[seqname '_ssim_rate'])
    grid on
    hold on 
	%mserdo
    rate=avg_mserdo(qp+1,10);
    ssimY=avg_mserdo(qp+1,1);
     plot(rate,ssimY,'*k-');
	 %ssimrdo
    rate=avg_ssimrdo(:,10);
    ssimY=avg_ssimrdo(:,1);
    plot(rate,ssimY,'*b-');
	%ssim3drdo
%    rate=avg_ssim3drdo(:,10);
%    ssimY=avg_ssim3drdo(:,1);
%    plot(rate,ssimY,'*y-');
	%stvssimrdo
    rate=avg_stvssimrdo(:,10);
    ssimY=avg_stvssimrdo(:,1);
    plot(rate,ssimY,'*r-');
   
    pause
end
%%
%3. thirdly, psnr is used. we are hoping the profermance
for i=1:len
    seqname=yuvname{i}(1:end-4);
    aa_mserdo=[seqname '_avg_mserdo'];
    eval(['avg_mserdo=' aa_mserdo ';']);
    aa_ssimrdo=[seqname '_avg_ssimrdo'];
    eval(['avg_ssimrdo=' aa_ssimrdo ';']);
%    aa_ssim3drdo=[seqname '_avg_ssim3drdo'];
%    eval(['avg_ssim3drdo=' aa_ssim3drdo ';']);
    aa_stvssimrdo=[seqname '_avg_stvssimrdo'];
    eval(['avg_stvssimrdo=' aa_stvssimrdo ';']);
    figure('name',[seqname '_ssim_rate'])
    grid on
    hold on 
	%mserdo
    rate=avg_mserdo(qp+1,10);
    psnrY=avg_mserdo(qp+1,11);
     plot(rate,psnrY,'*k-');
	 %ssimrdo
    rate=avg_ssimrdo(:,10);
    psnrY=avg_ssimrdo(:,11);
    plot(rate,psnrY,'*b-');
	%ssim3drdo
%    rate=avg_ssim3drdo(:,10);
%    psnrY=avg_ssim3drdo(:,11);
%    plot(rate,psnrY,'*y-');
	%stvssimrdo
    rate=avg_stvssimrdo(:,10);
    psnrY=avg_stvssimrdo(:,11);
    plot(rate,psnrY,'*r-');
   
    pause
end



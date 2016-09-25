%novelty-SSIMº∆À„
sender_mail = 'programrunning@163.com'; %
password = '123456program';  %
smtp='smtp.163.com';
receiver_mail = 'programrunning@163.com'; %
%
setpref('Internet','E_mail',sender_mail);
setpref('Internet','SMTP_Server',smtp);
setpref('Internet','SMTP_Username',sender_mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465')

%
fullPath=pwd;
[pathstr,kindname,ext]=fileparts(fullPath);


exfile='..\lencod.exe';
reso='cif';
path='e:/sequences/';
%yuvname={'akiyo_qcif.yuv','bridge_qcif.yuv','carphone_qcif.yuv',...
%         'coastguard_qcif.yuv','container_qcif.yuv','football_qcif.yuv',...
%         'foreman_qcif.yuv','grandma_qcif.yuv','hall_qcif.yuv',...
%         'highway_qcif.yuv','mobile_qcif.yuv','news_qcif.yuv',...
%         'salesman_qcif.yuv','silent_qcif.yuv','soccer_qcif.yuv','suzie_qcif.yuv'};
%yuvname={'akiyo_qcif.yuv','bridge_qcif.yuv','bus_qcif.yuv','carphone_qcif.yuv',...
%         'coastguard_qcif.yuv','container_qcif.yuv','football_qcif.yuv','foreman_qcif.yuv',...
%         'grandma_qcif.yuv','hall_qcif.yuv','harbour_qcif.yuv','highway_qcif.yuv',...
%         'ice_qcif.yuv','mobile_qcif.yuv','news_qcif.yuv','salesman_qcif.yuv',...
%         'silent_qcif.yuv','soccer_qcif.yuv','stefan_qcif.yuv','suzie_qcif.yuv'};
yuvname={'football','foreman','news'};

len=length(yuvname);
%FrameSkip             = [1 1 0 1  1 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1]; 
%FrameSkip             = [2 2 0 1  1 2 1 1  2 2 1 2  1 1 2 2  2 1 1 1]; 
FrameSkip=zeros(len,0);
seqList=1:len;

FramesToBeEncoded     = 100;
FrameRate             = 30.0;
SourceWidth           = 176*2; 
SourceHeight          = 144*2;
MDDistortion          = 2 ; 
NumberReferenceFrames = 5 ;
NumberBFrames         = 0;  
QPBSlice              = 30;   
SymbolMode            =  0;  
DistortionSSIM        =  1;  
SSIMOverlapSize       =  4;   
WeightY               =  1;  
WeightCb              =  1;   
WeightCr              =  1;  
LambdaFunc            = 2;
precmd=[exfile ' -f encoder.cfg' ' -p FramesToBeEncoded=' num2str(FramesToBeEncoded) ' -p FrameRate=' num2str(FrameRate) ' -p SourceWidth=' num2str(SourceWidth) ' -p SourceHeight=' num2str(SourceHeight) ' -p MDDistortion=' num2str(MDDistortion) ' -p NumberReferenceFrames=' num2str(NumberReferenceFrames) ' -p NumberBFrames=' num2str(NumberBFrames) ' -p QPBSlice=' num2str(QPBSlice) ' -p SymbolMode=' num2str(SymbolMode) ' -p DistortionSSIM=' num2str(DistortionSSIM) ' -p SSIMOverlapSize=' num2str(SSIMOverlapSize) ' -p WeightY=' num2str(WeightY) ' -p WeightCb=' num2str(WeightCb) ' -p WeightCr=' num2str(WeightCr) ' -p LambdaFunc=' num2str(LambdaFunc) ];
QP=28:42;
finalcontent='';
for i=seqList
    dirname=strcat(yuvname{i},'_',reso);
    if(7~=exist(dirname,'dir'))
        mkdir(dirname);
    end
    InputFile=strcat(path,reso,'/',dirname{i},'.yuv');   
    inputcmd=[' -p InputFile=' InputFile ' -p FrameSkip=' num2str(FrameSkip(i)) ];
     
    cd(dirname);
    for qp=QP
        if qp<10
            qpstr=['0' num2str(qp)];
        else
            qpstr=num2str(qp);
        end
       infix=[dirname '_',qpstr];
       TraceFile =['trace_' infix '_enc.txt'];
       ReconFile =['rec_' infix '.yuv'];
       OutputFile=['out_' infix '.264'];
       StatsFile =['stats_',infix '.dat'];
       InfoFile  =['info_',infix, '.txt'];
       AvgFile   =['avg_',infix,'.txt'];
       MbRDdataFile   =['mbrddata_',infix,'.txt'];
      midcmd=[' -p QPISlice=' qpstr ' -p QPPSlice=' qpstr ' -p TraceFile =' TraceFile  ' -p ReconFile =' ReconFile  ' -p OutputFile=' OutputFile ' -p StatsFile=' StatsFile ' -p InfoFile=' InfoFile ' -p AvgFile=' AvgFile ' -p MbRDdataFile=' MbRDdataFile ];

      cmd=[precmd inputcmd midcmd];
      disp(cmd); 
      %system(cmd); 
      
      if qp==33 || qp==42
      subject=[kindname ':' dirname '_' qpstr];
      content=cmd;
	  flag=0;count=0;
	  while flag~=1&&count<10
          count=count+1;
          try
              sendmail(receiver_mail,subject,content);
              flag=1;
          catch exception
              flag=0;
              disp(exception);
          end
      end
      end
             
    end
    
    cd ..
    finalcontent=[finalcontent ' ' yuvname{i}];
end
subject=['Complete:' kindname];
%finalcontent=[fullPath ':\n' finalcontent];
finalcontent=sprintf('%s:\n%s',fullPath,finalcontent);
flag=0;count=0;
while flag~=1&&count<10
      count=count+1;
      try
         sendmail(receiver_mail,subject,finalcontent);
         flag=1;
      catch exception
         flag=0;
         disp(exception);
      end
end
exit

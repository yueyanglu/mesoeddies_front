function y=randdf(size1,df1,flag)
% This function is to generate random numbers
% according to a user defined probability density function (pdf)
% or cumulative distribution function (cdf).
% pdf or cdf is described by a matrix,
% whose size is N-by-2.
% Sampling points of pdf or cdf form the second row.
% Function value of pdf or cdf form the first row.
% coded by halleyhit on Aug. 15th, 2018
% Email: halleyhit@sjtu.edu.cn or halleyhit@163.com
% Syntax
% y = randdf(S,D,F)
% S - Size of dimension, integer values.
% Example: S=10 creates a 10-by-1 array
% Example: S=[10,2] creates a 10-by-2 matrix
% D - Density function, numeric matrix
% F - Flag, 'pdf' or 'cdf'
% Example:
% x=[-1:0.01:1];% Sampling points
% y=2*(x<0).*(x>-0.1)+4*(x<0.5).*(x>0.3);% Function value of pdf
% plot(x,y,'black')
% r=randdf([10000],[y;x],'pdf'); % generate random numbers
% hold on
% h=histogram(r);
% h.Normalization='pdf';
% h.BinWidth=0.01;
% h.EdgeColor='none'; % view the pdf of the generated numbers
if numel(size1)==1
    n=1;m=size1(1);
elseif numel(size1)==2
    n=size1(1);m=size1(2);
else
    return
end
all=n*m;
if size(df1,1)~=2
    return
end
if sum(df1(1,:)<0)>0
    return
end
if size(df1,2)<2
    return
end
if strcmp(flag,'pdf')
    df1(1,:)=cumsum(df1(1,:))/sum(df1(1,:));
elseif strcmp(flag,'cdf')
    if sum(diff(df1(1,:))<0)>0
        return
    end
    df1(1,:)=df1(1,:)/df1(1,end);
else
    return
end
df1(1,:)=df1(1,:)+[1:size(df1,2)]*eps;
temp=rand(1,all);
temp=interp1(df1(1,:),df1(2,:),temp,'linear','extrap');
y=reshape(temp,n,m);
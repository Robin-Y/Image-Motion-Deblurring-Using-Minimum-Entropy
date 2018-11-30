% 2D dual-tree complex wavelet transform and adaptive svd image denoising (additive Gaussian noise);
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
sigma = 10;
B_1 = B_0 + sigma*randn(size(B_0));
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
w = cplxdual2D(B_1, 1, Faf, af);
% detail coefficients
CH_1_real_up = w{1}{1}{1}{1};
CV_1_real_up = w{1}{1}{1}{2};
CD_1_real_up = w{1}{1}{1}{3};
CH_2_real_down = w{1}{1}{2}{1};
CV_2_real_down = w{1}{1}{2}{2};
CD_2_real_down = w{1}{1}{2}{3};
CH_1_imag_up = w{1}{2}{1}{1};
CV_1_imag_up = w{1}{2}{1}{2};
CD_1_imag_up = w{1}{2}{1}{3};
CH_2_imag_down = w{1}{2}{2}{1};
CV_2_imag_down = w{1}{2}{2}{2};
CD_2_imag_down = w{1}{2}{2}{3};

p_1 = size(CH_1_real_up);
p_2 = size(CV_1_real_up);
p_3 = size(CD_1_real_up);
p_4 = size(CH_2_real_down);
p_5 = size(CV_2_real_down);
p_6 = size(CD_2_real_down);
p_7 = size(CH_1_imag_up);
p_8 = size(CV_1_imag_up);
p_9 = size(CD_1_imag_up);
p_10 = size(CH_2_imag_down);
p_11 = size(CV_2_imag_down);
p_12 = size(CD_2_imag_down);
% standard deviation estimation (using the real component of 45 degree direction subband wavelet coefficients)
CD_1_real_up_re = reshape(CD_1_real_up,1,p_3(1)*p_3(2));
sigma_1_HH = median(abs(CD_1_real_up_re))/0.6745;
% zero-padding of real and complex detail subbands
if  mod(p_1(1),8)~= 0
    CH_1_real_up(p_1(1)+1:p_1(1)+8-mod(p_1(1),8),:) = 0;
end

if  mod(p_1(2),8)~= 0
    CH_1_real_up(:,p_1(2)+1:p_1(2)+8-mod(p_1(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_2(1),8)~= 0
    CV_1_real_up(p_2(1)+1:p_2(1)+8-mod(p_2(1),8),:) = 0;
end

if  mod(p_2(2),8)~= 0
    CV_1_real_up(:,p_2(2)+1:p_2(2)+8-mod(p_2(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_3(1),8)~= 0
    CD_1_real_up(p_3(1)+1:p_3(1)+8-mod(p_3(1),8),:) = 0;
end

if  mod(p_3(2),8)~= 0
    CD_1_real_up(:,p_3(2)+1:p_3(2)+8-mod(p_3(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_4(1),8)~= 0
    CH_2_real_down(p_4(1)+1:p_4(1)+8-mod(p_4(1),8),:) = 0;
end

if  mod(p_4(2),8)~= 0
    CH_2_real_down(:,p_4(2)+1:p_4(2)+8-mod(p_4(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_5(1),8)~= 0
    CV_2_real_down(p_5(1)+1:p_5(1)+8-mod(p_5(1),8),:) = 0;
end

if  mod(p_5(2),8)~= 0
    CV_2_real_down(:,p_5(2)+1:p_5(2)+8-mod(p_5(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_6(1),8)~= 0
    CD_2_real_down(p_6(1)+1:p_6(1)+8-mod(p_6(1),8),:) = 0;
end

if  mod(p_6(2),8)~= 0
    CD_2_real_down(:,p_6(2)+1:p_6(2)+8-mod(p_6(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_7(1),8)~= 0
    CH_1_imag_up(p_7(1)+1:p_7(1)+8-mod(p_7(1),8),:) = 0;
end

if  mod(p_7(2),8)~= 0
    CH_1_imag_up(:,p_7(2)+1:p_7(2)+8-mod(p_7(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_8(1),8)~= 0
    CV_1_imag_up(p_8(1)+1:p_8(1)+8-mod(p_8(1),8),:) = 0;
end

if  mod(p_8(2),8)~= 0
    CV_1_imag_up(:,p_8(2)+1:p_8(2)+8-mod(p_8(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_9(1),8)~= 0
    CD_1_imag_up(p_9(1)+1:p_9(1)+8-mod(p_9(1),8),:) = 0;
end

if  mod(p_9(2),8)~= 0
    CD_1_imag_up(:,p_9(2)+1:p_9(2)+8-mod(p_9(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_10(1),8)~= 0
    CH_2_imag_down(p_10(1)+1:p_10(1)+8-mod(p_10(1),8),:) = 0;
end

if  mod(p_10(2),8)~= 0
    CH_2_imag_down(:,p_10(2)+1:p_10(2)+8-mod(p_10(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_11(1),8)~= 0
    CV_2_imag_down(p_11(1)+1:p_11(1)+8-mod(p_11(1),8),:) = 0;
end

if  mod(p_11(2),8)~= 0
    CV_2_imag_down(:,p_11(2)+1:p_11(2)+8-mod(p_11(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_12(1),8)~= 0
    CD_2_imag_down(p_12(1)+1:p_12(1)+8-mod(p_12(1),8),:) = 0;
end

if  mod(p_12(2),8)~= 0
    CD_2_imag_down(:,p_12(2)+1:p_12(2)+8-mod(p_12(2),8)) = 0;
end
% size of each real and complex detail subbands after zero-padding
s_1 = size(CH_1_real_up);
s_2 = size(CV_1_real_up);
s_3 = size(CD_1_real_up);
s_4 = size(CH_2_real_down);
s_5 = size(CV_2_real_down);
s_6 = size(CD_2_real_down);
s_7 = size(CH_1_imag_up);
s_8 = size(CV_1_imag_up);
s_9 = size(CD_1_imag_up);
s_10 = size(CH_2_imag_down);
s_11 = size(CV_2_imag_down);
s_12 = size(CD_2_imag_down);

%making 8*8 Blocks
CH_1_real_up_cell = cell((s_1(1)/8),(s_1(2)/8));
for i=0:s_1(1)/8-1
    for j=0:s_1(2)/8-1
        CH_1_real_up_cell{i+1,j+1}=CH_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_1_real_up_cell = cell((s_2(1)/8),(s_2(2)/8));
for i=0:s_2(1)/8-1
    for j=0:s_2(2)/8-1
        CV_1_real_up_cell{i+1,j+1}=CV_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_1_real_up_cell = cell((s_3(1)/8),(s_3(2)/8));
for i=0:s_3(1)/8-1
    for j=0:s_3(2)/8-1
        CD_1_real_up_cell{i+1,j+1}=CD_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_2_real_down_cell = cell((s_4(1)/8),(s_4(2)/8));
for i=0:s_4(1)/8-1
    for j=0:s_4(2)/8-1
        CH_2_real_down_cell{i+1,j+1}=CH_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_2_real_down_cell = cell((s_5(1)/8),(s_5(2)/8));
for i=0:s_5(1)/8-1
    for j=0:s_5(2)/8-1
        CV_2_real_down_cell{i+1,j+1}=CV_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_2_real_down_cell = cell((s_6(1)/8),(s_6(2)/8));
for i=0:s_6(1)/8-1
    for j=0:s_6(2)/8-1
        CD_2_real_down_cell{i+1,j+1}=CD_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_1_imag_up_cell = cell((s_7(1)/8),(s_7(2)/8));
for i=0:s_7(1)/8-1
    for j=0:s_7(2)/8-1
        CH_1_imag_up_cell{i+1,j+1}=CH_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_1_imag_up_cell = cell((s_8(1)/8),(s_8(2)/8));
for i=0:s_8(1)/8-1
    for j=0:s_8(2)/8-1
        CV_1_imag_up_cell{i+1,j+1}=CV_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_1_imag_up_cell = cell((s_9(1)/8),(s_9(2)/8));
for i=0:s_9(1)/8-1
    for j=0:s_9(2)/8-1
        CD_1_imag_up_cell{i+1,j+1}=CD_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_2_imag_down_cell = cell((s_10(1)/8),(s_10(2)/8));
for i=0:s_10(1)/8-1
    for j=0:s_10(2)/8-1
        CH_2_imag_down_cell{i+1,j+1}=CH_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_2_imag_down_cell = cell((s_11(1)/8),(s_11(2)/8));
for i=0:s_11(1)/8-1
    for j=0:s_11(2)/8-1
        CV_2_imag_down_cell{i+1,j+1}=CV_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_2_imag_down_cell = cell((s_12(1)/8),(s_12(2)/8));
for i=0:s_12(1)/8-1
    for j=0:s_12(2)/8-1
        CD_2_imag_down_cell{i+1,j+1}=CD_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end
% s, v, and d matrices for each block of each subband
ss_1 = cell((s_1(1)/8),(s_1(2)/8));
vv_1 = cell((s_1(1)/8),(s_1(2)/8));
dd_1 = cell((s_1(1)/8),(s_1(2)/8));

ss_2 = cell((s_2(1)/8),(s_2(2)/8));
vv_2 = cell((s_2(1)/8),(s_2(2)/8));
dd_2 = cell((s_2(1)/8),(s_2(2)/8));

ss_3 = cell((s_3(1)/8),(s_3(2)/8));
vv_3 = cell((s_3(1)/8),(s_3(2)/8));
dd_3 = cell((s_3(1)/8),(s_3(2)/8));

ss_4 = cell((s_4(1)/8),(s_4(2)/8));
vv_4 = cell((s_4(1)/8),(s_4(2)/8));
dd_4 = cell((s_4(1)/8),(s_4(2)/8));

ss_5 = cell((s_5(1)/8),(s_5(2)/8));
vv_5 = cell((s_5(1)/8),(s_5(2)/8));
dd_5 = cell((s_5(1)/8),(s_5(2)/8));

ss_6 = cell((s_6(1)/8),(s_6(2)/8));
vv_6 = cell((s_6(1)/8),(s_6(2)/8));
dd_6 = cell((s_6(1)/8),(s_6(2)/8));

ss_7 = cell((s_7(1)/8),(s_7(2)/8));
vv_7 = cell((s_7(1)/8),(s_7(2)/8));
dd_7 = cell((s_7(1)/8),(s_7(2)/8));

ss_8 = cell((s_8(1)/8),(s_8(2)/8));
vv_8 = cell((s_8(1)/8),(s_8(2)/8));
dd_8 = cell((s_8(1)/8),(s_8(2)/8));

ss_9 = cell((s_9(1)/8),(s_9(2)/8));
vv_9 = cell((s_9(1)/8),(s_9(2)/8));
dd_9 = cell((s_9(1)/8),(s_9(2)/8));

ss_10 = cell((s_10(1)/8),(s_10(2)/8));
vv_10 = cell((s_10(1)/8),(s_10(2)/8));
dd_10 = cell((s_10(1)/8),(s_10(2)/8));

ss_11 = cell((s_11(1)/8),(s_11(2)/8));
vv_11 = cell((s_11(1)/8),(s_11(2)/8));
dd_11 = cell((s_11(1)/8),(s_11(2)/8));

ss_12 = cell((s_12(1)/8),(s_12(2)/8));
vv_12 = cell((s_12(1)/8),(s_12(2)/8));
dd_12 = cell((s_12(1)/8),(s_12(2)/8));
% applying SVD to each block of each real and complex detail subbands.
for i=1:s_1(1)/8
    for j=1:s_1(2)/8
    [ss_1{i,j} vv_1{i,j} dd_1{i,j}] = svd(CH_1_real_up_cell{i,j});
    end
end

for i=1:s_2(1)/8
    for j=1:s_2(2)/8
    [ss_2{i,j} vv_2{i,j} dd_2{i,j}] = svd(CV_1_real_up_cell{i,j});
    end
end

for i=1:s_3(1)/8
    for j=1:s_3(2)/8
    [ss_3{i,j} vv_3{i,j} dd_3{i,j}] = svd(CD_1_real_up_cell{i,j});
    end
end

for i=1:s_4(1)/8
    for j=1:s_4(2)/8
    [ss_4{i,j} vv_4{i,j} dd_4{i,j}] = svd(CH_2_real_down_cell{i,j});
    end
end

for i=1:s_5(1)/8
    for j=1:s_5(2)/8
    [ss_5{i,j} vv_5{i,j} dd_5{i,j}] = svd(CV_2_real_down_cell{i,j});
    end
end

for i=1:s_6(1)/8
    for j=1:s_6(2)/8
    [ss_6{i,j} vv_6{i,j} dd_6{i,j}] = svd(CD_2_real_down_cell{i,j});
    end
end

for i=1:s_7(1)/8
    for j=1:s_7(2)/8
    [ss_7{i,j} vv_7{i,j} dd_7{i,j}] = svd(CH_1_imag_up_cell{i,j});
    end
end

for i=1:s_8(1)/8
    for j=1:s_8(2)/8
    [ss_8{i,j} vv_8{i,j} dd_8{i,j}] = svd(CV_1_imag_up_cell{i,j});
    end
end

for i=1:s_9(1)/8
    for j=1:s_9(2)/8
    [ss_9{i,j} vv_9{i,j} dd_9{i,j}] = svd(CD_1_imag_up_cell{i,j});
    end
end

for i=1:s_10(1)/8
    for j=1:s_10(2)/8
    [ss_10{i,j} vv_10{i,j} dd_10{i,j}] = svd(CH_2_imag_down_cell{i,j});
    end
end

for i=1:s_11(1)/8
    for j=1:s_11(2)/8
    [ss_11{i,j} vv_11{i,j} dd_11{i,j}] = svd(CV_2_imag_down_cell{i,j});
    end
end

for i=1:s_12(1)/8
    for j=1:s_12(2)/8
    [ss_12{i,j} vv_12{i,j} dd_12{i,j}] = svd(CD_2_imag_down_cell{i,j});
    end
end

%changing svd
vvv_1 = cell((s_1(1)/8),(s_1(2)/8));
sss_1 = cell((s_1(1)/8),(s_1(2)/8));
ddd_1 = cell((s_1(1)/8),(s_1(2)/8));
CH_1_real_up_cell_new = cell((s_1(1)/8),(s_1(2)/8));

vvv_2 = cell((s_2(1)/8),(s_2(2)/8));
sss_2 = cell((s_2(1)/8),(s_2(2)/8));
ddd_2 = cell((s_2(1)/8),(s_2(2)/8));
CV_1_real_up_cell_new = cell((s_2(1)/8),(s_2(2)/8));

vvv_3 = cell((s_3(1)/8),(s_3(2)/8));
sss_3 = cell((s_3(1)/8),(s_3(2)/8));
ddd_3 = cell((s_3(1)/8),(s_3(2)/8));
CD_1_real_up_cell_new = cell((s_3(1)/8),(s_3(2)/8));

vvv_4 = cell((s_4(1)/8),(s_4(2)/8));
sss_4 = cell((s_4(1)/8),(s_4(2)/8));
ddd_4 = cell((s_4(1)/8),(s_4(2)/8));
CH_2_real_down_cell_new = cell((s_4(1)/8),(s_4(2)/8));

vvv_5 = cell((s_5(1)/8),(s_5(2)/8));
sss_5 = cell((s_5(1)/8),(s_5(2)/8));
ddd_5 = cell((s_5(1)/8),(s_5(2)/8));
CV_2_real_down_cell_new = cell((s_5(1)/8),(s_5(2)/8));

vvv_6 = cell((s_6(1)/8),(s_6(2)/8));
sss_6 = cell((s_6(1)/8),(s_6(2)/8));
ddd_6 = cell((s_6(1)/8),(s_6(2)/8));
CD_2_real_down_cell_new = cell((s_6(1)/8),(s_6(2)/8));

vvv_7 = cell((s_7(1)/8),(s_7(2)/8));
sss_7 = cell((s_7(1)/8),(s_7(2)/8));
ddd_7 = cell((s_7(1)/8),(s_7(2)/8));
CH_1_imag_up_cell_new = cell((s_7(1)/8),(s_7(2)/8));

vvv_8 = cell((s_8(1)/8),(s_8(2)/8));
sss_8 = cell((s_8(1)/8),(s_8(2)/8));
ddd_8 = cell((s_8(1)/8),(s_8(2)/8));
CV_1_imag_up_cell_new = cell((s_8(1)/8),(s_8(2)/8));

vvv_9 = cell((s_9(1)/8),(s_9(2)/8));
sss_9 = cell((s_9(1)/8),(s_9(2)/8));
ddd_9 = cell((s_9(1)/8),(s_9(2)/8));
CD_1_imag_up_cell_new = cell((s_9(1)/8),(s_9(2)/8));

vvv_10 = cell((s_10(1)/8),(s_10(2)/8));
sss_10 = cell((s_10(1)/8),(s_10(2)/8));
ddd_10 = cell((s_10(1)/8),(s_10(2)/8));
CH_2_imag_down_cell_new = cell((s_10(1)/8),(s_10(2)/8));

vvv_11 = cell((s_11(1)/8),(s_11(2)/8));
sss_11 = cell((s_11(1)/8),(s_11(2)/8));
ddd_11 = cell((s_11(1)/8),(s_11(2)/8));
CV_2_imag_down_cell_new = cell((s_11(1)/8),(s_11(2)/8));

vvv_12 = cell((s_12(1)/8),(s_12(2)/8));
sss_12 = cell((s_12(1)/8),(s_12(2)/8));
ddd_12 = cell((s_12(1)/8),(s_12(2)/8));
CD_2_imag_down_cell_new = cell((s_12(1)/8),(s_12(2)/8));

% real detail subband denoising;
T_1 = zeros(8,8);
T_2 = zeros(8,8);
T_3 = zeros(8,8);
T_4 = zeros(1,8);
K = zeros(8,8);
beta = 7;
for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_2 = ss_1{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_3 = transpose(dd_1{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
        for j=1:s_1(2)/8
        CH_1_real_up_cell_new{i,j}=sss_1{i,j}*vvv_1{i,j}*ddd_1{i,j};
        end
    end
    CH_1_real = cell2mat(CH_1_real_up_cell_new);
    CH_1_real_up_new = CH_1_real(1:p_1(1),1:p_1(2));
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_2 = ss_2{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_3 = transpose(dd_2{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        CV_1_real_up_cell_new{i,j}=sss_2{i,j}*vvv_2{i,j}*ddd_2{i,j};
    end
    end
    CV_1_real = cell2mat(CV_1_real_up_cell_new); 
    CV_1_real_up_new = CV_1_real(1:p_2(1),1:p_2(2));
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_2 = ss_3{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_3 = transpose(dd_3{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        CD_1_real_up_cell_new{i,j}=sss_3{i,j}*vvv_3{i,j}*ddd_3{i,j};
    end
    end
    CD_1_real = cell2mat(CD_1_real_up_cell_new);
    CD_1_real_up_new = CD_1_real(1:p_3(1),1:p_3(2));
for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_2 = ss_4{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_3 = transpose(dd_4{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        CH_2_real_down_cell_new{i,j}=sss_4{i,j}*vvv_4{i,j}*ddd_4{i,j};
    end
    end
    CH_2_real = cell2mat(CH_2_real_down_cell_new);
    CH_2_real_down_new = CH_2_real(1:p_4(1),1:p_4(2));
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_2 = ss_5{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_3 = transpose(dd_5{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        CV_2_real_down_cell_new{i,j}=sss_5{i,j}*vvv_5{i,j}*ddd_5{i,j};
    end
    end
    CV_2_real = cell2mat(CV_2_real_down_cell_new); 
    CV_2_real_down_new = CV_2_real(1:p_5(1),1:p_5(2));
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_2 = ss_6{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_3 = transpose(dd_6{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        CD_2_real_down_cell_new{i,j}=sss_6{i,j}*vvv_6{i,j}*ddd_6{i,j};
    end
    end
    CD_2_real = cell2mat(CD_2_real_down_cell_new);
    CD_2_real_down_new = CD_2_real(1:p_6(1),1:p_6(2));
% imaginary detail subband denoising; 
for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_2 = ss_7{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_3 = transpose(dd_7{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
        for j=1:s_7(2)/8
        CH_1_imag_up_cell_new{i,j}=sss_7{i,j}*vvv_7{i,j}*ddd_7{i,j};
        end
    end
    CH_1_imag = cell2mat(CH_1_imag_up_cell_new);
    CH_1_imag_up_new = CH_1_imag(1:p_7(1),1:p_7(2));
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_2 = ss_8{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_3 = transpose(dd_8{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        CV_1_imag_up_cell_new{i,j}=sss_8{i,j}*vvv_8{i,j}*ddd_8{i,j};
    end
    end
    CV_1_imag = cell2mat(CV_1_imag_up_cell_new); 
    CV_1_imag_up_new = CV_1_imag(1:p_8(1),1:p_8(2));
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_2 = ss_9{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_3 = transpose(dd_9{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        CD_1_imag_up_cell_new{i,j}=sss_9{i,j}*vvv_9{i,j}*ddd_9{i,j};
    end
    end
    CD_1_imag = cell2mat(CD_1_imag_up_cell_new);
    CD_1_imag_up_new = CD_1_imag(1:p_9(1),1:p_9(2));
for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_2 = ss_10{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_3 = transpose(dd_10{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        CH_2_imag_down_cell_new{i,j}=sss_10{i,j}*vvv_10{i,j}*ddd_10{i,j};
    end
    end
    CH_2_imag = cell2mat(CH_2_imag_down_cell_new);
    CH_2_imag_down_new = CH_2_imag(1:p_10(1),1:p_10(2));
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_2 = ss_11{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_3 = transpose(dd_11{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        CV_2_imag_down_cell_new{i,j}=sss_11{i,j}*vvv_11{i,j}*ddd_11{i,j};
    end
    end
    CV_2_imag = cell2mat(CV_2_imag_down_cell_new); 
    CV_2_imag_down_new = CV_2_imag(1:p_11(1),1:p_11(2));
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_2 = ss_12{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_3 = transpose(dd_12{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        CD_2_imag_down_cell_new{i,j}=sss_12{i,j}*vvv_12{i,j}*ddd_12{i,j};
    end
    end
    CD_2_imag = cell2mat(CD_2_imag_down_cell_new);
    CD_2_imag_down_new = CD_2_imag(1:p_12(1),1:p_12(2));
% create the new array w
w{1}{1}{1}{1} = CH_1_real_up_new;
w{1}{1}{1}{2} = CV_1_real_up_new;
w{1}{1}{1}{3} = CD_1_real_up_new;
w{1}{1}{2}{1} = CH_2_real_down_new;
w{1}{1}{2}{2} = CV_2_real_down_new;
w{1}{1}{2}{3} = CD_2_real_down_new;
w{1}{2}{1}{1} = CH_1_imag_up_new;
w{1}{2}{1}{2} = CV_1_imag_up_new;
w{1}{2}{1}{3} = CD_1_imag_up_new;
w{1}{2}{2}{1} = CH_2_imag_down_new;
w{1}{2}{2}{2} = CV_2_imag_down_new;
w{1}{2}{2}{3} = CD_2_imag_down_new;
y = icplxdual2D(w, 1, Fsf, sf);
% calculate MSE for the denoised image;
MSE_1 = MSE(B_0,y);
% calculate MSE for the nosiy image;    
MSE_2 = MSE(B_0,B_1);    
M = [MSE_1 MSE_2];
subplot(1,3,1);
imagesc(B_0);
colormap gray;
title('original image');
subplot(1,3,2);
imagesc(B_1);
colormap gray;
title('noisy image');
subplot(1,3,3);
imagesc(y);
colormap gray;
title('denoised image');
% denoised observation
snr_1 = calSNR(B_0,y);
% noisy observation
snr_2 = calSNR(B_0,B_1);
snr = [snr_1 snr_2];

%%
% dual-tree complex wavelet transform and adaptive svd image despeckling;
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
B_1 = speck(B_0,0.5,3,3);
B_1 = log(B_1);
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
w = cplxdual2D(B_1, 1, Faf, af);

CH_1_real_up = w{1}{1}{1}{1};
CV_1_real_up = w{1}{1}{1}{2};
CD_1_real_up = w{1}{1}{1}{3};
CH_2_real_down = w{1}{1}{2}{1};
CV_2_real_down = w{1}{1}{2}{2};
CD_2_real_down = w{1}{1}{2}{3};
CH_1_imag_up = w{1}{2}{1}{1};
CV_1_imag_up = w{1}{2}{1}{2};
CD_1_imag_up = w{1}{2}{1}{3};
CH_2_imag_down = w{1}{2}{2}{1};
CV_2_imag_down = w{1}{2}{2}{2};
CD_2_imag_down = w{1}{2}{2}{3};

p_1 = size(CH_1_real_up);
p_2 = size(CV_1_real_up);
p_3 = size(CD_1_real_up);
p_4 = size(CH_2_real_down);
p_5 = size(CV_2_real_down);
p_6 = size(CD_2_real_down);
p_7 = size(CH_1_imag_up);
p_8 = size(CV_1_imag_up);
p_9 = size(CD_1_imag_up);
p_10 = size(CH_2_imag_down);
p_11 = size(CV_2_imag_down);
p_12 = size(CD_2_imag_down);

CD_1_real_up_re = reshape(CD_1_real_up,1,p_3(1)*p_3(2));
sigma_1_HH = median(abs(CD_1_real_up_re))/0.6745;

if  mod(p_1(1),8)~= 0
    CH_1_real_up(p_1(1)+1:p_1(1)+8-mod(p_1(1),8),:) = 0;
end

if  mod(p_1(2),8)~= 0
    CH_1_real_up(:,p_1(2)+1:p_1(2)+8-mod(p_1(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_2(1),8)~= 0
    CV_1_real_up(p_2(1)+1:p_2(1)+8-mod(p_2(1),8),:) = 0;
end

if  mod(p_2(2),8)~= 0
    CV_1_real_up(:,p_2(2)+1:p_2(2)+8-mod(p_2(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_3(1),8)~= 0
    CD_1_real_up(p_3(1)+1:p_3(1)+8-mod(p_3(1),8),:) = 0;
end

if  mod(p_3(2),8)~= 0
    CD_1_real_up(:,p_3(2)+1:p_3(2)+8-mod(p_3(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_4(1),8)~= 0
    CH_2_real_down(p_4(1)+1:p_4(1)+8-mod(p_4(1),8),:) = 0;
end

if  mod(p_4(2),8)~= 0
    CH_2_real_down(:,p_4(2)+1:p_4(2)+8-mod(p_4(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_5(1),8)~= 0
    CV_2_real_down(p_5(1)+1:p_5(1)+8-mod(p_5(1),8),:) = 0;
end

if  mod(p_5(2),8)~= 0
    CV_2_real_down(:,p_5(2)+1:p_5(2)+8-mod(p_5(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_6(1),8)~= 0
    CD_2_real_down(p_6(1)+1:p_6(1)+8-mod(p_6(1),8),:) = 0;
end

if  mod(p_6(2),8)~= 0
    CD_2_real_down(:,p_6(2)+1:p_6(2)+8-mod(p_6(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_7(1),8)~= 0
    CH_1_imag_up(p_7(1)+1:p_7(1)+8-mod(p_7(1),8),:) = 0;
end

if  mod(p_7(2),8)~= 0
    CH_1_imag_up(:,p_7(2)+1:p_7(2)+8-mod(p_7(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_8(1),8)~= 0
    CV_1_imag_up(p_8(1)+1:p_8(1)+8-mod(p_8(1),8),:) = 0;
end

if  mod(p_8(2),8)~= 0
    CV_1_imag_up(:,p_8(2)+1:p_8(2)+8-mod(p_8(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_9(1),8)~= 0
    CD_1_imag_up(p_9(1)+1:p_9(1)+8-mod(p_9(1),8),:) = 0;
end

if  mod(p_9(2),8)~= 0
    CD_1_imag_up(:,p_9(2)+1:p_9(2)+8-mod(p_9(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_10(1),8)~= 0
    CH_2_imag_down(p_10(1)+1:p_10(1)+8-mod(p_10(1),8),:) = 0;
end

if  mod(p_10(2),8)~= 0
    CH_2_imag_down(:,p_10(2)+1:p_10(2)+8-mod(p_10(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_11(1),8)~= 0
    CV_2_imag_down(p_11(1)+1:p_11(1)+8-mod(p_11(1),8),:) = 0;
end

if  mod(p_11(2),8)~= 0
    CV_2_imag_down(:,p_11(2)+1:p_11(2)+8-mod(p_11(2),8)) = 0;
end
%%%%%%%%%%
if  mod(p_12(1),8)~= 0
    CD_2_imag_down(p_12(1)+1:p_12(1)+8-mod(p_12(1),8),:) = 0;
end

if  mod(p_12(2),8)~= 0
    CD_2_imag_down(:,p_12(2)+1:p_12(2)+8-mod(p_12(2),8)) = 0;
end

s_1 = size(CH_1_real_up);
s_2 = size(CV_1_real_up);
s_3 = size(CD_1_real_up);
s_4 = size(CH_2_real_down);
s_5 = size(CV_2_real_down);
s_6 = size(CD_2_real_down);
s_7 = size(CH_1_imag_up);
s_8 = size(CV_1_imag_up);
s_9 = size(CD_1_imag_up);
s_10 = size(CH_2_imag_down);
s_11 = size(CV_2_imag_down);
s_12 = size(CD_2_imag_down);

%making 8*8 Blocks
CH_1_real_up_cell = cell((s_1(1)/8),(s_1(2)/8));
for i=0:s_1(1)/8-1
    for j=0:s_1(2)/8-1
        CH_1_real_up_cell{i+1,j+1}=CH_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_1_real_up_cell = cell((s_2(1)/8),(s_2(2)/8));
for i=0:s_2(1)/8-1
    for j=0:s_2(2)/8-1
        CV_1_real_up_cell{i+1,j+1}=CV_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_1_real_up_cell = cell((s_3(1)/8),(s_3(2)/8));
for i=0:s_3(1)/8-1
    for j=0:s_3(2)/8-1
        CD_1_real_up_cell{i+1,j+1}=CD_1_real_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_2_real_down_cell = cell((s_4(1)/8),(s_4(2)/8));
for i=0:s_4(1)/8-1
    for j=0:s_4(2)/8-1
        CH_2_real_down_cell{i+1,j+1}=CH_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_2_real_down_cell = cell((s_5(1)/8),(s_5(2)/8));
for i=0:s_5(1)/8-1
    for j=0:s_5(2)/8-1
        CV_2_real_down_cell{i+1,j+1}=CV_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_2_real_down_cell = cell((s_6(1)/8),(s_6(2)/8));
for i=0:s_6(1)/8-1
    for j=0:s_6(2)/8-1
        CD_2_real_down_cell{i+1,j+1}=CD_2_real_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_1_imag_up_cell = cell((s_7(1)/8),(s_7(2)/8));
for i=0:s_7(1)/8-1
    for j=0:s_7(2)/8-1
        CH_1_imag_up_cell{i+1,j+1}=CH_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_1_imag_up_cell = cell((s_8(1)/8),(s_8(2)/8));
for i=0:s_8(1)/8-1
    for j=0:s_8(2)/8-1
        CV_1_imag_up_cell{i+1,j+1}=CV_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_1_imag_up_cell = cell((s_9(1)/8),(s_9(2)/8));
for i=0:s_9(1)/8-1
    for j=0:s_9(2)/8-1
        CD_1_imag_up_cell{i+1,j+1}=CD_1_imag_up(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CH_2_imag_down_cell = cell((s_10(1)/8),(s_10(2)/8));
for i=0:s_10(1)/8-1
    for j=0:s_10(2)/8-1
        CH_2_imag_down_cell{i+1,j+1}=CH_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CV_2_imag_down_cell = cell((s_11(1)/8),(s_11(2)/8));
for i=0:s_11(1)/8-1
    for j=0:s_11(2)/8-1
        CV_2_imag_down_cell{i+1,j+1}=CV_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

CD_2_imag_down_cell = cell((s_12(1)/8),(s_12(2)/8));
for i=0:s_12(1)/8-1
    for j=0:s_12(2)/8-1
        CD_2_imag_down_cell{i+1,j+1}=CD_2_imag_down(8*i+1:8*(i+1),8*j+1:8*(j+1));
    end
end

ss_1 = cell((s_1(1)/8),(s_1(2)/8));
vv_1 = cell((s_1(1)/8),(s_1(2)/8));
dd_1 = cell((s_1(1)/8),(s_1(2)/8));

ss_2 = cell((s_2(1)/8),(s_2(2)/8));
vv_2 = cell((s_2(1)/8),(s_2(2)/8));
dd_2 = cell((s_2(1)/8),(s_2(2)/8));

ss_3 = cell((s_3(1)/8),(s_3(2)/8));
vv_3 = cell((s_3(1)/8),(s_3(2)/8));
dd_3 = cell((s_3(1)/8),(s_3(2)/8));

ss_4 = cell((s_4(1)/8),(s_4(2)/8));
vv_4 = cell((s_4(1)/8),(s_4(2)/8));
dd_4 = cell((s_4(1)/8),(s_4(2)/8));

ss_5 = cell((s_5(1)/8),(s_5(2)/8));
vv_5 = cell((s_5(1)/8),(s_5(2)/8));
dd_5 = cell((s_5(1)/8),(s_5(2)/8));

ss_6 = cell((s_6(1)/8),(s_6(2)/8));
vv_6 = cell((s_6(1)/8),(s_6(2)/8));
dd_6 = cell((s_6(1)/8),(s_6(2)/8));

ss_7 = cell((s_7(1)/8),(s_7(2)/8));
vv_7 = cell((s_7(1)/8),(s_7(2)/8));
dd_7 = cell((s_7(1)/8),(s_7(2)/8));

ss_8 = cell((s_8(1)/8),(s_8(2)/8));
vv_8 = cell((s_8(1)/8),(s_8(2)/8));
dd_8 = cell((s_8(1)/8),(s_8(2)/8));

ss_9 = cell((s_9(1)/8),(s_9(2)/8));
vv_9 = cell((s_9(1)/8),(s_9(2)/8));
dd_9 = cell((s_9(1)/8),(s_9(2)/8));

ss_10 = cell((s_10(1)/8),(s_10(2)/8));
vv_10 = cell((s_10(1)/8),(s_10(2)/8));
dd_10 = cell((s_10(1)/8),(s_10(2)/8));

ss_11 = cell((s_11(1)/8),(s_11(2)/8));
vv_11 = cell((s_11(1)/8),(s_11(2)/8));
dd_11 = cell((s_11(1)/8),(s_11(2)/8));

ss_12 = cell((s_12(1)/8),(s_12(2)/8));
vv_12 = cell((s_12(1)/8),(s_12(2)/8));
dd_12 = cell((s_12(1)/8),(s_12(2)/8));

for i=1:s_1(1)/8
    for j=1:s_1(2)/8
    [ss_1{i,j} vv_1{i,j} dd_1{i,j}] = svd(CH_1_real_up_cell{i,j});
    end
end

for i=1:s_2(1)/8
    for j=1:s_2(2)/8
    [ss_2{i,j} vv_2{i,j} dd_2{i,j}] = svd(CV_1_real_up_cell{i,j});
    end
end

for i=1:s_3(1)/8
    for j=1:s_3(2)/8
    [ss_3{i,j} vv_3{i,j} dd_3{i,j}] = svd(CD_1_real_up_cell{i,j});
    end
end

for i=1:s_4(1)/8
    for j=1:s_4(2)/8
    [ss_4{i,j} vv_4{i,j} dd_4{i,j}] = svd(CH_2_real_down_cell{i,j});
    end
end

for i=1:s_5(1)/8
    for j=1:s_5(2)/8
    [ss_5{i,j} vv_5{i,j} dd_5{i,j}] = svd(CV_2_real_down_cell{i,j});
    end
end

for i=1:s_6(1)/8
    for j=1:s_6(2)/8
    [ss_6{i,j} vv_6{i,j} dd_6{i,j}] = svd(CD_2_real_down_cell{i,j});
    end
end

for i=1:s_7(1)/8
    for j=1:s_7(2)/8
    [ss_7{i,j} vv_7{i,j} dd_7{i,j}] = svd(CH_1_imag_up_cell{i,j});
    end
end

for i=1:s_8(1)/8
    for j=1:s_8(2)/8
    [ss_8{i,j} vv_8{i,j} dd_8{i,j}] = svd(CV_1_imag_up_cell{i,j});
    end
end

for i=1:s_9(1)/8
    for j=1:s_9(2)/8
    [ss_9{i,j} vv_9{i,j} dd_9{i,j}] = svd(CD_1_imag_up_cell{i,j});
    end
end

for i=1:s_10(1)/8
    for j=1:s_10(2)/8
    [ss_10{i,j} vv_10{i,j} dd_10{i,j}] = svd(CH_2_imag_down_cell{i,j});
    end
end

for i=1:s_11(1)/8
    for j=1:s_11(2)/8
    [ss_11{i,j} vv_11{i,j} dd_11{i,j}] = svd(CV_2_imag_down_cell{i,j});
    end
end

for i=1:s_12(1)/8
    for j=1:s_12(2)/8
    [ss_12{i,j} vv_12{i,j} dd_12{i,j}] = svd(CD_2_imag_down_cell{i,j});
    end
end

%changing svd
vvv_1 = cell((s_1(1)/8),(s_1(2)/8));
sss_1 = cell((s_1(1)/8),(s_1(2)/8));
ddd_1 = cell((s_1(1)/8),(s_1(2)/8));
CH_1_real_up_cell_new = cell((s_1(1)/8),(s_1(2)/8));

vvv_2 = cell((s_2(1)/8),(s_2(2)/8));
sss_2 = cell((s_2(1)/8),(s_2(2)/8));
ddd_2 = cell((s_2(1)/8),(s_2(2)/8));
CV_1_real_up_cell_new = cell((s_2(1)/8),(s_2(2)/8));

vvv_3 = cell((s_3(1)/8),(s_3(2)/8));
sss_3 = cell((s_3(1)/8),(s_3(2)/8));
ddd_3 = cell((s_3(1)/8),(s_3(2)/8));
CD_1_real_up_cell_new = cell((s_3(1)/8),(s_3(2)/8));

vvv_4 = cell((s_4(1)/8),(s_4(2)/8));
sss_4 = cell((s_4(1)/8),(s_4(2)/8));
ddd_4 = cell((s_4(1)/8),(s_4(2)/8));
CH_2_real_down_cell_new = cell((s_4(1)/8),(s_4(2)/8));

vvv_5 = cell((s_5(1)/8),(s_5(2)/8));
sss_5 = cell((s_5(1)/8),(s_5(2)/8));
ddd_5 = cell((s_5(1)/8),(s_5(2)/8));
CV_2_real_down_cell_new = cell((s_5(1)/8),(s_5(2)/8));

vvv_6 = cell((s_6(1)/8),(s_6(2)/8));
sss_6 = cell((s_6(1)/8),(s_6(2)/8));
ddd_6 = cell((s_6(1)/8),(s_6(2)/8));
CD_2_real_down_cell_new = cell((s_6(1)/8),(s_6(2)/8));

vvv_7 = cell((s_7(1)/8),(s_7(2)/8));
sss_7 = cell((s_7(1)/8),(s_7(2)/8));
ddd_7 = cell((s_7(1)/8),(s_7(2)/8));
CH_1_imag_up_cell_new = cell((s_7(1)/8),(s_7(2)/8));

vvv_8 = cell((s_8(1)/8),(s_8(2)/8));
sss_8 = cell((s_8(1)/8),(s_8(2)/8));
ddd_8 = cell((s_8(1)/8),(s_8(2)/8));
CV_1_imag_up_cell_new = cell((s_8(1)/8),(s_8(2)/8));

vvv_9 = cell((s_9(1)/8),(s_9(2)/8));
sss_9 = cell((s_9(1)/8),(s_9(2)/8));
ddd_9 = cell((s_9(1)/8),(s_9(2)/8));
CD_1_imag_up_cell_new = cell((s_9(1)/8),(s_9(2)/8));

vvv_10 = cell((s_10(1)/8),(s_10(2)/8));
sss_10 = cell((s_10(1)/8),(s_10(2)/8));
ddd_10 = cell((s_10(1)/8),(s_10(2)/8));
CH_2_imag_down_cell_new = cell((s_10(1)/8),(s_10(2)/8));

vvv_11 = cell((s_11(1)/8),(s_11(2)/8));
sss_11 = cell((s_11(1)/8),(s_11(2)/8));
ddd_11 = cell((s_11(1)/8),(s_11(2)/8));
CV_2_imag_down_cell_new = cell((s_11(1)/8),(s_11(2)/8));

vvv_12 = cell((s_12(1)/8),(s_12(2)/8));
sss_12 = cell((s_12(1)/8),(s_12(2)/8));
ddd_12 = cell((s_12(1)/8),(s_12(2)/8));
CD_2_imag_down_cell_new = cell((s_12(1)/8),(s_12(2)/8));

% real detail subband denoising;
T_1 = zeros(8,8);
T_2 = zeros(8,8);
T_3 = zeros(8,8);
T_4 = zeros(1,8);
K = zeros(8,8);
beta = input('insert the value of beta = ');
for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_2 = ss_1{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
    for j=1:s_1(2)/8
        T_1 = vv_1{i,j};
        T_3 = transpose(dd_1{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_1{i,j}= K;
    end
    end
    for i=1:s_1(1)/8
        for j=1:s_1(2)/8
        CH_1_real_up_cell_new{i,j}=sss_1{i,j}*vvv_1{i,j}*ddd_1{i,j};
        end
    end
    CH_1_real = cell2mat(CH_1_real_up_cell_new);
    CH_1_real_up_new = CH_1_real(1:p_1(1),1:p_1(2));
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_2 = ss_2{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        T_1 = vv_2{i,j};
        T_3 = transpose(dd_2{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_2{i,j}= K;
    end
    end
    for i=1:s_2(1)/8
    for j=1:s_2(2)/8
        CV_1_real_up_cell_new{i,j}=sss_2{i,j}*vvv_2{i,j}*ddd_2{i,j};
    end
    end
    CV_1_real = cell2mat(CV_1_real_up_cell_new); 
    CV_1_real_up_new = CV_1_real(1:p_2(1),1:p_2(2));
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_2 = ss_3{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        T_1 = vv_3{i,j};
        T_3 = transpose(dd_3{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_3{i,j}= K;
    end
    end
    for i=1:s_3(1)/8
    for j=1:s_3(2)/8
        CD_1_real_up_cell_new{i,j}=sss_3{i,j}*vvv_3{i,j}*ddd_3{i,j};
    end
    end
    CD_1_real = cell2mat(CD_1_real_up_cell_new);
    CD_1_real_up_new = CD_1_real(1:p_3(1),1:p_3(2));
for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_2 = ss_4{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        T_1 = vv_4{i,j};
        T_3 = transpose(dd_4{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_4{i,j}= K;
    end
    end
    for i=1:s_4(1)/8
    for j=1:s_4(2)/8
        CH_2_real_down_cell_new{i,j}=sss_4{i,j}*vvv_4{i,j}*ddd_4{i,j};
    end
    end
    CH_2_real = cell2mat(CH_2_real_down_cell_new);
    CH_2_real_down_new = CH_2_real(1:p_4(1),1:p_4(2));
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_2 = ss_5{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        T_1 = vv_5{i,j};
        T_3 = transpose(dd_5{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_5{i,j}= K;
    end
    end
    for i=1:s_5(1)/8
    for j=1:s_5(2)/8
        CV_2_real_down_cell_new{i,j}=sss_5{i,j}*vvv_5{i,j}*ddd_5{i,j};
    end
    end
    CV_2_real = cell2mat(CV_2_real_down_cell_new); 
    CV_2_real_down_new = CV_2_real(1:p_5(1),1:p_5(2));
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_2 = ss_6{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        T_1 = vv_6{i,j};
        T_3 = transpose(dd_6{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_6{i,j}= K;
    end
    end
    for i=1:s_6(1)/8
    for j=1:s_6(2)/8
        CD_2_real_down_cell_new{i,j}=sss_6{i,j}*vvv_6{i,j}*ddd_6{i,j};
    end
    end
    CD_2_real = cell2mat(CD_2_real_down_cell_new);
    CD_2_real_down_new = CD_2_real(1:p_6(1),1:p_6(2));
% imaginary detail subband denoising; 
for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_2 = ss_7{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
    for j=1:s_7(2)/8
        T_1 = vv_7{i,j};
        T_3 = transpose(dd_7{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_7{i,j}= K;
    end
    end
    for i=1:s_7(1)/8
        for j=1:s_7(2)/8
        CH_1_imag_up_cell_new{i,j}=sss_7{i,j}*vvv_7{i,j}*ddd_7{i,j};
        end
    end
    CH_1_imag = cell2mat(CH_1_imag_up_cell_new);
    CH_1_imag_up_new = CH_1_imag(1:p_7(1),1:p_7(2));
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_2 = ss_8{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        T_1 = vv_8{i,j};
        T_3 = transpose(dd_8{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_8{i,j}= K;
    end
    end
    for i=1:s_8(1)/8
    for j=1:s_8(2)/8
        CV_1_imag_up_cell_new{i,j}=sss_8{i,j}*vvv_8{i,j}*ddd_8{i,j};
    end
    end
    CV_1_imag = cell2mat(CV_1_imag_up_cell_new); 
    CV_1_imag_up_new = CV_1_imag(1:p_8(1),1:p_8(2));
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_2 = ss_9{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        T_1 = vv_9{i,j};
        T_3 = transpose(dd_9{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_9{i,j}= K;
    end
    end
    for i=1:s_9(1)/8
    for j=1:s_9(2)/8
        CD_1_imag_up_cell_new{i,j}=sss_9{i,j}*vvv_9{i,j}*ddd_9{i,j};
    end
    end
    CD_1_imag = cell2mat(CD_1_imag_up_cell_new);
    CD_1_imag_up_new = CD_1_imag(1:p_9(1),1:p_9(2));
for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_2 = ss_10{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        T_1 = vv_10{i,j};
        T_3 = transpose(dd_10{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_10{i,j}= K;
    end
    end
    for i=1:s_10(1)/8
    for j=1:s_10(2)/8
        CH_2_imag_down_cell_new{i,j}=sss_10{i,j}*vvv_10{i,j}*ddd_10{i,j};
    end
    end
    CH_2_imag = cell2mat(CH_2_imag_down_cell_new);
    CH_2_imag_down_new = CH_2_imag(1:p_10(1),1:p_10(2));
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_2 = ss_11{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        T_1 = vv_11{i,j};
        T_3 = transpose(dd_11{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_11{i,j}= K;
    end
    end
    for i=1:s_11(1)/8
    for j=1:s_11(2)/8
        CV_2_imag_down_cell_new{i,j}=sss_11{i,j}*vvv_11{i,j}*ddd_11{i,j};
    end
    end
    CV_2_imag = cell2mat(CV_2_imag_down_cell_new); 
    CV_2_imag_down_new = CV_2_imag(1:p_11(1),1:p_11(2));
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_1(1:e,1:e);
        elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_1(1:e,1:e);
                    break
                elseif a==1
                    K = T_1(1:a,1:a);
                end
            end
        end
        vvv_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_2 = ss_12{i,j}; 
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_2(1:8,1:e);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_2(1:8,1:e);
                    break
                elseif a==1
                    K = T_2(1:8,1:a);
                end
            end
        end
        sss_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        T_1 = vv_12{i,j};
        T_3 = transpose(dd_12{i,j});
        T_4 = diag(T_1);
        if median(T_4)>beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>sigma_1_HH
                    e = a;
                    break
                end
            end
            K = T_3(1:e,1:8);
            elseif median(T_4)<=beta*sigma_1_HH
            for a=8:-1:1
                if T_4(a)>beta*sigma_1_HH
                    e = a;
                    K = T_3(1:e,1:8);
                    break
                elseif a==1
                    K = T_3(1:a,1:8);
                end
            end
        end
        ddd_12{i,j}= K;
    end
    end
    for i=1:s_12(1)/8
    for j=1:s_12(2)/8
        CD_2_imag_down_cell_new{i,j}=sss_12{i,j}*vvv_12{i,j}*ddd_12{i,j};
    end
    end
    CD_2_imag = cell2mat(CD_2_imag_down_cell_new);
    CD_2_imag_down_new = CD_2_imag(1:p_12(1),1:p_12(2));
% create the new array w;
w{1}{1}{1}{1} = CH_1_real_up_new;
w{1}{1}{1}{2} = CV_1_real_up_new;
w{1}{1}{1}{3} = CD_1_real_up_new;
w{1}{1}{2}{1} = CH_2_real_down_new;
w{1}{1}{2}{2} = CV_2_real_down_new;
w{1}{1}{2}{3} = CD_2_real_down_new;
w{1}{2}{1}{1} = CH_1_imag_up_new;
w{1}{2}{1}{2} = CV_1_imag_up_new;
w{1}{2}{1}{3} = CD_1_imag_up_new;
w{1}{2}{2}{1} = CH_2_imag_down_new;
w{1}{2}{2}{2} = CV_2_imag_down_new;
w{1}{2}{2}{3} = CD_2_imag_down_new;
y = icplxdual2D(w, 1, Fsf, sf);
y = exp(y);
y = round(y);
B_1 = exp(B_1);
% calculate MSE for the denoised image;
MSE_1=0;
    for i=1:l(1)
    for j=1:l(2)
        MSE_1 = MSE_1+((y(i,j)-B_0(i,j))^2)/(l(1)*l(2));
    end
    end
% calculate MSE for the nosiy image;    
MSE_2=0;
    for i=1:l(1)
    for j=1:l(2)
        MSE_2 = MSE_2+((B_1(i,j)-B_0(i,j))^2)/(l(1)*l(2));
    end
    end
    
M = [MSE_1 MSE_2];
subplot(1,3,1);
imagesc(B_0);
colormap gray;
title('original image');
subplot(1,3,2);
imagesc(B_1);
colormap gray;
title('noisy image');
subplot(1,3,3);
imagesc(y);
colormap gray;
title('denoised image');
snr_1 = calSNR(B_0,y);
snr_2 = calSNR(B_0,B_1);
snr = [snr_1 snr_2];


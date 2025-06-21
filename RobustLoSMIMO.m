clear all 
close all
clc

W_side_vect = [2];% [0,1,2,3];
Nwside = length(W_side_vect);
redund_vect = [2];%[1,2,3,4]; % how many extra 4 receivers we have
Nredund = length(redund_vect);
SISO_SNR_dB_vect =[[],-5:1:35];%[-6:3:18,inf];
Nsnr = length(SISO_SNR_dB_vect);
Nrun = 150; % # randomized runs for each setting

SINR_MonteCarlo = zeros(Nredund,Nwside,Nsnr);
BER_MonteCarlo = zeros(Nredund,Nwside,Nsnr);
SINR_analytical = zeros(Nredund,Nwside,Nsnr);
BER_analytical = zeros(Nredund,Nwside,Nsnr);
ratio_dB_EVM = zeros(Nredund,Nwside,Nsnr); % computed using EVM
ratio_lin_anal = zeros(Nredund,Nwside,Nsnr); % computed analytically

correct_flag = 1;

lambda = 0.0023;
speed_of_light = 3e8*1e-12;
fc = speed_of_light/lambda;
T_symb = 50; %pico seconds
max_antenna_time_offset = T_symb/8;
pulse_sidewidth = 5;
L_pulse = 2*pulse_sidewidth + 1;
Nsymb = 51000 + 50;
Ntraining_symb = length(10*(-2 - pulse_sidewidth - 11)*1 : 1*(11 + 2 + pulse_sidewidth) )+ length((-2 - pulse_sidewidth - 11)*1 : 1*(11 + 2 + pulse_sidewidth)) -1  + 50; 
backoff = 0.25;


R_vect = 50;%[[],0:1:200];%100
%R_vect = [[],0:1:200];
R_len = length(R_vect);
Czf_R = zeros(R_len,Nrun);
mean_ZF = zeros(R_len,Nrun);
H_corr12or34 = zeros(R_len,Nrun);
H_corr13or24 = zeros(R_len,Nrun);
H_corr23or14 = zeros(R_len,Nrun);
Czf_R_avg = zeros(R_len,1);
mean_ZF_avg = zeros(R_len,1);
H_corr12or34_avg  = zeros(R_len,1);
H_corr13or24_avg  = zeros(R_len,1);
H_corr23or14_avg  = zeros(R_len,1);
% MIMO link distance and apertures
%R  = 50;  % m
d_tx  = 0.34; % m
txloc = ([0,1; 1,1; 1,0; 0,0]-0.5)*d_tx; % TX locations
txloc(:,3) = 0;
%d_rx_list_0 = d_tx./[1,2,1,2];%abcd %Can add more Rx here (New addition for every 4 Rx elements)
d_rx_list_0 = d_tx./[1,2,2,1];%abdc
d_rx = max(d_rx_list_0);
% Misalignment angles
% vtilt_tx = 0.0499;%((rand()-0.5)*15*1 + 0*7.5)*pi/180;
% vtilt_rx = -0.0708;%((rand()-0.5)*15*1 - 0*4.5)*pi/180;
% htilt_tx =0.0506;%((rand()-0.5)*15*1 + 0*8)*pi/180;
% htilt_rx =-0.0629;%((rand()-0.5)*15*1 + 0*6)*pi/180;

allow_randomness = 1;
Xtx_fixed = exp(1i*pi/4 + 1i*pi/2*round(rand(4,Nsymb)*4));
Xtx_tr = exp(1i*pi/4 + 1i*pi/2*round(rand(4,Nsymb)*4));

numtot = Nredund*Nwside*Nsnr;

for r_val = 1:R_len
    R = R_vect(r_val);
for irun = 1:Nrun
vtilt_tx = ((rand()-0.5)*15*1 + 0*7.5)*pi/180;
vtilt_rx = ((rand()-0.5)*15*1 - 0*4.5)*pi/180;
htilt_tx =  ((rand()-0.5)*15*1 + 0*8)*pi/180;
htilt_rx = ((rand()-0.5)*15*1 + 0*6)*pi/180;
% % 
    % vtilt_tx = 0;
    % vtilt_rx = 0; 
    % htilt_tx =  0;
    % htilt_rx =0; 

    
    for iredund = 1:Nredund
        %Unblock this block of code if u want to check NE for ideal Rx
        %structure and not the og one
%         rxloc = []; %RX location
%         RxFactor = redund_vect(iredund);
%         NrxT = RxFactor*4;
%         for rxi = 1:RxFactor
%             rxloc = [rxloc; d_tx/2,d_tx/2+(rxi-1)*d_tx];
%             rxloc = [rxloc;-d_tx/2, d_tx/2+(rxi-1)*d_tx];
%             rxloc = [rxloc;-d_tx/2,-d_tx/2+(rxi-1)*d_tx];
%             rxloc = [rxloc;d_tx/2,-d_tx/2+(rxi-1)*d_tx];
%          end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%correct placement with abcd in that order of rx grouping (NEW 1)
        d_rx_list = d_rx_list_0(1:redund_vect(iredund));
        rxloc = []; %RX location
% Define the positions of RX antennas
% RX 1 to RX 4 (as per the original positions)
rx_positions = ([0,1; 1,1; 1,0; 0,0]-0.5)*d_tx;

% Calculate the positions for RX5 to RX10 based on midpoints
rx_5 = (rx_positions(1,:) + rx_positions(4,:)) / 2;  % Midpoint of RX1 and RX4
rx_6 = mean(rx_positions, 1);  % Center of the RX array (average of RX1, RX2, RX3, RX4)
rx_7 = (rx_positions(3,:) + rx_positions(4,:)) / 2;  % Midpoint of RX3 and RX4
% rx_8 = (rx_positions(4,:) + rx_5) / 2;  % Midpoint of RX4 and RX5
% rx_9 = (rx_7 + rx_positions(4,:)) / 2;  % Midpoint of RX7 and RX4
% rx_10 = (rx_positions(4,:) + rx_6) / 2;  % Midpoint of RX4 and RX6

% Combine all RX positions into one matrix
all_rx_positions = [
    rx_positions;
    rx_5;
    rx_6;
    rx_7;
    % rx_8;
    % rx_9;
    % rx_10
];
rxloc = all_rx_positions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correct placement with abcd in that order of rx grouping (NEW 2)
%       d_rx_list = d_rx_list_0(1:redund_vect(iredund));
% rxloc = []; % RX location
% 
% % Define the positions of RX1 to RX4 (original positions)
% rx_positions = ([0,1; 1,1; 1,0; 0,0] - 0.5) * d_tx;
% 
% % Calculate the positions for RX5 to RX7 based on midpoints
% rx_5 = (rx_positions(1,:) + rx_positions(4,:)) / 2;  % Midpoint of RX1 and RX4
% rx_6 = mean(rx_positions, 1);  % Center of the RX array (average of RX1, RX2, RX3, RX4)
% rx_7 = (rx_positions(3,:) + rx_positions(4,:)) / 2;  % Midpoint of RX3 and RX4
% 
% % % Calculate RX8 to RX10 positions
% % rx_8 = rx_positions(1,:) + 0.25 * (rx_positions(2,:) - rx_positions(1,:));  % 1/4th distance between RX1 and RX2
% % rx_9 = rx_positions(3,:) + 0.25 * (rx_positions(2,:) - rx_positions(3,:));  % 1/4th distance between RX3 and RX2
% % rx_10 = (rx_6 + rx_positions(2,:)) / 2;  % Midpoint of RX6 (center) and RX2 (diagonal placement)
% 
% % Combine all RX positions into one matrix
% all_rx_positions = [
%     rx_positions;
%     rx_5;
%     rx_6;
%     rx_7;
%     % rx_8;
%     % rx_9;
%     % rx_10
% ];
% 
% rxloc = all_rx_positions;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % d_rx_list = d_rx_list_0(1:redund_vect(iredund));
        % rxloc = []; %RX location
        % for ii1 = 1:length(d_rx_list)
        %     if ii1<=2
        %         rxloc = [rxloc; ([0,1; 1,1; 1,0; 0,0]-0.5)*d_rx_list(ii1)];
        %     end
        % if ii1>2
        %     rxloc = [rxloc; ([0,0.5; 0.5,0; 0,-0.5; -0.5,0])*d_rx_list(ii1)];
        % end
        % 
        % end
% % 
%         % 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % %removing rx10 and rx 11
%         % % Remove RX10 and RX11 (10th and 11th entries)
        % rxloc([10,11], :) = [];  
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%% this portion is for abe setting Nrx = 12
% % Define 1/8 d_rx shift
% shift = d_rx / 8;
% % 
% % Adjust RX9 - RX12 positions
% rxloc(end-3, :) = rxloc(end-3, :) + [-shift, -shift]; % RX9: down, left
% rxloc(end-2, :) = rxloc(end-2, :) + [-shift, shift];  % RX10: left, up
% rxloc(end-1, :) = rxloc(end-1, :) + [shift, shift];   % RX11: up, right
% rxloc(end, :)   = rxloc(end, :)   + [shift, -shift];  % RX12: right, down
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        rxloc(:,3) = R;
        Ntx = size(txloc,1);
        Nrx = size(rxloc,1);
        % Introducing the tilt
        txloc_tilted = tilt_array (txloc, vtilt_tx, htilt_tx);
        rxloc_tilted = tilt_array (rxloc, vtilt_rx, htilt_rx);
        
        % if we want to randomize this will have to do multiple runs for each
        % case
        TX_LO_phase_offsets = (rand(Ntx,1)-0.5)*2*pi*(allow_randomness~=0);
        RX_LO_phase_offsets = (rand(Nrx,1)-0.5)*2*pi*(allow_randomness~=0);
        TX_time_offsets = 0;%(rand(Ntx,1)-0.5)*max_antenna_time_offset*(allow_randomness~=0);
        RX_time_offsets = 0;%(rand(Nrx,1)-0.5)*max_antenna_time_offset*(allow_randomness~=0);
        
        %Calculating differential distance and delay in picoseconds
        dist_tilted = sqrt( (txloc_tilted(:,1).' - rxloc_tilted(:,1)).^2 ...
            + (txloc_tilted(:,2).' - rxloc_tilted(:,2)).^2 ...
            + (txloc_tilted(:,3).' - rxloc_tilted(:,3)).^2 ) ;...
            % Difference in distance calculated in meters
        channel_dist = dist_tilted; % in meters
        channel_delay = (dist_tilted / speed_of_light); %in picoseconds
        channel_delay = channel_delay ...
            + ones(Nrx,1)*TX_time_offsets(:)'...
            + RX_time_offsets(:)*ones(1,Ntx);
        
        channel_phase_response = exp(-1i*(...
            channel_delay*2*pi*fc +...
            ones(Nrx,1)*TX_LO_phase_offsets(:)' +...
            RX_LO_phase_offsets(:)*ones(1,Ntx)...
            ));
        
        channel_delay = channel_delay - min(min(channel_delay));
        
        delays_whole = round(channel_delay/T_symb);
        delays_offset = channel_delay/T_symb - delays_whole;
        max_del_whole = max(delays_whole(:));
        conv_chan = cell(1,Ntx);
        for itx = 1:Ntx
            conv_chan{itx} = raisedcosine_generator(...
                (-pulse_sidewidth:pulse_sidewidth) + delays_offset(:,itx)...
                ,1,backoff)...
                .*channel_phase_response(:,itx);
        end
          
        % QPSK transmit symbols
        Xtx = exp(1i*pi/4 + 1i*pi/2*round(rand(Ntx,Nsymb)*4));
%         Xtx = (allow_randomness~=0)*Xtx + (allow_randomness==0)*Xtx_fixed;
        % received samples
        Yrx = zeros(Nrx,Nsymb+L_pulse-1);
        Yrx_tr = zeros(Nrx,Nsymb+L_pulse-1);
        for itx = 1:Ntx
            resp_mat = conv_chan{itx};
            for irx = 1:Nrx
                xtemp = circshift(Xtx(itx,:),delays_whole(irx,itx)); % don't forget to disregard the first 50 symbols
                Yrx(irx,:) = Yrx(irx,:) + ...
                    conv(xtemp,(fliplr(resp_mat(irx,:))));
                
                xtemp = circshift(Xtx_tr(itx,:),delays_whole(irx,itx)); % don't forget to disregard the first 50 symbols
                Yrx_tr(irx,:) = Yrx_tr(irx,:) + ...
                    conv(xtemp,(fliplr(resp_mat(irx,:))));
                
            end
        end
        Yrx = Yrx(:,pulse_sidewidth+1:end-pulse_sidewidth); % cut off the transient parts
        Yrx_tr = Yrx_tr(:,pulse_sidewidth+1:end-pulse_sidewidth);
        Yrx_noiseless = Yrx;
        Yrx_tr_noiseless = Yrx_tr;
        ylen = size(Yrx,2);
        
        
        % compute the windowed channels
        for iwside = 1:Nwside
            W_side = W_side_vect(iwside);
            W = 1 + 2*W_side; % let's keep the window symmetrical
            
            H_model = cell(Ntx,1);
            Xgrid = (-2 - pulse_sidewidth - 11)*1 : 1*(11 + 2 + pulse_sidewidth);
            Xgridlen = length(Xgrid);
            for itx = 1:Ntx
                H = zeros(Nrx*W , Ntx*Xgridlen);
                for irx = 1:Nrx
                    del_y = delays_whole(irx,itx);
                    for iw = 1:W
                        w = -W_side + iw - 1;
                        for iitx = 1:Ntx
                            del_xy = del_y + w - delays_whole(irx,iitx) - delays_offset(irx,iitx);
                            resp_xy = raisedcosine_generator(Xgrid - del_xy , 1 , backoff)...
                                *(channel_phase_response(irx,iitx));
                            H( (iw-1)*Nrx + irx , (iitx-1)*Xgridlen + (1:Xgridlen) ) = resp_xy;
                        end
                    end
                end
                H_model{itx} = H;
                
            end
             H_est_model  = cell(Ntx,1);               
            %% Training
            % Only for the first run     

                % Defining the Zadoff sequence
                Xtrain_tr = cell(Ntx, 1);
                Xtrain = cell(Ntx, 1);
                desired_ind_tx = zeros(1,Ntx);
                real_ind = zeros(Ntx, Nrx); % What is recieved 
                Htrain = cell(Ntx, Nrx);
                
                grid_center_ind = find(Xgrid == 0);
                for iz = 1:Ntx           
                   z = exp(1i*pi/4 + 1i*pi/2*round(rand(1,Ntraining_symb)*4))';
                    %z = zadoffChuSeq(iz*2,Ntraining_symb);
                    %ztrain = [exp(1i*pi/4 + 1i*pi/2*round(rand(1,50)*4)) z.' exp(1i*pi/4 + 1i*pi/2*round(rand(1,50)*4))];
%                     ztrain = [exp(1i*pi/4 + 1i*pi/2*round(rand(1,50)*4)) z.' z(1:max_del_whole*5).' exp(1i*pi/4 + 1i*pi/2*round(rand(1,50)*4))];
                    Xtrain_tr{iz} = z.';
                    %Xtrain{iz} = ztrain;
                    desired_ind_tx(iz) = (iz - 1)*Xgridlen + grid_center_ind;
                end
                for itx = 1:1 % only Tx1 perspective 
                    H = H_model{itx};
                    for irx = 1:Nrx
                        for iitx = 1:Ntx
                            h1 = H(irx, (iitx - 1)*Xgridlen + 1 : (iitx - 1)*Xgridlen  + Xgridlen);
                            Htrain{iitx,irx} = h1;
                            [M_val, ind1] = max(abs(h1));
                            real_ind(iitx,irx) = (iitx - 1)*Xgridlen + ind1; % Only Tx1 perspective 
                        end
                    end
                end
                % Convolve corresponding sequence [Tx1 perspective still!!]
                Ytrain = cell(Ntx, Nrx);
                H_estimate = cell(Ntx,Nrx);
                y_tx_sum = 0;
                Data_train = cell(Ntx,1);
                D = zeros((Ntraining_symb + Xgridlen - 1) - Xgridlen, Xgridlen);
                for irx =1:Nrx
                    for itx = 1: Ntx
                        h_train = Htrain{itx,irx};
                        x_train = Xtrain_tr{itx};
                        cd =1;
                        for id = Ntraining_symb -Xgridlen  : -1 :1
                            D(id,:) = fliplr(x_train((cd-1)+1:(cd-1)+Xgridlen));
                            cd=cd+1;
                        end
                        noise_lin = 10^(25/10);
                        sigma_sq = 1/noise_lin;
                        y_train = D * h_train.' ;
                        noise_stuff = 1*sqrt(sigma_sq /2)*(randn(size(y_train)) + 1i*randn(size(y_train)));
                        y_train = y_train + noise_stuff;
                        % Applying 1-bit I/Q equal prob Quantizer to the YTrain output
                        
                        Ytrain{itx,irx} = y_train;                   
                        Data_train{itx} = D;

                       h_est = inv(D'*D) * D' *y_train;
                       [M_val, ind] = max(abs(h_est));
                       H_estimate{itx,irx} = circshift(h_est,real_ind(itx,irx) - ind);

                    end
                end
                
                %Forming the Tx1 perspective estimate matrix

                H_tx1 = zeros(Nrx*W, Ntx*Xgridlen);
                H_tx_est = zeros(Nrx*W, Ntx*Xgridlen);
                for irx = 1:Nrx
                    for iw =1:W
                        for itx = 1:Ntx
                            h_est1 = H_estimate{itx,irx};
                            H_tx1( (iw-1)*Nrx + irx , (itx-1)*(Xgridlen) + (1 + (iw-1):Xgridlen + (iw-1)) ) = h_est1;
                            
                        end
                    end
                end
                H_est_model{1} = H_tx1;
                del_txs = cell(Ntx-1,1);
                %Finding first 8 rows for the other 3 matrices 
                for irx =1:Nrx
                    del_source(irx,:) = real_ind(:,irx)' - desired_ind_tx + W_side;
                end
                del_mat = zeros(Nrx,Ntx);
                for itx = 2:Ntx
                    for iirx = 1:Nrx
                        del_tx = del_source(iirx,:) - del_source(iirx,itx) -W_side;
                        del_mat(iirx, :) = del_tx + desired_ind_tx;
                    end
                    del_txs{itx-1} = del_mat;
                end
                for itx = 2:Ntx
                    del_mat = del_txs{itx -1};
                    for irx = 1:Nrx
                        for iw =1:W
                            for iitx = 1:Ntx
                                sym_del = del_mat(irx,iitx) - real_ind(iitx,irx);
                                H_tx1_est = H_tx1( (iw-1)*Nrx + irx, (iitx-1)*(Xgridlen) + (1 + (iw-1):Xgridlen + (iw-1))  ); 
                                H_tx_est( (iw-1)*Nrx + irx , (iitx-1)*(Xgridlen) + (1+(iw-1) :Xgridlen + (iw-1)))...
                                    =circshift(H_tx1_est,sym_del) ;
                            end
                        end
                    end
                    H_est_model{itx} = H_tx_est;
                end
            

        end
    end

            % H1 = H_est_model{1};
            % H2 = H_est_model{2};
            % H3 = H_est_model{3};
            % H4 = H_est_model{4};
            H1 = H_model{1};
            H2 = H_model{2};
            H3 = H_model{3};
            H4 = H_model{4};
            H5 = [H1(:,19) H2(:,56) H3(:,93) H4(:,130)];
            
            H_ind = [[],19,56,93,130];
            
            H_corr12or34(r_val,irun) = abs(H5(:,1)'*H5(:,2))/(norm(H5(:,1))*norm(H5(:,2)));
            H_corr13or24(r_val,irun) = abs(H5(:,1)'*H5(:,3))/(norm(H5(:,1))*norm(H5(:,3)));
            H_corr23or14(r_val,irun) = abs(H5(:,3)'*H5(:,2))/(norm(H5(:,3))*norm(H5(:,2)));
            % if H_corr12or34(r_val,irun) > 0.93
            %     H_corr12or34(r_val,irun) = 1;
            % end
            % if H_corr13or24(r_val,irun) > 0.93
            %     H_corr13or24(r_val,irun) = 1;
            % end

%%
% %Calculating NE
        R_mat = [1 H_corr12or34(r_val,irun) H_corr13or24(r_val,irun) H_corr12or34(r_val,irun);H_corr12or34(r_val,irun) 1 H_corr12or34(r_val,irun) H_corr13or24(r_val,irun);H_corr13or24(r_val,irun) H_corr12or34(r_val,irun) 1 H_corr12or34(r_val,irun);H_corr12or34(r_val,irun) H_corr13or24(r_val,irun) H_corr12or34(r_val,irun) 1];
        Rinv=pinv(R_mat+0.00001*eye(4)); %zero-forcing, small diagonal loading to ensure invertible
         mean_ZF(r_val,irun) = 10*log10(Rinv(1,1));


    RI=R_mat(2:4,2:4); %interference correlation matrix
    rhoI=R_mat(1,2:4); %correlation of signal vector with interference vectors
    %proportional of signal energy orthogonal to interference subspace
    %this is the inverse of the noise enhancement, and is a more stable
    %computation
    signal_remaining=1-rhoI*pinv(RI)*transpose(rhoI); 
    %saturate noise enhancement at 60 dB
    if signal_remaining > 10^(-6)
        mean_ZF(r_val,irun) = 10*log10(1/signal_remaining);
    else
        mean_ZF(r_val,irun) =60;
    end


%%

        %% Checking the ZF Noise Enhancement 

                %C_zf = zeros(Nrx,Ntx);
                C_zf = H3'*pinv(H3*H3');
                %C_zf = H5'*pinv(H5*H5');
                %Need to chnage to append max val instead of a random val
                Czf_R(r_val,irun) =10*log10(Ntx*norm(C_zf(93,:))^2);%+(3*log2(Nrx) -3*log2(4));


                % % %Mean Noise Enhancement
                % [U,S,V] = svd(H5);
                % %[U,S,V] = svd([H1 H2 H3 H4], 'econ');
                % singVals = diag(S);
                % sumZF = 0;
                % for si = 1:length(singVals)
                %     sumZF = sumZF + (1/singVals(si)^2);
                % end
                % 
                % mean_ZF(r_val,irun) = 10*log10(sumZF/(Ntx));% + (3*log2(Nrx) -3*log2(4));

                % %Another way to compute the avg ZF NE
                % sumZF = 0;
                % for itxx = 1:Ntx
                %     h_val = H_est_model{itxx};
                %     %h_val = H_est_model{itxx}/ sqrt(Ntx);
                %     C_zf1 = h_val'*pinv(h_val*h_val');
                %    % sumZF =sumZF + (Ntx*norm(C_zf1(H_ind(itxx),:))^2);
                %    sumZF =sumZF + ((Ntx/2)*norm(C_zf1(H_ind(itxx),:))^2);
                % 
                % end
                % %mean_ZF(r_val,irun) = 10*log10(sumZF/(Nrx/Ntx));
                % mean_ZF(r_val,irun) = 10*log10(sumZF);

                %Alt way to fully capture isi + csi with normalization
% sumZF = 0;
% 
% for itx = 1:Ntx
%     H_blk = H_est_model{itx} / sqrt(Ntx);     % Normalize total power
%     W_zf = pinv(H_blk' * H_blk) * H_blk';     % Nt x (Nrx*W)
% 
%     for row = 1:size(W_zf,1)
%         sumZF = sumZF + norm(W_zf(row,:))^2;  % NE for this stream
%     end
% end
% 
% % Divide by total number of streams (Ntx streams in each of Ntx views)
% mean_ZF(r_val, irun) = 10 * log10(sumZF / (Ntx * Ntx));



           
            for isnr = 1:Nsnr
                SISO_SNR_dB = SISO_SNR_dB_vect(isnr);
                SISO_SNR_lin = 10^(SISO_SNR_dB/10);
                sigmasq = 1/SISO_SNR_lin;
                
                % let's tell the folks how we're doing
                % numdone = numdone + 1;
                % clc;
                % disp(['run ',num2str(irun),' of ',num2str(Nrun)])
                % disp([num2str(round(numdone/numtot*100)),'% done'])
                
                Yrx = Yrx_noiseless +  sqrt(sigmasq/2)*(randn(size(Yrx)) + 1i*randn(size(Yrx)));
                Yrx_tr = Yrx_tr_noiseless +   1*sqrt(sigmasq/2)*(randn(size(Yrx)) + 1i*randn(size(Yrx)));
 
                C_LMMSE = cell(Ntx,1);
                % Xgrid = (-W - L_pulse - max_del_whole)*1 : 1*(max_del_whole + W + L_pulse);
                % Xgridlen = length(Xgrid);
                grid_center_ind = find(Xgrid == 0);
                SINR_vect = zeros(Ntx,1);
                for itx = 1:Ntx
                     H = H_est_model{itx};
                    % H = H_model{itx};
                    desired_ind = (itx - 1)*Xgridlen + grid_center_ind;
                    desired_delta = zeros(size(H,2),1);
                    desired_delta(desired_ind) = 1;
                    C = pinv(H'*H + sigmasq*eye(size(H,2)))*H';
                    C = C(desired_ind,:);
                    correction_factor = abs(C*H*desired_delta);
                    C_LMMSE{itx} = C/(1-correct_flag + correct_flag*correction_factor);
                    signal_power = correction_factor^2;
                    noise_power = norm(C)^2*sigmasq;
                    int_power = 0;
                    for i_int = 1:size(H,2)
                        if i_int ~= desired_ind
                            int_delta = zeros(size(H,2),1);
                            int_delta(i_int) = 1;
                            int_power = int_power + abs(C*H*int_delta)^2;
                        end
                    end
                    SINR_vect(itx) = signal_power/(noise_power + int_power);
                end
                SINR_analytical(iredund,iwside,isnr) = SINR_analytical(iredund,iwside,isnr) + 1/Nrun*mean(SINR_vect);
                BER_analytical(iredund,iwside,isnr) = BER_analytical(iredund,iwside,isnr) + 1/Nrun*mean(qfunc(sqrt(SINR_vect)));
                ratio_lin_anal(iredund,iwside,isnr) = ratio_lin_anal(iredund,iwside,isnr) + 1/Nrun*mean(SINR_vect)/(Nrx/sigmasq);
                
                Y_win = cell(Ntx,1);
                Y_win_tr = cell(Ntx,1);
                for itx = 1:Ntx
                    ywin = zeros(Nrx*W,ylen);
                    ywin_tr = zeros(Nrx*W,ylen);
                    for iw = 1:W
                        w = -W_side + iw - 1;
                        for irx = 1:Nrx
                            ywin((iw-1)*Nrx + irx,:) = circshift(Yrx(irx,:),-(delays_whole(irx,itx)+w));
                            ywin_tr((iw-1)*Nrx + irx,:) = circshift(Yrx_tr(irx,:),-(delays_whole(irx,itx)+w));
                        end
                    end
                    Y_win{itx} = ywin;
                    Y_win_tr{itx} = ywin_tr;
                end
                
                % now let's receive with the computed receivers
                Xest = zeros(Ntx,size(Yrx,2));
                for itx = 1:Ntx
                    C = C_LMMSE{itx};
                    ywin = Y_win{itx};
                    Xest(itx,:) = C*ywin;
                end
                
                Xtx_clipped = Xtx(:,26:end-25);
                Xest_clipped = Xest(:,26:end-25);
                len1 = size(Xest_clipped,2); % should be 10000
                EVM_mat = abs(Xest_clipped - Xtx_clipped);
                Ph_err_mat = angle(Xest_clipped./Xtx_clipped);
                BER_LMMSE_MC = (sum(abs(Ph_err_mat)'>pi/4) + sum(abs(Ph_err_mat)'>3*pi/4)*0)/len1/2;
                
                SNR_SISO_dB = SISO_SNR_dB;
                SNR_beamformed_dB = SISO_SNR_dB + 10*log10(Nrx);
                SINR_LMMSE_dB = 10*log10(mean(mean(abs(Xtx_clipped.').^2)./mean((EVM_mat.').^2)));
                
                SINR_MonteCarlo(iredund,iwside,isnr) = SINR_MonteCarlo(iredund,iwside,isnr) + 1/Nrun*SINR_LMMSE_dB;
                BER_MonteCarlo(iredund,iwside,isnr) = BER_MonteCarlo(iredund,iwside,isnr) + 1/Nrun*mean(BER_LMMSE_MC);
                ratio_dB_EVM(iredund,iwside,isnr) = ratio_dB_EVM(iredund,iwside,isnr) + 1/Nrun*(SINR_LMMSE_dB - SNR_beamformed_dB);
                
            end
end
H_corr12or34_avg(r_val) = sum(H_corr12or34(r_val,:))/Nrun;
H_corr13or24_avg(r_val) = sum(H_corr13or24(r_val,:))/Nrun;
mean_ZF_avg(r_val) = sum(mean_ZF(r_val,:))/Nrun;
Czf_R_avg(r_val) = sum(Czf_R(r_val,:))/Nrun;

end


%  figure
% % plot(R_vect/100,Czf_R_avg);
% % hold on;
%  plot(R_vect/100,mean_ZF_avg- 10*log10(Nrx/Ntx),'--');
% 
% figure
% plot(R_vect/100,H_corr12or34_avg);
% 
% figure
% plot(R_vect/100,H_corr13or24_avg);
% 
% % figure
% % plot(R_vect/100, (Czf_R_avg - 10*log10(Nrx/Ntx) + 20*log10(R_vect/100)'))
% 
% figure
% plot(R_vect/100, ((mean_ZF_avg)- 10*log10(Nrx/Ntx) + 20*log10(R_vect/100)'))


% plot stuff 

for plotstuff=[1]
    close all

    marker_list = ['*'; 's'; 'o'; 'd'; 'x'; '^'; '>'; 'v' ; '<'; 'p'];
    color_list = [.85,.33,.10; .47,.67,.19;   .00,.45,.74; 1.0,.84,.00; 1.0,.60,.78;  .60,.20,.00;   .00,.75,.75; .87,.49,.00; .08,.17,.55;   
             ];
    fig1_ratio_LMMSE = figure();
    hold on; box on; grid on;
    fig2_BER_SER_LMMSE = figure();
    box on; grid on;
    leg1 = cell(0);
    leg2 = cell(0);
    xx0 = SISO_SNR_dB_vect;
    if xx0(end)==inf; xx0(end) = 100; end
    for iredund = 1:Nredund
        nrx = redund_vect(iredund)*4;
        linwid = iredund/4 + 0.5;
        color_choice = color_list(1+mod(iredund-1,size(color_list,1)),:);
        for iwside = 1:Nwside
            marker_choice = marker_list(1+mod(iwside-1,length(marker_list)));
            w = 1 + 2*W_side_vect(iwside);
            figure(fig1_ratio_LMMSE)
            temp = ratio_dB_EVM(iredund,iwside,:);
            plot(xx0,temp(:),['-',marker_choice],'LineWidth',linwid,'Color',color_choice);
            leg1{end+1} = ['N_{rx} = ',num2str(nrx),', W = ',num2str(w)];
            temp = 10*log10(ratio_lin_anal(iredund,iwside,:));
            plot(xx0,temp(:),'--','LineWidth',linwid,'Color',color_choice,'HandleVisibility','off');

            figure(fig2_BER_SER_LMMSE);
            little_bit = 0*1e-20;
            temp2 = little_bit + BER_MonteCarlo(iredund,iwside,:);
            temp1 = little_bit + BER_analytical(iredund,iwside,:);
            semilogy(xx0,temp1(:),['-',marker_choice],'LineWidth',linwid,'Color',color_choice); hold on;
            semilogy(xx0,temp2(:),'--','LineWidth',linwid,'Color',color_choice,'HandleVisibility','off');
            leg2{end+1} = ['N_{rx} = ',num2str(nrx),', W = ',num2str(w)];
            % leg2{end+1} = ['num''l, N_{rx} = ',num2str(nrx),', W = ',num2str(w)];
        end
    end

    figure(fig1_ratio_LMMSE)
    legend(leg1)
    xlabel('SISO SNR (dB)')
    ylabel('ratio SINR/SNR_{beamformed} (dB)')

    figure(fig2_BER_SER_LMMSE)
    legend(leg2)
    xlabel('SISO SNR (dB)')
    ylabel('BER')
    ylim([1e-7,1])
    % (Nredund,Nwside,Nsnr)


end

% 
% Visualization of Rx antenna locations
figure;
scatter(rxloc(:,1), rxloc(:,2), 100, 'bo', 'filled'); % Rx antennas in blue
hold on;
scatter(txloc(:,1), txloc(:,2), 100, 'rs', 'filled'); % Tx antennas in red

Annotate points
for i = 1:size(rxloc,1)
    text(rxloc(i,1), rxloc(i,2), sprintf('Rx%d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
end
for i = 1:size(txloc,1)
    text(txloc(i,1), txloc(i,2), sprintf('Tx%d', i), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 12);
end

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Antenna Locations (Rx)');
grid on;
axis equal;
legend({'Rx Antennas'});
hold off;


% % Define the positions of RX1 to RX7 (unchanged)
% rx_positions = [
%     0, 0;   % RX1
%     1, 0;   % RX2
%     1, 1;   % RX3
%     0, 1;   % RX4
%     0.5, 0.5; % RX5 (Midpoint of RX1 and RX4)
%     mean([0, 0; 1, 0; 1, 1; 0, 1], 1);  % RX6 (Center of RX array)
%     (rx_positions(3,:) + rx_positions(4,:)) / 2;  % RX7 (Midpoint of RX3 and RX4)
% ];
% 
% % Calculate new RX8 to RX10 positions
% rx_8 = rx_positions(1,:) + (rx_positions(2,:) - rx_positions(1,:)) * 0.25;  % 1/4th distance between RX1 and RX2
% rx_9 = rx_positions(3,:) + (rx_positions(2,:) - rx_positions(3,:)) * 0.25;  % 1/4th distance between RX3 and RX2
% rx_10 = (rx_positions(6,:) + rx_positions(2,:)) / 2;  % Midpoint of RX6 and RX2 (diagonal)
% 
% % Combine all RX positions into one matrix
% all_rx_positions = [
%     rx_positions;
%     rx_8;
%     rx_9;
%     rx_10
% ];
% 
% % Plot the antenna positions
% figure;
% scatter(all_rx_positions(:, 1), all_rx_positions(:, 2), 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'b');
% hold on;
% 
% % Annotate the points with their RX labels
% for i = 1:size(all_rx_positions, 1)
%     text(all_rx_positions(i, 1), all_rx_positions(i, 2), [' RX', num2str(i)], 'FontSize', 12, 'HorizontalAlignment', 'right');
% end
% 
% % Set axis labels and title
% xlabel('X Position');
% ylabel('Y Position');
% title('Receiver Antenna Positions');
% grid on;
% axis equal;
% hold off;



%% functions

% Introducing the misalignment
function loc = tilt_array (loc, vtilt, htilt)
D = mean(loc(:,3));
loc(:,3) = loc(:,3) - D;
vrot=[cos(vtilt),-sin(vtilt);sin(vtilt),cos(vtilt)];
hrot=[cos(htilt),-sin(htilt);sin(htilt),cos(htilt)];
loc(:,2:3) = loc(:,2:3)*vrot';
loc(:,[1,3]) = loc(:,[1,3])*hrot';
loc(:,3) = loc(:,3) + D;
end

% Raised cosine function
function z = raisedcosine_generator(t,Ts,alpha)
z = (cos(pi*alpha*((t)/Ts)) ./ (1-(2*alpha*((t)/Ts)).^2)) .* sinc(t/Ts) ; %( (sin(pi*((t)/Ts))) ./ (pi*((t)/Ts)) );
z(abs(t/Ts)==2) = 0;
end

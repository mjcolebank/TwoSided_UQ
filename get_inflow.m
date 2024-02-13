function get_inflow(amp_perturb,amp_ind,save_ind)
load Qin_orig.dat

temp = fft(Qin_orig);
%%
imag_num = sqrt(-1);

% amp_ind = [1 2 3]; % Which components to change
n_amps_sign = sign(real(temp));
n_pha_sign  = sign(imag(temp));
new_Q = temp;

% Now just change the amplitude
for j=1:length(amp_ind)%j=amp_ind
    if amp_ind(j)==1
        new_Q(amp_ind(j)) = amp_perturb(j);%unifrnd(0.8*n_amps(j),1.20*n_amps(j),1,1); % Unif
    else
        new_Q(amp_ind(j)) = n_amps_sign(amp_ind(j)).*amp_perturb(j) + sqrt(-1)*imag(new_Q(j)); % Uniform
        new_Q(end-(amp_ind(j)-2)) = real(new_Q(amp_ind(j))) + imag_num*imag(new_Q(end-(amp_ind(j)-2)));
    end
end

%% or change the phase?
%     for j=n_change
%         if j==1
%             new_LA(j) = new_LA(j).*normrnd(1,0.5,1,1);
%         else
%         new_LA(j) = real(new_LA(j)) + sqrt(-1)*imag(new_LA(j)).*normrnd(1,0.3,1,1);
%         new_LA(end-(j-2)) = real(new_LA(j)) - sqrt(-1)*imag(new_LA(j));
%         end
%     end
q_test = ifft(new_Q);



% smooth the data to ensure that we (1) circumvent numerical instabilities,
% and (2) have flows that start and end in the same place
Q_filter = [q_test; q_test];
Q_filter2 = smoothdata(Q_filter,1,'gaussian','SmoothingFactor',0.01);

Q_test = Q_filter2(1:8193,:);
figure(99); hold on; plot(Q_test)

fname = strcat('Qin_',num2str(save_ind),'.dat');
% dlmwrite(fname,Q_test);
writematrix(Q_test,fname);
end
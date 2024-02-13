function get_outletpressure(amp_perturb,pha_perturb,amp_ind,pha_ind,save_ind)
load LA_orig.dat
amp_perturb
temp = fft(LA_orig);
%%
imag_num = sqrt(-1);

% amp_ind = [1 2 3]; % Which components to change
n_amps = abs(real(temp));
n_amps_sign = sign(real(temp));
n_pha_sign  = sign(imag(temp));
new_LA = temp;

% Now just change the amplitude
for j=1:length(amp_ind)%j=amp_ind
    if amp_ind(j)==1
        new_LA(amp_ind(j)) = amp_perturb(j);%unifrnd(0.8*n_amps(j),1.20*n_amps(j),1,1); % Unif
    else
        new_LA(amp_ind(j)) = n_amps_sign(amp_ind(j)).*amp_perturb(j) + sqrt(-1)*imag(new_LA(amp_ind(j))); % Uniform
        new_LA(end-(amp_ind(j)-2)) = real(new_LA(amp_ind(j))) + imag_num*imag(new_LA(end-(amp_ind(j)-2)));
    end
end

%% or change the phase?
    for j=1:length(pha_ind)%j=amp_ind
        if j==1
            % Do nothing; zero phase
        else
        new_LA(pha_ind(j)) = real(new_LA(pha_ind(j))) + n_pha_sign(pha_ind(j)).*sqrt(-1)*pha_perturb(j);
        new_LA(end-(pha_ind(j)-2)) = real(new_LA(pha_ind(j))) - sqrt(-1)*imag(new_LA(pha_ind(j)));
        end
    end
LA_test = ifft(new_LA);



% smooth the data to ensure that we (1) circumvent numerical instabilities,
% and (2) have flows that start and end in the same place
LA_filter = [LA_test; LA_test];
LA_filter2 = smoothdata(LA_filter,1,'gaussian','SmoothingFactor',0.01);

LA_test = LA_filter2(1:8193,:);
figure(99); hold on; plot(LA_test); plot(LA_orig,'k')
% Test for negative pressure
min_LA = min(LA_test);
LA_test = LA_test-min(min_LA,0);
fname = strcat('LAout_',num2str(save_ind),'.dat');
% dlmwrite(fname,Q_test);
% writematrix(LA_test,fname);
end
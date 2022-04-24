%%
%%%%%%%%%%%%%%%%%%%% Entering the Input Signal Matrix and T (Highest signal period) %%%%%%%%%%%%%%%%

M = input("Enter the M signals in the matrix form: \n");
d_M = size(M);

T = input("Enter the Highest Signal Duration (T): \n"); %Highest Signal Duration (T)
n = d_M(1); %number of bits (time slots)
Tb = T/n;   %Time of one Bit

%%%%%%%%%%%%%%%%%%% Plot the input signals %%%%%%%%%%%%%%%%%%%%

M %To show the Input siganal matrix

for j = 1 : d_M(2)
subplot(1,d_M(2),j)

y1 = stairs(M(:,j));  %plotting Signals (S1 --> SN)
l= 0:1:length(M(:,j));
y1.XData = Tb.*[l l(length(M(:,j)))+1] ;
y1.YData = [(M(:,j))' M(d_M(1),j) 0];

ylim([min(min(M))-1 max(max(M))+1]); %adjusting X-axis and Y-axis
xlim([0,Tb.*d_M(1)+Tb])
xticks(0:Tb:Tb*d_M(1)+Tb)

title_txt1 = sprintf('Input signal s_%d(t)' , j);
title_txt2 = sprintf('s_%d(t)' , j);
title(title_txt1);
xlabel('time(sec)');
ylabel(title_txt2);
grid on;
end

%%
%%%%%%%%%%%%%%%%%%%%% Getting the basis functions %%%%%%%%%%%%%%%%%%%%

%getting phi(Basis functions) and Coeff form our implemented function
[phi,Coeff] = Gram_schmidt(M,T); 
d_N = size(phi); %size of phi matrix (basis functions)

%%%%%%%%%%%%%%%%%%% Plotting the basis functions %%%%%%%%%%%%%%%%%%%%

figure;
for i = 1 : d_N(2)
subplot(1,d_N(2),i)

y = stairs(phi(:,i));  %plotting Signals (phi1 --> phiM)
l= 0:1:length(phi(:,i));
y.XData = Tb.*[l l(length(phi(:,i)))+1] ;
y.YData = [(phi(:,i))' phi(d_N(1),i) 0];

ylim([min(min(phi))-1 max(max(phi))+1]); %adjusting X-axis and Y-axis
xlim([0,Tb.*d_N(1)+Tb])
xticks(0:Tb:Tb*d_N(1)+Tb)

title("Basis Function \phi_"+int2str(i));
xlabel('time(sec)');
ylabel("\phi_"+int2str(i));
grid on;
end

%%
%%%%%%%%%%%%%%%%%%%%% Check the basis functions %%%%%%%%%%%%%%%%%%%%

d_Matrix = size(Coeff); % Size of coefficients Matrix 
M_Reconstructed = zeros(d_N(1),d_Matrix(1)); %initializing size of reconstructed Signals Matrix

%reconstructing the signals Matrix form coefficients Matrix and Basis functions 
for j = 1 : d_Matrix(1)  
    
    signal_vec = Coeff(j,:);
    signal = (signal_vec(1)*phi(:,1));
    
    for i = 2:d_N(2)
       signal = signal+(signal_vec(i)*phi(:,i));
    end
    M_Reconstructed(:,j)=signal;
end

M_Reconstructed %to show Reconstructed Signals Matrix in command window


figure;    %Plotting the Reconstructed Signals which is the same as entered Matrix
for j = 1 : d_Matrix(1)
    
    signal_vec = Coeff(j,:);   %reconstructing the signals Matrix again to plot it
    signal = (signal_vec(1)*phi(:,1));
    
    for i = 2:d_N(2)
    signal = signal+(signal_vec(i)*phi(:,i));
    end
    
    subplot(1,d_Matrix(1),j) %plotting Reconstructed Signals (S1 --> SN)
    y1 = stairs(signal);
    l= 0:1:length(signal);
    y1.XData = Tb.*[l l(length(signal))+1] ;
    y1.YData = [signal' signal(length(signal)), 0];
    
    ylim([min(min(M))-1 max(max(M))+1]);
    xlim([0, Tb.*d_N(1)+Tb])
    xticks(0:Tb: Tb*d_N(1)+Tb)
    
    title_txt1 = sprintf('Reconstructed signal s_%d(t)' , j);
    title_txt2 = sprintf('s_%d(t)' , j);
    title(title_txt1);
    xlabel('time(sec)');
    ylabel(title_txt2);
    grid on;
end   
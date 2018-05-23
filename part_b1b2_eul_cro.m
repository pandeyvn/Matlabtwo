%% This is Euler Cromer algorithm. Molecular Dynamics 9(a) page 96
%% Both variants are implemented, variant = 1 or 2 selects the desired variant.
clear all;
close all;
format long
tau=[0.2 0.02 0.002];
N=100; % The number of seconds
xlim_low=00; xlim_high=xlim_low+100;

% The variable variant has to be either 1 or 2 to select the preferred variant.
% It is 1 when we use (a) implementation Molecular Dynamics 9(a)
% It is 2 when we use (b) implementation Molecular Dynamics 9(a)

variant = 1;  % or 2

str_var='ab';

x_init=0; v_init=1;
x(1)=x_init; % at time 0
v(1)=v_init; % at time 0
t(1)=0;
m=1; e =0.5;
k=1;p(1)=1; pe(1)=(x(1)*x(1))/2;ke(1)=0.5*1*v(1)*v(1);te(1)=pe(1)+ke(1);
for i = 1:length(tau)
    for j= 1:N/tau(i) %1:1000/tau(i) % j is the number of steps and nothing else.
        t(j+1)=j*tau(i); % present discrete time This does not hold for j=0
        % finding x and the new time
        % x(j) is the distance at the jth instance which is at time j*tau
        vdash(j)=k*x(j);
        
        switch variant
            case 1  %% This is variant a
                p(j+1)=p(j)-tau(i)*vdash(j);
                x(j+1)=x(j)+ tau(i)*p(j+1);
            case 2  %% This is variant b
                x(j+1)=x(j)+ tau(i)*p(j);
                p(j+1)=p(j)-tau(i)*k*x(j+1); % because vdash=kx
            otherwise % Means variant variable not set properly.
                disp('The variant value is neither 1 or 2 ');
                quit
        end
        
        v(j+1)=p(j+1);
        pe(j+1)=(x(j+1)*x(j+1))/2;
        ke(j+1)=0.5*1*v(j+1)*v(j+1);
        te(j+1)=pe(j+1)+ke(j+1);
    end
    
    exact_disp(1:length(t))=sin(t);
    exact_energy(1:length(t))=0.5;
    
    figure(1);
    subplot(length(tau),1,i);
    plot(t,x,'g.',t,exact_disp,'r');
    legend('x(t)','sin(t)')
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Displacement');
    tit_str=['Euler Cromer variant (' num2str(str_var(variant))   ')  and \tau = ' num2str(tau(i)) ];
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    max(x-exact_disp)
    min(x-exact_disp)
    
    
    figure(2);
    subplot(length(tau),1,i)
    plot(t,x-exact_disp,'k');
    legend('x(t)-sin(t)')
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Error in displacement');
    tit_str=['Euler Cromer variant (' num2str(str_var(variant))   ')  and \tau = ' num2str(tau(i)) ];
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    
    figure(3)
    subplot(length(tau),1,i);
    plot(t,te,'b',t,exact_energy,'k')
    legend('Actual Energy','Exact Energy')
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Energy' );
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    
    figure(4)
    subplot(length(tau),1,i);
    plot(t,te,'g.',t,exact_energy,'r',t,pe,'b')
    legend('Total Energy','True Total Energy','Potential Energy')
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Energy' );
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    
    figure(5)
    subplot(length(tau),1,i);
    plot(t,te-exact_energy,'k')
    legend('Total Energy - True Total Energy');
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Error in Total Energy' );
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    
    figure(6)
    subplot(length(tau),1,i);
    plot(t,v,'k')
    legend('velocity')
    lim_ax=axis;
    axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
    xlabel ('Time in seconds (t)');
    ylabel('Velocity');
    title(tit_str)
    grid on;
    remove_extra_space(gca);
    
    if(tau(i) == 0.02)
        figure(7)
        plot(t,x,'r+',t,sin(t),'b',t,v,'gx',t,te,'m');
        legend('x(t)','sin(t)', 'v(t)','E(t)');
        lim_ax=axis;
        axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
        xlabel('Time in seconds (t)');
        ylabel('x(t),v(t)');
        title(tit_str)
        grid on;
	remove_extra_space(gca);
        %%%%%%%%%%%%%
        figure(8)
        plot(t,exact_energy,'k',t,te,'m');
        legend('Exact','E(t)');
        lim_ax=axis;
        axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
        xlabel('Time in seconds (t)');
        ylabel('Energy');
        title(tit_str)
        grid on;
	remove_extra_space(gca);
        figure(9)
        plot(x,v,'r.');
        lim_ax=axis;
 %       axis([xlim_low xlim_high lim_ax(3) lim_ax(4)]);
        xlabel('Displacement X(t)');
	supl_str=['    Time 00s  to  ',num2str(xlim_high),'s'];
        ylabel('Momentum');
        title([tit_str supl_str])
        grid on;
	remove_extra_space(gca);
	
    end
    
    %pause
    
end

fig_max=9;
file_str=['Euler_Cromer_Variant_',str_var(variant),];
supl_str=['_limit_',num2str(xlim_low),'_',num2str(xlim_high),'s'];
for i=1:fig_max
    file_name=[file_str,supl_str,'_figure',num2str(i),'.jpg']
    figure(i);
%    pause
     saveas(gcf,file_name,'jpg')
end


function phase_evoltn()
% PHASE_EVOLUTION
%	A Matlab version of the GRASP example: phase_evolution.
%	Example run: phase_evolution
%
% Steve Drasco
% Summer 1998

m1 = 1.4;
m2 = 1.4;
spin1 = 0;
spin2 = 0;
n_phaseterms = 5;
Initial_Freq = 60;
Max_Freq_Rqst = 2000;
Sample_Time = 1/9868.4208984375;
err_cd_sprs = 0;
phaseterms = [1 0 1 1 1];

[Max_Freq_Actual,phase,frequency,steps_filld,clscnc_time]=mxPhase_frequency(m1,...
	m2,spin1,spin2,n_phaseterms,phaseterms,Initial_Freq,Max_Freq_Rqst,Sample_Time,[], err_cd_sprs);

time = (1:steps_filld) * Sample_Time;

subplot(2,1,1)
	plot(time, phase,'b');
	xlabel('time (s)');
	ylabel('phase');
subplot(2,1,2)
	plot(time, frequency,'r');
        xlabel('time (s)');
        ylabel('frequency (Hz)');

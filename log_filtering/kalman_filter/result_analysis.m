clc;
clear all;
close all;


zero_kg = load('0kg','eval_vf_time','fx','fy','fz');
KF2_0kg = load('KF2_0kg','eval_fvf_time','ffx','ffy','ffz');
three_kg = load('3kg','eval_vf_time','fx','fy','fz');
KF2_3kg = load('KF2_3kg','eval_fvf_time','ffx','ffy','ffz');

ts_lift_0kg = find(zero_kg.eval_vf_time>7.5);
te_lift_0kg = find(zero_kg.eval_vf_time>24);
ts_drop_0kg = find(zero_kg.eval_vf_time>38);
te_drop_0kg = find(zero_kg.eval_vf_time>58);

ts_lift_3kg = find(three_kg.eval_vf_time>9);
te_lift_3kg = find(three_kg.eval_vf_time>24);
ts_drop_3kg = find(three_kg.eval_vf_time>51);
te_drop_3kg = find(three_kg.eval_vf_time>69);


figure ('Name','Virtual Force')

hold on
plot(zero_kg.eval_vf_time(ts_lift_0kg(1):te_lift_0kg(1)),zero_kg.fx(ts_lift_0kg(1):te_lift_0kg(1)),'b')
plot(three_kg.eval_vf_time(ts_lift_3kg(1):te_lift_3kg(1)),three_kg.fx(ts_lift_3kg(1):te_lift_3kg(1)),'r')
plot(zero_kg.eval_vf_time(ts_drop_0kg(1):te_drop_0kg(1)),zero_kg.fx(ts_drop_0kg(1):te_drop_0kg(1)),'b')
plot(three_kg.eval_vf_time(ts_drop_3kg(1):te_drop_3kg(1)),three_kg.fx(ts_drop_3kg(1):te_drop_3kg(1)),'r')
xlabel('time')
ylabel('fx')
legend('0kg-fx','3kg-fx')

suptitle('Virtual Force')

figure ('Name','Virtual Force filtered')

hold on
plot(KF2_0kg.eval_fvf_time(ts_lift_0kg(1):te_lift_0kg(1)),KF2_0kg.ffx(ts_lift_0kg(1):te_lift_0kg(1)),'b')
plot(KF2_3kg.eval_fvf_time(ts_lift_3kg(1):te_lift_3kg(1)),KF2_3kg.ffx(ts_lift_3kg(1):te_lift_3kg(1)),'r')
plot(KF2_0kg.eval_fvf_time(ts_drop_0kg(1):te_drop_0kg(1)),KF2_0kg.ffx(ts_drop_0kg(1):te_drop_0kg(1)),'b')
plot(KF2_3kg.eval_fvf_time(ts_drop_3kg(1):te_drop_3kg(1)),KF2_3kg.ffx(ts_drop_3kg(1):te_drop_3kg(1)),'r')
xlabel('time')
ylabel('fx')
legend('0kg-fx','3kg-fx')

suptitle('Virtual Force filtered')

% % org  = load('org','eval_mjs_time','mq1','mDq1','mDDq1');
% % PID0 = load('PID0','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % PID1 = load('PID1','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % PID2 = load('PID2','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % PID3 = load('PID3','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % PID4 = load('PID4','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % PID5 = load('PID5','eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % P0   = load('P0',  'eval_mjs_time','mq1','mDq1','eval_fjs_time','fq1','fDq1','fDDq1');
% % 
% % %org
% % eval_mjs_time = org.eval_mjs_time;
% % mq1 = org.mq1;
% % mDq1 = org.mDq1;
% % mDDq1 = org.mDDq1;
% % 
% % %PID0
% % f0_time = PID0.eval_fjs_time;
% % f0q1    = PID0.fq1;
% % f0Dq1   = PID0.fDq1;
% % f0DDq1  = PID0.fDDq1;
% % 
% % %PID1
% % f1_time = PID1.eval_fjs_time;
% % f1q1    = PID1.fq1;
% % f1Dq1   = PID1.fDq1;
% % f1DDq1  = PID1.fDDq1;
% % 
% % %PID2
% % f2_time = PID2.eval_fjs_time;
% % f2q1    = PID2.fq1;
% % f2Dq1   = PID2.fDq1;
% % f2DDq1  = PID2.fDDq1;
% % 
% % %PID3
% % f3_time = PID3.eval_fjs_time;
% % f3q1    = PID3.fq1;
% % f3Dq1   = PID3.fDq1;
% % f3DDq1  = PID3.fDDq1;
% % 
% % %PID4
% % f4_time = PID4.eval_fjs_time;
% % f4q1    = PID4.fq1;
% % f4Dq1   = PID4.fDq1;
% % f4DDq1  = PID4.fDDq1;
% % 
% % %PID5
% % f5_time = PID5.eval_fjs_time;
% % f5q1    = PID5.fq1;
% % f5Dq1   = PID5.fDq1;
% % f5DDq1  = PID5.fDDq1;
% % 
% % %P0
% % fp0_time = P0.eval_fjs_time;
% % fp0q1    = P0.fq1;
% % fp0Dq1   = P0.fDq1;
% % fp0DDq1  = P0.fDDq1;
% % 
% % %%
% % 
% % figure ('Name','Joint 1 Position  with Filter')
% % 
% % plot(eval_mjs_time,mq1)
% % hold on
% % plot(fp0_time,fp0q1)
% % hold on
% % plot(f0_time,f0q1)
% % % hold on
% % % plot(f1_time,f1q1)
% % % hold on
% % % plot(f2_time,f2q1)
% % % hold on
% % % plot(f3_time,f3q1)
% % % hold on
% % % plot(f4_time,f4q1)
% % hold on
% % plot(f5_time,f5q1)
% % 
% % xlabel('time')
% % ylabel('q1')
% % %legend('org','P0','PID0','PID1','PID2','PID3','PID4','PID5')
% % legend('org','P0','PID0','PID5')
% % 
% % suptitle('Joint 1 Position  with Filter')
% % 
% % %%
% % 
% % figure ('Name','Joint 1 Velocity  with Filter')
% % 
% % plot(eval_mjs_time,mDq1)
% % hold on
% % plot(fp0_time,fp0Dq1)
% % hold on
% % plot(f0_time,f0Dq1)
% % % hold on
% % % plot(f1_time,f1Dq1)
% % % hold on
% % % plot(f2_time,f2Dq1)
% % % hold on
% % % plot(f3_time,f3Dq1)
% % % hold on
% % % plot(f4_time,f4Dq1)
% % hold on
% % plot(f5_time,f5Dq1)
% % 
% % 
% % xlabel('time')
% % ylabel('Dq1')
% % %legend('org','P0','PID0','PID1','PID2','PID3','PID4','PID5')
% % legend('org','P0','PID0','PID5')
% % 
% % suptitle('Joint 1 Velocity  with Filter')
% % 
% % %%
% % 
% % figure ('Name','Joint 1 Accelearation')
% % 
% % plot(eval_mjs_time,mDDq1)
% % hold on
% % plot(fp0_time,fp0DDq1)
% % hold on
% % plot(f0_time,f0DDq1)
% % % hold on
% % % plot(f1_time,f1DDq1)
% % % hold on
% % % plot(f2_time,f2DDq1)
% % % hold on
% % % plot(f3_time,f3DDq1)
% % % hold on
% % % plot(f4_time,f4DDq1)
% % hold on
% % plot(f5_time,f5DDq1)
% % 
% % 
% % xlabel('time')
% % ylabel('fDDq1')
% % %legend('org','P0','PID0','PID1','PID2','PID3','PID4','PID5')
% % legend('org','P0','PID0','PID5')
% % 
% % suptitle('Joint 1 Accelearation')
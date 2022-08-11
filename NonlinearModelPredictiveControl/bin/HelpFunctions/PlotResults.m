figure(1);clf;
subplot(2,2,1)
stairs(ops.t(1,1:kk),x(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(1)*ones(1,kk),'r--');
ylabel('$x_1(k)$','fontsize',14,'interpreter','latex')
subplot(2,2,2)
stairs(ops.t(1,1:kk),x(2,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(2)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(2)*ones(1,kk),'r--');
ylabel('$x_2(k)$','fontsize',14,'interpreter','latex')
subplot(2,2,3)
stairs(ops.t(1,1:kk),u(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.u_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.u_min(1)*ones(1,kk),'r--');
ylabel('$u(k)$','fontsize',14,'interpreter','latex')
xlabel('$t$ (s)','fontsize',14,'interpreter','latex')
subplot(2,2,4)
stairs(ops.t(1,1:kk),y(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),r(1,1:kk),'k--','linewidth',1.5);
ylabel('$y(k)$','fontsize',14,'interpreter','latex')
xlabel('$t$ (s)','fontsize',14,'interpreter','latex')

figure(2);clf;
subplot(1,2,1)
stairs(ops.t(1,1:kk),u(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.u_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.u_min(1)*ones(1,kk),'r--');
ylim([0 ops.u_max(1)*1.1])
ylabel('$u(k)$','fontsize',14,'interpreter','latex')
xlabel('$t$ (s)','fontsize',14,'interpreter','latex')
subplot(1,2,2)
stairs(ops.t(1,1:kk),y(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),r(1,1:kk),'k--','linewidth',1.5);
ylabel('$y(k)$','fontsize',14,'interpreter','latex')
xlabel('$t$ (s)','fontsize',14,'interpreter','latex')
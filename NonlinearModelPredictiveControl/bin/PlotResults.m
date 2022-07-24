figure(1);clf;
subplot(2,2,1)
stairs(ops.t(1,1:kk),x(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(1)*ones(1,kk),'r--');
ylabel('$x_1(k)$','fontsize',18,'interpreter','latex')
subplot(2,2,2)
stairs(ops.t(1,1:kk),x(2,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.x_max(2)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.x_min(2)*ones(1,kk),'r--');
ylabel('$x_2(k)$','fontsize',18,'interpreter','latex')
subplot(2,2,3)
stairs(ops.t(1,1:kk),u(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),ops.u_max(1)*ones(1,kk),'r--');
stairs(ops.t(1,1:kk),ops.u_min(1)*ones(1,kk),'r--');
ylabel('$u(k)$','fontsize',18,'interpreter','latex')
xlabel('$t$ (s)','fontsize',18,'interpreter','latex')
subplot(2,2,4)
stairs(ops.t(1,1:kk),y(1,1:kk),'linewidth',1.5);grid;hold on
stairs(ops.t(1,1:kk),r(1,1:kk),'k--','linewidth',1.5);
ylabel('$y(k)$','fontsize',18,'interpreter','latex')
xlabel('$t$ (s)','fontsize',18,'interpreter','latex')
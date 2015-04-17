function [] = plot_spy(paramf, out_str)

load(paramf);
figure(1);
spy(mrf_weights);
title('MRF param Matrix');
xlabel('visNodes');
ylabel('visNodes');
saveas(gcf, [out_str '_mrf.png'])
close(1);

load(paramf);
figure(2);
spy(L');
title('RBM param Matrix');
xlabel('visNodes');
ylabel('hidNodes');
saveas(gcf, [out_str '_rbm.png'])
close(2);





end


filein = 'net_test.out';

fid = fopen (filein,'r');
net_arr = [];
net_arr_cnt = 1;
tline = fgetl(fid);
while ischar(tline)
    [store rem] = strtok(tline);
    [store rem] = strtok(rem);
    % rem should be the wire length
    net_arr(net_arr_cnt) = str2num(rem);
    net_arr_cnt = net_arr_cnt + 1;
    tline = fgetl(fid);
end

fclose(fid);

net_arr_sorted = sortrows(net_arr);
histog_net = histogram(net_arr_sorted, 10, 'orientation', 'vertical');
hold on
grid on
ylabel('Net Count');
xlabel('Wire Length (um)');
% yticks([0 7 14 21 28 35 42 49 56 63 70 77 84 91 98 112 126]);
% xticks([0 400 800 1200 1]);
title('Wire Length Histogram (M3)');
% legend('Cortex M3','Location','east');
x0 = 0;
y0 = 0;
width = 400;
height = 250;
set(gcf, 'position', [x0,y0,width,height])

% area(net_arr_sorted, 'EdgeColor', "#7E2F8E", 'FaceColor', "#7E2F8E");
% title('Wire Length per Net');
% xlabel('Net Count');
% ylabel('Length (um)');
% % xlim([0 9074]);
% % yticks([0 200 400 600 800 1000 1200]);
% % ylim([0 1000]);
% grid on;
% hold on;
%multiple delays


Tmax = 6;
lags = [1; 2; 4];
N = 100;             %sample time points
Nl = length(lags);
Nc = 2*length(lags)+1; %number of components
spec = linspecer(Nc);

f = @(x) 1*(0.5 * cos(3*x) -  sin(5.5*x) + 0.3*sin(5*x) + sin(0.1*x) - cos(10*x));

%lag intervals
lag_rev = reshape(lags(end:-1:1), 1, []);

%all time lag intervals together
lag_all = [-lag_rev, 0, Tmax - lag_rev, Tmax];

%stack up the spans together
lag_span = [lag_all(1:end-1); lag_all(2:end)];

x = cell(Nc,1);
t = cell(Nc, 1);

for i = 1:Nc
    lag_curr = lag_span(:, i);
    t{i} = linspace(lag_curr(1), lag_curr(2), N*diff(lag_curr));
    x{i} = f(t{i});
end

%% Component Measures
figure(39)
clf
hold on
for i = 1:Nc
    plot(t{i}, x{i}, 'color', spec(i, :), 'LineWidth', 2)
end
% title('Component Decomposition', 'FontSize', FS)
plot([0, 0], ylim, ':k', 'LineWidth', 2)
ylabel('x(t)')
xlabel('t')

%% Component Support
figure(40)
clf
hold on
text_up = 0.25;
text_left_0 = 0.15;
for i = 1:Nc
    plot([lag_span(1, i), lag_span(2, i)], [0,0],'color', spec(i, :), 'LineWidth', 5) 
    
    if i <= length(lags)
        text_left = text_left_0 + 0.15;
    else
        text_left = text_left_0;
    end
    text(mean(lag_span(:, i))-text_left, text_up, ['$\nu_{', num2str(i-1-length(lags)), '}$'], ...
        'interpreter', 'latex', 'FontSize', 18)
end
ax1 = gca;                   % gca = get current axis
ax1.YAxis.Visible = 'off';  
ax1.XAxis.FontSize=12;
xlabel('time')
ylim([-0.5, 1])

%% Component Stack
figure(41)
clf
ax_list = cell(Nl+1, 1);
tiledlayout(Nl+1, 1);
for i = 0:Nl
%     subplot(Nl+1, 1, i+1); 
    ax_list{i+1} = nexttile;    
    hold on
    if i==0
        curr_lag = 0;        
    else
        curr_lag = lags(i);
    end
%     curr_lag = 12;
    for j = -i:(length(lags)-i)
        plot(t{j+Nl+1} + curr_lag, x{j+Nl+1}, 'color', spec(j+Nl+1, :), 'LineWidth', 2)
    end
        
       ylabel(['x(t-\tau_{', num2str(i), '})']);
                
end
xlabel('time')
linkaxes([ax_list{:}], 'xy')
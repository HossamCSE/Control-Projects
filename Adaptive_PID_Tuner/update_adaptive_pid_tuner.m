function update_adaptive_pid_tuner()
    %% Adaptive PID Tuner with GUI Input, Metrics, and Comparison Plot

    % Step 1: Prompt user for transfer function
    prompt = {'Enter Transfer Function (e.g., 1/((s+1)(s+2)) or (s+2)/((s+1)(s+5))*exp(-0.2*s)):'};
    dlgtitle = 'System Input';
    definput = {'1/((s+1)*(s+2))'}; 
    answer = inputdlg(prompt, dlgtitle, [1 70], definput);

    if isempty(answer)
        return;
    end
    tf_str = answer{1};

    % Step 2: Convert string to transfer function
    try
        G = parse_tf(tf_str);
    catch ME
        errordlg(['Invalid TF! ', ME.message], 'Error');
        return;
    end

    % Step 3: Initial PID using Skogestad method
    [Kp_init, Ki_init, Kd_init] = skogestad_tuning(G);

    % Step 4: Closed-loop response before tuning
    C_init = pid(Kp_init, Ki_init, Kd_init);
    sys_init = feedback(C_init * G, 1);
    [y_init, t_init] = step(sys_init, 50);

    % Step 5: Adaptive PID tuning
    alpha = 0.01;
    max_iter = 300;
    [Kp, Ki, Kd, y_opt, t_opt] = adaptive_pid(G, Kp_init, Ki_init, Kd_init, alpha, max_iter);

    % Step 6: Final optimized closed-loop system
    C_final = pid(Kp, Ki, Kd);
    sys_final = feedback(C_final * G, 1);

    % Step 7: Plot both responses
    figure;
    plot(t_init, y_init, '--', 'LineWidth', 1.5); hold on;
    plot(t_opt, y_opt, '-', 'LineWidth', 2);
    legend('Before Tuning', 'After Tuning');
    title('System Response: Before vs After PID Tuning');
    xlabel('Time (s)');
    ylabel('Output');
    grid on;

    % Step 8: Analyze metrics
    info = stepinfo(y_opt, t_opt);
    overshoot = info.Overshoot;
    settling_time = info.SettlingTime;
    ss_error = abs(1 - y_opt(end));

    % Step 9: Display results
    msg = sprintf(['Optimized PID Parameters:\nKp = %.4f\nKi = %.4f\nKd = %.4f\n\n' ...
                   'Overshoot: %.2f%%\nSettling Time: %.2f s\nSteady-State Error: %.4f'], ...
                   Kp, Ki, Kd, overshoot, settling_time, ss_error);
    msgbox(msg, 'Optimized Results');
end

%% --- Parse Transfer Function String with Delay Support ---
function sys = parse_tf(str)
    s = sym('s');
    str = strrep(str, ' ', '');
    str = regexprep(str, '\)\(', ')*(');
    str = regexprep(str, ')(exp', ')*(exp');
    str = regexprep(str, '\)\s*(s)', ')*$1');

    delay = 0;
    delay_match = regexp(str, 'exp\(-([0-9\.]+)\*s\)', 'tokens');
    if ~isempty(delay_match)
        delay = str2double(delay_match{1}{1});
        str = regexprep(str, '\*?exp\(-[0-9\.]+\*s\)', '');
    end

    try
        expr = evalin(symengine, str);
    catch
        error('Invalid TF! Expression could not be evaluated.');
    end

    [num_sym, den_sym] = numden(expr);
    num = sym2poly(num_sym);
    den = sym2poly(den_sym);

    if isempty(num) || isempty(den)
        error('Invalid TF! Could not extract numerator/denominator.');
    end

    sys = tf(num, den);
    if delay > 0
        sys.InputDelay = delay;
    end
end

%% --- Skogestad Initial PID Tuning ---
function [Kp, Ki, Kd] = skogestad_tuning(G)
    [y, t] = step(G);
    K = dcgain(G);

    if isinf(K) || K == 0
        error('DC gain is zero or infinite.');
    end

    y_norm = y / K;
    idx_2 = find(y_norm >= 0.02, 1);
    idx_63 = find(y_norm >= 0.63, 1);

    if isempty(idx_2) || isempty(idx_63)
        error('Cannot estimate delay and time constant.');
    end

    theta = t(idx_2);
    tau = t(idx_63) - theta;

    Kp = (1 / K) * (0.5 * tau) / (theta + 1e-3);
    Ti = min(tau, 4 * (theta + 1e-3));
    Ki = Kp / Ti;
    Kd = 0;
end

%% --- Adaptive PID Optimization ---
function [Kp, Ki, Kd, y, t] = adaptive_pid(G, Kp, Ki, Kd, alpha, max_iter)
    bounds = struct('Kp', [0.01, 100], 'Ki', [0.001, 50], 'Kd', [0, 20]);
    decay = 0.99;
    min_alpha = 1e-4;
    prev_cost = Inf;

    for iter = 1:max_iter
        [J_curr, grad] = gradients(Kp, Ki, Kd);
        Kp = constrain(Kp - alpha * grad.Kp, bounds.Kp);
        Ki = constrain(Ki - alpha * grad.Ki, bounds.Ki);
        Kd = constrain(Kd - alpha * grad.Kd, bounds.Kd);

        alpha = max(min_alpha, alpha * decay);
        if iter > 10 && abs(J_curr - prev_cost) < 1e-6
            break;
        end
        prev_cost = J_curr;
    end

    C = pid(Kp, Ki, Kd);
    sys_cl = feedback(C * G, 1);
    [y, t] = step(sys_cl, 50);

    % Nested cost and gradient functions
    function [J, grad] = gradients(Kp_, Ki_, Kd_)
        delta = 1e-3;
        grad.Kp = (cost(Kp_ + delta, Ki_, Kd_) - cost(Kp_ - delta, Ki_, Kd_)) / (2 * delta);
        grad.Ki = (cost(Kp_, Ki_ + delta, Kd_) - cost(Kp_, Ki_ - delta, Kd_)) / (2 * delta);
        grad.Kd = (cost(Kp_, Ki_, Kd_ + delta) - cost(Kp_, Ki_, Kd_ - delta)) / (2 * delta);
        J = cost(Kp_, Ki_, Kd_);
    end

    function J = cost(Kp_, Ki_, Kd_)
        C_ = pid(Kp_, Ki_, Kd_);
        sys_ = feedback(C_ * G, 1);
        [y_, t_] = step(sys_, 50);
        J = trapz(t_, t_ .* abs(1 - y_));
    end

    function val = constrain(x, range)
        val = max(min(x, range(2)), range(1));
    end
end
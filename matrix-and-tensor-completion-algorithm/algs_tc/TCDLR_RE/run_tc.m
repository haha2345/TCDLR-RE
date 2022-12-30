function [Xhat,  err, iter, rank_e] = run_tc(params_tc)
    disp(['#######FTCN########']);

    T = params_tc.T;
    Idx = params_tc.Idx;
    Omega=find(Idx==1);
    normalize = max(T(:));
%     X = T / normalize;
    X = T*255;
    
    [n1, n2, n3] = size(T);
%     opts.mu = 1e-5;
%     opts.rho = 1.3;
%     opts.max_iter = 50;
%     opts.lp = 0.8;
%     rate = 0.4;

    opts.mu = params_tc.mu;
    opts.rho = params_tc.rho;
    opts.max_iter = params_tc.max_iter;
    opts.lp = params_tc.lp;
    rate = params_tc.rate;

    opts.DEBUG = 1;
%     opts.veck = 30 * ones(n3);
        opts.veck = 2 * ones(n3);

%     opts.rank_min = 10 * ones(n3);
        opts.rank_min = 1 * ones(n3);

%     opts.rank_max = floor(rate * min(size(X, 1), size(X, 2))) * ones(n3);
    opts.rank_max = 3 * ones(n3);

%     maxP = max(X(:));

    data = X(Omega(:));
    Xn = zeros(size(X));
    Xn(Omega(:)) = data;
%     tic
    [Xhat, err, iter, rank_e] = lrtc_tnn_fastNEW(Xn, Omega, X, opts);
    Xhat(Omega(:))=data;
%     time = toc;
%     psnr_result = PSNR(X, Xhat, maxP);
%     ssim_result = ssim(Xhat, X);
%     fsim_result = FeatureSIM(X, Xhat);

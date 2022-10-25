function [ mi, phase ] = get_mi( phase_lf, amp_hf, method, varargin )
%GET_MI computes the modulation index that measures the phase-amplitude
%coulping between low and high frequency components.
%   mi= get_mi( phase_lf, amp_hf, method, parameters )
%   phase_lf-- phase time series
%   amp_hf -- amplitude time series
%   method -- 'entropy', 'entropy_stats', 'mean_vector', 'direct_pac'
%   parameters -- for the entropy case, you can specify as first parameter the
%   size of the phase bins. For the 'entropy_stats' you can specify the
%   number of surrogates for computing the modulation index.
%
%   Authors: David Escobar, method entropy_stats based on code by Allison
%   Connolly.

%%

mi=0;
phase=0;
n_args = length(varargin);
switch method
    
    case 'mean_vector'
        vec = amp_hf.*exp(1i*phase_lf ) ;
        vec_sum=  sum(vec')/size(vec,2) ;
        phase = angle(vec_sum);
        mi = abs(vec_sum)/max(amp_hf);
        
    case 'direct_pac'
        vec_sum = sum (amp_hf.*exp(1i.*phase_lf ) );
        phase = angle(vec_sum) ;
        mi = ( norm ( vec_sum ) / sqrt(sum(amp_hf.^2)) ) / sqrt(size(amp_hf,2)) ;
        
        
    case 'entropy'
        d_ph=20;
        if n_args == 1
            n_bins=varargin{1};
            d_ph=2*pi/n_bins;
        end
        
        %amp= amp_hf';
        pha = wrapTo2Pi(phase_lf);
        pha= pha';
                
        % calc distribution of mean amp for bins of phases
        pha_bin_indexes = ceil( (pha/(2*pi))*n_bins  );
        amp_bin = zeros(n_bins,1);
        for k=1:n_bins
            amp_bin(k) = mean( amp_hf( pha_bin_indexes == k)',1);
            %amp_bin(k) = ( amp_hf( pha_bin_indexes == k)',1);
        end
        
        P = amp_bin/sum(amp_bin) + 1e-12;
        N= length(P);
        H = -sum( P.*log(P)) ;
        mi = (log(N)  - H)/log(N); % modulation index
        [~,ind_ph] = max(P);
        phase = ( d_ph )*(ind_ph-1) + d_ph/2;
        phase = wrapToPi(phase);
        
    case 'entropy_stats' % Code adapted from A. Connolly
        
        
        amp= amp_hf';
        
        pha = wrapTo2Pi(phase_lf);
        pha= pha';
        
        pha = pha';
        if n_args >= 1
            numbins=varargin{1};
            d_ph=2*pi/numbins;
        end
        
        if n_args > 1
            NSurr=varargin{2};
        else
            NSurr= 200;
        end
        
        % calc distribution of mean amp for bins of phases
        pha_binned = ceil(pha/(2*pi)*numbins);
        distr = zeros(numbins,1);
        for kk=1:numbins
            distr(kk) = mean( amp(pha_binned == kk),1);
        end
        %entropy of the distribution
        p = distr/sum(distr);
        H = -sum(p.*log2(p));
        [~,ind_ph] = max(p);
        phase = ( d_ph )*(ind_ph-1) + d_ph/2;
        phase = wrapToPi(phase);
        
        Hmax=log2(length(distr));
        %mi=(Hmax-H)/Hmax;
        
        HSurr = zeros(NSurr,1);
        parfor SS = 1:NSurr
            distr = zeros(numbins,1);
            randshift = ceil(rand(1)*size(pha,1)) 
            temp = [pha_binned(randshift:end); pha_binned(1:randshift-1)];
            for kk = 1:numbins
                distr(kk) = mean(amp(temp == kk),1);
            end
            pSurr = distr/sum(distr);
            HSurr(SS) = -sum(pSurr.*log2(pSurr));
        end
        
        mi = (mean(HSurr)-H)/std(HSurr);
        %phase = max(HSUrr)
        %}
        
        
    case 'entropy2' % weighted entropy - 
        d_ph=20;
        kappa =0.8;
        if n_args >= 1
            n_bins=varargin{1};
            d_ph=2*pi/n_bins;
        end
        
        amp= amp_hf';
        pha = phase_lf + pi;
        pha = pha';
        
        % calc distribution of mean amp for bins of phases
        pha_bin_indexes = ceil( (pha/(2*pi))*n_bins  );
        amp_bin = zeros(n_bins,1);
        for k=1:n_bins
            amp_bin(k) = mean( amp( pha_bin_indexes == k),1);
        end
        
        P = amp_bin/sum(amp_bin);
        
        N= length(P);
        H = -sum( P.*log(P)) ;
        D_kl = log(N)  - H; % KL distance to uniform distributions
        mi_1 = D_kl/log(N); % modulation index
        
        mi_2 = max(amp_bin);
        
        mi = (mi_1^(kappa))*(mi_2^(1-kappa));
        
        [~,ind_ph] = max(P);
        phase = ( d_ph )*(ind_ph-1) + d_ph/2 - pi;
        
    otherwise
        disp('method not valid')
        
end

if isnan(mi)
    mi =0;
end

end

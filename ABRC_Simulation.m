function output = ABRC_Simulation(niter)
% usage:   output = ABRC_Simulation(niter)
% author:  Scott Schafer
% date:    10/16/2014
% purpose: This function is designed to perform a power calculation for
%          detecting a difference between patients and controls for
%          BOLD/RSA covariation, presuming N = 10 controls, N = 7 patients,
%          and four samples of ~600 images each.
%
%          Input variables:
%               niter = The number of iterations to simulate the data.
%                       Higher iterations will be more accurate but will take
%                       longer to run. 
%
%          Output variables:
%             output =  A struct that contains the pertinent data
%                       regarding the simulation
%                       
%                       output.pval is a vector of p-values for the linear slope X patient group interaction for each simulation
%                       output.sigper is the percentage of p-values less than or equal to .05
%                       output.meanb is a vector of slopes for patients and controls for each simulation
%
%          Example: output = ABRC_Simulation(3000)
%
%                   ====== Within Subjects ======
% 
%                   For a within-subject effect size of 0.300, between session variance of 0.0016, and between subjects variance of 0.0000, the t-stat was significant  76% of the time
% 
%                   ====== Between Subjects ======
% 
%                   For a within-subject effect size of 0.300, between session variance of 0.0016, and between subjects variance of 0.0030, the t-stat was significant  37% of the time


% rows are weeks (Wk0, Wk2, Wk6, Wk12), columns are MDD, Control. Values
% are the true levels of RSA-BOLD covariation within an ROI on average
% across subjects. The values will have random subject errors applied to
% them.

%if these were the true values within all subjects,
truval = [-.025 .05
              0 .05
           .025 .05
            .05 .05];

%use these values to estimate errors (Model should detect difference 5% of the time)
% truval = [.05 .05
%           .05 .05
%           .05 .05
%           .05 .05];

% This is the variance of the RSA-BOLD covariation across weeks
truvar = var(truval(:,1));

eta2 = 0.3; %this is the within ROI effect size (related to noise variance)
variance = truvar./eta2 - truvar; %this is the within-subject noise variance, added on top of the signal

nimgs = 600; %this is the approximate number of times the data in the ROI is sampled each session

week = 1:4;
N = [7 10];

bet_sess_var = .0016; %this is the session variance within subjects (estimated from controls)
bet_subj_var = .003; %this is the variance in scores between subjects (estimated from controls)


%First do the within-subject case, then the between subject case
fprintf('\n====== Within Subjects ======\n\n');
within_output = runmodel(bet_sess_var, 0,'Within Subjects');
fprintf('\n====== Between Subjects ======\n\n');
between_output = runmodel(bet_sess_var, bet_subj_var, 'Between Subjects');

output = [within_output between_output];



    function output = runmodel(bet_sess_var, bet_subj_var,titlename) 
        
        for vidx = 1:numel(bet_subj_var)
            for bidx = 1:numel(bet_sess_var) %over all between session variances
                for pidx = 1:numel(eta2) %over all within-subject power levels
                    clear pval meanb
                    for i = 1:niter %over all model iterations
                        clear mdddata condata
                        for j = week %over all weeks
                            clear ROIvals
                            for k = 1:N(1) % get ROI values for each MDD patient within each of the images
                                ROIvals(:,k) = (truval(j,1)+randn*sqrt(bet_sess_var(bidx))+randn*sqrt(bet_subj_var(vidx)))*ones(nimgs,1)+randn(nimgs,1)*sqrt(variance(pidx));
                            end
                            mdddata(j,:) = mean(ROIvals); %average values across the images
                            clear ROIvals
                            for k = 1:N(2) %get ROI values for each control within each of the images
                                ROIvals(:,k) = (truval(j,2)+randn*sqrt(bet_sess_var(bidx))+randn*sqrt(bet_subj_var(vidx)))*ones(nimgs,1)+randn(nimgs,1)*sqrt(variance(pidx));
                            end
                            condata(j,:) = mean(ROIvals); %average values across the images
                        end
                        if bet_subj_var == 0 %within subjects
                            clear b
                            for k = 1:N(1) %calculate the within-subject slopes across week
                                b(:,k) = regress(mdddata(:,k),[[1 1 1 1];[-3 -1 1 3];[-1 1 1 -1];[-1 3 -3 1]]');
                            end
                            mdd_b = b(2,:); %We only care about the linear slope
                            %                         keyboard
                            clear b
                            for k = 1:N(2) %calculate the within-subject slopes across week
                                b(:,k) = regress(condata(:,k),[[1 1 1 1];[-3 -1 1 3];[-1 1 1 -1];[-1 3 -3 1]]');
                            end
                            con_b = b(2,:); %We only care about the linear slope
                            [~,p] = ttest2(mdd_b,con_b); %get p for a between groups t-test on the within subject slopes
                            meanb(i,:) = [mean(mdd_b) mean(con_b)];
                        else %between subjects
                            mdd_val = reshape(mdddata,N(1)*4,1);
                            con_val = reshape(condata,N(2)*4,1);
                            mddvcon = [ones(N(1)*4,1); -1*ones(N(2)*4,1)];
                            contrasts = [-3 -1 -1;
                                         -1 1 3;
                                         1 1 -3;
                                         3 -1 1];
                            contrasts = repmat(contrasts,N(1)+N(2),1);
                            [b,~,stats] = glmfit([mddvcon contrasts contrasts.*repmat(mddvcon,1,3)],[mdd_val;con_val]);
                            %get p-value for linslope X mddvcon interaction
                            p = stats.p(6);
                            %slopes are mean linslope + linslopeXmddvcon for patients and linslope - linslopXmddvcon for controls
                            meanb(i,:) = [stats.beta(3)+stats.beta(6) stats.beta(3)-stats.beta(6)];
                        end
                        pval(i,1) = p;
                        
                    end
                    test = pval <= .05; %How many p's are less than or equal to .05?
                    final = sum(test)/numel(test)*100; % turn it into a percentage
                    fprintf('For a within-subject effect size of %4.3f, between session variance of %5.4f, and between subjects variance of %5.4f, the t-stat was significant %3.0f%% of the time\n',eta2(pidx),bet_sess_var(bidx),bet_subj_var(vidx),final)
                    output(bidx,pidx).pval = pval;
                    output(bidx,pidx).sigper = final;
                    output(bidx,pidx).meanb = meanb;
                    
                    %The following script generates a histogram
                    %to view each simulation
                    figure
                    hist(meanb(:,1),60);
                    hold on
                    hist(meanb(:,2),60);
                    h = findobj(gca,'Type','patch');
                    set(h(1),'FaceColor','r','EdgeColor','k')
                    title(titlename,'FontSize',14)
                    ylabel('Count','FontSize',14)
                    xlabel('Mean beta value','FontSize',14)
                    legend({'MDD' 'Control'},'Location','Best')
                    box off
                    set(gca,'LineWidth',2,'FontSize',12)
                end
            end
        end
        
    end

end






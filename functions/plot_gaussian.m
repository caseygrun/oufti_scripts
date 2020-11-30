
function [percent_on, model, CDF_figure, PDF_figure] = plot_gaussian(mfis, show_plot, n_bins)
% PLOT_GAUSSIAN(MFIS) fit a Gaussian mixture model to the MFI distribution  
%    and plot the results over a histogram
%
%    Arguments:
%        MFIS : row vector of MFIs of cells in a frame
%        SHOW_PLOT : true to show histogram, false to skip (default: true)
%        N_BINS : number of bins to use for histogram; defaults to
%        sqrt(length(mfis))
%        
%    Returns:
%        percent_on : % of cells in the list that are classified as "ON" by
%        the model
%        model : the Gaussian mixture mode, fit by

	if iscolumn(mfis)
		mfis = mfis';
    end
    CDF_figure = [];
    PDF_figure = [];
    
	mdl = fitgmdist(mfis',2,'CovarianceType','diagonal');
    model = mdl;
    
    % smallcomp is the low-mean component of the mixture (e.g. the OFF
    % population)
	smallcomp=find(mdl.mu==min(mdl.mu));
	largecomp=find(mdl.mu==max(mdl.mu));
    percent_on = mdl.ComponentProportion(largecomp);
    
    if ~exist('show_plot','var')
        show_plot = 1;
    end
    
    if ~show_plot
        return
    end
    
    
    if ~exist('n_bins','var')
        binnum=ceil(sqrt(length(mfis)));
    else
        binnum = n_bins
    end
	
	xaxis=linspace(min(mfis),max(mfis),binnum);

	% Plot CDF
	CDF_figure = figure();
	xlabel('Mean Fluorescence Intensity') %minus background
	ylabel('Cumulative Probability Function')
	title('CDF evaluated at each data point')

	hold on
	[f,x]=ecdf(mfis);
	cdfsum=cdf(mdl,x);
	plot(x,cdfsum,'-g')
	plot(x,f,'.k')
	xax=linspace(min(x),max(x),100);

	cdf1=mdl.ComponentProportion(smallcomp)*normcdf(xax,mdl.mu(smallcomp),sqrt(mdl.Sigma(smallcomp)));
	cdf2=mdl.ComponentProportion(largecomp)*normcdf(xax,mdl.mu(largecomp),sqrt(mdl.Sigma(largecomp)));
	plot(xax,cdf1,'--b')
	plot(xax,cdf2+mdl.ComponentProportion(smallcomp),'--r')
	legend({'GMM fit','data','CDF1','CDF2+w1'})

	% Plot histograms of GFP MFIs
	PDF_figure = figure();
	[n, c]=ecdfhist(f,x,binnum);
	
	hold on;
	plot(c,n,'ok')
	nfit=ecdfhist(cdfsum,x,binnum);
	plot(c,nfit,'-g')
	xlabel('Mean Fluorescence Intensity')
	ylabel('Probability Density')
	pdf1=ecdfhist(cdf1,xax,binnum);
	pdf2=ecdfhist(cdf2,xax,binnum);

	% plot mixed model fits
	plot(c,pdf1,'--b')
	plot(c,pdf2,'--r')
	legend({'combined data',...
        'EM Fit',...
        ['c_{1}=' num2str(mdl.ComponentProportion(smallcomp))...
            ' \mu_{1}=' num2str(mdl.mu(smallcomp)) ' \sigma_{1}=' num2str(sqrt(mdl.Sigma(smallcomp))) ],...
        ['c_{2}=' num2str(mdl.ComponentProportion(largecomp))... 
            ' \mu_{2}=' num2str(mdl.mu(largecomp)) ' \sigma_{2}=' num2str(sqrt(mdl.Sigma(largecomp))) ]...
    })
	title(['PDF evaluated with N = ' num2str(binnum) 'bins'])
    hold off;
end
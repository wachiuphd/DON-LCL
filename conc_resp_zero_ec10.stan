// Concentration-response downward to 0 (three parameter Hill model)
data {
	real scale_factor; // overall scaling factor
	int<lower=0> Ni; // Number of cell lines
	int<lower=0> Nj; // Number of data points
	vector[Nj] x;    // concentrations
	vector[Nj] ys;    // responses
	int<lower=0,upper=Ni> cell[Nj]; // cell line for each data point
	int<lower=1> Nquants;	// Number of quantiles of EC10 to calculate
	vector[Nquants] quants;	// Quantiles (e.g., c(0.01,0.025,0.5,0.975,0.99))
}
parameters { // all natural log-transformed
	// Population mean 
	real m_y0;	// background
	real m_x10;	// numerator log scale EC10
	real<lower=-2,upper=2> m_n;	// Hill exponent
	// Population standard deviation
	real<lower=0> sd_y0;
	real<lower=0> sd_x10; // SD log scale EC10 
	real<lower=0> sd_n;
	// Inter-individual variability
	real z_y0[Ni];
	real z_x10[Ni];
	real z_n[Ni];
	 // Residual error variance
	real<lower=0> sigma_y;
}
transformed parameters {
	real y0[Ni];	// untransformed background
	real x10[Ni];	// untransformed EC10 scale
	real n[Ni];		// untransformed Hill exponent
	real x0[Ni];	// untransformed numerator scale
	for (i in 1:Ni) { // un-log transform
		y0[i] = exp(m_y0 + sd_y0 * z_y0[i]);
		x10[i] = exp(m_x10 + sd_x10 * z_x10[i]);
		n[i] = exp(m_n + sd_n * z_n[i]);
		x0[i] = x10[i] * 9.0^(1/n[i]); // original scale paramter
	}
}
model {
	vector[Nj] yp;
	// prior distributions
	m_y0 ~ normal(0,5);
	m_x10 ~ normal(0,5);
	m_n ~ normal(0,1);
	sd_y0 ~ normal(0,1);
	sd_x10 ~ normal(0,1);
	sd_n ~ normal(0,0.2);
	z_y0 ~ normal(0,1);
	z_x10 ~ normal(0,1);
	z_n ~ normal(0,1);
	sigma_y ~ normal(0,0.2); 
	for (j in 1:Nj) {
		yp[j] = y0[cell[j]]/ (1 + (x[j] / x0[cell[j]])^n[cell[j]]);
	}
	ys ~ student_t(5,yp,sigma_y);
}
generated quantities {
	vector[Ni] ec99;		// Concentration at which the relative response is 99%
	vector[Ni] ec975;		// Concentration at which the relative response is 97.5%
	vector[Ni] ec95;		// Concentration at which the relative response is 95%
	vector[Ni] ec90;		// Concentration at which the relative response is 90%
	vector[Ni] ec50;		// Concentration at which the relative response is 50%
	vector[Ni] ec10;		// Concentration at which the relative response is 10% 
	vector[Ni] ec05; 		// Concentration at which the relative response is 5%
	vector[Ni] ec025;		// Concentration at which the relative response is 2.5%
	vector[Ni] ec01;		// Concentration at which the relative response is 1%
	
	// real ec99_median;			// Median ec99
	// real ec975_median;			// Median ec975
	// real ec95_median;			// Median ec95
	// real ec90_median;			// Median ec90
	// real ec50_median;			// Median ec50
	// real ec10_median;			// Median ec10
	// real ec05_median;			// Median ec05
	// real ec025_median;			// Median ec025
	// real ec01_median;			// Median ec01
	// 
	// vector[Nquants] ec99_quants;	// Quantiles ec99
	// vector[Nquants] ec99_quant_ratios;	// Ratios of median to quantiles ec99
	// vector[Nquants] ec975_quants;	// Quantiles ec975
	// vector[Nquants] ec975_quant_ratios; // Ratios of median to quantiles ec975
	// vector[Nquants] ec95_quants;	// Quantiles ec95
	// vector[Nquants] ec95_quant_ratios; // Ratios of median to quantiles ec95
	// vector[Nquants] ec90_quants;	// Quantiles ec90
	// vector[Nquants] ec90_quant_ratios; // Ratios of median to quantiles ec90
	// vector[Nquants] ec50_quants;	// Quantiles ec50
	// vector[Nquants] ec50_quant_ratios; // Ratios of median to quantiles ec50
	// vector[Nquants] ec10_quants;	// Quantiles 
	// vector[Nquants] ec10_quant_ratios;	// Ratios of median to quantiles
	// vector[Nquants] ec05_quants;	// Quantiles ec05
	// vector[Nquants] ec05_quant_ratios; // Ratios of median to quantiles ec05
	// vector[Nquants] ec025_quants;	// Quantiles ec02.5
	// vector[Nquants] ec025_quant_ratios;	// Ratios of median to quantiles ec02.5
	// vector[Nquants] ec01_quants;	// Quantiles ec01
	// vector[Nquants] ec01_quant_ratios; // Ratios of median to quantiles ec01
	
	
	for (i in 1:Ni) {
		ec99[i] = x0[i]*(0.99 / (1 - 0.99))^(1/n[i]);
		ec975[i] = x0[i]*(0.975 / (1 - 0.975))^(1/n[i]);
		ec95[i] = x0[i]*(0.95 / (1 - 0.95))^(1/n[i]);
		ec90[i] = x0[i]*(0.90 / (1 - 0.90))^(1/n[i]);
		ec50[i] = x0[i]*(0.5 / (1 - 0.5))^(1/n[i]);
		ec10[i] = x0[i]*(0.1 / (1 - 0.1))^(1/n[i]);
		ec05[i] = x0[i]*(0.05 / (1 - 0.05))^(1/n[i]);
		ec025[i] = x0[i]*(0.025 / (1 - 0.025))^(1/n[i]);
		ec01[i] = x0[i]*(0.01 / (1 - 0.01))^(1/n[i]);
	}
	// {
	// int ec99_indx[Ni];
	// int ec975_indx[Ni];
	// int ec95_indx[Ni];
	// int ec90_indx[Ni];
	// int ec50_indx[Ni];
	// int ec10_indx[Ni];
	// int ec05_indx[Ni];
	// int ec025_indx[Ni];
	// int ec01_indx[Ni];
	// real p;
	// real m;
	// real g;
	// int i;
	// int itmp;
	// ec99_indx= sort_indices_asc(ec99);
	// ec975_indx= sort_indices_asc(ec975);
	// ec95_indx= sort_indices_asc(ec95);
	// ec90_indx= sort_indices_asc(ec90);
	// ec50_indx= sort_indices_asc(ec50);
	// ec10_indx= sort_indices_asc(ec10);
	// ec05_indx= sort_indices_asc(ec05);
	// ec025_indx= sort_indices_asc(ec025);
	// ec01_indx= sort_indices_asc(ec01);
	// p= 0.5;
	// m=1-p;
	// itmp=1;
	// while (itmp <(Ni*p + m))
	// 	itmp=itmp+1;
	// i = itmp-1;
	// g = Ni*p + m - i;
	// ec99_median = (1-g)*ec99[ec99_indx[i]]+g*ec99[ec99_indx[i+1]];
	// ec975_median = (1-g)*ec975[ec975_indx[i]]+g*ec975[ec975_indx[i+1]];
	// ec95_median = (1-g)*ec95[ec95_indx[i]]+g*ec95[ec95_indx[i+1]];
	// ec90_median = (1-g)*ec90[ec90_indx[i]]+g*ec90[ec90_indx[i+1]];
	// ec50_median = (1-g)*ec50[ec50_indx[i]]+g*ec50[ec50_indx[i+1]];
	// ec10_median = (1-g)*ec10[ec10_indx[i]]+g*ec10[ec10_indx[i+1]];
	// ec05_median = (1-g)*ec05[ec05_indx[i]]+g*ec05[ec05_indx[i+1]];
	// ec025_median = (1-g)*ec025[ec025_indx[i]]+g*ec025[ec025_indx[i+1]];
	// ec01_median = (1-g)*ec01[ec01_indx[i]]+g*ec01[ec01_indx[i+1]];
	// for (q in 1:Nquants) {
	// 	p= quants[q];
	// 	m= 1-p;
	// 	itmp=1;
	// 	while (itmp <(Ni*p + m)) 
	// 		itmp=itmp+1;
	// 	i = itmp-1;
	// 	g = Ni*p + m - i;
	// 	ec99_quants[q] = (1-g)*ec99[ec99_indx[i]]+g*ec99[ec99_indx[i+1]];
	// 	ec975_quants[q] = (1-g)*ec975[ec975_indx[i]]+g*ec975[ec975_indx[i+1]];
	// 	ec95_quants[q] = (1-g)*ec95[ec95_indx[i]]+g*ec95[ec95_indx[i+1]];
	// 	ec90_quants[q] = (1-g)*ec90[ec90_indx[i]]+g*ec90[ec90_indx[i+1]];
	// 	ec50_quants[q] = (1-g)*ec50[ec50_indx[i]]+g*ec50[ec50_indx[i+1]];
	// 	ec10_quants[q] = (1-g)*ec10[ec10_indx[i]]+g*ec10[ec10_indx[i+1]];
	// 	ec05_quants[q] = (1-g)*ec05[ec05_indx[i]]+g*ec05[ec05_indx[i+1]];
	// 	ec025_quants[q] = (1-g)*ec025[ec025_indx[i]]+g*ec025[ec025_indx[i+1]];
	// 	ec01_quants[q] = (1-g)*ec01[ec01_indx[i]]+g*ec01[ec01_indx[i+1]];
	// }
	// for (q in 1:Nquants) { 
	// 	ec99_quant_ratios[q] = ec99_median/ec99_quants[q];
	// 	ec975_quant_ratios[q] = ec975_median/ec975_quants[q];
	// 	ec95_quant_ratios[q] = ec95_median/ec95_quants[q];
	// 	ec90_quant_ratios[q] = ec90_median/ec90_quants[q];
	// 	ec50_quant_ratios[q] = ec50_median/ec50_quants[q];
	// 	ec10_quant_ratios[q] = ec10_median/ec10_quants[q];
	// 	ec05_quant_ratios[q] = ec05_median/ec05_quants[q];
	// 	ec025_quant_ratios[q] = ec025_median/ec025_quants[q];
	// 	ec01_quant_ratios[q] = ec01_median/ec01_quants[q];
	// }
	// }
}

# gxp-sim

[Vignette](https://nickwisniewski.com/gxp-sim/vignette.html)

## Simulating gene microarrays simply

We set out to simulate gene expression arrays, in order to later investigate the effects of various biases on methods like clustering and PCA. Typical benchmarking or validation studies seek to generate “realistic” models of gene expression data, by simulating Michaelis-Mentin kinetics etc. (e.g. see the DREAM challenges). However, these simulations are much too sophisticated for our purposes, and their complexity may actually obfuscate some of the simple biases we would like to investigate. In contrast, we focused on creating a simple and understandable generative model that incorporates only the following idealized features:

* Mean expression levels (each gene has its own expression level and variance, and we can control the distribution of each)
* Sample “read coverage” or similar (each sample has its own mean expression level and variance, and we can control the distribution)

* Differential expression (mean differences between treatment groups)
* Gene co-expression (modular structure)
* Sample correlation (batch structure, treatment structure)

We especially note that the last bullet point concerning sample correlation often goes underappreciated. It is typically assumed that samples are independent. But correlation between samples reduces the effective sample size, while also inducing correlation between genes – both important details when using gene networks for statistical inference. For more details on this topic, see Efron (2009).

In this study, we designed a generative model for the purpose of investigating these features on gene network inference. 

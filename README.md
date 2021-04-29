<!--[![Build status](https://ci.appveyor.com/api/projects/status/7po303l6lv4fd18a?svg=true)](https://ci.appveyor.com/project/tkimhofer/jres)-->

# jres

R package for importing, pre-processing and visualising J-resolved (jres) Magnetic Resonance Spectroscopy (MRS) experiments. 

### functionality 
In the field of small molecule research (metabolomics/metabonomics) the interpretation of 1D MRS features can be tricky when there is a high degree of peak overlap/convolution. The visualiser is designed to address this issue by combining the standard 1D spectrum with a J-resolved spectrum (a fast 2D experiment that is standardly applied in high-throughput settings) acquired from the same sample. The 2D experiment is naturally less convoluted and the joint visualisation with the 1D experiment improves spectral interpretation and metabolite annotation.

### extensions
Functions also work for 2D MRS experiments other than jres.

### installation
```R
#install.packages('devtools')
devtools::install_github('tkimhofer/jres')
```






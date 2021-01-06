[![Build status](https://ci.appveyor.com/api/projects/status/7po303l6lv4fd18a?svg=true)](https://ci.appveyor.com/project/tkimhofer/jres)

# jres

R package for importing, pre-processing and visualising J-resolved (jres) Magnetic Resonance Spectroscopy (MRS) experiments. 

### functionality 
In the field of small molecule analysis (metabolomics/metabonomics) the annotation of 1D MRS features can be tricky in situations where there is a high degree of peak overlap. The visualiser is designed to to improve data interpretation and metaoblite annotaion by combining a J-resolved spectrum (a fast 2D experiment that is standardly applied in high-throughput settings) with the respective standard 1D spectrum acquired from the same sample.

### extensions
Fucntions can also be used for other 2D MRS experiments than jres.

##€ installation
```R
#install.packages('devtools')
devtools::install_github('tkimhofer/jres')
```






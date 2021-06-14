Rscript --vanilla 1_fit_variational.R --dir "decrease/vb/popdensity/naked" --intervention decrease --no_temporal --no_spatial --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "decrease/vb/popdensity/no_temporal" --intervention decrease --no_temporal --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "decrease/vb/popdensity/temporal05" --intervention decrease --autocor 0.5 --ar_scale 0.05 --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "decrease/vb/popdensity/temporal10" --intervention decrease --autocor 0.5 --ar_scale 0.1 --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "decrease/vb/popdensity/temporal15" --intervention decrease --autocor 0.5 --ar_scale 0.15 --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "stayhome/vb/popdensity/naked" --intervention stayhome --no_temporal --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "stayhome/vb/popdensity/no_temporal" --intervention stayhome --no_temporal --no_spatial --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "stayhome/vb/popdensity/temporal05" --intervention stayhome --autocor 0.5 --ar_scale 0.05 --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "stayhome/vb/popdensity/temporal10" --intervention stayhome --autocor 0.5 --ar_scale 0.1 --iter 50000  
Rscript --vanilla 1_fit_variational.R --dir "stayhome/vb/popdensity/temporal15" --intervention stayhome --autocor 0.5 --ar_scale 0.15 --iter 50000  

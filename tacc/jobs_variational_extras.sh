Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/stayhome/rand_lag" --bent_cable --intervention stayhome --no_temporal 
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/decrease/rand_lag" --bent_cable --intervention decrease --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/stayhome/bent_cable" --bent_cable --intervention stayhome --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/decrease/bent_cable" --bent_cable --intervention decrease --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/stayhome/bent_cable" --bent_cable --intervention stayhome --autocor 0.5 --ar_scale 0.1
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/decrease/bent_cable" --bent_cable --intervention stayhome --autocor 0.5 --ar_scale 0.1
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/stayhome/no_ny_no_temporal" --exclude_ny --intervention stayhome --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/decrease/no_ny_no_temporal" --exclude_ny --intervention decrease --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/stayhome/no_cities_temporal" --exclude_cities_post_inter --intervention stayhome --no_temporal
Rscript --vanilla 1_fit_variational.R --iter 20000 --dir "vb_extra/decrease/no_cities_temporal" --exclude_cities_post_inter --intervention decrease --no_temporal

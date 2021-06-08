Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention decrease --ar_relax_prior_scale --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention decrease --ar_relax_prior_scale --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention stayhome --ar_relax_prior_scale --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention stayhome --ar_relax_prior_scale --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention decrease --ar_relax_prior_scale --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention decrease --ar_relax_prior_scale --autocor 1.0 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention stayhome --ar_relax_prior_scale --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_relax/$LAUNCHER_TSK_ID" --intervention stayhome --ar_relax_prior_scale --autocor 1.0 --no_post_inter

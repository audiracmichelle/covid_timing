Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --intervention stayhome --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --intervention stayhome --autocor 1.0 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --exclude_post_only --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --exclude_post_only --intervention stayhome --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --exclude_post_only --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb_no_ny/$LAUNCHER_TSK_ID" --exclude_ny --exclude_post_only --intervention stayhome --autocor 1.0 --no_post_inter

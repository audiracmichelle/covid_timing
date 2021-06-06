Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 1.0 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 1.0 --no_post_inter


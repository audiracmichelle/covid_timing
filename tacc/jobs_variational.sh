Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type decrease --autocor 0.7 
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type decrease --autocor 0.9
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.7 
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.9
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type decrease --autocor 0.7 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type decrease --autocor 0.9 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.7 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.9 --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.9 --exclude_ny
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.9 --exclude_ny --no_post_inter
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.7 --exclude_ny
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --type stayhome --autocor 0.7 --exclude_ny --no_post_inter
s
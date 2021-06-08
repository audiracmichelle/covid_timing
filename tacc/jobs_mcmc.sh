Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention decrease --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention decrease --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention stayhome --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention decrease --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention decrease --autocor 1.0 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds" --intervention stayhome --autocor 1.0 --no_post_inter

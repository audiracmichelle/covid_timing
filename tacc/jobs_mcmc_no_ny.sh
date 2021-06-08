Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --intervention stayhome --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --intervention stayhome --autocor 1.0 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --exclude_post_only --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --exclude_post_only --intervention stayhome --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --exclude_post_only --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc_no_ny/$LAUNCHER_TSK_ID" --init "results/vb_no_ny/$LAUNCHER_TSK_ID/fit.rds" --exclude_ny --exclude_post_only --intervention stayhome --autocor 1.0 --no_post_inter

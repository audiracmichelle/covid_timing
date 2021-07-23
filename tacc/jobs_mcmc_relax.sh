Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention decrease --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention decrease --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention stayhome --autocor 0.8 
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention stayhome --autocor 1.0
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention decrease --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention decrease --autocor 1.0 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention stayhome --autocor 0.8 --no_post_inter
Rscript --vanilla fit_mcmc.R --dir "mcmc/$LAUNCHER_TSK_ID" --init "results/vb/$LAUNCHER_TSK_ID/fit.rds"  --ar_relax_prior_scale --intervention stayhome --autocor 1.0 --no_post_inter

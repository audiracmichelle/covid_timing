Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 0.7
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 0.7 
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 1.0
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 0.7 --spatial_scale 0.75
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 1.0 --spatial_scale 0.75
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 0.7 --spatial_scale 0.75
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 1.0 --spatial_scale 0.75
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 0.7 --spatial_scale 0.5
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention decrease --autocor 1.0 --spatial_scale 0.5
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 0.7 --spatial_scale 0.5
Rscript --vanilla fit_variational.R --dir "vb/$LAUNCHER_TSK_ID" --intervention stayhome --autocor 1.0 --spatial_scale 0.5

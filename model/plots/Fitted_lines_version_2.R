#MS FIG 1

	postscript("plots/fig3-2.eps", horizontal=FALSE, onefile=FALSE, height=6,width=6,pointsize=10, paper ="special" )
	#LOAD MAMMAL DATA
	X=obs.mam$Wa; Y=obs.mam$W0;
	Y_ax=seq.log(10^-9, 10^3, 10);	X_ax=seq.log(10^-7, 10^5, 10);
	Y_lab =c(expression(10^-9), expression(10^-8), expression(10^-7), expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3))
	X_lab =c(Y_lab[3:13], expression(10^4), expression(10^5))

	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="darkgray", pch = 19, cex=1.2, xaxs="i", yaxs="i", family="helvetica")
		points(X, Y, type="p",pch = 1, cex=1.3,lwd =0.25)
		
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext(expression(paste("Offspring size (kg) ", " [log scale]")), side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)
		mtext(expression(paste("Adult size (kg) ", " [log scale]")), side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)


	#------------------------------------------------------------
	#ADD FITTED CURVES - ver 2 --> use different values of D
	c_LRI = 0.927/2;

	#mammals
  D=1; b_LRI = 0.9519; s = 0.944;
  points(X, ESS_analytic(X, LRI(X)), type="l", lty="solid", lwd=2, col ="darkgray")

  #LOAD PLANT DATA
	X=obs.plt$Wa; Y=obs.plt$W0;
 	points(X, Y, type="p", pch=1, cex=1.3, lwd =0.25)

  D=2; b_LRI = 0.632; s = 0.00677/17;
  points(X, ESS_analytic(X, LRI(X)), type="l", lty="solid", lwd=2)

	#1:1
	points(c(10^-7, 10^5), c(10^-7, 10^5), type="l", lty="dashed", lwd=1.5);
	

	dev.off()

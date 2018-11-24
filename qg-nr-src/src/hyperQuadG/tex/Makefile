.IGNORE:

TARGET = paper

#all: dvi
all: klk
#all: dvi pdf

links:
 	#ln -s figures/optdepthBNS/output.ps  optdepthBNS.ps
	#ln -s figures/optdepthNS/output.ps   optdepthNS.ps
	#ln -s figures/bnshshenlums/output.ps bnshshenlums.ps
	#ln -s figures/hotdenserotatingmagnetizedstar/final5x2.png hotdense2D.png
	#convert figures/hotdenserotatingmagnetizedstar/final5x2.png eps3:hotdense2D.eps
	#convert figures/hotdenserotatingmagnetizedstar/oned.png     hotdense1D.eps
	#ln -s figures/hotdenserotatingmagnetizedstar/oned.eps        hotdense1D.eps
	#convert figures/bns2dplots/bns.png                          bns.eps
	#ln -s figures/bns2dplots/bns.eps                             bns.eps
	#ln -s figures/bns2dplots/bns2.eps                             bns2.eps
	#ln -s figures/bns2dplots/bnsleft.eps                             bnsleft.eps
	#ln -s figures/bns2dplots/bnsright.eps                             bnsright.eps
	#ln -s figures/hotdenserotatingmagnetizedstar/final5x2.png hotdense2D.png
	#ln -s figures/hotdenserotatingmagnetizedstar/oned.png hotdense1D.png

klk: 
	pdflatex     $(TARGET)
	bibtex    $(TARGET)
	pdflatex     $(TARGET)
	pdflatex     $(TARGET)

units:
	pdflatex units
	pdflatex units

dvi:  links
	latex     $(TARGET)
	bibtex    $(TARGET)
	latex     $(TARGET)
	latex     $(TARGET)

pdf: dvi
	dvipdf     -sPAPERSIZE=letter $(TARGET)
#	dvipdf     -sPAPERSIZE=a4 $(TARGET)
	#dvipdf     $(TARGET)
#	dvips -tletter -o $(TARGET).ps $(TARGET).dvi
#	ps2pdf $(TARGET).ps

ps: dvi
	dvips     $(TARGET)

clean:
	/bin/rm -f *.dvi *.log *.aux *.toc *.bbl *blg *.end  *~ paperNotes.bib
	/bin/rm -f $(TARGET).ps $(TARGET).pdf  $(TARGET).tex~ Makefile~ $(TARGET).out $(TARGET).bib~

arxiv: dvi
	egrep -i "epsfig|includegraph" $(TARGET).tex | grep -v "^%"
	tar czvhf arxiv.tgz aas_macros.tex	$(TARGET).tex $(TARGET).bbl optdepthNS.ps optdepthBNS.ps ./figures/rhoc_lsk220_cowling.eps ./figures/fourier_rhoc_lsk220_cowling.eps ./figures/rhoBc_lsk220_dynamical.eps ./figures/fourier_rhoBc_lsk220_dynamical.eps ./figures/rhoBc_gshen_hot_cowling.eps ./figures/fourier_rhoBLc_gshen_hot_cowling.eps ./figures/luminosities_gshen_hot_cowling.eps ./figures/luminosities_gshen_hot_cowling2.eps hotdense1D.eps hotdense2D.eps bnshshenlums.ps bnsleft.eps bnsright.eps



prdresubmit: dvi
	egrep -i "epsfig|includegraph" $(TARGET).tex | grep -v "^%"
	tar czvhf prdresubmit.tgz $(TARGET).tex $(TARGET).bbl ./figures/orbits_rAvg.ps ./figures/GW_psi4/c2p2.ps ./figures/Bflds_lines/UD_s17.ps ./figures/Bflds_Bcomp/bns03-Bflds_wcolours--06.ps ./figures/Bflds_Bcomp/bns03-Bflds_wcolours--11.ps ./figures/Bflds_Bcomp/bns03-Bflds_wcolours--16.ps ./figures/Bflds_Bcomp/bns03-Bflds_wcolours--20.ps ./figures/Bflds_xi/bns03-Bfields_xi--12.ps ./figures/Bflds_lines/UU_s17.ps ./figures/Bflds_Bcomp/bns01-Bflds_wcolours--06.ps ./figures/Bflds_Bcomp/bns01-Bflds_wcolours--11.ps ./figures/Bflds_Bcomp/bns01-Bflds_wcolours--16.ps ./figures/Bflds_Bcomp/bns01-Bflds_wcolours--20.ps ./figures/Bflds_xi/bns01-Bfields_xi--12.ps ./figures/Bflds_lines/U0_s17.ps ./figures/Bflds_Bcomp/bns02-Bflds_wcolours--06.ps ./figures/Bflds_Bcomp/bns02-Bflds_wcolours--12.ps ./figures/Bflds_Bcomp/bns02-Bflds_wcolours--16.ps ./figures/Bflds_Bcomp/bns02-Bflds_wcolours--20.ps ./figures/Bflds_xi/bns02-Bfields-xi_16.ps ./figures/currents/bns03-currents_06.ps ./figures/currents/bns01-currents_06.ps ./figures/currents/bns02-currents_06.ps ./figures/sqphi2/bns03_UD-sqphi2_3d-14.ps ./figures/sqphi2/sphere_sqphi2_UD.ps ./figures/sqphi2/bns01_UU-sqphi2_3d-14.ps ./figures/sqphi2/sphere_sqphi2_UU.ps ./figures/sqphi2/bns02_Uu-sqphi2_3d-14.ps ./figures/sqphi2/sphere_sqphi2_Uu0.ps ./figures/EM_phi2/luminosity.ps ./figures/UD-Bz-XYplane_1600--z_5.ps ./figures/UU-Bz-XYplane_1600--z_5.ps ./figures/U0-Bz-XYplane_1600--z_5.ps aas_macros.tex ./figures/Bflds_xi/bns03-Bfields-xi_12_onBb.ps ./figures/Bflds_xi/bns01-Bfields-xi_12_onBb.ps ./figures/Bflds_xi/bns02-Bfields-xi_16_onBb.ps


arxivresubmit: dvi
	egrep -i "epsfig|includegraph" $(TARGET).tex | grep -v "^%"
	tar czvhf arxiv2.tgz $(TARGET).tex $(TARGET).bbl  ./figures/Bflds_xi/bns03-Bfields-xi_12_onBb.ps ./figures/Bflds_xi/bns01-Bfields-xi_12_onBb.ps ./figures/Bflds_xi/bns02-Bfields-xi_16_onBb.ps

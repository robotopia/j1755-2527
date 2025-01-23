include ../.env
export

RMS = $(shell seq -970 0.5 -950)
YAMLS = $(shell ls *-I.yaml)

DM = 710
LIGHTCURVES = $(patsubst %.yaml,%_lightcurve_DM$(DM).txt,$(YAMLS))

lightcurves: $(LIGHTCURVES)

%_lightcurve_DM$(DM).txt: %.yaml
	python dedisperse_dynspec.py --yaml $< --lightcurve $@ --dms $(DM) --dynspec_image $*_DM$(DM).png

#1410773120_and_1410773416-RM_%-fullpol.png:
#	python polprof.py --rel_time --off_pulse_lims 0:200 --DM 1100 --RM $* --outdynspec 1410773120_and_1410773416-dynspec-RM_$*.png --outfile 1410773120_and_1410773416-RM_$*-fullpol.png 1410773120-?_and_1410773416-?.yaml


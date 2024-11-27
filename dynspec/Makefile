include ../.env
export

RMS = $(shell seq -970 0.5 -950)
PNGS = $(patsubst %,1410773120_and_1410773416-RM_%-fullpol.png,$(RMS))
YAMLS = 1410773120-I_and_1410773416-I.yaml \
				1410777264-I_and_1410777560-I.yaml \
				1410785848-I.yaml \
				1410789992-I.yaml \
				1412506568-I.yaml \
				1412510752-I.yaml \
				1412594488-I.yaml \
				1412598672-I.yaml \
				1412678224-I.yaml \
				1412682408-I.yaml \
				1412766144-I.yaml \
				1412770328-I.yaml \
				1412849880-I.yaml \
				1412854064-I.yaml \
				1412937800-I.yaml \
				1412941984-I.yaml \
				1413109456-I.yaml \
				1413113640-I.yaml \
				1413201560-I.yaml \
				1413369032-I.yaml \
				1413373216-I.yaml \
				1413381294-I.yaml

DM = 1100
LIGHTCURVES = $(patsubst %.yaml,%_lightcurve_DM$(DM).txt,$(YAMLS))

lightcurves: $(LIGHTCURVES)

%_lightcurve_DM$(DM).txt: %.yaml
	dedisperse_dynspec --yaml $< --lightcurve $@ --dms $(DM) --dynspec_image $*_DM$(DM).png

rm-attempt: $(PNGS)

1410773120_and_1410773416-RM_%-fullpol.png:
	python polprof.py --rel_time --off_pulse_lims 0:200 --DM 1100 --RM $* --outdynspec 1410773120_and_1410773416-dynspec-RM_$*.png --outfile 1410773120_and_1410773416-RM_$*-fullpol.png 1410773120-?_and_1410773416-?.yaml

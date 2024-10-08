RMS = $(shell seq 900 0.5 1100)
PNGS = $(patsubst %,1410773120_and_1410773416-RM_%-fullpol.png,$(RMS))


all: $(PNGS)

1410773120_and_1410773416-RM_%-fullpol.png:
	python polprof.py --rel_time --off_pulse_lims 0:200 --DM 1100 --RM $* --outdynspec 1410773120_and_1410773416-dynspec-RM_$*.png --outfile 1410773120_and_1410773416-RM_$*-fullpol.png 1410773120-?_and_1410773416-?.yaml

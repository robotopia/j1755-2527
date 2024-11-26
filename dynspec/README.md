# Data processing notes

Where (some of) the data came from (to the best of my knowledge):
## 1413381294

### Pull the data from Acacia

```
rclone copy mwasci:askapj1755/J1755-25_1729341386_scan18_images_models.tar.gz .
tar -xzf J1755-25_1729341386_scan18_images_models.tar.gz \
  J1755-25_1729341386_scan18-I.csv \
  J1755-25_1729341386_scan18-Q.csv \
  J1755-25_1729341386_scan18-U.csv \
  J1755-25_1729341386_scan18-V.csv \
  J1755-25_1729341386_scan18-I.yaml \
  J1755-25_1729341386_scan18-Q.yaml \
  J1755-25_1729341386_scan18-U.yaml \
  J1755-25_1729341386_scan18-V.yaml
```

### Rename the files to be consistent with the naming convention

```
mv J1755-25_1729341386_scan18-I.yaml 1413381294-I.yaml
mv J1755-25_1729341386_scan18-Q.yaml 1413381294-Q.yaml
mv J1755-25_1729341386_scan18-U.yaml 1413381294-U.yaml
mv J1755-25_1729341386_scan18-V.yaml 1413381294-V.yaml
mv J1755-25_1729341386_scan18-I.csv 1413381294-I.csv
mv J1755-25_1729341386_scan18-Q.csv 1413381294-Q.csv
mv J1755-25_1729341386_scan18-U.csv 1413381294-U.csv
mv J1755-25_1729341386_scan18-V.csv 1413381294-V.csv
```

I had to manually edit the YAML file to fix some typos, update it with the new file names, correct the "transpose" option for this data set, and also flag one of the time bins at the end.
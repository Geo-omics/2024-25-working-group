# Extract cyanobacterial index (CI) values from NOAA Lake Erie Harmful Algal Bloom Forecast

For times and locations of interest. Single and bulk input options.

## Extract Cyanobacteria Index (CI) Values Predicted by HABTracker
**Usage:**
- habs_grab_tracker_ci_v1.R -l LATITUDE -L LONGITUDE -d DATE [options]
- habs_grab_tracker_ci_v1.R --csv CSV [options]

## Options
- `-l, --lat=LATITUDE`  
  Latitude in decimal degrees.
  
- `-L, --lon=LONGITUDE`  
  Longitude in decimal degrees.
  
- `-d, --date=DATE`  
  Date of sample collection (format: YYYY-MM-DD).
  
- `-t, --time=TIME`  
  Time of sample collection in Eastern Daylight Time (EDT) (format: hh:mm).  
  DEFAULT = 12:00.
  
- `-r, --rad=DISTANCE`  
  Radial distance in meters to calculate average CI values.  
  DEFAULT = nearest CI value.
  
- `-c, --csv=CSV`  
  CSV input file path (for batch mode).
  
- `-o, --out=OUTPUT`  
  Output file.
  
- `-h, --help`  
  Show this help screen.

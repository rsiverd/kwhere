# kwhere
A small utility that tells you which KELT field(s) contain the specified RA, Dec. Also makes nifty plots (requires matplotlib's basemap).

Installation
============

Try the following using your normal Python environment (you may have 
dependencies satisfied already):
```
$ git clone https://github.com/rsiverd/kwhere.git

$ cd kwhere
$ ./kwhere.py 123.45 54.32 -P derp.png
Basemap not installed ... plots will be ugly!
Found 1 KELT field(s) containing the test point:
KN19 --- (120.750, +57.000)
Saving diagram to derp.png ... done.

```
If your default Python already has the necessary packages installed, you will
see something similar to the output above. If the script does not run because
some dependencies are missing, try the following:
```
$ pip install --upgrade pip
$ pip install -r requirements.txt
```

This should pull in everything needed to run the script. Tested working with
a Python2.7 virtual environment on CentOS 7.


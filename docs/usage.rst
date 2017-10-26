=====
Usage
=====

Two scripts are included in MoBamPlot in the mobamplot folder::

    1. getcountvectordata.py
    2. mobamplot.py
    
GetCountVectorData
------------------

The first script, ``getcountvectordata.py``, is used to create an .hdf5 file that stores counts read from a bam file. Each region of interest is stored in a dset that saves its information such as the ID label, chromosome, segment start coordinate, and segment end coordinate.

There are two required inputs to ``getcountvectordata.py``:

     --bedfile  The full path of the desired bed file. This does NOT use a flag.
     --bamfiles  The full path of the desired bam file. Multiple bam files are separated by spaces. The folder containing the bam file must also include an indexed file for that bam. This DOES require a flag. 

and some optional arguments:

     --outfolder  Indicates where the .hdf5 file will be saved. If none is specified, it is saved in the current working directory.
     --bamIDs  Desired name to identify each BAM.

This is an example call made from the folder containing the scripts in the terminal:

.. code-block:: shell

   $ python getcountvectordata.py /path/bedfile.bed --bamfiles /path/bam1.bam /path/to/bams/bam2.bam --outfolder /path/outputfolder

The output is a .hdf5 file with the countdata in the output folder. The .hdf5 file is currently named after the bed file name.


MoBamPlot
---------

The second script, ``mobamplot.py``, takes the .hdf5 counts file and creates a visualization of the data. The data can be graphed as an interactive histogram html or as a static line plot.

These are the following arguments for ``mobamplot.py``:

     --hdf5  Input hdf5 file containing the counts. Include path.
     --outfolder  Input desired folder to save figures to. If none is specified, figures will be saved to the current working directory.
     --binsize  Desired bin size. This value is the size of an individual bin. The default is 200.
     --line  Changes graph type from histogram, the default, to line.
     --viewplot  Once a graph is created, a web browser will automatically open to the html of the graph.
     --format  Indicates the file type for a static image to be saved as. Options are png, svg, jpeg, and webp.

This is an example call made from the folder containing the scripts in the terminal that creates an html figure of the data wit a bin size of 50 saved in outputfolder:

.. code-block:: shell

   $ python mobamplot.py --hdf5 /path/counts.hdf5 --outfolder /path/figurefolder --binsize 50

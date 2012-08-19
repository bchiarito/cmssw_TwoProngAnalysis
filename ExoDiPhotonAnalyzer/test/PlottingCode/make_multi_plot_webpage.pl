#!/usr/bin/perl -w

# A list of figures is passed as command line argument
# eg ./make_multi_plot_webpage.pl *.png
# output is html code for a webpage where these plots are displayed
# as a table, with thumbnails and links to the full image

# the output html is just printed to screen
# so it may be more useful to redirect output to a file on the command line
# ie ./make_multi_plot_webpage.pl *.png > index.html

# options are used to distinguish 'dataFlag', date and lumi of sample
# so usage is:
# ./make_multi_plot_webpage.pl --dataFlag="Data" --date="25 July 2012" --lumi=""5.4 fb-1" *.png

use Getopt::Long;

#dataFlag would be "MC, "Data", "Data+Bkg", "Signal", etc
#date is taken as a string to display literally, eg "25 July 2012"
#lumi is also treated as string, eg "5.4 fb-1"
GetOptions("dataFlag=s"=>\$dataFlag, 
	   "date=s"=>\$date,
	   "lumi=s"=>\$lumi
);

# boiler plate header for html page
print "<html> \n";
print "<head>\n";

print "</head>\n";


print "<body>\n";

#print any headings here 
print "<h1> \n";
print "Diphoton Plots";
print "</h1> \n";

# insert labels like Data/MC, date and lumi?
print "<h2>\n";
print "$dataFlag \t $date \t $lumi";
print "</h2>\n";



# now start table
print "<table border=\"0\" cellpadding=\"0\" cellspacing=\"8\">\n";

# list of files to process passed on command line
# note that the file format is restricted to png here
# if instead gif or jpg, need to edit below

# for the first row
print "<tr>";
$count = 0;

foreach (@ARGV) {
    if(/^(.*)\.png/) {
        $file = $_ ; # ie the current file name from list
        $thumb = "${1}_thumb.png";
#make the thumbnail
	`convert $file -resize 240x240 $thumb`;

#how many columns? 3?
	if($count % 3 == 0) {
          #start a new row
	    print "<tr>";

	}

# for each photo, display the thumbnail, with link to full photo
	print "<td><a href=\"$file\"><img src=\"$thumb\"></a></td>\n";

#increment count
	$count++;



    }
} # end of file list loop

print "</table>\n";

# boiler plate to close html document
print "</body>\n";
print "</html> \n";





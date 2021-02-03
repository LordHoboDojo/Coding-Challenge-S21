Libraries Used:
Image - Used to create images and was used to make the final product from the data
ImageDraw - Used to draw/write on the images. This library was used to draw the circles, make the legend, and map the features on the circular genome map
ImageFont - This part of the Image library was used to adjust the font of the text. 
SeqIO from biopython - This library was used to read the genbank file and map the features onto the genome map
Math - The math library was needed because I had to utilize pi as well as well as a few trigonometric functions.

Graph Information:
The graph consists of two circles and the areas shaded on the circles represent various features on the genome. The inner circle represents the complementary strand and the outer circle represents the template strand. The shaded regions are color coded and each color represents a different feature on the genome. The key shows which colors correspond to each feature. Each of the lines surrounding the circle are separated by 50 base pairs with the starting position at 6:00.

Sources:
https://biopython.org/wiki/SeqIO - This is the documentation for biopython, one of the libraries I am using. I also looked at their sample code in order to learn how to read a genbank file.

https://pillow.readthedocs.io/en/stable/reference/ImageDraw.html#PIL.ImageDraw.PIL.ImageDraw.ImageDraw.ellipse - Documentation for ImageDraw. I used the examples to model my solution

https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html - This site has a sample genbank file and explains the schema for the data. I used it to get an understanding of how the genbank file data is to be read.







from Bio import SeqIO
from PIL import Image, ImageDraw, ImageFont
import math


def get_line_endpoints(sequence_length: int) -> list:
    angle_increment = 2 * math.pi / sequence_length
    angle = 0
    line_points = []
    # Loops through points on the circle and stores them for future use
    for x in range(0, sequence_length):
        point1 = (1950 * math.sin(angle) + 2000, 1950 * math.cos(angle) + 2000)
        point2 = (1900 * math.sin(angle) + 2000, 1900 * math.cos(angle) + 2000)
        point3 = (1850 * math.sin(angle) + 2000, 1850 * math.cos(angle) + 2000)
        three_point_list = [point1, point2, point3]
        line_points.append(three_point_list)
        angle += angle_increment
    return line_points


# Initialize the image file and ImageDraw variable and get data from genbank file
records = list(SeqIO.parse("Genome.gb", "genbank"))
out = Image.new("RGB", (6000, 5000), (255, 255, 255))
draw = ImageDraw.Draw(out)

# Draw circles
twoPointList = [(100, 100), (3900, 3900)]
draw.ellipse(twoPointList, fill=None, outline=5, width=5)
twoPointList = [(50, 50), (3950, 3950)]
draw.ellipse(twoPointList, fill=None, outline=5, width=5)
twoPointList = [(150, 150), (3850, 3850)]
draw.ellipse(twoPointList, fill=None, outline=5, width=5)

# Draw main title
font = ImageFont.truetype("Lato.ttf", 100)
text = "Circular Genome Map of Tomato Stunt Curly Virus"
w, h = font.getsize(text)
draw.text(((4000 - w) / 2, (4000 - h) / 2), text=text, fill="black", font=font)

# Create key
key_text = "Key"
draw.text((5000, 2500), text=key_text, fill="black", font=font)
x = 600
draw.rectangle((4300 + x, 2900, 4400 + x, 3000), outline=(255, 0, 0), fill=(125, 100, 121))
draw.rectangle((4300 + x, 3100, 4400 + x, 3200), outline=(255, 0, 0), fill=(171, 151, 135))
draw.rectangle((4300 + x, 3300, 4400 + x, 3400), outline=(255, 0, 0), fill=(44, 28, 17))
draw.rectangle((4300 + x, 3500, 4400 + x, 3600), outline=(255, 0, 0), fill=(153, 121, 30))
draw.rectangle((4300 + x, 3700, 4400 + x, 3800), outline=(255, 0, 0), fill=(95, 105, 38))

draw.line([(4450 + x, 2950), (4500 + x, 2950)], fill="black", width=20)
draw.line([(4450 + x, 3150), (4500 + x, 3150)], fill="black", width=20)
draw.line([(4450 + x, 3350), (4500 + x, 3350)], fill="black", width=20)
draw.line([(4450 + x, 3550), (4500 + x, 3550)], fill="black", width=20)
draw.line([(4450 + x, 3750), (4500 + x, 3750)], fill="black", width=20)

draw.text((4550 + x, 2890), fill="black", font=font, text="v1")
draw.text((4550 + x, 3090), fill="black", font=font, text="v2")
draw.text((4550 + x, 3290), fill="black", font=font, text="c1")
draw.text((4550 + x, 3490), fill="black", font=font, text="c2")
draw.text((4550 + x, 3690), fill="black", font=font, text="c3")

# Get points that correspond to the same angle on each of the three circles created
line_pts = get_line_endpoints(sequence_length=len(records[0]))
# Draw lines to show position of the genes
for i in range(0, len(line_pts)):
    if i % 50 == 0:
        x_offset = 30 * math.sin(2 * math.pi * i / len(line_pts))
        y_offset = 30 * math.cos(2 * math.pi * i / len(line_pts))
        draw.line([line_pts[i][0], (line_pts[i][0][0] + x_offset, line_pts[i][0][1] + y_offset)], fill="black", width=8)

seq_length = len(records[0])
for feature in records[0].features:
    if feature.type == "gene":
        if feature.strand == 1:
            if feature.qualifiers['gene'][0] == "v1":
                for x in range(feature.location.start, feature.location.end):
                    draw.line((line_pts[x][0], line_pts[x][1]), fill=(125, 100, 121), width=7)
            if feature.qualifiers['gene'][0] == "v2":
                for x in range(feature.location.start, feature.location.end):
                    draw.line((line_pts[x][0], line_pts[x][1]), fill=(171, 151, 135), width=7)
        if feature.strand == -1:
            if feature.qualifiers['gene'][0] == "c1":
                for x in range(seq_length - feature.location.end, seq_length - feature.location.start):
                    draw.line((line_pts[x][1], line_pts[x][2]), fill=(44, 28, 17), width=7)
            if feature.qualifiers['gene'][0] == "c2":
                for x in range(seq_length - feature.location.end, seq_length - feature.location.start):
                    draw.line((line_pts[x][1], line_pts[x][2]), fill=(153, 121, 30), width=7)
            if feature.qualifiers['gene'][0] == "c3":
                for x in range(1, 600):
                    draw.line((line_pts[x][1], line_pts[x][2]), fill=(95, 105, 38), width=7)

out.save(fp="result.png")

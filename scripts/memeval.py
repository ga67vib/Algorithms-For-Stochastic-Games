import csv

def min_max_diff(row):
	min = 10000000
	max = 0
	for entry in row:
		if entry == "X":
			continue
		entry = int(entry)
		if entry < min:
			min = entry
		if entry > max:
			max = entry
	diff = max/min
	return (min,max,diff)


with open("mem1/mem.csv", "r") as memfile:
	with open("mem1/mem-relative.csv", "w") as outfile:
		memreader = csv.reader(memfile, delimiter=',')
		
		header = next(memreader)
		for entry in header:
			outfile.write(entry+ ",")
		outfile.write("\n")


		for row in memreader:
			min = min_max_diff(row[2:])[0]

			outfile.write(row[0] + "," + row[1] + ",")

			for entry in row[2:]:
				if entry == "X":
					outfile.write("X")
				else:
					outfile.write(str(int(int(entry)/min)))
				outfile.write(",")
			outfile.write("\n")

							

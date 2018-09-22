#strain_file = open("strainfile.txt","r").read()
#file = file.lower()
# file = file.replace(".","")
# file = file.replace(" ","")
# file = file.split("\n")

# names = ["coli", "subtilis", "crescentus", "fischeri", "cerevisiae", "pastoris", "bassiana", "tumefa"]
# counts = []
# for i in range(len(names)):
#   counts.append(0)
#
# for line in file:
#   for i in range(len(names)):
#     if names[i] in line:
#       counts[i] += 1

# print(counts)
# print(names)

strain_file = open("strainfile.txt","r").read()

words_file = open("words.txt", "r").read()
words = words_file.split("\n")

possible_strains = []
lines = strain_file.split("*****")
print(lines)
for line in lines:
        if line not in words:
            possible_strains += line
print(possible_strains)

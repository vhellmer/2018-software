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

from collections import defaultdict

strain_file = open("strainfile-Reformatted.txt","r").read()

words_file = open("words.txt", "r").read()
common_words = words_file.split("\n")

possible_strains = []
lines = strain_file.split("*****")
for line in lines:
    line = line.replace("\n", " ")
    line = line.replace(",", " ")
    line = line.replace(".", " ")
    line = line.replace("(", " ")
    line = line.replace(":", " ")
    line = line.replace(")", " ")

    words = line.split(" ")
    team_words = []
    for word in words:
        if word.lower() not in common_words:
            team_words.append(word.lower())
    possible_strains.append(team_words)

print(possible_strains)

common_names = ["coli", "subtilis", "crescentus", "fischeri", "cerevisiae", "pastoris", "bassiana", "tumefa"]
common_dict = defaultdict(int)
all_dict = defaultdict(int)

for team_list in possible_strains:
    for word in set(team_list):
        for name in common_names:
            if name in word:
                common_dict[name] += 1
        all_dict[word] += 1

print(common_dict)
for key, value in sorted(all_dict.items()):
    print(key, value)
print("NUMBER OF TEAMS: " + str(len(possible_strains)))
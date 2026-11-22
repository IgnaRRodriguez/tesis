import csv
import os
import datetime
import math
from collections import deque

class Node(object):
	"""docstring for Node"""
	def __init__(self,):
		self.families = {}
		self.nexts= [None] * 26
		
def getTrie():
	return Node()


def addProtein(node:Node,protein:str,family:str,affectedProteins:int):
	currentNode = node
	for letter in protein:
		position= ord(letter) - 65
		if currentNode.nexts[position] is None :
			currentNode.nexts[position]= Node()
		currentNode = currentNode.nexts[position]
		if family not in currentNode.families.keys(): 
			currentNode.families[family] = affectedProteins
		elif currentNode.families[family]< affectedProteins:
			currentNode.families[family] = affectedProteins

def getDistanceFamily(trie:Node,family1:str,family2:str,totalNodes:dict):
	nextTries = deque()
	similarty = 0
	for nexts in filter(None,trie.nexts):
		if ((family1 in nexts.families.keys()) and (family2 in nexts.families.keys())):
			nextTries.append(nexts)
	level = 0
	while nextTries:
		currentTries = nextTries 
		nextTries = deque()
		level += 1
		levelSimilarity = 0
		while currentTries:
			currentNode = currentTries.popleft()
			nodeScore = 1
			for family in currentNode.families.keys():
				nodeScore *= (totalNodes[level][family]/totalNodes[level]['total'])
			levelSimilarity += nodeScore
			for nexts in filter(None,currentNode.nexts):
				if ((family1 in nexts.families.keys()) and (family2 in nexts.families.keys())):
					nextTries.append(nexts)
		similarty += (level*levelSimilarity)
	return similarty
	


def pruneTrie(trie:Node,totalNodes:dict):
	level = 0
	for i in range(len(trie.nexts)):
		nextNode = trie.nexts[i]
		if nextNode is not None:
			if len(nextNode.families.keys()) <= 1:
				trie.nexts[i] = None

	nextLevel = deque(trie.nexts)
	while nextLevel:
		level += 1
		totalNodes[level] = {}
		totalNodes[level]['total'] = 0
		currentLevel = nextLevel
		nextLevel= deque()
		while currentLevel:
			currentNode = currentLevel.popleft()
			if currentNode is not None:
				totalNodes[level]['total'] += 1
				for i in range(len(currentNode.nexts)):
					if (currentNode.nexts[i] is not None):
						nextNode = currentNode.nexts[i]
						if len(nextNode.families.keys()) <= 1:
							currentNode.nexts[i] = None
				for family in currentNode.families.keys():
					if family in totalNodes[level]:
						totalNodes[level][family] += 1
					else:
						totalNodes[level][family] = 1
				nextLevel.extend(currentNode.nexts)


print("start")
trie= getTrie()
dirList = list(filter(lambda x: (x.endswith("NE.csv")),os.listdir("../trie/output1")))
families = []
totalNodesByLevel = {}
proteinsInFamily = {
"G_ABCtran" : 27442,
"C_ABCtranSCRAMBLED" : 27442,
"C_ABCtranUNIFORM" : 27442,
"R_Annexin" : 3424,
"R_Arm" : 19821,
"R_Collagen" : 14747,
"G_CorA" : 13398,
"R_CWbinding1" : 5403,
"R_Fer4" : 30523,
"R_Filamin" : 4137,
"G_GDCP" : 8371,
"G_Globin" : 6606,
"G_Glycohydro19" : 3363,
"G_GreAGreB" : 8643,
"R_HEAT" : 25339,
"G_HelicaseC" : 37928,
"R_HemolysinCabind" : 12515,
"R_Hexapep" : 37258,
"R_Kelch1" : 12700,
"R_Ldlrecepta" : 12312,
"R_Ldlreceptb" : 4233,
"R_LRR1" : 33060,
"G_MgtE" : 5518,
"G_MIP" : 11824,
"R_Mitocarr" : 26363,
"C_mix" : 70859,
"C_mixScrambled" : 70859,
"C_mixUniformDistributed" : 70866,
"R_MORN" : 8079,
"R_Nebulin" : 1062,
"R_NEWAnk" : 32169,
"C_NEWAnkSCRAMBLED" : 32169,
"C_NEWAnkTPR1" : 68542,
"C_NEWAnkUNIFORM" : 32169,
"C_NEWAnkWithoutP25963" : 32168,
"R_NEWWD40" : 30898,
"G_PBP" : 6734,
"R_PD40" : 32370,
"R_Pentapeptide" : 15325,
"G_PeptidaseC25" : 1392,
"G_PFLlike" : 3068,
"G_PIN" : 21705,
"G_Pkinase" : 27797,
"G_PPbinding" : 4230,
"R_PPR" : 38342,
"G_PUD" : 924,
"R_PUF" : 3956,
"G_Rhomboid" : 15532,
"R_Sel1" : 21424,
"G_TerB" : 6952,
"G_Thaumatin" : 2827,
"R_TPR1" : 36389,
"R_TSP1" : 11296,
"G_ValtRNAsyntC" : 4782,
"R_YadAhead" : 3933
} 
startTime=datetime.datetime.now()
print("get trie")
for family in dirList:
	with open('../trie/output1/'+family,mode='r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		name= family.split("_")[0]+"_"+family.split("_")[1]
		families.append(name)
		for row in csv_reader:
			#Pattern,PatternLength,Instances,AffectedProteins
			if ((int(row["PatternLength"])<25) and ((int(row["AffectedProteins"])/proteinsInFamily[name])>=0.05)):
				addProtein(trie,row["Pattern"],name,int(row["AffectedProteins"]))

pruneTrie(trie,totalNodesByLevel)
print(totalNodesByLevel)

distances= [([0]* len(families)) for i in range(len(families))]
for i in range(len(families)):
	family1 = families[i]
	print(family1)
	for j in range(i,len(families)):
		family2 = families[j]
		print(family2)
		val = getDistanceFamily(trie,family1,family2,totalNodesByLevel)
		distances[i][j] = val
		distances[j][i] = distances[i][j]
	print("-------------------------------------------------")


endTime=datetime.datetime.now()
print('time: ',abs(endTime-startTime).total_seconds() / 60.0)
with open('resultV3-proteinasAfectadas-05.csv','w',newline='') as csv_file:
	csvWriter = csv.writer(csv_file)
	csvWriter.writerow(families)
	csvWriter.writerows(distances)
